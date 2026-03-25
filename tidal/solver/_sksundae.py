"""Typed wrappers for sksundae CVODE and IDA solvers.

All sksundae imports are isolated here so that ``# pyright: ignore``
suppressions for ``reportMissingTypeStubs`` are confined to this module.
The rest of the codebase imports ``call_cvode`` / ``call_ida`` and receives
fully typed :class:`SundialsResult` objects.

Stepwise variants (``call_cvode_stepwise``, ``call_ida_stepwise``) use
the ``init_step()`` + ``step(t)`` API for progress reporting with zero
overhead on the solver's internal computation.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

# sksundae >= 1.1.2 emits a UserWarning when both ``sparsity`` and ``jacfn``
# are provided: "Sparse Jacobian approximation will be ignored in favor of
# 'jacfn'".  This is expected — we *want* the analytical jacfn to take
# precedence over FD.  Suppress globally for this module.
warnings.filterwarnings(
    "ignore",
    message="Sparse Jacobian approximation will be ignored",
    category=UserWarning,
    module=r"sksundae\.",
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.solver.progress import SimulationProgress


@dataclass(frozen=True)
class SundialsResult:
    """Typed result from ``CVODE.solve()`` or ``IDA.solve()``.

    Both CVODE and IDA return result objects with identical attribute names
    and the same array shape convention (``n_times`` x ``n_states``).

    IDA additionally provides ``yp`` (time-derivative vector) which is used
    to extract constraint field velocities for energy measurement.
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]  # shape (n_times, n_states) — TIDAL convention
    success: bool
    message: str
    yp: NDArray[np.float64] | None = None  # IDA only: shape (n_times, n_states)


def call_cvode(
    rhsfn: Any,  # noqa: ANN401
    t_eval: NDArray[np.float64],
    y0: NDArray[np.float64],
    **options: Any,  # noqa: ANN401
) -> SundialsResult:
    """Create a CVODE solver and run ``.solve(t_eval, y0)``.

    Parameters
    ----------
    rhsfn:
        RHS function with signature ``f(t, y, yp) -> None`` (writes into
        *yp* in-place).
    t_eval:
        Output time points.
    y0:
        Initial state vector.
    **options:
        Forwarded to ``CVODE.__init__``.  Recognized keys: ``method``,
        ``rtol``, ``atol``, ``max_num_steps``, ``max_step``, ``linsolver``,
        ``sparsity``.

    Returns
    -------
    SundialsResult
        Typed result container.
    """
    from sksundae.cvode import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        CVODE,
    )

    solver = CVODE(rhsfn, **options)
    raw = solver.solve(t_eval, y0)
    return SundialsResult(
        t=np.asarray(raw.t, dtype=np.float64),  # pyright: ignore[reportUnknownArgumentType,reportUnknownMemberType]
        y=np.asarray(raw.y, dtype=np.float64),  # pyright: ignore[reportUnknownArgumentType,reportUnknownMemberType]
        success=bool(raw.success),  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        message=str(raw.message),  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
    )


def call_ida(
    resfn: Any,  # noqa: ANN401
    t_eval: NDArray[np.float64],
    y0: NDArray[np.float64],
    yp0: NDArray[np.float64],
    **options: Any,  # noqa: ANN401
) -> SundialsResult:
    """Create an IDA solver and run ``.solve(t_eval, y0, yp0)``.

    Parameters
    ----------
    resfn:
        Residual function with signature ``f(t, y, yp, res) -> None``
        (writes into *res* in-place).
    t_eval:
        Output time points.
    y0:
        Initial state vector.
    yp0:
        Initial time-derivative vector.
    **options:
        Forwarded to ``IDA.__init__``.  Recognized keys: ``rtol``,
        ``atol``, ``max_num_steps``, ``algebraic_idx``,
        ``calc_initcond``, ``calc_init_dt``, ``linsolver``, ``sparsity``,
        ``jacfn``, ``jactimes``.

    Returns
    -------
    SundialsResult
        Typed result container.
    """
    from sksundae.ida import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        IDA,
    )

    solver = IDA(resfn, **options)
    raw = solver.solve(t_eval, y0, yp0)
    return SundialsResult(
        t=np.asarray(raw.t, dtype=np.float64),  # pyright: ignore[reportUnknownArgumentType,reportUnknownMemberType]
        y=np.asarray(raw.y, dtype=np.float64),  # pyright: ignore[reportUnknownArgumentType,reportUnknownMemberType]
        success=bool(raw.success),  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        message=str(raw.message),  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        yp=np.asarray(raw.yp, dtype=np.float64),  # pyright: ignore[reportUnknownArgumentType,reportUnknownMemberType]
    )


def call_cvode_stepwise(
    rhsfn: Any,  # noqa: ANN401
    t_eval: NDArray[np.float64],
    y0: NDArray[np.float64],
    progress: SimulationProgress,
    snapshot_callback: Callable[..., None] | None = None,
    **options: Any,  # noqa: ANN401
) -> SundialsResult:
    """CVODE step-by-step integration with progress reporting.

    Uses ``init_step()`` + ``step(t)`` instead of ``solve()`` so that
    progress can be reported between solver steps with zero overhead on
    the solver's internal computation.
    """
    from sksundae.cvode import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        CVODE,
    )

    solver = CVODE(rhsfn, **options)
    progress.set_phase("CVODE init")
    solver.init_step(float(t_eval[0]), y0)  # pyright: ignore[reportUnknownMemberType]
    progress.set_phase("CVODE")

    t_list: list[float] = [float(t_eval[0])]
    y_list: list[NDArray[np.float64]] = [y0.copy()]

    if snapshot_callback is not None:
        snapshot_callback(t_list[0], y_list[0])

    success = True
    message = "CVODE stepwise integration completed"

    for t_next in t_eval[1:]:
        raw = solver.step(float(t_next))  # pyright: ignore[reportUnknownMemberType]
        t_val = float(raw.t)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        y_val = np.asarray(raw.y, dtype=np.float64)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        t_list.append(t_val)
        y_list.append(y_val)
        progress.update(t_val)
        if snapshot_callback is not None:
            snapshot_callback(t_val, y_val)
        if not raw.success:  # pyright: ignore[reportUnknownMemberType]
            success = False
            message = str(raw.message)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
            break

    progress.finish()
    return SundialsResult(
        t=np.array(t_list, dtype=np.float64),
        y=np.array(y_list, dtype=np.float64),
        success=success,
        message=message,
    )


def call_ida_stepwise(  # noqa: PLR0913, PLR0917
    resfn: Any,  # noqa: ANN401
    t_eval: NDArray[np.float64],
    y0: NDArray[np.float64],
    yp0: NDArray[np.float64],
    progress: SimulationProgress,
    snapshot_callback: Callable[..., None] | None = None,
    **options: Any,  # noqa: ANN401
) -> SundialsResult:
    """IDA step-by-step integration with progress reporting.

    Uses ``init_step()`` + ``step(t)`` instead of ``solve()`` so that
    progress can be reported between solver steps with zero overhead on
    the solver's internal computation.
    """
    from sksundae.ida import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        IDA,
    )

    solver = IDA(resfn, **options)
    progress.set_phase("IDA init")
    solver.init_step(float(t_eval[0]), y0, yp0)  # pyright: ignore[reportUnknownMemberType]
    progress.set_phase("IDA")

    t_list: list[float] = [float(t_eval[0])]
    y_list: list[NDArray[np.float64]] = [y0.copy()]
    yp_list: list[NDArray[np.float64]] = [yp0.copy()]

    if snapshot_callback is not None:
        snapshot_callback(t_list[0], y_list[0], yp_list[0])

    success = True
    message = "IDA stepwise integration completed"

    for t_next in t_eval[1:]:
        raw = solver.step(float(t_next))  # pyright: ignore[reportUnknownMemberType]
        t_val = float(raw.t)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        y_val = np.asarray(raw.y, dtype=np.float64)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        yp_val = np.asarray(raw.yp, dtype=np.float64)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
        t_list.append(t_val)
        y_list.append(y_val)
        yp_list.append(yp_val)
        progress.update(t_val)
        if snapshot_callback is not None:
            snapshot_callback(t_val, y_val, yp_val)
        if not raw.success:  # pyright: ignore[reportUnknownMemberType]
            success = False
            message = str(raw.message)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
            break

    progress.finish()
    return SundialsResult(
        t=np.array(t_list, dtype=np.float64),
        y=np.array(y_list, dtype=np.float64),
        success=success,
        message=message,
        yp=np.array(yp_list, dtype=np.float64),
    )
