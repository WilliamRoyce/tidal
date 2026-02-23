"""Typed wrappers for sksundae CVODE and IDA solvers.

All sksundae imports are isolated here so that ``# pyright: ignore``
suppressions for ``reportMissingTypeStubs`` are confined to this module.
The rest of the codebase imports ``call_cvode`` / ``call_ida`` and receives
fully typed :class:`SundialsResult` objects.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


@dataclass(frozen=True)
class SundialsResult:
    """Typed result from ``CVODE.solve()`` or ``IDA.solve()``.

    Both CVODE and IDA return result objects with identical attribute names
    and the same array shape convention (``n_times`` x ``n_states``).
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]  # shape (n_times, n_states) — TIDAL convention
    success: bool
    message: str


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
        ``calc_initcond``, ``calc_init_dt``, ``linsolver``, ``sparsity``.

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
    )
