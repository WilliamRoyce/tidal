"""Typed wrappers and type aliases for scipy integration and sparse solvers.

All scipy imports for ``scipy.integrate`` and ``scipy.sparse`` are isolated
here so that ``# pyright: ignore`` suppressions are confined to this module.
The rest of the codebase imports :class:`IVPResult`, :func:`call_solve_ivp`,
:func:`sparse_solve`, ``SparseMatrix``, and the sparse matrix constructors
from this module instead of importing scipy directly.

scipy does not ship complete type stubs for these sub-modules as of 1.16.x,
hence the inline suppressions below.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, TypeAlias

import numpy as np
from scipy import sparse as _sparse  # pyright: ignore[reportMissingTypeStubs]
from scipy.sparse.linalg import (  # pyright: ignore[reportMissingTypeStubs]
    spsolve as _spsolve,  # pyright: ignore[reportUnknownVariableType]
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

# ---------------------------------------------------------------------------
# Sparse matrix constructors re-exported from this module.
# Callers use ``from tidal.solver._scipy_types import lil_matrix, csc_matrix``
# instead of importing scipy directly.  The types are ``Any`` because scipy
# sparse matrices are untyped, but the names document intent clearly.
# ---------------------------------------------------------------------------

lil_matrix: Any = _sparse.lil_matrix  # pyright: ignore[reportUnknownMemberType]
csc_matrix: Any = _sparse.csc_matrix  # pyright: ignore[reportUnknownMemberType]
coo_matrix: Any = _sparse.coo_matrix  # pyright: ignore[reportUnknownMemberType]
speye: Any = _sparse.eye  # pyright: ignore[reportUnknownMemberType,reportUnknownVariableType]
diags: Any = _sparse.diags  # pyright: ignore[reportUnknownMemberType,reportUnknownVariableType]

# Type alias for scipy sparse matrices.  Opaque to pyright (= Any) but
# communicates intent in function signatures.
SparseMatrix: TypeAlias = Any


# ---------------------------------------------------------------------------
# IVP solver wrapper
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IVPResult:
    """Typed result from ``scipy.integrate.solve_ivp``.

    ``y`` has shape ``(n_times, n_states)`` — transposed from scipy's native
    ``(n_states, n_times)`` convention to match the TIDAL / SUNDIALS
    convention used by :class:`SundialsResult`.
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]  # shape (n_times, n_states)
    success: bool
    message: str


def call_solve_ivp(  # noqa: PLR0913
    fun: Callable[[float, NDArray[np.float64]], NDArray[np.float64]],
    t_span: tuple[float, float],
    y0: NDArray[np.float64],
    *,
    method: str,
    t_eval: NDArray[np.float64],
    rtol: float,
    atol: float,
    max_step: float,
    jac_sparsity: SparseMatrix | None = None,
) -> IVPResult:
    """Typed wrapper around ``scipy.integrate.solve_ivp``.

    Handles the ``(n_states, n_times)`` → ``(n_times, n_states)``
    transposition internally, returning an :class:`IVPResult` with the
    TIDAL array convention.

    Parameters
    ----------
    fun:
        RHS function ``f(t, y) -> dy/dt``.
    t_span:
        ``(t_start, t_end)`` integration interval.
    y0:
        Initial state vector.
    method:
        Integrator name, e.g. ``"DOP853"``, ``"RK45"``, ``"BDF"``.
    t_eval:
        Output time points.
    rtol:
        Relative tolerance.
    atol:
        Absolute tolerance.
    max_step:
        Maximum allowed step size (``numpy.inf`` for no limit).
    jac_sparsity:
        Optional Jacobian sparsity pattern (for implicit methods).

    Returns
    -------
    IVPResult
        Typed result container with ``y`` in ``(n_times, n_states)`` order.
    """
    from scipy.integrate import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        solve_ivp,  # pyright: ignore[reportUnknownVariableType]
    )

    kwargs: dict[str, Any] = {
        "t_eval": t_eval,
        "rtol": rtol,
        "atol": atol,
        "max_step": max_step,
        "dense_output": False,
    }
    if jac_sparsity is not None:
        kwargs["jac_sparsity"] = jac_sparsity

    raw: Any = solve_ivp(fun, t_span, y0, method=method, **kwargs)  # pyright: ignore[reportUnknownVariableType]

    raw_y: Any = raw.y  # pyright: ignore[reportUnknownMemberType,reportUnknownVariableType]
    y_out: NDArray[np.float64] = (
        np.asarray(raw_y, dtype=np.float64).T  # pyright: ignore[reportUnknownArgumentType]
        if raw_y is not None
        else np.empty((0, len(y0)), dtype=np.float64)
    )
    raw_t: Any = raw.t  # pyright: ignore[reportUnknownMemberType,reportUnknownVariableType]
    t_out: NDArray[np.float64] = (
        np.asarray(raw_t, dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]
        if raw_t is not None
        else np.empty(0, dtype=np.float64)
    )
    raw_msg: Any = raw.message  # pyright: ignore[reportUnknownMemberType,reportUnknownVariableType]
    msg: str = (
        str(raw_msg)  # pyright: ignore[reportUnknownArgumentType]
        if raw_msg
        else (
            "scipy integration completed"
            if bool(raw.success)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]
            else "scipy integration failed"
        )
    )
    return IVPResult(t=t_out, y=y_out, success=bool(raw.success), message=msg)  # pyright: ignore[reportUnknownMemberType,reportUnknownArgumentType]


# ---------------------------------------------------------------------------
# Sparse direct solver wrapper
# ---------------------------------------------------------------------------


def sparse_solve(
    matrix: SparseMatrix,
    rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Solve ``matrix @ x = rhs`` via sparse direct factorization.

    Wraps ``scipy.sparse.linalg.spsolve`` with a typed signature.

    Parameters
    ----------
    matrix:
        Square CSC sparse matrix.
    rhs:
        Right-hand side vector (1-D).

    Returns
    -------
    NDArray[np.float64]
        Solution vector.
    """
    result: Any = _spsolve(matrix, rhs)  # pyright: ignore[reportUnknownVariableType]
    return np.asarray(result, dtype=np.float64)  # pyright: ignore[reportUnknownArgumentType]
