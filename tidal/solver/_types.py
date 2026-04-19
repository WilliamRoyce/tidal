"""Shared type definitions and constants for TIDAL solvers.

``SolverResult`` provides a typed contract for all solver return values.
Threshold constants for linear solver selection are shared by CVODE and IDA.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, TypedDict

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


class SolverResult(TypedDict):
    """Typed result dict returned by all TIDAL solver functions.

    All four solvers (CVODE, IDA, scipy, leapfrog) return dicts with
    these exact keys.  Using a TypedDict rather than ``dict[str, Any]``
    gives pyright full key/value type checking at every call site.
    """

    t: NDArray[np.float64]
    y: NDArray[np.float64]
    success: bool
    message: str


class PerturbativePass1Result(SolverResult):
    """Typed result dict for :func:`solve_modal_pass1`.

    Extends :class:`SolverResult` with the two additional keys required
    by the v6 perturbative driver's Gap C constraint recovery and #272
    silent-drop diagnostics. Using a dedicated TypedDict (rather than
    a `# type: ignore` on key injection into SolverResult) lets
    downstream callers typo-check field access and removes the silent
    `result.get("y_hat_dyn")` defensive-None path.
    """

    y_hat_dyn: NDArray[np.complex128]
    correction_drops: list[dict[str, object]]


# System size thresholds for linear solver selection.
# Dense Jacobian at N=2000: N^2 * 8 = 32 MB (safe for all environments).
# Sparse (SuperLU_MT) scales well up to ~200K state variables.
# Beyond that, matrix-free GMRES avoids storing any Jacobian.
DENSE_THRESHOLD: int = 2_000
SPARSE_THRESHOLD: int = 200_000

# Sparsity pattern nnz limit for SuperLU_MT. Above this, SuperLU's fill-in
# during LU factorization can exhaust memory. Fall back to FD GMRES instead.
# Calibrated against massive_gravity (32x32 fails at nnz~172K, 16x16 passes
# at nnz~43K). Value of 100_000 gives a conservative safety margin.
# Reference: Li & Demmel, "SuperLU_MT", ACM TOMS 29(2), 2003.
SUPERLU_NNZ_LIMIT: int = 100_000
