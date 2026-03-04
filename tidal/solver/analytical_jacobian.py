"""Analytical Jacobian for time-independent linear systems.

For time-independent linear systems (most TIDAL examples), the IDA
Jacobian ``J = dF/dy + cj * dF/dyp`` consists of two constant sparse
matrices.  Instead of finite-difference approximation (which requires
O(n_colors) residual evaluations per Newton step), we precompute ``dF/dy``
and ``dF/dyp`` once and supply them analytically.

Two delivery modes depending on system size:

- **Dense tier** (N <= DENSE_THRESHOLD):  ``jacfn`` fills a 2D numpy
  array with ``dF_dy + cj * dF_dyp``.
- **GMRES tier** (N > SPARSE_THRESHOLD):  ``jactimes`` provides an
  analytical Jacobian-vector product ``Jv = dF_dy @ v + cj * dF_dyp @ v``
  using sparse matrix-vector products, eliminating finite-difference
  residual evaluations per GMRES iteration.

For the sparse tier (DENSE_THRESHOLD < N <= SPARSE_THRESHOLD), the system
falls through to the normal FD-based ``configure_linear_solver`` path.
Colored finite-differences with SuperLU_MT direct factorisation outperform
unpreconditioned GMRES at these sizes.  A sparse analytical ``jacfn`` (1D
CSC data) is implemented here (``_create_sparse_jacfn``), but disabled
because sksundae v1.1.1 has a bug where ``_setup_memory()``
unconditionally overwrites ``aux.jacfn`` when a sparsity pattern is
provided (``_cy_ida.pyx:720``), preventing user-supplied ``jacfn`` from
being called.

Position-dependent (but time-independent) coefficients are supported:
the spatial grid is fixed, so the Jacobian is still constant.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._scipy_types import SparseMatrix, diags, lil_matrix, speye

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.state import SlotInfo, StateLayout
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        OperatorTerm,
    )

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Time-independence detection
# ---------------------------------------------------------------------------


def _is_system_time_independent(
    spec: EquationSystem,
    grid: GridInfo,
    parameters: dict[str, float],
) -> bool:
    """Check if the system Jacobian is constant (time-independent).

    Returns True when ``CoefficientEvaluator.all_time_independent()``
    holds.  This includes both fully constant coefficients AND position-
    dependent (but time-independent) coefficients, since the spatial grid
    is fixed and so the Jacobian doesn't vary over time.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    coeff_eval = CoefficientEvaluator(spec, grid, parameters)
    return coeff_eval.all_time_independent()


# ---------------------------------------------------------------------------
# Operator matrix construction
# ---------------------------------------------------------------------------


def build_operator_matrix(
    operator: str,
    grid: GridInfo,
    bc: BCSpec | None,
) -> SparseMatrix:
    """Build the N x N sparse matrix for a single spatial operator by probing.

    Probes ``apply_operator(op, e_j)`` with unit vectors to construct the
    matrix column-by-column.  Used by both the analytical Jacobian builder
    and the constraint IC solver.
    """
    from tidal.solver.operators import apply_operator  # noqa: PLC0415

    n = grid.num_points
    mat = lil_matrix((n, n))
    e_j: NDArray[np.float64] = np.zeros(grid.shape)

    for j in range(n):
        e_j.flat[:] = 0.0
        e_j.flat[j] = 1.0
        col = apply_operator(operator, e_j, grid, bc)
        col_flat = col.ravel()
        nz = np.nonzero(col_flat)[0]
        for row in nz:
            mat[row, j] = col_flat[row]

    return mat.tocsc()


class _OperatorCache:
    """Caches operator matrices for a fixed (grid, bc) pair."""

    def __init__(self, grid: GridInfo, bc: BCSpec | None) -> None:
        self._grid = grid
        self._bc = bc
        self._cache: dict[str, SparseMatrix] = {}
        self._identity: SparseMatrix | None = None

    def get_identity(self) -> SparseMatrix:
        """Return the identity matrix (cached)."""
        if self._identity is None:
            self._identity = speye(self._grid.num_points, format="csc")
        return self._identity

    def get(self, operator: str) -> SparseMatrix:
        """Return operator matrix, building on first access."""
        if operator == "identity":
            return self.get_identity()
        if operator not in self._cache:
            self._cache[operator] = build_operator_matrix(
                operator, self._grid, self._bc
            )
        return self._cache[operator]


# ---------------------------------------------------------------------------
# Term target resolution
# ---------------------------------------------------------------------------


def _resolve_term_target(  # noqa: PLR0911
    term: OperatorTerm,
    layout: StateLayout,
    constraint_fields: set[str],
    op_cache: _OperatorCache,
) -> tuple[int, SparseMatrix, bool] | None:
    """Determine where an RHS term couples in the Jacobian.

    Returns ``(col_slot_idx, operator_matrix, is_dyp)`` or ``None`` if the
    term's field reference cannot be resolved to a state slot.

    ``is_dyp=True`` means the coupling enters ``dF/dyp`` (the term depends
    on ``yp`` via constraint velocity injection); ``False`` means ``dF/dy``.
    """
    I_mat = op_cache.get_identity()  # noqa: N806

    # Mechanism A: first_derivative_t(X) resolves to velocity of X
    if term.operator == "first_derivative_t":
        if term.field in constraint_fields:
            # Constraint velocity: v_X = yp[X_field_slot]
            return (layout.field_slot_map[term.field], I_mat, True)
        # Dynamical velocity: v_X = y[velocity_slot]
        vel_slot = layout.velocity_slot_map.get(term.field)
        if vel_slot is not None:
            return (vel_slot, I_mat, False)
        return None

    # Mechanism B: explicit v_X field reference (e.g. gradient_x(v_A_0))
    if term.field.startswith("v_"):
        base_field = term.field[2:]
        if base_field in constraint_fields:
            # Constraint velocity: v_X = yp[X_field_slot]
            op_mat = op_cache.get(term.operator)
            return (layout.field_slot_map[base_field], op_mat, True)
        # Dynamical velocity: v_X = y[velocity_slot]
        vel_slot = layout.velocity_slot_map.get(base_field)
        if vel_slot is not None:
            op_mat = op_cache.get(term.operator)
            return (vel_slot, op_mat, False)
        # Fall through: field might literally be named "v_something"

    # Normal field reference
    field_slot = layout.field_slot_map.get(term.field)
    if field_slot is None:
        return None
    op_mat = op_cache.get(term.operator)
    return (field_slot, op_mat, False)


# ---------------------------------------------------------------------------
# Jacobian matrix builder
# ---------------------------------------------------------------------------

# Reuse the same detection logic as ida.py for gauge/no-self-term cases.
LAPLACIAN_OPS = frozenset({
    "laplacian", "laplacian_x", "laplacian_y", "laplacian_z",
})


def _is_pure_laplacian(eq: ComponentEquation) -> bool:
    """Check if equation has Laplacian self-terms but no identity/mass."""
    has_lap = False
    for term in eq.rhs_terms:
        if term.field != eq.field_name:
            continue
        if term.operator in LAPLACIAN_OPS:
            has_lap = True
        elif term.operator == "identity":
            return False
    return has_lap


def _detect_no_self_term_fields(spec: EquationSystem) -> set[str]:
    """Fields whose constraint equations have no self-referencing terms."""
    result: set[str] = set()
    for eq in spec.equations:
        if eq.time_derivative_order != 0:
            continue
        if eq.constraint_solver.enabled:
            continue
        has_self = any(t.field == eq.field_name for t in eq.rhs_terms)
        if not has_self:
            result.add(eq.field_name)
    return result


def _detect_gauge_fix_fields(
    spec: EquationSystem,
    grid: GridInfo,
    bc: BCSpec | None,
) -> set[str]:
    """Fields needing gauge regularization (pure Laplacian + periodic BCs)."""
    from tidal.solver.operators import is_periodic_bc  # noqa: PLC0415

    # Check all-periodic
    if bc is not None:
        bcs = (bc,) * grid.ndim if isinstance(bc, str) else tuple(bc)
        all_periodic = all(is_periodic_bc(b) for b in bcs)
    else:
        all_periodic = all(grid.periodic)
    if not all_periodic:
        return set()

    return {
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order == 0
        and eq.constraint_solver.enabled
        and _is_pure_laplacian(eq)
    }


def build_jacobian_matrices(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    parameters: dict[str, float],
) -> tuple[SparseMatrix, SparseMatrix]:
    """Build the analytical dF/dy and dF/dyp sparse matrices.

    Mirrors the residual handler structure in ``ida.py`` exactly:
    constraint (3 sub-cases), velocity, dynamical_field, first_order.

    Parameters
    ----------
    spec : EquationSystem
        Parsed equation specification.
    layout : StateLayout
        State vector layout.
    grid : GridInfo
        Spatial grid.
    bc : BCSpec or None
        Boundary conditions.
    parameters : dict[str, float]
        Resolved parameter values.

    Returns
    -------
    (dF_dy, dF_dyp) : tuple[SparseMatrix, SparseMatrix]
        Sparse Jacobian component matrices.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    n = grid.num_points
    n_total = layout.total_size
    eq_map: dict[str, int] = spec.equation_map

    coeff_eval = CoefficientEvaluator(spec, grid, parameters)
    op_cache = _OperatorCache(grid, bc)
    I_mat = op_cache.get_identity()  # noqa: N806

    constraint_fields = {
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order == 0
    }
    no_self_term_fields = _detect_no_self_term_fields(spec)
    gauge_fix_fields = _detect_gauge_fix_fields(spec, grid, bc)

    # Build in LIL format for efficient incremental construction
    dF_dy = lil_matrix((n_total, n_total))  # noqa: N806
    dF_dyp = lil_matrix((n_total, n_total))  # noqa: N806

    for slot_idx, slot in enumerate(layout.slots):
        if slot.time_order == 0:
            _build_constraint_block(
                slot_idx, slot, layout, spec, eq_map,
                coeff_eval, op_cache, constraint_fields,
                no_self_term_fields, gauge_fix_fields,
                n, dF_dy, dF_dyp,
            )

        elif slot.kind == "velocity":
            # res = yp - RHS → dF/dyp[slot, slot] = I
            _set_diagonal_identity(dF_dyp, slot_idx, n, I_mat)
            eq_idx = eq_map.get(slot.field_name)
            if eq_idx is not None:
                _add_rhs_terms(
                    slot_idx, eq_idx, spec, layout, coeff_eval,
                    op_cache, constraint_fields, n, dF_dy, dF_dyp,
                    negate=True,
                )

        elif slot.time_order >= 2 and slot.kind == "field":  # noqa: PLR2004
            # res = yp - v → dF/dyp[slot, slot] = I, dF/dy[slot, vel] = -I
            _set_diagonal_identity(dF_dyp, slot_idx, n, I_mat)
            vel_slot = layout.velocity_slot_map.get(slot.field_name)
            if vel_slot is not None:
                _set_block(dF_dy, slot_idx, vel_slot, n, -I_mat)

        elif slot.time_order == 1:
            # res = yp - RHS → dF/dyp[slot, slot] = I
            _set_diagonal_identity(dF_dyp, slot_idx, n, I_mat)
            eq_idx = eq_map.get(slot.field_name)
            if eq_idx is not None:
                _add_rhs_terms(
                    slot_idx, eq_idx, spec, layout, coeff_eval,
                    op_cache, constraint_fields, n, dF_dy, dF_dyp,
                    negate=True,
                )

    return dF_dy.tocsc(), dF_dyp.tocsc()


# ---------------------------------------------------------------------------
# Block-assembly helpers
# ---------------------------------------------------------------------------


def _set_diagonal_identity(
    mat: SparseMatrix,
    slot_idx: int,
    n: int,
    I_mat: SparseMatrix,  # noqa: N803
) -> None:
    """Set diagonal block mat[slot, slot] = I."""
    s = slice(slot_idx * n, (slot_idx + 1) * n)
    mat[s, s] = I_mat


def _set_block(
    mat: SparseMatrix,
    row_slot: int,
    col_slot: int,
    n: int,
    block: SparseMatrix,
) -> None:
    """Set block mat[row_slot, col_slot] = block."""
    rs = slice(row_slot * n, (row_slot + 1) * n)
    cs = slice(col_slot * n, (col_slot + 1) * n)
    mat[rs, cs] = block


def _add_block(
    mat: SparseMatrix,
    row_slot: int,
    col_slot: int,
    n: int,
    block: SparseMatrix,
) -> None:
    """Add block: mat[row_slot, col_slot] += block."""
    rs = slice(row_slot * n, (row_slot + 1) * n)
    cs = slice(col_slot * n, (col_slot + 1) * n)
    mat[rs, cs] += block


def _add_rhs_terms(  # noqa: PLR0913, PLR0917
    row_slot: int,
    eq_idx: int,
    spec: EquationSystem,
    layout: StateLayout,
    coeff_eval: CoefficientEvaluator,
    op_cache: _OperatorCache,
    constraint_fields: set[str],
    n: int,
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
    *,
    negate: bool,
) -> None:
    """Add all RHS term contributions for an equation.

    When ``negate=True``, contributions are negated (for ``res = yp - RHS``).
    When ``negate=False``, contributions are added as-is (for ``res = RHS``).
    """
    eq = spec.equations[eq_idx]
    sign = -1.0 if negate else 1.0

    for term_idx, term in enumerate(eq.rhs_terms):
        resolved = _resolve_term_target(
            term, layout, constraint_fields, op_cache,
        )
        if resolved is None:
            continue

        col_slot, op_mat, is_dyp = resolved
        coeff = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
        if isinstance(coeff, np.ndarray):
            # Position-dependent coefficient: diag(c(x)) @ Op
            scaled = sign * diags(coeff.ravel()) @ op_mat
        else:
            scaled = sign * float(coeff) * op_mat

        target_mat = dF_dyp if is_dyp else dF_dy
        _add_block(target_mat, row_slot, col_slot, n, scaled)


def _build_constraint_block(  # noqa: PLR0913, PLR0917
    slot_idx: int,
    slot: SlotInfo,
    layout: StateLayout,
    spec: EquationSystem,
    eq_map: dict[str, int],
    coeff_eval: CoefficientEvaluator,
    op_cache: _OperatorCache,
    constraint_fields: set[str],
    no_self_term_fields: set[str],
    gauge_fix_fields: set[str],
    n: int,
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
) -> None:
    """Build Jacobian blocks for a constraint slot (time_order=0).

    Three sub-cases matching ida.py handle_constraint():
    1. No-self-term: res = y[field] → dF/dy = I
    2. Gauge-fix: normal RHS, then pin row 0 in both matrices
    3. Normal: res = RHS
    """
    eq_idx = eq_map.get(slot.field_name)
    if eq_idx is None:
        return

    # Case 1: no self-terms — field frozen at zero, res = y[field]
    if slot.field_name in no_self_term_fields:
        field_slot = layout.field_slot_map[slot.field_name]
        I_mat = op_cache.get_identity()  # noqa: N806
        _set_block(dF_dy, slot_idx, field_slot, n, I_mat)
        return

    # Case 3 (and setup for case 2): normal RHS
    # res = RHS → dF/dy += coeff * op_mat (no negation)
    _add_rhs_terms(
        slot_idx, eq_idx, spec, layout, coeff_eval,
        op_cache, constraint_fields, n, dF_dy, dF_dyp,
        negate=False,
    )

    # Case 2: gauge-fix — overwrite row 0 with y[field_slot*n] = 0
    if slot.field_name in gauge_fix_fields:
        field_slot = layout.field_slot_map[slot.field_name]
        global_row = slot_idx * n
        # Zero out entire row in both matrices
        dF_dy[global_row, :] = 0
        dF_dyp[global_row, :] = 0
        # Set pinning entry: dF/dy[row, field_col] = 1
        dF_dy[global_row, field_slot * n] = 1.0


# ---------------------------------------------------------------------------
# Jacobian delivery: dense jacfn
# ---------------------------------------------------------------------------


def _create_jacfn(
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
) -> Callable[..., None]:
    """Create a dense ``jacfn`` callback for IDA.

    Signature: ``jacfn(t, y, yp, res, cj, JJ)`` — fills the pre-allocated
    2D numpy array ``JJ`` with ``dF_dy + cj * dF_dyp``.
    """
    jac_y = dF_dy.toarray()
    jac_yp = dF_dyp.toarray()
    jac_yp_scaled = np.empty_like(jac_yp)

    def jacfn(  # noqa: PLR0913, PLR0917
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        res: NDArray[np.float64],  # noqa: ARG001
        cj: float,
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        np.multiply(jac_yp, cj, out=jac_yp_scaled)
        np.add(jac_y, jac_yp_scaled, out=JJ)

    return jacfn


# ---------------------------------------------------------------------------
# CVODE-specific delivery (ODE Jacobian df/dy = -dF_dy)
# ---------------------------------------------------------------------------


def _create_cvode_jacfn(
    dF_dy: SparseMatrix,  # noqa: N803
) -> Callable[..., None]:
    """Create a dense ``jacfn`` callback for CVODE.

    Signature: ``jacfn(t, y, yp, JJ)`` — fills ``JJ`` with the ODE
    Jacobian ``df/dy = -dF_dy`` (negated because ``dF_dy`` stores the
    IDA residual derivatives where ``F = yp - RHS``).
    """
    ode_jac = (-dF_dy).toarray()

    def jacfn(
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        JJ[:] = ode_jac

    return jacfn


# ---------------------------------------------------------------------------
# Sparse tier delivery: 1D jacfn (CSC data order)
# ---------------------------------------------------------------------------


def _prepare_sparse_data(  # pyright: ignore[reportUnusedFunction]  # reserved for sparse tier
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
) -> tuple[NDArray[np.float64], NDArray[np.float64], SparseMatrix]:
    """Align two CSC matrices to their union sparsity pattern.

    Returns ``(dy_data, dyp_data, pattern)`` where ``pattern`` is a CSC
    matrix encoding the union sparsity, and ``dy_data``/``dyp_data`` are
    1D arrays aligned to ``pattern``'s CSC data ordering.

    The sparse ``jacfn`` then computes ``JJ[:] = dy_data + cj * dyp_data``
    — a single O(nnz) vectorized operation with zero allocation.
    """
    # Build union sparsity via structural OR.
    # Replace data with ones to get structural patterns, then add.
    ones_dy = dF_dy.copy()
    ones_dy.data[:] = 1.0
    ones_dyp = dF_dyp.copy()
    ones_dyp.data[:] = 1.0
    pattern = (ones_dy + ones_dyp).tocsc()
    pattern.sort_indices()
    pattern.data[:] = 1.0  # Normalize to binary

    nnz = pattern.nnz
    dy_data = np.zeros(nnz, dtype=np.float64)
    dyp_data = np.zeros(nnz, dtype=np.float64)

    dy_csc = dF_dy.tocsc()
    dy_csc.sort_indices()
    dyp_csc = dF_dyp.tocsc()
    dyp_csc.sort_indices()

    n_cols = pattern.shape[1]
    for col in range(n_cols):
        p_start = pattern.indptr[col]
        p_rows = pattern.indices[p_start:pattern.indptr[col + 1]]

        # Project dF_dy values onto union pattern for this column
        d_start = dy_csc.indptr[col]
        d_end = dy_csc.indptr[col + 1]
        if d_start < d_end:
            d_rows = dy_csc.indices[d_start:d_end]
            idxs = np.searchsorted(p_rows, d_rows)  # pyright: ignore[reportUnknownVariableType]
            dy_data[p_start + idxs] = dy_csc.data[d_start:d_end]

        # Project dF_dyp values onto union pattern for this column
        d_start = dyp_csc.indptr[col]
        d_end = dyp_csc.indptr[col + 1]
        if d_start < d_end:
            d_rows = dyp_csc.indices[d_start:d_end]
            idxs = np.searchsorted(p_rows, d_rows)  # pyright: ignore[reportUnknownVariableType]
            dyp_data[p_start + idxs] = dyp_csc.data[d_start:d_end]

    return dy_data, dyp_data, pattern


def _create_sparse_jacfn(  # pyright: ignore[reportUnusedFunction]  # reserved for sparse tier
    dy_data: NDArray[np.float64],
    dyp_data: NDArray[np.float64],
) -> Callable[..., None]:
    """Create a sparse ``jacfn`` callback for IDA.

    Signature: ``jacfn(t, y, yp, res, cj, JJ)`` where ``JJ`` is a 1D
    numpy array of size ``nnz`` in CSC data order.  sksundae's
    ``np2smat_sparse1D`` transfers this directly to SUNDIALS.
    """
    def jacfn(  # noqa: PLR0913, PLR0917
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        res: NDArray[np.float64],  # noqa: ARG001
        cj: float,
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        np.multiply(dyp_data, cj, out=JJ)
        np.add(dy_data, JJ, out=JJ)

    return jacfn


def _create_cvode_sparse_jacfn(  # pyright: ignore[reportUnusedFunction]  # reserved for sparse tier
    neg_dy_data: NDArray[np.float64],
) -> Callable[..., None]:
    """Create a sparse ``jacfn`` callback for CVODE.

    Signature: ``jacfn(t, y, yp, JJ)`` where ``JJ`` is 1D (CSC data).
    ODE Jacobian is ``-dF_dy`` (negation of DAE residual derivative).
    """
    def jacfn(
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        JJ[:] = neg_dy_data

    return jacfn


# ---------------------------------------------------------------------------
# GMRES tier delivery: Jacobian-vector product (jactimes)
# ---------------------------------------------------------------------------


def _create_ida_jactimes(
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
) -> Callable[..., None]:
    """Create a ``jactimes.solvefn`` callback for IDA GMRES.

    Signature: ``solvefn(t, y, yp, res, v, Jv, cj)`` — fills ``Jv``
    with ``(dF_dy + cj * dF_dyp) @ v`` using two sparse mat-vec products.

    Avoids finite-difference residual evaluations per GMRES iteration.
    """
    dy_csc = dF_dy.tocsc()
    dyp_csc = dF_dyp.tocsc()

    def solvefn(  # noqa: PLR0913, PLR0917
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        res: NDArray[np.float64],  # noqa: ARG001
        v: NDArray[np.float64],
        jv: NDArray[np.float64],
        cj: float,
    ) -> None:
        # Two sparse mat-vec products: Jv = dF_dy @ v + cj * dF_dyp @ v
        jv[:] = dy_csc @ v
        jv += cj * (dyp_csc @ v)

    return solvefn


def _create_cvode_jactimes(
    dF_dy: SparseMatrix,  # noqa: N803
) -> Callable[..., None]:
    """Create a ``jactimes.solvefn`` callback for CVODE GMRES.

    Signature: ``solvefn(t, y, yp, v, Jv)`` — fills ``Jv`` with
    ``(-dF_dy) @ v`` (the ODE Jacobian-vector product).
    """
    neg_dy_csc = (-dF_dy).tocsc()

    def solvefn(
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        v: NDArray[np.float64],
        jv: NDArray[np.float64],
    ) -> None:
        jv[:] = neg_dy_csc @ v

    return solvefn


# ---------------------------------------------------------------------------
# Public API: integration into configure_linear_solver
# ---------------------------------------------------------------------------


def try_analytical_jacobian(  # noqa: PLR0913, PLR0917
    options: dict[str, Any],
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    parameters: dict[str, float],
    *,
    solver: str = "ida",
) -> bool:
    """Try to configure analytical Jacobian for time-independent systems.

    Mutates *options* in-place if successful.  Returns True on success,
    False if the system has time-dependent coefficients (caller should
    fall back to the finite-difference tier system).

    Two delivery modes by system size:

    - **Dense** (N <= DENSE_THRESHOLD): 2D ``jacfn`` fills dense array.
    - **GMRES** (N > SPARSE_THRESHOLD): ``jactimes`` provides analytical
      Jacobian-vector product for iterative GMRES.

    Systems in the sparse tier (DENSE_THRESHOLD < N <= SPARSE_THRESHOLD)
    return False, falling through to FD sparse (colored FD + SuperLU_MT)
    which outperforms unpreconditioned GMRES at these sizes.

    Parameters
    ----------
    solver : str
        ``"ida"`` or ``"cvode"``.  Controls the callback signature.
    """
    from tidal.solver._types import (  # noqa: PLC0415
        DENSE_THRESHOLD,
        SPARSE_THRESHOLD,
    )

    n_state = layout.total_size

    if not _is_system_time_independent(spec, grid, parameters):
        return False

    # Sparse tier: fall through to FD colored-Jacobian + SuperLU_MT.
    # This outperforms unpreconditioned GMRES with analytical jactimes
    # at moderate sizes.  A sparse analytical jacfn is implemented
    # (_create_sparse_jacfn) but disabled because sksundae v1.1.1 has
    # a bug where _setup_memory() overwrites aux.jacfn when sparsity
    # is provided (_cy_ida.pyx:720).
    if DENSE_THRESHOLD < n_state <= SPARSE_THRESHOLD:
        return False

    if n_state <= DENSE_THRESHOLD:
        # Dense tier: 2D jacfn
        jac_y, jac_yp = build_jacobian_matrices(
            spec, layout, grid, bc, parameters,
        )
        options["linsolver"] = "dense"
        if solver == "cvode":
            options["jacfn"] = _create_cvode_jacfn(jac_y)
        else:
            options["jacfn"] = _create_jacfn(jac_y, jac_yp)
        logger.info(
            "Analytical Jacobian (dense %s jacfn) for %d-state system",
            solver, n_state,
        )

    else:
        # GMRES tier (N > SPARSE_THRESHOLD): analytical Jacobian-vector
        # product.  Eliminates O(n_colors) residual evaluations per GMRES
        # iteration compared to the FD GMRES path.
        jac_y, jac_yp = build_jacobian_matrices(
            spec, layout, grid, bc, parameters,
        )
        if solver == "cvode":
            from sksundae.cvode import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
                CVODEJacTimes,  # pyright: ignore[reportUnknownVariableType]
            )

            solvefn = _create_cvode_jactimes(jac_y)
            options["linsolver"] = "gmres"
            options["jactimes"] = CVODEJacTimes(setupfn=None, solvefn=solvefn)  # pyright: ignore[reportUnknownArgumentType]
        else:
            from sksundae.ida import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
                IDAJacTimes,  # pyright: ignore[reportUnknownVariableType]
            )

            solvefn = _create_ida_jactimes(jac_y, jac_yp)
            options["linsolver"] = "gmres"
            options["jactimes"] = IDAJacTimes(setupfn=None, solvefn=solvefn)  # pyright: ignore[reportUnknownArgumentType]
        logger.info(
            "Analytical Jacobian (GMRES %s jactimes) for %d-state system",
            solver, n_state,
        )

    return True
