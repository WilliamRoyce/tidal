"""Analytical Jacobian for time-independent linear systems.

For time-independent linear systems (most TIDAL examples), the IDA
Jacobian ``J = dF/dy + cj * dF/dyp`` consists of two constant sparse
matrices.  Instead of finite-difference approximation (which requires
O(n_colors) residual evaluations per Newton step), we precompute ``dF/dy``
and ``dF/dyp`` once and supply them analytically.

Three delivery modes depending on system size:

- **Dense tier** (N <= DENSE_THRESHOLD):  ``jacfn`` fills a 2D numpy
  array with ``dF_dy + cj * dF_dyp``.  ~5.3x speedup vs FD.
- **Sparse tier** (DENSE_THRESHOLD < N <= SPARSE_THRESHOLD):  ``jacfn``
  fills a 1D CSC data array with ``dF_dy.data + cj * dF_dyp.data``,
  paired with a sparsity pattern for SuperLU_MT direct factorisation.
  Requires sksundae >= 1.1.2.  IDA: ~2.5x, CVODE: ~1.2-1.4x.
- **GMRES tier** (N > SPARSE_THRESHOLD):  ``jactimes`` provides an
  analytical Jacobian-vector product ``Jv = dF_dy @ v + cj * dF_dyp @ v``
  using sparse matrix-vector products, eliminating finite-difference
  residual evaluations per GMRES iteration.

Performance optimizations:

- **COO accumulation**: ``build_jacobian_matrices()`` uses
  ``_COOAccumulator`` to append block triples and do a single CSC
  conversion, avoiding O(N²)-per-block LIL slice assignment.
- **Circulant operators**: For all-periodic BCs,
  ``build_operator_matrix()`` probes a single grid point and tiles the
  stencil via modular arithmetic — O(nnz) instead of O(N²) probing.

Position-dependent (but time-independent) coefficients are supported:
the spatial grid is fixed, so the Jacobian is still constant.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._scipy_types import (
    SparseMatrix,
    diags,
    lil_matrix,
    speye,
)

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


def _build_operator_matrix_circulant(
    operator: str,
    grid: GridInfo,
    bc: BCSpec | None,
) -> SparseMatrix:
    """Build operator matrix for all-periodic BCs using translation invariance.

    For periodic BCs on a uniform grid, the operator is translation-invariant:
    probing a single center point gives the full stencil, which is then tiled
    across all grid points.  This is O(nnz) instead of O(N²).
    """
    from tidal.solver.operators import apply_operator  # noqa: PLC0415

    n = grid.num_points
    shape = grid.shape

    # Probe center point to get stencil
    center_nd = tuple(s // 2 for s in shape)
    e0: NDArray[np.float64] = np.zeros(shape)
    e0[center_nd] = 1.0
    col0 = apply_operator(operator, e0, grid, bc)

    nz_nd = np.nonzero(col0)
    vals = col0[nz_nd]
    # Compute N-D offsets relative to center
    offsets_nd = [nz_nd[d] - center_nd[d] for d in range(len(shape))]

    # Build COO data: for each grid point, apply the stencil with wrapping
    all_coords = np.array(np.unravel_index(np.arange(n), shape))  # (ndim, n)
    rows_list: list[NDArray[np.intp]] = []
    cols_list: list[NDArray[np.intp]] = []
    data_list: list[NDArray[np.float64]] = []
    arange_n = np.arange(n)

    for k in range(len(vals)):
        # For column j (all grid points), the output row is at j + offset
        row_coords = tuple(
            (all_coords[d] + offsets_nd[d][k]) % shape[d] for d in range(len(shape))
        )
        row_flat: NDArray[np.intp] = np.ravel_multi_index(row_coords, shape)  # type: ignore[assignment]
        rows_list.append(row_flat)
        cols_list.append(arange_n)
        data_list.append(np.full(n, vals[k]))

    from tidal.solver._scipy_types import csc_matrix as _csc  # noqa: PLC0415

    return _csc(
        (
            np.concatenate(data_list),
            (np.concatenate(rows_list), np.concatenate(cols_list)),
        ),
        shape=(n, n),
    )


def build_operator_matrix(
    operator: str,
    grid: GridInfo,
    bc: BCSpec | None,
) -> SparseMatrix:
    """Build the N x N sparse matrix for a single spatial operator.

    For all-periodic BCs, uses a fast O(nnz) circulant construction
    (single-probe + tiling).  Otherwise falls back to O(N²) column-by-column
    probing.
    """
    from tidal.solver.operators import apply_operator, is_periodic_bc  # noqa: PLC0415

    # Fast path: all-periodic BCs → circulant (translation-invariant)
    if bc is not None and operator != "identity":
        if isinstance(bc, str):
            all_periodic = is_periodic_bc(bc)
        else:
            all_periodic = all(is_periodic_bc(b) for b in bc)  # type: ignore[union-attr]
        if all_periodic:
            return _build_operator_matrix_circulant(operator, grid, bc)

    # Fallback: probe each unit vector
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
                operator, self._grid, self._bc,
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
LAPLACIAN_OPS = frozenset(
    {
        "laplacian",
        "laplacian_x",
        "laplacian_y",
        "laplacian_z",
    },
)


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


def build_jacobian_matrices(  # noqa: C901, PLR0914
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
        eq.field_name for eq in spec.equations if eq.time_derivative_order == 0
    }
    no_self_term_fields = _detect_no_self_term_fields(spec)
    gauge_fix_fields = _detect_gauge_fix_fields(spec, grid, bc)

    # COO accumulation: append blocks as triples, single CSC conversion
    dF_dy = _COOAccumulator((n_total, n_total))  # noqa: N806
    dF_dyp = _COOAccumulator((n_total, n_total))  # noqa: N806

    for slot_idx, slot in enumerate(layout.slots):
        if slot.time_order == 0:
            _build_constraint_block(
                slot_idx,
                slot,
                layout,
                spec,
                eq_map,
                coeff_eval,
                op_cache,
                constraint_fields,
                no_self_term_fields,
                n,
                dF_dy,
                dF_dyp,
            )

        elif slot.kind == "velocity":
            # res = yp - RHS → dF/dyp[slot, slot] = I
            dF_dyp.add_block(slot_idx, slot_idx, n, I_mat)
            eq_idx = eq_map.get(slot.field_name)
            if eq_idx is not None:
                _add_rhs_terms(
                    slot_idx,
                    eq_idx,
                    spec,
                    layout,
                    coeff_eval,
                    op_cache,
                    constraint_fields,
                    n,
                    dF_dy,
                    dF_dyp,
                    negate=True,
                )

        elif slot.time_order >= 2 and slot.kind == "field":  # noqa: PLR2004
            # res = yp - v → dF/dyp[slot, slot] = I, dF/dy[slot, vel] = -I
            dF_dyp.add_block(slot_idx, slot_idx, n, I_mat)
            vel_slot = layout.velocity_slot_map.get(slot.field_name)
            if vel_slot is not None:
                dF_dy.add_block(slot_idx, vel_slot, n, -I_mat)

        elif slot.time_order == 1:
            # res = yp - RHS → dF/dyp[slot, slot] = I
            dF_dyp.add_block(slot_idx, slot_idx, n, I_mat)
            eq_idx = eq_map.get(slot.field_name)
            if eq_idx is not None:
                _add_rhs_terms(
                    slot_idx,
                    eq_idx,
                    spec,
                    layout,
                    coeff_eval,
                    op_cache,
                    constraint_fields,
                    n,
                    dF_dy,
                    dF_dyp,
                    negate=True,
                )

    result_dy = dF_dy.to_csc()
    result_dyp = dF_dyp.to_csc()

    # Post-process gauge-fix rows: zero out entire row, set pinning entry.
    # This is done on the final CSC because COO can't "overwrite" earlier
    # entries — it can only sum.
    for slot_idx, slot in enumerate(layout.slots):
        if slot.time_order == 0 and slot.field_name in gauge_fix_fields:
            field_slot = layout.field_slot_map[slot.field_name]
            global_row = slot_idx * n
            # Zero out row in both matrices
            result_dy[global_row, :] = 0  # type: ignore[index]
            result_dyp[global_row, :] = 0  # type: ignore[index]
            # Pinning entry: dF/dy[row, field_col] = 1
            result_dy[global_row, field_slot * n] = 1.0  # type: ignore[index]
            # Re-convert after row modification
    if gauge_fix_fields:
        result_dy = result_dy.tocsc()
        result_dyp = result_dyp.tocsc()

    return result_dy, result_dyp


# ---------------------------------------------------------------------------
# Block-assembly helpers
# ---------------------------------------------------------------------------


class _COOAccumulator:
    """Accumulate sparse blocks as COO triples for fast final assembly.

    Instead of LIL slice assignment (O(N²) per block due to dense
    intermediate), this appends COO triples and does a single CSC
    conversion at the end.  Duplicate entries are summed automatically
    by scipy's CSC constructor, so ``+=`` semantics come for free.
    """

    def __init__(self, shape: tuple[int, int]) -> None:
        self._shape = shape
        self._rows: list[NDArray[np.intp]] = []
        self._cols: list[NDArray[np.intp]] = []
        self._data: list[NDArray[np.float64]] = []

    def add_block(
        self,
        row_slot: int,
        col_slot: int,
        n: int,
        block: SparseMatrix,
    ) -> None:
        """Append a block at the ``(row_slot, col_slot)`` position."""
        coo = block.tocoo()
        self._rows.append(coo.row.astype(np.intp) + row_slot * n)
        self._cols.append(coo.col.astype(np.intp) + col_slot * n)
        self._data.append(np.asarray(coo.data, dtype=np.float64))

    def to_csc(self) -> SparseMatrix:
        """Convert accumulated triples to a CSC matrix."""
        from tidal.solver._scipy_types import csc_matrix as _csc  # noqa: PLC0415

        if not self._rows:
            return _csc(self._shape)
        rows = np.concatenate(self._rows)
        cols = np.concatenate(self._cols)
        data = np.concatenate(self._data)
        return _csc((data, (rows, cols)), shape=self._shape)


def _add_rhs_terms(  # noqa: PLR0913, PLR0917
    row_slot: int,
    eq_idx: int,
    spec: EquationSystem,
    layout: StateLayout,
    coeff_eval: CoefficientEvaluator,
    op_cache: _OperatorCache,
    constraint_fields: set[str],
    n: int,
    dF_dy: _COOAccumulator,  # noqa: N803
    dF_dyp: _COOAccumulator,  # noqa: N803
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
            term,
            layout,
            constraint_fields,
            op_cache,
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
        target_mat.add_block(row_slot, col_slot, n, scaled)


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
    n: int,
    dF_dy: _COOAccumulator,  # noqa: N803
    dF_dyp: _COOAccumulator,  # noqa: N803
) -> None:
    """Build Jacobian blocks for a constraint slot (time_order=0).

    Three sub-cases matching ida.py handle_constraint():
    1. No-self-term: res = y[field] → dF/dy = I
    2. Gauge-fix: normal RHS (gauge-fix row pinning applied post-hoc)
    3. Normal: res = RHS
    """
    eq_idx = eq_map.get(slot.field_name)
    if eq_idx is None:
        return

    # Case 1: no self-terms — field frozen at zero, res = y[field]
    if slot.field_name in no_self_term_fields:
        field_slot = layout.field_slot_map[slot.field_name]
        I_mat = op_cache.get_identity()  # noqa: N806
        dF_dy.add_block(slot_idx, field_slot, n, I_mat)
        return

    # Case 3 (and setup for case 2): normal RHS
    # res = RHS → dF/dy += coeff * op_mat (no negation)
    # Note: gauge-fix row pinning (case 2) is applied post-hoc in
    # build_jacobian_matrices() after COO→CSC conversion.
    _add_rhs_terms(
        slot_idx,
        eq_idx,
        spec,
        layout,
        coeff_eval,
        op_cache,
        constraint_fields,
        n,
        dF_dy,
        dF_dyp,
        negate=False,
    )


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
# Jacobian delivery: sparse jacfn (1D CSC data array)
# ---------------------------------------------------------------------------


def _extract_aligned_data(
    matrix: SparseMatrix,
    union: SparseMatrix,
) -> NDArray[np.float64]:
    """Extract *matrix* values aligned to *union*'s CSC data layout.

    Returns a 1D array of length ``union.nnz`` where each entry
    corresponds to the same ``(row, col)`` position as ``union.data[i]``.
    Positions present in *union* but absent in *matrix* are zero.

    This is necessary because sksundae's sparse ``jacfn`` callback
    receives a 1D array whose positions correspond to the sparsity
    pattern's CSC data layout.
    """
    mat_csc = matrix.tocsc()
    n_cols = union.shape[1]
    result = np.zeros(union.nnz, dtype=np.float64)

    for j in range(n_cols):
        # Union column j entries
        u_start, u_end = int(union.indptr[j]), int(union.indptr[j + 1])
        if u_start == u_end:
            continue

        # Matrix column j entries
        m_start, m_end = int(mat_csc.indptr[j]), int(mat_csc.indptr[j + 1])
        if m_start == m_end:
            continue
        m_rows = mat_csc.indices[m_start:m_end]
        m_vals = mat_csc.data[m_start:m_end]

        # Both are sorted within each column (CSC guarantee), use searchsorted
        u_rows = union.indices[u_start:u_end]
        idx: NDArray[np.intp] = np.searchsorted(u_rows, m_rows)  # type: ignore[assignment]
        result[u_start + idx] = m_vals

    return result


def _create_sparse_jacfn(
    dF_dy: SparseMatrix,  # noqa: N803
    dF_dyp: SparseMatrix,  # noqa: N803
) -> tuple[Callable[..., None], SparseMatrix]:
    """Create a sparse ``jacfn`` callback for IDA.

    Returns ``(jacfn, sparsity_pattern)`` where ``jacfn`` fills a 1D
    array ``JJ`` of length ``nnz`` in CSC column-compressed order, and
    ``sparsity_pattern`` is the CSC binary pattern to pass to IDA via
    ``options["sparsity"]``.

    The callback computes ``JJ[:] = dF_dy_data + cj * dF_dyp_data``
    using pre-extracted, structurally aligned data arrays.

    Requires sksundae >= 1.1.2 (fixes the ``aux.jacfn`` overwrite bug).
    """
    # Union nonzero structure — superset of both matrices
    union_csc = (abs(dF_dy) + abs(dF_dyp)).tocsc()
    union_csc.eliminate_zeros()

    # Binary sparsity pattern (same CSC structure, all 1.0)
    sparsity = union_csc.copy()
    sparsity.data[:] = 1.0

    # Extract data arrays aligned to union CSC layout
    dy_data = _extract_aligned_data(dF_dy, union_csc)
    dyp_data = _extract_aligned_data(dF_dyp, union_csc)
    dyp_scaled = np.empty_like(dyp_data)

    def jacfn(  # noqa: PLR0913, PLR0917
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        res: NDArray[np.float64],  # noqa: ARG001
        cj: float,
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        np.multiply(dyp_data, cj, out=dyp_scaled)
        np.add(dy_data, dyp_scaled, out=JJ)

    return jacfn, sparsity


def _create_cvode_sparse_jacfn(
    dF_dy: SparseMatrix,  # noqa: N803
) -> tuple[Callable[..., None], SparseMatrix]:
    """Create a sparse ``jacfn`` callback for CVODE.

    Returns ``(jacfn, sparsity_pattern)`` where ``jacfn`` fills a 1D
    CSC data array with the ODE Jacobian ``-dF_dy``.
    """
    neg_dy = (-dF_dy).tocsc()
    neg_dy.eliminate_zeros()

    sparsity = neg_dy.copy()
    sparsity.data[:] = 1.0

    ode_jac_data = neg_dy.data.copy()

    def jacfn(
        t: float,  # noqa: ARG001
        y: NDArray[np.float64],  # noqa: ARG001
        yp: NDArray[np.float64],  # noqa: ARG001
        JJ: NDArray[np.float64],  # noqa: N803
    ) -> None:
        JJ[:] = ode_jac_data

    return jacfn, sparsity


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


def _configure_gmres_tier(
    options: dict[str, Any],
    jac_y: SparseMatrix,
    jac_yp: SparseMatrix,
    solver: str,
    n_state: int,
) -> None:
    """Configure GMRES with analytical Jacobian-vector product."""
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
        solver,
        n_state,
    )


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

    Three delivery modes by system size:

    - **Dense** (N <= DENSE_THRESHOLD): 2D ``jacfn`` fills dense array.
    - **Sparse** (DENSE_THRESHOLD < N <= SPARSE_THRESHOLD): 1D CSC
      ``jacfn`` + SuperLU_MT direct factorisation.  Falls through to
      GMRES if nnz exceeds ``SUPERLU_NNZ_LIMIT``.
    - **GMRES** (N > SPARSE_THRESHOLD): ``jactimes`` provides analytical
      Jacobian-vector product for iterative GMRES.

    Parameters
    ----------
    solver : str
        ``"ida"`` or ``"cvode"``.  Controls the callback signature.
    """
    from tidal.solver._types import (  # noqa: PLC0415
        DENSE_THRESHOLD,
        SPARSE_THRESHOLD,
        SUPERLU_NNZ_LIMIT,
    )

    n_state = layout.total_size

    if not _is_system_time_independent(spec, grid, parameters):
        return False

    if n_state <= DENSE_THRESHOLD:
        # Dense tier: 2D jacfn
        jac_y, jac_yp = build_jacobian_matrices(
            spec,
            layout,
            grid,
            bc,
            parameters,
        )
        options["linsolver"] = "dense"
        if solver == "cvode":
            options["jacfn"] = _create_cvode_jacfn(jac_y)
        else:
            options["jacfn"] = _create_jacfn(jac_y, jac_yp)
        logger.info(
            "Analytical Jacobian (dense %s jacfn) for %d-state system",
            solver,
            n_state,
        )

    elif n_state <= SPARSE_THRESHOLD:
        # Sparse tier: 1D CSC jacfn + SuperLU_MT direct factorisation.
        # Eliminates O(n_colors) FD residual evaluations per Newton step.
        # Requires sksundae >= 1.1.2 (fixes aux.jacfn overwrite bug).

        # Cheap nnz pre-check: build_jacobian_sparsity is O(nnz) and
        # returns a conservative superset of the actual pattern.  If it
        # already exceeds SUPERLU_NNZ_LIMIT, skip straight to GMRES.
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        pattern_est = build_jacobian_sparsity(spec, layout, grid, bc)
        if pattern_est.nnz > SUPERLU_NNZ_LIMIT:
            jac_y, jac_yp = build_jacobian_matrices(
                spec,
                layout,
                grid,
                bc,
                parameters,
            )
            _configure_gmres_tier(options, jac_y, jac_yp, solver, n_state)
        else:
            jac_y, jac_yp = build_jacobian_matrices(
                spec,
                layout,
                grid,
                bc,
                parameters,
            )
            if solver == "cvode":
                jacfn_cb, sparsity = _create_cvode_sparse_jacfn(jac_y)
            else:
                jacfn_cb, sparsity = _create_sparse_jacfn(jac_y, jac_yp)
            options["linsolver"] = "sparse"
            options["sparsity"] = sparsity
            options["jacfn"] = jacfn_cb
            logger.info(
                "Analytical Jacobian (sparse %s jacfn, nnz=%d) for %d-state system",
                solver,
                sparsity.nnz,
                n_state,
            )

    else:
        # GMRES tier (N > SPARSE_THRESHOLD): analytical Jacobian-vector
        # product.  Eliminates O(n_colors) residual evaluations per GMRES
        # iteration compared to the FD GMRES path.
        jac_y, jac_yp = build_jacobian_matrices(
            spec,
            layout,
            grid,
            bc,
            parameters,
        )
        _configure_gmres_tier(options, jac_y, jac_yp, solver, n_state)

    return True
