"""Analytical Jacobian sparsity pattern for TIDAL IDA solver.

Builds the sparsity structure of the IDA Jacobian ``J = dF/dy + cj * dF/dyp``
directly from the equation spec and operator stencils — no numerical probing.
This avoids the O(N^2) cost of ``sksundae.jacband.j_pattern()`` (which
allocates a dense N x N array) and gives exact structure in milliseconds.

The sparsity pattern is returned as a ``scipy.sparse.csc_matrix`` suitable
for passing to ``sksundae.ida.IDA(linsolver='sparse', sparsity=pattern)``,
which uses SuperLU_MT for direct factorisation.

Reference: SuperLU_MT, Li & Demmel, ACM TOMS 29(2), 2003.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np
from scipy import sparse

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm

# Second-order threshold (same constant as ida.py / leapfrog.py)
_SECOND_ORDER = 2

# ---------------------------------------------------------------------------
# Operator -> stencil offsets
# ---------------------------------------------------------------------------

# Each operator maps to a set of relative grid-coordinate offsets that
# describe which neighbor grid points contribute to the output at a given
# point.  These are used to populate the Jacobian sparsity blocks.
#
# Convention: offset tuple length = spatial_dimension (ndim).
# For operators that work per-axis, ``_axis_stencil`` builds the full
# offset tuple from the axis index.

_AXIS_MAP = {"x": 0, "y": 1, "z": 2}


def _axis_offsets(axis: int, ndim: int, deltas: list[int]) -> list[tuple[int, ...]]:
    """Build N-D offset tuples for *deltas* along a single *axis*."""
    offsets: list[tuple[int, ...]] = []
    for d in deltas:
        o = [0] * ndim
        o[axis] = d
        offsets.append(tuple(o))
    return offsets


def _cross_offsets(
    ax1: int, ax2: int, ndim: int,
) -> list[tuple[int, ...]]:
    """Build offsets for cross-derivative d^2/(dx_i dx_j).

    Stencil: (f[i+1,j+1] - f[i+1,j-1] - f[i-1,j+1] + f[i-1,j-1]) / (4 dx dy)
    Plus self at origin (conservative -- covers the cj diagonal contribution).
    """
    offsets: list[tuple[int, ...]] = []
    for d1 in (-1, 1):
        for d2 in (-1, 1):
            o = [0] * ndim
            o[ax1] = d1
            o[ax2] = d2
            offsets.append(tuple(o))
    # Include self (origin) for robustness
    offsets.append((0,) * ndim)
    return offsets


def _biharmonic_offsets(ndim: int) -> list[tuple[int, ...]]:
    r"""Build offsets for biharmonic nabla^4 f = nabla^2(nabla^2 f).

    The biharmonic stencil is the convolution of the Laplacian stencil with
    itself.  For 1D: 5-point stencil {-2, -1, 0, 1, 2}.  For 2D, all
    (dx, dy) with \|dx\| + \|dy\| <= 2.
    """
    # Conservative: all offsets within L1-radius 2
    offsets_set: set[tuple[int, ...]] = set()
    for combo in itertools.product(range(-2, 3), repeat=ndim):
        if sum(abs(c) for c in combo) <= 2:  # noqa: PLR2004
            offsets_set.add(combo)
    return list(offsets_set)


def _gradient_or_laplacian_offsets(name: str, ndim: int) -> list[tuple[int, ...]] | None:
    """Try to resolve gradient_X or laplacian_X operators."""
    if name.startswith("gradient_"):
        ax = _AXIS_MAP.get(name[-1])
        if ax is not None:
            return _axis_offsets(ax, ndim, [-1, 0, 1])

    if name.startswith("laplacian_"):
        parts = name.split("_")
        ax = _AXIS_MAP.get(parts[1]) if len(parts) > 1 else None
        if ax is not None:
            return _axis_offsets(ax, ndim, [-1, 0, 1])

    return None


def _full_laplacian_offsets(ndim: int) -> list[tuple[int, ...]]:
    """Full Laplacian = union of per-axis stencils."""
    zero = (0,) * ndim
    offsets_set: set[tuple[int, ...]] = {zero}
    for ax in range(ndim):
        offsets_set.update(_axis_offsets(ax, ndim, [-1, 1]))
    return list(offsets_set)


def _cross_derivative_offsets(name: str, ndim: int) -> list[tuple[int, ...]] | None:
    """Try to resolve cross_derivative_XY operators."""
    if not name.startswith("cross_derivative_"):
        return None
    suffix = name.rsplit("_", maxsplit=1)[-1]  # "xy", "xz", "yz"
    ax1 = _AXIS_MAP.get(suffix[0])
    ax2 = _AXIS_MAP.get(suffix[1]) if len(suffix) > 1 else None
    if ax1 is not None and ax2 is not None:
        return _cross_offsets(ax1, ax2, ndim)
    return None


def operator_stencil_offsets(name: str, ndim: int) -> list[tuple[int, ...]]:
    """Return the grid-coordinate offsets for operator *name*.

    Parameters
    ----------
    name : str
        Operator name (e.g. ``"laplacian"``, ``"gradient_x"``).
    ndim : int
        Number of spatial dimensions.

    Returns
    -------
    list[tuple[int, ...]]
        Offset tuples, each of length *ndim*.
    """
    zero = (0,) * ndim

    # Identity and first_derivative_t: self-coupling only
    if name in {"identity", "first_derivative_t"}:
        return [zero]

    # Directional gradient or Laplacian
    result = _gradient_or_laplacian_offsets(name, ndim)
    if result is not None:
        return result

    # Full Laplacian
    if name == "laplacian":
        return _full_laplacian_offsets(ndim)

    # Cross derivative
    result = _cross_derivative_offsets(name, ndim)
    if result is not None:
        return result

    # Biharmonic
    if name == "biharmonic":
        return _biharmonic_offsets(ndim)

    # Unknown operator -- conservative: assume full local coupling (radius 1)
    return list(itertools.product(range(-1, 2), repeat=ndim))


# ---------------------------------------------------------------------------
# Stencil -> sparse block construction (vectorized)
# ---------------------------------------------------------------------------


def _stencil_block(
    n_points: int,
    grid_shape: tuple[int, ...],
    offsets: list[tuple[int, ...]],
    periodic: tuple[bool, ...],
) -> tuple[np.ndarray, np.ndarray]:
    """Build (row, col) index arrays for a single operator stencil block.

    Parameters
    ----------
    n_points : int
        Total grid points (product of grid_shape).
    grid_shape : tuple[int, ...]
        Shape of the spatial grid.
    offsets : list[tuple[int, ...]]
        Grid-coordinate offsets from :func:`operator_stencil_offsets`.
    periodic : tuple[bool, ...]
        Per-axis periodic flags.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (rows, cols) -- flat grid indices (0-based, within a single slot).
    """
    ndim = len(grid_shape)
    flat = np.arange(n_points)
    # Unravel to multi-index: shape (ndim, n_points)
    multi = np.array(np.unravel_index(flat, grid_shape))

    all_rows: list[np.ndarray] = []
    all_cols: list[np.ndarray] = []

    for offset in offsets:
        neighbor = multi.copy()
        valid = np.ones(n_points, dtype=bool)

        for ax in range(ndim):
            shifted = neighbor[ax] + offset[ax]
            if periodic[ax]:
                shifted %= grid_shape[ax]
            else:
                # Non-periodic: clip and mark out-of-bounds
                oob = (shifted < 0) | (shifted >= grid_shape[ax])
                valid &= ~oob
                shifted = np.clip(shifted, 0, grid_shape[ax] - 1)
            neighbor[ax] = shifted

        col_flat = np.asarray(np.ravel_multi_index(tuple(neighbor), grid_shape))

        if not valid.all():
            all_rows.append(flat[valid])
            all_cols.append(col_flat[valid])
        else:
            all_rows.append(flat)
            all_cols.append(col_flat)

    if not all_rows:
        return np.array([], dtype=np.intp), np.array([], dtype=np.intp)
    return np.concatenate(all_rows), np.concatenate(all_cols)


# ---------------------------------------------------------------------------
# Slot-coupling analysis
# ---------------------------------------------------------------------------


def _find_slot_for_field(layout: StateLayout, field_ref: str) -> int | None:
    """Map a field reference (e.g. ``"A_0"`` or ``"pi_1"``) to a slot index.

    The field_ref may be:
    - A field name (``"A_0"``) -> field slot
    - A momentum name (``"pi_A_1"`` or ``"pi_1"``) -> momentum slot
    """
    # Direct match on slot name
    for i, slot in enumerate(layout.slots):
        if slot.name == field_ref:
            return i

    # Try as field name in field_slot_map
    if field_ref in layout.field_slot_map:
        return layout.field_slot_map[field_ref]

    # Try as momentum reference: "pi_X" -> momentum slot for field X
    if field_ref.startswith("pi_"):
        base_field = field_ref[3:]  # Strip "pi_" prefix
        if base_field in layout.momentum_slot_map:
            return layout.momentum_slot_map[base_field]
        # Also try finding a field whose momentum slot name matches
        for slot_idx in layout.momentum_slot_map.values():
            mom_slot = layout.slots[slot_idx]
            if mom_slot.name == field_ref:
                return slot_idx

    return None


def _rhs_couplings(
    terms: tuple[OperatorTerm, ...],
    layout: StateLayout,
    ndim: int,
) -> list[tuple[int, list[tuple[int, ...]]]]:
    """Extract (col_slot_idx, stencil_offsets) for each RHS term.

    For ``first_derivative_t(field)``, the coupling target is the momentum
    slot ``pi_{field}`` with identity stencil.
    """
    couplings: list[tuple[int, list[tuple[int, ...]]]] = []
    for term in terms:
        if term.operator == "first_derivative_t":
            # Couples to momentum of the referenced field
            pi_name = f"pi_{term.field}"
            col_slot = _find_slot_for_field(layout, pi_name)
            if col_slot is not None:
                couplings.append((col_slot, [((0,) * ndim)]))
        else:
            col_slot = _find_slot_for_field(layout, term.field)
            if col_slot is not None:
                offsets = operator_stencil_offsets(term.operator, ndim)
                couplings.append((col_slot, offsets))
    return couplings


# ---------------------------------------------------------------------------
# Per-slot-type sparsity helpers
# ---------------------------------------------------------------------------


class _SparsityBuilder:
    """Accumulates (row, col) pairs for the full Jacobian sparsity pattern."""

    def __init__(
        self,
        spec: EquationSystem,
        layout: StateLayout,
        grid: GridInfo,
    ) -> None:
        self.spec = spec
        self.layout = layout
        self.grid = grid
        self.n = grid.num_points
        self.ndim = grid.ndim
        self.canonical = spec.canonical
        self.eq_map: dict[str, int] = {
            eq.field_name: i for i, eq in enumerate(spec.equations)
        }
        self.all_rows: list[np.ndarray] = []
        self.all_cols: list[np.ndarray] = []

    def add_block(
        self,
        row_slot: int,
        col_slot: int,
        offsets: list[tuple[int, ...]],
    ) -> None:
        """Add a stencil coupling block between two slots."""
        r, c = _stencil_block(self.n, self.grid.shape, offsets, self.grid.periodic)
        self.all_rows.append(row_slot * self.n + r)
        self.all_cols.append(col_slot * self.n + c)

    def add_diagonal(self, slot_idx: int) -> None:
        """Add identity (diagonal) block for a slot."""
        idx = np.arange(slot_idx * self.n, (slot_idx + 1) * self.n)
        self.all_rows.append(idx)
        self.all_cols.append(idx)

    def add_rhs_couplings(self, row_slot: int, eq_idx: int) -> None:
        """Add coupling from all RHS terms of equation *eq_idx*."""
        eq = self.spec.equations[eq_idx]
        for col_slot, offsets in _rhs_couplings(eq.rhs_terms, self.layout, self.ndim):
            self.add_block(row_slot, col_slot, offsets)

    def handle_constraint(self, slot_idx: int, eq_idx: int | None) -> None:
        """Constraint: res = RHS(y) -- no yp coupling, no cj*I diagonal.

        If the equation has no self-referencing terms, the IDA residual
        is ``res = y[field]`` (freeze at zero) → diagonal Jacobian block.
        """
        if eq_idx is not None:
            eq = self.spec.equations[eq_idx]
            has_self = any(t.field == eq.field_name for t in eq.rhs_terms)
            if not has_self:
                # No self-terms: residual is y[field] → identity coupling
                self.add_diagonal(slot_idx)
                return
            self.add_rhs_couplings(slot_idx, eq_idx)

    def handle_momentum(self, slot_idx: int, eq_idx: int | None) -> None:
        """Momentum: res = yp - RHS(y) -- diagonal from cj*I + RHS coupling."""
        self.add_diagonal(slot_idx)
        if eq_idx is not None:
            self.add_rhs_couplings(slot_idx, eq_idx)

    def handle_dynamical_field(self, slot_idx: int, field_name: str, dyn_i: int | None) -> None:
        """Dynamical field: res = K*yp - (pi - S)."""
        zero = (0,) * self.ndim
        canonical = self.canonical

        # yp coupling from kinetic matrix
        self._add_kinetic_coupling(slot_idx, dyn_i)

        # y coupling: momentum slot (pi_i) via identity
        if field_name in self.layout.momentum_slot_map:
            mom_slot = self.layout.momentum_slot_map[field_name]
            self.add_block(slot_idx, mom_slot, [zero])

        # y coupling: spatial_momenta S_i
        if canonical and canonical.spatial_momenta:
            sm_terms = canonical.spatial_momenta.get(field_name)
            if sm_terms:
                for col_slot, offsets in _rhs_couplings(sm_terms, self.layout, self.ndim):
                    self.add_block(slot_idx, col_slot, offsets)

        # Fallback: field_rates coupling
        if canonical and field_name in canonical.field_rates:
            fr_terms = canonical.field_rates[field_name]
            for col_slot, offsets in _rhs_couplings(fr_terms, self.layout, self.ndim):
                self.add_block(slot_idx, col_slot, offsets)

    def _add_kinetic_coupling(self, slot_idx: int, dyn_i: int | None) -> None:
        """Add yp coupling from kinetic matrix K_{ij} * yp_j."""
        zero = (0,) * self.ndim
        canonical = self.canonical
        if dyn_i is not None and canonical and canonical.kinetic_matrix:
            k_dense = canonical.kinetic_matrix.to_dense()
            for j, dyn_field in enumerate(self.layout.dynamical_fields):
                if k_dense[dyn_i][j] != 0:
                    fs = self.layout.field_slot_map[dyn_field]
                    self.add_diagonal(slot_idx)  # Self from cj*I
                    if fs != slot_idx:
                        self.add_block(slot_idx, fs, [zero])
        else:
            self.add_diagonal(slot_idx)  # cj*I from yp

    def handle_first_order(self, slot_idx: int, eq_idx: int | None) -> None:
        """First-order: res = yp - RHS(y) -- diagonal + RHS."""
        self.add_diagonal(slot_idx)
        if eq_idx is not None:
            self.add_rhs_couplings(slot_idx, eq_idx)

    def assemble(self) -> sparse.csc_matrix:
        """Assemble accumulated entries into a CSC sparse matrix."""
        n_total = self.layout.total_size
        if not self.all_rows:
            return sparse.csc_matrix((n_total, n_total))

        rows = np.concatenate(self.all_rows)
        cols = np.concatenate(self.all_cols)
        data = np.ones(len(rows), dtype=np.float64)

        # Deduplicate via COO -> CSC conversion (scipy handles duplicates)
        pattern = sparse.coo_matrix((data, (rows, cols)), shape=(n_total, n_total))
        return sparse.csc_matrix(pattern)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------


def build_jacobian_sparsity(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    _bc: str | tuple[str, ...] | None,
) -> sparse.csc_matrix:
    """Build the analytical Jacobian sparsity pattern for IDA.

    Constructs the nonzero structure of ``J = dF/dy + cj * dF/dyp`` from
    the equation system's operator stencils and slot-coupling structure.
    No numerical probing -- the pattern is exact and built in O(nnz) time.

    Parameters
    ----------
    spec : EquationSystem
        Parsed equation specification.
    layout : StateLayout
        State vector layout.
    grid : GridInfo
        Spatial grid (periodic flags used for wrapping).
    _bc : str or tuple of str, optional
        Boundary conditions string (kept for API compatibility with
        ``ida.py`` call site; periodicity is read from ``grid.periodic``).

    Returns
    -------
    scipy.sparse.csc_matrix
        Binary sparsity pattern of shape ``(N_total, N_total)``.
    """
    builder = _SparsityBuilder(spec, layout, grid)

    for slot_idx, slot in enumerate(layout.slots):
        eq_idx = builder.eq_map.get(slot.field_name)

        if slot.time_order == 0:
            builder.handle_constraint(slot_idx, eq_idx)
        elif slot.kind == "momentum":
            builder.handle_momentum(slot_idx, eq_idx)
        elif slot.time_order >= _SECOND_ORDER and slot.kind == "field":
            builder.handle_dynamical_field(slot_idx, slot.field_name, slot.dynamical_index)
        elif slot.time_order == 1:
            builder.handle_first_order(slot_idx, eq_idx)

    return builder.assemble()
