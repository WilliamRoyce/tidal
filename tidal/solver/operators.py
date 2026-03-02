"""Spatial finite-difference operators on plain numpy arrays.

All operators accept an ``np.ndarray`` (field data) and a ``GridInfo``
descriptor, returning an ``np.ndarray`` of the same shape.  No py-pde
ScalarField / VectorField wrappers — pure numpy with explicit ghost-cell
boundary handling.

Boundary condition types
------------------------
- ``"periodic"`` — wraparound via ``np.roll``
- ``"neumann"``  — zero normal derivative (mirror ghost cell)
- ``"dirichlet"``— zero value (antisymmetric ghost cell)
- ``"robin"``    -- mixed: d_n f + gamma*f = beta

All non-periodic BCs are unified via the ghost-cell formula
``ghost = const + factor * interior_cell`` where ``(const, factor)``
are determined by BC type, value, and grid spacing.  See ``SideBCSpec``
for the derivation.

The ``bc`` parameter may be:
- A single string applied to all axes (e.g. ``"periodic"``)
- A tuple of strings, one per axis (e.g. ``("neumann", "periodic")``)
- A tuple of ``AxisBCSpec`` objects for per-side / non-zero / Robin BCs

Reference: py-pde's ghost-cell virtual-point approach
(David Zwicker, J. Open Source Software 5(48), 2020).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo

# ---------------------------------------------------------------------------
# BC data model
# ---------------------------------------------------------------------------

_VALID_BC_KINDS = frozenset({"dirichlet", "neumann", "robin"})
_VALID_BC = frozenset({"periodic", "neumann", "dirichlet", "robin"})


@dataclass(frozen=True)
class SideBCSpec:
    """Boundary condition for one side (low or high) of one axis.

    All first-order BCs map to a ghost-cell formula::

        ghost = const + factor * interior_cell

    where ``(const, factor)`` depend on the BC kind:

    - **Neumann** (df/dn = D):  const = dx*D,  factor = +1
    - **Dirichlet** (f = V):    const = 2*V,    factor = -1
    - **Robin** (d_n f + gamma*f = beta):
      const = 2*dx*beta / (gamma*dx + 2),
      factor = (2 - gamma*dx) / (gamma*dx + 2)

    Setting D=0, V=0 recovers the homogeneous cases.
    """

    kind: str
    """BC type: ``"dirichlet"``, ``"neumann"``, or ``"robin"``."""
    value: float = 0.0
    """Dirichlet boundary value V, or Robin inhomogeneity β."""
    derivative: float = 0.0
    """Neumann ghost-cell offset D.

    Ghost cell is placed at ``interior + dx * D``.  For homogeneous Neumann
    (zero normal flux), set D=0.  The sign is the same at both low and high
    boundaries (py-pde convention; see Zwicker, JOSS 2020).
    """
    gamma: float = 0.0
    """Robin coefficient gamma in d_n f + gamma*f = beta. Must be >= 0."""

    def __post_init__(self) -> None:
        if self.kind not in _VALID_BC_KINDS:
            msg = (
                f"Invalid BC kind: {self.kind!r}. "
                f"Must be one of {sorted(_VALID_BC_KINDS)}"
            )
            raise ValueError(msg)
        if self.kind == "robin" and self.gamma < 0:
            msg = f"Robin gamma must be ≥ 0, got {self.gamma}"
            raise ValueError(msg)

    def ghost_params(self, dx: float) -> tuple[float, float]:
        """Return ``(const, factor)`` for ``ghost = const + factor * interior``.

        Parameters
        ----------
        dx : float
            Grid spacing in the normal direction.
        """
        if self.kind == "neumann":
            return (dx * self.derivative, 1.0)
        if self.kind == "dirichlet":
            return (2.0 * self.value, -1.0)
        # Robin: d_n f + gamma*f = beta
        # When gamma*dx >= 2, factor <= 0 (unstable); see check_robin_stability().
        denom = self.gamma * dx + 2.0
        return (
            2.0 * dx * self.value / denom,
            (2.0 - self.gamma * dx) / denom,
        )


@dataclass(frozen=True)
class AxisBCSpec:
    """Boundary condition for one axis: periodic or distinct low/high sides.

    For periodic axes, ``low`` and ``high`` must be ``None``.
    For non-periodic axes, both must be ``SideBCSpec`` instances.
    """

    periodic: bool = False
    low: SideBCSpec | None = None
    high: SideBCSpec | None = None

    def __post_init__(self) -> None:
        if self.periodic:
            if self.low is not None or self.high is not None:
                msg = "Periodic axis cannot have low/high BC specs"
                raise ValueError(msg)
        elif self.low is None or self.high is None:
            msg = "Non-periodic axis requires both low and high BC specs"
            raise ValueError(msg)


# Singleton cache for _str_to_axis_bc — avoids repeated AxisBCSpec/SideBCSpec
# allocation for the same BC string.  Safe because AxisBCSpec is frozen.
_AXIS_BC_CACHE: dict[str, AxisBCSpec] = {}


def _str_to_axis_bc(bc_str: str) -> AxisBCSpec:
    """Convert a legacy BC string to an ``AxisBCSpec``.

    ``"periodic"`` maps to a periodic axis.  ``"neumann"``/``"dirichlet"``
    map to symmetric homogeneous BC on both sides.  ``"robin"`` maps to
    zero Robin on both sides (gamma=0, beta=0, equivalent to Neumann).

    Results are cached (singleton pattern) since ``AxisBCSpec`` is frozen.

    Raises
    ------
    ValueError
        If *bc_str* is not a recognized BC type.
    """
    cached = _AXIS_BC_CACHE.get(bc_str)
    if cached is not None:
        return cached
    if bc_str == "periodic":
        result = AxisBCSpec(periodic=True)
    elif bc_str not in _VALID_BC:
        msg = f"Unknown BC type: {bc_str!r}; valid: {sorted(_VALID_BC)}"
        raise ValueError(msg)
    else:
        side = SideBCSpec(kind=bc_str)
        result = AxisBCSpec(periodic=False, low=side, high=side)
    _AXIS_BC_CACHE[bc_str] = result
    return result


def is_periodic_bc(bc_entry: str | AxisBCSpec) -> bool:
    """Check whether a single axis BC entry is periodic.

    Handles both legacy string BCs and structured ``AxisBCSpec`` objects,
    so callers don't need to type-check manually.
    """
    if isinstance(bc_entry, str):
        return bc_entry == "periodic"
    return bc_entry.periodic


# Legacy type alias — widened to accept structured specs
BCSpec = str | tuple[str, ...] | tuple[AxisBCSpec, ...]
"""Single BC string, per-axis strings, or per-axis ``AxisBCSpec`` objects."""


def _normalize_bc(bc: BCSpec, grid: GridInfo) -> tuple[str | AxisBCSpec, ...]:
    """Expand a BC spec to a per-axis tuple, validating each entry.

    Returns a tuple of either strings or ``AxisBCSpec`` objects (never mixed
    within one call — the type depends on what was passed in).

    Pre-normalized tuples (from a previous ``_normalize_bc`` call) pass
    through with only a fast length check, avoiding per-element validation.

    Raises
    ------
    ValueError
        If the number of BC entries doesn't match ``grid.ndim`` or an
        unknown BC type is encountered.
    TypeError
        If a BC entry is neither a string nor an ``AxisBCSpec``.
    """
    # Fast path: already a tuple of correct length (pre-normalized)
    if isinstance(bc, tuple) and len(bc) == grid.ndim:
        return bc

    if isinstance(bc, str):
        bcs: tuple[str | AxisBCSpec, ...] = (bc,) * grid.ndim
    else:
        bcs = tuple(bc)

    if len(bcs) != grid.ndim:
        msg = f"Expected {grid.ndim} BC entries, got {len(bcs)}"
        raise ValueError(msg)

    for i, b in enumerate(bcs):
        if isinstance(b, str):
            if b not in _VALID_BC:
                msg = f"Unknown BC type {b!r} for axis {i}; valid: {sorted(_VALID_BC)}"
                raise ValueError(msg)
        elif not isinstance(b, AxisBCSpec):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"BC entry {i} must be str or AxisBCSpec, got {type(b).__name__}"
            raise TypeError(msg)

    return bcs


def _bc_from_grid(grid: GridInfo) -> tuple[str, ...]:
    """Infer BCs from grid periodicity (periodic or neumann)."""
    return tuple("periodic" if p else "neumann" for p in grid.periodic)


def _resolve_axis_bc(bc_entry: str | AxisBCSpec) -> AxisBCSpec:
    """Convert a per-axis BC entry to an ``AxisBCSpec``.

    Strings are converted via ``_str_to_axis_bc``.  ``AxisBCSpec`` objects
    pass through unchanged.
    """
    if isinstance(bc_entry, AxisBCSpec):
        return bc_entry
    return _str_to_axis_bc(bc_entry)


# ---------------------------------------------------------------------------
# Cached slice-tuple helpers (avoid per-call list + tuple allocation)
# ---------------------------------------------------------------------------

# Module-level constants for common slices
_SLC_LO = slice(0, -2)   # ghost-cell frame: left neighbor
_SLC_MID = slice(1, -1)  # ghost-cell frame: center
_SLC_HI = slice(2, None)  # ghost-cell frame: right neighbor
_SLC_GHOST_L = slice(0, 1)       # left ghost cell
_SLC_INTERIOR = slice(1, -1)     # interior (for writing)


def _slice_tuple(ndim: int, axis: int, axis_slice: slice) -> tuple[slice, ...]:
    """Build a tuple of slices with *axis_slice* at position *axis*."""
    s = [slice(None)] * ndim
    s[axis] = axis_slice
    return tuple(s)


# Pre-computed cache for the most common (ndim, axis, slice) combos.
# Populated lazily on first use per combo.  Thread-safe for reads (dict).
_SLICE_CACHE: dict[tuple[int, int, int], tuple[slice, ...]] = {}

# Map slice objects to small ints for cache keys
_SLC_ID = {id(_SLC_LO): 0, id(_SLC_MID): 1, id(_SLC_HI): 2}


def _cached_slice(ndim: int, axis: int, slc: slice) -> tuple[slice, ...]:
    """Return cached slice tuple for (ndim, axis, slc).

    Falls back to ``_slice_tuple`` for non-standard slices.
    """
    slc_id = _SLC_ID.get(id(slc))
    if slc_id is None:
        return _slice_tuple(ndim, axis, slc)
    key = (ndim, axis, slc_id)
    cached = _SLICE_CACHE.get(key)
    if cached is not None:
        return cached
    result = _slice_tuple(ndim, axis, slc)
    _SLICE_CACHE[key] = result
    return result


def _cached_slice_2(
    ndim: int, axis1: int, slc1: slice, axis2: int, slc2: slice,
) -> tuple[slice, ...]:
    """Build a 2-axis slice tuple for cross_derivative corners."""
    s = [slice(None)] * ndim
    s[axis1] = slc1
    s[axis2] = slc2
    return tuple(s)


# Cross-derivative corner cache: (ndim, axis1, axis2, slc1_id, slc2_id)
_CORNER_CACHE: dict[tuple[int, int, int, int, int], tuple[slice, ...]] = {}


def _cached_corner(
    ndim: int, axis1: int, slc1: slice, axis2: int, slc2: slice,
) -> tuple[slice, ...]:
    """Return cached corner slice for cross_derivative (2-axis)."""
    s1 = _SLC_ID.get(id(slc1), -1)
    s2 = _SLC_ID.get(id(slc2), -1)
    if s1 < 0 or s2 < 0:
        return _cached_slice_2(ndim, axis1, slc1, axis2, slc2)
    key = (ndim, axis1, axis2, s1, s2)
    cached = _CORNER_CACHE.get(key)
    if cached is not None:
        return cached
    result = _cached_slice_2(ndim, axis1, slc1, axis2, slc2)
    _CORNER_CACHE[key] = result
    return result


# ---------------------------------------------------------------------------
# Ghost-cell padding helpers
# ---------------------------------------------------------------------------


class _PadEntry:
    """Pre-allocated buffer and cached slices for one (shape, axis) combo."""

    __slots__ = ("buf", "interior_slc", "left_slc", "right_slc", "src_hi", "src_lo")

    def __init__(self, data_shape: tuple[int, ...], axis: int) -> None:
        n = data_shape[axis]
        ndim = len(data_shape)
        padded_shape = list(data_shape)
        padded_shape[axis] += 2
        self.buf = np.empty(tuple(padded_shape))
        # Cache all slice tuples (avoids per-call list+tuple creation)
        self.interior_slc = _slice_tuple(ndim, axis, slice(1, n + 1))
        self.left_slc = _slice_tuple(ndim, axis, slice(0, 1))
        self.right_slc = _slice_tuple(ndim, axis, slice(n + 1, n + 2))
        self.src_lo = _slice_tuple(ndim, axis, slice(0, 1))
        self.src_hi = _slice_tuple(ndim, axis, slice(n - 1, n))


class _PadBufferCache:
    """Lazily allocate and reuse padded buffers for ghost-cell operations.

    Keyed by ``(data.shape, axis)`` — since TIDAL operates on a fixed grid,
    each key is allocated only once.  Also caches the slice tuples needed
    for ghost-cell writes.  Single-threaded assumption (no locking).
    """

    __slots__ = ("_cache",)

    def __init__(self) -> None:
        self._cache: dict[tuple[tuple[int, ...], int], _PadEntry] = {}

    def get(self, data_shape: tuple[int, ...], axis: int) -> _PadEntry:
        key = (data_shape, axis)
        entry = self._cache.get(key)
        if entry is None:
            entry = _PadEntry(data_shape, axis)
            self._cache[key] = entry
        return entry


_pad_cache = _PadBufferCache()


def _pad_axis(
    data: np.ndarray,
    axis: int,
    bc: str | AxisBCSpec,
    dx: float = 1.0,
) -> np.ndarray:
    """Pad *data* with one ghost cell on each side along *axis*.

    Uses pre-allocated buffers and cached slice tuples (via ``_pad_cache``)
    for zero-alloc ghost-cell writes in the hot path.

    Parameters
    ----------
    data : np.ndarray
        Field data to pad.
    axis : int
        Axis along which to pad.
    bc : str or AxisBCSpec
        Boundary condition for this axis.
    dx : float
        Grid spacing along this axis (needed for Neumann derivative and
        Robin formula; default 1.0 for backward compat).
    """
    if isinstance(bc, str):
        bc = _str_to_axis_bc(bc)

    entry = _pad_cache.get(data.shape, axis)
    buf = entry.buf

    # Copy interior into buffer
    buf[entry.interior_slc] = data

    if bc.periodic:
        buf[entry.left_slc] = data[entry.src_hi]
        buf[entry.right_slc] = data[entry.src_lo]
    else:
        assert bc.low is not None  # guaranteed by AxisBCSpec.__post_init__
        assert bc.high is not None
        c_lo, f_lo = bc.low.ghost_params(dx)
        c_hi, f_hi = bc.high.ghost_params(dx)

        # Write ghost cells in-place: ghost = const + factor * interior
        np.multiply(f_lo, data[entry.src_lo], out=buf[entry.left_slc])
        if c_lo != 0.0:
            buf[entry.left_slc] += c_lo
        np.multiply(f_hi, data[entry.src_hi], out=buf[entry.right_slc])
        if c_hi != 0.0:
            buf[entry.right_slc] += c_hi

    return buf


# ---------------------------------------------------------------------------
# Core stencil operations
# ---------------------------------------------------------------------------


def gradient(
    data: np.ndarray,
    axis: int,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Central-difference gradient ∂f/∂x_i.

    Stencil: ``(f[i+1] - f[i-1]) / (2·dx)``

    Parameters
    ----------
    data : np.ndarray
        Field values, shape must match ``grid.shape``.
    axis : int
        Spatial axis along which to differentiate.
    grid : GridInfo
        Grid descriptor.
    bc : str or tuple of str, optional
        Boundary conditions.  Defaults to grid periodicity inference.

    Returns
    -------
    np.ndarray
        Gradient array, same shape as *data*.
    """
    bcs = _normalize_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    dx = grid.dx[axis]
    bc_axis = bcs[axis]

    # Inline _resolve_axis_bc to avoid function call overhead
    axis_bc = bc_axis if isinstance(bc_axis, AxisBCSpec) else _str_to_axis_bc(bc_axis)
    padded = _pad_axis(data, axis, axis_bc, dx)
    slc_left = _cached_slice(data.ndim, axis, _SLC_LO)
    slc_right = _cached_slice(data.ndim, axis, _SLC_HI)
    result = np.subtract(padded[slc_right], padded[slc_left])
    result *= 1.0 / (2.0 * dx)
    return result


def _directional_laplacian_raw(
    data: np.ndarray,
    axis: int,
    grid: GridInfo,
    bc_axis: str | AxisBCSpec,
) -> np.ndarray:
    """Core 3-point second derivative with pre-resolved per-axis BC.

    Used by ``directional_laplacian`` and ``laplacian`` to avoid
    redundant BC normalization when computing multiple axes.

    Uses ghost-cell padding with cached slice tuples and in-place arithmetic.
    """
    dx = grid.dx[axis]
    inv_dx2 = 1.0 / (dx * dx)

    # Inline _resolve_axis_bc to avoid function call overhead
    axis_bc = bc_axis if isinstance(bc_axis, AxisBCSpec) else _str_to_axis_bc(bc_axis)
    padded = _pad_axis(data, axis, axis_bc, dx)
    slc_left = _cached_slice(data.ndim, axis, _SLC_LO)
    slc_center = _cached_slice(data.ndim, axis, _SLC_MID)
    slc_right = _cached_slice(data.ndim, axis, _SLC_HI)
    # In-place arithmetic: result = (right - 2*center + left) * inv_dx2
    result = np.subtract(padded[slc_right], padded[slc_center])
    result -= padded[slc_center]
    result += padded[slc_left]
    result *= inv_dx2
    return result


def directional_laplacian(
    data: np.ndarray,
    axis: int,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """3-point second derivative ∂²f/∂x_i².

    Stencil: ``(f[i+1] - 2·f[i] + f[i-1]) / dx²``
    """
    bcs = _normalize_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    return _directional_laplacian_raw(data, axis, grid, bcs[axis])


def laplacian(
    data: np.ndarray,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Full Laplacian ∇²f = Σ_i ∂²f/∂x_i².

    Normalizes BCs once and fuses per-axis directional Laplacians
    to avoid redundant BC validation and resolution.
    """
    bcs = _normalize_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    result = _directional_laplacian_raw(data, 0, grid, bcs[0])
    for ax in range(1, grid.ndim):
        result += _directional_laplacian_raw(data, ax, grid, bcs[ax])
    return result


def cross_derivative(
    data: np.ndarray,
    axis1: int,
    axis2: int,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Mixed second derivative ∂²f/(∂x_i ∂x_j).

    Uses a fused 4-point stencil with ghost-cell padding and cached slice
    tuples for all BC types.
    """
    bcs = _normalize_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    # Inline _resolve_axis_bc
    bc_e1 = bcs[axis1]
    bc_e2 = bcs[axis2]
    bc1 = bc_e1 if isinstance(bc_e1, AxisBCSpec) else _str_to_axis_bc(bc_e1)
    bc2 = bc_e2 if isinstance(bc_e2, AxisBCSpec) else _str_to_axis_bc(bc_e2)

    dx = grid.dx[axis1]
    dy = grid.dx[axis2]

    # Pad along axis1, then axis2 on the padded result
    padded1 = _pad_axis(data, axis1, bc1, dx)
    padded = _pad_axis(padded1, axis2, bc2, dy)

    # Cached 4-corner slices: f[i±1, j±1] on the doubly-padded array
    ndim = data.ndim
    pp = padded[_cached_corner(ndim, axis1, _SLC_HI, axis2, _SLC_HI)]
    pm = padded[_cached_corner(ndim, axis1, _SLC_HI, axis2, _SLC_LO)]
    mp = padded[_cached_corner(ndim, axis1, _SLC_LO, axis2, _SLC_HI)]
    mm = padded[_cached_corner(ndim, axis1, _SLC_LO, axis2, _SLC_LO)]

    # In-place: (pp - pm - mp + mm) / (4*dx*dy)
    result = np.subtract(pp, pm)
    result -= mp
    result += mm
    result *= 1.0 / (4.0 * dx * dy)
    return result


def biharmonic(
    data: np.ndarray,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Biharmonic operator ∇⁴f = ∇²(∇²f)."""
    return laplacian(laplacian(data, grid, bc), grid, bc)


def identity(
    data: np.ndarray,
    grid: GridInfo,  # noqa: ARG001
    bc: BCSpec | None = None,  # noqa: ARG001
) -> np.ndarray:
    """Identity operator — returns data unchanged (no copy)."""
    return data


# ---------------------------------------------------------------------------
# Operator registry
# ---------------------------------------------------------------------------

# Each entry maps operator name → callable(data, grid, bc) → ndarray.
# Partials bind axis arguments for directional operators.


def _make_directional_laplacian(ax: int):  # noqa: ANN202
    def _op(data: np.ndarray, grid: GridInfo, bc: BCSpec | None = None) -> np.ndarray:
        # Bypass directional_laplacian wrapper — normalize BCs once, call raw
        bcs = _normalize_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
        return _directional_laplacian_raw(data, ax, grid, bcs[ax])

    return _op


def _make_gradient(ax: int):  # noqa: ANN202
    def _op(data: np.ndarray, grid: GridInfo, bc: BCSpec | None = None) -> np.ndarray:
        return gradient(data, ax, grid, bc)

    return _op


def _make_cross_derivative(ax1: int, ax2: int):  # noqa: ANN202
    def _op(data: np.ndarray, grid: GridInfo, bc: BCSpec | None = None) -> np.ndarray:
        return cross_derivative(data, ax1, ax2, grid, bc)

    return _op


OPERATOR_REGISTRY: dict[str, Any] = {
    "identity": identity,
    "laplacian": laplacian,
    "laplacian_x": _make_directional_laplacian(0),
    "laplacian_y": _make_directional_laplacian(1),
    "laplacian_z": _make_directional_laplacian(2),
    "gradient_x": _make_gradient(0),
    "gradient_y": _make_gradient(1),
    "gradient_z": _make_gradient(2),
    "cross_derivative_xy": _make_cross_derivative(0, 1),
    "cross_derivative_xz": _make_cross_derivative(0, 2),
    "cross_derivative_yz": _make_cross_derivative(1, 2),
    "biharmonic": biharmonic,
}


# Minimum spatial dimension required by each operator
_OPERATOR_MIN_DIM: dict[str, int] = {
    "identity": 1,
    "laplacian": 1,
    "laplacian_x": 1,
    "laplacian_y": 2,
    "laplacian_z": 3,
    "gradient_x": 1,
    "gradient_y": 2,
    "gradient_z": 3,
    "cross_derivative_xy": 2,
    "cross_derivative_xz": 3,
    "cross_derivative_yz": 3,
    "biharmonic": 1,
    "first_derivative_t": 1,
}


def operator_min_dim(name: str) -> int:
    """Return the minimum spatial dimension required by *name*.

    Raises
    ------
    ValueError
        If *name* is not a recognized operator.
    """
    dim = _OPERATOR_MIN_DIM.get(name)
    if dim is not None:
        return dim
    msg = f"Unknown operator {name!r}; known: {sorted(_OPERATOR_MIN_DIM)}"
    raise ValueError(msg)


def apply_operator(
    name: str,
    data: np.ndarray,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Look up *name* in the registry and apply it to *data*.

    Raises
    ------
    ValueError
        If *name* is not a known operator.
    """
    fn = OPERATOR_REGISTRY.get(name)
    if fn is None:
        msg = f"Unknown operator {name!r}; known: {sorted(OPERATOR_REGISTRY)}"
        raise ValueError(msg)
    result: np.ndarray = fn(data, grid, bc)
    return result
