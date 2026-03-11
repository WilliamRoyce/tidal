"""Spatial finite-difference operators on plain numpy arrays.

All operators accept an ``np.ndarray`` (field data) and a ``GridInfo``
descriptor, returning an ``np.ndarray`` of the same shape.  No py-pde
ScalarField / VectorField wrappers — pure numpy with explicit ghost-cell
boundary handling.

Finite-difference accuracy orders
----------------------------------
Supports 2nd, 4th, and 6th-order central-difference stencils, selectable
at runtime via :func:`set_fd_order`.  Higher orders use wider stencils
(more ghost cells) for dramatically improved spatial accuracy:

- **Order 2** (default): 3-point stencil, 1 ghost cell/side.
- **Order 4**: 5-point stencil, 2 ghost cells/side.
- **Order 6**: 7-point stencil, 3 ghost cells/side.

Stencil coefficients from Fornberg (1988), "Generation of Finite Difference
Formulas on Arbitrarily Spaced Grids", Mathematics of Computation 51(184),
pp. 699-706, Table 1 (central differences).

Boundary condition types
------------------------
- ``"periodic"`` — wraparound via ``np.roll``
- ``"neumann"``  — zero normal derivative (mirror ghost cell)
- ``"dirichlet"``— zero value (antisymmetric ghost cell)
- ``"robin"``    -- mixed: d_n f + gamma*f = beta

All non-periodic BCs are unified via the ghost-cell formula
``ghost = const + factor * interior_cell`` where ``(const, factor)``
are determined by BC type, value, and grid spacing.  See ``SideBCSpec``
for the derivation.  For higher-order stencils (order > 2), ghost cells
are filled layer-by-layer from interior outward (recursive ghost fill;
LeVeque, "Finite Difference Methods for ODEs and PDEs", SIAM, 2007, ch. 2).

The ``bc`` parameter may be:
- A single string applied to all axes (e.g. ``"periodic"``)
- A tuple of strings, one per axis (e.g. ``("neumann", "periodic")``)
- A tuple of ``AxisBCSpec`` objects for per-side / non-zero / Robin BCs

References
----------
- Fornberg (1988), Math. Comp. 51(184), pp. 699-706 — stencil coefficients.
- LeVeque (2007), "Finite Difference Methods for ODEs and PDEs", SIAM — convergence theory.
- Zwicker (2020), J. Open Source Software 5(48) — ghost-cell virtual-point approach.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo

# ---------------------------------------------------------------------------
# Finite-difference order (module-level state)
# ---------------------------------------------------------------------------
# Set once at simulation startup via set_fd_order().  Determines stencil
# width and ghost-cell count for all spatial operators.
#
# Reference: Fornberg, "Generation of Finite Difference Formulas on
# Arbitrarily Spaced Grids", Mathematics of Computation 51(184), 1988.
# Stencil coefficients from Table 1 (central differences).

_fd_order: int = 2
"""Current FD accuracy order: 2 (default), 4, or 6."""

_n_ghost: int = 1
"""Number of ghost cells per side: fd_order // 2."""


def set_fd_order(order: int) -> None:
    """Set the finite-difference accuracy order for all spatial operators.

    Must be called before any operator evaluation (typically at simulation
    startup).  Clears the ghost-cell padding cache since buffer shapes change.

    Parameters
    ----------
    order : int
        Accuracy order: 2 (3-point stencil), 4 (5-point), or 6 (7-point).

    Raises
    ------
    ValueError
        If *order* is not 2, 4, or 6.
    """
    global _fd_order, _n_ghost  # noqa: PLW0603
    if order not in {2, 4, 6}:
        msg = f"FD order must be 2, 4, or 6; got {order}"
        raise ValueError(msg)
    _fd_order = order
    _n_ghost = order // 2
    # Invalidate cached padding buffers — shapes depend on ghost count
    _pad_cache.clear()


def get_fd_order() -> int:
    """Return the current FD accuracy order."""
    return _fd_order


def get_n_ghost() -> int:
    """Return the current number of ghost cells per side."""
    return _n_ghost


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


def _resolve_axis_bc(bc_entry: str | AxisBCSpec) -> AxisBCSpec:  # pyright: ignore[reportUnusedFunction]  # used by rhs.py
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
    """Pre-allocated buffer and cached slices for one (shape, axis, n_ghost) combo.

    ALL slice tuples are pre-computed at construction time so the hot-path
    ``_pad_axis()`` does zero per-call tuple/slice construction.

    Supports variable ghost-cell width for higher-order FD stencils:
    - Order 2: 1 ghost cell per side (3-point stencil)
    - Order 4: 2 ghost cells per side (5-point stencil)
    - Order 6: 3 ghost cells per side (7-point stencil)
    """

    __slots__ = (
        "buf",
        "interior_slc",
        "left_slcs",
        "ng",
        "periodic_left_src",
        "periodic_right_src",
        "recursive_left_src",
        "recursive_right_src",
        "right_slcs",
        "src_hi_slcs",
        "src_lo_slcs",
    )

    def __init__(self, data_shape: tuple[int, ...], axis: int, ng: int) -> None:
        n = data_shape[axis]
        ndim = len(data_shape)
        self.ng = ng
        padded_shape = list(data_shape)
        padded_shape[axis] += 2 * ng
        self.buf = np.empty(tuple(padded_shape))
        # Interior: data sits at [ng .. n+ng)
        self.interior_slc = _slice_tuple(ndim, axis, slice(ng, n + ng))

        # Per-layer ghost destination slices in padded buffer
        self.left_slcs: list[tuple[slice, ...]] = []
        self.right_slcs: list[tuple[slice, ...]] = []
        # Per-layer source slices from original data (for non-periodic k=0)
        self.src_lo_slcs: list[tuple[slice, ...]] = []
        self.src_hi_slcs: list[tuple[slice, ...]] = []
        # Per-layer source slices from data for periodic wrapping
        self.periodic_left_src: list[tuple[slice, ...]] = []
        self.periodic_right_src: list[tuple[slice, ...]] = []
        # Per-layer source slices from padded buffer for recursive ghost fill (k>0)
        self.recursive_left_src: list[tuple[slice, ...]] = []
        self.recursive_right_src: list[tuple[slice, ...]] = []

        for k in range(ng):
            # Ghost cell destination: left_slcs[k] = layer k (0=outermost)
            self.left_slcs.append(_slice_tuple(ndim, axis, slice(k, k + 1)))
            self.right_slcs.append(
                _slice_tuple(ndim, axis, slice(n + 2 * ng - 1 - k, n + 2 * ng - k))
            )
            # Non-periodic source from data boundary (innermost=index 0)
            self.src_lo_slcs.append(_slice_tuple(ndim, axis, slice(ng - 1 - k, ng - k)))
            self.src_hi_slcs.append(
                _slice_tuple(ndim, axis, slice(n + k, n + k + 1))
            )
            # Periodic source from data: left ghost copies from right end
            self.periodic_left_src.append(
                _slice_tuple(ndim, axis, slice(n - ng + k, n - ng + k + 1))
            )
            # Periodic source from data: right ghost copies from left end
            self.periodic_right_src.append(
                _slice_tuple(ndim, axis, slice(ng - 1 - k, ng - k))
            )
            # Recursive ghost fill source (from padded buffer, for k>0)
            if k > 0:
                # Left: source from previously-filled ghost at index (ng - k)
                src_l = ng - k
                self.recursive_left_src.append(
                    _slice_tuple(ndim, axis, slice(src_l, src_l + 1))
                )
                # Right: source from previously-filled ghost
                src_r = n + ng + k - 1
                self.recursive_right_src.append(
                    _slice_tuple(ndim, axis, slice(src_r, src_r + 1))
                )
            else:
                # Placeholders for k=0 (use data directly, not buffer)
                self.recursive_left_src.append(
                    _slice_tuple(ndim, axis, slice(0, 1))
                )
                self.recursive_right_src.append(
                    _slice_tuple(ndim, axis, slice(n - 1, n))
                )


class _PadBufferCache:
    """Lazily allocate and reuse padded buffers for ghost-cell operations.

    Keyed by ``(data.shape, axis, n_ghost)`` — since TIDAL operates on a
    fixed grid with a fixed FD order per run, each key is allocated only
    once.  Single-threaded assumption (no locking).
    """

    __slots__ = ("_cache",)

    def __init__(self) -> None:
        self._cache: dict[tuple[tuple[int, ...], int, int], _PadEntry] = {}

    def clear(self) -> None:
        """Invalidate all cached padding buffers."""
        self._cache.clear()

    def get(self, data_shape: tuple[int, ...], axis: int, ng: int | None = None) -> _PadEntry:
        if ng is None:
            ng = _n_ghost
        key = (data_shape, axis, ng)
        entry = self._cache.get(key)
        if entry is None:
            entry = _PadEntry(data_shape, axis, ng)
            self._cache[key] = entry
        return entry


_pad_cache = _PadBufferCache()


def _pad_axis(
    data: np.ndarray,
    axis: int,
    bc: str | AxisBCSpec,
    dx: float = 1.0,
    ng: int | None = None,
) -> np.ndarray:
    """Pad *data* with ghost cells on each side along *axis*.

    Uses pre-allocated buffers and cached slice tuples (via ``_pad_cache``)
    for zero-alloc ghost-cell writes in the hot path.

    The number of ghost cells is determined by the current FD order
    (``_n_ghost``), or overridden by *ng*.

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
    ng : int, optional
        Override ghost-cell count (default: ``_n_ghost`` from FD order).
    """
    if isinstance(bc, str):
        bc = _str_to_axis_bc(bc)

    if ng is None:
        ng = _n_ghost
    entry = _pad_cache.get(data.shape, axis, ng)
    buf = entry.buf

    # Copy interior into buffer
    buf[entry.interior_slc] = data

    if bc.periodic:
        # All source slices pre-cached in _PadEntry — zero per-call construction.
        for k in range(ng):
            buf[entry.left_slcs[k]] = data[entry.periodic_left_src[k]]
            buf[entry.right_slcs[k]] = data[entry.periodic_right_src[k]]
    else:
        assert bc.low is not None  # guaranteed by AxisBCSpec.__post_init__
        assert bc.high is not None
        c_lo, f_lo = bc.low.ghost_params(dx)
        c_hi, f_hi = bc.high.ghost_params(dx)

        # Fill ghost cells layer-by-layer from interior outward.
        # Layer 0 (innermost, adjacent to data boundary) uses data edge.
        # Layer k>0 uses the previously-filled ghost layer as source,
        # implementing the recursive ghost-cell fill standard in
        # computational physics (LeVeque, "Finite Difference Methods for
        # Ordinary and Partial Differential Equations", SIAM, 2007, §2.12).
        #
        # All slice tuples pre-cached in _PadEntry for zero per-call overhead.
        for k in range(ng):
            left_ghost_idx = ng - 1 - k  # innermost first
            left_ghost_slc = entry.left_slcs[left_ghost_idx]
            if k == 0:
                np.multiply(f_lo, data[entry.recursive_left_src[0]], out=buf[left_ghost_slc])
            else:
                np.multiply(f_lo, buf[entry.recursive_left_src[k]], out=buf[left_ghost_slc])
            if c_lo != 0.0:
                buf[left_ghost_slc] += c_lo

            right_ghost_idx = ng - 1 - k
            right_ghost_slc = entry.right_slcs[right_ghost_idx]
            if k == 0:
                np.multiply(f_hi, data[entry.recursive_right_src[0]], out=buf[right_ghost_slc])
            else:
                np.multiply(f_hi, buf[entry.recursive_right_src[k]], out=buf[right_ghost_slc])
            if c_hi != 0.0:
                buf[right_ghost_slc] += c_hi

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

    Stencil width adapts to the current FD order (set via ``set_fd_order``):

    - **Order 2** (3-point): ``(f[i+1] - f[i-1]) / (2·dx)``
    - **Order 4** (5-point): ``(f[i-2]/12 - 2f[i-1]/3 + 2f[i+1]/3 - f[i+2]/12) / dx``
    - **Order 6** (7-point): ``(-f[i-3]/60 + 3f[i-2]/20 - 3f[i-1]/4 + 3f[i+1]/4
      - 3f[i+2]/20 + f[i+3]/60) / dx``

    Reference: Fornberg (1988), Mathematics of Computation 51(184), Table 1.

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
    ng = _n_ghost
    padded = _pad_axis(data, axis, axis_bc, dx, ng=ng)
    ndim = data.ndim
    inv_dx = 1.0 / dx

    if _fd_order == 2:  # noqa: PLR2004
        slc_left = _cached_slice(ndim, axis, _SLC_LO)
        slc_right = _cached_slice(ndim, axis, _SLC_HI)
        result = np.subtract(padded[slc_right], padded[slc_left])
        result *= 0.5 * inv_dx
    elif _fd_order == 4:  # noqa: PLR2004
        result = _gradient_o4(padded, axis, ndim, ng, inv_dx)
    else:  # order 6
        result = _gradient_o6(padded, axis, ndim, ng, inv_dx)
    return result


def _gradient_o4(
    padded: np.ndarray, axis: int, ndim: int, ng: int, inv_dx: float,
) -> np.ndarray:
    """4th-order central gradient on a padded array (5-point stencil).

    Coefficients: {1/12, -2/3, 0, 2/3, -1/12} / dx

    Reference: Fornberg (1988), Table 1, d=1, N=2.
    """
    n = padded.shape[axis] - 2 * ng

    def _slc(offset: int) -> tuple[slice, ...]:
        start = ng + offset
        return _slice_tuple(ndim, axis, slice(start, start + n))

    # Fornberg (1988) Table 1, d=1, N=2 coefficients
    result = np.multiply(1.0 / 12.0, padded[_slc(-2)])
    result -= np.multiply(2.0 / 3.0, padded[_slc(-1)])
    result += np.multiply(2.0 / 3.0, padded[_slc(1)])
    result -= np.multiply(1.0 / 12.0, padded[_slc(2)])
    result *= inv_dx
    return result


def _gradient_o6(
    padded: np.ndarray, axis: int, ndim: int, ng: int, inv_dx: float,
) -> np.ndarray:
    """6th-order central gradient on a padded array (7-point stencil).

    Coefficients: {-1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60} / dx

    Reference: Fornberg (1988), Table 1, d=1, N=3.
    """
    n = padded.shape[axis] - 2 * ng

    def _slc(offset: int) -> tuple[slice, ...]:
        start = ng + offset
        return _slice_tuple(ndim, axis, slice(start, start + n))

    result = np.multiply(-1.0 / 60.0, padded[_slc(-3)])
    result += np.multiply(3.0 / 20.0, padded[_slc(-2)])
    result -= np.multiply(3.0 / 4.0, padded[_slc(-1)])
    result += np.multiply(3.0 / 4.0, padded[_slc(1)])
    result -= np.multiply(3.0 / 20.0, padded[_slc(2)])
    result += np.multiply(1.0 / 60.0, padded[_slc(3)])
    result *= inv_dx
    return result


def _directional_laplacian_raw(
    data: np.ndarray,
    axis: int,
    grid: GridInfo,
    bc_axis: str | AxisBCSpec,
) -> np.ndarray:
    """Central second derivative ∂²f/∂x_i² with pre-resolved per-axis BC.

    Stencil width adapts to the current FD order:

    - **Order 2** (3-point): ``(f[i+1] - 2f[i] + f[i-1]) / dx²``
    - **Order 4** (5-point): ``(-f[i-2]/12 + 4f[i-1]/3 - 5f[i]/2
      + 4f[i+1]/3 - f[i+2]/12) / dx²``
    - **Order 6** (7-point): ``(f[i-3]/90 - 3f[i-2]/20 + 3f[i-1]/2
      - 49f[i]/18 + 3f[i+1]/2 - 3f[i+2]/20 + f[i+3]/90) / dx²``

    Reference: Fornberg (1988), Mathematics of Computation 51(184), Table 1.

    Used by ``directional_laplacian`` and ``laplacian`` to avoid
    redundant BC normalization when computing multiple axes.
    """
    dx = grid.dx[axis]
    inv_dx2 = 1.0 / (dx * dx)

    # Inline _resolve_axis_bc to avoid function call overhead
    axis_bc = bc_axis if isinstance(bc_axis, AxisBCSpec) else _str_to_axis_bc(bc_axis)
    ng = _n_ghost
    padded = _pad_axis(data, axis, axis_bc, dx, ng=ng)
    ndim = data.ndim

    if _fd_order == 2:  # noqa: PLR2004
        slc_left = _cached_slice(ndim, axis, _SLC_LO)
        slc_center = _cached_slice(ndim, axis, _SLC_MID)
        slc_right = _cached_slice(ndim, axis, _SLC_HI)
        # In-place arithmetic: result = (right - 2*center + left) * inv_dx2
        result = np.subtract(padded[slc_right], padded[slc_center])
        result -= padded[slc_center]
        result += padded[slc_left]
        result *= inv_dx2
    elif _fd_order == 4:  # noqa: PLR2004
        result = _laplacian_raw_o4(padded, axis, ndim, ng, inv_dx2)
    else:  # order 6
        result = _laplacian_raw_o6(padded, axis, ndim, ng, inv_dx2)
    return result


def _laplacian_raw_o4(
    padded: np.ndarray, axis: int, ndim: int, ng: int, inv_dx2: float,
) -> np.ndarray:
    """4th-order central Laplacian on a padded array (5-point stencil).

    Coefficients: {-1/12, 4/3, -5/2, 4/3, -1/12} / dx²

    Reference: Fornberg (1988), Table 1, d=2, N=2.
    """
    n = padded.shape[axis] - 2 * ng

    def _slc(offset: int) -> tuple[slice, ...]:
        start = ng + offset
        return _slice_tuple(ndim, axis, slice(start, start + n))

    result = np.multiply(-1.0 / 12.0, padded[_slc(-2)])
    result += np.multiply(4.0 / 3.0, padded[_slc(-1)])
    result -= np.multiply(5.0 / 2.0, padded[_slc(0)])
    result += np.multiply(4.0 / 3.0, padded[_slc(1)])
    result -= np.multiply(1.0 / 12.0, padded[_slc(2)])
    result *= inv_dx2
    return result


def _laplacian_raw_o6(
    padded: np.ndarray, axis: int, ndim: int, ng: int, inv_dx2: float,
) -> np.ndarray:
    """6th-order central Laplacian on a padded array (7-point stencil).

    Coefficients: {1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90} / dx²

    Reference: Fornberg (1988), Table 1, d=2, N=3.
    """
    n = padded.shape[axis] - 2 * ng

    def _slc(offset: int) -> tuple[slice, ...]:
        start = ng + offset
        return _slice_tuple(ndim, axis, slice(start, start + n))

    result = np.multiply(1.0 / 90.0, padded[_slc(-3)])
    result -= np.multiply(3.0 / 20.0, padded[_slc(-2)])
    result += np.multiply(3.0 / 2.0, padded[_slc(-1)])
    result -= np.multiply(49.0 / 18.0, padded[_slc(0)])
    result += np.multiply(3.0 / 2.0, padded[_slc(1)])
    result -= np.multiply(3.0 / 20.0, padded[_slc(2)])
    result += np.multiply(1.0 / 90.0, padded[_slc(3)])
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

    Implemented as successive application of the 1D gradient operator along
    each axis: ``gradient(gradient(f, axis1), axis2)``.  This naturally
    inherits the current FD order from :func:`set_fd_order`:

    - **Order 2**: 4-point stencil (f[±1,±1])
    - **Order 4**: product of 5-point 1D stencils (5x5 = 25-point)
    - **Order 6**: product of 7-point 1D stencils (7x7 = 49-point)

    Reference: LeVeque, "Finite Difference Methods for Ordinary and Partial
    Differential Equations", SIAM, 2007, §1.4.
    """
    # Apply gradient along axis1 first, then along axis2
    # This is equivalent to the fused stencil but automatically inherits
    # the correct FD order from the gradient function
    tmp = gradient(data, axis1, grid, bc)
    return gradient(tmp, axis2, grid, bc)


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
