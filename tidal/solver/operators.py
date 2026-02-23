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
    """Neumann normal derivative D (outward convention)."""
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


def _str_to_axis_bc(bc_str: str) -> AxisBCSpec:
    """Convert a legacy BC string to an ``AxisBCSpec``.

    ``"periodic"`` maps to a periodic axis.  ``"neumann"``/``"dirichlet"``
    map to symmetric homogeneous BC on both sides.  ``"robin"`` maps to
    zero Robin on both sides (gamma=0, beta=0, equivalent to Neumann).

    Raises
    ------
    ValueError
        If *bc_str* is not a recognized BC type.
    """
    if bc_str == "periodic":
        return AxisBCSpec(periodic=True)
    if bc_str not in _VALID_BC:
        msg = f"Unknown BC type: {bc_str!r}; valid: {sorted(_VALID_BC)}"
        raise ValueError(msg)
    side = SideBCSpec(kind=bc_str)
    return AxisBCSpec(periodic=False, low=side, high=side)


# Legacy type alias — widened to accept structured specs
BCSpec = str | tuple[str, ...] | tuple[AxisBCSpec, ...]
"""Single BC string, per-axis strings, or per-axis ``AxisBCSpec`` objects."""


def _normalize_bc(
    bc: BCSpec, grid: GridInfo
) -> tuple[str | AxisBCSpec, ...]:
    """Expand a BC spec to a per-axis tuple, validating each entry.

    Returns a tuple of either strings or ``AxisBCSpec`` objects (never mixed
    within one call — the type depends on what was passed in).

    Raises
    ------
    ValueError
        If the number of BC entries doesn't match ``grid.ndim`` or an
        unknown BC type is encountered.
    TypeError
        If a BC entry is neither a string nor an ``AxisBCSpec``.
    """
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
        elif not isinstance(b, AxisBCSpec):
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
# Ghost-cell padding helpers
# ---------------------------------------------------------------------------


def _pad_axis(
    data: np.ndarray,
    axis: int,
    bc: str | AxisBCSpec,
    dx: float = 1.0,
) -> np.ndarray:
    """Pad *data* with one ghost cell on each side along *axis*.

    Uses the unified ghost-cell formula: ``ghost = const + factor * interior``
    where ``(const, factor)`` are determined by the BC type.

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

    n = data.shape[axis]

    if bc.periodic:
        left = np.take(data, [n - 1], axis=axis)
        right = np.take(data, [0], axis=axis)
    else:
        assert bc.low is not None  # guaranteed by AxisBCSpec.__post_init__
        assert bc.high is not None
        c_lo, f_lo = bc.low.ghost_params(dx)
        c_hi, f_hi = bc.high.ghost_params(dx)

        interior_lo = np.take(data, [0], axis=axis)
        interior_hi = np.take(data, [n - 1], axis=axis)

        # Optimised path: skip addition when const == 0 (homogeneous BCs)
        left = f_lo * interior_lo if c_lo == 0.0 else c_lo + f_lo * interior_lo
        right = f_hi * interior_hi if c_hi == 0.0 else c_hi + f_hi * interior_hi

    return np.concatenate([left, data, right], axis=axis)


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

    # Fast path for periodic: np.roll avoids ghost-cell allocation
    axis_bc = _resolve_axis_bc(bc_axis)
    if axis_bc.periodic:
        return (np.roll(data, -1, axis=axis) - np.roll(data, 1, axis=axis)) / (2.0 * dx)

    # Ghost-cell path for non-periodic BCs
    padded = _pad_axis(data, axis, axis_bc, dx)
    slc_left = [slice(None)] * data.ndim
    slc_right = [slice(None)] * data.ndim
    slc_left[axis] = slice(0, -2)  # f[i-1]
    slc_right[axis] = slice(2, None)  # f[i+1]
    return (padded[tuple(slc_right)] - padded[tuple(slc_left)]) / (2.0 * dx)


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
    dx = grid.dx[axis]
    inv_dx2 = 1.0 / (dx * dx)
    bc_axis = bcs[axis]

    axis_bc = _resolve_axis_bc(bc_axis)
    if axis_bc.periodic:
        return (
            np.roll(data, -1, axis=axis) - 2.0 * data + np.roll(data, 1, axis=axis)
        ) * inv_dx2

    padded = _pad_axis(data, axis, axis_bc, dx)
    slc_left = [slice(None)] * data.ndim
    slc_center = [slice(None)] * data.ndim
    slc_right = [slice(None)] * data.ndim
    slc_left[axis] = slice(0, -2)
    slc_center[axis] = slice(1, -1)
    slc_right[axis] = slice(2, None)
    return (
        padded[tuple(slc_right)]
        - 2.0 * padded[tuple(slc_center)]
        + padded[tuple(slc_left)]
    ) * inv_dx2


def laplacian(
    data: np.ndarray,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Full Laplacian ∇²f = Σ_i ∂²f/∂x_i²."""
    result = np.zeros_like(data)
    for ax in range(grid.ndim):
        result += directional_laplacian(data, ax, grid, bc)
    return result


def cross_derivative(
    data: np.ndarray,
    axis1: int,
    axis2: int,
    grid: GridInfo,
    bc: BCSpec | None = None,
) -> np.ndarray:
    """Mixed second derivative ∂²f/(∂x_i ∂x_j).

    Applied as gradient along axis2 first, then gradient along axis1.
    """
    d1 = gradient(data, axis2, grid, bc)
    return gradient(d1, axis1, grid, bc)


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
        return directional_laplacian(data, ax, grid, bc)

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
