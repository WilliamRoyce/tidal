"""Spatial finite-difference operators on plain numpy arrays.

All operators accept an ``np.ndarray`` (field data) and a ``GridInfo``
descriptor, returning an ``np.ndarray`` of the same shape.  No py-pde
ScalarField / VectorField wrappers — pure numpy with explicit ghost-cell
boundary handling.

Boundary condition types
------------------------
- ``"periodic"`` — wraparound via ``np.roll``
- ``"neumann"`` — zero normal derivative (mirror ghost cell)
- ``"dirichlet"`` — zero value (antisymmetric ghost cell)

The ``bc`` parameter may be:
- A single string applied to all axes
- A tuple of strings, one per axis
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo

# ---------------------------------------------------------------------------
# BC types
# ---------------------------------------------------------------------------

BCSpec = str | tuple[str, ...]
"""Single BC string (all axes) or per-axis tuple."""

_VALID_BC = frozenset({"periodic", "neumann", "dirichlet"})


def _normalise_bc(bc: BCSpec, grid: GridInfo) -> tuple[str, ...]:
    """Expand a BC spec to a per-axis tuple, validating each entry.

    Raises
    ------
    ValueError
        If the number of BC entries doesn't match ``grid.ndim`` or an
        unknown BC type is encountered.
    """
    bcs = (bc,) * grid.ndim if isinstance(bc, str) else tuple(bc)
    if len(bcs) != grid.ndim:
        msg = f"Expected {grid.ndim} BC entries, got {len(bcs)}"
        raise ValueError(msg)
    for i, b in enumerate(bcs):
        if b not in _VALID_BC:
            msg = f"Unknown BC type {b!r} for axis {i}; valid: {sorted(_VALID_BC)}"
            raise ValueError(msg)
    return bcs


def _bc_from_grid(grid: GridInfo) -> tuple[str, ...]:
    """Infer BCs from grid periodicity (periodic or neumann)."""
    return tuple("periodic" if p else "neumann" for p in grid.periodic)


# ---------------------------------------------------------------------------
# Ghost-cell padding helpers
# ---------------------------------------------------------------------------


def _pad_axis(data: np.ndarray, axis: int, bc: str) -> np.ndarray:
    """Pad *data* with one ghost cell on each side along *axis*.

    For periodic axes, wraps around.  For neumann, mirrors the boundary
    cell.  For dirichlet, negates the boundary cell (antisymmetric
    reflection about the zero-value boundary).

    Raises
    ------
    ValueError
        If *bc* is not one of ``periodic``, ``neumann``, ``dirichlet``.
    """
    n = data.shape[axis]
    if bc == "periodic":
        # Wrap: left ghost = last cell, right ghost = first cell
        left = np.take(data, [n - 1], axis=axis)
        right = np.take(data, [0], axis=axis)
    elif bc == "neumann":
        # Mirror: ∂f/∂n = 0 ⟹ ghost = boundary cell
        left = np.take(data, [0], axis=axis)
        right = np.take(data, [n - 1], axis=axis)
    elif bc == "dirichlet":
        # Antisymmetric: f(boundary) = 0 ⟹ ghost = -boundary cell
        left = -np.take(data, [0], axis=axis)
        right = -np.take(data, [n - 1], axis=axis)
    else:
        msg = f"Unknown BC type: {bc!r}"
        raise ValueError(msg)
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
    bcs = _normalise_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    dx = grid.dx[axis]
    bc_axis = bcs[axis]

    if bc_axis == "periodic":
        return (
            np.roll(data, -1, axis=axis) - np.roll(data, 1, axis=axis)
        ) / (2.0 * dx)

    # Ghost-cell path for non-periodic BCs
    padded = _pad_axis(data, axis, bc_axis)
    # Slice out interior after central difference
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
    bcs = _normalise_bc(bc, grid) if bc is not None else _bc_from_grid(grid)
    dx = grid.dx[axis]
    inv_dx2 = 1.0 / (dx * dx)
    bc_axis = bcs[axis]

    if bc_axis == "periodic":
        return (
            np.roll(data, -1, axis=axis)
            - 2.0 * data
            + np.roll(data, 1, axis=axis)
        ) * inv_dx2

    padded = _pad_axis(data, axis, bc_axis)
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
