"""Scalar field profile generators for spatial coefficient fields.

This module provides utilities for creating spatially-varying mass and potential
fields used with InhomogeneousKGPDE.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import cast

import numpy as np
from pde import CartesianGrid, ScalarField

ArrayLikeFunc = Callable[[np.ndarray], np.ndarray]
# coords is (N, dim) array of cell centers


def constant_field(grid: CartesianGrid, value: float) -> ScalarField:
    """Create a constant scalar field on the given grid.

    Parameters
    ----------
    grid : CartesianGrid
        The computational grid on which to define the field.
    value : float
        The constant value to assign to all grid points.

    Returns
    -------
    ScalarField
        A scalar field with uniform value across the entire grid.

    Examples
    --------
    >>> from pde import CartesianGrid
    >>> grid = CartesianGrid([[0, 10]], 100)
    >>> field = constant_field(grid, 1.5)
    >>> field.data.min() == field.data.max() == 1.5
    True
    """
    return ScalarField(grid, data=float(value))


def step_region_1d(
    grid: CartesianGrid,
    x0: float,
    x1: float,
    inside_value: float,
    outside_value: float,
) -> ScalarField:
    """Create a 1D scalar field that is a step function between x0 and x1.

    Parameters
    ----------
    grid : CartesianGrid
        The computational grid (must be 1D).
    x0 : float
        Left boundary (inclusive) of the inside region.
    x1 : float
        Right boundary (inclusive) of the inside region.
    inside_value : float
        Value assigned to cells whose centers lie in [x0, x1].
    outside_value : float
        Value assigned to cells outside [x0, x1].

    Returns
    -------
    ScalarField
        A scalar field defined on `grid` with the step region applied.

    Raises
    ------
    ValueError
        If `grid` is not one-dimensional.
    """
    if grid.dim != 1:
        msg = "step_region_1d requires a 1D grid"
        raise ValueError(msg)
    coordinates = cast("np.ndarray", grid.cell_coords)  # (N, dim)
    x = coordinates[:, 0].astype(float)
    mask = (x >= float(x0)) & (x <= float(x1))
    data = np.where(mask, inside_value, outside_value)
    return ScalarField(grid, data=data.reshape(grid.shape))


def from_callable(grid: CartesianGrid, fn: ArrayLikeFunc) -> ScalarField:
    """Create a scalar field by evaluating a callable at grid cell centers.

    This is a flexible factory for creating arbitrary field profiles from
    user-defined functions that operate on coordinate arrays.

    Parameters
    ----------
    grid : CartesianGrid
        The computational grid on which to define the field.
    fn : Callable[[np.ndarray], np.ndarray]
        A function that takes cell center coordinates as input and returns
        field values. The input array has shape (N, dim) where N is the total
        number of grid cells and dim is the spatial dimension. The function
        should return a 1D array of shape (N,) with the field value at each
        cell center.

    Returns
    -------
    ScalarField
        A scalar field with values computed from the callable, reshaped to
        match the grid structure.

    Examples
    --------
    >>> from pde import CartesianGrid
    >>> import numpy as np
    >>> grid = CartesianGrid([[0, 10]], 100)
    >>> # Create a linear gradient field
    >>> field = from_callable(grid, lambda coords: coords[:, 0])
    >>> # Create a radial field in 2D
    >>> grid_2d = CartesianGrid([[0, 1], [0, 1]], [10, 10])
    >>> radial = from_callable(grid_2d, lambda c: np.sqrt(c[:, 0]**2 + c[:, 1]**2))
    """
    coordinates = cast("np.ndarray", grid.cell_coords)  # (N, dim)
    vals = fn(coordinates)  # shape (N,)
    return ScalarField(grid, data=np.asarray(vals).reshape(grid.shape))
