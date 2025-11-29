from __future__ import annotations

from collections.abc import Callable
from typing import cast

import numpy as np
from pde import CartesianGrid, ScalarField

ArrayLikeFunc = Callable[[np.ndarray], np.ndarray]
# coords is (N, dim) array of cell centers


def constant_field(grid: CartesianGrid, value: float) -> ScalarField:
    """Create a constant scalar field on the given grid."""
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
    """fn(coords) -> values; coords has shape (N, dim)."""
    coordinates = cast("np.ndarray", grid.cell_coords)  # (N, dim)
    vals = fn(coordinates)  # shape (N,)
    return ScalarField(grid, data=np.asarray(vals).reshape(grid.shape))
