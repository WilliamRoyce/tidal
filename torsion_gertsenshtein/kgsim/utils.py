from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from pde import CartesianGrid, ScalarField

if TYPE_CHECKING:
    # type-checker only import (stubs may be missing at runtime)
    import numpy as np
    from pde.fields.datafield_base import DataFieldBase  # type: ignore[import]

# Type alias for boundary condition data compatible with py-pde
# Can be: string for auto BC, dict for explicit BC mapping, or None
BCDescriptor = str | dict[str, dict[str, str | float]] | None


def natural_center(bounds: Sequence[tuple[float, float]]) -> list[float]:
    """
    Return the midpoint (natural center) for each interval in `bounds`.

    Parameters
    ----------
    bounds : Sequence[tuple[float, float]]
        An iterable of 2-tuples (or 2-length sequences) representing numeric intervals
        (a, b). Each tuple's elements should be numbers (ints or floats). The function
        computes the midpoint for each interval as (a + b) / 2 and does not require
        that a <= b.

    Returns
    -------
    list[float]
        A list of midpoints corresponding to each input interval, in the same order.
    """
    return [(a + b) / 2 for (a, b) in bounds]


def extract_grid_coordinates(
    grid: CartesianGrid, *, flatten: bool = True
) -> np.ndarray:
    """Extract and normalize grid cell coordinates from CartesianGrid.

    Handles different py-pde versions that may return cell coordinates in
    different formats (grid-shaped vs flattened).

    Parameters
    ----------
    grid : CartesianGrid
        The computational grid.
    flatten : bool, optional
        If True, return coordinates in shape (N_cells, dim).
        If False, return coordinates in shape (*grid.shape, dim).
        Default is True.

    Returns
    -------
    np.ndarray
        Normalized coordinate array:
        - If flatten=True: shape (N_cells, dim)
        - If flatten=False: shape (*grid.shape, dim)

    Raises
    ------
    ValueError
        If cell_coords has an unexpected shape that cannot be normalized.

    Notes
    -----
    Different py-pde versions may provide cell_coords in different formats:
    - Grid-shaped: (*grid.shape, dim) - newer versions
    - Flattened: (N_cells, dim) - older versions
    This function provides a unified interface regardless of format.
    """
    coordinates = cast("np.ndarray", grid.cell_coords)

    # Handle grid-shaped format: (*grid.shape, dim)
    if coordinates.ndim == grid.dim + 1 and coordinates.shape[-1] == grid.dim:
        if flatten:
            return coordinates.reshape(-1, grid.dim)
        return coordinates

    # Handle flattened format: (N_cells, dim)
    if coordinates.ndim == 2 and coordinates.shape[1] == grid.dim:  # noqa: PLR2004
        if not flatten:
            return coordinates.reshape(*grid.shape, grid.dim)
        return coordinates

    msg = f"Unexpected cell_coords shape: {coordinates.shape} for grid with dim={grid.dim}"
    raise ValueError(msg)


def sub_scalar_fields(
    a: ScalarField | DataFieldBase, b: ScalarField | DataFieldBase
) -> ScalarField:
    """
    Subtract two scalar fields and return a ScalarField.

    Raises
    ------
    TypeError
        If either argument is not a ScalarField at runtime.
    """
    if not isinstance(a, ScalarField) or not isinstance(b, ScalarField):
        msg = f"Cannot subtract {type(a)} and {type(b)}"
        raise TypeError(msg)
    # result has runtime type ScalarField, cast so the type checker knows it
    result = a - b
    return cast("ScalarField", result)


def mul_scalar_field(
    scalar: complex, field: ScalarField | DataFieldBase
) -> ScalarField:
    """Multiply a numeric scalar and a ScalarField, returning a ScalarField (runtime-checked).

    Raises
    ------
    TypeError
        If the field argument is not a ScalarField at runtime.
    """
    if not isinstance(field, ScalarField):
        msg = f"Cannot multiply {type(field)} with {type(scalar)}"
        raise TypeError(msg)
    # rely on field.__rmul__ / __mul__ at runtime; cast so the type checker knows the result type
    return cast("ScalarField", scalar * field)


def infer_bc_from_grid(
    grid: object,
) -> BCDescriptor:
    """
    Infer boundary-condition descriptor from a grid-like object.

    For periodic grids, returns 'auto_periodic_neumann' which is required for
    gradient chaining to work correctly in py-pde.
    For non-periodic grids, returns an explicit BC mapping with Neumann (derivative=0).

    Returns
    -------
    BCDescriptor
        Boundary condition descriptor compatible with py-pde's BoundariesData type.
        Can be: str (e.g., 'auto_periodic_neumann'), dict[str, dict[str, str | float]], or None.

    Notes
    -----
    When using gradient().gradient() for second derivatives, py-pde requires
    an explicit boundary condition string like 'auto_periodic_neumann' rather
    than None, even for periodic grids.
    """
    # Try explicit boundaries first
    bd = getattr(grid, "boundaries", None)
    periodic = getattr(grid, "periodic", None)

    result: BCDescriptor
    # explicit boundaries override periodic inference
    if bd is not None:
        result = bd  # type: ignore[assignment]
    # handle periodic being None/True/False/sequence/other in branches
    elif periodic is None or periodic is True:
        # Use explicit BC string for periodic grids (required for gradient chaining)
        result = "auto_periodic_neumann"
    elif periodic is False:
        # Non-periodic: use Neumann (derivative=0) boundary conditions
        result = "derivative"
    elif isinstance(periodic, Sequence):
        periodic_seq = cast("Sequence[bool]", periodic)
        # At least one periodic direction - use auto BC, otherwise derivative
        result = "auto_periodic_neumann" if any(periodic_seq) else "derivative"
    else:
        result = "auto_periodic_neumann"

    return result
