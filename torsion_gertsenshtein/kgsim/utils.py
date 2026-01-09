from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from pde import ScalarField

if TYPE_CHECKING:
    # type-checker only import (stubs may be missing at runtime)
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
        result = {"all": {"type": "derivative", "value": 0.0}}
    elif isinstance(periodic, Sequence):
        periodic_seq = cast("Sequence[bool]", periodic)
        if any(periodic_seq):
            # At least one periodic direction - use auto BC
            result = "auto_periodic_neumann"
        else:
            # No periodic directions
            result = {"all": {"type": "derivative", "value": 0.0}}
    else:
        result = "auto_periodic_neumann"

    return result
