from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, cast

from pde import ScalarField

if TYPE_CHECKING:
    # type-checker only import (stubs may be missing at runtime)
    from collections.abc import Mapping

    from pde.fields.datafield_base import DataFieldBase  # type: ignore[import]


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


def infer_bc_from_grid(grid: object) -> Mapping[str, Mapping[str, object]] | None:
    """
    Infer boundary-condition descriptor from a grid-like object.

    Returns None to indicate "let py-pde use the grid's (periodic) defaults".
    Only returns an explicit BC mapping for confirmed non-periodic grids.
    """
    # Try explicit boundaries first
    bd = getattr(grid, "boundaries", None)
    periodic = getattr(grid, "periodic", None)

    result: Mapping[str, Mapping[str, object]] | None
    # explicit boundaries override periodic inference
    if bd is not None:
        result = bd
    # handle periodic being None/True/False/sequence/other in branches
    elif periodic is None or periodic is True:
        result = None
    elif periodic is False:
        result = {"all": {"type": "derivative", "value": 0.0}}
    elif isinstance(periodic, Sequence):
        periodic_seq = cast("Sequence[bool]", periodic)
        if any(periodic_seq):
            result = None
        else:
            result = {"all": {"type": "derivative", "value": 0.0}}
    else:
        result = None

    return result
