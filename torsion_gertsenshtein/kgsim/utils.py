from __future__ import annotations

from typing import TYPE_CHECKING, cast

from pde import ScalarField

if TYPE_CHECKING:
    # type-checker only import (stubs may be missing at runtime)
    from collections.abc import Mapping, Sequence

    from pde.fields.datafield_base import DataFieldBase  # type: ignore[import]


def natural_center(bounds: Sequence[tuple[float, float]]) -> list[float]:
    """
    Return the midpoint (natural center) for each interval in `bounds`.

    Parameters.
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
    Infer a boundary condition argument suitable for py-pde functions.

    Returns
    -------
      - None: let py-pde use the grid's own (e.g. periodic) boundary handling.
      - dict: explicit boundary condition (e.g. homogeneous Neumann).
    """
    # If the grid exposes a boundaries mapping (already explicit), reuse it.
    bd = getattr(grid, "boundaries", None)
    if bd is not None:
        return bd  # assume it's already the right mapping at runtime

    # If the grid is periodic on any axis, prefer letting py-pde handle BCs implicitly.
    periodic = getattr(grid, "periodic", None)
    if periodic is None:
        return None

    if isinstance(periodic, bool):
        # fully periodic -> let py-pde handle it
        if periodic:
            return None
        # fully non-periodic -> explicit homogeneous Neumann
        return {"all": {"value": "natural"}}

    # periodic is a Sequence[bool]: if any axis is periodic, let py-pde determine BCs
    if any(bool(p) for p in periodic):
        return None

    # no periodic axes -> explicit homogeneous Neumann
    return {"all": {"value": "natural"}}
