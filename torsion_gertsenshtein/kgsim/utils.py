from __future__ import annotations

from typing import TYPE_CHECKING, cast

from pde import ScalarField

if TYPE_CHECKING:
    # type-checker only import (stubs may be missing at runtime)
    from pde.fields.datafield_base import DataFieldBase  # type: ignore[import]


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
