"""Utility functions."""

from __future__ import annotations

from typing import Any

from pde import FieldCollection

# Type alias for boundary condition data compatible with py-pde
# Can be: string for auto BC, dict for explicit BC mapping, or None
BCDescriptor = str | dict[str, dict[str, str | float]] | None


def normalize_solve_result(
    result: FieldCollection | tuple[FieldCollection | None, Any] | None,
) -> FieldCollection:
    """Normalize solve results to a FieldCollection.

    Accepts:
    - FieldCollection
    - (FieldCollection | None, info)
    - None

    Raises
    ------
    TypeError
        If result is None or not a FieldCollection.
    """
    if isinstance(result, tuple):
        sol, _info = result
    else:
        sol = result

    if sol is None:
        msg = "Expected FieldCollection from solve"
        raise TypeError(msg)
    if not isinstance(sol, FieldCollection):  # type: ignore[arg-type]
        msg = "Expected FieldCollection from solve"
        raise TypeError(msg)
    return sol
