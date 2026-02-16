"""Utility functions."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, cast

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


def infer_bc_from_grid(
    grid: object,
) -> BCDescriptor:
    """Infer boundary-condition descriptor from a grid-like object.

    For periodic grids, returns 'auto_periodic_neumann' which is required for
    gradient chaining to work correctly in py-pde.
    For non-periodic grids, returns an explicit BC mapping with Neumann (derivative=0).

    Parameters
    ----------
    grid : object
        A grid-like object, typically a py-pde CartesianGrid. May have
        'boundaries' and/or 'periodic' attributes.

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
        result = bd
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
