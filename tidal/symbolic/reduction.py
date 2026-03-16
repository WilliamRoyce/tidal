"""Plane-wave dimensional reduction for JSON equation specs.

After the Wolfram pipeline exports a JSON spec with the original higher-dimensional
metadata (dimension, coordinates, operator names), this module transforms it into a
clean reduced-dimension spec suitable for efficient simulation.

For example, a 3+1D theory reduced along z produces a 1+1D spec with:
- ``dimension: 2``, ``coordinates: ["t", "x"]``, ``signature: [-1, 1]``
- Operators remapped: ``laplacian_z → laplacian_x``, ``gradient_z → gradient_x``
- Killed-axis operators (``laplacian_x``, ``gradient_y``, etc.) removed
- ``coordinate_dependent`` arrays and ``coefficient_symbolic`` strings updated
- Provenance metadata recording the reduction

Curved-coordinate support:
- Coefficients depending only on the surviving coordinate are preserved and remapped
- Volume element is kept if it depends only on the surviving coordinate
- If any surviving coefficient or volume element references a killed coordinate,
  ``ValueError`` is raised (incompatible reduction)
"""

from __future__ import annotations

import copy
import re
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from collections.abc import Mapping

# ---------------------------------------------------------------------------
# Plane-wave dimensional reduction
# ---------------------------------------------------------------------------

# Coordinate reference pattern in Wolfram expressions: x[], y[], z[]
_COORD_REF_RE = re.compile(r"\b([xyzwvut])\s*\[\s*\]")

# Axis-specific operator patterns
_AXIS_OP_PREFIXES = ("laplacian_", "gradient_", "cross_derivative_")


def _build_operator_remap(
    propagation_axis: str,
    spatial_coords: list[str],
) -> dict[str, str | None]:
    """Build a mapping from original operator names to reduced operator names.

    The propagation axis operators are remapped to the first spatial axis ("x").
    Killed-axis operators map to ``None`` (remove the term).
    Dimension-agnostic operators (``identity``, ``laplacian``, ``first_derivative_t``,
    ``biharmonic``) are preserved.

    Parameters
    ----------
    propagation_axis : str
        The surviving spatial axis (e.g., "z").
    spatial_coords : list[str]
        All spatial coordinates (e.g., ["x", "y", "z"]).

    Returns
    -------
    dict[str, str | None]
        Mapping from original operator name to new name (or None to remove).
    """
    remap: dict[str, str | None] = {}

    for axis in spatial_coords:
        if axis == propagation_axis:
            remap[f"laplacian_{axis}"] = "laplacian_x"
            remap[f"gradient_{axis}"] = "gradient_x"
        else:
            remap[f"laplacian_{axis}"] = None
            remap[f"gradient_{axis}"] = None

    # All cross derivatives involve at least one killed axis → remove
    for i, a in enumerate(spatial_coords):
        for b in spatial_coords[i + 1 :]:
            remap[f"cross_derivative_{a}{b}"] = None

    # Dimension-agnostic operators preserved as-is
    remap.update(
        {
            op: op
            for op in (
                "identity",
                "laplacian",
                "first_derivative_t",
                "biharmonic",
                "time_derivative",
            )
        }
    )

    return remap


def _remap_coord_string(expr: str, propagation_axis: str) -> str:
    """Remap coordinate references in a Wolfram expression string.

    Replaces ``{prop_axis}[]`` with ``x[]`` if the propagation axis
    is not already "x".

    Parameters
    ----------
    expr : str
        Wolfram expression string (e.g., ``"2/z[]"``).
    propagation_axis : str
        The surviving spatial axis.

    Returns
    -------
    str
        Expression with coordinate references remapped.
    """
    if propagation_axis == "x":
        return expr
    return re.sub(
        rf"\b{re.escape(propagation_axis)}\s*\[\s*\]",
        "x[]",
        expr,
    )


def _check_killed_coord_refs(
    expr: str,
    killed: set[str],
    context: str,
) -> None:
    """Raise if expression references any killed coordinate.

    Parameters
    ----------
    expr : str
        Wolfram expression string to check.
    killed : set[str]
        Set of killed coordinate names.
    context : str
        Description of where this expression appears (for error messages).

    Raises
    ------
    ValueError
        If ``expr`` references any coordinate in ``killed``.
    """
    for match in _COORD_REF_RE.finditer(expr):
        coord = match.group(1)
        if coord in killed:
            msg = (
                f"Incompatible plane-wave reduction: {context} "
                f"references killed coordinate '{coord}[]' in expression: {expr}"
            )
            raise ValueError(msg)


def _remap_operator(op: str, remap: Mapping[str, str | None]) -> str | None:
    """Look up an operator in the remap table.

    Returns the new name, ``None`` (remove term), or the original name
    if not in the table (passthrough for unknown operators).
    """
    if op in remap:
        return remap[op]
    for prefix in _AXIS_OP_PREFIXES:
        if op.startswith(prefix):
            return None
    return op


def _remap_coord_deps(
    deps: list[str],
    propagation_axis: str,
    killed: set[str],
    context: str,
) -> list[str]:
    """Remap a ``coordinate_dependent`` array.

    Raises
    ------
    ValueError
        If any entry in ``deps`` is a killed coordinate.
    """
    for dep in deps:
        if dep in killed:
            msg = (
                f"Incompatible plane-wave reduction: {context} "
                f"has coordinate_dependent entry '{dep}' which is "
                f"a killed coordinate"
            )
            raise ValueError(msg)
    return [
        "x" if (dep == propagation_axis and propagation_axis != "x") else dep
        for dep in deps
    ]


def _transform_term(
    term: dict[str, Any],
    remap: dict[str, str | None],
    propagation_axis: str,
    killed: set[str],
) -> dict[str, Any] | None:
    """Transform a single RHS term, returning None if it should be removed."""
    op = term.get("operator", "")
    new_op = _remap_operator(op, remap)
    if new_op is None:
        return None

    result = copy.deepcopy(term)
    result["operator"] = new_op

    if "coefficient_symbolic" in result:
        sym = result["coefficient_symbolic"]
        _check_killed_coord_refs(sym, killed, f"term operator={op}")
        result["coefficient_symbolic"] = _remap_coord_string(sym, propagation_axis)

    if "coordinate_dependent" in result:
        result["coordinate_dependent"] = _remap_coord_deps(
            result["coordinate_dependent"],
            propagation_axis,
            killed,
            f"term operator={op}",
        )

    return result


def _transform_hamiltonian_term(
    term: dict[str, Any],
    remap: dict[str, str | None],
    propagation_axis: str,
    killed: set[str],
) -> dict[str, Any] | None:
    """Transform a single Hamiltonian bilinear term."""
    for factor_key in ("factor_a", "factor_b"):
        factor = term.get(factor_key, {})
        op = factor.get("operator", "")
        if _remap_operator(op, remap) is None:
            return None

    result = copy.deepcopy(term)

    for factor_key in ("factor_a", "factor_b"):
        factor = result[factor_key]
        op = factor["operator"]
        new_op = _remap_operator(op, remap)
        factor["operator"] = new_op  # type: ignore[assignment]

    if "coefficient_symbolic" in result:
        sym = result["coefficient_symbolic"]
        _check_killed_coord_refs(sym, killed, "hamiltonian_term coefficient")
        result["coefficient_symbolic"] = _remap_coord_string(sym, propagation_axis)

    if "coordinate_dependent" in result:
        result["coordinate_dependent"] = _remap_coord_deps(
            result["coordinate_dependent"],
            propagation_axis,
            killed,
            "hamiltonian_term",
        )

    return result


def _reduce_equations(
    result: dict[str, Any],
    remap: dict[str, str | None],
    propagation_axis: str,
    killed: set[str],
) -> tuple[set[str], list[str]]:
    """Transform equations and return (surviving_fields, eliminated_fields)."""
    original_fields: set[str] = {eq["field"] for eq in result.get("equations", [])}
    new_equations: list[dict[str, Any]] = []

    for eq in result.get("equations", []):
        terms = eq.get("rhs", {}).get("terms", [])
        new_terms = [
            t
            for term in terms
            if (t := _transform_term(term, remap, propagation_axis, killed)) is not None
        ]
        if not new_terms:
            continue
        new_eq = copy.deepcopy(eq)
        new_eq["rhs"]["terms"] = new_terms
        new_equations.append(new_eq)

    surviving: set[str] = {eq["field"] for eq in new_equations}
    result["equations"] = new_equations

    # Filter and reindex fields
    result["fields"] = [f for f in result.get("fields", []) if f["name"] in surviving]
    for idx, field in enumerate(result["fields"]):
        field["index"] = idx

    return surviving, sorted(original_fields - surviving)


def _reduce_hamiltonian(
    canonical: dict[str, Any],
    surviving_fields: set[str],
    remap: dict[str, str | None],
    propagation_axis: str,
    killed: set[str],
) -> None:
    """Filter and remap hamiltonian_terms in-place."""
    if "hamiltonian_terms" not in canonical:
        return
    new_terms: list[dict[str, Any]] = []
    for hterm in canonical["hamiltonian_terms"]:
        fa_field = hterm.get("factor_a", {}).get("field", "")
        fb_field = hterm.get("factor_b", {}).get("field", "")
        if fa_field not in surviving_fields or fb_field not in surviving_fields:
            continue
        transformed = _transform_hamiltonian_term(
            hterm,
            remap,
            propagation_axis,
            killed,
        )
        if transformed is not None:
            new_terms.append(transformed)
    canonical["hamiltonian_terms"] = new_terms


def _reduce_volume_element(
    canonical: dict[str, Any],
    propagation_axis: str,
    killed: set[str],
) -> None:
    """Handle volume element: remap or error if incompatible.

    Raises
    ------
    ValueError
        If volume element references killed coordinates.
    """
    if "volume_element" not in canonical:
        return

    vol_expr: str = canonical["volume_element"]
    for match in _COORD_REF_RE.finditer(vol_expr):
        if match.group(1) in killed:
            msg = (
                f"Volume element '{vol_expr}' references killed coordinate(s). "
                f"The Wolfram pipeline's factored volume element computation "
                f"should have handled this. The reduction may be incompatible "
                f"with this coordinate system for energy measurement."
            )
            raise ValueError(msg)

    new_vol = _remap_coord_string(vol_expr, propagation_axis)
    if new_vol.strip() in {"1", "1.", "1.0"}:
        del canonical["volume_element"]
    else:
        canonical["volume_element"] = new_vol


def _filter_coupling_matrices(
    result: dict[str, Any],
    surviving_indices: list[int],
) -> None:
    """Remove rows/columns for eliminated fields from coupling matrices.

    Parameters
    ----------
    result : dict
        Spec dict containing a ``"coupling"`` section.
    surviving_indices : list[int]
        Original field indices that survive (in order).
    """
    coupling: dict[str, Any] = result.get("coupling", {})
    if not coupling:
        return

    for matrix_key in (
        "mass_matrix",
        "coupling_matrix",
        "mass_matrix_symbolic",
        "coupling_matrix_symbolic",
    ):
        if matrix_key not in coupling:
            continue
        raw = coupling[matrix_key]
        if not isinstance(raw, list) or not raw:
            continue
        matrix = cast("list[list[Any]]", raw)
        new_matrix: list[list[Any]] = []
        for i in surviving_indices:
            if i < len(matrix):
                row = matrix[i]
                new_row: list[Any] = [row[j] for j in surviving_indices if j < len(row)]
                new_matrix.append(new_row)
        coupling[matrix_key] = new_matrix


def reduce_spec(
    spec_data: dict[str, Any],
    reduction_config: dict[str, Any],
) -> dict[str, Any]:
    """Apply plane-wave dimensional reduction to a JSON equation spec.

    Transforms a higher-dimensional spec into a clean 1+1D spec by:

    1. Remapping operator names (``laplacian_z → laplacian_x``)
    2. Removing terms with killed-axis operators
    3. Updating spacetime metadata (dimension, signature, coordinates)
    4. Remapping coordinate references in coefficient expressions
    5. Handling volume element (keep if surviving-coord-only, error if not)
    6. Adding provenance metadata

    Parameters
    ----------
    spec_data : dict
        Parsed JSON spec (as loaded by ``json.loads``).
    reduction_config : dict
        The ``[reduction]`` section from the TOML config.

    Returns
    -------
    dict
        Transformed spec with reduced dimension.

    """
    result = copy.deepcopy(spec_data)

    propagation_axis: str = reduction_config["propagation_axis"]
    spacetime: dict[str, Any] = result["spacetime"]
    coords: list[str] = spacetime["coordinates"]
    spatial_coords = [c for c in coords if c != "t"]
    killed = {c for c in spatial_coords if c != propagation_axis}

    remap = _build_operator_remap(propagation_axis, spatial_coords)

    surviving_fields, eliminated_fields = _reduce_equations(
        result,
        remap,
        propagation_axis,
        killed,
    )

    canonical: dict[str, Any] = result.get("canonical", {})
    _reduce_hamiltonian(canonical, surviving_fields, remap, propagation_axis, killed)
    _reduce_volume_element(canonical, propagation_axis, killed)

    spacetime["dimension"] = 2
    spacetime["signature"] = [-1, 1]
    spacetime["coordinates"] = ["t", "x"]

    surviving_indices = [
        i
        for i, eq in enumerate(spec_data.get("equations", []))
        if str(eq.get("field", "")) in surviving_fields
    ]
    _filter_coupling_matrices(result, surviving_indices)

    metadata: dict[str, Any] = result.get("metadata", {})
    metadata["reduction"] = {
        "type": reduction_config["type"],
        "original_dimension": spec_data["spacetime"]["dimension"],
        "propagation_axis": propagation_axis,
        "eliminated_fields": eliminated_fields,
    }
    result["metadata"] = metadata

    return result
