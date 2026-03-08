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
import warnings
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


# ---------------------------------------------------------------------------
# Degenerate algebraic constraint elimination
# ---------------------------------------------------------------------------
#
# After gauge fixing (TT, Lorenz, etc.) and plane-wave reduction, the derived
# equation system may contain "degenerate algebraic constraints" — algebraic
# equations where the constraint field cancels from its own RHS, leaving an
# implicit relation between dynamical fields.
#
# Example (TT gauge traceless condition):
#   h_9 = h_4 + h_7 + h_9   →   0 = h_4 + h_7   →   h_7 = -h_4
#
# The DAE solver (IDA) treats this as constraining h_9 = -(h_4 + h_7), but
# does NOT enforce h_4 + h_7 = 0.  If the dynamical equations for h_4 and h_7
# have cross-field Laplacians (common in gauge-fixed gravity), the mode
# u = h_4 + h_7 satisfies an elliptic equation d²t(u) = -∂²x(u) and grows
# exponentially from any numerical perturbation.
#
# The fix: substitute the constraint (h_7 → -h_4) throughout all equations,
# eliminating the dependent field entirely.  This removes the unstable mode,
# reduces system size, and produces a well-posed system for any solver backend.
#
# The algorithm handles:
# - Multi-round elimination (constraint A may expose constraint B)
# - N-field constraints (0 = c₁f₁ + c₂f₂ + ... + cₙfₙ), not just 2-field
# - Term merging after substitution (combines like terms on same operator+field)
# - Trivial equation cleanup (equations that become 0=0 after substitution)
# - Coupling matrix updates

#: Tolerance for detecting self_coeff ≈ 1.0 in degenerate constraints.
_DEGENERATE_TOL = 1e-12
#: Near-zero threshold for coefficients and term filtering.
_ZERO_TOL = 1e-15


def _find_degenerate_constraint(  # noqa: C901
    equations: list[dict[str, Any]],
) -> tuple[str, str, dict[str, float]] | None:
    """Find one degenerate algebraic constraint.

    Returns ``(constraint_eq_field, dep_field, {target_field: coeff, ...})``
    where ``dep_field = sum(coeff * target_field for all targets)``,
    or ``None`` if no degenerate constraint exists.

    The substitution means: replace ``dep_field`` with the linear combination.
    """
    for eq in equations:
        lhs = eq.get("lhs", {})
        order = lhs.get("order", {})

        # Must be algebraic (time=0, space=0)
        if order.get("time", 0) != 0 or order.get("space", 0) != 0:
            continue

        constraint_field = eq["field"]
        terms = eq.get("rhs", {}).get("terms", [])

        # All terms must use identity operator (pure algebraic, no derivatives)
        if not all(t.get("operator") == "identity" for t in terms):
            continue

        # Sum coefficient of constraint_field on RHS
        self_coeff = sum(
            t["coefficient"] for t in terms if t["field"] == constraint_field
        )

        # Degenerate: self_coeff ≈ 1.0 (cancels with LHS)
        if abs(self_coeff - 1.0) > _DEGENERATE_TOL:
            continue

        # Collect coefficients of other fields
        other: dict[str, float] = {}
        for t in terms:
            if t["field"] != constraint_field:
                other[t["field"]] = other.get(t["field"], 0.0) + t["coefficient"]

        if not other:
            continue  # Trivial identity h_9 = h_9, nothing to eliminate

        # The constraint is: 0 = Σ other[f] * f
        # Solve for the field with largest |coefficient| (for numerical stability)
        dep_field = max(other, key=lambda f: abs(other[f]))
        dep_coeff = other[dep_field]

        if abs(dep_coeff) < _ZERO_TOL:
            continue  # Degenerate coefficient

        # dep_field = -Σ (other[f]/dep_coeff) * f  for f ≠ dep_field
        substitution: dict[str, float] = {}
        for f, c in other.items():
            if f != dep_field:
                substitution[f] = -c / dep_coeff

        return constraint_field, dep_field, substitution

    return None


#: Spatial derivative operators whose null space (with periodic BCs) is
#: the constant function.  ``op(f) = 0`` with periodic BCs → ``f = const``.
#: With zero-mean initial data this gives ``f = 0``.
_DERIVATIVE_ZERO_OPS: frozenset[str] = frozenset(
    {
        "gradient_x",
        "gradient_y",
        "gradient_z",
        "laplacian",
        "laplacian_x",
        "laplacian_y",
        "laplacian_z",
    }
)


def _find_gradient_zero_constraint(
    equations: list[dict[str, Any]],
) -> str | None:
    """Find a field constrained to zero by a derivative self-constraint.

    Detects algebraic equations of the form ``op(f) = 0`` where ``op`` is a
    spatial derivative operator and ``f`` is the equation's own field.  With
    periodic BCs and zero-mean initial data, ``gradient(f) = 0`` or
    ``laplacian(f) = 0`` implies ``f = 0``.

    These arise from TT gauge transverse constraints after plane-wave
    reduction simplifies multi-field divergence conditions to single-term
    self-referencing derivative equations.

    Returns the field name, or ``None`` if no such constraint exists.
    """
    for eq in equations:
        lhs = eq.get("lhs", {})
        order = lhs.get("order", {})

        # Must be algebraic (time=0, space=0)
        if order.get("time", 0) != 0 or order.get("space", 0) != 0:
            continue

        terms = eq.get("rhs", {}).get("terms", [])

        # Single self-referencing derivative term
        if len(terms) != 1:
            continue
        term = terms[0]
        if term["field"] != eq["field"]:
            continue
        if term.get("operator") not in _DERIVATIVE_ZERO_OPS:
            continue

        # op(field) = 0 with periodic BCs → field = 0
        return eq["field"]

    return None


def _negate_symbolic(expr: str) -> str:
    """Negate a symbolic coefficient expression: ``expr`` → ``-(expr)``."""
    stripped = expr.strip()
    if stripped.startswith("-") and "+" not in stripped and "-" not in stripped[1:]:
        return stripped[1:]  # Simple case: "-X" → "X"
    return f"-({stripped})"


def _scale_symbolic(expr: str, factor: float) -> str:
    """Scale a symbolic expression by a numeric factor."""
    if factor == -1.0:
        return _negate_symbolic(expr)
    if factor == 1.0:
        return expr
    return f"({factor})*({expr})"


def _sum_symbolic(a: str, b: str) -> str:
    """Sum two symbolic coefficient expressions."""
    return f"({a}) + ({b})"


def _substitute_in_terms(  # noqa: C901
    terms: list[dict[str, Any]],
    dep_field: str,
    substitution: dict[str, float],
) -> list[dict[str, Any]]:
    """Replace references to ``dep_field`` in a list of RHS terms.

    For each term referencing ``dep_field`` with operator ``op`` and
    coefficient ``c``: replace with one term per substitution entry,
    ``c * sub_coeff * op(target_field)``.

    Then merge terms with identical ``(operator, field)`` keys.
    Symbolic coefficients are preserved through transformation and merging.
    """
    expanded: list[dict[str, Any]] = []

    for term in terms:
        if term["field"] != dep_field:
            expanded.append(copy.deepcopy(term))
            continue

        # Expand: dep_field → Σ sub_coeff * target
        orig_coeff = term["coefficient"]
        orig_symbolic = term.get("coefficient_symbolic")

        if not substitution:
            # dep_field = 0 (constraint says field is zero)
            continue

        for target_field, sub_coeff in substitution.items():
            new_term = copy.deepcopy(term)
            new_term["field"] = target_field
            new_term["coefficient"] = orig_coeff * sub_coeff
            # Transform symbolic coefficient
            if orig_symbolic is not None:
                new_term["coefficient_symbolic"] = _scale_symbolic(
                    orig_symbolic,
                    sub_coeff,
                )
            else:
                new_term.pop("coefficient_symbolic", None)
            expanded.append(new_term)

    # Merge terms with same (operator, field)
    merged: dict[tuple[str, str], dict[str, Any]] = {}
    for term in expanded:
        key = (term["operator"], term["field"])
        if key in merged:
            existing = merged[key]
            existing["coefficient"] += term["coefficient"]
            # Merge symbolic coefficients
            sym_a = existing.get("coefficient_symbolic")
            sym_b = term.get("coefficient_symbolic")
            if sym_a is not None and sym_b is not None:
                existing["coefficient_symbolic"] = _sum_symbolic(sym_a, sym_b)
            elif sym_b is not None:
                # Only new term has symbolic — use numeric from existing
                existing["coefficient_symbolic"] = _sum_symbolic(
                    str(existing["coefficient"] - term["coefficient"]),
                    sym_b,
                )
            elif sym_a is not None:
                # Only existing has symbolic — add numeric from new
                existing["coefficient_symbolic"] = _sum_symbolic(
                    sym_a,
                    str(term["coefficient"]),
                )
        else:
            merged[key] = term

    # Remove near-zero terms
    return [t for t in merged.values() if abs(t["coefficient"]) > _ZERO_TOL]


def eliminate_degenerate_constraints(  # noqa: C901, PLR0912, PLR0914, PLR0915
    spec_data: dict[str, Any],
) -> dict[str, Any]:
    """Eliminate fields determined by degenerate algebraic constraints.

    Scans for algebraic equations where the constraint field cancels from its
    own RHS, leaving an implicit relation between dynamical fields.  Substitutes
    the dependent field throughout all equations and removes it from the system.

    This is applied iteratively — eliminating one constraint may expose another.

    Parameters
    ----------
    spec_data : dict
        Parsed JSON equation spec.

    Returns
    -------
    dict
        Spec with degenerate constraints eliminated.  Returns the original
        spec unchanged if no degenerate constraints are found.
    """
    result = copy.deepcopy(spec_data)
    equations: list[dict[str, Any]] = result.get("equations", [])
    eliminated: list[str] = []

    # Build original field name → index map for coupling matrix filtering
    original_field_names = [f["name"] for f in result.get("fields", [])]

    while True:
        found = _find_degenerate_constraint(equations)
        if found is None:
            break

        constraint_eq_field, dep_field, substitution = found

        # Build velocity substitution: v_{dep_field} → Σ coeff * v_{target}
        vel_dep = f"v_{dep_field}"
        vel_substitution: dict[str, float] = {
            f"v_{target}": coeff for target, coeff in substitution.items()
        }

        # Phase 1: Substitute dep_field AND v_{dep_field} in all equation terms
        new_equations: list[dict[str, Any]] = []
        for eq in equations:
            field = eq["field"]

            # Remove the degenerate constraint equation itself
            if field == constraint_eq_field:
                continue

            # Remove the equation *for* the dependent field
            if field == dep_field:
                continue

            # Substitute in RHS terms (both field and velocity)
            terms = eq.get("rhs", {}).get("terms", [])
            new_terms = _substitute_in_terms(terms, dep_field, substitution)
            new_terms = _substitute_in_terms(new_terms, vel_dep, vel_substitution)

            # If no terms remain, equation is trivially 0=0 → remove
            if not new_terms:
                continue

            eq["rhs"]["terms"] = new_terms
            new_equations.append(eq)

        equations = new_equations
        eliminated.append(dep_field)

        # Also track the constraint equation field if it disappeared
        if constraint_eq_field != dep_field and constraint_eq_field not in {
            eq["field"] for eq in equations
        }:
            eliminated.append(constraint_eq_field)

    # --- Gradient-zero elimination ---
    # After degenerate pairs are gone, detect derivative self-constraints
    # like ``gradient_z(h_6) = 0`` that imply ``h_6 = 0`` with periodic BCs.
    # Eliminating such fields may expose further degenerate constraints or
    # trivial (0=0) equations, so we iterate.
    while True:
        zero_field = _find_gradient_zero_constraint(equations)
        if zero_field is None:
            break

        # Remove the zero-field's equation and substitute zero everywhere
        new_equations = []
        for eq in equations:
            if eq["field"] == zero_field:
                continue

            terms = eq.get("rhs", {}).get("terms", [])
            new_terms = _substitute_in_terms(terms, zero_field, {})
            # Also substitute the velocity form
            new_terms = _substitute_in_terms(
                new_terms, f"v_{zero_field}", {}
            )

            if not new_terms:
                # Equation became 0=0 → remove
                continue

            eq["rhs"]["terms"] = new_terms
            new_equations.append(eq)

        equations = new_equations
        eliminated.append(zero_field)

    if not eliminated:
        return spec_data  # No changes

    # Validate: non-trivial coefficients should retain symbolic expressions.
    # This catches bugs in the symbolic merging logic.
    for eq in equations:
        for term in eq.get("rhs", {}).get("terms", []):
            coeff = term.get("coefficient", 0.0)
            if abs(coeff) not in {0.0, 1.0} and "coefficient_symbolic" not in term:
                warnings.warn(
                    f"Constraint elimination: term {coeff}*{term['operator']}"
                    f"({term['field']}) in equation for '{eq['field']}' lost its "
                    f"coefficient_symbolic during elimination. This may cause "
                    f"incorrect results when parameters are overridden at runtime.",
                    UserWarning,
                    stacklevel=2,
                )

    result["equations"] = equations

    # Phase 2: Update fields list and reindex
    surviving_names = {eq["field"] for eq in equations}
    result["fields"] = [
        f for f in result.get("fields", []) if f["name"] in surviving_names
    ]
    for idx, field in enumerate(result["fields"]):
        field["index"] = idx

    # Phase 3: Filter coupling matrices
    elim_set = set(eliminated)
    surviving_original_indices = [
        i
        for i, name in enumerate(original_field_names)
        if name not in elim_set and name in surviving_names
    ]
    _filter_coupling_matrices(result, surviving_original_indices)

    # Phase 4: Filter hamiltonian terms
    canonical = result.get("canonical", {})
    if "hamiltonian_terms" in canonical:
        canonical["hamiltonian_terms"] = [
            ht
            for ht in canonical["hamiltonian_terms"]
            if ht.get("factor_a", {}).get("field") not in elim_set
            and ht.get("factor_b", {}).get("field") not in elim_set
        ]

    # Phase 5: Record metadata
    metadata = result.get("metadata", {})
    metadata["constraint_elimination"] = {
        "eliminated_fields": eliminated,
    }

    return result
