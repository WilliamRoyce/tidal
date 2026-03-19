"""TOML configuration validation for ``tidal derive``.

Validates the structure and consistency of theory TOML files before
Wolfram code generation.  Each ``_validate_*`` function checks one
TOML section; ``_validate_config`` orchestrates all checks.

This module is pure validation — it does not emit Wolfram code or
use ``_WlsContext``.  Separated from ``_derive.py`` for readability
and single-responsibility.
"""

from __future__ import annotations

import re
from typing import Any, cast

# ---------------------------------------------------------------------------
# Constants and presets
# ---------------------------------------------------------------------------

_VALID_FIELD_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")

_MIN_DIM = 2
_MAX_DIM = 7
_MIN_PREFIX_LEN = 2

# Names reserved for built-in operators and xAct auto-created tensors.
# Field, constant, and perturbation names must not collide with these.
# The substitution logic in _substitute_field_names replaces these globally,
# so a field named "CD" would be double-prefixed or corrupted.
_RESERVED_NAMES: frozenset[str] = frozenset(
    {
        # Covariant derivatives
        "CD",
        "CDT",
        # Metric and chart
        "eta",
        "bg",
        "chart",
        # xAct auto-created curvature/torsion tensors (case-insensitive first letter)
        "Riemann",
        "Ricci",
        "RicciScalar",
        "Einstein",
        "Weyl",
        "Schouten",
        "Cotton",
        "Kretschner",
        "Torsion",
        "TorsionCDT",
        # Other reserved
        "PD",
        "Christoffel",
    }
)

# Gauge-fixing preset registry.  Each entry maps a preset name to its
# mechanism ("lagrangian_term" or "constraint"), and the required field type.
# Type A (lagrangian_term) presets also have a "builder" key naming the
# Wolfram function in GaugeFix.wl.  Type B (constraint) presets are
# implemented at the component level in Python (_wls_gauge_fixing_type_b).
_GAUGE_PRESETS: dict[str, dict[str, str]] = {
    "lorenz": {
        "mechanism": "lagrangian_term",
        "builder": "BuildLorenzGaugeTerm",
        "requires": "vector",
    },
    "de_donder": {
        "mechanism": "lagrangian_term",
        "builder": "BuildDeDonderGaugeTerm",
        "requires": "tensor",
    },
    "temporal": {
        "mechanism": "constraint",
        "requires": "vector",
    },
    "coulomb": {
        "mechanism": "constraint",
        "requires": "vector",
    },
    "axial": {
        "mechanism": "constraint",
        "requires": "vector",
    },
    "tt": {
        "mechanism": "constraint",
        "requires": "tensor",
    },
}

# Spatial coordinates per dimension (excluding time) — used by _validate_reduction
_SPATIAL_COORDS: dict[int, list[str]] = {
    2: ["x"],
    3: ["x", "y"],
    4: ["x", "y", "z"],
}

# Must match _INDEX_LETTERS in _derive.py — used here only for max-rank check
_INDEX_LETTERS = list("abcdefghijklmnop")  # 16 letters, enough for dim 7 + 4 = 11


# ---------------------------------------------------------------------------
# Section validators
# ---------------------------------------------------------------------------


def _validate_spacetime(config: dict[str, Any]) -> None:
    """Validate [spacetime] section of TOML config.

    Raises
    ------
    ValueError
        If dimension is missing, non-integer, or out of range.
    """
    dim = config["spacetime"].get("dimension")
    if dim is None:
        msg = "[spacetime] must specify 'dimension'"
        raise ValueError(msg)
    if not isinstance(dim, int) or dim < _MIN_DIM or dim > _MAX_DIM:
        msg = f"[spacetime].dimension must be integer {_MIN_DIM}-{_MAX_DIM}, got: {dim}"
        raise ValueError(msg)


def _check_name_valid(name: str, context: str) -> None:
    """Check that a name is a valid identifier and not reserved.

    Raises
    ------
    ValueError
        If the name is not alphanumeric or is reserved.
    """
    if not _VALID_FIELD_NAME.match(name):
        msg = f"{context} '{name}' must be alphanumeric starting with a letter"
        raise ValueError(msg)
    if name in _RESERVED_NAMES:
        msg = (
            f"{context} '{name}' is reserved (collides with built-in operator). "
            f"Reserved names: {sorted(_RESERVED_NAMES)}"
        )
        raise ValueError(msg)


def _validate_single_field(field: dict[str, Any], index: int, dim: int) -> None:
    """Validate a single [[fields]] entry.

    Raises
    ------
    ValueError
        If field definition is missing required keys or has invalid values.
    """
    if "name" not in field:
        msg = f"[[fields]] entry {index} missing 'name'"
        raise ValueError(msg)
    fname = field["name"]
    _check_name_valid(fname, "Field name")
    if "type" not in field:
        msg = f"[[fields]] entry {index} ('{fname}') missing 'type'"
        raise ValueError(msg)
    ftype = field["type"]
    if ftype not in {"scalar", "vector", "tensor"}:
        msg = f"Unknown field type '{ftype}' for '{fname}'. Use: scalar, vector, tensor"
        raise ValueError(msg)
    if ftype == "tensor":
        if "rank" not in field:
            msg = f"Tensor field '{fname}' must specify 'rank'"
            raise ValueError(msg)
        if "symmetry" not in field:
            msg = f"Tensor field '{fname}' must specify 'symmetry' (e.g. 'antisymmetric', 'symmetric', 'none')"
            raise ValueError(msg)
        if field["rank"] > dim:
            msg = f"Tensor field '{fname}' rank {field['rank']} exceeds spacetime dimension {dim}"
            raise ValueError(msg)
        max_rank = len(_INDEX_LETTERS)
        if field["rank"] > max_rank:
            msg = f"Tensor field '{fname}' rank {field['rank']} exceeds maximum supported rank {max_rank}"
            raise ValueError(msg)


def _validate_fields(config: dict[str, Any]) -> None:
    """Validate [[fields]] entries of TOML config.

    Raises
    ------
    ValueError
        If fields are missing, empty, or have invalid entries.
    """
    if "fields" not in config or not config["fields"]:
        msg = "Must define at least one [[fields]] entry"
        raise ValueError(msg)

    dim = config["spacetime"]["dimension"]
    for i, field in enumerate(config["fields"]):
        _validate_single_field(field, i, dim)


def _validate_constants(config: dict[str, Any]) -> None:
    """Validate optional [constants] section of TOML config.

    Constant names flow directly into Wolfram ``DefConstantSymbol[name]``
    calls, so they must be valid Wolfram identifiers (alphanumeric,
    starting with a letter, no underscores).

    """
    names: list[str] = config.get("constants", {}).get("names", [])
    for name in names:
        _check_name_valid(name, "Constant name")


def _validate_lagrangian(config: dict[str, Any]) -> None:
    """Validate [lagrangian] section of TOML config.

    Raises
    ------
    ValueError
        If expression is missing, empty, or contains covariant derivatives
        of background fields (which are not supported).
    """
    expr = config["lagrangian"].get("expression")
    if not expr or not expr.strip():
        msg = "[lagrangian].expression must be a non-empty string"
        raise ValueError(msg)

    # Check for covariant derivatives of background fields.
    # CD[-a][G[]] or CD[a][G[-b]] etc. would silently produce wrong results
    # because the scalar ReplaceAll substitution doesn't handle derivatives.
    bg_names = [f["name"] for f in config.get("background_fields", [])]
    for bg_name in bg_names:
        # Match CD[...][<bg_name>[...]] — any index pattern around the BG field
        pattern = rf"CD\[[^\]]*\]\s*\[\s*{re.escape(bg_name)}\b"
        if re.search(pattern, expr):
            msg = (
                f"Lagrangian contains a covariant derivative of background "
                f"field '{bg_name}' (e.g., CD[-a][{bg_name}[...]]). "
                f"Derivatives of background fields are not yet supported — "
                f"they would produce incorrect equations of motion. "
                f"Use the background field directly (without derivatives) "
                f"in the coupling term."
            )
            raise ValueError(msg)


def _validate_parameters(config: dict[str, Any]) -> None:
    """Validate optional [parameters] section of TOML config.

    Raises
    ------
    TypeError
        If parameters is not a dict or contains non-numeric values.
    ValueError
        If keys aren't declared constants.
    """
    params: dict[str, object] | None = config.get("parameters")
    if params is None:
        return
    if not isinstance(params, dict):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = "[parameters] must be a table of key = value pairs"
        raise TypeError(msg)
    const_names = set(config.get("constants", {}).get("names", []))
    for key, val in params.items():
        if not isinstance(val, (int, float)):
            msg = f"[parameters].{key} must be numeric, got {type(val).__name__}"
            raise TypeError(msg)
        if const_names and key not in const_names:
            msg = f"[parameters].{key} is not declared in [constants].names"
            raise ValueError(msg)


def _validate_derived_fields(config: dict[str, Any]) -> None:
    """Validate optional [[derived_fields]] entries.

    Raises
    ------
    ValueError
        If a derived field has invalid type/rank/symmetry, missing definition,
        or name collision with a fundamental field.
    """
    derived = config.get("derived_fields", [])
    if not derived:
        return

    dim = config["spacetime"]["dimension"]
    fundamental_names = {f["name"] for f in config.get("fields", [])}

    for i, field in enumerate(derived):
        _validate_single_field(field, i, dim)

        # Require definition
        defn = field.get("definition")
        if not defn or not isinstance(defn, str) or not defn.strip():
            msg = f"[[derived_fields]] entry {i} ('{field.get('name', '?')}') must have a non-empty 'definition'"
            raise ValueError(msg)

        # Name collision check
        if field["name"] in fundamental_names:
            msg = f"Derived field '{field['name']}' conflicts with a fundamental [[fields]] entry"
            raise ValueError(msg)


def _expected_component_count(field: dict[str, Any], dim: int) -> int:
    """Return the expected number of components for a field given spacetime dimension."""
    ftype = field["type"]
    if ftype == "scalar":
        return 1
    if ftype == "vector":
        return dim
    # Tensor: dim^rank (full, no symmetry reduction for component values)
    rank = field.get("rank", 2)
    return dim**rank


def _validate_background_fields(config: dict[str, Any]) -> None:
    """Validate optional ``[[background_fields]]`` entries.

    Background fields are non-dynamical tensors that appear in the Lagrangian
    but are NOT varied in the Euler-Lagrange derivation.  They survive as
    (possibly position-dependent) coefficients in the equations of motion.

    Each entry must have a ``components`` key specifying component values:
    - scalar: ``components = ["B0val"]`` (1 element)
    - vector: ``components = [0, 0, "B0val"]`` (dim elements)

    Component values can be numbers (0, 1.0) or symbolic Wolfram expressions
    (``"B0val"``, ``"B0val * Sin[2*Pi*x[]/L]"``).  Constants used in
    expressions must be declared in ``[constants]``.

    Raises
    ------
    ValueError
        If a background field has invalid type/rank/symmetry, missing
        components, wrong component count, or name collision.
    TypeError
        If ``components`` is not a list or contains non-numeric/non-string values.
    """
    bg_fields = config.get("background_fields", [])
    if not bg_fields:
        return

    dim = config["spacetime"]["dimension"]
    dynamical_names = {f["name"] for f in config.get("fields", [])}
    derived_names = {f["name"] for f in config.get("derived_fields", [])}

    for i, field in enumerate(bg_fields):
        _validate_single_field(field, i, dim)

        fname = field["name"]

        # Name collision check
        if fname in dynamical_names:
            msg = f"Background field '{fname}' conflicts with a dynamical [[fields]] entry"
            raise ValueError(msg)
        if fname in derived_names:
            msg = (
                f"Background field '{fname}' conflicts with a [[derived_fields]] entry"
            )
            raise ValueError(msg)

        # Require components
        comps: list[int | float | str] | None = field.get("components")
        if comps is None:
            msg = (
                f"[[background_fields]] entry {i} ('{fname}') must have "
                f"'components' specifying component values"
            )
            raise ValueError(msg)
        if not isinstance(comps, list):  # pyright: ignore[reportUnnecessaryIsInstance]
            msg = f"[[background_fields]] entry {i} ('{fname}'): 'components' must be a list"
            raise TypeError(msg)

        expected = _expected_component_count(field, dim)
        if len(comps) != expected:
            msg = (
                f"[[background_fields]] entry {i} ('{fname}'): expected "
                f"{expected} components for {field['type']} in {dim}D, got {len(comps)}"
            )
            raise ValueError(msg)

        # Validate component types (numbers or strings)
        for j, comp in enumerate(comps):
            if not isinstance(comp, (int, float, str)):  # pyright: ignore[reportUnnecessaryIsInstance]
                msg = (
                    f"[[background_fields]] entry {i} ('{fname}'): "
                    f"component {j} must be a number or string expression, "
                    f"got {type(comp).__name__}"
                )
                raise TypeError(msg)


def _validate_matter_perturbations(
    config: dict[str, Any],
    matter_perts: list[dict[str, Any]],
    *,
    has_lagrangian: bool = False,
) -> None:
    """Validate ``[[linearization.matter_perturbations]]`` entries.

    Each entry defines a field split: field = background + ε·perturbation.
    Requires ``field`` (in ``[[fields]]``), ``perturbation_name`` (unique),
    and ``background`` (in ``[[background_fields]]``).

    Raises
    ------
    ValueError
        If entries are invalid, reference non-existent fields/backgrounds,
        or have duplicate names.
    """
    if not has_lagrangian:
        msg = (
            "[[linearization.matter_perturbations]] requires [lagrangian] "
            "(Lagrangian-first linearization path)"
        )
        raise ValueError(msg)

    field_names = {f["name"] for f in config.get("fields", [])}
    bg_names = {f["name"] for f in config.get("background_fields", [])}
    seen_fields: set[str] = set()
    seen_pnames: set[str] = set()

    for i, mp in enumerate(matter_perts):
        _validate_single_matter_perturbation(
            i,
            mp,
            field_names,
            bg_names,
            (seen_fields, seen_pnames),
        )


def _validate_single_matter_perturbation(
    i: int,
    mp: dict[str, Any],
    field_names: set[str],
    bg_names: set[str],
    seen: tuple[set[str], set[str]],
) -> None:
    """Validate a single ``[[linearization.matter_perturbations]]`` entry.

    Raises
    ------
    ValueError
        If the entry is missing required keys or references invalid names.
    """
    seen_fields, seen_pnames = seen
    tag = f"[[linearization.matter_perturbations]] entry {i}"

    # Require 'field'
    mf = mp.get("field")
    if not mf or not isinstance(mf, str):
        msg = f"{tag}: 'field' is required"
        raise ValueError(msg)
    if mf not in field_names:
        msg = f"{tag}: field '{mf}' not found in [[fields]]"
        raise ValueError(msg)
    if mf in seen_fields:
        msg = f"{tag}: field '{mf}' already has a perturbation defined"
        raise ValueError(msg)
    seen_fields.add(mf)

    # Require 'perturbation_name'
    pname = mp.get("perturbation_name")
    if not pname or not isinstance(pname, str):
        msg = f"{tag}: 'perturbation_name' is required"
        raise ValueError(msg)
    if pname in seen_pnames:
        msg = f"{tag}: perturbation_name '{pname}' already used"
        raise ValueError(msg)
    seen_pnames.add(pname)

    # Require 'background'
    bg = mp.get("background")
    if not bg or not isinstance(bg, str):
        msg = (
            f"{tag}: 'background' is required "
            f"(name of matching [[background_fields]] entry)"
        )
        raise ValueError(msg)
    if bg not in bg_names:
        msg = f"{tag}: background '{bg}' not found in [[background_fields]]"
        raise ValueError(msg)


def _validate_linearization(
    config: dict[str, Any],
    *,
    has_lagrangian: bool = False,
) -> None:
    """Validate optional ``[linearization]`` section.

    Parameters
    ----------
    config : dict
        Full TOML config.
    has_lagrangian : bool
        Whether ``[lagrangian]`` is also present.  When *True*,
        ``[linearization].expression`` is not required (the Lagrangian
        provides the expression to linearize).

    Raises
    ------
    TypeError
        If the section is not a table.
    ValueError
        If required keys are missing or ``perturbation_field`` is not declared.
    """
    if "linearization" not in config:
        return
    lin: dict[str, Any] = config["linearization"]
    if not isinstance(lin, dict):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = "[linearization] must be a table"
        raise TypeError(msg)

    # expression is required only for legacy path (no [lagrangian])
    if not has_lagrangian:
        expr = lin.get("expression")
        if not expr or not isinstance(expr, str) or not expr.strip():
            msg = "[linearization].expression must be a non-empty string"
            raise ValueError(msg)

    pf = lin.get("perturbation_field")
    matter_perts: list[dict[str, Any]] = lin.get("matter_perturbations", [])

    # perturbation_field is required unless matter_perturbations are present
    # (matter-only perturbation without metric perturbation)
    if not matter_perts and (not pf or not isinstance(pf, str)):
        msg = "[linearization].perturbation_field is required"
        raise ValueError(msg)
    if pf is not None and isinstance(pf, str):
        field_names = {f["name"] for f in config.get("fields", [])}
        if pf not in field_names:
            msg = f"[linearization].perturbation_field '{pf}' not found in [[fields]]"
            raise ValueError(msg)

    if matter_perts:
        _validate_matter_perturbations(
            config, matter_perts, has_lagrangian=has_lagrangian
        )


def _validate_gauge_entry_preset(
    i: int,
    entry: dict[str, Any],
    field_info: dict[str, Any],
) -> None:
    """Validate a named gauge preset against the target field type.

    Raises
    ------
    ValueError
        If the field type is incompatible with the preset.
    """
    gauge_type = entry["type"]
    preset = _GAUGE_PRESETS[gauge_type]
    required_type = preset["requires"]
    if field_info["type"] == "scalar":
        msg = f"[[gauge]] entry {i}: scalars have no gauge freedom"
        raise ValueError(msg)
    if required_type == "vector" and field_info["type"] != "vector":
        msg = f"[[gauge]] entry {i}: '{gauge_type}' requires a vector field, got '{field_info['type']}'"
        raise ValueError(msg)
    if required_type == "tensor" and field_info["type"] != "tensor":
        msg = f"[[gauge]] entry {i}: '{gauge_type}' requires a tensor field, got '{field_info['type']}'"
        raise ValueError(msg)
    if gauge_type == "tt" and field_info.get("symmetry") != "symmetric":
        msg = f"[[gauge]] entry {i}: 'tt' gauge requires a symmetric tensor, got symmetry='{field_info.get('symmetry', 'none')}'"
        raise ValueError(msg)


def _validate_gauge_entry_custom(i: int, entry: dict[str, Any]) -> None:
    """Validate a custom gauge entry has ``mechanism`` and ``expression``.

    Raises
    ------
    ValueError
        If ``mechanism`` or ``expression`` is missing or invalid.
    """
    mechanism = entry.get("mechanism")
    if mechanism not in {"lagrangian_term", "constraint"}:
        msg = (
            f"[[gauge]] entry {i}: custom gauge requires "
            f"'mechanism' = 'lagrangian_term' or 'constraint'"
        )
        raise ValueError(msg)
    if mechanism == "constraint":
        msg = (
            f"[[gauge]] entry {i}: custom constraint gauges not yet supported; "
            f"use a named preset (temporal, coulomb, axial) or mechanism='lagrangian_term'"
        )
        raise ValueError(msg)
    expr = entry.get("expression")
    if not expr or not isinstance(expr, str) or not expr.strip():
        msg = f"[[gauge]] entry {i}: custom gauge requires non-empty 'expression'"
        raise ValueError(msg)


def _validate_gauge_entry_common(
    i: int,
    entry: dict[str, Any],
    field_map: dict[str, dict[str, Any]],
    seen_fields: set[str],
) -> tuple[str, str]:
    """Validate keys common to all ``[[gauge]]`` entries.

    Returns the ``(field_name, gauge_type)`` pair on success.

    Raises
    ------
    ValueError
        If required keys are missing, field is undeclared or duplicated,
        or ``xi`` is not a positive number.
    """
    if "field" not in entry or not isinstance(entry.get("field"), str):
        msg = f"[[gauge]] entry {i}: missing or invalid 'field'"
        raise ValueError(msg)
    if "type" not in entry or not isinstance(entry.get("type"), str):
        msg = f"[[gauge]] entry {i}: missing or invalid 'type'"
        raise ValueError(msg)

    field_name: str = entry["field"]
    gauge_type: str = entry["type"]

    if field_name not in field_map:
        msg = f"[[gauge]] entry {i}: field '{field_name}' not in [[fields]]"
        raise ValueError(msg)
    if field_name in seen_fields:
        msg = f"[[gauge]] entry {i}: duplicate gauge for field '{field_name}'"
        raise ValueError(msg)
    seen_fields.add(field_name)

    if "xi" in entry:
        xi = entry["xi"]
        if not isinstance(xi, (int, float)) or xi <= 0:
            msg = f"[[gauge]] entry {i}: 'xi' must be a positive number"
            raise ValueError(msg)

    return field_name, gauge_type


def _validate_gauge(config: dict[str, Any]) -> None:
    """Validate optional ``[[gauge]]`` entries.

    Each entry specifies a per-field gauge-fixing choice.  Named presets
    (``lorenz``, ``de_donder``, etc.) are validated against field type;
    ``type = "custom"`` requires ``mechanism`` and ``expression``.

    Raises
    ------
    TypeError
        If ``[[gauge]]`` is not a list.
    ValueError
        If a gauge entry references a non-existent field, uses an unknown
        type, targets a scalar (no gauge freedom), duplicates a field, or
        is missing required keys for custom gauges.
    """
    gauge_list: list[dict[str, Any]] = config.get("gauge", [])
    if not gauge_list:
        return
    if not isinstance(gauge_list, list):  # pyright: ignore[reportUnnecessaryIsInstance]
        msg = "[[gauge]] must be an array of tables"
        raise TypeError(msg)

    field_map = {f["name"]: f for f in config.get("fields", [])}

    # Also include matter perturbation names as valid gauge targets.
    # When [[linearization.matter_perturbations]] splits A → Ā + εa,
    # the gauge references the perturbation name "a" (the dynamical variable),
    # with the same type/rank/symmetry as the original field.
    lin = config.get("linearization", {})
    for mp in lin.get("matter_perturbations", []) if isinstance(lin, dict) else []:  # type: ignore[reportUnknownVariableType]
        pname: str = mp.get("perturbation_name", "")  # type: ignore[reportUnknownVariableType]
        mf_name: str = mp.get("field", "")  # type: ignore[reportUnknownVariableType]
        if pname and mf_name and mf_name in field_map:
            # Perturbation inherits type/rank/symmetry from the original field
            field_map[pname] = dict(field_map[mf_name], name=pname)  # type: ignore[reportUnknownArgumentType]

    seen_fields: set[str] = set()

    for i, entry in enumerate(gauge_list):
        field_name, gauge_type = _validate_gauge_entry_common(
            i,
            entry,
            field_map,
            seen_fields,
        )

        if gauge_type in _GAUGE_PRESETS:
            _validate_gauge_entry_preset(i, entry, field_map[field_name])
        elif gauge_type == "custom":
            _validate_gauge_entry_custom(i, entry)
        else:
            known = ", ".join([*_GAUGE_PRESETS, "custom"])
            msg = f"[[gauge]] entry {i}: unknown gauge type '{gauge_type}' (known: {known})"
            raise ValueError(msg)


def _validate_reduction(config: dict[str, Any]) -> None:
    """Validate optional ``[reduction]`` section.

    The plane-wave reduction zeroes all transverse spatial derivatives,
    reducing the spacetime dimension to 1+1D (or 2+1D for partial
    reductions in future).

    Raises
    ------
    ValueError
        If reduction config is invalid (bad type, axis, or dimension).
    TypeError
        If ``coordinate_values`` is not a dict.
    """
    reduction = config.get("reduction")
    if reduction is None:
        return

    rtype = reduction.get("type")
    if rtype != "plane_wave":
        msg = (
            f"[reduction] type must be 'plane_wave', got '{rtype}'. "
            "Other reduction types may be added in future."
        )
        raise ValueError(msg)

    dim: int = config["spacetime"]["dimension"]
    min_reducible_dim = 3  # need at least 2 spatial dimensions
    if dim < min_reducible_dim:
        msg = (
            f"[reduction] cannot reduce a {dim - 1}+1D theory further. "
            "Plane-wave reduction requires at least 2 spatial dimensions."
        )
        raise ValueError(msg)

    prop_axis = reduction.get("propagation_axis")
    if not prop_axis:
        msg = "[reduction] requires 'propagation_axis' (e.g., 'x', 'y', or 'z')"
        raise ValueError(msg)

    spatial = _SPATIAL_COORDS.get(dim)
    if spatial is None:
        msg = f"[reduction] unsupported dimension {dim}"
        raise ValueError(msg)

    if prop_axis not in spatial:
        msg = (
            f"[reduction] propagation_axis '{prop_axis}' is not a valid "
            f"spatial coordinate for dimension {dim}. "
            f"Valid axes: {spatial}"
        )
        raise ValueError(msg)

    # Optional coordinate_values: map of killed-axis names → Wolfram value strings
    coord_values_raw = reduction.get("coordinate_values", {})
    if not isinstance(coord_values_raw, dict):
        msg = "[reduction] coordinate_values must be a table (dict), e.g. {y = 'Pi/2'}"
        raise TypeError(msg)
    coord_values = dict(cast("dict[str, Any]", coord_values_raw))
    killed = [c for c in (spatial or []) if c != prop_axis]
    for key in coord_values:
        if key not in killed:
            msg = (
                f"[reduction] coordinate_values key '{key}' is not a killed "
                f"coordinate for propagation_axis '{prop_axis}'. "
                f"Valid killed coordinates: {killed}"
            )
            raise ValueError(msg)


# ---------------------------------------------------------------------------
# Orchestrator
# ---------------------------------------------------------------------------


def _validate_config(config: dict[str, Any]) -> None:  # type: ignore[reportUnusedFunction]
    """Validate TOML config structure.

    Raises
    ------
    ValueError
        If any required section or field is missing or invalid.
    """
    if "spacetime" not in config:
        msg = "Missing required section: [spacetime]"
        raise ValueError(msg)

    has_lagrangian = "lagrangian" in config
    has_linearization = "linearization" in config

    if not has_lagrangian and not has_linearization:
        msg = "Config must have either [lagrangian] or [linearization]"
        raise ValueError(msg)
    if has_linearization and not has_lagrangian:
        import warnings

        warnings.warn(
            "[linearization] without [lagrangian] is deprecated. "
            "Provide the Lagrangian in [lagrangian] and use "
            "[linearization] as a modifier.",
            DeprecationWarning,
            stacklevel=2,
        )
    if has_linearization and not has_lagrangian and config.get("gauge"):
        msg = (
            "[[gauge]] with [linearization] requires [lagrangian] "
            "(Lagrangian-first linearization path). "
            "Add a [lagrangian] section or remove [[gauge]]."
        )
        raise ValueError(msg)

    _validate_spacetime(config)
    _validate_fields(config)
    _validate_constants(config)
    _validate_derived_fields(config)
    _validate_background_fields(config)
    if has_lagrangian:
        _validate_lagrangian(config)
    if has_linearization:
        _validate_linearization(config, has_lagrangian=has_lagrangian)
    _validate_gauge(config)
    _validate_reduction(config)
    _validate_parameters(config)
    if "torsion" in config:
        _validate_torsion(config)


def _validate_torsion(config: dict[str, Any]) -> None:
    """Validate optional [torsion] section.

    The [torsion] section extends the spacetime connection to include torsion
    (Poincaré gauge theory).  It requires a ``perturbation_name`` for the
    linearized torsion field, which must be a valid identifier and not collide
    with reserved names or existing field/constant names.

    Raises
    ------
    ValueError
        If perturbation_name is missing, invalid, or collides with other names.
    """
    torsion = config["torsion"]
    pert_name = torsion.get("perturbation_name")
    if not pert_name:
        msg = "[torsion] requires 'perturbation_name' (e.g., perturbation_name = \"t\")"
        raise ValueError(msg)
    if not _VALID_FIELD_NAME.match(pert_name):
        msg = (
            f"[torsion].perturbation_name '{pert_name}' must be alphanumeric "
            f"starting with a letter"
        )
        raise ValueError(msg)
    if pert_name in _RESERVED_NAMES:
        msg = (
            f"[torsion].perturbation_name '{pert_name}' is reserved. "
            f"Reserved names: {sorted(_RESERVED_NAMES)}"
        )
        raise ValueError(msg)

    # Check for collision with field or constant names
    field_names = {f["name"] for f in config.get("fields", [])}
    const_names = set(config.get("constants", {}).get("names", []))
    if pert_name in field_names:
        msg = (
            f"[torsion].perturbation_name '{pert_name}' collides with "
            f"a [[fields]] entry. Use a different name."
        )
        raise ValueError(msg)
    if pert_name in const_names:
        msg = (
            f"[torsion].perturbation_name '{pert_name}' collides with "
            f"a [constants] entry. Use a different name."
        )
        raise ValueError(msg)
