"""``tidal derive`` — Derive equations from Lagrangian via Wolfram/xAct."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

_VALID_FIELD_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")

_MIN_DIM = 2
_MAX_DIM = 7
_MIN_PREFIX_LEN = 2

if TYPE_CHECKING:
    from argparse import Namespace

# Coordinates per spacetime dimension
_COORDS: dict[int, list[str]] = {
    2: ["t", "x"],
    3: ["t", "x", "y"],
    4: ["t", "x", "y", "z"],
    5: ["t", "x", "y", "z", "w"],
    6: ["t", "x", "y", "z", "w", "v"],
    7: ["t", "x", "y", "z", "w", "v", "u"],
}

_INDEX_LETTERS = list("abcdefghijklmnop")  # 16 letters, enough for dim 7 + 4 = 11

# Minkowski signatures per spacetime dimension
_MINKOWSKI_SIGNATURES: dict[int, list[int]] = {
    2: [-1, 1],
    3: [-1, 1, 1],
    4: [-1, 1, 1, 1],
}


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


# --- Validation ---


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
    if not _VALID_FIELD_NAME.match(fname):
        msg = f"Field name '{fname}' must be alphanumeric starting with a letter"
        raise ValueError(msg)
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
    if not pf or not isinstance(pf, str):
        msg = "[linearization].perturbation_field is required"
        raise ValueError(msg)
    field_names = {f["name"] for f in config.get("fields", [])}
    if pf not in field_names:
        msg = f"[linearization].perturbation_field '{pf}' not found in [[fields]]"
        raise ValueError(msg)


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


def _validate_config(config: dict[str, Any]) -> None:
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
    _validate_derived_fields(config)
    _validate_background_fields(config)
    if has_lagrangian:
        _validate_lagrangian(config)
    if has_linearization:
        _validate_linearization(config, has_lagrangian=has_lagrangian)
    _validate_gauge(config)
    _validate_parameters(config)


# --- Code generation helpers ---


def _make_prefix(config: dict[str, Any]) -> str:
    """Generate a 2-3 letter symbol prefix from the theory name."""
    name = config.get("theory", {}).get("name", "")
    if not name:
        return "tidal"
    # Take initials of first 2-3 words (alpha only), lowercase
    words = name.split()
    initials = [w[0].lower() for w in words if w[0].isalpha()]
    prefix = "".join(initials[:3])
    return prefix if len(prefix) >= _MIN_PREFIX_LEN else "tidal"


def _generate_metric_code(config: dict[str, Any], prefix: str) -> str:
    """Generate Wolfram code for metric definition.

    Raises
    ------
    ValueError
        If metric type is unknown or dimension mismatch in diagonal/matrix.
    """
    dim = config["spacetime"]["dimension"]
    metric_type = config["spacetime"].get("metric", "minkowski")

    if metric_type == "minkowski":
        sig = _MINKOWSKI_SIGNATURES.get(dim)
        if sig is None:
            sig = [-1] + [1] * (dim - 1)
        diag_str = ", ".join(str(s) for s in sig)
        return f"{prefix}MetricMatrix = DiagonalMatrix[{{{diag_str}}}];"

    if metric_type == "diagonal":
        entries = config["spacetime"].get("diagonal", [])
        if len(entries) != dim:
            msg = f"[spacetime].diagonal must have {dim} entries, got {len(entries)}"
            raise ValueError(msg)
        # Entries can be numbers or strings (for coordinate-dependent)
        parts = [str(e) for e in entries]
        return f"{prefix}MetricMatrix = DiagonalMatrix[{{{', '.join(parts)}}}];"

    if metric_type == "matrix":
        matrix = config["spacetime"].get("matrix", [])
        if len(matrix) != dim:
            msg = f"[spacetime].matrix must be {dim}x{dim}"
            raise ValueError(msg)
        rows: list[str] = []
        for row in matrix:
            if len(row) != dim:
                msg = f"[spacetime].matrix rows must have {dim} entries"
                raise ValueError(msg)
            rows.append("{" + ", ".join(str(e) for e in row) + "}")
        return f"{prefix}MetricMatrix = {{{', '.join(rows)}}};"

    msg = f"Unknown metric type: '{metric_type}'. Use: minkowski, diagonal, matrix"
    raise ValueError(msg)


def _generate_field_def(field: dict[str, Any], prefix: str, manifold: str) -> str:
    """Generate DefTensor code for a field."""
    name = field["name"]
    ftype = field["type"]
    prefixed = f"{prefix}{name.capitalize()}"

    if ftype == "scalar":
        return f"If[!xTensorQ[{prefixed}],\n  DefTensor[{prefixed}[], {manifold}]\n];"

    if ftype == "vector":
        return f"If[!xTensorQ[{prefixed}],\n  DefTensor[{prefixed}[-a], {manifold}]\n];"

    # Tensor with rank and symmetry
    rank = field["rank"]
    symmetry = field.get("symmetry", "none")
    indices = ", ".join(f"-{_INDEX_LETTERS[i]}" for i in range(rank))

    if symmetry == "antisymmetric":
        sym_spec = f", Antisymmetric[{{{indices}}}]"
    elif symmetry == "symmetric":
        sym_spec = f", Symmetric[{{{indices}}}]"
    else:
        sym_spec = ""

    return (
        f"If[!xTensorQ[{prefixed}],\n"
        f"  DefTensor[{prefixed}[{indices}], {manifold}{sym_spec},\n"
        f'    PrintAs -> "{name}"]\n'
        f"];"
    )


def _field_expression(field: dict[str, Any], prefix: str) -> str:
    """Return the xAct expression for a field reference (e.g., 'phi[]' or 'C[-a,-b,-c]')."""
    name = field["name"]
    ftype = field["type"]
    prefixed = f"{prefix}{name.capitalize()}"

    if ftype == "scalar":
        return f"{prefixed}[]"
    if ftype == "vector":
        return f"{prefixed}[-a]"
    rank = field["rank"]
    indices = ", ".join(f"-{_INDEX_LETTERS[i]}" for i in range(rank))
    return f"{prefixed}[{indices}]"


def _substitute_field_names(
    expression: str,
    fields: list[dict[str, Any]],
    prefix: str,
    *,
    derived_fields: list[dict[str, Any]] | None = None,
    background_fields: list[dict[str, Any]] | None = None,
) -> str:
    """Replace user field names with prefixed xAct names in the Lagrangian."""
    result = expression

    # Merge fundamental, derived, and background fields for substitution
    all_fields = list(fields)
    if derived_fields:
        all_fields.extend(derived_fields)
    if background_fields:
        all_fields.extend(background_fields)

    # Sort by name length descending to avoid partial replacements
    sorted_fields = sorted(all_fields, key=lambda f: len(f["name"]), reverse=True)

    for field in sorted_fields:
        name = field["name"]
        prefixed = f"{prefix}{name.capitalize()}"
        # Replace field name references (e.g., phi → {prefix}Phi, but not inside other words)
        # Handle: CD[-a][phi[]] → CD[-a][{prefix}Phi[]]
        # Handle: phi[] → {prefix}Phi[]
        # Handle: C[-a, -b, -c] → {prefix}C[-a, -b, -c]
        result = result.replace(f"{name}[", f"{prefixed}[")
        result = result.replace(f"{name} ", f"{prefixed} ")

    # Also substitute eta, CD, and bg (background/reference metric) with prefixed versions
    result = result.replace("eta[", f"{prefix}Eta[")
    result = result.replace("bg[", f"{prefix}Bg[")
    result = result.replace("CD[", f"{prefix}CD[")
    result = result.replace("CD]", f"{prefix}CD]")
    result = result.replace("CD ]", f"{prefix}CD ]")
    # Substitute chart placeholder for component-derivative notation
    # e.g., CD[{0, -chart}][ux[]] → {prefix}CD[{0, -{prefix}Cart}][...]
    return result.replace("-chart}", f"-{prefix}Cart}}")


# --- WLS script generation ---


@dataclass(frozen=True)
class _WlsContext:
    """Shared context for WLS script generation helpers."""

    prefix: str
    dim: int
    fields: list[dict[str, Any]]
    constants: list[str]
    coords: list[str]
    manifold: str
    metric: str
    cd: str
    chart: str
    theory_name: str
    output_path: str
    lagrangian_expr: str
    is_multi: bool
    pipeline_path: str
    parameters: dict[str, float]
    derived_fields: list[dict[str, Any]]
    background_fields: list[dict[str, Any]]
    linearization: dict[str, Any] | None
    constraint_solver: dict[str, Any] | None
    gauge: list[dict[str, Any]]


def _wls_header(ctx: _WlsContext) -> list[str]:
    """Generate script header lines."""
    return [
        "#!/usr/bin/env wolframscript",
        f"(* Auto-generated by tidal derive: {ctx.theory_name} *)",
        "",
        f'Print["=== {ctx.theory_name} ==="];',
        'Print[""];',
        "",
    ]


def _wls_packages(
    pipeline_path: str,
    *,
    load_xpert: bool = False,
    load_gauge: bool = False,
) -> list[str]:
    """Generate xAct package loading and pipeline import lines.

    Parameters
    ----------
    pipeline_path : str
        Absolute path to the ``tidal/wolfram/`` directory.
    load_xpert : bool
        If *True*, also load ``xAct`xPert`` and ``Linearize.wl``.
    load_gauge : bool
        If *True*, also load ``GaugeFix.wl``.
    """
    escaped = pipeline_path.replace("\\", "\\\\").replace('"', '\\"')
    lines = [
        "(* Load xAct packages *)",
        "<< xAct`xTensor`;",
        "<< xAct`xCoba`;",
    ]
    if load_xpert:
        lines.append("<< xAct`xPert`;")
    lines.extend(
        (
            "",
            "(* Load pipeline modules *)",
            f'pipelinePath = "{escaped}";',
            'Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];',
            'Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];',
        )
    )
    if load_xpert:
        lines.append('Get[FileNameJoin[{pipelinePath, "Linearize.wl"}]];')
    if load_gauge:
        lines.append('Get[FileNameJoin[{pipelinePath, "GaugeFix.wl"}]];')
    lines.extend(("", "$DefInfoQ = False;", ""))
    return lines


def _wls_spacetime(config: dict[str, Any], ctx: _WlsContext) -> list[str]:
    """Generate spacetime definition lines (manifold, metric, chart)."""
    coord_funcs = ", ".join(f"{c}[]" for c in ctx.coords)
    indices = ", ".join(str(i) for i in range(ctx.dim))
    idx_str = ", ".join(_INDEX_LETTERS[: min(ctx.dim + 4, 8)])

    return [
        f"(* Step 1: Define {ctx.dim}D spacetime *)",
        f"If[!xTensorQ[{ctx.manifold}],",
        f"  DefManifold[{ctx.manifold}, {ctx.dim}, {{{idx_str}}}]",
        "];",
        "",
        f"If[!MetricQ[{ctx.metric}],",
        f"  DefMetric[-1, {ctx.metric}[-a, -b], {ctx.cd},",
        '    SymbolOfCovD -> {";", "\\[Del]"},',
        '    PrintAs -> "\\[Eta]"]',
        "];",
        "",
        f"If[!ChartQ[{ctx.chart}],",
        f"  DefChart[{ctx.chart}, {ctx.manifold}, {{{indices}}}, {{{coord_funcs}}}]",
        "];",
        "",
        _generate_metric_code(config, ctx.prefix),
        f"MetricInBasis[{ctx.metric}, -{ctx.chart}, {ctx.prefix}MetricMatrix];",
        "",
    ]


def _needs_bg_tensor(config: dict[str, Any]) -> bool:
    """Check if any expression references the background metric ``bg[...]``.

    The background metric is a non-dynamical reference tensor used in theories
    like massive gravity (Fierz-Pauli) and bimetric gravity.  It is defined
    via ``DefTensor`` so that xPert treats it as unperturbed.
    """
    exprs: list[str] = []
    if "lagrangian" in config:
        exprs.append(config["lagrangian"].get("expression", ""))
    if "linearization" in config:
        exprs.append(config["linearization"].get("expression", ""))
    exprs.extend(df.get("definition", "") for df in config.get("derived_fields", []))
    return any("bg[" in expr for expr in exprs)


def _wls_background_component_values(
    field: dict[str, Any],
    prefix: str,
    chart: str,
    dim: int,
) -> list[str]:
    """Generate ``ComponentValue`` lines for a single background field."""
    comps: list[int | float | str] = field.get("components", [])
    prefixed = f"{prefix}{field['name'].capitalize()}"
    ftype = field["type"]
    lines: list[str] = []

    if ftype == "scalar":
        lines.append(f"ComponentValue[{prefixed}[], {comps[0]}];")
    elif ftype == "vector":
        for idx, val in enumerate(comps):
            lines.append(f"ComponentValue[{prefixed}[{{{idx}, -{chart}}}], {val}];")
    else:
        # Tensor rank 2+: iterate over all index tuples
        rank = field.get("rank", 2)
        for flat_idx, val in enumerate(comps):
            multi_idx: list[int] = []
            remaining = flat_idx
            for _ in range(rank):
                multi_idx.append(remaining % dim)
                remaining //= dim
            multi_idx.reverse()
            idx_str = ", ".join(f"{{{k}, -{chart}}}" for k in multi_idx)
            lines.append(f"ComponentValue[{prefixed}[{idx_str}], {val}];")

    return lines


def _wls_validate_backgrounds_after_decompose(
    ctx: _WlsContext,
    comp_var: str,
) -> list[str]:
    """Generate validation that vector/tensor backgrounds resolved after ToBasis.

    After DecomposeToComponents, all background tensor symbols should be
    fully resolved to coordinate expressions via ComponentValue + ToBasis.
    This emits a check that catches silent substitution failures.

    Only emitted when there are non-scalar background fields (scalars use
    explicit ReplaceAll, so they're always resolved).
    """
    non_scalar_bgs = [f for f in ctx.background_fields if f["type"] != "scalar"]
    if not non_scalar_bgs:
        return []

    bg_heads = ", ".join(
        f"{ctx.prefix}{f['name'].capitalize()}" for f in non_scalar_bgs
    )
    return [
        "(* Validate background field substitution succeeded *)",
        f"ValidateNoUnresolvedBackgrounds[{comp_var}, {{{bg_heads}}}];",
    ]


def _wls_scalar_background_substitution(
    ctx: _WlsContext,
    eom_var: str,
) -> list[str]:
    """Generate explicit ``ReplaceAll`` for scalar background fields.

    Scalar backgrounds (rank 0) have no indices, so ``ToBasis`` inside
    ``DecomposeToComponents`` does **not** trigger ``ComponentValue``
    substitution.  We must substitute explicitly before decomposition.

    Vector/tensor backgrounds need explicit substitution AFTER
    decomposition — see ``_wls_vector_background_substitution``.
    """
    lines: list[str] = []
    for field in ctx.background_fields:
        if field["type"] != "scalar":
            continue
        prefixed = f"{ctx.prefix}{field['name'].capitalize()}"
        comps = field.get("components", [])
        if not comps:
            continue
        value = comps[0]
        lines.extend(
            (
                f"(* Substitute scalar background {field['name']} -> {value} *)",
                f"{eom_var} = {eom_var} /. {{{prefixed}[] -> {value}}};",
            )
        )
    return lines


def _wls_vector_background_substitution(
    ctx: _WlsContext,
    comp_var: str,
) -> list[str]:
    """Generate explicit ``ReplaceAll`` for vector/tensor background fields.

    Unlike scalar backgrounds (substituted BEFORE decomposition because
    ``ToBasis`` doesn't trigger ``ComponentValue`` for rank-0 tensors),
    vector/tensor backgrounds are substituted AFTER decomposition so that
    xAct handles the index algebra (contractions, metric raising/lowering)
    before we inject numeric component values.

    Both covariant (``{i, -chart}``) and contravariant (``{i, chart}``)
    index orientations get the same component value.  This is correct for
    diagonal metrics (Minkowski); non-diagonal metrics would need metric
    factors for index raising/lowering.
    """
    non_scalar_bgs = [f for f in ctx.background_fields if f["type"] != "scalar"]
    if not non_scalar_bgs:
        return []

    lines: list[str] = []
    for field in non_scalar_bgs:
        prefixed = f"{ctx.prefix}{field['name'].capitalize()}"
        comps: list[int | float | str] = field.get("components", [])
        if not comps:
            continue

        rules: list[str] = []
        if field["type"] == "vector":
            for idx, val in enumerate(comps):
                rules.extend(
                    (
                        f"{prefixed}[{{{idx}, -{ctx.chart}}}] -> {val}",
                        f"{prefixed}[{{{idx}, {ctx.chart}}}] -> {val}",
                        # Component function form (after ReplaceTensorFieldComponents
                        # converts vbdB[{i, -chart}] -> vbdBi[t, x, y])
                        f"{prefixed}{idx}[__] -> {val}",
                    )
                )
        else:
            # Tensor rank 2+: iterate over all index tuples
            rank = field.get("rank", 2)
            for flat_idx, val in enumerate(comps):
                multi_idx: list[int] = []
                remaining = flat_idx
                for _ in range(rank):
                    multi_idx.append(remaining % ctx.dim)
                    remaining //= ctx.dim
                multi_idx.reverse()
                idx_down = ", ".join(f"{{{k}, -{ctx.chart}}}" for k in multi_idx)
                idx_up = ", ".join(f"{{{k}, {ctx.chart}}}" for k in multi_idx)
                # Component function name: head + concatenated index digits
                comp_name = "".join(str(k) for k in multi_idx)
                rules.extend(
                    (
                        f"{prefixed}[{idx_down}] -> {val}",
                        f"{prefixed}[{idx_up}] -> {val}",
                        # Component function form (after ReplaceTensorFieldComponents)
                        f"{prefixed}{comp_name}[__] -> {val}",
                    )
                )

        rules_str = ", ".join(rules)
        lines.extend(
            [
                f"(* Substitute vector/tensor background {field['name']} *)",
                f"{comp_var} = {comp_var} /. {{{rules_str}}};",
            ]
        )

    return lines


def _wls_fields(ctx: _WlsContext, *, include_bg: bool = False) -> list[str]:
    """Generate field definitions and constant symbol declarations."""
    lines: list[str] = ["(* Step 2: Define fields *)"]
    for field in ctx.fields:
        lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))

    if ctx.constants:
        lines.append("(* Define constant symbols *)")
        lines.extend(
            f"If[!ConstantSymbolQ[{c}], DefConstantSymbol[{c}]];" for c in ctx.constants
        )
        lines.append("")

    # Background fields — non-dynamical tensors (not varied by VarD)
    if ctx.background_fields:
        lines.append(
            "(* Background fields — non-dynamical (not varied in Euler-Lagrange) *)"
        )
        for field in ctx.background_fields:
            lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))
            # Set component values via xCoba's ComponentValue mechanism.
            # After ToBasis in ComponentDecompose, xCoba replaces the tensor
            # with these values — no pipeline changes needed.
            lines.extend(
                _wls_background_component_values(field, ctx.prefix, ctx.chart, ctx.dim)
            )
            lines.append("")
        lines.append("")

    if include_bg:
        bg_name = f"{ctx.prefix}Bg"
        lines.extend(
            [
                "(* Background/reference metric — not perturbed by xPert *)",
                f"If[!xTensorQ[{bg_name}],",
                f'  DefTensor[{bg_name}[-a, -b], {ctx.manifold}, Symmetric[{{-a, -b}}], PrintAs -> "bg"]',
                "];",
                "(* Explicit zero perturbation: bg is non-dynamical *)",
                "Unprotect[Perturbation];",
                f"Perturbation[{bg_name}[__], ___] := 0;",
                "Protect[Perturbation];",
                "",
            ]
        )

    return lines


def _wls_derived_fields(ctx: _WlsContext) -> list[str]:
    """Generate DefTensor and MakeRule for each derived field."""
    if not ctx.derived_fields:
        return []

    lines: list[str] = ["(* Derived field definitions *)"]
    for field in ctx.derived_fields:
        # DefTensor (reuse existing helper)
        lines.extend((_generate_field_def(field, ctx.prefix, ctx.manifold), ""))

        # Build the MakeRule LHS: prefixed tensor with lowered indices
        fexpr = _field_expression(field, ctx.prefix)
        # Substitute field names in the definition
        defn = _substitute_field_names(
            field["definition"].strip(),
            ctx.fields,
            ctx.prefix,
            derived_fields=ctx.derived_fields,
            background_fields=ctx.background_fields,
        )
        rule_var = f"{ctx.prefix}{field['name'].capitalize()}Rules"
        lines.extend(
            (
                f"{rule_var} = MakeRule[{{{fexpr}, {defn}}}, MetricOn -> All, ContractMetrics -> True];",
                "",
            )
        )

    return lines


def _wls_lagrangian(ctx: _WlsContext) -> list[str]:
    """Generate Lagrangian definition lines."""
    prefixed = _substitute_field_names(
        ctx.lagrangian_expr,
        ctx.fields,
        ctx.prefix,
        derived_fields=ctx.derived_fields,
        background_fields=ctx.background_fields,
    )
    lines = [
        "(* Step 3: Lagrangian *)",
        f"{ctx.prefix}Lagrangian = (",
        f"  {prefixed}",
        ");",
        "",
    ]

    if ctx.derived_fields:
        lines.append("(* Expand derived field definitions *)")
        for field in ctx.derived_fields:
            rule_var = f"{ctx.prefix}{field['name'].capitalize()}Rules"
            lines.append(
                f"{ctx.prefix}Lagrangian = {ctx.prefix}Lagrangian /. {rule_var};"
            )
        lines.extend(
            (
                f"{ctx.prefix}Lagrangian = ToCanonical[{ctx.prefix}Lagrangian];",
                f"{ctx.prefix}Lagrangian = ContractMetric[{ctx.prefix}Lagrangian, {ctx.metric}];",
                "",
                f'Print["Lagrangian (expanded): ", {ctx.prefix}Lagrangian];',
                "",
            )
        )
    else:
        lines.extend(
            (
                f'Print["Lagrangian: ", {ctx.prefix}Lagrangian];',
                "",
            )
        )

    return lines


def _wls_linearize_from_lagrangian(
    ctx: _WlsContext,
    *,
    include_bg: bool = False,
) -> list[str]:
    """Lagrangian-first linearization: L → L^(2) → EOM + canonical from L^(2).

    Single-path approach using xPert's 2nd-order perturbation:

    1. ``Perturbation[L, 2] / 2`` → L^(2) (quadratic Lagrangian).
       xPert correctly perturbs all metric-dependent objects (Christoffels,
       Ricci tensor, etc.) *before* evaluating on the flat background, so
       L^(2) retains the full linearized Einstein tensor contribution even
       though R₀ = 0 on Minkowski.

    2. Expand ``Scalar[x]^n`` → ``∏ Scalar[RenameDummies[x]]`` so that VarD
       can vary through each copy independently (fixes the index-collision
       problem from Fierz-Pauli trace-squared terms).

    3. ``VarD[H[-a,-b], CD][L^(2)]`` → correct linearized EOM.

    4. Same L^(2) serves the canonical pipeline (momenta π, Hamiltonian H).

    Valid for flat Minkowski background where ``√(-g₀) = 1`` and ``L₀ = 0``
    (Brizuela, Martín-García, Mena Marugán 2009; Carroll 2004, Ch. 7).

    Raises
    ------
    ValueError
        If ``ctx.linearization`` is ``None``.
    """
    if ctx.linearization is None:
        msg = "_wls_linearize_from_lagrangian called without linearization config"
        raise ValueError(msg)
    lin = ctx.linearization
    pert_field_name = lin["perturbation_field"]
    pert_field = next(f for f in ctx.fields if f["name"] == pert_field_name)
    fexpr = _field_expression(pert_field, ctx.prefix)

    pert_sym = f"{ctx.prefix}hpert"
    eps_sym = f"{ctx.prefix}Epsilon"
    field_head = f"{ctx.prefix}{pert_field_name.capitalize()}"
    bg_name = f"{ctx.prefix}Bg"
    ricci_sym = f"Ricci{ctx.cd}"
    ricci_scalar_sym = f"RicciScalar{ctx.cd}"

    p = ctx.prefix

    # ------------------------------------------------------------------
    # Step 1: L^(2) = Perturbation[L, 2] / 2 via xPert
    # ------------------------------------------------------------------
    lines: list[str] = [
        "",
        "(* ============================================================ *)",
        "(* Lagrangian-first linearization (single-path via L^(2))       *)",
        "(* L -> L^(2) = Perturbation[L, 2] / 2  (xPert 2nd-order)     *)",
        "(* Then: VarD[H, CD][L^(2)] -> linearized EOM                   *)",
        "(* Same L^(2) feeds canonical pipeline (pi, H)                  *)",
        "(* ============================================================ *)",
        "",
        "(* Save original nonlinear Lagrangian *)",
        f"lOriginal = {p}Lagrangian;",
        "",
        "(* Set up metric perturbation: g = eta + epsilon * h  via xPert *)",
        f"{pert_sym}Tensor = SetupMetricPerturbation[{ctx.metric}, {pert_sym}, {eps_sym}];",
        f'Print["Perturbation: {ctx.metric} -> {ctx.metric} + {eps_sym} * {pert_field_name}"];',
        "",
        "(* 2nd-order perturbation of Lagrangian *)",
        "(* xPert perturbs curvature tensors symbolically BEFORE evaluating *)",
        "(* on flat background, so L^(2) retains linearized Einstein tensor *)",
        "l2Raw = Perturbation[lOriginal, 2];",
        "l2Raw = ExpandPerturbation[l2Raw];",
        'Print["L^(2) raw (expanded): ", Short[l2Raw, 3]];',
        "",
        "(* Validate that xPert fully expanded *)",
        "If[!FreeQ[l2Raw, Perturbation],",
        '  Throw["Linearization: ExpandPerturbation did not fully expand L^(2)."]',
        "];",
        "",
        "(* Drop 2nd-order metric perturbation h^(2) -- keep h^(1)*h^(1) only *)",
        f"l2Raw = l2Raw /. {pert_sym}[LI[2], idx__] :> 0;",
        "",
        "(* Replace xPert notation with declared field tensor *)",
        f"l2Raw = l2Raw /. {pert_sym}[LI[1], idx__] :> {field_head}[idx];",
    ]

    # Replace bg → metric if bg tensor is used
    if include_bg:
        lines.extend(
            [
                "",
                "(* Replace reference metric bg -> background metric *)",
                f"l2Raw = l2Raw /. {bg_name} -> {ctx.metric};",
            ]
        )

    # Scalar BG substitutions
    lines.extend(_wls_scalar_background_substitution(ctx, "l2Raw"))

    lines.extend(
        [
            "",
            "(* Set background curvature to zero (flat Minkowski) *)",
            f"l2Raw = l2Raw /. {{{ricci_sym}[__] :> 0, {ricci_scalar_sym}[] :> 0}};",
            "",
            "(* Canonical simplifications *)",
            "l2Raw = ToCanonical[l2Raw];",
            f"l2Raw = ContractMetric[l2Raw, {ctx.metric}];",
            "",
            "(* L^(2) = delta^2 L / 2 *)",
            f"{p}Lagrangian = l2Raw / 2;",
            f'Print["L^(2) set: ", Short[{p}Lagrangian, 5]];',
            "",
        ]
    )

    # ------------------------------------------------------------------
    # Step 1b: Type A gauge fixing (if any) — add L_gf to L^(2)
    # Perturbation-level gauge terms (current presets) are already
    # quadratic in the perturbation field, so they are added directly
    # to L^(2) without further xPert expansion.
    # ------------------------------------------------------------------
    if ctx.gauge and any(
        _resolve_gauge_mechanism(g) == "lagrangian_term" for g in ctx.gauge
    ):
        lines.extend(_wls_gauge_fixing_type_a(ctx))

    # ------------------------------------------------------------------
    # Step 2: EOM via VarD[H[-a,-b], CD][L^(2)]
    # ------------------------------------------------------------------
    lines.extend(
        [
            "(* ------------------------------------------------------------ *)",
            "(* EOM: VarD[H[-a,-b], CD][L^(2)]                              *)",
            "(* ------------------------------------------------------------ *)",
            "",
            "(* Expand Scalar[x]^n into products with renamed dummies.       *)",
            "(* This fixes VarD index collision from Fierz-Pauli (tr h)^2:   *)",
            "(*   Scalar[eta^ab H_ab]^2 -> Scalar[eta^ab H_ab]*Scalar[eta^cd H_cd] *)",
            f"l2ForVarD = {p}Lagrangian;",
            "l2ForVarD = l2ForVarD //. Scalar[x_]^n_Integer?Positive :>",
            "  Times @@ Table[Scalar[RenameDummies[x]], {n}];",
            'Print["L^(2) for VarD (Scalar expanded): ", Short[l2ForVarD, 5]];',
            "",
            "(* Vary L^(2) with respect to perturbation field H *)",
            f"eomLin = VarD[{field_head}[-a, -b], {ctx.cd}][l2ForVarD];",
            "eomLin = ToCanonical[eomLin];",
            f"eomLin = ContractMetric[eomLin, {ctx.metric}];",
            'Print["Linearized EOM: ", Short[eomLin, 5]];',
            "",
        ]
    )

    lines.extend(
        [
            "(* Decompose to components *)",
            f'componentEqs = DecomposeToComponents[eomLin, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {p}MetricMatrix];',
            'Print["Components: ", Length[componentEqs]];',
        ]
    )
    lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
    lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))

    # Build fieldEquations table
    lines.extend(
        [
            "",
            "fieldEquations = Table[",
            f'  {{"{pert_field_name}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
            "  {k, Length[componentEqs]}",
            "];",
            "",
        ]
    )

    # ------------------------------------------------------------------
    # Step 3: Type B gauge fixing (if any) — constraints on fieldEquations
    # ------------------------------------------------------------------
    if ctx.gauge and any(
        _resolve_gauge_mechanism(g) == "constraint" for g in ctx.gauge
    ):
        lines.extend(_wls_gauge_fixing_type_b(ctx))

    return lines


def _resolve_gauge_mechanism(entry: dict[str, Any]) -> str:
    """Return ``'lagrangian_term'`` or ``'constraint'`` for a gauge entry."""
    if entry["type"] == "custom":
        return str(entry["mechanism"])
    return _GAUGE_PRESETS[entry["type"]]["mechanism"]


def _wls_gauge_fixing_type_a(ctx: _WlsContext) -> list[str]:
    """Generate Type A gauge-fixing terms (added to Lagrangian before EL).

    For named presets, emits a call to the corresponding ``Build*GaugeTerm``
    function in ``GaugeFix.wl``.  For custom expressions, emits the user's
    Wolfram expression directly (after field-name substitution).
    """
    lines: list[str] = ["(* Gauge fixing: Lagrangian terms *)"]
    for entry in ctx.gauge:
        if _resolve_gauge_mechanism(entry) != "lagrangian_term":
            continue
        field_name: str = entry["field"]
        xi = entry.get("xi", 1)
        pfx_field = f"{ctx.prefix}{field_name.capitalize()}"

        if entry["type"] == "custom":
            expr = _substitute_field_names(
                entry["expression"],
                ctx.fields,
                ctx.prefix,
                derived_fields=ctx.derived_fields,
                background_fields=ctx.background_fields,
            )
            lines.append(f"{ctx.prefix}GaugeTerm = ({expr});")
        else:
            builder = _GAUGE_PRESETS[entry["type"]]["builder"]
            lines.append(
                f"{ctx.prefix}GaugeTerm = {builder}[{pfx_field}, {ctx.metric}, {ctx.cd}, {xi}];"
            )
        lines.extend(
            (
                f"{ctx.prefix}Lagrangian = AddGaugeFixingTerm["
                f"{ctx.prefix}Lagrangian, {ctx.prefix}GaugeTerm];",
                f"{ctx.prefix}Lagrangian = ToCanonical[{ctx.prefix}Lagrangian];",
                f"{ctx.prefix}Lagrangian = ContractMetric[{ctx.prefix}Lagrangian, {ctx.metric}];",
                f'Print["Gauge-fixed Lagrangian '
                f'({entry["type"]} on {field_name}): ", {ctx.prefix}Lagrangian];',
                "",
            )
        )
    return lines


def _type_b_zero_component(
    comp_name: str, field_name: str, gauge_type: str
) -> list[str]:
    """Generate WLS to substitute a component and its derivatives with zero.

    Used by temporal gauge (``A_0 = 0``) and axial gauge (``A_n = 0``).
    Applied to the ``fieldEquations`` variable after component decomposition.
    """
    return [
        f"(* {gauge_type.capitalize()} gauge: {field_name} component {comp_name} = 0 *)",
        f"fieldEquations = fieldEquations /. {{{comp_name}[args___] :> 0, Derivative[ders__][{comp_name}][args___] :> 0}};",
        f'Print["Applied {gauge_type} gauge: {comp_name} = 0"];',
        "",
    ]


def _type_b_coulomb_constraint(
    ctx: _WlsContext,
    field_name: str,
    comp_pfx: str,
) -> list[str]:
    """Generate WLS to add a Coulomb gauge constraint (spatial divergence = 0).

    Appends an additional constraint equation ``div A_spatial = 0`` to
    ``fieldEquations`` after component decomposition.
    """
    coords = _COORDS[ctx.dim]
    coord_args = ", ".join(f"{c}[]" for c in coords)

    # Build sum of spatial first derivatives: d_1 A_1 + d_2 A_2 + ...
    deriv_terms: list[str] = []
    for i in range(1, ctx.dim):
        deriv_indices = ", ".join("1" if j == i else "0" for j in range(ctx.dim))
        comp = f"{comp_pfx}{i}"
        deriv_terms.append(f"Derivative[{deriv_indices}][{comp}][{coord_args}]")

    constraint_expr = " + ".join(deriv_terms)
    constraint_name = f"{field_name}_coulomb"

    return [
        f"(* Coulomb gauge: spatial divergence of {field_name} = 0 *)",
        f'AppendTo[fieldEquations, {{"{constraint_name}", {constraint_expr}}}];',
        f'Print["Added Coulomb constraint: div {field_name} = 0"];',
        "",
    ]


def _sym_flat_index(a: int, b: int, dim: int) -> int:
    """Flat index for symmetric pair ``(a, b)`` in ``dim``-dimensional space.

    Matches the lexicographic ordering used by ``ReplaceRank2FieldComponents``
    in ``ComponentDecompose.wl``: pairs ``(i, j)`` with ``i <= j``, flattened
    to ``0, 1, 2, ...``
    """
    if a > b:
        a, b = b, a
    return a * (2 * dim - 1 - a) // 2 + b


def _tt_transverse_constraints(
    dim: int,
    comp_pfx: str,
    field_name: str,
    coord_args: str,
) -> list[str]:
    """Generate WLS for TT transverse constraints: d^i h_{i,j} = 0 per spatial j."""
    coords = _COORDS[dim]
    lines: list[str] = []
    for j in range(1, dim):
        deriv_terms: list[str] = []
        for i in range(1, dim):
            flat_idx = _sym_flat_index(i, j, dim)
            comp = f"{comp_pfx}{flat_idx}"
            deriv_indices = ", ".join("1" if k == i else "0" for k in range(dim))
            deriv_terms.append(f"Derivative[{deriv_indices}][{comp}][{coord_args}]")

        constraint_expr = " + ".join(deriv_terms)
        coord_label = coords[j] if j < len(coords) else str(j)
        constraint_name = f"{field_name}_transverse_{coord_label}"
        lines.extend(
            [
                f"(* TT transverse: d^i h_{{i,{j}}} = 0 *)",
                f'AppendTo[fieldEquations, {{"{constraint_name}", {constraint_expr}}}];',
                f'Print["Added TT transverse constraint: d^i h_{{i,{j}}} = 0"];',
                "",
            ]
        )
    return lines


def _tt_traceless_substitution(
    dim: int,
    comp_pfx: str,
    field_name: str,
    coord_args: str,
) -> list[str]:
    """Generate WLS for TT traceless substitution: h_{d-1,d-1} → -(other diags)."""
    spatial_diag_indices = [_sym_flat_index(i, i, dim) for i in range(1, dim)]
    last_diag_idx = spatial_diag_indices[-1]
    other_diag_indices = spatial_diag_indices[:-1]
    last_comp = f"{comp_pfx}{last_diag_idx}"

    repl_sum = " + ".join(f"{comp_pfx}{idx}[args]" for idx in other_diag_indices)
    deriv_repl_sum = " + ".join(
        f"Derivative[d][{comp_pfx}{idx}][args]" for idx in other_diag_indices
    )
    trace_terms = " + ".join(
        f"{comp_pfx}{idx}[{coord_args}]" for idx in spatial_diag_indices
    )

    return [
        f"(* TT traceless: substitute {last_comp} → -(sum of other spatial diags) *)",
        f"fieldEquations = fieldEquations /. {{"
        f"{last_comp}[args___] :> -({repl_sum}), "
        f"Derivative[d__][{last_comp}][args___] :> -({deriv_repl_sum})}};",
        "",
        "(* Expand all equations to simplify kinetic terms after traceless sub *)",
        "fieldEquations = Table[{fieldEquations[[k, 1]], "
        "Expand[fieldEquations[[k, 2]]]},"
        " {k, Length[fieldEquations]}];",
        "",
        f"(* Replace {last_comp} equation with algebraic traceless constraint *)",
        f'Do[If[fieldEquations[[k, 1]] === "{field_name}_{last_diag_idx}",'
        f' fieldEquations[[k]] = {{"{field_name}_{last_diag_idx}", {trace_terms}}}],'
        f" {{k, Length[fieldEquations]}}];",
        f'Print["Applied TT traceless: {last_comp} → -(spatial diag sum), '
        f'eq replaced with constraint"];',
        "",
    ]


def _type_b_tt_gauge(
    ctx: _WlsContext,
    field_name: str,
    comp_pfx: str,
) -> list[str]:
    """Generate WLS for TT (transverse-traceless) gauge on a symmetric rank-2 tensor.

    Applies three conditions (dimension-aware):

    1. **Temporal** (``h_{0,mu} = 0``): zero the first ``dim`` components
       (all pairs with one time index) via ``ReplaceAll``.
    2. **Transverse** (``partial^i h_{ij} = 0``): append ``dim - 1`` spatial
       divergence constraint equations.
    3. **Traceless** (``eta^{ij} h_{ij} = 0``): substitute last spatial diagonal
       ``h_{d-1,d-1} → -(sum of other spatial diags)`` throughout all equations,
       ``Expand`` to simplify, then replace h_{d-1,d-1}'s equation with the
       algebraic traceless constraint.

    Ordering matters: transverse constraints are appended *before* the traceless
    substitution so that ``h_{d-1,d-1}`` references inside the transverse
    constraints also get substituted.  Each equation is ``Expand``-ed after
    substitution to ensure the kinetic matrix simplifies correctly (analytically
    diagonal after TT substitution).
    """
    dim = ctx.dim
    coord_args = ", ".join(f"{c}[]" for c in _COORDS[dim])
    lines: list[str] = []

    # --- 1. Temporal: zero h_{0,mu} for mu = 0 .. dim-1 ---
    for mu in range(dim):
        idx = _sym_flat_index(0, mu, dim)
        lines.extend(
            _type_b_zero_component(f"{comp_pfx}{idx}", field_name, "TT-temporal")
        )

    # --- 2. Transverse: d^i h_{i,j} = 0 per spatial j ---
    lines.extend(_tt_transverse_constraints(dim, comp_pfx, field_name, coord_args))

    # --- 3. Traceless: h_{d-1,d-1} → -(sum of other spatial diags) ---
    lines.extend(_tt_traceless_substitution(dim, comp_pfx, field_name, coord_args))

    return lines


def _wls_gauge_fixing_type_b(ctx: _WlsContext) -> list[str]:
    """Generate Type B gauge-fixing code (constraints applied after decomposition).

    Type B constraints modify ``fieldEquations`` after component decomposition:

    - **temporal**: substitute ``A_0 → 0`` (and all derivatives) everywhere
    - **axial**: substitute last spatial component → 0
    - **coulomb**: append ``div A_spatial = 0`` constraint equation
    - **tt**: temporal zeroing + traceless substitution + transverse constraints
    """
    lines: list[str] = ["(* Gauge fixing: post-decomposition constraints *)"]
    for entry in ctx.gauge:
        if _resolve_gauge_mechanism(entry) != "constraint":
            continue
        gauge_type: str = entry["type"]
        field_name: str = entry["field"]
        comp_pfx = f"{ctx.prefix}{field_name.capitalize()}"

        if gauge_type == "temporal":
            lines.extend(_type_b_zero_component(f"{comp_pfx}0", field_name, "temporal"))
        elif gauge_type == "axial":
            last_spatial = ctx.dim - 1
            lines.extend(
                _type_b_zero_component(
                    f"{comp_pfx}{last_spatial}", field_name, "axial"
                ),
            )
        elif gauge_type == "coulomb":
            lines.extend(_type_b_coulomb_constraint(ctx, field_name, comp_pfx))
        elif gauge_type == "tt":
            lines.extend(_type_b_tt_gauge(ctx, field_name, comp_pfx))
    return lines


def _wls_euler_lagrange_multi(ctx: _WlsContext) -> list[str]:
    """Generate Euler-Lagrange, decomposition, and export lines for multi-field."""
    lines: list[str] = ["(* Step 4: Euler-Lagrange equations *)"]

    for field in ctx.fields:
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        eom_var = f"eom{fname.capitalize()}"
        lines.extend(
            (
                f"{eom_var} = VarD[{fexpr}, {ctx.cd}][{ctx.prefix}Lagrangian];",
                f"{eom_var} = ToCanonical[{eom_var}];",
                f"{eom_var} = ContractMetric[{eom_var}, {ctx.metric}];",
            )
        )
        # Scalar background fields need explicit substitution (ToBasis won't touch them)
        lines.extend(_wls_scalar_background_substitution(ctx, eom_var))
        lines.extend(
            (
                f'Print["EOM {fname}: ", {eom_var}];',
                "",
            )
        )

    # Step 5: Decompose
    lines.append("(* Step 5: Decompose to components *)")
    for i, field in enumerate(ctx.fields):
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        eom_var = f"eom{fname.capitalize()}"
        comp_var = f"comp{fname.capitalize()}"

        other_exprs = [
            _field_expression(f, ctx.prefix) for j, f in enumerate(ctx.fields) if j != i
        ]
        others_str = ", ".join(other_exprs) if other_exprs else ""

        lines.extend(
            (
                f'{comp_var} = DecomposeToComponents[{eom_var}, {fexpr}, {ctx.chart}, {{{others_str}}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
                f'Print["{fname} components: ", Length[{comp_var}]];',
            )
        )
        # Substitute vector/tensor backgrounds explicitly (ToBasis unreliable)
        lines.extend(_wls_vector_background_substitution(ctx, comp_var))
        # Validate all backgrounds resolved
        lines.extend(_wls_validate_backgrounds_after_decompose(ctx, comp_var))
        lines.append("")

    # Step 6: Export
    lines.extend(("(* Step 6: Build and export JSON *)", "fieldEquations = Flatten[{"))
    for i, field in enumerate(ctx.fields):
        fname = field["name"]
        comp_var = f"comp{fname.capitalize()}"
        comma = "," if i < len(ctx.fields) - 1 else ""
        lines.append(
            f'  Table[{{"{fname}_" <> ToString[{comp_var}[[k, 1]]], {comp_var}[[k, 2]]}}, {{k, Length[{comp_var}]}}]{comma}'
        )
    lines.extend(("}, 1];", ""))

    return lines


def _wls_euler_lagrange_single(ctx: _WlsContext) -> list[str]:
    """Generate Euler-Lagrange, decomposition, and export lines for single field."""
    field = ctx.fields[0]
    fname = field["name"]
    fexpr = _field_expression(field, ctx.prefix)

    lines = [
        "(* Step 4: Euler-Lagrange equations *)",
        f"eom = EulerLagrangeEquation[{ctx.prefix}Lagrangian, {fexpr}, {ctx.cd}];",
        'Print["EOM: ", eom];',
        "",
    ]

    # Scalar background fields need explicit substitution (ToBasis won't touch them)
    lines.extend(_wls_scalar_background_substitution(ctx, "eom"))

    lines.extend(
        [
            "(* Step 5: Decompose to components *)",
            f'componentEqs = DecomposeToComponents[eom, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
            'Print["Components: ", Length[componentEqs]];',
        ]
    )
    lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
    lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))
    lines.extend(
        [
            "",
            "fieldEquations = Table[",
            f'  {{"{fname}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
            "  {k, Length[componentEqs]}",
            "];",
            "",
        ]
    )

    return lines


def _wls_linearization(ctx: _WlsContext, *, include_bg: bool = False) -> list[str]:
    """Generate xPert linearization, decomposition, and export lines.

    Raises
    ------
    ValueError
        If ``ctx.linearization`` is ``None``.
    """
    if ctx.linearization is None:
        msg = "_wls_linearization called without linearization config"
        raise ValueError(msg)
    lin = ctx.linearization
    pert_field_name = lin["perturbation_field"]
    pert_field = next(f for f in ctx.fields if f["name"] == pert_field_name)
    fexpr = _field_expression(pert_field, ctx.prefix)

    # Auto-generated internal symbols
    pert_sym = f"{ctx.prefix}hpert"
    eps_sym = f"{ctx.prefix}Epsilon"

    # Substitute CD and field names in the expression
    prefixed_expr = _substitute_field_names(
        lin["expression"].strip(),
        ctx.fields,
        ctx.prefix,
        derived_fields=ctx.derived_fields,
        background_fields=ctx.background_fields,
    )

    bg_name = f"{ctx.prefix}Bg"

    lines: list[str] = [
        "(* Step 3: xPert metric perturbation setup *)",
        f"{pert_sym}Tensor = SetupMetricPerturbation[{ctx.metric}, {pert_sym}, {eps_sym}];",
        f'Print["Perturbation: metric = background + {eps_sym} * {pert_field_name}"];',
        "",
        "(* Step 4: Linearize tensor expression *)",
        f"linExpr = LinearizeTensorExpression[{prefixed_expr}];",
        'Print["Linearized expression: ", Short[linExpr, 3]];',
        "",
        "(* Convert xPert notation to plain tensor for pipeline *)",
        f"linExprPlain = linExpr /. {pert_sym}[LI[1], idx__] :> {ctx.prefix}{pert_field_name.capitalize()}[idx];",
    ]

    # Replace bg tensor with metric (bg = background metric by construction)
    if include_bg:
        lines.extend(
            (
                "(* Replace bg with metric — bg is the background by construction *)",
                f"linExprPlain = linExprPlain /. {bg_name} -> {ctx.metric};",
            )
        )

    lines.append("linExprPlain = Simplify[linExprPlain];")

    # Scalar background fields need explicit substitution (ToBasis won't touch them)
    lines.extend(_wls_scalar_background_substitution(ctx, "linExprPlain"))

    lines.extend(
        [
            'Print["Converted to plain tensor: ", Short[linExprPlain, 3]];',
            "",
            "(* Step 5: Decompose to components *)",
            f'componentEqs = DecomposeToComponents[linExprPlain, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
            'Print["Components: ", Length[componentEqs]];',
        ]
    )
    lines.extend(_wls_vector_background_substitution(ctx, "componentEqs"))
    lines.extend(_wls_validate_backgrounds_after_decompose(ctx, "componentEqs"))
    lines.extend(
        [
            "",
            "fieldEquations = Table[",
            f'  {{"{pert_field_name}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
            "  {k, Length[componentEqs]}",
            "];",
            "",
        ]
    )

    return lines


def _wls_constraint_metadata(
    cs_config: dict[str, Any], spatial_coords: list[str]
) -> list[str]:
    """Generate Wolfram metadata lines for constraint solver configuration.

    Translates the ``[constraint_solver]`` TOML section into metadata keys
    that ``ConstraintSolverHints`` in ExportJSON.wl reads at export time.
    """
    method = cs_config.get("method", "auto")
    max_iter = cs_config.get("max_iterations", 30)
    tol = cs_config.get("tolerance", 1e-10)

    # Format tolerance in Wolfram scientific notation: 1e-10 → 1*^-10
    tol_str = f"{tol:.0e}"
    tol_str = tol_str.replace("e+0", "*^").replace("e-0", "*^-")
    tol_str = tol_str.replace("e+", "*^").replace("e-", "*^-")

    lines = [
        '  "solve_constraints" -> True,',
        f'  "constraint_method" -> "{method}",',
        f'  "constraint_max_iterations" -> {max_iter},',
        f'  "constraint_tolerance" -> {tol_str},',
    ]

    # Build boundary conditions Association
    bcs = cs_config.get("boundary_conditions", {})
    if bcs:
        bc_parts: list[str] = []
        for coord in spatial_coords:
            bc_info = bcs.get(coord, {})
            bc_type = bc_info.get("type", "periodic")
            if "value" in bc_info:
                bc_parts.append(
                    f'    "{coord}" -> <|"type" -> "{bc_type}", '
                    f'"value" -> {bc_info["value"]}|>'
                )
            else:
                bc_parts.append(f'    "{coord}" -> <|"type" -> "{bc_type}"|>')
        bc_str = ",\n".join(bc_parts)
        lines.append(f'  "constraint_boundary_conditions" -> <|\n{bc_str}\n  |>,')

    return lines


def _canonical_field_heads(ctx: _WlsContext) -> tuple[str, str]:
    """Return (heads_str, all_heads_str) for canonical pipeline WLS generation."""
    p = ctx.prefix
    field_heads = [f"{p}{f['name'].capitalize()}" for f in ctx.fields]
    all_heads = list(field_heads)
    all_heads.extend(f"{p}{bg['name'].capitalize()}" for bg in ctx.background_fields)
    return ", ".join(field_heads), ", ".join(all_heads)


def _field_component_count(field: dict[str, Any], dim: int) -> int:
    """Return the number of independent components for *field* in *dim* spacetime."""
    ftype = field["type"]
    if ftype == "scalar":
        return 1
    if ftype == "vector":
        return dim
    rank = field.get("rank", 2)
    sym = field.get("symmetry")
    if sym == "symmetric":
        # C(dim + rank - 1, rank) = (dim+rank-1)! / (rank! * (dim-1)!)
        from math import comb

        return comb(dim + rank - 1, rank)
    if sym == "antisymmetric":
        # C(dim, rank) = dim! / (rank! * (dim-rank)!)
        from math import comb

        return comb(dim, rank)
    return dim**rank


def _wls_canonical_hamiltonian(ctx: _WlsContext, all_heads_str: str) -> list[str]:
    """Generate WLS code to compute H via component-level Legendre transform.

    Computes canonical momenta using ``D[Lcomp, velocity_i]`` at the
    component level, bypassing the abstract-index ``CanonicalMomentum``
    (which fails when the Lagrangian uses abstract metric contractions).

    Sets up WLS variables: ``lagComp``, ``piCompList``, ``canonicalH``,
    ``compToFunc``, ``velOrders``, ``coordSyms``, ``allCompNames``,
    ``hamiltonianTerms``.
    """
    p = ctx.prefix

    lines: list[str] = [
        "",
        "(* === Phase K: Canonical Momentum & Hamiltonian (component-level) === *)",
        'Print[""];',
        'Print["Computing canonical momenta and Hamiltonian..."];',
        "",
        "(* Copy Lagrangian for canonical analysis *)",
        f"lagForCanon = {p}Lagrangian;",
        "",
        "(* Re-introduce explicit metric tensors for correct sign handling.       *)",
        "(* When derived fields are expanded, ContractMetric absorbs the metric   *)",
        "(* into index positions (raised/lowered). DecomposeScalarExpression needs *)",
        "(* explicit metric tensors to correctly apply the Minkowski signature     *)",
        "(* (or any user-supplied metric) during component evaluation.             *)",
        f"lagForCanon = SeparateMetric[{ctx.metric}][lagForCanon];",
    ]

    # Apply scalar BG substitutions before decomposition
    lines.extend(_wls_scalar_background_substitution(ctx, "lagForCanon"))

    lines.extend(
        [
            "",
            "(* Decompose Lagrangian to component form *)",
            f"lagComp = DecomposeScalarExpression[lagForCanon, {ctx.chart}, {{{all_heads_str}}}, "
            f'"MetricMatrix" -> {p}MetricMatrix];',
        ]
    )

    # Apply vector BG substitutions after decomposition
    lines.extend(_wls_vector_background_substitution(ctx, "lagComp"))

    lines.extend(
        [
            'Print["L (components): ", Short[lagComp, 5]];',
            "",
            "(* Build velocity-order pattern: {1, 0, ...} for d_t *)",
            f"coordSyms = ScalarsOfChart[{ctx.chart}];",
            "nCoords = Length[coordSyms];",
            "velOrders = Table[If[i == 1, 1, 0], {i, nCoords}];",
            "",
            "(* Normalize all Derivative arities to full dimension.  *)",
            "(* ConvertCDToDerivatives may produce mixed arities     *)",
            "(* (e.g. Derivative[1,0] and Derivative[0,0,1] in 2+1D) *)",
            "(* which causes D[] to fail matching. Pad to nCoords.   *)",
            "lagComp = lagComp /. Derivative[orders__][g_][args__] /;",
            "  Length[{orders}] < nCoords :>",
            "  Derivative[Sequence @@ PadRight[{orders}, nCoords, 0]][g][args];",
            "",
            "(* Map component names to Wolfram function symbols *)",
            "compToFunc = <||>;",
        ]
    )

    # Build component-to-function mapping from Python field definitions
    for field in ctx.fields:
        fname = field["name"]
        head = f"{p}{fname.capitalize()}"
        n_comps = _field_component_count(field, ctx.dim)
        lines.extend(f'compToFunc["{fname}_{j}"] = {head}{j};' for j in range(n_comps))

    lines.extend(
        [
            "",
            "(* Compute canonical momenta: pi_i = dL/d(d_t q_i) *)",
            "allCompNames = fieldEquations[[All, 1]];",
            "piCompList = {};",
            "canonicalH = 0;",
            "Do[",
            "  Module[{compName, compFunc, vel, piComp},",
            "    compName = allCompNames[[k]];",
            "    compFunc = compToFunc[compName];",
            "    vel = Derivative[Sequence @@ velOrders][compFunc][Sequence @@ coordSyms];",
            "    piComp = D[lagComp, vel];",
            "    AppendTo[piCompList, {compName, piComp}];",
            "    canonicalH += piComp * vel;",
            '    Print["pi(", compName, "): ", piComp];',
            "  ],",
            "  {k, Length[allCompNames]}",
            "];",
            "",
            "(* Legendre transform: H = Sigma pi_i * vel_i - L *)",
            "canonicalH = Expand[canonicalH - lagComp];",
            'Print["H (components): ", Short[canonicalH, 5]];',
            "",
            "(* Parse H into structured quadratic terms *)",
            "hamiltonianTerms = ParseHamiltonianExpression[canonicalH, allCompNames];",
            'Print["Hamiltonian terms: ", Length[hamiltonianTerms]];',
            "",
        ]
    )
    return lines


def _total_raw_component_count(ctx: _WlsContext) -> int:
    """Total *raw* (pre-symmetry-reduction) component count for all fields.

    Used to decide whether `DecomposeScalarExpression` on the abstract
    Lagrangian is feasible.  For rank-3 antisymmetric tensors in 4D the
    raw count is 4^3 = 64 even though only C(4,3) = 4 are independent;
    the Lagrangian still uses the full tensor head C[-a,-b,-c] so the
    decomposition cost scales with the *raw* count.
    """
    total = 0
    for f in ctx.fields:
        ftype = f["type"]
        if ftype == "scalar":
            total += 1
        elif ftype == "vector":
            total += ctx.dim
        else:
            total += ctx.dim ** f.get("rank", 2)
    return total


# Maximum raw component count for which we attempt full Lagrangian
# decomposition.  Beyond this threshold the EOM-based fast path is used.
_LAGRANGIAN_DECOMPOSE_THRESHOLD = 30


def _wls_canonical_from_eom(ctx: _WlsContext) -> list[str]:
    """Generate canonical structure directly from the already-decomposed EOM.

    This is the **fast path** for high-rank tensor fields (rank >= 3) where
    ``DecomposeScalarExpression`` on the abstract Lagrangian is prohibitively
    slow.

    E-L velocity form: equations are preserved as-is.  Only hamiltonian_terms
    are injected (empty for fast path — H reconstruction is optional).

    ``fieldEquations`` must already exist in the WLS script context.
    """
    p = ctx.prefix
    return [
        "",
        "(* === Canonical Structure (EOM-based fast path) === *)",
        "(* Lagrangian decomposition skipped: high raw component count. *)",
        "(* E-L equations preserved as-is. hamiltonian_terms left empty. *)",
        'Print[""];',
        'Print["Building canonical structure from EOM (fast path)..."];',
        "",
        "(* Compute spatial volume element sqrt|det(g_spatial)| *)",
        f"sqrtDetGSpatial = Simplify[Sqrt[Abs[Det[{p}MetricMatrix[[2;;, 2;;]]]]]];",
        'Print["sqrt|g_spatial|: ", sqrtDetGSpatial];',
        "",
        "(* Inject canonical structure — hamiltonian_terms empty (fast path). *)",
        "canonicalSection = <|",
        '  "hamiltonian_terms" -> {}',
        "|>;",
        "If[sqrtDetGSpatial =!= 1,",
        '  canonicalSection["volume_element"] = ToString[sqrtDetGSpatial, InputForm]',
        "];",
        'jsonStructure["canonical"] = canonicalSection;',
        "",
        'Print["Canonical structure (EOM-based, hamiltonian_terms empty) injected."];',
        'Print["E-L equations preserved (no Hamilton equation injection)."];',
        'Print[""];',
        "",
    ]


def _wls_canonical_pipeline(ctx: _WlsContext) -> list[str]:
    """Generate canonical momentum + Hamiltonian computation and JSON injection.

    Two paths:

    **Full path** (component count <= threshold): Decomposes the Lagrangian
    to component form, computes momenta, Hamiltonian, and inverts K.

    **EOM-based fast path** (component count > threshold): Constructs the
    canonical structure directly from the already-decomposed EOM, assuming
    K = I.  Used for high-rank tensors (e.g. rank-3 antisymmetric) where
    ``DecomposeScalarExpression`` would be prohibitively slow.

    The Lagrangian variable ``{prefix}Lagrangian`` and ``fieldEquations`` must
    already exist in the WLS script context (set by EL/linearization steps).
    """
    raw_count = _total_raw_component_count(ctx)
    if raw_count > _LAGRANGIAN_DECOMPOSE_THRESHOLD:
        return _wls_canonical_from_eom(ctx)

    _, all_heads_str = _canonical_field_heads(ctx)

    lines: list[str] = _wls_canonical_hamiltonian(ctx, all_heads_str)

    # E-L velocity form: keep original E-L equations, only inject
    # Hamiltonian terms for energy measurement.
    p = ctx.prefix
    lines.extend(
        [
            "(* === E-L Velocity Form: Inject Canonical Structure === *)",
            "(* E-L equations are preserved as-is in equations[] array. *)",
            "(* Only hamiltonian_terms are injected for energy measurement. *)",
            "",
            "(* Compute spatial volume element sqrt|det(g_spatial)| for energy integration *)",
            f"sqrtDetGSpatial = Simplify[Sqrt[Abs[Det[{p}MetricMatrix[[2;;, 2;;]]]]]];",
            'Print["sqrt|g_spatial|: ", sqrtDetGSpatial];',
            "",
            "(* Inject canonical structure into JSON *)",
            "canonicalSection = <|",
            '  "hamiltonian_terms" -> hamiltonianTerms',
            "|>;",
            "(* Only include volume_element when non-trivial (curved coordinates) *)",
            "If[sqrtDetGSpatial =!= 1,",
            '  canonicalSection["volume_element"] = ToString[sqrtDetGSpatial, InputForm]',
            "];",
            'jsonStructure["canonical"] = canonicalSection;',
            "",
            'Print["Canonical structure (hamiltonian_terms only) injected into JSON."];',
            'Print["E-L equations preserved (no Hamilton equation injection)."];',
            'Print[""];',
            "",
        ]
    )

    return lines


def _wls_metadata_and_export(config: dict[str, Any], ctx: _WlsContext) -> list[str]:
    """Generate metadata and JSON export lines.

    Raises
    ------
    ValueError
        If dimension > 4 and no explicit signature provided.
    """
    sig_default = _MINKOWSKI_SIGNATURES.get(ctx.dim)
    if sig_default is None and "signature" not in config.get("spacetime", {}):
        msg = f"Dimension {ctx.dim}: must specify [spacetime].signature (no default for dim > 4)"
        raise ValueError(msg)
    sig_str = ", ".join(
        str(s) for s in config["spacetime"].get("signature", sig_default or [])
    )
    coord_str = ", ".join(f'"{c}"' for c in ctx.coords)

    is_linearization = ctx.linearization is not None

    if is_linearization and ctx.lagrangian_expr:
        # Lagrangian-first linearization: store the Lagrangian (not the EOM)
        raw_expr: str = ctx.lagrangian_expr
        linearized_str = "True"
    elif is_linearization and ctx.linearization is not None:
        # Legacy: direct EOM linearization
        raw_expr = str(ctx.linearization["expression"])
        linearized_str = "True"
    else:
        raw_expr = ctx.lagrangian_expr
        linearized_str = "False"

    escaped_expr = (
        raw_expr.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()
    )

    # Build gauge metadata string
    if ctx.gauge:
        gauge_parts: list[str] = []
        for g in ctx.gauge:
            if g["type"] == "custom":
                gauge_parts.append(f"custom({g['field']})")
            else:
                gauge_parts.append(f"{g['type']}({g['field']})")
        gauge_val = "+".join(gauge_parts)
    else:
        gauge_val = "none"

    lines: list[str] = [
        "metadata = <|",
        '  "source" -> "xAct",',
        f'  "lagrangian_expr" -> "{escaped_expr}",',
        '  "derived_from" -> "Euler-Lagrange",',
        f'  "gauge" -> "{gauge_val}",',
        f'  "linearized" -> {linearized_str},',
        f'  "dimension" -> {ctx.dim},',
        f'  "signature" -> {{{sig_str}}},',
        f'  "coordinates" -> {{{coord_str}}}',
    ]

    # Add constraint solver metadata if configured in TOML
    if ctx.constraint_solver is not None:
        # Spatial coordinates are all coords except "t"
        spatial_coords = [c for c in ctx.coords if c != "t"]
        cs_lines = _wls_constraint_metadata(ctx.constraint_solver, spatial_coords)
        # The last base metadata line needs a trailing comma
        lines[-1] += ","
        lines.extend(cs_lines)
        # Strip trailing comma from last line to avoid Null in Association
        lines[-1] = lines[-1].removesuffix(",")

    lines.extend(["|>;", ""])

    # Build JSON — always use multi-field builder since fieldEquations
    # is constructed with proper labels by both single and multi-field paths
    lines.extend(
        ("jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];", "")
    )

    # Inject runtime parameter defaults into JSON metadata
    if ctx.parameters:
        param_entries = ", ".join(f'"{k}" -> {v}' for k, v in ctx.parameters.items())
        lines.extend(
            (
                f'jsonStructure["metadata", "parameters"] = <|{param_entries}|>;',
                "",
            )
        )

    # Canonical momentum + Hamiltonian pipeline (Phase K)
    # Must run after jsonStructure is built (needs allCompNames from fieldEquations)
    # and before JSON export. Runs for ALL theories that have a Lagrangian
    # (including linearization, where {prefix}Lagrangian = L^(2)).
    # Only skip for legacy linearization (no Lagrangian available).
    if ctx.lagrangian_expr:
        lines.extend(_wls_canonical_pipeline(ctx))

    # Export
    escaped_output = str(ctx.output_path).replace("\\", "\\\\").replace('"', '\\"')
    lines.extend(
        (
            f'outputPath = "{escaped_output}";',
            "outputDir = DirectoryName[outputPath];",
            'If[outputDir =!= "" && !DirectoryQ[outputDir], CreateDirectory[outputDir]];',
            "",
            'Print["JSON Output:"];',
            'Print[ExportString[jsonStructure, "JSON"]];',
            "",
            'Export[outputPath, jsonStructure, "JSON"];',
            'Print[""];',
            'Print["Exported to: ", outputPath];',
            "",
            f'Print["*** {ctx.theory_name} derivation complete! ***"];',
        )
    )

    return lines


def generate_wls(
    config: dict[str, Any],
    output_override: str | None = None,
    *,
    config_dir: Path | None = None,
) -> str:
    """Generate a complete .wls script from a TOML config.

    Parameters
    ----------
    config : dict
        Parsed TOML configuration.
    output_override : str | None
        Override output JSON path.
    config_dir : Path | None
        Directory of the TOML config file.  Relative output paths are resolved
        against this directory (falls back to CWD if *None*).

    Returns
    -------
    str
        Complete Wolfram Language script.

    """
    _validate_config(config)

    prefix = _make_prefix(config)
    dim = config["spacetime"]["dimension"]

    # Resolve pipeline path to absolute so the WLS script works from any location
    wolfram_dir = Path(__file__).resolve().parent.parent / "wolfram"

    # Resolve output path to absolute — relative paths resolve against config_dir
    raw_output = output_override or config.get("output", {}).get("path", "output.json")
    resolved_output = Path(raw_output)
    if not resolved_output.is_absolute():
        base = config_dir if config_dir is not None else Path.cwd()
        resolved_output = (base / resolved_output).resolve()

    ctx = _WlsContext(
        prefix=prefix,
        dim=dim,
        fields=config["fields"],
        constants=config.get("constants", {}).get("names", []),
        coords=_COORDS[dim],
        manifold=f"{prefix}M{dim}",
        metric=f"{prefix}Eta",
        cd=f"{prefix}CD",
        chart=f"{prefix}Cart",
        theory_name=config.get("theory", {}).get("name", "Custom Theory"),
        output_path=str(resolved_output),
        lagrangian_expr=config.get("lagrangian", {}).get("expression", "").strip(),
        is_multi=len(config["fields"]) > 1,
        pipeline_path=str(wolfram_dir),
        parameters={k: float(v) for k, v in config.get("parameters", {}).items()},
        derived_fields=config.get("derived_fields", []),
        background_fields=config.get("background_fields", []),
        linearization=config.get("linearization"),
        constraint_solver=config.get("constraint_solver"),
        gauge=config.get("gauge", []),
    )

    is_linearization = ctx.linearization is not None
    has_gauge = bool(ctx.gauge)
    has_type_a = has_gauge and any(
        _resolve_gauge_mechanism(g) == "lagrangian_term" for g in ctx.gauge
    )
    has_type_b = has_gauge and any(
        _resolve_gauge_mechanism(g) == "constraint" for g in ctx.gauge
    )

    lines: list[str] = []
    lines.extend(_wls_header(ctx))
    lines.extend(
        _wls_packages(
            ctx.pipeline_path,
            load_xpert=is_linearization,
            load_gauge=has_type_a,  # Only load GaugeFix.wl for Type A
        )
    )
    lines.extend(_wls_spacetime(config, ctx))
    lines.extend(_wls_fields(ctx, include_bg=_needs_bg_tensor(config)))
    lines.extend(_wls_derived_fields(ctx))

    if is_linearization and ctx.lagrangian_expr:
        # Lagrangian-first linearization: single-path via L^(2)
        # L -> L^(2) = Perturbation[L,2]/2 -> VarD[H,CD][L^(2)] -> EOM
        # Same L^(2) feeds canonical pipeline (pi, H)
        lines.extend(_wls_lagrangian(ctx))
        lines.extend(
            _wls_linearize_from_lagrangian(
                ctx,
                include_bg=_needs_bg_tensor(config),
            )
        )
        # EOM computed inside _wls_linearize_from_lagrangian
    elif is_linearization:
        # Legacy: direct EOM linearization (deprecated — no [lagrangian])
        lines.extend(_wls_linearization(ctx, include_bg=_needs_bg_tensor(config)))
    else:
        lines.extend(_wls_lagrangian(ctx))
        # Type A gauge fixing: modify Lagrangian before EL derivation
        if has_type_a:
            lines.extend(_wls_gauge_fixing_type_a(ctx))
        if ctx.is_multi:
            lines.extend(_wls_euler_lagrange_multi(ctx))
        else:
            lines.extend(_wls_euler_lagrange_single(ctx))
        # Type B gauge fixing: constraints applied after decomposition
        if has_type_b:
            lines.extend(_wls_gauge_fixing_type_b(ctx))

    lines.extend(_wls_metadata_and_export(config, ctx))

    return "\n".join(lines) + "\n"


# --- Execution ---


def _run_wolframscript(script_path: Path) -> int:
    """Run wolframscript on a .wls file.

    Returns
    -------
    int
        Exit code from wolframscript.
    """
    if shutil.which("wolframscript") is None:
        print("Error: 'wolframscript' not found on PATH.", file=sys.stderr)
        print(file=sys.stderr)
        print("Install Wolfram Engine (free for development):", file=sys.stderr)
        print("  https://www.wolfram.com/engine/", file=sys.stderr)
        print(file=sys.stderr)
        print(
            "Or use --dry-run to see the generated script without execution.",
            file=sys.stderr,
        )
        return 1

    print(f"Running: wolframscript -file {script_path}")
    print()

    result = subprocess.run(
        ["wolframscript", "-file", str(script_path)],
        capture_output=False,
        check=False,
    )

    return result.returncode


def _derive_from_toml(config_path: Path, args: Namespace) -> int:
    """Run derivation from a TOML config file.

    Parameters
    ----------
    config_path : Path
        Path to the TOML config file.
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    with config_path.open("rb") as f:
        config = tomllib.load(f)

    script_content = generate_wls(
        config, output_override=args.output, config_dir=config_path.parent.resolve()
    )

    if args.dry_run:
        print(script_content)
        return 0

    if args.save_script:
        save_path = Path(args.save_script)
        save_path.write_text(script_content, encoding="utf-8")
        print(f"Saved script to: {save_path.resolve()}")
        if shutil.which("wolframscript") is None:
            print(
                "Note: wolframscript not found. Run the script manually when available."
            )
            return 0
        return _run_wolframscript(save_path)

    # Use temp file
    with tempfile.NamedTemporaryFile(
        encoding="utf-8", mode="w", suffix=".wls", delete=False, prefix="tidal_derive_"
    ) as tmp:
        tmp.write(script_content)
        tmp_path = Path(tmp.name)

    try:
        ret = _run_wolframscript(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)

    # Post-validate output JSON if wolframscript succeeded
    if ret == 0:
        raw_output = args.output or config.get("output", {}).get("path", "output.json")
        resolved = Path(raw_output)
        if not resolved.is_absolute():
            resolved = (config_path.parent.resolve() / resolved).resolve()
        if resolved.exists():
            try:
                from tidal.symbolic.json_loader import (
                    load_equation_system,
                )

                spec = load_equation_system(resolved)
                print()
                print(
                    f"Validation: JSON loaded successfully ({spec.n_components} components)"
                )
            except Exception as exc:  # noqa: BLE001
                print(f"\nWarning: JSON validation failed: {exc}", file=sys.stderr)
                ret = 1

    return ret


def derive_command(args: Namespace) -> int:
    """Execute the derive command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"Error: file not found: {config_path}", file=sys.stderr)
        return 1

    ext = config_path.suffix.lower()

    if ext == ".wls":
        return _run_wolframscript(config_path)

    if ext in {".toml", ".tml"}:
        return _derive_from_toml(config_path, args)

    print(
        f"Error: unsupported file extension '{ext}'. Use .toml for config or .wls for script.",
        file=sys.stderr,
    )
    return 1
