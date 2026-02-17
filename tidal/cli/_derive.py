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
        If expression is missing or empty.
    """
    expr = config["lagrangian"].get("expression")
    if not expr or not expr.strip():
        msg = "[lagrangian].expression must be a non-empty string"
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
            msg = f"Background field '{fname}' conflicts with a [[derived_fields]] entry"
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


def _validate_linearization(config: dict[str, Any]) -> None:
    """Validate optional ``[linearization]`` section.

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
    if has_lagrangian and has_linearization:
        msg = "[lagrangian] and [linearization] are mutually exclusive"
        raise ValueError(msg)

    _validate_spacetime(config)
    _validate_fields(config)
    _validate_derived_fields(config)
    _validate_background_fields(config)
    if has_lagrangian:
        _validate_lagrangian(config)
    else:
        _validate_linearization(config)
    _validate_parameters(config)


# --- Code generation helpers ---


def _make_prefix(config: dict[str, Any]) -> str:
    """Generate a 2-3 letter symbol prefix from the theory name."""
    name = config.get("theory", {}).get("name", "")
    if not name:
        return "tidal"
    # Take initials of first 2-3 words, lowercase
    words = name.split()
    prefix = "".join(w[0].lower() for w in words[:3])
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


def _wls_packages(pipeline_path: str, *, load_xpert: bool = False) -> list[str]:
    """Generate xAct package loading and pipeline import lines.

    Parameters
    ----------
    pipeline_path : str
        Absolute path to the ``tidal/wolfram/`` directory.
    load_xpert : bool
        If *True*, also load ``xAct`xPert`` and ``Linearize.wl``.
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
    exprs.extend(df.get("rule", "") for df in config.get("derived_fields", []))
    return any("bg[" in expr for expr in exprs)


def _wls_background_component_values(
    field: dict[str, Any], prefix: str, chart: str, dim: int,
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
            lines.append(
                f"ComponentValue[{prefixed}[{{{idx}, -{chart}}}], {val}];"
            )
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


def _wls_scalar_background_substitution(
    ctx: _WlsContext, eom_var: str,
) -> list[str]:
    """Generate explicit ``ReplaceAll`` for scalar background fields.

    Scalar backgrounds (rank 0) have no indices, so ``ToBasis`` inside
    ``DecomposeToComponents`` does **not** trigger ``ComponentValue``
    substitution.  We must substitute explicitly before decomposition.

    Vector/tensor backgrounds are handled correctly by ``ToBasis``, so
    this only targets ``type == "scalar"``.
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
        lines.extend((
            f"(* Substitute scalar background {field['name']} -> {value} *)",
            f"{eom_var} = {eom_var} /. {{{prefixed}[] -> {value}}};",
        ))
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
        lines.append("(* Background fields — non-dynamical (not varied in Euler-Lagrange) *)")
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
        lines.extend([
            "(* Background/reference metric — not perturbed by xPert *)",
            f"If[!xTensorQ[{bg_name}],",
            f'  DefTensor[{bg_name}[-a, -b], {ctx.manifold}, Symmetric[{{-a, -b}}], PrintAs -> "bg"]',
            "];",
            "(* Explicit zero perturbation: bg is non-dynamical *)",
            "Unprotect[Perturbation];",
            f"Perturbation[{bg_name}[__], ___] := 0;",
            "Protect[Perturbation];",
            "",
        ])

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
                "",
            )
        )

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

    lines.extend([
        "(* Step 5: Decompose to components *)",
        f'componentEqs = DecomposeToComponents[eom, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
        'Print["Components: ", Length[componentEqs]];',
        "",
        "fieldEquations = Table[",
        f'  {{"{fname}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
        "  {k, Length[componentEqs]}",
        "];",
        "",
    ])

    return lines


def _wls_linearization(
    ctx: _WlsContext, *, include_bg: bool = False
) -> list[str]:
    """Generate xPert linearization, decomposition, and export lines."""
    assert ctx.linearization is not None
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
        lines.extend((
            "(* Replace bg with metric — bg is the background by construction *)",
            f"linExprPlain = linExprPlain /. {bg_name} -> {ctx.metric};",
        ))

    lines.append("linExprPlain = Simplify[linExprPlain];")

    # Scalar background fields need explicit substitution (ToBasis won't touch them)
    lines.extend(_wls_scalar_background_substitution(ctx, "linExprPlain"))

    lines.extend([
        'Print["Converted to plain tensor: ", Short[linExprPlain, 3]];',
        "",
        "(* Step 5: Decompose to components *)",
        f'componentEqs = DecomposeToComponents[linExprPlain, {fexpr}, {ctx.chart}, {{}}, "MetricMatrix" -> {ctx.prefix}MetricMatrix];',
        'Print["Components: ", Length[componentEqs]];',
        "",
        "fieldEquations = Table[",
        f'  {{"{pert_field_name}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
        "  {k, Length[componentEqs]}",
        "];",
        "",
    ])

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

    if is_linearization and ctx.linearization is not None:
        raw_expr: str = str(ctx.linearization["expression"])
        linearized_str = "True"
    else:
        raw_expr = ctx.lagrangian_expr
        linearized_str = "False"

    escaped_expr = (
        raw_expr.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()
    )

    lines: list[str] = [
        "metadata = <|",
        '  "source" -> "xAct",',
        f'  "lagrangian_expr" -> "{escaped_expr}",',
        '  "derived_from" -> "Euler-Lagrange",',
        '  "gauge" -> "none",',
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
    )

    is_linearization = ctx.linearization is not None

    lines: list[str] = []
    lines.extend(_wls_header(ctx))
    lines.extend(_wls_packages(ctx.pipeline_path, load_xpert=is_linearization))
    lines.extend(_wls_spacetime(config, ctx))
    lines.extend(_wls_fields(ctx, include_bg=_needs_bg_tensor(config)))
    lines.extend(_wls_derived_fields(ctx))

    if is_linearization:
        lines.extend(_wls_linearization(ctx, include_bg=_needs_bg_tensor(config)))
    else:
        lines.extend(_wls_lagrangian(ctx))
        if ctx.is_multi:
            lines.extend(_wls_euler_lagrange_multi(ctx))
        else:
            lines.extend(_wls_euler_lagrange_single(ctx))

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
        print(f"Saved script to: {save_path}")
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
