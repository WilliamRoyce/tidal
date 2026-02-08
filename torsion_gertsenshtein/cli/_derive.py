"""``tg derive`` — Derive equations from Lagrangian via Wolfram/xAct."""

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


def _validate_config(config: dict[str, Any]) -> None:
    """Validate TOML config structure.

    Raises
    ------
    ValueError
        If any required section or field is missing or invalid.
    """
    required_sections = ["spacetime", "lagrangian"]
    for section in required_sections:
        if section not in config:
            msg = f"Missing required section: [{section}]"
            raise ValueError(msg)

    _validate_spacetime(config)
    _validate_fields(config)
    _validate_lagrangian(config)


# --- Code generation helpers ---


def _make_prefix(config: dict[str, Any]) -> str:
    """Generate a 2-3 letter symbol prefix from the theory name."""
    name = config.get("theory", {}).get("name", "")
    if not name:
        return "tg"
    # Take initials of first 2-3 words, lowercase
    words = name.split()
    prefix = "".join(w[0].lower() for w in words[:3])
    return prefix if len(prefix) >= _MIN_PREFIX_LEN else "tg"


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
    expression: str, fields: list[dict[str, Any]], prefix: str
) -> str:
    """Replace user field names with prefixed xAct names in the Lagrangian."""
    result = expression

    # Sort by name length descending to avoid partial replacements
    sorted_fields = sorted(fields, key=lambda f: len(f["name"]), reverse=True)

    for field in sorted_fields:
        name = field["name"]
        prefixed = f"{prefix}{name.capitalize()}"
        # Replace field name references (e.g., phi → tgPhi, but not inside other words)
        # Handle: CD[-a][phi[]] → CD[-a][tgPhi[]]
        # Handle: phi[] → tgPhi[]
        # Handle: C[-a, -b, -c] → tgC[-a, -b, -c] (but with proper prefix)
        result = result.replace(f"{name}[", f"{prefixed}[")
        result = result.replace(f"{name} ", f"{prefixed} ")

    # Also substitute eta and CD with prefixed versions
    result = result.replace("eta[", f"{prefix}Eta[")
    result = result.replace("CD[", f"{prefix}CD[")
    return result.replace("CD ]", f"{prefix}CD ]")


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


def _wls_header(ctx: _WlsContext) -> list[str]:
    """Generate script header lines."""
    return [
        "#!/usr/bin/env wolframscript",
        f"(* Auto-generated by tg derive: {ctx.theory_name} *)",
        "",
        f'Print["=== {ctx.theory_name} ==="];',
        'Print[""];',
        "",
    ]


def _wls_packages(pipeline_path: str) -> list[str]:
    """Generate xAct package loading and pipeline import lines.

    Parameters
    ----------
    pipeline_path : str
        Absolute path to the ``torsion_gertsenshtein/wolfram/`` directory.
    """
    escaped = pipeline_path.replace("\\", "\\\\").replace('"', '\\"')
    return [
        "(* Load xAct packages *)",
        "<< xAct`xTensor`;",
        "<< xAct`xCoba`;",
        "",
        "(* Load pipeline modules *)",
        f'pipelinePath = "{escaped}";',
        'Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];',
        'Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];',
        'Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];',
        'Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];',
        "",
        "$DefInfoQ = False;",
        "",
    ]


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


def _wls_fields(ctx: _WlsContext) -> list[str]:
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

    return lines


def _wls_lagrangian(ctx: _WlsContext) -> list[str]:
    """Generate Lagrangian definition lines."""
    prefixed = _substitute_field_names(ctx.lagrangian_expr, ctx.fields, ctx.prefix)
    return [
        "(* Step 3: Lagrangian *)",
        f"{ctx.prefix}Lagrangian = (",
        f"  {prefixed}",
        ");",
        "",
        f'Print["Lagrangian: ", {ctx.prefix}Lagrangian];',
        "",
    ]


def _wls_euler_lagrange_multi(ctx: _WlsContext) -> list[str]:
    """Generate Euler-Lagrange, decomposition, and export lines for multi-field."""
    lines: list[str] = ["(* Step 4: Euler-Lagrange equations *)"]

    for field in ctx.fields:
        fname = field["name"]
        fexpr = _field_expression(field, ctx.prefix)
        eom_var = f"eom{fname.capitalize()}"
        lines.extend((
            f"{eom_var} = VarD[{fexpr}, {ctx.cd}][{ctx.prefix}Lagrangian];",
            f"{eom_var} = ToCanonical[{eom_var}];",
            f"{eom_var} = ContractMetric[{eom_var}, {ctx.metric}];",
            f'Print["EOM {fname}: ", {eom_var}];',
            "",
        ))

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

        lines.extend((
            f"{comp_var} = DecomposeToComponents[{eom_var}, {fexpr}, {ctx.chart}, {{{others_str}}}];",
            f'Print["{fname} components: ", Length[{comp_var}]];',
            "",
        ))

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

    return [
        "(* Step 4: Euler-Lagrange equations *)",
        f"eom = EulerLagrangeEquation[{ctx.prefix}Lagrangian, {fexpr}, {ctx.cd}];",
        'Print["EOM: ", eom];',
        "",
        "(* Step 5: Decompose to components *)",
        f"componentEqs = DecomposeToComponents[eom, {fexpr}, {ctx.chart}];",
        'Print["Components: ", Length[componentEqs]];',
        "",
        "fieldEquations = Table[",
        f'  {{"{fname}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},',
        "  {k, Length[componentEqs]}",
        "];",
        "",
    ]


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

    escaped_lagrangian = (
        ctx.lagrangian_expr.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", " ")
        .strip()
    )

    lines: list[str] = [
        "metadata = <|",
        '  "source" -> "xAct",',
        f'  "lagrangian_expr" -> "{escaped_lagrangian}",',
        '  "derived_from" -> "Euler-Lagrange",',
        '  "gauge" -> "none",',
        '  "linearized" -> False,',
        f'  "dimension" -> {ctx.dim},',
        f'  "signature" -> {{{sig_str}}},',
        f'  "coordinates" -> {{{coord_str}}}',
        "|>;",
        "",
    ]

    # Build JSON — always use multi-field builder since fieldEquations
    # is constructed with proper labels by both single and multi-field paths
    lines.extend(("jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];", ""))

    # Export
    escaped_output = str(ctx.output_path).replace("\\", "\\\\").replace('"', '\\"')
    lines.extend((
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
    ))

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
        lagrangian_expr=config["lagrangian"]["expression"].strip(),
        is_multi=len(config["fields"]) > 1,
        pipeline_path=str(wolfram_dir),
    )

    lines: list[str] = []
    lines.extend(_wls_header(ctx))
    lines.extend(_wls_packages(ctx.pipeline_path))
    lines.extend(_wls_spacetime(config, ctx))
    lines.extend(_wls_fields(ctx))
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
        print("Or use --dry-run to see the generated script without execution.", file=sys.stderr)
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
        encoding="utf-8", mode="w", suffix=".wls", delete=False, prefix="tg_derive_"
    ) as tmp:
        tmp.write(script_content)
        tmp_path = Path(tmp.name)

    try:
        ret = _run_wolframscript(tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)

    # Post-validate output JSON if wolframscript succeeded
    if ret == 0:
        raw_output = args.output or config.get("output", {}).get(
            "path", "output.json"
        )
        resolved = Path(raw_output)
        if not resolved.is_absolute():
            resolved = (config_path.parent.resolve() / resolved).resolve()
        if resolved.exists():
            try:
                from torsion_gertsenshtein.symbolic.json_loader import (
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
