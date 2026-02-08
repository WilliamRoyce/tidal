"""``tg derive`` — Derive equations from Lagrangian via Wolfram/xAct."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import tomllib
from pathlib import Path
from typing import TYPE_CHECKING, Any

_VALID_FIELD_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9]*$")

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


def _validate_config(config: dict[str, Any]) -> None:
    """Validate TOML config structure."""
    required_sections = ["spacetime", "lagrangian"]
    for section in required_sections:
        if section not in config:
            msg = f"Missing required section: [{section}]"
            raise ValueError(msg)

    dim = config["spacetime"].get("dimension")
    if dim is None:
        msg = "[spacetime] must specify 'dimension'"
        raise ValueError(msg)
    if not isinstance(dim, int) or dim < 2 or dim > 7:
        msg = f"[spacetime].dimension must be integer 2-7, got: {dim}"
        raise ValueError(msg)

    if "fields" not in config or not config["fields"]:
        msg = "Must define at least one [[fields]] entry"
        raise ValueError(msg)

    for i, field in enumerate(config["fields"]):
        if "name" not in field:
            msg = f"[[fields]] entry {i} missing 'name'"
            raise ValueError(msg)
        fname = field["name"]
        if not _VALID_FIELD_NAME.match(fname):
            msg = f"Field name '{fname}' must be alphanumeric starting with a letter"
            raise ValueError(msg)
        if "type" not in field:
            msg = f"[[fields]] entry {i} ('{fname}') missing 'type'"
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

    expr = config["lagrangian"].get("expression")
    if not expr or not expr.strip():
        msg = "[lagrangian].expression must be a non-empty string"
        raise ValueError(msg)


def _make_prefix(config: dict[str, Any]) -> str:
    """Generate a 2-3 letter symbol prefix from the theory name."""
    name = config.get("theory", {}).get("name", "")
    if not name:
        return "tg"
    # Take initials of first 2-3 words, lowercase
    words = name.split()
    prefix = "".join(w[0].lower() for w in words[:3])
    return prefix if len(prefix) >= 2 else "tg"


def _generate_metric_code(config: dict[str, Any], prefix: str) -> str:
    """Generate Wolfram code for metric definition."""
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
        rows = []
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


def generate_wls(config: dict[str, Any], output_override: str | None = None) -> str:
    """Generate a complete .wls script from a TOML config.

    Parameters
    ----------
    config : dict
        Parsed TOML configuration.
    output_override : str | None
        Override output JSON path.

    Returns
    -------
    str
        Complete Wolfram Language script.
    """
    _validate_config(config)

    prefix = _make_prefix(config)
    dim = config["spacetime"]["dimension"]
    fields = config["fields"]
    constants = config.get("constants", {}).get("names", [])
    lagrangian_expr = config["lagrangian"]["expression"].strip()
    theory_name = config.get("theory", {}).get("name", "Custom Theory")

    output_path = output_override or config.get("output", {}).get("path", "output.json")

    coords = _COORDS[dim]
    coord_funcs = ", ".join(f"{c}[]" for c in coords)
    indices = ", ".join(str(i) for i in range(dim))
    manifold = f"{prefix}M{dim}"
    metric = f"{prefix}Eta"
    cd = f"{prefix}CD"
    chart = f"{prefix}Cart"

    lines: list[str] = []

    # Header
    lines.append("#!/usr/bin/env wolframscript")
    lines.append(f"(* Auto-generated by tg derive: {theory_name} *)")
    lines.append("")

    lines.append(f'Print["=== {theory_name} ==="];')
    lines.append('Print[""];')
    lines.append("")

    # Load packages
    lines.append("(* Load xAct packages *)")
    lines.append("<< xAct`xTensor`;")
    lines.append("<< xAct`xCoba`;")
    lines.append("")

    # Load pipeline
    lines.append("(* Load pipeline modules *)")
    lines.append(
        'pipelinePath = FileNameJoin[{DirectoryName[$InputFileName], "..", "..", "torsion_gertsenshtein", "wolfram"}];'
    )
    lines.append('Get[FileNameJoin[{pipelinePath, "CommonUtilities.wl"}]];')
    lines.append('Get[FileNameJoin[{pipelinePath, "EulerLagrange.wl"}]];')
    lines.append('Get[FileNameJoin[{pipelinePath, "ComponentDecompose.wl"}]];')
    lines.append('Get[FileNameJoin[{pipelinePath, "ExportJSON.wl"}]];')
    lines.append("")
    lines.append("$DefInfoQ = False;")
    lines.append("")

    # Step 1: Spacetime
    lines.append(f"(* Step 1: Define {dim}D spacetime *)")
    idx_str = ", ".join(_INDEX_LETTERS[: min(dim + 4, 8)])
    lines.extend(
        (
            f"If[!xTensorQ[{manifold}],",
            f"  DefManifold[{manifold}, {dim}, {{{idx_str}}}]",
            "];",
            "",
            f"If[!MetricQ[{metric}],",
            f"  DefMetric[-1, {metric}[-a, -b], {cd},",
            '    SymbolOfCovD -> {";", "\\[Del]"},',
            '    PrintAs -> "\\[Eta]"]',
            "];",
            "",
            f"If[!ChartQ[{chart}],",
            f"  DefChart[{chart}, {manifold}, {{{indices}}}, {{{coord_funcs}}}]",
            "];",
            "",
        )
    )

    metric_code = _generate_metric_code(config, prefix)
    lines.append(metric_code)
    lines.append(f"MetricInBasis[{metric}, -{chart}, {prefix}MetricMatrix];")
    lines.append("")

    # Step 2: Fields
    lines.append("(* Step 2: Define fields *)")
    for field in fields:
        lines.extend((_generate_field_def(field, prefix, manifold), ""))

    # Constants
    if constants:
        lines.append("(* Define constant symbols *)")
        lines.extend(
            f"If[!ConstantSymbolQ[{c}], DefConstantSymbol[{c}]];" for c in constants
        )
        lines.append("")

    # Step 3: Lagrangian
    prefixed_lagrangian = _substitute_field_names(lagrangian_expr, fields, prefix)
    lines.append("(* Step 3: Lagrangian *)")
    lines.append(f"{prefix}Lagrangian = (")
    lines.append(f"  {prefixed_lagrangian}")
    lines.append(");")
    lines.append("")
    lines.append(f'Print["Lagrangian: ", {prefix}Lagrangian];')
    lines.append("")

    # Step 4: Euler-Lagrange
    lines.append("(* Step 4: Euler-Lagrange equations *)")

    is_multi = len(fields) > 1

    if is_multi:
        # Multi-field: VarD for each, decompose with cross-field refs
        for i, field in enumerate(fields):
            fname = field["name"]
            fexpr = _field_expression(field, prefix)
            eom_var = f"eom{fname.capitalize()}"

            lines.extend(
                (
                    f"{eom_var} = VarD[{fexpr}, {cd}][{prefix}Lagrangian];",
                    f"{eom_var} = ToCanonical[{eom_var}];",
                    f"{eom_var} = ContractMetric[{eom_var}, {metric}];",
                    f'Print["EOM {fname}: ", {eom_var}];',
                    "",
                )
            )

        # Step 5: Decompose
        lines.append("(* Step 5: Decompose to components *)")
        for i, field in enumerate(fields):
            fname = field["name"]
            fexpr = _field_expression(field, prefix)
            eom_var = f"eom{fname.capitalize()}"
            comp_var = f"comp{fname.capitalize()}"

            # Additional fields = all other fields
            other_exprs = [
                _field_expression(f, prefix) for j, f in enumerate(fields) if j != i
            ]
            others_str = ", ".join(other_exprs) if other_exprs else ""

            lines.extend(
                (
                    f"{comp_var} = DecomposeToComponents[{eom_var}, {fexpr}, {chart}, {{{others_str}}}];",
                    f'Print["{fname} components: ", Length[{comp_var}]];',
                    "",
                )
            )

        # Step 6: Export
        lines.extend(
            ("(* Step 6: Build and export JSON *)", "fieldEquations = Flatten[{")
        )
        for i, field in enumerate(fields):
            fname = field["name"]
            comp_var = f"comp{fname.capitalize()}"
            comma = "," if i < len(fields) - 1 else ""
            lines.append(
                f'  Table[{{"{fname}_" <> ToString[{comp_var}[[k, 1]]], {comp_var}[[k, 2]]}}, {{k, Length[{comp_var}]}}]{comma}'
            )
        lines.extend(("}, 1];", ""))

    else:
        # Single field
        field = fields[0]
        fname = field["name"]
        fexpr = _field_expression(field, prefix)

        lines.append(f"eom = EulerLagrangeEquation[{prefix}Lagrangian, {fexpr}, {cd}];")
        lines.append('Print["EOM: ", eom];')
        lines.append("")

        # Step 5: Decompose
        lines.append("(* Step 5: Decompose to components *)")
        lines.append(f"componentEqs = DecomposeToComponents[eom, {fexpr}, {chart}];")
        lines.append('Print["Components: ", Length[componentEqs]];')
        lines.append("")

        # Format for export
        lines.append("fieldEquations = Table[")
        lines.append(
            f'  {{"{fname}_" <> ToString[componentEqs[[k, 1]]], componentEqs[[k, 2]]}},'
        )
        lines.append("  {k, Length[componentEqs]}")
        lines.append("];")
        lines.append("")

    # Metadata
    sig_default = _MINKOWSKI_SIGNATURES.get(dim)
    if sig_default is None and "signature" not in config.get("spacetime", {}):
        msg = f"Dimension {dim}: must specify [spacetime].signature (no default for dim > 4)"
        raise ValueError(msg)
    sig_str = ", ".join(
        str(s) for s in config["spacetime"].get("signature", sig_default or [])
    )
    coord_str = ", ".join(f'"{c}"' for c in coords)

    escaped_lagrangian = (
        lagrangian_expr.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()
    )

    lines.extend(
        (
            "metadata = <|",
            '  "source" -> "xAct",',
            f'  "lagrangian_expr" -> "{escaped_lagrangian}",',
            '  "derived_from" -> "Euler-Lagrange",',
            '  "gauge" -> "none",',
            '  "linearized" -> False,',
            f'  "dimension" -> {dim},',
            f'  "signature" -> {{{sig_str}}},',
            f'  "coordinates" -> {{{coord_str}}}',
            "|>;",
            "",
        )
    )

    # Build JSON
    if is_multi or any(f["type"] != "scalar" or f.get("rank", 0) > 0 for f in fields):
        lines.append(
            "jsonStructure = BuildMultiFieldJSONStructure[fieldEquations, metadata];"
        )
    else:
        lines.append("jsonStructure = BuildJSONStructure[componentEqs, metadata];")
    lines.append("")

    # Export
    escaped_output = str(output_path).replace("\\", "\\\\").replace('"', '\\"')
    lines.append(f'outputPath = "{escaped_output}";')
    lines.append("outputDir = DirectoryName[outputPath];")
    lines.append(
        'If[outputDir =!= "" && !DirectoryQ[outputDir], CreateDirectory[outputDir]];'
    )
    lines.append("")
    lines.append('Print["JSON Output:"];')
    lines.append('Print[ExportString[jsonStructure, "JSON"]];')
    lines.append("")
    lines.append('Export[outputPath, jsonStructure, "JSON"];')
    lines.append('Print[""];')
    lines.append('Print["Exported to: ", outputPath];')
    lines.append("")
    lines.append(f'Print["*** {theory_name} derivation complete! ***"];')

    return "\n".join(lines) + "\n"


def _run_wolframscript(script_path: Path) -> int:
    """Run wolframscript on a .wls file.

    Returns
    -------
    int
        Exit code from wolframscript.
    """
    if shutil.which("wolframscript") is None:
        print("Error: 'wolframscript' not found on PATH.")
        print()
        print("Install Wolfram Engine (free for development):")
        print("  https://www.wolfram.com/engine/")
        print()
        print("Or use --dry-run to see the generated script without execution.")
        return 1

    print(f"Running: wolframscript -file {script_path}")
    print()

    result = subprocess.run(
        ["wolframscript", "-file", str(script_path)],
        capture_output=False,
        check=False,
    )

    return result.returncode


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
        print(f"Error: file not found: {config_path}")
        return 1

    # Detect mode from extension
    ext = config_path.suffix.lower()

    if ext == ".wls":
        # Mode B: Script pass-through
        return _run_wolframscript(config_path)

    if ext in {".toml", ".tml"}:
        # Mode A: TOML config → generate .wls → run
        with config_path.open("rb") as f:
            config = tomllib.load(f)

        script_content = generate_wls(config, output_override=args.output)

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
            output_path = args.output or config.get("output", {}).get(
                "path", "output.json"
            )
            if Path(output_path).exists():
                try:
                    from torsion_gertsenshtein.symbolic.json_loader import (
                        load_equation_system,
                    )

                    spec = load_equation_system(output_path)
                    print()
                    print(
                        f"Validation: JSON loaded successfully ({spec.n_components} components)"
                    )
                except Exception as exc:  # noqa: BLE001
                    print(f"\nWarning: JSON validation failed: {exc}")
                    ret = 1

        return ret

    print(
        f"Error: unsupported file extension '{ext}'. Use .toml for config or .wls for script."
    )
    return 1
