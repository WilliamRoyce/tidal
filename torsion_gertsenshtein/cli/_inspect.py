"""``tg inspect`` — Display equation system information from a JSON spec."""

from __future__ import annotations

import re
from pathlib import Path
from typing import TYPE_CHECKING

from torsion_gertsenshtein.symbolic.json_loader import load_equation_system

if TYPE_CHECKING:
    from argparse import Namespace


def _discover_parameters(spec: object) -> dict[str, list[str]]:
    """Scan all terms for symbolic coefficients and map parameter → field names.

    Returns
    -------
    dict[str, list[str]]
        Mapping from parameter name to list of field names that reference it.
    """
    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        msg = f"Expected EquationSystem, got {type(spec).__name__}"
        raise TypeError(msg)
    param_map: dict[str, list[str]] = {}
    # Match bare identifiers (possibly preceded by '-'), ignoring math functions/operators
    ident_re = re.compile(r"[A-Za-z_]\w*")

    # Known math functions/constants to exclude
    math_names = {
        "E",
        "Pi",
        "Sin",
        "Cos",
        "Tan",
        "Exp",
        "Log",
        "Sqrt",
        "Power",
        "Abs",
        "Sign",
        "Floor",
        "Ceiling",
        "Round",
        "Mod",
        "sin",
        "cos",
        "tan",
        "exp",
        "log",
        "sqrt",
        "abs",
        "t",
        "x",
        "y",
        "z",
        "w",
        "v",
        "u",  # coordinates
    }

    for eq in spec.equations:
        for term in eq.rhs_terms:
            sym = term.coefficient_symbolic
            if sym is None:
                continue
            # Extract identifiers from the symbolic expression
            for match in ident_re.findall(sym):
                if match in math_names:
                    continue
                # Skip field names (they appear in field references, not as params)
                if match in spec.component_names:
                    continue
                if match.startswith("pi_"):
                    continue
                param_map.setdefault(match, [])
                if eq.field_name not in param_map[match]:
                    param_map[match].append(eq.field_name)

    return param_map


def _format_matrix(
    matrix: tuple[tuple[float, ...], ...],
    symbolic: tuple[tuple[object, ...], ...] | None = None,
) -> str:
    """Format a matrix for display, preferring symbolic form if available."""
    n = len(matrix)
    if n == 0:
        return "  (empty)"

    # Check if symbolic is available and non-trivial
    use_symbolic = symbolic and len(symbolic) == n

    # Check if diagonal
    is_diag = all(
        (matrix[i][j] == 0.0 if i != j else True) for i in range(n) for j in range(n)
    )
    # Check if zero
    is_zero = all(matrix[i][j] == 0.0 for i in range(n) for j in range(n))

    if is_zero:
        return f"zeros({n}x{n})"

    if is_diag:
        if use_symbolic:
            diag_entries = []
            for i in range(n):
                s = symbolic[i][i]  # type: ignore[index]
                diag_entries.append(str(s) if s is not None else str(matrix[i][i]))
            return f"diag({', '.join(diag_entries)})"
        diag_vals = [str(matrix[i][i]) for i in range(n)]
        return f"diag({', '.join(diag_vals)})"

    # Full matrix
    rows = []
    for i in range(n):
        if use_symbolic:
            entries = []
            for j in range(n):
                s = symbolic[i][j]  # type: ignore[index]
                entries.append(str(s) if s is not None else str(matrix[i][j]))
            rows.append("  [" + ", ".join(entries) + "]")
        else:
            rows.append("  [" + ", ".join(str(matrix[i][j]) for j in range(n)) + "]")
    return "\n".join(rows)


def inspect_command(args: Namespace) -> int:
    """Execute the inspect command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    json_path = Path(args.json_path)
    if not json_path.exists():
        print(f"Error: file not found: {json_path}")
        return 1

    spec = load_equation_system(json_path)

    # Header
    lagrangian = spec.metadata.get("lagrangian_expr", "")
    source = spec.metadata.get("source", "unknown")
    derived_from = spec.metadata.get("derived_from", "")

    print(
        f"  Lagrangian: {lagrangian}"
        if lagrangian
        else "  (no Lagrangian expression in metadata)"
    )
    source_parts = [s for s in [source, derived_from] if s]
    if source_parts:
        print(f"  Source: {', '.join(source_parts)}")
    print()

    # Spacetime
    dim_label = f"{spec.spatial_dimension}+1D"
    print("Spacetime:")
    print(f"  Dimension: {spec.dimension} ({dim_label})")
    print(f"  Coordinates: {spec.effective_coordinates}")
    if spec.metadata.get("signature"):
        sig = spec.metadata["signature"]
        sig_str = ", ".join(f"{'+' if s > 0 else ''}{s}" for s in sig)
        print(f"  Signature: ({sig_str})")
    print()

    # Fields
    print(
        f"Fields ({spec.n_components} component{'s' if spec.n_components != 1 else ''}):"
    )
    for eq in spec.equations:
        t_order = eq.time_derivative_order
        dynamical = "dynamical" if t_order > 0 else "constraint"
        print(f"  {eq.field_name:<12s} {dynamical:<12s} time_order={t_order}")
    print()

    # Equations
    print("Equations:")
    for eq in spec.equations:
        t_order = eq.time_derivative_order
        lhs_label = f"d{t_order}_t" if t_order > 0 else "constraint"
        terms_strs = []
        for term in eq.rhs_terms:
            coeff = term.coefficient
            sym = term.coefficient_symbolic
            coeff_str = f"[{sym}]" if sym is not None else f"{coeff:+.4g}"
            terms_strs.append(f"{coeff_str} {term.operator}({term.field})")
        print(f"  {lhs_label}({eq.field_name}) = {' '.join(terms_strs)}")
    print()

    # Parameters
    param_map = _discover_parameters(spec)
    if param_map:
        print("Required parameters:")
        for param, fields in sorted(param_map.items()):
            print(f"  {param}  (in: {', '.join(fields)})")
    else:
        print("Required parameters: none")

    if args.params and spec.metadata.get("parameters"):
        print()
        print("Default parameter values (from metadata):")
        for k, v in spec.metadata["parameters"].items():
            print(f"  {k} = {v}")
    print()

    # Matrices
    print(
        f"Mass matrix:    {_format_matrix(spec.mass_matrix, spec.mass_matrix_symbolic or None)}"
    )
    print(
        f"Coupling matrix: {_format_matrix(spec.coupling_matrix, spec.coupling_matrix_symbolic or None)}"
    )

    return 0
