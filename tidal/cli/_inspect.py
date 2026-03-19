"""``tidal inspect`` — Display equation system information from a JSON spec."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Namespace

# Known math functions/constants to exclude from parameter discovery
# (Mathematica + Python names)
_MATH_NAMES = {
    # Constants
    "E",
    "Pi",
    "I",
    "Infinity",
    # Basic trig (Mathematica)
    "Sin",
    "Cos",
    "Tan",
    "Cot",
    "Sec",
    "Csc",
    # Inverse trig (Mathematica)
    "ArcSin",
    "ArcCos",
    "ArcTan",
    "ArcCot",
    "ArcSec",
    "ArcCsc",
    # Hyperbolic
    "Sinh",
    "Cosh",
    "Tanh",
    "Coth",
    "Sech",
    "Csch",
    "ArcSinh",
    "ArcCosh",
    "ArcTanh",
    # Exponential / logarithm
    "Exp",
    "Log",
    "Log2",
    "Log10",
    # Algebraic
    "Sqrt",
    "Power",
    "Abs",
    "Sign",
    "Re",
    "Im",
    "Conjugate",
    # Rounding / integer
    "Floor",
    "Ceiling",
    "Round",
    "Mod",
    "Quotient",
    "IntegerPart",
    "FractionalPart",
    # Calculus / special
    "Derivative",
    "D",
    "Integrate",
    "Sum",
    "Product",
    "Max",
    "Min",
    # Python names for the same
    "sin",
    "cos",
    "tan",
    "cot",
    "sec",
    "csc",
    "exp",
    "log",
    "sqrt",
    "abs",
    "np",
    # Coordinate variable names
    "t",
    "x",
    "y",
    "z",
    "w",
    "v",
    "u",
    "r",
    "theta",
    "phi",
}


def discover_parameters(spec: object) -> dict[str, list[str]]:
    """Scan all terms for symbolic coefficients and map parameter → field names.

    Returns
    -------
    dict[str, list[str]]
        Mapping from parameter name to list of field names that reference it.

    Raises
    ------
    TypeError
        If *spec* is not an ``EquationSystem`` instance.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        msg = f"Expected EquationSystem, got {type(spec).__name__}"
        raise TypeError(msg)
    param_map: dict[str, list[str]] = {}
    # Match bare identifiers (possibly preceded by '-'), ignoring math functions/operators
    ident_re = re.compile(r"[A-Za-z_]\w*")

    for eq in spec.equations:
        for term in eq.rhs_terms:
            sym = term.coefficient_symbolic
            if sym is None:
                continue
            # Extract identifiers from the symbolic expression
            for match in ident_re.findall(sym):
                if match in _MATH_NAMES:
                    continue
                # Skip field names (they appear in field references, not as params)
                if match in spec.component_names:
                    continue
                if match.startswith("v_"):
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
        if use_symbolic and symbolic is not None:
            diag_entries: list[str] = []
            for i in range(n):
                s = symbolic[i][i]
                diag_entries.append(str(s) if s is not None else str(matrix[i][i]))
            return f"diag({', '.join(diag_entries)})"
        diag_vals = [str(matrix[i][i]) for i in range(n)]
        return f"diag({', '.join(diag_vals)})"

    # Full matrix
    rows: list[str] = []
    for i in range(n):
        if use_symbolic and symbolic is not None:
            entries: list[str] = []
            for j in range(n):
                s = symbolic[i][j]
                entries.append(str(s) if s is not None else str(matrix[i][j]))
            rows.append("  [" + ", ".join(entries) + "]")
        else:
            rows.append("  [" + ", ".join(str(matrix[i][j]) for j in range(n)) + "]")
    return "\n".join(rows)


def _print_header(spec: object) -> None:
    """Print Lagrangian and source information."""
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return
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


def _print_spacetime(spec: object) -> None:
    """Print spacetime dimension and coordinate information."""
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return
    dim_label = f"{spec.spatial_dimension}+1D"
    print("Spacetime:")
    print(f"  Dimension: {spec.dimension} ({dim_label})")
    print(f"  Coordinates: {spec.effective_coordinates}")
    if spec.metadata.get("signature"):
        sig = spec.metadata["signature"]
        sig_str = ", ".join(f"{'+' if s > 0 else ''}{s}" for s in sig)
        print(f"  Signature: ({sig_str})")
    print()


def _print_equations(spec: object) -> None:
    """Print field summary and detailed equation terms."""
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return
    print(
        f"Fields ({spec.n_components} component{'s' if spec.n_components != 1 else ''}):"
    )
    for eq in spec.equations:
        t_order = eq.time_derivative_order
        dynamical = "dynamical" if t_order > 0 else "constraint"
        print(f"  {eq.field_name:<12s} {dynamical:<12s} time_order={t_order}")
    print()

    print("Equations:")
    for eq in spec.equations:
        t_order = eq.time_derivative_order
        lhs_label = f"d{t_order}_t" if t_order > 0 else "constraint"
        terms_strs: list[str] = []
        for term in eq.rhs_terms:
            sym = term.coefficient_symbolic
            coeff_str = f"[{sym}]" if sym is not None else f"{term.coefficient:+.4g}"
            terms_strs.append(f"{coeff_str} {term.operator}({term.field})")
        print(f"  {lhs_label}({eq.field_name}) = {' '.join(terms_strs)}")
    print()


def _build_json_output(spec: object, *, show_params: bool) -> dict[str, Any]:
    """Build a JSON-serializable dict from an EquationSystem.

    Parameters
    ----------
    spec : EquationSystem
        The loaded equation system.
    show_params : bool
        Whether to include default parameter values from metadata.

    Returns
    -------
    dict[str, Any]
        JSON-serializable representation of the equation system.

    Raises
    ------
    TypeError
        If *spec* is not an ``EquationSystem`` instance.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        msg = f"Expected EquationSystem, got {type(spec).__name__}"
        raise TypeError(msg)

    equations: list[dict[str, Any]] = []
    for eq in spec.equations:
        terms = [
            {
                "coefficient": term.coefficient,
                "coefficient_symbolic": term.coefficient_symbolic,
                "operator": term.operator,
                "field": term.field,
            }
            for term in eq.rhs_terms
        ]
        equations.append(
            {
                "field_name": eq.field_name,
                "time_derivative_order": eq.time_derivative_order,
                "terms": terms,
            }
        )

    result: dict[str, Any] = {
        "spacetime": {
            "dimension": spec.dimension,
            "spatial_dimension": spec.spatial_dimension,
            "coordinates": list(spec.effective_coordinates),
        },
        "fields": list(spec.component_names),
        "equations": equations,
        "parameters": {
            "required": dict(sorted(discover_parameters(spec).items())),
        },
        "mass_matrix": [list(row) for row in spec.mass_matrix],
        "coupling_matrix": [list(row) for row in spec.coupling_matrix],
        "mass_matrix_symbolic": (
            [list(row) for row in spec.mass_matrix_symbolic]
            if spec.mass_matrix_symbolic
            else None
        ),
        "coupling_matrix_symbolic": (
            [list(row) for row in spec.coupling_matrix_symbolic]
            if spec.coupling_matrix_symbolic
            else None
        ),
    }

    if spec.metadata.get("signature"):
        result["spacetime"]["signature"] = list(spec.metadata["signature"])
    if spec.metadata.get("lagrangian_expr"):
        result["lagrangian"] = spec.metadata["lagrangian_expr"]
    if show_params and spec.metadata.get("parameters"):
        result["parameters"]["defaults"] = dict(spec.metadata["parameters"])

    return result


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
    from tidal.cli._console import error as _cerror
    from tidal.symbolic.json_loader import load_equation_system

    json_path = Path(args.json_path)
    if not json_path.exists():
        _cerror(f"file not found: {json_path}")
        return 1

    spec = load_equation_system(json_path)

    if args.json_output:
        data = _build_json_output(spec, show_params=args.params)
        print(json.dumps(data, indent=2))
        return 0

    if args.latex:
        from tidal.symbolic.latex import system_to_latex

        print(system_to_latex(spec, output_format=args.latex_format))
        return 0

    _print_header(spec)
    _print_spacetime(spec)
    _print_equations(spec)

    # Parameters
    param_map = discover_parameters(spec)
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
