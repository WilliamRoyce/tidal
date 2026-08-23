"""``tidal inspect`` — Display equation system information from a JSON spec."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.symbolic.json_loader import ComponentEquation

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


def _show(value: object) -> str:
    """Render a matrix entry, normalizing xCoba coordinate calls to bare names."""
    return _COORD_CALL.sub(r"\1", str(value))


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
                diag_entries.append(_show(s) if s is not None else _show(matrix[i][i]))
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
                entries.append(_show(s) if s is not None else _show(matrix[i][j]))
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
        else "  (no Lagrangian expression in metadata)",
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


# Coordinate calls are stored as ``x[]`` (xCoba nullary form). Rendered as a
# bare ``x`` so the brackets do not nest inside the ``[coefficient]`` display
# convention, which would make the coefficient boundary ambiguous.
_COORD_CALL = re.compile(r"\b([A-Za-z]\w*)\[\]")

# Proven-sign markers. Spelled as inequalities rather than ``+``/``-`` so they
# cannot be confused with the ``+`` that separates terms, and prefixed ``eff``
# because they describe the coefficient *after* division by the LHS kinetic
# coefficient -- not the numerator printed in the bracket. Those differ in sign
# whenever the kinetic coefficient is negative: ``h_0`` in coupled_scalars shows
# ``[-kappa^(-2)]`` yet is ``eff>0``, because the LHS carries ``-kappa^(-2)``
# too and the two cancel.
_SIGN_MARKER = {
    "positive": "eff>0",
    "negative": "eff<0",
    "nonnegative": "eff>=0",
    "nonpositive": "eff<=0",
    "zero": "eff=0",
}


def _render_equation(
    equation: ComponentEquation,
    *,
    summary: bool = False,
    promoted: frozenset[str] = frozenset(),
) -> list[str]:
    """Render one equation as lines, kinetic coefficient on the left-hand side.

    The printed line *is* the equation: ``[K] d2_t(phi) = [c] op(field) + ...``.
    Carrying ``K`` on the LHS is what makes the output complete — without it a
    reader cannot tell that a bare ``-1`` coefficient in one component and a
    bare ``+1`` in another may describe the same physics, which is the
    misreading recorded as row 1 of GH #401.

    Each ``(operator, field)`` key appears exactly once with its terms summed,
    via :func:`tidal.symbolic.spec_query.effective_coefficient`, so no key can
    be read while a duplicate goes unnoticed. The individual contributions stay
    visible inside the sum, preserving base-versus-correction provenance.

    Two annotations appear only when they carry information:

    * ``eps:0,1`` when a key spans perturbative orders — merging would
      otherwise hide the split ``filter_by_order`` depends on.
    * ``eff>0`` / ``eff<0`` / ``eff>=0`` when the sign is *proven*. This is the
      sign of the coefficient **after** dividing by the LHS kinetic
      coefficient, which is the physically meaningful quantity and is not
      generally the sign of the bracket beside it: a negative kinetic
      coefficient flips the two apart. Most coefficients are free parameters
      whose sign is genuinely unknowable, and those stay unmarked rather than
      guessed at.

    Parameters
    ----------
    equation : ComponentEquation
        The equation to render.
    summary : bool
        Emit one compact line per key with the sign only, omitting coefficients.

    Returns
    -------
    list[str]
        Rendered lines, without trailing newlines.
    """
    from tidal.symbolic.spec_query import effective_coefficient

    order = equation.time_derivative_order
    name = equation.field_name
    kinetic = equation.kinetic_coefficient_symbolic
    if order:
        head = f"d{order}_t({name})"
    elif name in promoted:
        # Algebraic LHS, but the row carries (or is targeted by) inter-
        # constraint time derivatives: it belongs to the second-order
        # sector and cannot be Schur-eliminated as algebraic (GH #457).
        head = f"{name} (algebraic LHS — promoted to second-order sector, GH #457)"
    else:
        head = f"{name} (constraint)"
    kinetic_shown = _COORD_CALL.sub(r"\1", kinetic) if kinetic else None
    lhs = f"[{kinetic_shown}] {head}" if kinetic_shown else head

    keys = sorted({(t.operator, t.field) for t in equation.rhs_terms})
    if summary:
        out: list[str] = []
        for operator, field in keys:
            eff = effective_coefficient(equation, field, operator)
            marker = _SIGN_MARKER.get(eff.sign().sign.value, "?")
            out.append(f"  {name} {operator}({field}) {marker}")
        return out

    lines = [f"  {lhs} ="]
    for operator, field in keys:
        eff = effective_coefficient(equation, field, operator)
        coeff = _COORD_CALL.sub(r"\1", eff.numerator)
        tags: list[str] = []
        orders = {t.order_in_eps for t in eff.terms}
        if len(orders) > 1:
            tags.append("eps:" + ",".join(str(o) for o in sorted(orders)))
        marker = _SIGN_MARKER.get(eff.sign().sign.value)
        if marker is not None:
            tags.append(marker)
        suffix = "   " + "  ".join(tags) if tags else ""
        lines.append(f"    + [{coeff}] {operator}({field}){suffix}")
    return lines


def _print_equations(spec: object, *, detail: str = "full") -> None:
    """Print the field summary and the equations.

    Parameters
    ----------
    spec : EquationSystem
        The loaded equation system.
    detail : str
        ``"full"`` for the equation form, ``"summary"`` for one line per key
        carrying only the proven sign.
    """
    from tidal.symbolic.json_loader import EquationSystem

    if not isinstance(spec, EquationSystem):
        return

    promoted = spec.second_order_sector.promoted
    if detail != "summary":
        print(
            f"Fields ({spec.n_components} "
            f"component{'s' if spec.n_components != 1 else ''}):",
        )
        for eq in spec.equations:
            if eq.time_derivative_order > 0:
                kind = "dynamical"
            elif eq.field_name in promoted:
                kind = "promoted"  # second-order sector, GH #457
            else:
                kind = "constraint"
            print(
                f"  {eq.field_name:<12s} {kind:<12s} time_order={eq.time_derivative_order}"
            )
        print()

    print(
        "Equations:" if detail != "summary" else "Signs (proven only; ? = undecided):"
    )
    for eq in spec.equations:
        for line in _render_equation(
            eq, summary=detail == "summary", promoted=promoted
        ):
            print(line)
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
            },
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


def inspect_command(args: Namespace) -> int:  # noqa: C901, PLR0911, PLR0912
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
    from tidal.symbolic.json_loader import load_equation_system

    json_path = Path(args.json_path)
    if not json_path.exists():
        from tidal.cli._console import error_with_hint

        error_with_hint(
            f"file not found: {json_path}",
            hints=["Use `tidal list` to find specs, or `tidal derive` to generate"],
        )
        return 1

    # `tidal inspect` is read-only — never evolves the system — so the v6
    # higher-time-order guard is relaxed to a warning when ``--latex`` (or
    # other read-only output flags) is requested. This lets us render
    # equations for theories whose [perturbation] block was intentionally
    # disabled to extract the exact (non-LPS) EOMs.
    # The semantic query flags (#401) are read-only in the same sense as
    # --latex, so they relax the guard too: analysis never evolves the system.
    query_flags = ("equation", "coefficient", "families", "diff")
    is_query = any(getattr(args, name, None) for name in query_flags)
    strict_v6 = (
        not bool(getattr(args, "latex", False))
        and not bool(getattr(args, "json_output", False))
        and not is_query
    )
    spec = load_equation_system(json_path, strict_v6=strict_v6)

    if is_query:
        from tidal.cli._inspect_query import run_query

        result = run_query(args, spec)
        if result is not None:
            return result

    if args.json_output:
        data = _build_json_output(spec, show_params=args.params)
        print(json.dumps(data, indent=2))
        return 0

    if args.latex:
        from tidal.symbolic.latex import (
            kinetic_matrix_to_latex,
            load_symbol_overrides,
            system_to_latex,
        )

        if getattr(args, "symbols", None) is not None:
            load_symbol_overrides(args.symbols)

        if args.latex_format == "kinetic-matrix":
            from tidal.symbolic.kinetic_matrix import build_kinetic_matrix

            # system_to_latex sets the metric symbol as a side effect;
            # call it just to prime the module state, but we discard
            # the output and emit only the kinetic matrix.
            _ = system_to_latex(spec, output_format="raw")
            print(kinetic_matrix_to_latex(build_kinetic_matrix(spec), spec))
            return 0

        print(system_to_latex(spec, output_format=args.latex_format))
        return 0

    detail = getattr(args, "detail", "full")
    if detail != "summary":
        _print_header(spec)
        _print_spacetime(spec)
    _print_equations(spec, detail=detail)

    if detail == "summary":
        # summary exists to be cheap to scan: the parameter list and the
        # n x n matrices are the two largest remaining blocks and neither
        # carries a sign verdict, so they are omitted rather than doubling
        # the cost of the mode.
        return 0

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
        f"Mass matrix:    {_format_matrix(spec.mass_matrix, spec.mass_matrix_symbolic or None)}",
    )
    print(
        f"Coupling matrix: {_format_matrix(spec.coupling_matrix, spec.coupling_matrix_symbolic or None)}",
    )

    return 0
