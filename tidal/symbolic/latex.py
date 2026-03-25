r"""Convert JSON equation specifications to LaTeX math notation.

Provides functions to render TIDAL equation systems as publication-ready
LaTeX, including:

- Component PDEs with proper operator notation
- Lagrangian expressions with tensor index notation (``\\tensor{}`` package)
- Hamiltonian density terms
- Symbolic coefficients (Mathematica InputForm → LaTeX)

Primary public entry point:

- ``system_to_latex(spec, ...)`` — full equation system
"""

from __future__ import annotations

import re
from fractions import Fraction
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        HamiltonianTerm,
    )

# ---------------------------------------------------------------------------
# Greek letter mapping (lowercase + uppercase used in physics)
# ---------------------------------------------------------------------------

_GREEK_MAP: dict[str, str] = {
    "alpha": r"\alpha",
    "beta": r"\beta",
    "gamma": r"\gamma",
    "delta": r"\delta",
    "epsilon": r"\epsilon",
    "zeta": r"\zeta",
    "eta": r"\eta",
    "theta": r"\theta",
    "iota": r"\iota",
    "kappa": r"\kappa",
    "lambda": r"\lambda",
    "lam": r"\lambda",
    "mu": r"\mu",
    "nu": r"\nu",
    "xi": r"\xi",
    "pi": r"\pi",
    "rho": r"\rho",
    "sigma": r"\sigma",
    "tau": r"\tau",
    "upsilon": r"\upsilon",
    "phi": r"\phi",
    "chi": r"\chi",
    "psi": r"\psi",
    "omega": r"\omega",
    "Gamma": r"\Gamma",
    "Delta": r"\Delta",
    "Theta": r"\Theta",
    "Lambda": r"\Lambda",
    "Sigma": r"\Sigma",
    "Phi": r"\Phi",
    "Omega": r"\Omega",
}

# Pre-compiled patterns for Greek substitution (longest match first to avoid
# partial replacement, e.g. "epsilon" before "eps")
_GREEK_RE = re.compile(
    r"\b(" + "|".join(sorted(_GREEK_MAP, key=len, reverse=True)) + r")\b"
)

# ---------------------------------------------------------------------------
# Mathematica function → LaTeX mapping
# (mirrors _FUNCTION_MAP in _eval_utils.py; Sqrt/Abs/Exp handled separately)
# ---------------------------------------------------------------------------

_MATH_FUNC_LATEX: dict[str, str] = {
    # Trig
    "Sin": r"\sin",
    "Cos": r"\cos",
    "Tan": r"\tan",
    "Cot": r"\cot",
    "Sec": r"\sec",
    "Csc": r"\csc",
    # Inverse trig
    "ArcSin": r"\arcsin",
    "ArcCos": r"\arccos",
    "ArcTan": r"\arctan",
    # Hyperbolic
    "Sinh": r"\sinh",
    "Cosh": r"\cosh",
    "Tanh": r"\tanh",
    # Inverse hyperbolic
    "ArcSinh": r"\operatorname{arcsinh}",
    "ArcCosh": r"\operatorname{arccosh}",
    "ArcTanh": r"\operatorname{arctanh}",
    # Logarithmic
    "Log": r"\ln",
    # Special
    "Erf": r"\operatorname{erf}",
    "Sign": r"\operatorname{sgn}",
    "UnitStep": r"\Theta",
    "HeavisideTheta": r"\Theta",
}

# Regex: FuncName[expr] → \funcname(expr)   (sorted longest-first)
_RE_MATH_FUNC = re.compile(
    r"\b("
    + "|".join(sorted(_MATH_FUNC_LATEX, key=len, reverse=True))
    + r")\[([^\[\]]+)\]"
)

# Sqrt[expr] → \sqrt{expr}  (separate: uses braces not parens)
_RE_SQRT = re.compile(r"\bSqrt\[([^\[\]]+)\]")

# Abs[expr] → \left| expr \right|
_RE_ABS = re.compile(r"\bAbs\[([^\[\]]+)\]")


def _convert_math_functions(s: str) -> str:
    """Convert Mathematica math functions to LaTeX equivalents."""
    # Sqrt → \sqrt{} (must come before general func replacement)
    s = _RE_SQRT.sub(r"\\sqrt{\1}", s)
    # Abs → |...|
    s = _RE_ABS.sub(r"\\left| \1 \\right|", s)
    # General functions: FuncName[expr] → \funcname(expr)
    return _RE_MATH_FUNC.sub(
        lambda m: rf"{_MATH_FUNC_LATEX[m.group(1)]}({m.group(2)})", s
    )


# ---------------------------------------------------------------------------
# Operator → LaTeX mapping
# ---------------------------------------------------------------------------

_OPERATOR_LATEX: dict[str, str] = {
    "identity": "",
    "laplacian": r"\nabla^2",
    "laplacian_x": r"\partial_x^2",
    "laplacian_y": r"\partial_y^2",
    "laplacian_z": r"\partial_z^2",
    "gradient_x": r"\partial_x",
    "gradient_y": r"\partial_y",
    "gradient_z": r"\partial_z",
    "cross_derivative_xy": r"\partial_x \partial_y",
    "cross_derivative_xz": r"\partial_x \partial_z",
    "cross_derivative_yz": r"\partial_y \partial_z",
    "first_derivative_t": r"\partial_t",
    "biharmonic": r"\nabla^4",
    "time_derivative": "dot",  # sentinel for Hamiltonian factors
    # Ostrogradsky higher time derivatives
    "d2_t": r"\partial_t^2",
    "d3_t": r"\partial_t^{3}",
    "d4_t": r"\partial_t^{4}",
}

# Dynamic operator patterns (from json_loader.py)
_RE_SINGLE_AXIS = re.compile(r"^derivative_(\d+)_([xyzwvu])$")
_RE_MULTI_AXIS = re.compile(r"^derivative_((?:\d+[xyzwvu])+)$")
# New format from rhs.py: mixed_T{n}_S{m}{axis}
_RE_MIXED_NEW = re.compile(r"^mixed_T(\d+)_S(\d+)([xyzwvu])$")
# Old format from energy.py: mixed_{T}_{S1}_{S2}_... (positional spatial orders)
_RE_MIXED_OLD = re.compile(r"^mixed_(\d+(?:_\d+)+)$")

# ---------------------------------------------------------------------------
# Coefficient rendering helpers
# ---------------------------------------------------------------------------

_RE_E_POWER = re.compile(r"\bE\^")
_RE_RATIONAL = re.compile(r"Rational\[([^,\]]+),\s*([^,\]]+)\]")
_RE_PI = re.compile(r"\bPi\b")
_RE_SQRT = re.compile(r"Sqrt\[([^\[\]]+)\]")
_RE_COORD_CALL = re.compile(r"\b([xyzwvut])\s*\[\s*\]")
_RE_DIGIT_SUFFIX = re.compile(r"^([A-Za-z]+?)(\d+)$")
_RE_POWER_PAREN = re.compile(r"\^\(([^)]+)\)")
_RE_POWER_SIMPLE = re.compile(r"\^(\d+)")


def _find_top_level_division(s: str) -> int | None:
    """Find the index of the top-level '/' in *s*, respecting parentheses.

    Returns None if no top-level division is found.
    """
    depth = 0
    for i, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "/" and depth == 0 and i > 0:
            return i
    return None


def _strip_outer_parens(s: str) -> str:
    """Remove one layer of matching outer parentheses if present."""
    s = s.strip()
    if s.startswith("(") and s.endswith(")"):
        inner = s[1:-1]
        # Verify the parens are actually matching (not "(a)*(b)")
        depth = 0
        for ch in inner:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
            if depth < 0:
                return s  # Not matching outer parens
        if depth == 0:
            return inner
    return s


def coefficient_to_latex(expr: str) -> str:
    r"""Convert a Mathematica-style symbolic coefficient to LaTeX.

    Examples
    --------
    >>> coefficient_to_latex("-(B0^2*kappa^2)")
    '-B_0^{2} \\\\kappa^{2}'
    >>> coefficient_to_latex("1/2")
    '\\\\frac{1}{2}'
    """
    if not expr:
        return ""

    s = expr.strip()

    # Step 1: Rational[p, q] → \frac{p}{q}
    s = _RE_RATIONAL.sub(
        lambda m: rf"\frac{{{m.group(1).strip()}}}{{{m.group(2).strip()}}}", s
    )

    # Step 2: Sqrt[expr] → \sqrt{expr}
    s = _RE_SQRT.sub(lambda m: rf"\sqrt{{{m.group(1)}}}", s)

    # Step 3: E^ → e^
    s = _RE_E_POWER.sub(r"e^", s)

    # Step 4: Pi → \pi
    s = _RE_PI.sub(r"\\pi", s)

    # Step 5: Coordinate calls x[] → x
    s = _RE_COORD_CALL.sub(r"\1", s)

    # Step 5b: Math functions: Tanh[x] → \tanh(x), Abs[x] → |x|, etc.
    s = _RE_ABS.sub(r"\\left| \1 \\right|", s)
    s = _RE_MATH_FUNC.sub(lambda m: rf"{_MATH_FUNC_LATEX[m.group(1)]}({m.group(2)})", s)

    # Step 6: Numeric prefix fraction: -1/2*(rest) or 1/2*(rest) → -\frac{1}{2}(rest)
    prefix_frac = re.match(r"^(-?)(\d+)/(\d+)\*(.+)$", s)
    if prefix_frac:
        sign = prefix_frac.group(1)
        num = prefix_frac.group(2)
        den = prefix_frac.group(3)
        rest = prefix_frac.group(4)
        rest_tex = _coefficient_inner(rest)
        return rf"{sign}\frac{{{num}}}{{{den}}} {rest_tex}"

    # Step 7: Top-level fraction detection A/B → \frac{A}{B}
    div_idx = _find_top_level_division(s)
    if div_idx is not None:
        numer = s[:div_idx].strip()
        denom = s[div_idx + 1 :].strip()
        # Handle sign: -(A)/B → -\frac{A}{B}
        sign = ""
        if numer.startswith("-") and numer[1:].strip().startswith("("):
            sign = "-"
            numer = numer[1:].strip()
        numer = _strip_outer_parens(numer)
        denom = _strip_outer_parens(denom)
        numer_tex = _coefficient_inner(numer)
        denom_tex = _coefficient_inner(denom)
        return rf"{sign}\frac{{{numer_tex}}}{{{denom_tex}}}"

    return _coefficient_inner(s)


def _coefficient_inner(s: str) -> str:
    """Apply Greek, subscript, and power transforms to an expression fragment."""
    # Strip outer parens for cleaner output
    s = _strip_outer_parens(s)

    # Convert Mathematica functions before any other processing
    s = _convert_math_functions(s)

    # Greek prefix extraction: omegaP2 → omega P2 (before Greek substitution)
    for greek in sorted(_GREEK_MAP, key=len, reverse=True):
        s = re.sub(rf"\b{greek}([A-Z])", rf"{greek} \1", s)

    # Greek letter substitution
    s = _GREEK_RE.sub(lambda m: _GREEK_MAP[m.group(1)], s)

    # Parameter name subscripts: B0 → B_0, mA2 → m_{A2}
    # Process word tokens individually
    def _subscript_token(m: re.Match[str]) -> str:
        token = m.group(0)
        # Skip if it's a Greek command
        if token.startswith("\\"):
            return token
        dm = _RE_DIGIT_SUFFIX.match(token)
        if dm:
            base, digits = dm.group(1), dm.group(2)
            # Apply Greek to base if applicable
            base_tex = _GREEK_RE.sub(lambda gm: _GREEK_MAP[gm.group(1)], base)
            return rf"{base_tex}_{{{digits}}}"
        return token

    s = re.sub(r"\\?[A-Za-z]+\d*", _subscript_token, s)

    # Powers with parentheses: ^(expr) → ^{expr}
    s = _RE_POWER_PAREN.sub(lambda m: f"^{{{m.group(1)}}}", s)

    # Simple powers: ^N (multi-digit) → ^{N}
    s = _RE_POWER_SIMPLE.sub(
        lambda m: f"^{{{m.group(1)}}}" if len(m.group(1)) > 1 else f"^{m.group(1)}", s
    )

    # Multiplication: * → \, (thin space, implicit multiplication)
    return s.replace("*", r" \, ")


# ---------------------------------------------------------------------------
# Field name rendering
# ---------------------------------------------------------------------------

_FIELD_GREEK: frozenset[str] = frozenset(
    {"phi", "chi", "psi", "alpha", "beta", "gamma", "sigma", "omega", "theta"}
)


def _calligraphic_head(head_tex: str) -> str:
    r"""Wrap non-Greek Roman field heads in calligraphic font.

    Greek commands (starting with ``\\``) are left unchanged.
    Roman single letters are promoted to uppercase and wrapped
    in ``\\mathcal{}``.
    """
    if head_tex.startswith("\\"):
        return head_tex  # Already a LaTeX command (Greek)
    return rf"\mathcal{{{head_tex.upper()}}}"


def field_to_latex(
    name: str,
    *,
    tensor_meta: dict[str, list[int] | int | str] | None = None,
    coordinates: tuple[str, ...] = (),
) -> str:
    """Convert a field component name to LaTeX.

    Parameters
    ----------
    name : str
        Field name (e.g., "h_5", "phi_0", "v_phi_0").
    tensor_meta : dict, optional
        Tensor metadata from enriched JSON: ``{"tensor_head": "h",
        "tensor_rank": 2, "tensor_indices": [2, 2]}``.
    coordinates : tuple[str, ...], optional
        Coordinate names for resolving index labels (e.g., ("t", "x", "y", "z")).

    Returns
    -------
    str
        LaTeX string for the field.
    """
    # Velocity prefix: v_phi_0 → \dot{base}
    if name.startswith("v_"):
        base = field_to_latex(
            name[2:], tensor_meta=tensor_meta, coordinates=coordinates
        )
        return rf"\dot{{{base}}}"

    # With tensor metadata: use proper index labels
    if tensor_meta is not None:
        head = str(tensor_meta.get("tensor_head", name.split("_", maxsplit=1)[0]))
        raw_rank = tensor_meta.get("tensor_rank", 0)
        rank = int(raw_rank) if isinstance(raw_rank, (int, str)) else 0
        raw_indices = tensor_meta.get("tensor_indices", [])
        indices: list[int] = list(raw_indices) if isinstance(raw_indices, list) else []

        # Render head with Greek if applicable
        head_tex = _GREEK_MAP.get(head, head)

        # Non-Greek Roman heads: calligraphic font (uppercase promoted)
        head_tex = _calligraphic_head(head_tex)

        if rank == 0 or not indices:
            return head_tex

        # Map numeric indices to coordinate labels
        if coordinates and len(coordinates) > max(indices):
            idx_labels = "".join(coordinates[i] for i in indices)
        else:
            idx_labels = "".join(str(i) for i in indices)

        return rf"{head_tex}_{{{idx_labels}}}"

    # Without tensor metadata: best-effort from name pattern
    parts = name.split("_", 1)
    base = parts[0]
    base_tex = _GREEK_MAP.get(base, base)

    # Non-Greek Roman heads: calligraphic font
    base_tex = _calligraphic_head(base_tex)

    if len(parts) > 1:
        return rf"{base_tex}_{{{parts[1]}}}"
    return base_tex


# ---------------------------------------------------------------------------
# Operator rendering
# ---------------------------------------------------------------------------


def _partial_order(axis: str, order: str) -> str:
    r"""Render \\partial_axis^order, omitting ^1."""
    if order == "1":
        return rf"\partial_{axis}"
    return rf"\partial_{axis}^{{{order}}}"


def _operator_dynamic(operator: str, field_latex: str) -> str | None:
    """Try to match dynamic operator patterns. Returns None if no match."""
    # derivative_N_x → \partial_x^N
    m = _RE_SINGLE_AXIS.match(operator)
    if m:
        return _partial_order(m.group(2), m.group(1)) + f" {field_latex}"

    # derivative_2x_1y → \partial_x^2 \partial_y
    m = _RE_MULTI_AXIS.match(operator)
    if m:
        parts = [
            _partial_order(dm.group(2), dm.group(1))
            for dm in re.finditer(r"(\d+)([xyzwvu])", m.group(1))
        ]
        return " ".join(parts) + f" {field_latex}"

    # New format: mixed_T{n}_S{m}{axis} → \partial_t^n \partial_{axis}^m
    m = _RE_MIXED_NEW.match(operator)
    if m:
        t_order, s_order, axis = m.group(1), m.group(2), m.group(3)
        parts_list: list[str] = []
        if t_order != "0":
            parts_list.append(_partial_order("t", t_order))
        if s_order != "0":
            parts_list.append(_partial_order(axis, s_order))
        return " ".join(parts_list) + f" {field_latex}"

    # Old format: mixed_{T}_{S1}_{S2}_{S3} → \partial_t^T \partial_x^S1 ...
    m = _RE_MIXED_OLD.match(operator)
    if m:
        nums = operator.split("_")[1:]  # strip "mixed" prefix
        t_order = nums[0]
        parts_list = []
        if t_order != "0":
            parts_list.append(_partial_order("t", t_order))
        axes = "xyzwvu"
        for i, s_order in enumerate(nums[1:]):
            if s_order != "0" and i < len(axes):
                parts_list.append(_partial_order(axes[i], s_order))
        return " ".join(parts_list) + f" {field_latex}"

    return None


def operator_to_latex(operator: str, field_latex: str) -> str:
    """Render an operator applied to a field in LaTeX.

    Parameters
    ----------
    operator : str
        Operator name (e.g., "laplacian_x", "gradient_y", "identity").
    field_latex : str
        Already-rendered LaTeX for the field.

    Returns
    -------
    str
        LaTeX expression for the term.
    """
    # Static operators
    if operator in _OPERATOR_LATEX:
        op_tex = _OPERATOR_LATEX[operator]
        if not op_tex:
            return field_latex
        if op_tex == "dot":
            return rf"\dot{{{field_latex}}}"
        return rf"{op_tex} {field_latex}"

    # Dynamic operator patterns
    result = _operator_dynamic(operator, field_latex)
    if result is not None:
        return result

    # Unknown operator: render verbatim
    return rf"\mathrm{{{operator}}}({field_latex})"


# ---------------------------------------------------------------------------
# Equation rendering
# ---------------------------------------------------------------------------


def equation_to_latex(
    eq: ComponentEquation,
    spec: EquationSystem,
) -> str:
    """Convert a single component equation to LaTeX.

    Parameters
    ----------
    eq : ComponentEquation
        The equation to render.
    spec : EquationSystem
        The parent equation system (for coordinates and tensor metadata).

    Returns
    -------
    str
        LaTeX string (without environment wrapping).
    """
    field_meta = _get_field_meta(eq.field_name, spec)
    coords = _tensor_coordinates(spec)
    axis_remap = _operator_axis_remap(spec)

    # LHS
    field_tex = field_to_latex(
        eq.field_name, tensor_meta=field_meta, coordinates=coords
    )
    t_order = eq.time_derivative_order
    if t_order == 0:
        lhs = "0"
    elif t_order == 1:
        lhs = rf"\partial_t {field_tex}"
    else:
        lhs = rf"\partial_t^{{{t_order}}} {field_tex}"

    # RHS
    rhs_parts: list[str] = []
    for i, term in enumerate(eq.rhs_terms):
        term_field_meta = _get_field_meta(term.field, spec)
        tf_tex = field_to_latex(
            term.field, tensor_meta=term_field_meta, coordinates=coords
        )
        op_tex = operator_to_latex(term.operator, tf_tex)

        # Coefficient
        coeff_str = _render_term_coefficient(
            term.coefficient, term.coefficient_symbolic, is_first=(i == 0)
        )
        if coeff_str:
            rhs_parts.append(f"{coeff_str} {op_tex}")
        else:
            rhs_parts.append(op_tex)

    rhs = " ".join(rhs_parts) if rhs_parts else "0"
    result = rf"{lhs} &= {rhs}"

    # Apply axis remapping for plane-wave reduced specs: e.g., \partial_x → \partial_z
    for from_axis, to_axis in axis_remap.items():
        result = result.replace(rf"\partial_{from_axis}", rf"\partial_{to_axis}")

    return result


_COEFF_TOL = 1e-12


_FRAC_TOL = 1e-12  # tolerance for fraction approximation


def _format_numeric_coeff(value: float) -> str:
    """Format a numeric coefficient, using fractions when exact."""
    frac = Fraction(value).limit_denominator(10000)
    # Check if the fraction is a good approximation
    if abs(float(frac) - value) < _FRAC_TOL:
        if frac.denominator == 1:
            return str(frac.numerator)
        sign = "-" if frac.numerator < 0 else ""
        return rf"{sign}\frac{{{abs(frac.numerator)}}}{{{frac.denominator}}}"
    return f"{value:g}"


def _render_term_coefficient(
    numeric: float, symbolic: str | None, *, is_first: bool
) -> str:
    """Render a term's coefficient for display in an equation.

    Returns empty string if the coefficient is effectively +1 (for first term)
    or "+" sign-only otherwise.
    """
    if symbolic is not None:
        # If symbolic contains unresolvable Mathematica (e.g., Derivative[...]),
        # fall back to the numeric value rendered as a fraction.
        if re.search(r"Derivative\[|PD\w+\[", symbolic):
            symbolic = None  # fall through to numeric path below
        else:
            tex = coefficient_to_latex(symbolic)
            if not is_first and not tex.lstrip().startswith("-"):
                return f"+ {tex}"
            return tex

    # Numeric only — try to render as fraction if possible
    if abs(numeric - 1.0) < _COEFF_TOL:
        return "" if is_first else "+"
    if abs(numeric + 1.0) < _COEFF_TOL:
        return "-"
    formatted = _format_numeric_coeff(numeric)
    if not is_first and numeric > 0:
        formatted = f"+ {formatted}"
    return formatted


# ---------------------------------------------------------------------------
# Hamiltonian rendering
# ---------------------------------------------------------------------------


def hamiltonian_to_latex(
    terms: list[HamiltonianTerm],
    spec: EquationSystem,
) -> str:
    r"""Render the Hamiltonian density as a LaTeX equation.

    Returns
    -------
    str
        LaTeX for ``\\mathcal{H} = ...``.
    """
    coords = _tensor_coordinates(spec)
    axis_remap = _operator_axis_remap(spec)
    parts: list[str] = []
    for i, term in enumerate(terms):
        fa_meta = _get_field_meta(term.factor_a.field, spec)
        fb_meta = _get_field_meta(term.factor_b.field, spec)
        fa_tex = field_to_latex(
            term.factor_a.field, tensor_meta=fa_meta, coordinates=coords
        )
        fb_tex = field_to_latex(
            term.factor_b.field, tensor_meta=fb_meta, coordinates=coords
        )
        fa_op = operator_to_latex(term.factor_a.operator, fa_tex)
        fb_op = operator_to_latex(term.factor_b.operator, fb_tex)

        coeff = _render_term_coefficient(
            term.coefficient, term.coefficient_symbolic, is_first=(i == 0)
        )
        # Quadratic form: coeff * factor_a * factor_b
        term_tex = rf"{fa_op}^2" if fa_op == fb_op else rf"{fa_op} \, {fb_op}"

        if coeff:
            parts.append(f"{coeff} {term_tex}")
        else:
            parts.append(term_tex)

    rhs = " ".join(parts) if parts else "0"
    result = rf"\mathscr{{H}} &= {rhs}"

    # Apply axis remapping for plane-wave reduced specs
    for from_axis, to_axis in axis_remap.items():
        result = result.replace(rf"\partial_{from_axis}", rf"\partial_{to_axis}")

    return result


# ---------------------------------------------------------------------------
# Lagrangian rendering (xAct notation → LaTeX with \tensor{})
# ---------------------------------------------------------------------------

# Special tensor name mappings for LaTeX.
# NOTE: "eta" is NOT included here — it is handled conditionally
# by _metric_symbol (set per-render based on flat/curved/linearized).
_TENSOR_NAME_MAP: dict[str, str] = {
    "bg": r"\bar{g}",
    "epsiloneta": r"\epsilon",
    "TorsionCDT": "T",
    "RicciScalarCDT": r"\tilde{\mathcal{R}}",
    "RicciScalarCD": r"\mathcal{R}",
}

# Module-level metric symbol for the current render pass.
# Set by system_to_latex() before calling lagrangian_to_latex().
# \eta for flat unperturbed Minkowski, g for curved or linearized.
_metric_symbol: str = r"\eta"

# Covariant derivative pattern: CD[-a][expr] or CD[{N, -chart}][expr]
_RE_CD_ABSTRACT = re.compile(r"CD\[(-?\w+)\]\[([^\[\]]+(?:\[[^\[\]]*\])*)\]")
_RE_CD_BASIS = re.compile(r"CD\[\{(\d+),\s*-\w+\}\]\[([^\[\]]+(?:\[[^\[\]]*\])*)\]")

# Tensor object with indices: name[-a, b, -c] or name[a, b]
_RE_TENSOR_INDICES = re.compile(r"(\w+)\[((?:-?\w+(?:\s*,\s*-?\w+)*))\]")

# Scalar field with empty brackets: phi[]
_RE_SCALAR_FIELD = re.compile(r"(\w+)\[\]")


def _replace_cd_basis(m: re.Match[str]) -> str:
    """Replace basis covariant derivative CD[{N, -chart}][expr]."""
    idx = int(m.group(1))
    inner = _lagrangian_inner(m.group(2))
    coord_labels = ("t", "x", "y", "z", "w", "v")
    label = coord_labels[idx] if idx < len(coord_labels) else str(idx)
    return rf"\partial_{{{label}}} {inner}"


def _replace_cd_abstract(m: re.Match[str]) -> str:
    """Replace abstract covariant derivative CD[-a][expr]."""
    idx = m.group(1)
    inner = _lagrangian_inner(m.group(2))
    label = idx.lstrip("-")
    return rf"\nabla_{{{label}}} {inner}"


_SKIP_TENSOR_NAMES = frozenset({"CD", "Sqrt", "Rational", "Exp", "Log", "Sin", "Cos"})


def _indices_to_tensor_spec(indices_str: str) -> str:
    r"""Parse comma-separated xAct indices into \\tensor{} index specification."""
    raw_indices = [idx.strip() for idx in indices_str.split(",")]
    parts: list[str] = []
    for idx in raw_indices:
        if idx.startswith("-"):
            parts.append(f"_{idx[1:]}")  # lowered
        else:
            parts.append(f"^{idx}")  # raised
    return "".join(parts)


def _tensor_head_to_latex(name: str) -> str:
    """Map a tensor head name to its LaTeX representation."""
    if name == "eta":
        return _metric_symbol
    if name in _TENSOR_NAME_MAP:
        return _TENSOR_NAME_MAP[name]
    if name in _GREEK_MAP:
        return _GREEK_MAP[name]
    return name


def _replace_tensor_match(m: re.Match[str]) -> str:
    r"""Replace a tensor object with indices: name[idx1, idx2, ...] → \\tensor{}."""
    name = m.group(1)
    if name in _SKIP_TENSOR_NAMES:
        return m.group(0)
    index_spec = _indices_to_tensor_spec(m.group(2))
    name_tex = _tensor_head_to_latex(name)
    return rf"\tensor{{{name_tex}}}{{{index_spec}}}"


def _replace_scalar_field(m: re.Match[str]) -> str:
    r"""Replace scalar field with empty brackets: phi[] → \\phi or \\mathcal{H}."""
    name = m.group(1)
    if name in {"G", "V"}:
        return rf"{name}(\mathbf{{x}})"
    greek = _GREEK_MAP.get(name)
    if greek:
        return greek
    # Non-Greek field head: calligraphic
    return _calligraphic_head(name)


def _paren_frac(m: re.Match[str]) -> str:
    r"""Convert parenthesized fraction (A/B) → \\frac{A}{B}."""
    inner = m.group(1)
    slash = inner.find("/")
    if slash > 0:
        return rf"\frac{{{inner[:slash].strip()}}}{{{inner[slash + 1 :].strip()}}}"
    return m.group(0)


# Pre-compiled pattern for Greek in Lagrangian cleanup (negative lookbehind)
_RE_GREEK_NO_BACKSLASH = re.compile(
    r"(?<!\\)\b(" + "|".join(sorted(_GREEK_MAP, key=len, reverse=True)) + r")\b"
)


def _lagrangian_cleanup(s: str) -> str:
    """Apply final cleanup to Lagrangian LaTeX output."""
    # Convert Mathematica functions before any other processing
    s = _convert_math_functions(s)
    s = s.replace("*", r" \, ")
    # Parenthesized fractions: (A/B) → \frac{A}{B} (before Greek, so names stay intact)
    s = re.sub(r"\(([^()]+/[^()]+)\)", _paren_frac, s)
    # Simple fractions: A/B where A,B are word tokens (before Greek substitution)
    s = re.sub(
        r"(?<![\\{])(\w+)/((?:\w+(?:\^[{\d]+}?)?))",
        lambda m: rf"\frac{{{m.group(1)}}}{{{m.group(2)}}}",
        s,
    )
    # Greek prefix extraction: omegaP2 → omega P2, so Greek + subscript work
    # Inserts a space between a Greek name and a trailing uppercase letter.
    for greek in sorted(_GREEK_MAP, key=len, reverse=True):
        s = re.sub(rf"\b{greek}([A-Z])", rf"{greek} \1", s)
    # Subscript splitting BEFORE Greek so that alpha1 → alpha_{1} → \alpha_{1}
    # (Greek regex uses \b which fails on alpha1 since 1 is a word char)
    s = re.sub(
        r"(?<!\\)\b([A-Za-z]+?)(\d+)\b",
        lambda m: rf"{m.group(1)}_{{{m.group(2)}}}",
        s,
    )
    # Greek substitution for remaining parameter names (skip already-escaped).
    # After subscript splitting, alpha_{1} has _ after alpha, which is a word char,
    # so the standard \b boundary fails. Use lookahead for \b OR _ OR {.
    s = re.sub(
        r"(?<!\\)\b("
        + "|".join(sorted(_GREEK_MAP, key=len, reverse=True))
        + r")(?=\b|[_{])",
        lambda m: _GREEK_MAP[m.group(1)],
        s,
    )
    # Powers: ^(expr) → ^{expr}
    s = _RE_POWER_PAREN.sub(lambda m: f"^{{{m.group(1)}}}", s)
    # Clean up double spaces
    return re.sub(r"\s+", " ", s).strip()


def lagrangian_to_latex(expr: str) -> str:
    r"""Convert a Lagrangian expression from xAct notation to LaTeX.

    This is a best-effort conversion. The xAct abstract index notation is
    rich and idiosyncratic; this handles the patterns found in the 33
    example JSONs in this project.

    Parameters
    ----------
    expr : str
        The ``lagrangian_expr`` string from JSON metadata.

    Returns
    -------
    str
        LaTeX representation using ``\\tensor{}`` for index placement.
    """
    if not expr:
        return ""

    s = expr.strip()
    s = _strip_outer_parens(s)

    # Pass 1: Bracket functions
    s = _RE_RATIONAL.sub(
        lambda m: rf"\frac{{{m.group(1).strip()}}}{{{m.group(2).strip()}}}", s
    )
    s = _RE_SQRT.sub(lambda m: rf"\sqrt{{{m.group(1)}}}", s)
    s = _RE_E_POWER.sub(r"e^", s)
    s = _RE_PI.sub(r"\\pi", s)

    # Pass 2: Covariant derivatives
    s = _RE_CD_BASIS.sub(_replace_cd_basis, s)
    prev = ""
    while prev != s:
        prev = s
        s = _RE_CD_ABSTRACT.sub(_replace_cd_abstract, s)

    # Pass 3: Named special objects (longer names first to avoid prefix match)
    s = re.sub(r"\bRicciScalarCDT\[\]", r"\\tilde{\\mathcal{R}}", s)
    s = re.sub(r"\bRicciScalarCD\[\]", r"\\mathcal{R}", s)

    # Pass 4: Tensor objects with indices
    s = _RE_TENSOR_INDICES.sub(_replace_tensor_match, s)

    # Pass 5: Scalar fields (empty brackets)
    s = _RE_SCALAR_FIELD.sub(_replace_scalar_field, s)

    # Pass 6: Cleanup
    return _lagrangian_cleanup(s)


def _lagrangian_inner(fragment: str) -> str:
    """Recursively process a Lagrangian sub-expression.

    Used for the inner part of CD[-a][...] to handle nested tensors.
    """
    result = _RE_SCALAR_FIELD.sub(_replace_scalar_field, fragment)
    return _RE_TENSOR_INDICES.sub(_replace_tensor_match, result)


# ---------------------------------------------------------------------------
# System-level rendering
# ---------------------------------------------------------------------------

_DOCUMENT_PREAMBLE = r"""\documentclass{article}
\usepackage{amsmath}
\usepackage{tensor}
\usepackage{mathrsfs}

\begin{document}
"""

_DOCUMENT_POSTAMBLE = r"""
\end{document}
"""


def system_to_latex(
    spec: EquationSystem,
    *,
    output_format: str = "align",
    include_hamiltonian: bool = True,
    include_lagrangian: bool = True,
) -> str:
    """Convert a full equation system to LaTeX.

    Parameters
    ----------
    spec : EquationSystem
        The equation system to render.
    output_format : {"align", "document", "raw"}
        Output format.
    include_hamiltonian : bool
        Whether to include the Hamiltonian density.
    include_lagrangian : bool
        Whether to include the Lagrangian expression.

    Returns
    -------
    str
        LaTeX output.
    """
    global _metric_symbol  # noqa: PLW0603

    # Determine metric symbol: η for flat unperturbed Minkowski, g otherwise
    is_linearized = spec.metadata.get("linearized", False)
    metric_type = spec.metadata.get("metric_type", "minkowski")
    _metric_symbol = (
        r"\eta" if (not is_linearized and metric_type == "minkowski") else "g"
    )

    sections: list[str] = []

    # Lagrangian
    lagrangian_expr = spec.metadata.get("lagrangian_expr", "")
    if include_lagrangian and lagrangian_expr:
        lag_tex = lagrangian_to_latex(lagrangian_expr)
        sections.append(rf"\mathcal{{L}} &= {lag_tex}")

    # Equations of motion
    sections.extend(equation_to_latex(eq, spec) for eq in spec.equations)

    # Hamiltonian
    if include_hamiltonian and spec.canonical and spec.canonical.hamiltonian_terms:
        sections.append(
            hamiltonian_to_latex(list(spec.canonical.hamiltonian_terms), spec)
        )

    if output_format == "raw":
        return "\n".join(sections)

    # Build align environment
    body = " \\\\\n  ".join(sections)
    align_block = f"\\begin{{align}}\n  {body}\n\\end{{align}}"

    if output_format == "document":
        # Extract title from metadata
        source = spec.metadata.get("source", "TIDAL")
        gauge = spec.metadata.get("gauge", "")
        title_parts = [f"Equations from {source}"]
        if gauge and gauge != "none":
            title_parts.append(f"(gauge: {gauge})")
        title = " ".join(title_parts)

        return (
            _DOCUMENT_PREAMBLE
            + f"\\title{{{title}}}\n\\maketitle\n\n"
            + "% Requires: \\usepackage{tensor} for index notation\n"
            + align_block
            + _DOCUMENT_POSTAMBLE
        )

    # align format: include a package reminder comment
    return "% Requires: \\usepackage{amsmath, tensor, mathrsfs}\n" + align_block


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_field_meta(
    field_name: str, spec: EquationSystem
) -> dict[str, list[int] | int | str] | None:
    """Look up tensor metadata for a field from the EquationSystem.

    Returns None if metadata is not available (backward compat with old JSONs).
    """
    tensor_metadata = spec.metadata.get("tensor_metadata")
    if tensor_metadata is None:
        return None
    return tensor_metadata.get(field_name)  # type: ignore[no-any-return]


#: Default axis letters for coordinate labelling (matching AXIS_LETTERS in json_loader).
_DEFAULT_COORD_LABELS: tuple[str, ...] = ("t", "x", "y", "z", "w", "v")


def _tensor_coordinates(spec: EquationSystem) -> tuple[str, ...]:
    """Return coordinate labels suitable for tensor index rendering.

    For plane-wave reduced specs, the tensor indices refer to the *original*
    higher-dimensional coordinate system, not the reduced 1+1D coordinates.
    This function infers the original coordinates from the ``reduction``
    metadata when available.
    """
    reduction = spec.metadata.get("reduction")
    if reduction:
        orig_dim = reduction.get("original_dimension")
        if orig_dim and isinstance(orig_dim, int):
            return _DEFAULT_COORD_LABELS[:orig_dim]
    return spec.effective_coordinates


def _operator_axis_remap(spec: EquationSystem) -> dict[str, str]:
    r"""Build a mapping from reduced axis labels to original axis labels.

    In plane-wave reduced specs, the single spatial coordinate ``x`` in the
    1+1D PDE actually represents the original propagation axis (e.g., ``z``).
    Operators like ``\\partial_x`` should render as ``\\partial_z`` to avoid
    confusion with the tensor index label ``x`` (which refers to the original
    transverse x-coordinate).

    Returns an empty dict when no remapping is needed.
    """
    reduction = spec.metadata.get("reduction")
    if not reduction:
        return {}
    prop_axis = reduction.get("propagation_axis")
    if not prop_axis or prop_axis == "x":
        # No remapping needed: propagation is already along x
        return {}
    # The reduced PDE uses "x" for the single spatial axis, but it
    # physically represents the propagation axis (e.g., "z")
    return {"x": prop_axis}
