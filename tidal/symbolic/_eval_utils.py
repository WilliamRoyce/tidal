"""Standalone evaluation of Mathematica symbolic expressions on numpy grids.

This module provides functions to convert Mathematica InputForm expressions
to Python and evaluate them on coordinate arrays.  It is used by both the
PDE builder (``PDEFromSpec``) and the measurement module (``_energy``) to
resolve position-dependent coefficients.

The conversion pipeline:
    Mathematica InputForm → Python string → eval() with numpy namespace

Primary public entry point:

- ``evaluate_coefficient(expr, parameters, coordinates, ...)`` — full eval

Building blocks (also used directly by ``coefficients.CoefficientEvaluator``):

- ``mathematica_to_python(expr, coordinates)`` — string-level conversion
- ``build_eval_namespace(parameters)`` — static evaluation namespace
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import special  # type: ignore[import-untyped]

if TYPE_CHECKING:
    from numpy.typing import NDArray


# ------------------------------------------------------------------
# Step 1: Mathematica InputForm → Python expression string
# ------------------------------------------------------------------

_FUNCTION_MAP: list[tuple[str, str]] = [
    # Basic trig
    ("Sin", "sin"),
    ("Cos", "cos"),
    ("Tan", "tan"),
    # Reciprocal trig
    ("Cot", "cot"),
    ("Sec", "sec"),
    ("Csc", "csc"),
    # Inverse trig (1-arg)
    ("ArcSin", "arcsin"),
    ("ArcCos", "arccos"),
    ("ArcTan", "arctan"),
    # Hyperbolic
    ("Sinh", "sinh"),
    ("Cosh", "cosh"),
    ("Tanh", "tanh"),
    # Inverse hyperbolic
    ("ArcSinh", "arcsinh"),
    ("ArcCosh", "arccosh"),
    ("ArcTanh", "arctanh"),
    # Other
    ("Exp", "exp"),
    ("Log", "log"),
    ("Sqrt", "sqrt"),
    ("Abs", "abs"),
    ("Sign", "sign"),
    ("Max", "maximum"),
    ("Min", "minimum"),
    # Step functions
    ("UnitStep", "heaviside"),
    ("HeavisideTheta", "heaviside"),
    # Special functions
    ("Erf", "erf"),
    ("BesselJ", "jv"),
    ("BesselY", "yv"),
]

# Pre-compiled regex patterns for mathematica_to_python (avoid re-compilation per call)
_RE_E_POWER = re.compile(r"\bE\^")
_RE_RATIONAL = re.compile(r"Rational\[([^,\]]+),\s*([^,\]]+)\]")
_RE_PI = re.compile(r"\bPi\b")
_COMPILED_FUNCTION_PATTERNS: list[tuple[re.Pattern[str], str]] = [
    (re.compile(rf"\b{mma}\b"), py) for mma, py in _FUNCTION_MAP
]


def _find_matching_paren(s: str, start: int) -> int:
    """Return index *past* the closing paren matching the open-paren at *start*."""
    depth = 1
    j = start + 1
    n = len(s)
    while j < n and depth > 0:
        if s[j] == "(":
            depth += 1
        elif s[j] == ")":
            depth -= 1
        j += 1
    return j


def _replace_compound_denom(denom: str) -> str | None:
    """If *denom* contains ``exp(arg)``, return the numerator replacement.

    Returns ``None`` if no ``exp(`` is found inside *denom*.
    """
    exp_pos = denom.find("exp(")
    if exp_pos < 0:
        return None
    ej = _find_matching_paren(denom, exp_pos + 3)
    exp_arg = denom[exp_pos + 4 : ej - 1]
    pre = denom[:exp_pos].rstrip("*")
    post = denom[ej:].lstrip("*")
    remaining = f"{pre}*{post}" if (pre and post) else (pre or post)
    if remaining:
        return f"*exp(-({exp_arg}))/({remaining})"
    return f"*exp(-({exp_arg}))"


def _invert_exp_denominator(expr: str) -> str:
    """Convert exp-in-denominator patterns to ``*exp(-(arg))`` to avoid overflow.

    Wolfram's InputForm may serialize ``Exp[-x^2/R^2]`` in several forms that
    all produce positive-exponent overflow when evaluated naively:

    1. ``k/E^(x^2/R^2)`` after ``E^`` to ``exp``: ``k/exp(x^2/R^2)``
       Handled as: ``/exp(arg)`` to ``*exp(-(arg))``

    2. ``k/(A*E^(x^2/R^2)*B)`` after ``E^`` to ``exp``: ``k/(A*exp(x^2/R^2)*B)``
       Handled as: ``/(A*exp(arg)*B)`` to ``*exp(-(arg))/(A*B)``

    Both patterns prevent evaluating ``exp(+large)``, eliminating numpy
    ``RuntimeWarning: overflow encountered in exp``.
    """
    result: list[str] = []
    i = 0
    n = len(expr)
    while i < n:
        if expr[i : i + 5] == "/exp(":
            # Pattern 1: /exp(arg) → *exp(-(arg))
            j = _find_matching_paren(expr, i + 4)
            inner = expr[i + 5 : j - 1]
            result.append(f"*exp(-({inner}))")
            i = j
        elif expr[i : i + 2] == "/(":
            # Pattern 2: /(denom) where denom may contain exp(arg) anywhere
            k = _find_matching_paren(expr, i + 1)
            denom = expr[i + 2 : k - 1]
            replacement = _replace_compound_denom(denom)
            if replacement is not None:
                result.append(replacement)
                i = k
            else:
                result.append(expr[i])
                i += 1
        else:
            result.append(expr[i])
            i += 1
    return "".join(result)


_COMPARISON_OPS: dict[str, str] = {
    "LessEqual": "<=",
    "Less": "<",
    "GreaterEqual": ">=",
    "Greater": ">",
    "Equal": "==",
}


# ------------------------------------------------------------------
# Bracket-aware parsing helpers
# ------------------------------------------------------------------


def _split_bracket_aware(s: str) -> list[str]:
    """Split string on commas, respecting ``[...]`` bracket nesting.

    Parameters
    ----------
    s : str
        Content string (typically inside a Mathematica function call).

    Returns
    -------
    list[str]
        Comma-separated parts, preserving nested bracket structure.
    """
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    for ch in s:
        if ch == "[":
            depth += 1
            current.append(ch)
        elif ch == "]":
            depth -= 1
            current.append(ch)
        elif ch == "," and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(ch)
    if current:
        parts.append("".join(current))
    return parts


def _find_top_level_comma(s: str) -> int:
    """Find index of first comma not inside ``()[]{}``."""
    depth = 0
    for i, ch in enumerate(s):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            return i
    return -1


def _find_matching_brace(s: str) -> int:
    """Find index of the closing ``}`` that matches the opening ``{`` at s[0]."""
    depth = 0
    for i, ch in enumerate(s):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return i
    return -1


def _extract_brace_pairs(cases_str: str) -> list[tuple[str, str]]:
    """Extract ``{value, condition}`` pairs from the inner list string."""
    cases: list[tuple[str, str]] = []
    i = 0
    while i < len(cases_str):
        if cases_str[i] == "{":
            end = _find_matching_brace(cases_str[i:])
            if end < 0:
                break
            case_content = cases_str[i + 1 : i + end]
            comma_idx = _find_top_level_comma(case_content)
            if comma_idx >= 0:
                val = case_content[:comma_idx].strip()
                cond = case_content[comma_idx + 1 :].strip()
                cases.append((val, cond))
            i += end + 1
        else:
            i += 1
    return cases


def _parse_piecewise_content(
    content: str,
) -> tuple[list[tuple[str, str]], str]:
    """Parse the inner content of a ``Piecewise[...]`` expression.

    Parameters
    ----------
    content : str
        Everything between ``Piecewise[`` and the closing ``]``.
        Expected form: ``{{val1, cond1}, {val2, cond2}, ...}, default``

    Returns
    -------
    tuple[list[tuple[str, str]], str]
        List of ``(value, condition)`` pairs and the default value string.
    """
    content = content.strip()
    if not content.startswith("{{"):
        return [], content

    end_idx = _find_matching_brace(content)
    if end_idx < 0:
        return [], content  # Malformed

    cases_str = content[1:end_idx]  # Strip outer { ... }
    default_str = content[end_idx + 1 :].strip()
    if default_str.startswith(","):
        default_str = default_str[1:].strip()
    if not default_str:
        default_str = "0"

    return _extract_brace_pairs(cases_str), default_str


# ------------------------------------------------------------------
# Conversion functions
# ------------------------------------------------------------------


def _convert_power(expr: str) -> str:
    """Convert ``Power[base, exp]`` to ``(base)**(exp)``.

    Uses manual bracket counting to handle arbitrary nesting depth.
    """
    pattern = r"Power\["

    def _replace(expr: str) -> str:
        match = re.search(pattern, expr)
        if match is None:
            return expr

        start = match.end()
        depth = 1
        comma_pos = -1
        i = start
        while i < len(expr) and depth > 0:
            if expr[i] == "[":
                depth += 1
            elif expr[i] == "]":
                depth -= 1
            elif expr[i] == "," and depth == 1 and comma_pos < 0:
                comma_pos = i
            i += 1

        if comma_pos < 0 or depth != 0:
            return expr

        base = expr[start:comma_pos].strip()
        exponent = expr[comma_pos + 1 : i - 1].strip()
        return expr[: match.start()] + f"({base})**({exponent})" + expr[i:]

    prev = None
    result = expr
    while prev != result:
        prev = result
        result = _replace(result)
    return result


def _convert_arctan2(expr: str) -> str:
    """Convert ``ArcTan[x, y]`` to ``arctan2(y, x)`` (note arg swap).

    Uses manual bracket counting to handle arbitrary nesting depth.
    """
    pattern = r"ArcTan\["

    def _replace(expr: str) -> str:
        match = re.search(pattern, expr)
        if match is None:
            return expr

        start = match.end()
        depth = 1
        comma_pos = -1
        i = start
        while i < len(expr) and depth > 0:
            if expr[i] == "[":
                depth += 1
            elif expr[i] == "]":
                depth -= 1
            elif expr[i] == "," and depth == 1 and comma_pos < 0:
                comma_pos = i
            i += 1

        if comma_pos < 0 or depth != 0:
            return expr

        x_arg = expr[start:comma_pos].strip()
        y_arg = expr[comma_pos + 1 : i - 1].strip()
        return expr[: match.start()] + f"arctan2({y_arg}, {x_arg})" + expr[i:]

    prev = None
    result = expr
    while prev != result:
        prev = result
        result = _replace(result)
    return result


def _convert_inequality(expr: str) -> str:
    """Convert ``Inequality[a, op, b, op, c]`` to comparison chains.

    Uses bracket-aware splitting to handle nested bracket expressions
    inside the arguments (e.g. ``x[]``).
    """
    pattern = r"Inequality\[((?:[^[\]]|\[[^\]]*\])*)\]"

    def replacer(match: re.Match[str]) -> str:
        inner = match.group(1)
        args = _split_bracket_aware(inner)
        min_args = 3  # Inequality[a, op, b] minimum
        if len(args) < min_args or len(args) % 2 == 0:
            return match.group(0)  # Malformed, leave unchanged
        parts: list[str] = []
        for i in range(0, len(args) - 2, 2):
            left = args[i].strip()
            op_name = args[i + 1].strip()
            op_sym = _COMPARISON_OPS.get(op_name, op_name)
            right = args[i + 2].strip()
            parts.append(f"(({left}) {op_sym} ({right}))")
        return " & ".join(parts)

    prev = None
    result = expr
    while prev != result:
        prev = result
        result = re.sub(pattern, replacer, result)
    return result


def _convert_piecewise(expr: str) -> str:
    """Convert ``Piecewise[{{val, cond}, ...}, default]`` to ``piecewise()`` calls.

    Uses structured parsing with ``_parse_piecewise_content`` to correctly
    handle nested brackets and multiple case branches.
    """
    pattern = r"Piecewise\[((?:[^[\]]|\[[^\]]*\])*)\]"

    def replacer(match: re.Match[str]) -> str:
        content = match.group(1)
        cases, default = _parse_piecewise_content(content)
        if not cases:
            return default
        result = default
        for val, cond in reversed(cases):
            result = f"piecewise({cond}, {val}, {result})"
        return result

    prev = None
    result = expr
    while prev != result:
        prev = result
        result = re.sub(pattern, replacer, result)
    return result


# ------------------------------------------------------------------
# Main conversion entry point
# ------------------------------------------------------------------


# Per-process cache of Mathematica→Python translations.  In nested
# sampling the same ~14 coefficient strings are re-translated ~13k
# times per likelihood evaluation (see #291 profiling), purely because
# the translator is stateless and gets called fresh each time.  The
# cache is keyed on (expression string, coordinate tuple) — deterministic
# in both inputs — so it's safe to share across calls within a process.
# Per-process (not shared) so MPI ranks / multiprocessing workers each
# keep their own copy, matching the #143 BLAS-pinning convention.
_PARSED_MATH_CACHE: dict[tuple[str, tuple[str, ...]], str] = {}


def mathematica_to_python(
    expr: str,
    coordinates: tuple[str, ...] = ("t", "x", "y"),
) -> str:
    """Convert a Mathematica InputForm expression to evaluable Python.

    Parameters
    ----------
    expr : str
        Mathematica InputForm expression string.
    coordinates : tuple[str, ...]
        Coordinate names (e.g. ``("t", "x", "y")``).  Used to strip
        xCoba ``()`` from coordinate symbols.

    Returns
    -------
    str
        Python-evaluable expression string.
    """
    # Coerce to hashable form and look up.
    cache_coords: tuple[str, ...] = tuple(coordinates)
    cached = _PARSED_MATH_CACHE.get((expr, cache_coords))
    if cached is not None:
        return cached

    result = expr

    # E^(...) → exp(...)
    result = _RE_E_POWER.sub("exp", result)

    # /exp(arg) → *exp(-(arg)) to prevent overflow from positive exponents in denominators
    # (Wolfram may serialize Exp[-x^2/R^2] as 1/E^(x^2/R^2))
    result = _invert_exp_denominator(result)

    # Power[x, y] → (x)**(y)
    result = _convert_power(result)

    # ArcTan[x, y] → arctan2(y, x) (2-arg special case)
    result = _convert_arctan2(result)

    # Rational[p, q] → (p)/(q)
    result = _RE_RATIONAL.sub(r"(\1)/(\2)", result)

    # Inequality → comparison chains
    result = _convert_inequality(result)

    # Piecewise → piecewise()
    result = _convert_piecewise(result)

    # Function name conversions
    for pat, py_func in _COMPILED_FUNCTION_PATTERNS:
        result = pat.sub(py_func, result)

    # Pi → np.pi
    result = _RE_PI.sub("np.pi", result)

    # Mathematica brackets → Python parens
    result = result.replace("[", "(").replace("]", ")")

    # Mathematica ^ → Python **
    result = result.replace("^", "**")

    # xCoba coordinate symbols: t() → t, x() → x, etc.
    for coord in coordinates:
        result = result.replace(f"{coord}()", coord)

    _PARSED_MATH_CACHE[expr, cache_coords] = result
    return result


# ------------------------------------------------------------------
# Step 2: Evaluation namespace
# ------------------------------------------------------------------


# Static math-function namespace, built once.  build_eval_namespace
# used to allocate this 30-key dict on every call — in nested sampling
# that's ~400 allocations per likelihood evaluation (see #291).  The
# static entries don't depend on parameters, so we build them once and
# shallow-copy per call.
_STATIC_EVAL_NAMESPACE: dict[str, object] = {
    "exp": np.exp,
    "sin": np.sin,
    "cos": np.cos,
    "tan": np.tan,
    "cot": lambda x: np.cos(x) / np.sin(x),
    "sec": lambda x: 1.0 / np.cos(x),
    "csc": lambda x: 1.0 / np.sin(x),
    "arcsin": np.arcsin,
    "arccos": np.arccos,
    "arctan": np.arctan,
    "arctan2": np.arctan2,
    "sinh": np.sinh,
    "cosh": np.cosh,
    "tanh": np.tanh,
    "arcsinh": np.arcsinh,
    "arccosh": np.arccosh,
    "arctanh": np.arctanh,
    "log": np.log,
    "sqrt": np.sqrt,
    "abs": np.abs,
    "sign": np.sign,
    "maximum": np.maximum,
    "minimum": np.minimum,
    "heaviside": lambda x: np.heaviside(x, 0.5),
    "piecewise": np.where,
    "erf": special.erf,
    "jv": special.jv,
    "yv": special.yv,
    "np": np,
    "True": True,
    "False": False,
}


def build_eval_namespace(parameters: dict[str, float]) -> dict[str, object]:
    """Build the evaluation namespace for ``eval()`` of converted expressions.

    Contains numpy/scipy math functions and user-provided parameters.
    Math functions come from a module-level static dict; only the
    parameter values are merged per call.

    Parameters
    ----------
    parameters : dict[str, float]
        User-provided parameter values (e.g. ``{"g0": 1.0, "R": 8.0}``).
    """
    ns = _STATIC_EVAL_NAMESPACE.copy()
    ns.update(parameters)
    return ns


# ------------------------------------------------------------------
# Step 3: Full evaluation pipeline
# ------------------------------------------------------------------


def evaluate_coefficient(
    symbolic_expr: str,
    parameters: dict[str, float],
    coordinates: tuple[str, ...],
    coord_arrays: dict[str, NDArray[np.float64]] | None = None,
    t: float = 0.0,
) -> float | NDArray[np.float64]:
    """Evaluate a Mathematica symbolic expression to a numeric value or grid array.

    Parameters
    ----------
    symbolic_expr : str
        Mathematica InputForm expression (e.g. ``"g0 * Exp[-(x[]^2) / R^2]"``).
    parameters : dict[str, float]
        User-provided parameter values.
    coordinates : tuple[str, ...]
        Coordinate names including time (e.g. ``("t", "x", "y")``).
    coord_arrays : dict[str, NDArray] | None
        Spatial coordinate arrays (e.g. ``{"x": x_grid, "y": y_grid}``).
        If provided, the result is an ndarray; otherwise scalar.
    t : float
        Current simulation time (for time-dependent expressions).

    Returns
    -------
    float | NDArray[np.float64]
        Scalar for constant/time-only coefficients, grid array for
        position-dependent ones.

    Raises
    ------
    ValueError
        If the expression cannot be evaluated or produces NaN/Inf.
    TypeError
        If the result is complex.
    """
    py_expr = mathematica_to_python(symbolic_expr, coordinates)
    namespace = build_eval_namespace(parameters)
    namespace["t"] = t

    if coord_arrays is not None:
        namespace.update(coord_arrays)

    try:
        result = eval(py_expr, {"__builtins__": {}}, namespace)  # noqa: S307
    except Exception as e:
        msg = (
            f"Cannot evaluate symbolic coefficient '{symbolic_expr}' "
            f"(Python form: '{py_expr}'): {e}"
        )
        raise ValueError(msg) from e

    # Validate result
    if isinstance(result, complex):
        msg = (
            f"Coefficient '{symbolic_expr}' evaluated to complex number {result}. "
            f"Only real-valued coefficients are supported."
        )
        raise TypeError(msg)

    if isinstance(result, np.ndarray):
        arr = np.asarray(result, dtype=np.float64)
        if np.isnan(arr).any():
            msg = f"Coefficient '{symbolic_expr}' produced NaN values."
            raise ValueError(msg)
        if np.isinf(arr).any():
            msg = f"Coefficient '{symbolic_expr}' produced Inf values."
            raise ValueError(msg)
        return arr

    scalar = float(result)
    if np.isnan(scalar):
        msg = f"Coefficient '{symbolic_expr}' evaluated to NaN."
        raise ValueError(msg)
    if np.isinf(scalar):
        msg = f"Coefficient '{symbolic_expr}' evaluated to Inf."
        raise ValueError(msg)
    return scalar
