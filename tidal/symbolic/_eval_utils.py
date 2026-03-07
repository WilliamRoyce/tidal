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
    result = expr

    # E^(...) → exp(...)
    result = _RE_E_POWER.sub("exp", result)

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

    return result


# ------------------------------------------------------------------
# Step 2: Evaluation namespace
# ------------------------------------------------------------------


def build_eval_namespace(parameters: dict[str, float]) -> dict[str, object]:
    """Build the evaluation namespace for ``eval()`` of converted expressions.

    Contains numpy/scipy math functions and user-provided parameters.

    Parameters
    ----------
    parameters : dict[str, float]
        User-provided parameter values (e.g. ``{"g0": 1.0, "R": 8.0}``).
    """
    ns: dict[str, object] = dict(parameters)
    ns["exp"] = np.exp
    ns["sin"] = np.sin
    ns["cos"] = np.cos
    ns["tan"] = np.tan
    ns["cot"] = lambda x: np.cos(x) / np.sin(x)  # type: ignore[reportUnknownLambdaType]
    ns["sec"] = lambda x: 1.0 / np.cos(x)  # type: ignore[reportUnknownLambdaType]
    ns["csc"] = lambda x: 1.0 / np.sin(x)  # type: ignore[reportUnknownLambdaType]
    ns["arcsin"] = np.arcsin
    ns["arccos"] = np.arccos
    ns["arctan"] = np.arctan
    ns["arctan2"] = np.arctan2
    ns["sinh"] = np.sinh
    ns["cosh"] = np.cosh
    ns["tanh"] = np.tanh
    ns["arcsinh"] = np.arcsinh
    ns["arccosh"] = np.arccosh
    ns["arctanh"] = np.arctanh
    ns["log"] = np.log
    ns["sqrt"] = np.sqrt
    ns["abs"] = np.abs
    ns["sign"] = np.sign
    ns["maximum"] = np.maximum
    ns["minimum"] = np.minimum
    ns["heaviside"] = lambda x: np.heaviside(x, 0.5)  # type: ignore[reportUnknownLambdaType]
    ns["piecewise"] = np.where
    ns["erf"] = special.erf
    ns["jv"] = special.jv
    ns["yv"] = special.yv
    ns["np"] = np
    ns["True"] = True
    ns["False"] = False
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
