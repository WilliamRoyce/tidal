"""Safe symbolic evaluation of kinetic-coefficient strings at ε=0.

Used by :meth:`tidal.symbolic.json_loader.EquationSystem.base_spec` (v6
perturbative reduction, Gap B) to detect when the ``lhs`` of a
higher-derivative equation vanishes at zero small-parameter value and
must therefore be demoted to an algebraic constraint.

The :meth:`evaluate_at_zero` function restricts itself to a small
subset of :mod:`ast` nodes (literals, arithmetic, pow, unary) so only
polynomials and rational expressions in named parameters can be parsed.
Arbitrary code execution is impossible by construction.

Wolfram-emitted coefficient strings use ``^`` for exponentiation;
:meth:`evaluate_at_zero` substitutes ``^`` → ``**`` before parsing.

Reference: mirrors the restricted-AST pattern in
``tidal/inference/_constraints.py``.
"""

from __future__ import annotations

import ast
import operator
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable, Sequence

# Safe binary / unary operators.
_BIN_OPS: dict[type, Callable[[float, float], float]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
}

_UNARY_OPS: dict[type, Callable[[float], float]] = {
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
}


class KineticEvalError(ValueError):
    """Raised when a kinetic coefficient string cannot be parsed safely."""


def evaluate_at_zero(expr: str, zero_names: set[str] | frozenset[str]) -> float | None:
    """Return the literal value of ``expr`` after substituting every name in
    ``zero_names`` with ``0.0``.

    Returns
    -------
    float
        The evaluated literal (possibly ``0.0``) if the substitution
        leaves no un-evaluated symbol behind.
    None
        If, after substitution, at least one ``Name`` node survives
        (i.e. the expression is still symbolic in some non-small
        parameter). The caller distinguishes "literal zero — demote"
        (result == 0.0) from "still symbolic — keep" (result is None).

    Parameters
    ----------
    expr : str
        The coefficient expression as emitted by ExportJSON.wl (e.g.,
        ``"2*b5"``, ``"(-25*b5)/2"``, ``"-kappa^(-2)"``).  Caret
        exponents are converted to Python ``**`` form before parsing.
    zero_names : set[str]
        The names to substitute with zero. Typically the small
        parameters declared in the theory's ``[perturbation]`` section.

    Raises
    ------
    KineticEvalError
        If ``expr`` contains a node kind outside the safe subset
        (literals, names, unary +/-, binary +/-/*/////**).
    """
    if not isinstance(expr, str):
        msg = f"evaluate_at_zero: expected a str, got {type(expr).__name__}"
        raise KineticEvalError(msg)

    # Wolfram uses ``^`` for exponentiation; Python expects ``**``.
    normalised = expr.replace("^", "**")

    try:
        tree = ast.parse(normalised, mode="eval")
    except SyntaxError as e:
        msg = (
            f"evaluate_at_zero: cannot parse {expr!r} (after ^→**: {normalised!r}): {e}"
        )
        raise KineticEvalError(msg) from e

    try:
        return _eval_node(tree.body, zero_names)
    except _StillSymbolic:
        return None


class _StillSymbolic(Exception):
    """Signal raised when evaluation encounters a surviving Name."""


def _eval_node(node: ast.AST, zero_names: set[str] | frozenset[str]) -> float:
    """Recursively evaluate a safe AST node; raise _StillSymbolic if a
    Name outside ``zero_names`` is encountered.
    """
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return float(node.value)
        msg = (
            f"evaluate_at_zero: unsupported literal type {type(node.value).__name__!r}"
        )
        raise KineticEvalError(msg)

    if isinstance(node, ast.Name):
        if node.id in zero_names:
            return 0.0
        raise _StillSymbolic

    if isinstance(node, ast.UnaryOp):
        op_fn = _UNARY_OPS.get(type(node.op))
        if op_fn is None:
            msg = f"evaluate_at_zero: unsupported unary op {type(node.op).__name__}"
            raise KineticEvalError(msg)
        return op_fn(_eval_node(node.operand, zero_names))

    if isinstance(node, ast.BinOp):
        op_fn_bin = _BIN_OPS.get(type(node.op))
        if op_fn_bin is None:
            msg = f"evaluate_at_zero: unsupported binary op {type(node.op).__name__}"
            raise KineticEvalError(msg)
        left = _eval_node(node.left, zero_names)
        right = _eval_node(node.right, zero_names)
        # 0**0 is mathematically ambiguous; treat as 1 per IEEE, which is
        # also Python's convention. 0**negative is ZeroDivisionError from
        # operator.pow — propagate with a clearer message.
        if isinstance(node.op, ast.Pow) and left == 0.0 and right < 0:
            msg = (
                f"evaluate_at_zero: divide-by-zero in expression "
                f"(0 ** {right} undefined)"
            )
            raise KineticEvalError(msg)
        if isinstance(node.op, ast.Div) and right == 0.0:
            msg = "evaluate_at_zero: divide-by-zero in expression"
            raise KineticEvalError(msg)
        return op_fn_bin(left, right)

    msg = f"evaluate_at_zero: unsupported AST node {type(node).__name__}"
    raise KineticEvalError(msg)


def lhs_collapses_to_zero(
    kinetic_symbolic: str | None, small_parameters: Sequence[str]
) -> bool:
    """Return True iff the LHS kinetic coefficient evaluates to literal 0
    when every small parameter is set to zero.

    ``None`` kinetic coefficients are implicitly non-perturbative
    (e.g. ``kinetic == 1`` by convention); they never trigger demotion.
    """
    if kinetic_symbolic is None or not small_parameters:
        return False
    zero_names = frozenset(small_parameters)
    value = evaluate_at_zero(kinetic_symbolic, zero_names)
    return value == 0.0
