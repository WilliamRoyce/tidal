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

    Thin wrapper over :func:`evaluate_with_substitutions`; kept for its
    existing callers in ``json_loader.base_spec``.
    """
    return evaluate_with_substitutions(expr, dict.fromkeys(zero_names, 0.0))


def evaluate_at_one(expr: str, one_names: set[str] | frozenset[str]) -> float | None:
    """Return the literal value of ``expr`` after substituting every name in
    ``one_names`` with ``1.0``.

    v6 R8.2 / #290: used by :class:`PerturbativeSolver` to extract the
    numeric kinetic coefficient for the augmented constraint recovery.
    For example, evaluating ``"2*b5"`` with ``{"b5"}`` returns ``2.0``.

    Returns ``None`` if the expression remains symbolic after substitution
    (e.g., ``"kappa*b5"`` still has the un-substituted ``kappa``).
    """
    return evaluate_with_substitutions(expr, dict.fromkeys(one_names, 1.0))


def evaluate_with_substitutions(
    expr: str, substitutions: dict[str, float]
) -> float | None:
    """Return the literal value of ``expr`` after replacing each
    ``name`` → ``substitutions[name]`` in the AST.

    Returns
    -------
    float
        The evaluated literal if every ``Name`` node either matches a
        substitution key or is a numeric literal after recursion.
    None
        If at least one ``Name`` outside ``substitutions`` survives
        — the expression is still symbolic in some other parameter.

    Parameters
    ----------
    expr : str
        Coefficient expression as emitted by ExportJSON.wl (e.g.
        ``"2*b5"``, ``"(-25*b5)/2"``, ``"-kappa^(-2)"``). Caret
        exponents convert to Python ``**`` before parsing.
    substitutions : dict[str, float]
        Name → float map. Only these symbols are replaced; others
        surviving after reduction cause a ``None`` return.

    Raises
    ------
    KineticEvalError
        If ``expr`` contains an AST node kind outside the safe
        subset (literals, names, unary +/-, binary +/-/*/////**).
    """
    if not isinstance(expr, str):
        msg = f"evaluate_with_substitutions: expected a str, got {type(expr).__name__}"
        raise KineticEvalError(msg)

    normalised = expr.replace("^", "**")

    try:
        tree = ast.parse(normalised, mode="eval")
    except SyntaxError as e:
        msg = (
            f"evaluate_with_substitutions: cannot parse {expr!r} "
            f"(after ^→**: {normalised!r}): {e}"
        )
        raise KineticEvalError(msg) from e

    try:
        return _eval_node(tree.body, substitutions)
    except _StillSymbolic:
        return None


class _StillSymbolic(Exception):
    """Signal raised when evaluation encounters a surviving Name."""


def _eval_node(node: ast.AST, substitutions: dict[str, float]) -> float:
    """Recursively evaluate a safe AST node; raise _StillSymbolic if a
    Name outside ``substitutions`` is encountered.
    """
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float)):
            return float(node.value)
        msg = (
            f"evaluate_with_substitutions: unsupported literal type "
            f"{type(node.value).__name__!r}"
        )
        raise KineticEvalError(msg)

    if isinstance(node, ast.Name):
        if node.id in substitutions:
            return substitutions[node.id]
        raise _StillSymbolic

    if isinstance(node, ast.UnaryOp):
        op_fn = _UNARY_OPS.get(type(node.op))
        if op_fn is None:
            msg = (
                f"evaluate_with_substitutions: unsupported unary op "
                f"{type(node.op).__name__}"
            )
            raise KineticEvalError(msg)
        return op_fn(_eval_node(node.operand, substitutions))

    if isinstance(node, ast.BinOp):
        op_fn_bin = _BIN_OPS.get(type(node.op))
        if op_fn_bin is None:
            msg = (
                f"evaluate_with_substitutions: unsupported binary op "
                f"{type(node.op).__name__}"
            )
            raise KineticEvalError(msg)
        left = _eval_node(node.left, substitutions)
        right = _eval_node(node.right, substitutions)
        # 0**0 is mathematically ambiguous; treat as 1 per IEEE, which is
        # also Python's convention. 0**negative is ZeroDivisionError from
        # operator.pow — propagate with a clearer message.
        if isinstance(node.op, ast.Pow) and left == 0.0 and right < 0:
            msg = (
                f"evaluate_with_substitutions: divide-by-zero in expression "
                f"(0 ** {right} undefined)"
            )
            raise KineticEvalError(msg)
        if isinstance(node.op, ast.Div) and right == 0.0:
            msg = "evaluate_with_substitutions: divide-by-zero in expression"
            raise KineticEvalError(msg)
        return op_fn_bin(left, right)

    msg = f"evaluate_with_substitutions: unsupported AST node {type(node).__name__}"
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
