"""User-specifiable parameter constraints for inference.

Constraints are hard priors that reject parameter combinations without
running a simulation.  They are specified as inequality strings::

    --constraint "xi > 0"
    --constraint "alpha > 0"
    --constraint "deltam**2 < 2*alpha*xi"

The parser uses Python's :mod:`ast` module for safe evaluation —
only arithmetic operators and a small set of math functions are allowed.
No arbitrary code execution.
"""

from __future__ import annotations

import ast
import math
import operator
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable

# Safe operators for constraint evaluation
_SAFE_OPS: dict[type, Callable[..., Any]] = {
    ast.Add: operator.add,
    ast.Sub: operator.sub,
    ast.Mult: operator.mul,
    ast.Div: operator.truediv,
    ast.Pow: operator.pow,
    ast.USub: operator.neg,
    ast.UAdd: operator.pos,
}

# Safe comparison operators
_SAFE_CMPS: dict[type, Callable[..., bool]] = {
    ast.Gt: operator.gt,
    ast.GtE: operator.ge,
    ast.Lt: operator.lt,
    ast.LtE: operator.le,
    ast.Eq: operator.eq,
}

# Safe functions available in constraint expressions
_SAFE_FUNCS: dict[str, Callable[..., Any]] = {
    "abs": abs,
    "sqrt": math.sqrt,
    "log": math.log,
    "exp": math.exp,
    "sin": math.sin,
    "cos": math.cos,
}


class ConstraintError(ValueError):
    """Raised when a constraint expression cannot be parsed."""


def _eval_node(node: ast.AST, params: dict[str, float]) -> float:
    """Recursively evaluate an AST node with the given parameters.

    Raises
    ------
    ConstraintError
        If the node cannot be evaluated safely.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, int | float):
        return float(node.value)
    if isinstance(node, ast.Name):
        if node.id in params:
            return params[node.id]
        if node.id in _SAFE_FUNCS:
            msg = f"'{node.id}' is a function, not a parameter — use '{node.id}(...)'"
            raise ConstraintError(msg)
        msg = (
            f"Unknown parameter '{node.id}' in constraint. Available: {sorted(params)}"
        )
        raise ConstraintError(msg)
    if isinstance(node, ast.BinOp):
        bin_fn = _SAFE_OPS.get(type(node.op))
        if bin_fn is None:
            msg = f"Unsupported operator: {type(node.op).__name__}"
            raise ConstraintError(msg)
        left = _eval_node(node.left, params)
        right = _eval_node(node.right, params)
        return float(bin_fn(left, right))
    if isinstance(node, ast.UnaryOp):
        unary_fn = _SAFE_OPS.get(type(node.op))
        if unary_fn is None:
            msg = f"Unsupported unary operator: {type(node.op).__name__}"
            raise ConstraintError(msg)
        return float(unary_fn(_eval_node(node.operand, params)))
    if isinstance(node, ast.Call):
        if not isinstance(node.func, ast.Name):
            msg = "Only simple function calls supported (e.g. sqrt(x))"
            raise ConstraintError(msg)
        func_name = node.func.id
        if func_name not in _SAFE_FUNCS:
            msg = f"Unknown function '{func_name}'. Available: {sorted(_SAFE_FUNCS)}"
            raise ConstraintError(msg)
        if len(node.args) != 1 or node.keywords:
            msg = f"Function '{func_name}' takes exactly 1 positional argument"
            raise ConstraintError(msg)
        arg = _eval_node(node.args[0], params)
        func = _SAFE_FUNCS[func_name]
        return float(func(arg))
    msg = f"Unsupported expression node: {type(node).__name__}"
    raise ConstraintError(msg)


def _eval_comparison(node: ast.Compare, params: dict[str, float]) -> bool:
    """Evaluate a comparison expression like ``xi > 0``.

    Raises
    ------
    ConstraintError
        If the comparison contains unsupported operators.
    """
    left = _eval_node(node.left, params)
    for cmp_op, comparator in zip(node.ops, node.comparators, strict=False):
        compare_fn = _SAFE_CMPS.get(type(cmp_op))
        if compare_fn is None:
            msg = f"Unsupported comparison: {type(cmp_op).__name__}"
            raise ConstraintError(msg)
        right = _eval_node(comparator, params)
        if not compare_fn(left, right):
            return False
        left = right
    return True


def parse_constraint(expr: str) -> Callable[[dict[str, float]], bool]:
    """Parse a constraint expression string into a callable.

    Parameters
    ----------
    expr : str
        A comparison expression, e.g. ``"xi > 0"`` or
        ``"deltam**2 < 2*alpha*xi"``.

    Returns
    -------
    callable
        A function ``(params: dict[str, float]) -> bool`` that returns
        True if the constraint is satisfied.

    Raises
    ------
    ConstraintError
        If the expression cannot be parsed safely.
    """
    try:
        tree = ast.parse(expr, mode="eval")
    except SyntaxError as e:
        msg = f"Invalid constraint syntax: '{expr}': {e}"
        raise ConstraintError(msg) from e

    body = tree.body
    if not isinstance(body, ast.Compare):
        msg = (
            f"Constraint must be a comparison (e.g. 'xi > 0'), "
            f"got {type(body).__name__}: '{expr}'"
        )
        raise ConstraintError(msg)

    def check(params: dict[str, float]) -> bool:
        return _eval_comparison(body, params)

    # Store the original expression for repr
    check.__doc__ = expr
    return check


@dataclass
class ConstraintSet:
    """A collection of parameter constraints combined with logical AND.

    Parameters
    ----------
    constraints : list of callables
        Each callable takes ``dict[str, float]`` and returns ``bool``.
    expressions : list of str
        The original expression strings (for display/serialization).
    """

    constraints: list[Callable[[dict[str, float]], bool]] = field(  # pyright: ignore[reportUnknownVariableType]
        default_factory=list
    )
    expressions: list[str] = field(default_factory=list)  # pyright: ignore[reportUnknownVariableType]

    @classmethod
    def from_strings(cls, exprs: list[str]) -> ConstraintSet:
        """Create from a list of constraint expression strings."""
        constraints = [parse_constraint(expr) for expr in exprs]
        return cls(constraints=constraints, expressions=list(exprs))

    def check(self, params: dict[str, float]) -> bool:
        """Return True if all constraints are satisfied."""
        return all(c(params) for c in self.constraints)

    def __len__(self) -> int:
        return len(self.constraints)

    def __bool__(self) -> bool:
        return len(self.constraints) > 0
