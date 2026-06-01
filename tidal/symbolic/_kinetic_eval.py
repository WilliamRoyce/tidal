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
import re
import sys
from typing import TYPE_CHECKING

if sys.version_info >= (3, 12):
    from typing import override
else:
    from typing_extensions import override

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
    """Return ``expr`` evaluated after substituting every name in ``zero_names`` with ``0.0``.

    Thin wrapper over :func:`evaluate_with_substitutions`; kept for its
    existing callers in ``json_loader.base_spec``.
    """
    return evaluate_with_substitutions(expr, dict.fromkeys(zero_names, 0.0))


def evaluate_at_one(expr: str, one_names: set[str] | frozenset[str]) -> float | None:
    """Return ``expr`` evaluated after substituting every name in ``one_names`` with ``1.0``.

    See #290: used by :class:`PerturbativeSolver` to extract the
    numeric kinetic coefficient for the augmented constraint recovery.
    For example, evaluating ``"2*b5"`` with ``{"b5"}`` returns ``2.0``.

    Returns ``None`` if the expression remains symbolic after substitution
    (e.g., ``"kappa*b5"`` still has the un-substituted ``kappa``).
    """
    return evaluate_with_substitutions(expr, dict.fromkeys(one_names, 1.0))


def evaluate_with_substitutions(
    expr: str,
    substitutions: dict[str, float],
) -> float | None:
    """Return ``expr`` evaluated after replacing each ``name`` → ``substitutions[name]`` in the AST.

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
    if not isinstance(expr, str):  # type: ignore[reportUnnecessaryIsInstance]
        msg = (
            f"evaluate_with_substitutions: expected a str, got {type(expr).__name__!r}"
        )
        raise KineticEvalError(msg)
    # Wolfram InputForm to Python-AST: ^ → **, and strip nullary xCoba
    # coordinate calls like `x[]` (ast.parse rejects empty subscripts).
    # Coords are then bare names; if a caller provides a substitution for
    # them they're evaluated, otherwise the term stays symbolic. See #380.
    normalized = expr.replace("^", "**")
    normalized = re.sub(r"(\b\w+)\s*\[\s*\]", r"\1", normalized)

    try:
        tree = ast.parse(normalized, mode="eval")
    except SyntaxError as e:
        msg = (
            f"evaluate_with_substitutions: cannot parse {expr!r} "
            f"(after ^→**: {normalized!r}): {e}"
        )
        raise KineticEvalError(msg) from e

    try:
        return _eval_node(tree.body, substitutions)
    except _StillSymbolicError:
        return None


class _StillSymbolicError(Exception):
    """Signal raised when evaluation encounters a surviving Name."""


def _eval_node(node: ast.AST, substitutions: dict[str, float]) -> float:  # noqa: C901
    """Recursively evaluate a safe AST node.

    Raises
    ------
    KineticEvalError
        If an unsupported AST node kind is encountered.
    _StillSymbolicError
        If a Name outside ``substitutions`` is encountered.
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
        raise _StillSymbolicError

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


def split_small_parameter_kinetic(
    expr: str,
    small_parameters: Sequence[str],
) -> tuple[str, dict[str, str]]:
    """Split a kinetic coefficient into base (M₀) and per-small-parameter corrections.

    Given an expression like ``"1 - 2*B0^2*rho + 16*B0^2*sigma"`` and
    ``small_parameters = ["rho", "sigma"]``, returns::

        ("1", {"rho": "-2*B0**2", "sigma": "16*B0**2"})

    Used by :meth:`EquationSystem.canonicalize_kinetic_for_perturbation`
    (#301 Phase 3 / #303) to move small-parameter kinetic dependence into
    Pass-1 RHS corrections via the identity
    ``(M₀ + εM₁)⁻¹(K₀ + εK₁) ≈ M₀⁻¹K₀ + ε M₀⁻¹(K₁ − M₁·M₀⁻¹·K₀)``.

    Contract
    --------
    - Input must be a sum of monomials (post-``Expand`` on the Wolfram side;
      :func:`tidal.wolfram.ExportJSON.EquationToJSONMultiField` uses
      ``Expand`` for this reason). Parenthesized sub-sums raise.
    - Each summand contains each small parameter *at most linearly* (exponent
      exactly 1, in the numerator of any division, and not raised to a Power).
      Bilinear terms (``rho*sigma``), quadratic terms (``rho^2``), and
      in-denominator appearances (``1/(1-rho)``) all raise
      :class:`KineticEvalError`. These are gated by the TIDAL perturbative
      contract (#273: ``order_in_eps > 1`` not implemented).
    - A summand containing no small parameter goes into ``M₀``.
    - A summand containing exactly one small parameter ``p`` contributes its
      "coefficient when ``p = 1``" to ``c_p``.

    Parameters
    ----------
    expr
        Kinetic coefficient string as emitted by ExportJSON.wl.
    small_parameters
        The names in ``[perturbation].small_parameters``.

    Returns
    -------
    ``(M0_expr, {param: c_expr})`` where all expressions are Python-syntax
    strings (``**`` for exponentiation; ``^`` is converted internally).
    ``M0_expr`` is ``"0"`` if every summand contained a small parameter.
    The returned dict maps *every* small parameter, including those with a
    zero coefficient (``"0"``) — callers can skip those if desired.

    Raises
    ------
    KineticEvalError
        On unsupported structure as described in the contract.
    """
    if not isinstance(expr, str):  # type: ignore[reportUnnecessaryIsInstance]
        msg = (
            f"split_small_parameter_kinetic: expected a str, "
            f"got {type(expr).__name__!r}"
        )
        raise KineticEvalError(msg)

    small_set = frozenset(small_parameters)
    # Wolfram InputForm conventions to Python-AST: caret to **, and
    # nullary xCoba coordinate calls like `x[]` to bare `x` (ast.parse
    # rejects empty subscripts). The renamed `x` flows through to the
    # synthesized RHS coefficient strings and downstream evaluator;
    # mathematica_to_python is idempotent on bare coord names. See GH #380.
    normalized = expr.replace("^", "**")
    normalized = re.sub(r"(\b\w+)\s*\[\s*\]", r"\1", normalized)
    try:
        tree = ast.parse(normalized, mode="eval")
    except SyntaxError as e:
        msg = (
            f"split_small_parameter_kinetic: cannot parse {expr!r} "
            f"(after ^→**: {normalized!r}): {e}"
        )
        raise KineticEvalError(msg) from e

    summands: list[tuple[int, ast.AST]] = []  # (sign, AST)
    _flatten_summands(tree.body, +1, summands)

    m0_terms: list[ast.AST] = []
    coefficients: dict[str, list[ast.AST]] = {name: [] for name in small_set}

    for sign, summand in summands:
        small_params_seen = _count_small_parameters(summand, small_set)
        if sum(small_params_seen.values()) == 0:
            m0_terms.append(_apply_sign(summand, sign))
            continue
        if len(small_params_seen) > 1 or max(small_params_seen.values()) > 1:
            msg = (
                f"split_small_parameter_kinetic: summand in {expr!r} has "
                f"multiple or higher-order small-parameter dependence "
                f"{dict(small_params_seen)!r}; order_in_eps > 1 is gated by "
                "#273 (not implemented). Emit a linear Lagrangian or split "
                "into lower-order pieces."
            )
            raise KineticEvalError(msg)
        (param_name,) = small_params_seen.keys()
        coeff_ast = _substitute_names(summand, {param_name: 1.0})
        coefficients[param_name].append(_apply_sign(coeff_ast, sign))

    m0_expr = _terms_to_expr(m0_terms) if m0_terms else "0"
    coefficient_exprs = {
        name: _terms_to_expr(coefficients[name]) if coefficients[name] else "0"
        for name in small_set
    }
    return m0_expr, coefficient_exprs


def _flatten_summands(
    node: ast.AST,
    sign: int,
    out: list[tuple[int, ast.AST]],
) -> None:
    """Walk top-level Add/Sub tree, emit each monomial with a ±1 sign."""
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Add):
            _flatten_summands(node.left, sign, out)
            _flatten_summands(node.right, sign, out)
            return
        if isinstance(node.op, ast.Sub):
            _flatten_summands(node.left, sign, out)
            _flatten_summands(node.right, -sign, out)
            return
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            _flatten_summands(node.operand, -sign, out)
            return
        if isinstance(node.op, ast.UAdd):
            _flatten_summands(node.operand, sign, out)
            return
    # Leaf or multiplicative/pow/div subtree — one summand.
    _reject_inner_addition(node)
    out.append((sign, node))


def _reject_inner_addition(node: ast.AST) -> None:
    """Raise if any Add/Sub BinOp appears in *multiplicative* position inside a summand.

    ``Expand`` on the Wolfram side distributes ``(a+b)*c`` and ``(a+b)/c`` and
    ``(a+b)^2`` into monomials, so a surviving Add/Sub in those positions
    means input was not properly Expand'd. But Add/Sub inside a ``Pow.right``
    (e.g. ``E**((zc1-x)**2/sigB**2)``) or a ``Call`` argument is transcendental
    — Expand doesn't reduce it — and represents a coord-dependent atomic
    factor. Walk multiplicative structure only; treat Pow exponents and Call
    args as opaque. See GH #380.
    """

    def walk(n: ast.AST) -> None:
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Add, ast.Sub)):
                msg = (
                    "split_small_parameter_kinetic: parenthesized sub-sum "
                    f"{ast.unparse(n)!r} inside a monomial is not supported "
                    "(expected Expand'd input)."
                )
                raise KineticEvalError(msg)
            if isinstance(n.op, ast.Pow):
                # Base is in multiplicative position; exponent is opaque.
                walk(n.left)
                return
            # Mult / Div / etc.: both sides multiplicative.
            walk(n.left)
            walk(n.right)
            return
        if isinstance(n, ast.UnaryOp):
            walk(n.operand)
            return
        # Call args, Name, Constant: opaque/leaf — no further descent.
        return

    walk(node)


def _count_small_parameters(
    node: ast.AST,
    small_set: frozenset[str],
) -> dict[str, int]:
    """Return the occurrence count of each small parameter in *node*.

    Also enforces exponent==1 and numerator-only constraints: any small
    parameter inside ``Pow(base=p, exp != 1)`` or ``Div(num, denom)`` where
    ``p`` is in ``denom`` raises ``KineticEvalError``.
    """
    counts: dict[str, int] = {}
    for child in ast.walk(node):
        if not isinstance(child, ast.Name):
            continue
        if child.id not in small_set:
            continue
        counts[child.id] = counts.get(child.id, 0) + 1

    # Validate: each small-parameter Name must be a direct factor, not a Pow
    # base with exponent != 1 nor inside a Div denominator.
    for child in ast.walk(node):
        if isinstance(child, ast.BinOp):
            if isinstance(child.op, ast.Pow):
                _reject_small_param_in_pow_base(child, small_set)
            elif isinstance(child.op, ast.Div):
                _reject_small_param_in_denominator(child.right, small_set)
    return counts


def _reject_small_param_in_pow_base(
    pow_node: ast.BinOp,
    small_set: frozenset[str],
) -> None:
    """Reject ``p^n`` for n != 1 when ``p`` is a small parameter."""
    for child in ast.walk(pow_node.left):
        if isinstance(child, ast.Name) and child.id in small_set:
            exp_is_one = (
                isinstance(pow_node.right, ast.Constant) and pow_node.right.value == 1
            )
            if not exp_is_one:
                msg = (
                    f"split_small_parameter_kinetic: small parameter "
                    f"{child.id!r} appears with non-unit exponent "
                    f"{ast.unparse(pow_node.right)!r} (order_in_eps > 1, #273)."
                )
                raise KineticEvalError(msg)


def _reject_small_param_in_denominator(
    denom: ast.AST,
    small_set: frozenset[str],
) -> None:
    """Reject rational forms with a small parameter in the denominator."""
    for child in ast.walk(denom):
        if isinstance(child, ast.Name) and child.id in small_set:
            msg = (
                f"split_small_parameter_kinetic: small parameter "
                f"{child.id!r} appears in a denominator — rational dependence "
                "on a small parameter is outside the perturbative contract."
            )
            raise KineticEvalError(msg)


def _substitute_names(node: ast.AST, subs: dict[str, float]) -> ast.AST:
    """Return a new AST with ``Name(id=k)`` replaced by ``Constant(v)`` for ``k, v in subs``.

    Used to extract a small parameter's coefficient by setting the parameter
    to 1.0 in a monomial.
    """

    class _Sub(ast.NodeTransformer):
        @override
        def visit_Name(self, node: ast.Name) -> ast.AST:
            if node.id in subs:
                return ast.copy_location(ast.Constant(value=subs[node.id]), node)
            return node

    return _Sub().visit(
        ast.fix_missing_locations(ast.parse(ast.unparse(node), mode="eval").body)
    )


def _apply_sign(node: ast.AST, sign: int) -> ast.AST:
    """Return ``node`` if ``sign == +1`` else ``-node`` at the AST level."""
    if sign == 1:
        return node
    return ast.copy_location(
        ast.UnaryOp(op=ast.USub(), operand=node),  # pyright: ignore[reportArgumentType]
        node,
    )


def _terms_to_expr(terms: list[ast.AST]) -> str:
    """Join a list of AST nodes into a single sum expression string.

    Uses :func:`ast.unparse`; emits Python syntax (``**``, ``-x``). The
    consumer :func:`evaluate_coefficient` accepts both ``^`` and ``**``.
    """
    if not terms:
        return "0"
    if len(terms) == 1:
        return ast.unparse(terms[0])
    expr: ast.AST = terms[0]
    for t in terms[1:]:
        # Collapse ``+ (-x)`` into ``- x`` for readability.
        if isinstance(t, ast.UnaryOp) and isinstance(t.op, ast.USub):
            expr = ast.BinOp(left=expr, op=ast.Sub(), right=t.operand)  # pyright: ignore[reportArgumentType]
        else:
            expr = ast.BinOp(left=expr, op=ast.Add(), right=t)  # pyright: ignore[reportArgumentType]
    return ast.unparse(expr)


def lhs_collapses_to_zero(
    kinetic_symbolic: str | None,
    small_parameters: Sequence[str],
) -> bool:
    """Return True iff the LHS kinetic coefficient evaluates to literal 0 when every small parameter is set to zero.

    ``None`` kinetic coefficients are implicitly non-perturbative
    (e.g. ``kinetic == 1`` by convention); they never trigger demotion.
    """
    if kinetic_symbolic is None or not small_parameters:
        return False
    zero_names = frozenset(small_parameters)
    value = evaluate_at_zero(kinetic_symbolic, zero_names)
    return value == 0.0
