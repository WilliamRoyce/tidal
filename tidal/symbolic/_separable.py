"""Detect BSM-multiplicative separability of RHS coefficient expressions.

For GH #384 Phase A′ (convolution-block cache): determines whether a term's
``coefficient_symbolic`` factors as ``c_BSM(α) · c_geom(coords, geometry)``
with the sampled BSM symbols (α) appearing only as top-level multiplicative
Names. Such terms produce convolution-matrix blocks that are LINEAR in
``c_BSM(α)``; the geometric block is identical across PolyChord likelihood
calls and can be cached at chain start.

Public API
----------
``extract_separable_bsm_factors(coeff_expr, bsm_symbols)``
    Returns ``(c_BSM_expr, c_geom_expr)`` if separable, else ``(None, coeff_expr)``.

Non-separable patterns this rejects:
- BSM inside an exponential or other transcendental: ``Exp[α*x]``, ``Sin[α]``.
- BSM with non-unit exponent: ``α^2``, ``α**n``.
- BSM in a denominator: ``1/(1+α)``.
- BSM inside one summand but not another: ``α + x[]``.
- BSM mixed with coord/geometry inside ``Pow.base`` with non-integer exponent.

Mathematica-form input is normalised via ``^→**`` and ``x[]→x`` before
``ast.parse``, mirroring the convention in ``tidal/symbolic/_kinetic_eval.py``.
"""

from __future__ import annotations

import ast
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


# Strip nullary xCoba coordinate calls like ``x[]`` → ``x`` so ast.parse
# accepts them. Same pattern used in tidal/symbolic/_kinetic_eval.py (#380).
_NULLARY_CALL_RE = re.compile(r"(\b\w+)\s*\[\s*\]")


def _normalize(coeff_expr: str) -> str:
    """Wolfram InputForm → Python-parseable: ``^→**``, strip nullary brackets."""
    out = coeff_expr.replace("^", "**")
    return _NULLARY_CALL_RE.sub(r"\1", out)


def _collect_names(node: ast.AST) -> set[str]:
    """Return every ``Name.id`` appearing anywhere in ``node``."""
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def _split_top_level_mult(
    node: ast.AST,
    bsm_symbols: frozenset[str],
) -> tuple[list[ast.AST], list[ast.AST]] | None:
    """Walk down a top-level Mult/Div chain, partition into (BSM factors, geom factors).

    Returns ``None`` if any BSM symbol appears in a non-separable position
    (inside Pow.exponent, Call.args, Pow.base with non-integer exponent,
    or Div.denominator).

    Tracks an overall sign accumulator: every USub flips it. If the final
    sign is negative, a ``Constant(-1)`` is prepended to the geom factors
    (NOT the BSM factors — a minus sign is not a BSM scalar). Without this,
    expressions like ``-(α·x)`` lose the leading minus and the reconstructed
    block has the wrong sign vs the original.
    """
    bsm_factors: list[ast.AST] = []
    geom_factors: list[ast.AST] = []
    sign_state = [1]  # accumulator flipped by every USub during descent

    def visit(n: ast.AST, in_numerator: bool) -> bool:  # noqa: PLR0912
        """Return True if separable; populate the lists by side-effect."""
        # Mult / Div: descend into both operands. Div numerator: in_numerator;
        # Div denominator: BSM forbidden.
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, ast.Mult):
                return visit(n.left, in_numerator) and visit(n.right, in_numerator)
            if isinstance(n.op, ast.Div):
                # numerator continues; denominator must be BSM-free.
                if not visit(n.left, in_numerator):
                    return False
                denom_names = _collect_names(n.right)
                if denom_names & bsm_symbols:
                    return False
                # Whole denominator is a geom factor (as 1/d → just keep n.right
                # to be re-emitted in geom expression preserving the Div).
                # We treat the entire Div as a single geom factor at this level:
                # NOTE: this path is reached only when left has been visited
                # (and its BSM factors collected). The right side is geom.
                geom_factors.append(n.right)
                # mark to indicate the original was Div, not Mult — encoded by
                # synthesizing 1/right after extraction; handle in reconstruction.
                # For now, the geom expression is rebuilt as
                # "<remaining_left> / <right>". We tag this with a sentinel:
                geom_factors.append(_DIV_MARKER)
                return True
            if isinstance(n.op, ast.Pow):
                # Base may contain BSM only if exponent is integer-literal 1.
                # Pow.right (exponent) MUST be BSM-free.
                exp_names = _collect_names(n.right)
                if exp_names & bsm_symbols:
                    return False
                # Pow.base: if exponent is 1, descend; else BSM in base is
                # non-separable.
                base_names = _collect_names(n.left)
                if base_names & bsm_symbols:
                    if (
                        isinstance(n.right, ast.Constant)
                        and isinstance(n.right.value, int)
                        and n.right.value == 1
                    ):
                        return visit(n.left, in_numerator)
                    # BSM under a non-unit power: non-separable.
                    return False
                # BSM-free base: treat the whole Pow as a single geom factor.
                geom_factors.append(n)
                return True
            # Other BinOp (Add, Sub): handled at top by the caller. If we
            # reach here, the input wasn't a pure Mult/Div tree.
            if isinstance(n.op, (ast.Add, ast.Sub)):
                # Whole sub-expression is a geom factor only if BSM-free.
                if _collect_names(n) & bsm_symbols:
                    return False
                geom_factors.append(n)
                return True
            return False
        if isinstance(n, ast.UnaryOp):
            if isinstance(n.op, ast.USub):
                # -X: flip the overall sign accumulator (resolved at the end
                # by prepending Constant(-1) to geom_factors), then descend.
                # Without this, -(α·x) loses its minus and the reconstructed
                # convolution block has the wrong sign.
                sign_state[0] = -sign_state[0]
                return visit(n.operand, in_numerator)
            if isinstance(n.op, ast.UAdd):
                # +X: no-op; descend.
                return visit(n.operand, in_numerator)
            return False
        if isinstance(n, ast.Name):
            if n.id in bsm_symbols:
                bsm_factors.append(n)
            else:
                geom_factors.append(n)
            return True
        if isinstance(n, ast.Constant):
            geom_factors.append(n)
            return True
        if isinstance(n, ast.Call):
            # Call args MUST be BSM-free for separability.
            for arg in n.args:
                if _collect_names(arg) & bsm_symbols:
                    return False
            geom_factors.append(n)
            return True
        return False

    if not visit(node, in_numerator=True):
        return None
    if sign_state[0] == -1:
        # Prepend -1 as a numeric geom factor so the reconstructed product
        # carries the overall sign accumulated through USub descents.
        geom_factors.insert(0, ast.Constant(value=-1))
    return bsm_factors, geom_factors


# Sentinel marker: when present in geom_factors, the preceding factor was the
# denominator of a Div node and should be emitted as "/ <factor>" rather than
# "* <factor>". This encodes the Mult/Div distinction in the flat factor list.
class _DivMarker:
    __slots__ = ()

    def __repr__(self) -> str:
        return "<DIV>"


_DIV_MARKER = _DivMarker()


def _emit_factor_product(factors: list[ast.AST]) -> str:
    """Reconstruct a multiplicative product from flat factor list.

    Each entry is either an AST node (multiplied in) or ``_DIV_MARKER`` which
    converts the previously-appended factor into a division.

    Output is in Mathematica InputForm-compatible style (``^`` for exponents,
    not Python's ``**``) so the downstream ``evaluate_coefficient`` ->
    ``mathematica_to_python`` pipeline can apply ``_invert_exp_denominator``
    (which keys on ``1/E^(arg)`` patterns to avoid overflow). If we emit
    Python ``**`` directly the InputForm-detection regex misses, and
    ``1/E**(positive)`` overflows at the box edges. See GH #384 Phase A′.
    """
    if not factors:
        return "1"
    parts: list[str] = []
    i = 0
    while i < len(factors):
        f = factors[i]
        # Look ahead for div marker
        op = "*"
        is_div = i + 1 < len(factors) and factors[i + 1] is _DIV_MARKER
        if is_div:
            op = "/"
            i += 1
        if parts:
            parts.append(op)
        rendered = ast.unparse(f)
        # ast.unparse omits outer parens. For a / b where b is itself a
        # Mult/Add/etc., the result "a / b * c" parses as (a/b)*c — wrong
        # precedence. Wrap the denominator if it's not already atomic.
        if is_div and isinstance(f, ast.BinOp):
            rendered = f"({rendered})"
        parts.append(rendered)
        i += 1
    # Trim leading "*"
    if parts and parts[0] == "*":
        parts = parts[1:]
    joined = " ".join(parts) if parts else "1"
    # Convert Python ` ** ` → Mathematica `^` so downstream parsers route
    # through _invert_exp_denominator (key step for E^(positive) overflows).
    # Strip surrounding spaces — mathematica_to_python's `\bE\^` regex requires
    # NO space between E and ^ to convert `E^X` to `exp(X)`. With spaces left
    # in, ast.unparse-style "E ** X" → "E ^ X" wouldn't match and would
    # fall through to the generic ^→** replacement, producing literal `E**X`
    # via numpy.power which overflows for large X.
    return re.sub(r"\s*\*\*\s*", "^", joined)


def extract_separable_bsm_factors(
    coeff_expr: str,
    bsm_symbols: Iterable[str],
) -> tuple[str | None, str]:
    """Detect BSM-multiplicative separability.

    Returns ``(c_BSM_expr, c_geom_expr)`` if the input factors as a top-level
    product of BSM Names and a BSM-free residual; else ``(None, original)``.

    Parameters
    ----------
    coeff_expr
        Mathematica InputForm coefficient string (e.g. from
        ``OperatorTerm.coefficient_symbolic``).
    bsm_symbols
        Names that are sampled per PolyChord call (e.g. ``{"alpha1", "delta1"}``).

    Examples
    --------
    >>> extract_separable_bsm_factors("alpha1*x[]", {"alpha1"})
    ('alpha1', 'x')
    >>> extract_separable_bsm_factors("alpha1**2*x[]", {"alpha1"})  # non-unit power
    (None, 'alpha1**2*x[]')
    >>> extract_separable_bsm_factors("Exp[alpha1*x[]]", {"alpha1"})  # BSM in Call
    (None, 'Exp[alpha1*x[]]')
    >>> extract_separable_bsm_factors("5*Bpeak", {"alpha1"})  # no BSM present
    ('', '5 * Bpeak')
    """
    if not coeff_expr:
        return None, coeff_expr
    bsm_set = frozenset(bsm_symbols)
    normalized = _normalize(coeff_expr)
    try:
        tree = ast.parse(normalized, mode="eval").body
    except SyntaxError:
        return None, coeff_expr

    # Handle Add/Sub at top level: each summand must factor with the SAME BSM
    # product (modulo sign). Otherwise non-separable.
    if isinstance(tree, ast.BinOp) and isinstance(tree.op, (ast.Add, ast.Sub)):
        summands: list[tuple[int, ast.AST]] = []
        _flatten_addsub(tree, +1, summands)
        per_summand: list[tuple[set[str], list[ast.AST], int]] = []
        for sign, summand in summands:
            split = _split_top_level_mult(summand, bsm_set)
            if split is None:
                return None, coeff_expr
            bsm_factors, geom_factors = split
            bsm_names = {f.id for f in bsm_factors if isinstance(f, ast.Name)}
            per_summand.append((bsm_names, geom_factors, sign))
        # All summands must share the same set of BSM Names
        reference_bsm = per_summand[0][0]
        for bsm_names, _, _ in per_summand[1:]:
            if bsm_names != reference_bsm:
                return None, coeff_expr
        bsm_expr = " * ".join(sorted(reference_bsm)) if reference_bsm else ""
        geom_parts: list[str] = []
        for _, geom_factors, sign in per_summand:
            piece = _emit_factor_product(geom_factors)
            geom_parts.append(("-" if sign < 0 else "+") + " (" + piece + ")")
        geom_expr = " ".join(geom_parts).lstrip("+ ").strip()
        return bsm_expr or "1", geom_expr

    # Pure multiplicative tree (no top-level Add/Sub)
    split = _split_top_level_mult(tree, bsm_set)
    if split is None:
        return None, coeff_expr
    bsm_factors, geom_factors = split
    bsm_names = sorted({f.id for f in bsm_factors if isinstance(f, ast.Name)})
    bsm_expr = " * ".join(bsm_names) if bsm_names else ""
    geom_expr = _emit_factor_product(geom_factors)
    return bsm_expr or "1", geom_expr


def _flatten_addsub(
    node: ast.AST,
    sign: int,
    out: list[tuple[int, ast.AST]],
) -> None:
    """Decompose top-level Add/Sub tree into signed summands."""
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, ast.Add):
            _flatten_addsub(node.left, sign, out)
            _flatten_addsub(node.right, sign, out)
            return
        if isinstance(node.op, ast.Sub):
            _flatten_addsub(node.left, sign, out)
            _flatten_addsub(node.right, -sign, out)
            return
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            _flatten_addsub(node.operand, -sign, out)
            return
        if isinstance(node.op, ast.UAdd):
            _flatten_addsub(node.operand, sign, out)
            return
    out.append((sign, node))
