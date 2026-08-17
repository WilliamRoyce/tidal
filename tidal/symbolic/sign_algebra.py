"""Sound sign and ratio decisions for symbolic coefficient strings.

Coefficients exported by ``ExportJSON.wl`` are strings like ``-xi``,
``-kappa^(-2)``, ``B0^2/2`` or ``-1 + 2*B0^2*rho``.  Asking a physics
question about a spec — "is this component's laplacian the same sign as
its siblings?", "did re-derivation change this equation?" — means deciding
the **sign** or the **ratio** of such expressions, usually with no numeric
values for the free parameters.

The governing requirement is **soundness**: a definite answer is returned
only when it is *proven*.  Anything else is :data:`Sign.UNKNOWN`.  Callers
may then escalate (supply parameters, declare assumptions) or report that
they cannot tell — but they never receive a guess.  Every confident-but-wrong
diagnosis recorded in GH #401 was a case where a proof or an honest UNKNOWN
was available; see the module tests, which pin one case per misreading.

Two independent decision procedures cooperate:

1. **Rational normal form** (:class:`_Ratio`) — expressions are expanded
   into a quotient of Laurent polynomials over :class:`~fractions.Fraction`
   coefficients, with non-polynomial subexpressions (``E**u``, coordinate
   symbols) held as opaque atoms.  Two expressions whose quotient is a
   rational *constant* are then decided exactly, by cross-multiplication —
   never by division.  This is what makes an overall rescaling of an
   equation provably invisible.
2. **Sign lattice** (:class:`Sign`) — the classical sign domain from
   abstract interpretation, ``{0, +, −, 0+, 0−, ⊤}``.  It decides the
   residue that normal form cannot: ``kappa^(-2) > 0``, ``E^u > 0``,
   sums of same-signed summands.

``NONNEGATIVE`` is deliberately distinct from ``POSITIVE``.  ``kappa^2`` is
zero at ``kappa = 0``, and a zero kinetic coefficient means a field is
*constrained rather than dynamical* (see
:attr:`~tidal.symbolic.json_loader.LHSStructure.kinetic_coefficient_symbolic`),
so collapsing the two would let a caller claim a field is dynamical when it
may not be.

Parsing reuses :func:`tidal.symbolic._kinetic_eval.normalize_inputform` and
admits only the restricted node set that module already allows — literals,
names, unary ``±``, and ``+ - * / **``.  Nothing is ever ``eval``-ed.

References
----------
Cousot & Cousot (1977), *Abstract interpretation: a unified lattice model*,
POPL — the sign domain and the role of ⊤ as "unknown".
"""

from __future__ import annotations

import ast
import math
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from enum import Enum
from fractions import Fraction
from typing import TYPE_CHECKING

from tidal.symbolic._kinetic_eval import (
    KineticEvalError,
    evaluate_with_substitutions,
    normalize_inputform,
)

if TYPE_CHECKING:
    from collections.abc import Iterable, Mapping

__all__ = [
    "Sign",
    "SignResult",
    "are_equal",
    "canonical_form",
    "constant_ratio",
    "evaluate_numeric",
    "free_names",
    "ratio_sign",
    "sign_of",
]


# --- Sign lattice ---


class Sign(Enum):
    """Element of the sign domain ``{0, +, −, 0+, 0−, ⊤}``.

    ``NONNEGATIVE``/``NONPOSITIVE`` are distinct from ``POSITIVE``/``NEGATIVE``
    because the zero case is physically meaningful here: a vanishing kinetic
    coefficient turns a dynamical field into a constrained one.
    """

    ZERO = "zero"
    POSITIVE = "positive"
    NEGATIVE = "negative"
    NONNEGATIVE = "nonnegative"
    NONPOSITIVE = "nonpositive"
    UNKNOWN = "unknown"

    @property
    def is_definite(self) -> bool:
        """Whether this element pins the sign exactly (``0``, ``+`` or ``−``)."""
        return self in {Sign.ZERO, Sign.POSITIVE, Sign.NEGATIVE}

    @property
    def symbol(self) -> str:
        """Short display form, e.g. ``"+"``, ``"0+"``, ``"?"``."""
        return _SIGN_SYMBOLS[self]


_SIGN_SYMBOLS: dict[Sign, str] = {
    Sign.ZERO: "0",
    Sign.POSITIVE: "+",
    Sign.NEGATIVE: "-",
    Sign.NONNEGATIVE: "0+",
    Sign.NONPOSITIVE: "0-",
    Sign.UNKNOWN: "?",
}

_NEGATE: dict[Sign, Sign] = {
    Sign.ZERO: Sign.ZERO,
    Sign.POSITIVE: Sign.NEGATIVE,
    Sign.NEGATIVE: Sign.POSITIVE,
    Sign.NONNEGATIVE: Sign.NONPOSITIVE,
    Sign.NONPOSITIVE: Sign.NONNEGATIVE,
    Sign.UNKNOWN: Sign.UNKNOWN,
}

# Strictness of each element, used to build the multiplication table.
_STRICT = {Sign.POSITIVE: 1, Sign.NEGATIVE: -1}
_WEAK = {Sign.NONNEGATIVE: 1, Sign.NONPOSITIVE: -1}


def _sign_from_number(value: float | Fraction) -> Sign:
    """Return the exact lattice element for a numeric literal."""
    if value > 0:
        return Sign.POSITIVE
    if value < 0:
        return Sign.NEGATIVE
    return Sign.ZERO


def _sign_mul(left: Sign, right: Sign) -> Sign:
    """Return the lattice product ``left * right`` (sound, not necessarily tight)."""
    if left is Sign.ZERO or right is Sign.ZERO:
        return Sign.ZERO
    if left is Sign.UNKNOWN or right is Sign.UNKNOWN:
        return Sign.UNKNOWN
    # Both are strictly or weakly signed: multiply directions, and the result
    # is strict only when both inputs are strict.
    left_dir = _STRICT.get(left) or _WEAK[left]
    right_dir = _STRICT.get(right) or _WEAK[right]
    direction = left_dir * right_dir
    strict = left in _STRICT and right in _STRICT
    if strict:
        return Sign.POSITIVE if direction > 0 else Sign.NEGATIVE
    return Sign.NONNEGATIVE if direction > 0 else Sign.NONPOSITIVE


def _sign_add(left: Sign, right: Sign) -> Sign:
    """Return the lattice sum ``left + right``.

    Opposite directions admit cancellation, so they yield ``UNKNOWN`` — this
    is why ``-1 + 2*B0^2*rho`` is correctly undecided: it flips sign at
    ``B0^2*rho = 1/2``.
    """
    if left is Sign.ZERO:
        return right
    if right is Sign.ZERO:
        return left
    if left is Sign.UNKNOWN or right is Sign.UNKNOWN:
        return Sign.UNKNOWN
    left_dir = _STRICT.get(left) or _WEAK[left]
    right_dir = _STRICT.get(right) or _WEAK[right]
    if left_dir != right_dir:
        return Sign.UNKNOWN
    # Same direction: strict if either summand is strict.
    if left in _STRICT or right in _STRICT:
        return Sign.POSITIVE if left_dir > 0 else Sign.NEGATIVE
    return Sign.NONNEGATIVE if left_dir > 0 else Sign.NONPOSITIVE


def _constant_number(node: ast.Constant) -> float:
    """Return the numeric value of a constant node.

    ``_reject_unsafe`` already guarantees numeric literals; this narrows the
    type for the checker and fails loudly if that invariant is ever broken.

    Raises
    ------
    KineticEvalError
        If the literal is not numeric.
    """
    value = node.value
    if isinstance(value, (int, float)):
        return float(value)
    msg = f"sign_algebra: unsupported literal {value!r}"
    raise KineticEvalError(msg)


def _refine_nonzero(sign: Sign) -> Sign:
    """Sharpen *sign* given that the value is known to be non-zero.

    Used for denominators and negative exponents: for the expression to be
    defined at all, the base cannot vanish.
    """
    if sign is Sign.NONNEGATIVE:
        return Sign.POSITIVE
    if sign is Sign.NONPOSITIVE:
        return Sign.NEGATIVE
    if sign is Sign.ZERO:
        # A zero denominator is a malformed expression; report UNKNOWN rather
        # than inventing a value.
        return Sign.UNKNOWN
    return sign


# --- Parsing ---


def _parse(expr: str | float | None) -> ast.expr:
    """Parse a coefficient into a restricted AST.

    Parameters
    ----------
    expr : str | float | None
        Wolfram InputForm string, a plain number, or ``None`` (treated as ``1``,
        matching the convention that an absent kinetic coefficient is unity).

    Returns
    -------
    ast.expr
        The parsed expression node.

    Raises
    ------
    KineticEvalError
        If the source cannot be parsed or uses a node outside the safe subset.
    """
    if expr is None:
        return ast.parse("1", mode="eval").body
    if isinstance(expr, (int, float, Fraction)):
        return ast.parse(repr(float(expr)), mode="eval").body
    normalized = normalize_inputform(expr)
    try:
        tree = ast.parse(normalized, mode="eval")
    except SyntaxError as exc:
        msg = f"sign_algebra: cannot parse {expr!r} (normalized: {normalized!r}): {exc}"
        raise KineticEvalError(msg) from exc
    _reject_unsafe(tree.body)
    return tree.body


_SAFE_BINOPS = (ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow)
_SAFE_UNARYOPS = (ast.UAdd, ast.USub)


def _reject_unsafe(node: ast.expr) -> None:
    """Raise if *node* uses anything outside the restricted subset.

    Raises
    ------
    KineticEvalError
        On any unsupported node or operator kind.
    """
    for child in ast.walk(node):
        if isinstance(child, ast.Constant):
            if not isinstance(child.value, (int, float)):
                msg = f"sign_algebra: unsupported literal {child.value!r}"
                raise KineticEvalError(msg)
        elif isinstance(child, ast.BinOp):
            if not isinstance(child.op, _SAFE_BINOPS):
                msg = f"sign_algebra: unsupported operator {type(child.op).__name__}"
                raise KineticEvalError(msg)
        elif isinstance(child, ast.UnaryOp):
            if not isinstance(child.op, _SAFE_UNARYOPS):
                msg = f"sign_algebra: unsupported unary {type(child.op).__name__}"
                raise KineticEvalError(msg)
        elif not isinstance(
            child,
            # ast.walk yields the operator nodes themselves; the BinOp/UnaryOp
            # branches above already validated which operators are allowed.
            (ast.Expression, ast.Name, ast.Load, ast.operator, ast.unaryop),
        ):
            msg = f"sign_algebra: unsupported node {type(child).__name__}"
            raise KineticEvalError(msg)


def free_names(expr: str | float | None) -> tuple[str, ...]:
    """Return the sorted free symbol names appearing in *expr*.

    ``E`` is excluded — it denotes Euler's number, not a parameter.

    Parameters
    ----------
    expr : str | float | None
        Coefficient expression.

    Returns
    -------
    tuple[str, ...]
        Sorted distinct parameter names.
    """
    node = _parse(expr)
    return tuple(
        sorted({n.id for n in ast.walk(node) if isinstance(n, ast.Name)} - {"E"}),
    )


# --- Rational normal form ---

# A monomial maps atom name -> integer exponent (possibly negative).  Atoms are
# either bare parameter names or the canonical string of an opaque subexpression
# such as ``E**(-x**2)``.
_Monomial = tuple[tuple[str, int], ...]
# A Laurent polynomial maps monomial -> rational coefficient.
_Poly = dict[_Monomial, Fraction]

_ONE_MONOMIAL: _Monomial = ()


def _poly_const(value: Fraction) -> _Poly:
    """Return the constant polynomial *value*."""
    return {} if value == 0 else {_ONE_MONOMIAL: value}


def _poly_atom(name: str, exponent: int = 1) -> _Poly:
    """Return the single-atom monomial ``name**exponent``."""
    if exponent == 0:
        return _poly_const(Fraction(1))
    return {((name, exponent),): Fraction(1)}


def _mono_mul(left: _Monomial, right: _Monomial) -> _Monomial:
    """Multiply two monomials, dropping atoms whose exponents cancel."""
    merged: dict[str, int] = dict(left)
    for name, exponent in right:
        total = merged.get(name, 0) + exponent
        if total:
            merged[name] = total
        else:
            merged.pop(name, None)
    return tuple(sorted(merged.items()))


def _poly_add(left: _Poly, right: _Poly) -> _Poly:
    """Return ``left + right``."""
    out: _Poly = dict(left)
    for mono, coeff in right.items():
        total = out.get(mono, Fraction(0)) + coeff
        if total:
            out[mono] = total
        else:
            out.pop(mono, None)
    return out


def _poly_neg(poly: _Poly) -> _Poly:
    """Return ``-poly``."""
    return {mono: -coeff for mono, coeff in poly.items()}


def _poly_mul(left: _Poly, right: _Poly) -> _Poly:
    """Return ``left * right`` by distributing."""
    out: _Poly = {}
    for lmono, lcoeff in left.items():
        for rmono, rcoeff in right.items():
            mono = _mono_mul(lmono, rmono)
            total = out.get(mono, Fraction(0)) + lcoeff * rcoeff
            if total:
                out[mono] = total
            else:
                out.pop(mono, None)
    return out


def _poly_is_monomial(poly: _Poly) -> bool:
    """Whether *poly* has exactly one term (so it can be inverted exactly)."""
    return len(poly) == 1


def _poly_invert_monomial(poly: _Poly) -> _Poly:
    """Return ``1 / poly`` for a single-term *poly*."""
    ((mono, coeff),) = poly.items()
    return {tuple((name, -exp) for name, exp in mono): 1 / coeff}


def _poly_pow(poly: _Poly, exponent: int) -> _Poly:
    """Raise *poly* to an integer power, inverting monomials when negative.

    Raises
    ------
    _NotPolynomialError
        If a negative power is applied to a multi-term polynomial, since
        ``1/(a+b)`` is not a Laurent polynomial.
    """
    if exponent == 0:
        return _poly_const(Fraction(1))
    if exponent < 0:
        if not _poly_is_monomial(poly):
            # 1/(a+b) is not a Laurent polynomial; caller falls back to _Ratio.
            raise _NotPolynomialError
        return _poly_pow(_poly_invert_monomial(poly), -exponent)
    out = _poly_const(Fraction(1))
    for _ in range(exponent):
        out = _poly_mul(out, poly)
    return out


class _NotPolynomialError(Exception):
    """Signal that a subexpression cannot be expanded as a Laurent polynomial."""


@dataclass(frozen=True)
class _Ratio:
    """An expression as ``num / den``, both Laurent polynomials."""

    num: _Poly
    den: _Poly

    def is_zero(self) -> bool:
        """Whether the numerator is identically zero."""
        return not self.num


def _to_ratio(node: ast.expr) -> _Ratio:
    """Expand *node* into a quotient of Laurent polynomials.

    Subexpressions that are not rational in the parameters (``E**u``,
    non-integer powers) become opaque atoms keyed by their canonical string,
    which keeps the representation exact: two occurrences of the same opaque
    subexpression share an atom and therefore cancel against each other.

    Raises
    ------
    _NotPolynomialError
        If the node cannot be expanded as a quotient of Laurent polynomials.
    """
    if isinstance(node, ast.Constant):
        return _Ratio(
            _poly_const(Fraction(_constant_number(node)).limit_denominator()),
            _poly_const(Fraction(1)),
        )

    if isinstance(node, ast.Name):
        return _Ratio(_poly_atom(node.id), _poly_const(Fraction(1)))

    if isinstance(node, ast.UnaryOp):
        inner = _to_ratio(node.operand)
        if isinstance(node.op, ast.USub):
            return _Ratio(_poly_neg(inner.num), inner.den)
        return inner

    if isinstance(node, ast.BinOp):
        return _binop_to_ratio(node)

    raise _NotPolynomialError


def _binop_to_ratio(node: ast.BinOp) -> _Ratio:
    """Expand a binary operation into a :class:`_Ratio`.

    Raises
    ------
    _NotPolynomialError
        If the node cannot be expanded as a quotient of Laurent polynomials.
    """
    if isinstance(node.op, ast.Pow):
        return _pow_to_ratio(node)

    left = _to_ratio(node.left)
    right = _to_ratio(node.right)

    if isinstance(node.op, ast.Add):
        return _Ratio(
            _poly_add(_poly_mul(left.num, right.den), _poly_mul(right.num, left.den)),
            _poly_mul(left.den, right.den),
        )
    if isinstance(node.op, ast.Sub):
        return _Ratio(
            _poly_add(
                _poly_mul(left.num, right.den),
                _poly_neg(_poly_mul(right.num, left.den)),
            ),
            _poly_mul(left.den, right.den),
        )
    if isinstance(node.op, ast.Mult):
        return _Ratio(
            _poly_mul(left.num, right.num),
            _poly_mul(left.den, right.den),
        )
    # Div
    if right.is_zero():
        raise _NotPolynomialError
    return _Ratio(
        _poly_mul(left.num, right.den),
        _poly_mul(left.den, right.num),
    )


def _pow_to_ratio(node: ast.BinOp) -> _Ratio:
    """Expand ``base ** exponent`` when the exponent is an integer literal."""
    exponent = _integer_literal(node.right)
    if exponent is None:
        # Symbolic or fractional exponent (e.g. `E**(-x**2/s**2)`): keep the
        # whole node as one opaque atom. The key must come from the structural
        # renderer, not canonical_form -- the latter re-enters _to_ratio here
        # and would recurse forever.
        return _Ratio(_poly_atom(_structural_form(node)), _poly_const(Fraction(1)))
    base = _to_ratio(node.left)
    if exponent >= 0:
        return _Ratio(_poly_pow(base.num, exponent), _poly_pow(base.den, exponent))
    return _Ratio(_poly_pow(base.den, -exponent), _poly_pow(base.num, -exponent))


def _integer_literal(node: ast.expr) -> int | None:
    """Return the integer value of *node* if it is an integer literal, else None."""
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return node.value
    if isinstance(node, ast.Constant) and isinstance(node.value, float):
        return int(node.value) if node.value.is_integer() else None
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        inner = _integer_literal(node.operand)
        return None if inner is None else -inner
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
        return _integer_literal(node.operand)
    return None


def _mono_str(mono: _Monomial) -> str:
    """Render a monomial in canonical form."""
    if not mono:
        return "1"
    return "*".join(f"{name}^{exp}" if exp != 1 else name for name, exp in mono)


def _poly_str(poly: _Poly) -> str:
    """Render a polynomial in canonical form (sorted, explicit coefficients)."""
    if not poly:
        return "0"
    parts = [
        f"{coeff}*{_mono_str(mono)}" for mono, coeff in sorted(poly.items(), key=repr)
    ]
    return " + ".join(parts)


def canonical_form(expr: str | float | ast.expr | None) -> str:
    """Return a canonical string for *expr*, stable under rewriting.

    Two expressions that differ only in how they were written — ``3*chi - xi``
    versus ``(3*chi) + (-xi)``, or a reordered product — produce the same
    string.  Expressions that cannot be expanded rationally fall back to a
    structural rendering, which is still stable but only decides syntactic
    equality.

    Parameters
    ----------
    expr : str | float | ast.expr | None
        Coefficient expression, or an already-parsed node.

    Returns
    -------
    str
        Canonical representation.
    """
    node = expr if isinstance(expr, ast.expr) else _parse(expr)
    try:
        ratio = _to_ratio(node)
    except (_NotPolynomialError, KineticEvalError):
        return _structural_form(node)
    return f"({_poly_str(ratio.num)}) / ({_poly_str(ratio.den)})"


def _structural_form(node: ast.expr) -> str:  # noqa: PLR0911
    """Render *node* structurally, sorting commutative operands."""
    if isinstance(node, ast.Constant):
        return repr(_constant_number(node))
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.UnaryOp):
        inner = _structural_form(node.operand)
        return f"(neg {inner})" if isinstance(node.op, ast.USub) else inner
    if isinstance(node, ast.BinOp):
        left = _structural_form(node.left)
        right = _structural_form(node.right)
        if isinstance(node.op, ast.Add):
            return "(+ " + " ".join(sorted([left, right])) + ")"
        if isinstance(node.op, ast.Sub):
            return "(+ " + " ".join(sorted([left, f"(neg {right})"])) + ")"
        if isinstance(node.op, ast.Mult):
            return "(* " + " ".join(sorted([left, right])) + ")"
        if isinstance(node.op, ast.Div):
            return f"(/ {left} {right})"
        return f"(^ {left} {right})"
    return ast.dump(node)


# --- Sign evaluation over the lattice ---


@dataclass(frozen=True)
class _Assumptions:
    """Caller-declared facts about parameters. Empty by default — nothing is implicit.

    Two independent declarations, because they are genuinely different claims:

    ``positive``
        The parameter is strictly greater than zero.
    ``nonzero``
        The parameter never vanishes, but its sign is not claimed.  This is the
        common physical case: ``kappa`` cannot be zero (that would delete the
        Einstein-Hilbert term altogether), so ``kappa^2 > 0`` strictly — while
        parameters that genuinely *do* reach zero in this project, such as
        ``xi`` (at which torsion becomes constrained) or ``b5`` (see the
        ``torsion_gertsenshtein_b5_zero`` example), must keep the cautious
        ``0+`` reading unless the caller says otherwise.
    """

    positive: frozenset[str] = frozenset()
    nonzero: frozenset[str] = frozenset()

    def names_used(self, names: Iterable[str]) -> tuple[str, ...]:
        """Return the declared names that actually occur in *names*."""
        pool = set(names)
        return tuple(sorted((self.positive | self.nonzero) & pool))

    def is_empty(self) -> bool:
        """Whether no assumption was declared at all."""
        return not (self.positive or self.nonzero)


_NO_ASSUMPTIONS = _Assumptions()


def _is_nonzero(node: ast.expr, assumptions: _Assumptions) -> bool:
    """Whether *node* is provably non-zero.

    Either its sign is strictly known, or it is built multiplicatively from
    factors the caller declared non-zero.  Sums are never assumed non-zero:
    ``a - a`` vanishes.
    """
    if _sign_node(node, assumptions) in {Sign.POSITIVE, Sign.NEGATIVE}:
        return True
    if isinstance(node, ast.Name):
        return node.id in assumptions.nonzero or node.id in assumptions.positive
    if isinstance(node, ast.UnaryOp):
        return _is_nonzero(node.operand, assumptions)
    if isinstance(node, ast.BinOp):
        if isinstance(node.op, (ast.Mult, ast.Div)):
            return _is_nonzero(node.left, assumptions) and _is_nonzero(
                node.right,
                assumptions,
            )
        if isinstance(node.op, ast.Pow):
            return _is_nonzero(node.left, assumptions)
    return False


def _sign_node(node: ast.expr, assumptions: _Assumptions) -> Sign:  # noqa: PLR0911
    """Return the lattice sign of *node*, proven only.

    A flat dispatch over AST node kinds; one return per kind reads better here
    than threading a result variable through the branches.
    """
    if isinstance(node, ast.Constant):
        return _sign_from_number(_constant_number(node))

    if isinstance(node, ast.Name):
        if node.id == "E":  # Euler's number
            return Sign.POSITIVE
        if node.id in assumptions.positive:
            return Sign.POSITIVE
        return Sign.UNKNOWN

    if isinstance(node, ast.UnaryOp):
        inner = _sign_node(node.operand, assumptions)
        return _NEGATE[inner] if isinstance(node.op, ast.USub) else inner

    if isinstance(node, ast.BinOp):
        return _sign_binop(node, assumptions)

    return Sign.UNKNOWN


def _sign_binop(node: ast.BinOp, assumptions: _Assumptions) -> Sign:
    """Return the lattice sign of a binary operation."""
    if isinstance(node.op, ast.Pow):
        return _sign_pow(node, assumptions)

    left = _sign_node(node.left, assumptions)
    right = _sign_node(node.right, assumptions)

    if isinstance(node.op, ast.Add):
        return _sign_add(left, right)
    if isinstance(node.op, ast.Sub):
        return _sign_add(left, _NEGATE[right])
    if isinstance(node.op, ast.Mult):
        return _sign_mul(left, right)
    # Div: the expression being defined implies a non-zero denominator.
    return _sign_mul(left, _refine_nonzero(right))


def _sign_pow(node: ast.BinOp, assumptions: _Assumptions) -> Sign:
    """Return the lattice sign of ``base ** exponent``.

    A strictly positive base stays positive under any real exponent.  An even
    exponent yields a non-negative value; it sharpens to strictly positive when
    the base cannot vanish — either because the exponent is negative (so the
    expression is only defined for a non-zero base) or because the caller
    declared the base non-zero.
    """
    base = _sign_node(node.left, assumptions)
    if base is Sign.POSITIVE:
        return Sign.POSITIVE

    exponent = _integer_literal(node.right)
    if exponent is None:
        return Sign.UNKNOWN

    # A negative exponent is only defined for a non-zero base, so it proves
    # non-vanishing on its own; otherwise the caller must have declared it.
    base_nonzero = exponent < 0 or _is_nonzero(node.left, assumptions)

    if exponent % 2 == 0:
        return Sign.POSITIVE if base_nonzero else Sign.NONNEGATIVE
    # Odd power preserves the base's sign.
    return _refine_nonzero(base) if base_nonzero else base


# --- Public results ---


@dataclass(frozen=True)
class SignResult:
    """Outcome of a sign or ratio query, with the reasoning that produced it.

    Attributes
    ----------
    sign : Sign
        The lattice verdict.  ``Sign.UNKNOWN`` means *not proven*, never
        "probably zero" or "probably positive".
    tactic : str
        Which rung of the ladder decided it (``"literal"``, ``"normal-form"``,
        ``"lattice"``, ``"assumption"``, ``"numeric"``) or ``"undecided"``.
    free_names : tuple[str, ...]
        Unresolved parameter names, for reporting when undecided.
    assumptions : tuple[str, ...]
        Caller-declared assumptions that were actually used.  Never implicit,
        and always surfaced so a reader can audit what the verdict rests on.
    value : Fraction | None
        The exact rational value when the query reduced to a constant.
    numeric : float | None
        Corroborating numeric evaluation, reported separately and never used
        to justify ``sign``.
    """

    sign: Sign
    tactic: str
    free_names: tuple[str, ...] = ()
    assumptions: tuple[str, ...] = dataclass_field(default_factory=tuple)
    value: Fraction | None = None
    numeric: float | None = None

    @property
    def is_definite(self) -> bool:
        """Whether the sign was proven exactly."""
        return self.sign.is_definite

    def describe(self) -> str:
        """Return a one-line human summary including the deciding tactic."""
        head = f"{self.sign.symbol} ({self.sign.value})"
        if self.value is not None:
            head = f"{head}  value={self.value}"
        parts = [head, f"via {self.tactic}"]
        if self.assumptions:
            parts.append(f"assuming {', '.join(self.assumptions)}")
        if self.sign is Sign.UNKNOWN and self.free_names:
            parts.append(f"free: {', '.join(self.free_names)}")
        return "  ".join(parts)


def _normalize_assumptions(
    assume_positive: Iterable[str] | None,
    assume_nonzero: Iterable[str] | None = None,
) -> _Assumptions:
    """Build an :class:`_Assumptions` record from the caller's declarations."""
    positive = frozenset(assume_positive or ())
    # A positive parameter is non-zero by definition; fold that in so callers
    # need not state both.
    nonzero = frozenset(assume_nonzero or ()) | positive
    return _Assumptions(positive=positive, nonzero=nonzero)


def sign_of(
    expr: str | float | None,
    *,
    assume_positive: Iterable[str] | None = None,
    assume_nonzero: Iterable[str] | None = None,
    parameters: Mapping[str, float] | None = None,
) -> SignResult:
    """Decide the sign of a single coefficient expression.

    Parameters
    ----------
    expr : str | float | None
        Coefficient in Wolfram InputForm, a number, or ``None`` (= 1).
    assume_positive : Iterable[str] | None
        Parameter names the caller declares strictly positive.  Only these are
        assumed; nothing is positive by default.  Any assumption actually used
        is recorded on the result.
    assume_nonzero : Iterable[str] | None
        Parameter names the caller declares merely non-vanishing, without
        claiming a sign.  This is the usual physical situation for a coupling
        such as ``kappa``, whose vanishing would remove the Einstein-Hilbert
        term entirely, and it sharpens ``kappa^2`` from ``0+`` to ``+``.  Leave
        it empty for parameters that genuinely reach zero (``xi``, ``b5``).
    parameters : Mapping[str, float] | None
        Optional numeric values.  These populate :attr:`SignResult.numeric`
        for corroboration but never upgrade an unproven structural verdict.

    Returns
    -------
    SignResult
        The verdict, the deciding tactic, and the supporting context.
    """
    assumptions = _normalize_assumptions(assume_positive, assume_nonzero)
    node = _parse(expr)
    names = free_names(expr)

    numeric = evaluate_numeric(expr, parameters)

    # Rung 1: the expression reduces to a rational constant.
    constant = _as_constant(node)
    if constant is not None:
        return SignResult(
            sign=_sign_from_number(constant),
            tactic="literal" if not names else "normal-form",
            free_names=names,
            value=constant,
            numeric=numeric,
        )

    # Rung 2: the sign lattice, optionally using declared assumptions.
    sign = _sign_node(node, assumptions)
    used = assumptions.names_used(names)
    if sign is not Sign.UNKNOWN:
        return SignResult(
            sign=sign,
            tactic="assumption" if used else "lattice",
            free_names=names,
            assumptions=used,
            numeric=numeric,
        )

    return SignResult(
        sign=Sign.UNKNOWN,
        tactic="undecided",
        free_names=names,
        assumptions=used,
        numeric=numeric,
    )


def _as_constant(node: ast.expr) -> Fraction | None:
    """Return the exact rational value of *node* if it reduces to a constant.

    The numerator and denominator are compared as polynomials, so a symbolic
    expression that cancels — ``(-xi)/(-xi)``, the shape a kinetic-normalised
    self-term routinely takes — reduces to a constant even though neither side
    is one on its own.
    """
    try:
        ratio = _to_ratio(node)
    except (_NotPolynomialError, KineticEvalError):
        return None
    if not ratio.den:
        return None
    if not ratio.num:
        return Fraction(0)
    return _poly_quotient_constant(ratio.num, ratio.den)


def evaluate_numeric(
    expr: str | float | None,
    parameters: Mapping[str, float] | None,
) -> float | None:
    """Evaluate *expr* numerically, or return ``None`` if anything is unresolved.

    Delegates to :func:`tidal.symbolic._kinetic_eval.evaluate_with_substitutions`,
    the repository's existing restricted-AST evaluator, rather than reimplementing
    evaluation here.  ``E`` is bound to Euler's number so that ``E**u`` terms
    from localised-background coefficients evaluate.

    This result is only ever *corroboration*: it is reported alongside a
    structural verdict and never used to justify one, because a value at one
    point in parameter space proves nothing about the sign in general.

    Parameters
    ----------
    expr : str | float | None
        Coefficient expression; ``None`` means ``1``.
    parameters : Mapping[str, float] | None
        Symbol values. ``None`` disables evaluation.

    Returns
    -------
    float | None
        The value, or ``None`` when unresolved, non-finite, or malformed.
    """
    if parameters is None:
        return None
    if expr is None:
        return 1.0
    if isinstance(expr, (int, float)):
        return float(expr)
    substitutions = {"E": math.e, **{k: float(v) for k, v in parameters.items()}}
    try:
        value = evaluate_with_substitutions(expr, substitutions)
    except (KineticEvalError, ZeroDivisionError, OverflowError, ValueError):
        return None
    if value is None or not math.isfinite(value):
        return None
    return value


def constant_ratio(
    numerator: str | float | None,
    denominator: str | float | None,
) -> Fraction | None:
    """Return ``numerator / denominator`` when it is exactly a rational constant.

    Decided by cross-multiplication on the rational normal form, so no division
    is performed and no parameter values are needed.  This is what makes an
    overall rescaling of an equation provably invisible: if both sides were
    multiplied by the same factor, the ratio is exactly ``1``.

    Parameters
    ----------
    numerator, denominator : str | float | None
        Coefficient expressions; ``None`` means ``1``.

    Returns
    -------
    Fraction | None
        The exact ratio, or ``None`` when it is not a constant (or cannot be
        expanded rationally).
    """
    try:
        left = _to_ratio(_parse(numerator))
        right = _to_ratio(_parse(denominator))
    except (_NotPolynomialError, KineticEvalError):
        return None
    if right.is_zero():
        return None
    if left.is_zero():
        return Fraction(0)
    # left/right == c  <=>  left.num * right.den == c * (right.num * left.den)
    lhs = _poly_mul(left.num, right.den)
    rhs = _poly_mul(right.num, left.den)
    return _poly_quotient_constant(lhs, rhs)


def _poly_quotient_constant(lhs: _Poly, rhs: _Poly) -> Fraction | None:
    """Return ``c`` if ``lhs == c * rhs`` identically, else ``None``."""
    if not rhs or len(lhs) != len(rhs):
        return None
    ratio: Fraction | None = None
    for mono, coeff in rhs.items():
        if mono not in lhs:
            return None
        candidate = lhs[mono] / coeff
        if ratio is None:
            ratio = candidate
        elif candidate != ratio:
            return None
    return ratio


def _poly_quotient_monomial(
    lhs: _Poly,
    rhs: _Poly,
) -> tuple[Fraction, _Monomial] | None:
    """Return ``(c, m)`` if ``lhs == c * m * rhs`` for a monomial *m*, else ``None``.

    This generalises :func:`_poly_quotient_constant` to the case where the two
    polynomials differ by a common monomial factor as well as a rational one —
    ``-4*beta1^2*P`` over ``2*P`` is ``-2 * beta1^2``.  A candidate ``(c, m)``
    is *guessed* from one pair of terms and then **verified** by reconstructing
    ``c * m * rhs`` and comparing exactly, so the result is sound however the
    guess was formed.
    """
    if not rhs or not lhs or len(lhs) != len(rhs):
        return None
    lhs_key = min(lhs, key=repr)
    rhs_key = min(rhs, key=repr)
    coeff = lhs[lhs_key] / rhs[rhs_key]
    # Candidate monomial m = lhs_key / rhs_key.
    mono = _mono_mul(lhs_key, tuple((name, -exp) for name, exp in rhs_key))
    # Verify exactly.
    if _poly_mul({mono: coeff}, rhs) != lhs:
        return None
    return coeff, mono


def _monomial_sign(mono: _Monomial) -> Sign:
    """Return the proven sign of a bare monomial.

    An all-even monomial is non-negative.  It is strictly positive only when
    *every* atom carries a negative exponent, since then each base must be
    non-zero for the expression to be defined; a single positive even exponent
    (``a^2``) still admits ``a = 0``.
    """
    if not mono:
        return Sign.POSITIVE
    if any(exp % 2 for _, exp in mono):
        return Sign.UNKNOWN
    if all(exp < 0 for _, exp in mono):
        return Sign.POSITIVE
    return Sign.NONNEGATIVE


def _monomial_quotient(
    numerator: str | float | None,
    denominator: str | float | None,
) -> tuple[Fraction, _Monomial] | None:
    """Return ``(c, m)`` with ``numerator/denominator == c * m``, else ``None``."""
    try:
        left = _to_ratio(_parse(numerator))
        right = _to_ratio(_parse(denominator))
    except (_NotPolynomialError, KineticEvalError):
        return None
    if right.is_zero() or left.is_zero():
        return None
    return _poly_quotient_monomial(
        _poly_mul(left.num, right.den),
        _poly_mul(right.num, left.den),
    )


def ratio_sign(
    numerator: str | float | None,
    denominator: str | float | None,
    *,
    assume_positive: Iterable[str] | None = None,
    assume_nonzero: Iterable[str] | None = None,
    parameters: Mapping[str, float] | None = None,
) -> SignResult:
    """Decide the sign of ``numerator / denominator``.

    This is the primitive behind every relative question asked of a spec:
    comparing sibling components, or comparing an equation before and after
    re-derivation.  Ratios *sometimes* cancel unknown parameters — ``-xi/-xi``
    is exactly ``1`` regardless of ``xi`` — but not always, and sums over
    distinct parameters generally do not cancel.  Those cases return
    ``Sign.UNKNOWN`` rather than a guess.

    Parameters
    ----------
    numerator, denominator : str | float | None
        Coefficient expressions; ``None`` means ``1``.
    assume_positive : Iterable[str] | None
        Parameter names the caller declares strictly positive.  Only these are
        assumed; nothing is positive by default.  Any assumption actually used
        is recorded on the result.
    assume_nonzero : Iterable[str] | None
        Parameter names the caller declares merely non-vanishing, without
        claiming a sign.  This is the usual physical situation for a coupling
        such as ``kappa``, whose vanishing would remove the Einstein-Hilbert
        term entirely, and it sharpens ``kappa^2`` from ``0+`` to ``+``.  Leave
        it empty for parameters that genuinely reach zero (``xi``, ``b5``).
    parameters : Mapping[str, float] | None
        Optional numeric values, used only for corroboration.

    Returns
    -------
    SignResult
        The verdict, with :attr:`SignResult.value` set when the ratio reduced
        to an exact rational constant.
    """
    assumptions = _normalize_assumptions(assume_positive, assume_nonzero)
    names = tuple(sorted(set(free_names(numerator)) | set(free_names(denominator))))

    # Rung 1-3: exact constant ratio via rational normal form.
    exact = constant_ratio(numerator, denominator)
    if exact is not None:
        return SignResult(
            sign=_sign_from_number(exact),
            tactic="normal-form",
            free_names=names,
            value=exact,
        )

    # Rung 4: the ratio is a rational constant times a monomial, whose sign
    # may still be provable even though the monomial's value is not.
    monomial = _monomial_quotient(numerator, denominator)
    if monomial is not None:
        coeff, mono = monomial
        sign = _sign_mul(_sign_from_number(coeff), _monomial_sign(mono))
        if sign is not Sign.UNKNOWN:
            return SignResult(
                sign=sign,
                tactic="monomial-quotient",
                free_names=names,
            )

    # Rung 5: independent lattice verdicts on each side.
    num_sign = _sign_node(_parse(numerator), assumptions)
    den_sign = _refine_nonzero(_sign_node(_parse(denominator), assumptions))
    sign = _sign_mul(num_sign, den_sign)
    used = assumptions.names_used(names)

    numeric: float | None = None
    if parameters:
        num_val = evaluate_numeric(numerator, parameters)
        den_val = evaluate_numeric(denominator, parameters)
        if num_val is not None and den_val not in {None, 0.0}:
            numeric = num_val / den_val  # type: ignore[operator]

    if sign is not Sign.UNKNOWN:
        return SignResult(
            sign=sign,
            tactic="assumption" if used else "lattice",
            free_names=names,
            assumptions=used,
            numeric=numeric,
        )

    return SignResult(
        sign=Sign.UNKNOWN,
        tactic="undecided",
        free_names=names,
        assumptions=used,
        numeric=numeric,
    )


def are_equal(
    left: str | float | None,
    right: str | float | None,
) -> bool | None:
    """Decide whether two coefficient expressions are identically equal.

    Three-valued: ``True`` when proven equal, ``False`` when proven different,
    and ``None`` when neither could be established.

    Parameters
    ----------
    left, right : str | float | None
        Coefficient expressions; ``None`` means ``1``.

    Returns
    -------
    bool | None
        Proven equality, proven inequality, or undecided.
    """
    ratio = constant_ratio(left, right)
    if ratio is not None:
        return ratio == 1
    # Fall back to canonical strings, which also cover the non-rational cases.
    if canonical_form(left) == canonical_form(right):
        return True
    return None
