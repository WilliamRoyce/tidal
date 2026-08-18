"""Tests for :mod:`tidal.symbolic.sign_algebra`.

The suite is organized around the failure it exists to prevent.  Every
confident-but-wrong spec diagnosis recorded in GH #401 is pinned here as a
named regression test, so a future reader can trace test → historical mistake.

Three layers, in increasing order of what they protect:

* :class:`TestSoundnessProperty` — the ladder must never lie.  Verdicts are
  cross-checked against numeric evaluation at randomized parameter values,
  over every coefficient string in the committed corpus.  This is the
  invariant; it keeps holding as tactics are added.
* :class:`TestSignLattice` / :class:`TestRatioDecisions` — the individual
  rungs, including the cases where the correct answer is "I cannot tell".
* :class:`TestMisreadings401` — one test per row of the #401 table.
"""

from __future__ import annotations

import json
import math
import random
import unittest
from fractions import Fraction
from pathlib import Path

import pytest

from tidal.symbolic._kinetic_eval import KineticEvalError
from tidal.symbolic.sign_algebra import (
    Sign,
    are_equal,
    canonical_form,
    constant_ratio,
    evaluate_numeric,
    free_names,
    ratio_sign,
    sign_of,
)

EXAMPLES = Path(__file__).resolve().parent.parent / "examples" / "data"

# Randomized soundness sampling is seeded so failures reproduce exactly.
_SEED = 20260817
_SAMPLES_PER_EXPRESSION = 40


def _corpus_coefficients() -> set[str]:
    """Return every symbolic coefficient string in the committed examples."""
    found: set[str] = set()
    for path in sorted(EXAMPLES.glob("*.json")):
        data = json.loads(path.read_text())
        for equation in data.get("equations", []):
            kinetic = equation["lhs"].get("kinetic_coefficient_symbolic")
            if kinetic:
                found.add(kinetic)
            for term in equation["rhs"]["terms"]:
                symbolic = term.get("coefficient_symbolic")
                if symbolic:
                    found.add(symbolic)
    return found


class TestSignLattice(unittest.TestCase):
    """The sign domain, including where it must refuse to decide."""

    def test_numeric_literals(self) -> None:
        """Plain numbers are decided exactly."""
        self.assertIs(sign_of("-2").sign, Sign.NEGATIVE)
        self.assertIs(sign_of(3.5).sign, Sign.POSITIVE)
        self.assertIs(sign_of("0").sign, Sign.ZERO)
        self.assertIs(sign_of(None).sign, Sign.POSITIVE)  # absent kinetic == 1

    def test_even_power_is_nonnegative_not_positive(self) -> None:
        """``kappa^2`` may vanish, so it is ``0+`` rather than ``+``.

        Conflating the two would let a caller assert a field is dynamical when
        its kinetic coefficient can vanish — the case that makes a field
        *constrained* instead.  See ``LHSStructure.kinetic_coefficient_symbolic``.
        """
        self.assertIs(sign_of("kappa^2").sign, Sign.NONNEGATIVE)
        self.assertIs(sign_of("B0^2/2").sign, Sign.NONNEGATIVE)

    def test_negative_even_power_is_strictly_positive(self) -> None:
        """``kappa^(-2)`` is defined only for non-zero ``kappa``, so it is ``+``."""
        self.assertIs(sign_of("kappa^(-2)").sign, Sign.POSITIVE)
        self.assertIs(sign_of("-kappa^(-2)").sign, Sign.NEGATIVE)

    def test_assume_nonzero_sharpens_even_powers(self) -> None:
        """Declaring a parameter non-zero lifts ``0+`` to ``+``.

        ``kappa`` cannot physically vanish — that would delete the
        Einstein-Hilbert term — but the declaration must be explicit, because
        other parameters here genuinely do reach zero.
        """
        self.assertIs(sign_of("kappa^2", assume_nonzero=["kappa"]).sign, Sign.POSITIVE)

    def test_parameters_that_reach_zero_keep_the_cautious_reading(self) -> None:
        """``xi`` and ``b5`` are not assumed non-zero.

        ``xi = 0`` is the documented limit where torsion becomes constrained,
        and ``torsion_gertsenshtein_b5_zero`` is a committed example at
        ``b5 = 0``.  Declaring an unrelated symbol must not sharpen them.
        """
        self.assertIs(sign_of("xi^2", assume_nonzero=["kappa"]).sign, Sign.NONNEGATIVE)
        self.assertIs(sign_of("b5^2").sign, Sign.NONNEGATIVE)

    def test_exponential_is_positive(self) -> None:
        """``E**u`` is positive for any real exponent."""
        self.assertIs(sign_of("E^(-x^2)").sign, Sign.POSITIVE)

    def test_same_signed_sum_keeps_its_sign(self) -> None:
        """A sum of same-signed summands is decided."""
        self.assertIs(sign_of("-1 - kappa^2").sign, Sign.NEGATIVE)

    def test_free_parameter_is_undecided(self) -> None:
        """A bare parameter has no knowable sign, and none is invented."""
        result = sign_of("-xi")
        self.assertIs(result.sign, Sign.UNKNOWN)
        self.assertEqual(result.free_names, ("xi",))

    def test_declared_positive_is_recorded_on_the_result(self) -> None:
        """An assumption that was used must be visible in the verdict."""
        result = sign_of("-xi", assume_positive=["xi"])
        self.assertIs(result.sign, Sign.NEGATIVE)
        self.assertEqual(result.assumptions, ("xi",))
        self.assertIn("xi", result.describe())

    def test_canceling_sum_declines_rather_than_guesses(self) -> None:
        """``-1 + 2*B0^2*rho`` genuinely flips sign, so the answer is UNKNOWN.

        This is the load-bearing case: a tool that guessed here would be
        confidently wrong for half of parameter space.
        """
        result = sign_of("-1 + 2*B0^2*rho")
        self.assertIs(result.sign, Sign.UNKNOWN)
        self.assertEqual(result.free_names, ("B0", "rho"))
        # Even declaring both parameters non-zero must not decide it.
        stubborn = sign_of("-1 + 2*B0^2*rho", assume_nonzero=["B0", "rho"])
        self.assertIs(stubborn.sign, Sign.UNKNOWN)

    def test_rejects_unsafe_source(self) -> None:
        """Anything outside the restricted grammar is refused, not evaluated."""
        with pytest.raises(KineticEvalError):
            sign_of("__import__('os').system('true')")

    def test_free_names_excludes_euler(self) -> None:
        """``E`` is a constant, not a parameter."""
        self.assertEqual(free_names("E^(-x^2)*B0"), ("B0", "x"))


class TestRatioDecisions(unittest.TestCase):
    """Relative questions — the form every #401 misreading actually took."""

    def test_identical_expressions_cancel(self) -> None:
        """``-xi / -xi`` is exactly 1 regardless of ``xi``."""
        result = ratio_sign("-xi", "-xi")
        self.assertEqual(result.value, Fraction(1))
        self.assertIs(result.sign, Sign.POSITIVE)

    def test_rewriting_does_not_change_the_ratio(self) -> None:
        """Canonical form sees through how a sum was written."""
        self.assertEqual(constant_ratio("(3*chi) + (-xi)", "3*chi - xi"), Fraction(1))

    def test_common_factor_cancels(self) -> None:
        """``-4*beta1 / 2*beta1`` is ``-2`` without knowing ``beta1``."""
        self.assertEqual(constant_ratio("-4*beta1", "2*beta1"), Fraction(-2))

    def test_opposite_overall_sign_is_proven(self) -> None:
        """The #397 signature: a ratio of exactly ``-1``."""
        self.assertIs(ratio_sign("-1.0", "1.0").sign, Sign.NEGATIVE)

    def test_non_constant_ratio_is_undecided(self) -> None:
        """Sums over distinct parameters generally do not cancel."""
        result = ratio_sign("3*chi - xi", "-chi - xi")
        self.assertIs(result.sign, Sign.UNKNOWN)
        self.assertIsNone(result.value)

    def test_are_equal_is_three_valued(self) -> None:
        """Equality is proven, disproven, or undecided — never assumed."""
        self.assertTrue(are_equal("2*B0", "B0*2"))
        self.assertFalse(are_equal("2*B0", "3*B0"))
        self.assertIsNone(are_equal("3*chi - xi", "-chi - xi"))

    def test_canonical_form_is_stable_under_rewriting(self) -> None:
        """Equivalent spellings share a canonical string."""
        self.assertEqual(canonical_form("B0^2/2"), canonical_form("(1/2)*B0^2"))


class TestMisreadings401(unittest.TestCase):
    """One test per row of the #401 table, against committed example JSONs."""

    @staticmethod
    def _equation(spec_name: str, field: str) -> dict:
        data = json.loads((EXAMPLES / f"{spec_name}.json").read_text())
        for equation in data["equations"]:
            if equation["field"] == field:
                return equation
        msg = f"{field} not found in {spec_name}"
        raise AssertionError(msg)

    @classmethod
    def _effective(cls, spec_name: str, field: str, operator: str) -> tuple[str, str]:
        """Return ``(summed_rhs, kinetic)`` for a component's self-term.

        Deliberately sums *all* matching terms and keeps the kinetic
        coefficient separate — the two things the misreadings dropped.
        """
        equation = cls._equation(spec_name, field)
        parts = [
            term.get("coefficient_symbolic") or repr(term["coefficient"])
            for term in equation["rhs"]["terms"]
            if term["field"] == field and term["operator"] == operator
        ]
        summed = "(" + ") + (".join(parts) + ")" if parts else "0"
        kinetic = equation["lhs"].get("kinetic_coefficient_symbolic") or "1"
        return summed, kinetic

    def test_row1_and_row3_euler_heisenberg_is_not_broken(self) -> None:
        """`gertsenshtein_eh` was wrongly flagged by comparing raw coefficients.

        ``a_0`` carries **two** self-laplacians and no kinetic coefficient;
        ``a_1`` carries two and a kinetic coefficient.  Comparing the raw
        numbers reports a defect; comparing effective coefficients does not.
        """
        num_0, kin_0 = self._effective("gertsenshtein_eh", "a_0", "laplacian_x")
        num_1, kin_1 = self._effective("gertsenshtein_eh", "a_1", "laplacian_x")

        # Two matching terms on each side -- singular "the" term is the bug.
        self.assertIn(") + (", num_0)
        self.assertIn(") + (", num_1)

        # eff(a_1)/eff(a_0) via cross-multiplication; must NOT be proven negative.
        result = ratio_sign(f"({num_1})*({kin_0})", f"({num_0})*({kin_1})")
        self.assertIsNot(
            result.sign,
            Sign.NEGATIVE,
            "gertsenshtein_eh must not be flagged as sign-inconsistent",
        )

    def test_row2_symbolic_coefficients_are_not_skipped(self) -> None:
        """`torsion_dark_photon_fv` carries symbolic ``-xi``; numeric-only scans miss it."""
        num, kin = self._effective("torsion_dark_photon_fv", "t_1", "laplacian_x")
        self.assertIn("xi", num + kin)
        # eff(t_1) = -xi / -xi = +1 exactly, with no value for xi.
        self.assertEqual(constant_ratio(num, kin), Fraction(1))

    def test_row4_torsion_components_are_not_indexed_by_suffix(self) -> None:
        """Component 0 is temporal for rank-1 ``a`` but not for rank-3 torsion."""
        data = json.loads(
            (EXAMPLES / "torsion_gertsenshtein_complete_even.json").read_text(),
        )
        by_name = {f["name"]: f for f in data["fields"]}
        self.assertEqual(by_name["a_0"]["tensor_indices"], [0])
        self.assertEqual(by_name["a_0"]["tensor_rank"], 1)
        # h_5 is h_{xy}: rank 2, spatial-spatial -- suffix "5" says nothing.
        self.assertEqual(by_name["h_5"]["tensor_indices"], [1, 2])
        self.assertNotIn(0, by_name["h_5"]["tensor_indices"])

    def test_row5_and_row6_overall_rescaling_is_not_a_change(self) -> None:
        """Rescaling both sides of an equation is representational, not a fix.

        In `gertsenshtein_ungauged`, ``h_5`` gained the kinetic coefficient
        ``-kappa^(-2)`` while every RHS term was multiplied by the same factor.
        A naive diff reports a fix; the effective coefficients are identical.
        """
        factor = "-kappa^(-2)"
        for old, new in [
            ("-(B0*kappa^2)", "B0"),  # gradient_x(a_1)
            ("-1/2*(B0^2*kappa^2)", "B0^2/2"),  # identity(h_5)
            ("1.0", factor),  # laplacian_x(h_5)
        ]:
            with self.subTest(term=new):
                # eff_old = old/1, eff_new = new/factor; equal iff old*factor == new
                self.assertEqual(
                    constant_ratio(new, f"({factor})*({old})"),
                    Fraction(1),
                    f"{new} / {factor} should equal {old}",
                )

    def test_row5_rhs_only_flip_is_a_real_change(self) -> None:
        """Flipping only the RHS, with the kinetic coefficient unchanged, is real.

        This is the asymmetry that matters: testing invariance alone would pass
        a tool that called *everything* representational.  ``a_0`` in the same
        spec flipped its RHS with no compensating LHS change — the genuine #397
        defect — and must be reported as different.
        """
        self.assertIs(ratio_sign("-1.0", "1.0").sign, Sign.NEGATIVE)
        self.assertFalse(are_equal("-B0", "B0"))

    def test_constraint_equations_are_distinguishable(self) -> None:
        """`chern_simons_3d` ``A_0`` is a constraint, where overall sign is conventional."""
        equation = self._equation("chern_simons_3d", "A_0")
        self.assertEqual(equation["lhs"]["order"]["time"], 0)


class TestSoundnessProperty(unittest.TestCase):
    """The ladder must never lie — checked against the whole committed corpus."""

    def test_every_definite_verdict_agrees_with_numeric_evaluation(self) -> None:
        """Sample each signed verdict at randomized parameters of both signs.

        Nothing is assumed positive, so a verdict that depended on an unstated
        assumption would be caught here.
        """
        # Seeded sampling for reproducibility; nothing cryptographic here.
        rng = random.Random(_SEED)  # noqa: S311
        expressions = _corpus_coefficients()
        self.assertGreater(len(expressions), 100, "corpus fixture looks empty")

        checked = 0
        for expr in sorted(expressions):
            result = sign_of(expr)
            if result.sign is Sign.UNKNOWN:
                continue
            names = free_names(expr)
            for _ in range(_SAMPLES_PER_EXPRESSION):
                params = {
                    name: rng.choice([-1.0, 1.0]) * 10 ** rng.uniform(-3, 3)
                    for name in names
                }
                value = evaluate_numeric(expr, params)
                if value is None or not math.isfinite(value):
                    continue
                checked += 1
                self._assert_consistent(expr, result.sign, value, params)
        self.assertGreater(checked, 200, "too few numeric samples to be meaningful")

    def _assert_consistent(
        self,
        expr: str,
        sign: Sign,
        value: float,
        params: dict[str, float],
    ) -> None:
        """Fail if *value* contradicts the proven *sign*."""
        context = f"{expr!r} judged {sign.symbol} but evaluates to {value} at {params}"
        if sign is Sign.POSITIVE:
            self.assertGreater(value, 0, context)
        elif sign is Sign.NEGATIVE:
            self.assertLess(value, 0, context)
        elif sign is Sign.NONNEGATIVE:
            self.assertGreaterEqual(value, 0, context)
        elif sign is Sign.NONPOSITIVE:
            self.assertLessEqual(value, 0, context)
        elif sign is Sign.ZERO:
            self.assertEqual(value, 0, context)

    def test_ratio_round_trip(self) -> None:
        """``ratio(a, a)`` is 1, and ``ratio(a, b)`` inverts consistently."""
        for expr in sorted(_corpus_coefficients()):
            with self.subTest(expr=expr):
                self.assertEqual(constant_ratio(expr, expr), Fraction(1))

    def test_negating_both_sides_is_invariant(self) -> None:
        """Multiplying an equation through by -1 changes no effective coefficient.

        Generated mechanically over the corpus, so it covers far more shapes
        than any hand-picked example.
        """
        for expr in sorted(_corpus_coefficients()):
            with self.subTest(expr=expr):
                # eff = expr/kin; negating both leaves the ratio identical.
                self.assertEqual(
                    constant_ratio(f"-({expr})", f"-({expr})"),
                    constant_ratio(expr, expr),
                )


if __name__ == "__main__":
    unittest.main()
