"""Tests for tidal.symbolic._kinetic_eval (Gap B of v6 plan).

Covers the small set of coefficient shapes that the Wolfram pipeline
actually emits:
  * Simple monomials in small parameters: ``"2*b5"``, ``"-b5"``
  * Rational expressions: ``"(-25*b5)/2"``, ``"B0^2/8"``
  * Negative exponents: ``"-kappa^(-2)"`` (caret → ``**``)
  * Mixed small + non-small: must return None (still symbolic, not
    literal zero)
  * Edge cases: literal zero, self-cancelling ``"b5 - b5"``, malformed
    inputs

Safety: no parser node outside the restricted subset (Constant, Name,
UnaryOp, BinOp over {Add, Sub, Mult, Div, Pow}) is tolerated —
confirmed by the malformed-input tests.
"""

from __future__ import annotations

import pytest

from tidal.symbolic._kinetic_eval import (
    KineticEvalError,
    evaluate_at_one,
    evaluate_at_zero,
    evaluate_with_substitutions,
    lhs_collapses_to_zero,
    split_small_parameter_kinetic,
)


class TestEvaluateAtZero:
    # --- Literal zero / self-cancelling ---

    def test_literal_zero(self) -> None:
        assert evaluate_at_zero("0", frozenset()) == 0.0
        assert evaluate_at_zero("0.0", frozenset()) == 0.0

    def test_self_cancelling_small_param(self) -> None:
        # b5 - b5 with b5 substituted → 0 - 0 = 0
        assert evaluate_at_zero("b5 - b5", {"b5"}) == 0.0

    # --- Single small parameter ---

    def test_linear_small(self) -> None:
        assert evaluate_at_zero("2*b5", {"b5"}) == 0.0
        assert evaluate_at_zero("b5", {"b5"}) == 0.0
        assert evaluate_at_zero("-b5", {"b5"}) == 0.0

    def test_rational_small_numerator(self) -> None:
        # Real shipped shape from torsion_gertsenshtein.json
        assert evaluate_at_zero("(-25*b5)/2", {"b5"}) == 0.0
        assert evaluate_at_zero("(25*b5)/2", {"b5"}) == 0.0

    # --- Multi-parameter polynomials ---

    def test_sum_of_small_params(self) -> None:
        assert evaluate_at_zero("2*b5 + 3*b1", {"b5", "b1"}) == 0.0

    def test_product_of_small_params(self) -> None:
        assert evaluate_at_zero("b5 * b1", {"b5", "b1"}) == 0.0

    # --- Still-symbolic cases (return None) ---

    def test_non_small_parameter_only(self) -> None:
        # Shipped shapes: kappa, m2, xi — none are small parameters
        assert evaluate_at_zero("-kappa^(-2)", {"b5"}) is None
        assert evaluate_at_zero("m2", {"b5"}) is None
        assert evaluate_at_zero("-xi", {"b5"}) is None
        assert evaluate_at_zero("B0^2/8", {"b5"}) is None

    def test_mixed_small_and_non_small(self) -> None:
        # b5 → 0 leaves m2 surviving: still symbolic, not demoted
        assert evaluate_at_zero("b5 + m2", {"b5"}) is None
        # Shipped shape: "-1/2*1/kappa^2"
        assert evaluate_at_zero("-1/2*1/kappa^2", {"b5"}) is None

    def test_b5_in_denominator_symbolic_numerator(self) -> None:
        # Unlikely in practice but worth verifying: numerator stays
        # symbolic, so returns None.
        assert evaluate_at_zero("m2 + b5", {"b5"}) is None

    # --- Caret normalization ---

    def test_caret_converted_to_power(self) -> None:
        # Bare exponent on a small parameter: b5^2 at b5=0 → 0
        assert evaluate_at_zero("b5^2", {"b5"}) == 0.0
        # Exponent on a non-small mixed with small: still symbolic (m2 survives)
        assert evaluate_at_zero("2*b5^2 + m2", {"b5"}) is None
        # Pure non-small exponent expression stays symbolic
        assert evaluate_at_zero("kappa^2 * 3", {"b5"}) is None

    def test_zero_to_negative_power_raises(self) -> None:
        # Contrived but mathematically well-defined: 0^(-1) is undefined.
        # The safeguard should refuse to return a nonsense value rather
        # than silently propagate NaN/inf.
        with pytest.raises(KineticEvalError, match="divide-by-zero"):
            evaluate_at_zero("b5^(-1)", {"b5"})

    def test_real_shipped_kappa_neg_two(self) -> None:
        # Exact shape from torsion_gertsenshtein.json
        assert evaluate_at_zero("-kappa^(-2)", {"b5"}) is None

    # --- Empty small-parameter set ---

    def test_empty_zero_names_returns_none_for_any_symbol(self) -> None:
        assert evaluate_at_zero("b5", frozenset()) is None
        assert evaluate_at_zero("2*b5 + 3", frozenset()) is None

    def test_empty_zero_names_evaluates_pure_literals(self) -> None:
        assert evaluate_at_zero("3 + 4", frozenset()) == 7.0
        assert evaluate_at_zero("2*5", frozenset()) == 10.0

    # --- Malformed / unsafe input ---

    def test_unbalanced_parens_raises(self) -> None:
        with pytest.raises(KineticEvalError, match="cannot parse"):
            evaluate_at_zero("(2*b5", {"b5"})

    def test_function_call_rejected(self) -> None:
        # sin(b5) uses ast.Call which is NOT in the safe subset
        with pytest.raises(KineticEvalError, match="unsupported AST node"):
            evaluate_at_zero("sin(b5)", {"b5"})

    def test_string_literal_rejected(self) -> None:
        with pytest.raises(KineticEvalError, match="unsupported literal"):
            evaluate_at_zero("'hello'", set())

    def test_non_string_input_rejected(self) -> None:
        with pytest.raises(KineticEvalError, match="expected a str"):
            evaluate_at_zero(42, {"b5"})  # type: ignore[arg-type]


class TestLhsCollapsesToZero:
    def test_none_kinetic_never_collapses(self) -> None:
        assert lhs_collapses_to_zero(None, ["b5"]) is False

    def test_empty_small_params_never_collapses(self) -> None:
        assert lhs_collapses_to_zero("2*b5", []) is False

    def test_collapses_with_matching_small(self) -> None:
        assert lhs_collapses_to_zero("2*b5", ["b5"]) is True
        assert lhs_collapses_to_zero("(-25*b5)/2", ["b5"]) is True

    def test_does_not_collapse_when_independent(self) -> None:
        assert lhs_collapses_to_zero("-kappa^(-2)", ["b5"]) is False
        assert lhs_collapses_to_zero("m2", ["b5"]) is False

    def test_multi_small_param_collapses_only_when_all_vanish(self) -> None:
        # "2*b5 + 3*b1": both are small → 0+0 = 0
        assert lhs_collapses_to_zero("2*b5 + 3*b1", ["b5", "b1"]) is True
        # Only b5 is small → 0 + 3*b1 = symbolic, does not collapse
        assert lhs_collapses_to_zero("2*b5 + 3*b1", ["b5"]) is False


class TestEvaluateAtOne:
    """v6 R8.2 / #290: extract numeric kinetic coefficients by
    substituting every small parameter with 1.0.

    Used by PerturbativeSolver.solve to recover ``K`` in the augmented
    constraint recovery formula ``h_c¹ += S_cc_inv · K · d^n_t(h_c⁰)``.
    """

    def test_simple_monomial_2_times_b5(self) -> None:
        assert evaluate_at_one("2*b5", {"b5"}) == 2.0

    def test_signed_rational_minus_25_times_b5_over_2(self) -> None:
        assert evaluate_at_one("(-25*b5)/2", {"b5"}) == -12.5

    def test_bare_b5_evaluates_to_one(self) -> None:
        assert evaluate_at_one("b5", {"b5"}) == 1.0

    def test_negated_b5(self) -> None:
        assert evaluate_at_one("-b5", {"b5"}) == -1.0

    def test_still_symbolic_returns_none(self) -> None:
        # kappa is not in the substitution set — still symbolic.
        assert evaluate_at_one("kappa*b5", {"b5"}) is None

    def test_pure_literal_unaffected(self) -> None:
        # No substitution needed: literal 3.0 evaluates as-is.
        assert evaluate_at_one("3.0", {"b5"}) == 3.0

    def test_b5_squared_is_order_epsilon_squared(self) -> None:
        # evaluate_at_one doesn't distinguish between O(ε) and O(ε²)
        # kinetics; the caller must gate on this separately. But we
        # document the raw numeric result here.
        assert evaluate_at_one("b5**2", {"b5"}) == 1.0


class TestEvaluateWithSubstitutions:
    """Generalized evaluator used by the two thin wrappers above."""

    def test_custom_value_map(self) -> None:
        # Evaluating "2*b5 + 3*b1" at b5=0.5, b1=2 gives 2*0.5 + 3*2 = 7.
        assert evaluate_with_substitutions("2*b5 + 3*b1", {"b5": 0.5, "b1": 2.0}) == 7.0

    def test_unsubstituted_name_returns_none(self) -> None:
        # Only b5 substituted; kappa still symbolic.
        assert evaluate_with_substitutions("b5 + kappa", {"b5": 1.0}) is None

    def test_caret_is_normalized_to_double_star(self) -> None:
        assert evaluate_with_substitutions("b5^3", {"b5": 2.0}) == 8.0


class TestSplitSmallParameterKinetic:
    """Split of a kinetic expression into (M₀, {param: c_p}) for #301 Phase 3.

    Must handle the forms ExportJSON.wl actually emits after ``Expand``:
    sums of monomials. Rejects parenthesised sub-sums, bilinear /
    quadratic small-parameter terms, and denominators containing small
    parameters — all gated by the #273 perturbative contract.
    """

    def _eq_val(self, expr: str, subs: dict[str, float]) -> float | None:
        """Evaluate an AST-unparsed expression to a float for sign-independent
        comparison (so `"-2*B0 ** 2"` and `"-(2*B0**2)"` both work).
        """
        return evaluate_with_substitutions(expr, subs)

    def test_euler_heisenberg_compound(self) -> None:
        """The exact form that Phase 1 emits for a_1's kinetic."""
        m0, coeffs = split_small_parameter_kinetic(
            "-1 + 2*B0^2*rho - 16*B0^2*sigma",
            ["rho", "sigma"],
        )
        assert self._eq_val(m0, {}) == -1.0
        assert self._eq_val(coeffs["rho"], {"B0": 1.0}) == 2.0
        assert self._eq_val(coeffs["sigma"], {"B0": 1.0}) == -16.0

    def test_single_small_parameter(self) -> None:
        m0, coeffs = split_small_parameter_kinetic(
            "-1 + 2*B0^2*rho",
            ["rho"],
        )
        assert self._eq_val(m0, {}) == -1.0
        assert self._eq_val(coeffs["rho"], {"B0": 1.0}) == 2.0

    def test_no_small_parameters_is_noop(self) -> None:
        """Empty small_parameters list → entire expression is M₀."""
        m0, coeffs = split_small_parameter_kinetic(
            "-1 + 2*B0^2*rho",
            [],
        )
        assert self._eq_val(m0, {"B0": 1.0, "rho": 0.5}) == 0.0  # -1 + 2*1*0.5
        assert coeffs == {}

    def test_small_param_not_present_coefficient_is_zero(self) -> None:
        """If a declared small parameter doesn't appear, its coefficient is ``"0"``."""
        m0, coeffs = split_small_parameter_kinetic(
            "-1 + 2*B0^2*rho",
            ["rho", "sigma"],
        )
        assert self._eq_val(m0, {}) == -1.0
        assert self._eq_val(coeffs["rho"], {"B0": 1.0}) == 2.0
        assert self._eq_val(coeffs["sigma"], {}) == 0.0

    def test_single_negative_monomial(self) -> None:
        m0, _coeffs = split_small_parameter_kinetic("-kappa^(-2)", [])
        assert self._eq_val(m0, {"kappa": 2.0}) == -0.25

    def test_pure_small_parameter(self) -> None:
        """Expression that's ONLY small parameters leaves M₀ = 0."""
        m0, coeffs = split_small_parameter_kinetic("rho", ["rho"])
        assert m0 == "0"
        assert self._eq_val(coeffs["rho"], {}) == 1.0

    def test_rejects_bilinear(self) -> None:
        """rho*sigma is order_in_eps=2 and gated by #273."""
        with pytest.raises(KineticEvalError, match="order_in_eps > 1"):
            split_small_parameter_kinetic(
                "1 + rho*sigma",
                ["rho", "sigma"],
            )

    def test_rejects_quadratic(self) -> None:
        """rho**2 is order_in_eps=2 via a single parameter."""
        with pytest.raises(KineticEvalError, match="#273"):
            split_small_parameter_kinetic("1 + rho**2", ["rho"])

    def test_rejects_small_param_in_denominator(self) -> None:
        """1/rho and 1/(1-rho) are out of the perturbative contract."""
        with pytest.raises(KineticEvalError, match="denominator"):
            split_small_parameter_kinetic("1 + 1/rho", ["rho"])

    def test_rejects_parenthesized_sum(self) -> None:
        """Mult over a sum like ``2*B0^2*(rho - 8*sigma)`` must be ``Expand``'d
        on the Wolfram side before it reaches this helper.
        """
        with pytest.raises(KineticEvalError, match="parenthesized sub-sum"):
            split_small_parameter_kinetic(
                "-1 + 2*B0^2*(rho - 8*sigma)",
                ["rho", "sigma"],
            )

    def test_rejects_non_str_input(self) -> None:
        with pytest.raises(KineticEvalError, match="expected a str"):
            split_small_parameter_kinetic(1.0, ["rho"])  # type: ignore[arg-type]

    def test_roundtrips_via_evaluate(self) -> None:
        """Sanity: M₀ + Σ c_p · p evaluated at some params equals the original."""
        params = {"B0": 1.5, "rho": 0.05, "sigma": 0.03}
        expr = "-1 + 2*B0^2*rho - 16*B0^2*sigma"
        m0, coeffs = split_small_parameter_kinetic(expr, ["rho", "sigma"])

        m0_val = evaluate_with_substitutions(m0, params)
        c_rho = evaluate_with_substitutions(coeffs["rho"], params)
        c_sigma = evaluate_with_substitutions(coeffs["sigma"], params)
        assert m0_val is not None
        assert c_rho is not None
        assert c_sigma is not None
        reconstructed = m0_val + c_rho * params["rho"] + c_sigma * params["sigma"]

        original = evaluate_with_substitutions(expr, params)
        assert original is not None
        assert abs(reconstructed - original) < 1e-12
