"""Tests for tidal/symbolic/_eval_utils.py — Mathematica→Python conversion and evaluation."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from tidal.symbolic._eval_utils import (
    _invert_exp_denominator,
    evaluate_coefficient,
    mathematica_to_python,
)


class TestInvertExpDenominator:
    """Tests for _invert_exp_denominator() covering both handled patterns."""

    # --- Pattern 1: /exp(arg) — simple inverted exp (pre-existing) ---

    def test_simple_slash_exp_constant(self) -> None:
        assert _invert_exp_denominator("1/exp(x)") == "1*exp(-(x))"

    def test_simple_slash_exp_preserves_prefix(self) -> None:
        assert _invert_exp_denominator("k/exp(x**2/R**2)") == "k*exp(-(x**2/R**2))"

    def test_simple_slash_exp_nested_parens(self) -> None:
        # Argument itself contains parens
        assert _invert_exp_denominator("1/exp(x**2/(2*R**2))") == "1*exp(-(x**2/(2*R**2)))"

    # --- Pattern 2: /(exp(arg)) — compound denominator, exp only ---

    def test_compound_exp_alone(self) -> None:
        assert _invert_exp_denominator("k/(exp(x**2))") == "k*exp(-(x**2))"

    def test_compound_exp_alone_nested(self) -> None:
        assert _invert_exp_denominator("-(k)/(exp(x**2/(2*R**2)))") == "-(k)*exp(-(x**2/(2*R**2)))"

    # --- Pattern 3: /(exp(arg)*rest) — compound denominator with trailing factor ---

    def test_compound_exp_with_rest(self) -> None:
        # exp at start of denominator, rest follows
        result = _invert_exp_denominator("-(Bpeak*x)/(exp(x**2/(2*R**2))*R**2)")
        assert result == "-(Bpeak*x)*exp(-(x**2/(2*R**2)))/(R**2)"

    def test_compound_exp_with_scalar_rest(self) -> None:
        result = _invert_exp_denominator("a/(exp(t)*2)")
        assert result == "a*exp(-(t))/(2)"

    # --- Pattern 4: /(N*exp(arg)) — numeric factor precedes exp ---

    def test_compound_exp_with_preceding_factor(self) -> None:
        # The actual failing pattern from expression 3 in gertsenshtein_localized.json
        result = _invert_exp_denominator("(Bpeak**2*kappa**2)/(2*exp(x**2/R**2))")
        assert result == "(Bpeak**2*kappa**2)*exp(-(x**2/R**2))/(2)"

    def test_compound_exp_with_preceding_and_trailing_factor(self) -> None:
        result = _invert_exp_denominator("k/(A*exp(x**2)*B)")
        assert result == "k*exp(-(x**2))/(A*B)"

    # --- Unchanged cases ---

    def test_no_exp_denominator_unchanged(self) -> None:
        expr = "Bpeak*exp(-(x**2/(2*R**2)))"
        assert _invert_exp_denominator(expr) == expr

    def test_empty_string(self) -> None:
        assert _invert_exp_denominator("") == ""

    def test_unrelated_expression_unchanged(self) -> None:
        expr = "sin(x) + cos(y)"
        assert _invert_exp_denominator(expr) == expr


class TestNoOverflowLocalized:
    """Regression guard: exact coefficient strings from gertsenshtein_localized.json
    must evaluate without RuntimeWarning overflow."""

    @pytest.mark.parametrize(
        "expr",
        [
            # B'(x) coupling coefficient: -(Bpeak*x)/(exp(x^2/(2R^2))*R^2)
            "-((Bpeak*x[])/(E^(x[]^2/(2*R^2))*R^2))",
            # B(x) coupling coefficient: Bpeak/exp(x^2/(2R^2))
            "Bpeak/E^(x[]^2/(2*R^2))",
            # h_7 mass coefficient: -1/2*(Bpeak^2*kappa^2)/exp(x^2/R^2)
            "-1/2*(Bpeak^2*kappa^2)/E^(x[]^2/R^2)",
            # Expression 3 — compound denominator with preceding factor /(2*exp(...))
            "(-((Bpeak^2*kappa^2)/(2*E^(x[]^2/R^2)))) + (-1/2*(Bpeak^2*kappa^2)/E^(x[]^2/R^2))",
        ],
    )
    def test_no_overflow_warning(self, expr: str) -> None:
        """Coefficients must not trigger RuntimeWarning at grid boundaries (|x|=100)."""
        x = np.linspace(-100.0, 100.0, 201)
        params = {"Bpeak": 0.1, "R": 2.0, "kappa": 1.0}
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result = evaluate_coefficient(expr, params, ("t", "x"), {"x": x})
        assert np.all(np.isfinite(result)), f"Non-finite values in result for expr: {expr}"

    def test_compound_denominator_value_correctness(self) -> None:
        """-(Bpeak*x)/(exp(x^2/(2R^2))*R^2) = -Bpeak*x*exp(-x^2/(2R^2))/R^2."""
        x = np.array([-10.0, -1.0, 0.0, 1.0, 10.0])
        Bpeak, R = 0.1, 2.0
        params = {"Bpeak": Bpeak, "R": R, "kappa": 1.0}
        expr = "-((Bpeak*x[])/(E^(x[]^2/(2*R^2))*R^2))"
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            result = evaluate_coefficient(expr, params, ("t", "x"), {"x": x})
        expected = -Bpeak * x * np.exp(-(x**2) / (2 * R**2)) / R**2
        np.testing.assert_allclose(result, expected, rtol=1e-12)


class TestMathematicaToPython:
    """Basic sanity checks for mathematica_to_python() (non-overflow paths)."""

    def test_exp_conversion(self) -> None:
        result = mathematica_to_python("Exp[-x[]^2]", ("t", "x"))
        # Wolfram Exp[-x^2] → exp(-x**2) (valid Python)
        assert "exp(" in result and "x" in result

    def test_e_power_conversion(self) -> None:
        result = mathematica_to_python("1/E^(x[]^2/R^2)", ("t", "x"))
        # Should become *exp(-(arg)) not /exp(+arg)
        assert "exp(-(" in result
        assert "/exp(" not in result
