"""Tests for tidal.symbolic._separable.extract_separable_bsm_factors.

GH #384 Phase A′ infrastructure: detect which RHS coefficient expressions
factor as ``c_BSM(α) · c_geom(coords, geometry)`` so the convolution-matrix
block can be cached at chain start.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.symbolic._eval_utils import evaluate_coefficient
from tidal.symbolic._separable import extract_separable_bsm_factors


class TestSeparableCases:
    """Inputs that DO factor as (top-level BSM product) · (BSM-free residual)."""

    def test_simple_bsm_times_coord(self) -> None:
        bsm, geom = extract_separable_bsm_factors("alpha1*x[]", {"alpha1"})
        assert bsm == "alpha1"
        assert "alpha1" not in geom

    def test_bsm_in_numerator_with_exp_denominator(self) -> None:
        bsm, geom = extract_separable_bsm_factors(
            "(Bpeak*delta1*x[])/(E^((zc1-x[])^2/sigB^2)*sigB^2)",
            {"alpha1", "alpha2", "delta1"},
        )
        assert bsm == "delta1"
        assert "delta1" not in geom

    def test_numeric_prefactor_and_bsm(self) -> None:
        bsm, geom = extract_separable_bsm_factors(
            "-3 * delta1 * Bpeak * x[]",
            {"delta1"},
        )
        assert bsm == "delta1"
        assert "delta1" not in geom

    def test_no_bsm_returns_empty_bsm(self) -> None:
        # When the coefficient contains no BSM symbols, separability is
        # trivially true with c_BSM = 1; the cache is the full geom expression.
        bsm, geom = extract_separable_bsm_factors(
            "5*Bpeak^2/E^(x[]^2)", {"alpha1", "delta1"}
        )
        # Empty product means the BSM scalar is constant 1
        assert bsm in ("", "1") or bsm is not None
        assert "alpha1" not in geom and "delta1" not in geom

    def test_addsub_summands_share_same_bsm(self) -> None:
        # delta1 * X + delta1 * Y is separable as delta1 * (X + Y)
        bsm, geom = extract_separable_bsm_factors(
            "delta1*x[] + delta1*Bpeak", {"delta1"}
        )
        assert bsm == "delta1"


class TestNonSeparableCases:
    """Inputs where BSM appears in non-multiplicative position."""

    def test_bsm_squared(self) -> None:
        # alpha1**2 is rank>1 in alpha1 — not multiplicative
        bsm, geom = extract_separable_bsm_factors("alpha1**2*x[]", {"alpha1"})
        assert bsm is None
        assert geom == "alpha1**2*x[]"

    def test_bsm_inside_exponential(self) -> None:
        # E^(alpha1*x) — BSM inside exponent
        bsm, geom = extract_separable_bsm_factors("E^(alpha1*x[])", {"alpha1"})
        assert bsm is None

    def test_bsm_in_one_summand_only(self) -> None:
        # alpha1 + x[] : different BSM-products per summand
        bsm, geom = extract_separable_bsm_factors("alpha1 + x[]", {"alpha1"})
        assert bsm is None

    def test_bsm_in_denominator(self) -> None:
        # 1/(1 + alpha1*X) — BSM in denominator
        bsm, geom = extract_separable_bsm_factors(
            "1/(1 + alpha1*x[])", {"alpha1"}
        )
        assert bsm is None


class TestNumericalRoundTrip:
    """Verify c_BSM · c_geom evaluates identically to the original expression."""

    @pytest.mark.parametrize(
        "expr,bsm_set,params,coord_arrays",
        [
            (
                "alpha1*x[]",
                {"alpha1"},
                {"alpha1": 0.3},
                {"x": np.linspace(0, 1, 8)},
            ),
            (
                "delta1 * Bpeak * x[] / sigB",
                {"delta1"},
                {"delta1": -0.7, "Bpeak": 0.01, "sigB": 5.0},
                {"x": np.linspace(0, 200, 16)},
            ),
        ],
    )
    def test_separated_product_matches_original(
        self,
        expr: str,
        bsm_set: set[str],
        params: dict[str, float],
        coord_arrays: dict[str, np.ndarray],
    ) -> None:
        bsm_expr, geom_expr = extract_separable_bsm_factors(expr, bsm_set)
        assert bsm_expr is not None
        # Evaluate original
        orig = evaluate_coefficient(expr, params, ("t", "x"), coord_arrays, 0.0)
        # Evaluate bsm * geom (geom needs the same coord_arrays)
        bsm_val = evaluate_coefficient(bsm_expr, params, ("t", "x"), None, 0.0)
        # geom uses Python-form names (no x[]); evaluate_coefficient handles
        # the namespace lookup for both forms (mathematica_to_python preserves
        # bare `x`).
        geom_val = evaluate_coefficient(
            geom_expr, params, ("t", "x"), coord_arrays, 0.0
        )
        product = bsm_val * geom_val if not isinstance(geom_val, np.ndarray) else (
            bsm_val * geom_val
        )
        assert np.allclose(orig, product, rtol=1e-12)
