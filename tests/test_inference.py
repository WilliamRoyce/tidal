"""Tests for ``tidal.inference`` — Bayesian inference infrastructure."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
    from pathlib import Path

    from tidal.inference._results import InferenceResult

# ===================================================================
# Prior tests
# ===================================================================


class TestPrior:
    """Test Prior dataclass and distributions."""

    def test_uniform_sample(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "uniform", 0.0, 1.0)
        rng = np.random.default_rng(42)
        samples = p.sample(rng, 1000)
        assert samples.shape == (1000,)
        assert np.all(samples >= 0.0)
        assert np.all(samples <= 1.0)

    def test_uniform_log_prob(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "uniform", 0.0, 10.0)
        assert p.log_prob(5.0) == pytest.approx(-math.log(10.0))
        assert p.log_prob(-1.0) == -math.inf
        assert p.log_prob(11.0) == -math.inf

    def test_uniform_transform(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "uniform", 2.0, 8.0)
        assert p.transform(0.0) == pytest.approx(2.0)
        assert p.transform(1.0) == pytest.approx(8.0)
        assert p.transform(0.5) == pytest.approx(5.0)

    def test_log_uniform_sample(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "log_uniform", 0.01, 10.0)
        rng = np.random.default_rng(42)
        samples = p.sample(rng, 1000)
        assert np.all(samples >= 0.01)
        assert np.all(samples <= 10.0)

    def test_log_uniform_transform(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "log_uniform", 0.01, 100.0)
        assert p.transform(0.0) == pytest.approx(0.01)
        assert p.transform(1.0) == pytest.approx(100.0)

    def test_log_uniform_requires_positive(self) -> None:
        from tidal.inference._prior import Prior

        with pytest.raises(ValueError, match="positive bounds"):
            Prior("x", "log_uniform", -1.0, 10.0)

    def test_normal_sample(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "normal", 0.0, 1.0)
        rng = np.random.default_rng(42)
        samples = p.sample(rng, 10000)
        assert abs(np.mean(samples)) < 0.1
        assert abs(np.std(samples) - 1.0) < 0.1

    def test_normal_log_prob(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "normal", 0.0, 1.0)
        # At the mean, log_prob = -0.5 * log(2*pi)
        expected = -0.5 * math.log(2 * math.pi)
        assert p.log_prob(0.0) == pytest.approx(expected)

    def test_normal_requires_positive_std(self) -> None:
        from tidal.inference._prior import Prior

        with pytest.raises(ValueError, match="positive"):
            Prior("x", "normal", 0.0, -1.0)

    def test_arctan_uniform_sample(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "arctan_uniform", -30.0, 30.0)
        rng = np.random.default_rng(42)
        samples = p.sample(rng, 1000)
        # Should have samples in a wide range including negative
        assert np.any(samples < -1.0)
        assert np.any(samples > 1.0)

    def test_arctan_uniform_transform(self) -> None:
        from tidal.inference._prior import Prior

        p = Prior("x", "arctan_uniform", -30.0, 30.0)
        # u=0.5 maps to theta=0, tan(0)=0
        assert p.transform(0.5) == pytest.approx(0.0)
        # u=0 and u=1 should be at the extremes
        assert p.transform(0.0) < -10
        assert p.transform(1.0) > 10

    def test_invalid_distribution(self) -> None:
        from tidal.inference._prior import Prior

        with pytest.raises(ValueError, match="Unknown distribution"):
            Prior("x", "cauchy", 0.0, 1.0)

    def test_bounds_order(self) -> None:
        from tidal.inference._prior import Prior

        with pytest.raises(ValueError, match="low < high"):
            Prior("x", "uniform", 10.0, 5.0)


class TestParsePrior:
    """Test CLI prior string parsing."""

    def test_basic(self) -> None:
        from tidal.inference._prior import parse_prior

        p = parse_prior("alpha=uniform:0.01:10")
        assert p.name == "alpha"
        assert p.distribution == "uniform"
        assert p.low == 0.01
        assert p.high == 10.0

    def test_log_uniform(self) -> None:
        from tidal.inference._prior import parse_prior

        p = parse_prior("xi=log_uniform:0.01:10")
        assert p.distribution == "log_uniform"

    def test_missing_equals(self) -> None:
        from tidal.inference._prior import parse_prior

        with pytest.raises(ValueError, match="must contain"):
            parse_prior("no_equals_here")

    def test_too_few_args(self) -> None:
        from tidal.inference._prior import parse_prior

        with pytest.raises(ValueError, match="DIST:ARG1:ARG2"):
            parse_prior("x=uniform:0.01")


class TestBuildPriorTransform:
    """Test joint prior transform for nested sampling."""

    def test_two_params(self) -> None:
        from tidal.inference._prior import Prior, build_prior_transform

        priors = [
            Prior("a", "uniform", 0.0, 10.0),
            Prior("b", "uniform", -1.0, 1.0),
        ]
        pt = build_prior_transform(priors)
        theta = pt(np.array([0.5, 0.5]))
        assert theta[0] == pytest.approx(5.0)
        assert theta[1] == pytest.approx(0.0)


# ===================================================================
# Constraint tests
# ===================================================================


class TestConstraints:
    """Test user-specifiable parameter constraints."""

    def test_simple_inequality(self) -> None:
        from tidal.inference._constraints import parse_constraint

        check = parse_constraint("xi > 0")
        assert check({"xi": 1.0})
        assert not check({"xi": -1.0})

    def test_compound_expression(self) -> None:
        from tidal.inference._constraints import parse_constraint

        check = parse_constraint("deltam**2 < 2*alpha*xi")
        assert check({"deltam": 0.5, "alpha": 1.0, "xi": 1.0})
        assert not check({"deltam": 2.0, "alpha": 1.0, "xi": 1.0})

    def test_function_call(self) -> None:
        from tidal.inference._constraints import parse_constraint

        check = parse_constraint("sqrt(alpha) > 1")
        assert check({"alpha": 4.0})
        assert not check({"alpha": 0.25})

    def test_unknown_parameter(self) -> None:
        from tidal.inference._constraints import ConstraintError, parse_constraint

        check = parse_constraint("xi > 0")
        with pytest.raises(ConstraintError, match="Unknown parameter"):
            check({"alpha": 1.0})  # xi not provided

    def test_not_a_comparison(self) -> None:
        from tidal.inference._constraints import ConstraintError, parse_constraint

        with pytest.raises(ConstraintError, match="comparison"):
            parse_constraint("xi + 1")

    def test_invalid_syntax(self) -> None:
        from tidal.inference._constraints import ConstraintError, parse_constraint

        with pytest.raises(ConstraintError, match="Invalid constraint syntax"):
            parse_constraint("xi >>>> 0")

    def test_constraint_set(self) -> None:
        from tidal.inference._constraints import ConstraintSet

        cs = ConstraintSet.from_strings(["xi > 0", "alpha > 0"])
        assert cs.check({"xi": 1.0, "alpha": 1.0})
        assert not cs.check({"xi": -1.0, "alpha": 1.0})
        assert not cs.check({"xi": 1.0, "alpha": -1.0})
        assert len(cs) == 2

    def test_empty_constraint_set(self) -> None:
        from tidal.inference._constraints import ConstraintSet

        cs = ConstraintSet.from_strings([])
        assert cs.check({"anything": 42.0})
        assert len(cs) == 0


# ===================================================================
# Likelihood config tests
# ===================================================================


class TestLikelihoodConfig:
    """Test likelihood configuration and parsing."""

    def test_parse_maximize(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        lc = parse_likelihood("P_max:maximize")
        assert lc.metric == "P_max"
        assert lc.likelihood_type == "maximize"

    def test_parse_gaussian(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        lc = parse_likelihood("P_max:gaussian:0.5:0.1")
        assert lc.metric == "P_max"
        assert lc.likelihood_type == "gaussian"
        assert lc.target == 0.5
        assert lc.sigma == 0.1

    def test_parse_threshold(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        lc = parse_likelihood("P_max:threshold:0.01")
        assert lc.metric == "P_max"
        assert lc.likelihood_type == "threshold"
        assert lc.min_value == 0.01

    def test_parse_too_short(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        with pytest.raises(ValueError, match="METRIC:TYPE"):
            parse_likelihood("P_max")

    def test_parse_gaussian_missing_args(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        with pytest.raises(ValueError, match="gaussian"):
            parse_likelihood("P_max:gaussian:0.5")

    def test_compute_gaussian(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(
            metric="P_max", likelihood_type="gaussian", target=1.0, sigma=0.1
        )
        # At target, log_L = 0
        assert compute_log_likelihood(1.0, lc) == pytest.approx(0.0)
        # 1 sigma away, log_L = -0.5
        assert compute_log_likelihood(1.1, lc) == pytest.approx(-0.5)

    def test_compute_threshold(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(
            metric="P_max", likelihood_type="threshold", min_value=0.1
        )
        assert compute_log_likelihood(0.5, lc) == 0.0
        assert compute_log_likelihood(0.05, lc) == -math.inf

    def test_compute_maximize(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert compute_log_likelihood(0.5, lc) == 0.5
        assert compute_log_likelihood(-1.0, lc) == -1.0

    def test_nan_returns_neg_inf(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert compute_log_likelihood(float("nan"), lc) == -math.inf


# ===================================================================
# InferenceResult tests
# ===================================================================


class TestInferenceResult:
    """Test inference results storage and analysis."""

    @pytest.fixture
    def mc_result(self) -> InferenceResult:
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(42)
        n = 50
        samples = rng.uniform(0, 10, size=(n, 2))
        log_l = -np.sum((samples - 5) ** 2, axis=1) / 10
        log_p = np.zeros(n)
        return InferenceResult(
            samples=samples,
            log_likelihood=log_l,
            log_prior=log_p,
            param_names=["alpha", "beta"],
            method="mc",
        )

    def test_n_samples(self, mc_result: InferenceResult) -> None:
        assert mc_result.n_samples == 50

    def test_n_params(self, mc_result: InferenceResult) -> None:
        assert mc_result.n_params == 2

    def test_log_posterior(self, mc_result: InferenceResult) -> None:
        lp = mc_result.log_posterior
        assert lp.shape == (50,)
        np.testing.assert_allclose(lp, mc_result.log_prior + mc_result.log_likelihood)

    def test_ess_unweighted(self, mc_result: InferenceResult) -> None:
        assert mc_result.effective_sample_size() == 50.0

    def test_ess_weighted(self) -> None:
        from tidal.inference._results import InferenceResult

        r = InferenceResult(
            samples=np.ones((4, 1)),
            log_likelihood=np.zeros(4),
            log_prior=np.zeros(4),
            param_names=["x"],
            method="nested",
            weights=np.array([1.0, 1.0, 1.0, 1.0]),
        )
        assert r.effective_sample_size() == pytest.approx(4.0)

        # Non-uniform weights: ESS < n
        r2 = InferenceResult(
            samples=np.ones((4, 1)),
            log_likelihood=np.zeros(4),
            log_prior=np.zeros(4),
            param_names=["x"],
            method="nested",
            weights=np.array([10.0, 0.0, 0.0, 0.0]),
        )
        assert r2.effective_sample_size() < 4.0

    def test_best(self, mc_result: InferenceResult) -> None:
        best = mc_result.best()
        assert "alpha" in best
        assert "beta" in best

    def test_posterior_mean(self, mc_result: InferenceResult) -> None:
        mean = mc_result.posterior_mean()
        assert "alpha" in mean
        assert "beta" in mean

    def test_credible_interval(self, mc_result: InferenceResult) -> None:
        ci = mc_result.credible_interval(0.95)
        assert "alpha" in ci
        lo, hi = ci["alpha"]
        assert lo < hi

    def test_to_sweep_results(self, mc_result: InferenceResult) -> None:
        sr = mc_result.to_sweep_results()
        assert len(sr.rows) == 50
        assert "log_likelihood" in sr.rows[0]
        assert "log_prior" in sr.rows[0]
        assert "alpha" in sr.rows[0]

    def test_save(self, mc_result: InferenceResult, tmp_path: Path) -> None:
        mc_result.save(tmp_path / "test_out")
        assert (tmp_path / "test_out" / "inference.json").exists()
        assert (tmp_path / "test_out" / "results.csv").exists()


# ===================================================================
# CLI help test
# ===================================================================


class TestSampleCLIHelp:
    """Test that tidal sample --help works."""

    def test_subparser_registered(self) -> None:
        """Verify the sample subparser is registered."""
        from tidal.cli import (
            _build_parser as build_parser,  # pyright: ignore[reportPrivateUsage]
        )

        parser = build_parser()
        # --help causes sys.exit, so just parse a minimal valid command
        args = parser.parse_args(["sample", "spec.json", "--output", "/tmp/test"])
        assert args.command == "sample"

    def test_parse_basic_args(self) -> None:
        """Verify basic argument parsing."""
        from tidal.cli import (
            _build_parser as build_parser,  # pyright: ignore[reportPrivateUsage]
        )

        parser = build_parser()
        args = parser.parse_args(
            [
                "sample",
                "spec.json",
                "--prior",
                "g0=uniform:0.01:0.5",
                "--likelihood",
                "P_max:maximize",
                "--method",
                "mc",
                "--n-samples",
                "10",
                "--output",
                "/tmp/test",
            ]
        )
        assert args.command == "sample"
        assert args.json_path == "spec.json"
        assert args.prior == ["g0=uniform:0.01:0.5"]
        assert args.likelihood == "P_max:maximize"
        assert args.method == "mc"
        assert args.n_samples == 10
