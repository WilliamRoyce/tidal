"""Tests for ``tidal.inference`` — Bayesian inference infrastructure."""

from __future__ import annotations

import logging
import math
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

if TYPE_CHECKING:
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

        with pytest.warns(UserWarning, match="ignores its bounds"):
            p = Prior("x", "arctan_uniform", -30.0, 30.0)
        rng = np.random.default_rng(42)
        samples = p.sample(rng, 1000)
        # Should have samples in a wide range including negative
        assert np.any(samples < -1.0)
        assert np.any(samples > 1.0)

    def test_arctan_uniform_transform(self) -> None:
        from tidal.inference._prior import Prior

        with pytest.warns(UserWarning, match="ignores its bounds"):
            p = Prior("x", "arctan_uniform", -30.0, 30.0)
        # u=0.5 maps to theta=0, tan(0)=0
        assert p.transform(0.5) == pytest.approx(0.0)
        # u=0 and u=1 should be at the extremes
        assert p.transform(0.0) < -10
        assert p.transform(1.0) > 10

    def test_arctan_uniform_warns_bounds_ignored(self) -> None:
        """GH #425: arctan_uniform ignores low/high — construction must say
        so instead of letting the bounds masquerade as the support.
        """
        import math as _math

        from tidal.inference._prior import _ARCTAN_EPS, Prior

        with pytest.warns(UserWarning, match="ignores its bounds"):
            Prior("x", "arctan_uniform", -89.0, 89.0)
        # Bounds matching the implied support do not warn.
        support = _math.tan(_math.pi / 2 - _ARCTAN_EPS)
        import warnings as _warnings

        with _warnings.catch_warnings():
            _warnings.simplefilter("error")
            Prior("x", "arctan_uniform", -support, support)
            Prior("x", "arctan_uniform", 0.0, 0.0)  # sanctioned sentinel
            Prior("x", "uniform", -89.0, 89.0)  # other kinds never warn

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

    def test_parse_likelihood_propagates_baseline_for_maximize(self) -> None:
        """Regression for #319: ``--baseline-formula`` must reach
        ``LikelihoodConfig`` for ``maximize`` (was silently dropped).
        """
        from tidal.inference._likelihood import parse_likelihood

        formula = "sin(kappa*B0*t_end/2)**2"
        lc = parse_likelihood("P_max:maximize", baseline_formula=formula)
        assert lc.likelihood_type == "maximize"
        assert lc.baseline_formula == formula

    def test_parse_likelihood_propagates_baseline_for_minimize(self) -> None:
        """Regression for #319: same propagation for ``minimize``."""
        from tidal.inference._likelihood import parse_likelihood

        formula = "sin(kappa*B0*t_end/2)**2"
        lc = parse_likelihood("P_max:minimize", baseline_formula=formula)
        assert lc.likelihood_type == "minimize"
        assert lc.baseline_formula == formula

    def test_parse_likelihood_extremize_baseline_unchanged(self) -> None:
        """Sanity: ``extremize`` already worked before #319, ensure the fix
        didn't regress it.
        """
        from tidal.inference._likelihood import parse_likelihood

        formula = "0.05"
        lc = parse_likelihood("P_max:extremize", baseline_formula=formula)
        assert lc.likelihood_type == "extremize"
        assert lc.baseline_formula == formula

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
            metric="P_max",
            likelihood_type="gaussian",
            target=1.0,
            sigma=0.1,
        )
        # At target, log_L = 0
        assert compute_log_likelihood(1.0, lc) == pytest.approx(0.0)
        # 1 sigma away, log_L = -0.5
        assert compute_log_likelihood(1.1, lc) == pytest.approx(-0.5)

    def test_compute_threshold(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="threshold",
            min_value=0.1,
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
            ],
        )
        assert args.command == "sample"
        assert args.json_path == "spec.json"
        assert args.prior == ["g0=uniform:0.01:0.5"]
        assert args.likelihood == "P_max:maximize"
        assert args.method == "mc"
        assert args.n_samples == 10

    def test_parse_new_flags(self) -> None:
        """Verify new CLI flags (analyze, nlive-auto, importance, extremize)."""
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
                "P_max:extremize",
                "--baseline-formula",
                "sin(kappa * B0 * t_end / 2)**2",
                "--method",
                "nested",
                "--sampler",
                "polychord",
                "--nlive-auto",
                "production",
                "--analyze",
                "--importance",
                "--output",
                "/tmp/test",
            ],
        )
        assert args.likelihood == "P_max:extremize"
        assert args.baseline_formula == "sin(kappa * B0 * t_end / 2)**2"
        assert args.sampler == "polychord"
        assert args.nlive_auto == "production"
        assert args.analyze is True
        assert args.importance is True


# ===================================================================
# New likelihood types
# ===================================================================


class TestNewLikelihoodTypes:
    """Test minimize and extremize likelihood types."""

    def test_parse_minimize(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        lc = parse_likelihood("P_max:minimize")
        assert lc.likelihood_type == "minimize"

    def test_parse_extremize(self) -> None:
        from tidal.inference._likelihood import parse_likelihood

        lc = parse_likelihood(
            "P_max:extremize",
            baseline_formula="sin(kappa * B0 * t_end / 2)**2",
        )
        assert lc.likelihood_type == "extremize"
        assert lc.baseline_formula is not None

    def test_extremize_requires_formula(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig

        with pytest.raises(ValueError, match="baseline-formula"):
            LikelihoodConfig(metric="P_max", likelihood_type="extremize")

    def test_compute_minimize(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(metric="P_max", likelihood_type="minimize")
        assert compute_log_likelihood(0.5, lc) == -0.5
        assert compute_log_likelihood(2.0, lc) == -2.0

    def test_compute_extremize(self) -> None:
        from tidal.inference._likelihood import LikelihoodConfig, compute_log_likelihood

        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="extremize",
            baseline_formula="sin(kappa * B0 * t_end / 2)**2",
        )
        params = {"kappa": 1.0, "B0": 0.01, "t_end": 50.0}
        baseline = math.sin(1.0 * 0.01 * 50.0 / 2) ** 2

        # At baseline: logL = |log(1)| = 0
        assert compute_log_likelihood(baseline, lc, params) == pytest.approx(0.0)

        # 2x amplification: logL = log(2) ≈ 0.693
        assert compute_log_likelihood(2 * baseline, lc, params) == pytest.approx(
            math.log(2),
        )

        # 0.5x suppression: logL = |log(0.5)| = log(2) (symmetric)
        assert compute_log_likelihood(0.5 * baseline, lc, params) == pytest.approx(
            math.log(2),
        )


class TestBaselineFormulaSpecParams:
    """Regression for #270: baseline-formula must see spec-level params.

    Before the fix, ``_evaluate_likelihood`` built ``eval_params`` only
    from the swept theta + t_end/dt.  A formula like
    ``sin(kappa * B0 * t_end / 2)**2`` then raised NameError on
    ``kappa``, ``_eval_baseline`` silently returned ``None``, and every
    likelihood collapsed to ``-inf``.  PolyChord then spun forever.
    """

    def test_parse_params_merges_cli_overrides(self) -> None:
        # `_parse_params` is the canonical merge used by both the
        # simulator and (after #270) the likelihood.  Verify it surfaces
        # CLI --param values so baseline formulas can resolve them.
        from tidal.cli._simulate import (  # pyright: ignore[reportPrivateUsage]
            _parse_params,
        )
        from tidal.symbolic import load_equation_system

        repo_root = Path(__file__).resolve().parents[1]
        spec = load_equation_system(
            repo_root / "examples" / "data" / "dark_photon_plasma.json",
        )
        params = _parse_params(["kappa=1.0", "B0=0.1"], spec)
        assert params["kappa"] == 1.0
        assert params["B0"] == 0.1

    def test_baseline_formula_resolves_with_spec_params(self) -> None:
        # Direct test of the formula-eval namespace the fix assembles.
        from tidal.inference._likelihood import (
            LikelihoodConfig,
            _eval_baseline,  # pyright: ignore[reportPrivateUsage]
            compute_log_likelihood,
        )

        cfg = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            baseline_formula="sin(kappa * B0 * t_end / 2)**2",
        )
        # After #270: eval_params contains spec/CLI params + swept theta.
        eval_params = {
            "kappa": 1.0,
            "B0": 0.1,
            "t_end": 10.0,
            "mA2": 0.5,
            "deltam": 0.0,
        }
        baseline = _eval_baseline(cfg.baseline_formula, eval_params)
        assert baseline is not None
        assert baseline == pytest.approx(math.sin(0.5) ** 2)

        # logL should be finite (not -inf) — this is the regression check.
        logl = compute_log_likelihood(0.01, cfg, eval_params)
        assert math.isfinite(logl)


class TestPolyChordGuard:
    """Regression for #270: PolyChord must fail fast on all-(-inf) likelihoods.

    Without the pre-flight guard, PolyChord spins indefinitely trying to
    build its initial live-point pool.  Job 28002101 hit the 1-hour wall
    with zero chain output because of this.
    """

    def test_all_inf_raises_runtime_error(self) -> None:
        from tidal.inference._nested import (
            _run_polychord,  # pyright: ignore[reportPrivateUsage]
        )

        def always_inf(_theta: list[float]) -> float:
            return float("-inf")

        def unit_prior(hypercube: list[float]) -> list[float]:
            return list(hypercube)

        with pytest.raises(RuntimeError, match="non-finite logL"):
            _run_polychord(
                log_likelihood=always_inf,
                prior_transform=unit_prior,
                ndim=2,
                param_names=["a", "b"],
                nlive=25,
                n_workers=None,
                quiet=True,
                n_probe=5,
            )


class TestLikelihoodBackendEquivalence:
    """Memory and disk backends must produce identical metrics (#269).

    If these ever diverge, the in-memory fast path has drifted from the
    disk-writer reference semantics and callers could get wrong science.
    """

    def test_cache_warm_vs_cold_identical_metrics(self, tmp_path: Path) -> None:
        """Regression for #291: solver-layer caches (parsed-Mathematica,
        eval namespace) must NOT change numerical results whether the
        cache is cold (first call) or warm (subsequent).
        """
        import sys

        import tidal.cli

        spec_path = (
            Path(__file__).resolve().parents[1]
            / "examples"
            / "data"
            / "dark_photon_plasma.json"
        )
        if not spec_path.exists():
            pytest.skip(f"dark_photon_plasma spec not found at {spec_path}")

        argv = [
            "tidal",
            "sample",
            str(spec_path),
            "--param",
            "kappa=1.0",
            "--param",
            "B0=0.01",
            "--param",
            "alpha3=0.25",
            "--param",
            "xi=1.0",
            "--param",
            "deltam=0.3",
            "--prior",
            "mA2=uniform:0.4:0.5",
            "--likelihood",
            "P_max:maximize",
            "--method",
            "nested",
            "--sampler",
            "polychord",
            "--grid-shape",
            "32",
            "--bounds",
            "0:50",
            "--periodic",
            "--ic",
            "plane-wave",
            "--ic-component",
            "h_5",
            "--ic-wavevector",
            "2.0",
            "--ic-amplitude",
            "1e-2",
            "--source",
            "h_5",
            "--target",
            "a_1",
            "--snapshots",
            "3",
            "--t-end",
            "10.0",
            "--output",
            str(tmp_path / "ns"),
        ]
        old_argv = sys.argv
        try:
            sys.argv = argv
            parser = tidal.cli._build_parser()  # pyright: ignore[reportPrivateUsage]
            args = parser.parse_args(argv[1:])
        finally:
            sys.argv = old_argv

        from tidal.cli._sweep import (
            _measure_from_sim_data,  # pyright: ignore[reportPrivateUsage]
            run_inference_step,
        )
        from tidal.symbolic import load_equation_system
        from tidal.symbolic._eval_utils import (
            _PARSED_MATH_CACHE,  # pyright: ignore[reportPrivateUsage]
        )

        spec = load_equation_system(spec_path)

        # Cold run — explicitly clear the parsed-Mathematica cache so
        # the first call exercises the miss path.
        _PARSED_MATH_CACHE.clear()
        overrides = {"mA2": 0.42}
        sim_cold = run_inference_step(args, spec_path, overrides, spec=spec)
        m_cold = _measure_from_sim_data(
            sim_cold,
            {"conversion", "peak_conversion"},
            ("h_5",),
            ("a_1",),
            0.5,
        )
        assert len(_PARSED_MATH_CACHE) > 0, "cold run should populate the cache"

        # Warm run at same theta — must produce identical metrics.
        sim_warm = run_inference_step(args, spec_path, overrides, spec=spec)
        m_warm = _measure_from_sim_data(
            sim_warm,
            {"conversion", "peak_conversion"},
            ("h_5",),
            ("a_1",),
            0.5,
        )

        for key in ("P_max", "P_final"):
            assert m_cold.get(key) is not None, f"cold run missing {key}"
            assert m_warm.get(key) is not None, f"warm run missing {key}"
            assert m_cold[key] == pytest.approx(m_warm[key], rel=0, abs=0), (
                f"cache warmth changed {key}: cold={m_cold[key]}, warm={m_warm[key]}"
            )

    @pytest.mark.parametrize("metric", ["P_max", "P_final"])
    def test_backends_agree(self, metric: str, tmp_path: Path) -> None:
        import sys

        spec_path = (
            Path(__file__).resolve().parents[1]
            / "examples"
            / "data"
            / "dark_photon_plasma.json"
        )
        if not spec_path.exists():
            pytest.skip(f"dark_photon_plasma spec not found at {spec_path}")

        # Build a real Namespace via the tidal argument parser so every
        # simulate-facing attribute gets the correct type/default.
        import tidal.cli

        argv = [
            "tidal",
            "sample",
            str(spec_path),
            "--param",
            "kappa=1.0",
            "--param",
            "B0=0.01",
            "--param",
            "alpha3=0.25",
            "--param",
            "xi=1.0",
            "--param",
            "deltam=0.3",
            "--prior",
            "mA2=uniform:0.4:0.5",  # single point we evaluate manually
            "--likelihood",
            "P_max:maximize",
            "--method",
            "nested",
            "--sampler",
            "polychord",
            "--grid-shape",
            "32",
            "--bounds",
            "0:50",
            "--periodic",
            "--ic",
            "plane-wave",
            "--ic-component",
            "h_5",
            "--ic-wavevector",
            "2.0",
            "--ic-amplitude",
            "1e-2",
            "--source",
            "h_5",
            "--target",
            "a_1",
            "--snapshots",
            "3",
            "--t-end",
            "10.0",
            "--output",
            str(tmp_path / "ns"),
        ]
        old_argv = sys.argv
        try:
            sys.argv = argv
            parser = tidal.cli._build_parser()  # pyright: ignore[reportPrivateUsage]
            args = parser.parse_args(argv[1:])
        finally:
            sys.argv = old_argv

        from tidal.cli._sweep import (
            _measure_from_sim_data,  # pyright: ignore[reportPrivateUsage]
            _measure_run,  # pyright: ignore[reportPrivateUsage]
            run_inference_step,
            simulate_run,
        )
        from tidal.symbolic import load_equation_system

        spec = load_equation_system(spec_path)
        overrides = {"mA2": 0.45}

        # Memory path
        sim_data = run_inference_step(args, spec_path, overrides, spec=spec)
        mem_metrics = _measure_from_sim_data(
            sim_data,
            {"conversion", "peak_conversion"},
            ("h_5",),
            ("a_1",),
            0.5,
        )

        # Disk path — uses the same fixed config.
        disk_dir = tmp_path / "disk_run"
        disk_dir.mkdir()
        exit_code, _, _ = simulate_run(args, spec_path, overrides, disk_dir, spec=spec)
        assert exit_code == 0
        disk_metrics = _measure_run(
            disk_dir,
            spec_path,
            {"conversion", "peak_conversion"},
            ("h_5",),
            ("a_1",),
            0.5,
            spec=spec,
        )

        mem_val = mem_metrics.get(metric)
        disk_val = disk_metrics.get(metric)
        assert mem_val is not None
        assert disk_val is not None
        assert mem_val == pytest.approx(disk_val, rel=1e-12, abs=1e-15)


# ===================================================================
# recommend_nlive
# ===================================================================


class TestRecommendNlive:
    """Test nlive auto-scaling."""

    def test_fast(self) -> None:
        from tidal.inference._nested import recommend_nlive

        assert recommend_nlive(4, "fast") == 100
        assert recommend_nlive(1, "fast") == 50

    def test_standard(self) -> None:
        from tidal.inference._nested import recommend_nlive

        assert recommend_nlive(4, "standard") == 100
        assert recommend_nlive(10, "standard") == 250

    def test_production(self) -> None:
        from tidal.inference._nested import recommend_nlive

        assert recommend_nlive(4, "production") == 250
        assert recommend_nlive(10, "production") == 500


# ===================================================================
# Parameter importance
# ===================================================================


class TestParameterImportance:
    """Test parameter importance analysis via anesthetic."""

    @pytest.fixture
    def nested_result(self) -> InferenceResult:
        """Synthetic nested sampling result with one constrained and one
        unconstrained parameter.
        """
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(42)
        n = 300
        nlive = 50

        # x is tightly constrained (narrow posterior), y is broad
        x = rng.normal(0, 0.1, n)
        y = rng.uniform(-5, 5, n)
        samples = np.column_stack([x, y])
        log_l = -0.5 * (x / 0.1) ** 2  # only depends on x

        # Build approximate logL_birth for anesthetic
        sorted_idx = np.argsort(log_l)
        logl_birth = np.full(n, -np.inf)
        sorted_logl = log_l[sorted_idx]
        for i in range(nlive, n):
            logl_birth[sorted_idx[i]] = sorted_logl[i - nlive]

        return InferenceResult(
            samples=samples,
            log_likelihood=log_l,
            log_prior=np.zeros(n),
            param_names=["x_constrained", "y_unconstrained"],
            method="nested",
            log_evidence=-1.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={"sampler": "test", "nlive": nlive},
        )

    def test_to_anesthetic(self, nested_result: InferenceResult) -> None:
        pytest.importorskip("anesthetic")
        ns = nested_result.to_anesthetic()
        assert hasattr(ns, "D_KL")
        assert hasattr(ns, "d_G")
        assert hasattr(ns, "logZ")

    def test_importance_computes(self, nested_result: InferenceResult) -> None:
        """Parameter importance metrics are finite and sensible."""
        pytest.importorskip("anesthetic")
        imp = nested_result.parameter_importance(n_bootstrap=10)

        assert math.isfinite(imp.d_kl)
        assert imp.d_kl >= 0
        assert math.isfinite(imp.d_g)
        assert imp.d_g >= 0
        assert len(imp.marginal_d_kl) == 2
        for dkl in imp.marginal_d_kl.values():
            assert math.isfinite(dkl)

    def test_format_table(self, nested_result: InferenceResult) -> None:
        pytest.importorskip("anesthetic")
        from tidal.inference._importance import format_importance_table

        imp = nested_result.parameter_importance(n_bootstrap=10)
        table = format_importance_table(imp)
        assert "x_constrained" in table
        assert "D_KL" in table
        assert "d_G" in table

    def test_cross_kl(self, nested_result: InferenceResult) -> None:
        """Cross-KL is zero between a chain and itself; positive between
        different chains.  Identifies which parameters distinguish the two
        posteriors.
        """
        pytest.importorskip("anesthetic")
        import numpy as np

        from tidal.inference._importance import compute_cross_kl

        ns = nested_result.to_anesthetic()
        params = list(nested_result.param_names)

        # Self-cross-KL must be 0 (modulo histogram quantization noise).
        self_cross = compute_cross_kl(ns, ns, params, n_bins=40)
        for name in params:
            assert math.isfinite(self_cross[name])
            assert self_cross[name] < 0.05, f"self cross-KL for {name} should be ~0"

        # Cross-KL against a shifted copy of the same chain is positive.
        # Build a shifted version by adding a constant to one column.
        ns_shift = ns.copy()
        shift_param = params[0]
        shift = 2.0 * float(np.std(np.asarray(ns[shift_param])))
        ns_shift[shift_param] = np.asarray(ns_shift[shift_param]) + shift
        cross = compute_cross_kl(ns, ns_shift, params, n_bins=40)
        assert cross[shift_param] > 0.1, (
            f"shifted cross-KL for {shift_param} should be > 0.1, got {cross[shift_param]}"
        )

    def test_from_directory_roundtrip(
        self,
        nested_result: InferenceResult,
        tmp_path: Path,
    ) -> None:
        from tidal.inference._results import InferenceResult

        out = tmp_path / "ns_out"
        nested_result.save(out)

        loaded = InferenceResult.from_directory(out)
        assert loaded.n_samples == nested_result.n_samples
        assert loaded.param_names == nested_result.param_names
        assert loaded.log_evidence == nested_result.log_evidence
        np.testing.assert_allclose(loaded.samples, nested_result.samples)

    def test_marginal_dkl_with_integer_indexed_columns(
        self,
        nested_result: InferenceResult,
    ) -> None:
        """Regression for #287: anesthetic's read_chains returns integer-
        indexed parameter columns in v2.0+ (``[0, 1, 'logL', ...]``), but
        our importance loop previously indexed by parameter name, silently
        NaN'ing the marginal D_KL for every parameter.  Verify both
        integer- and name-indexed anesthetic samples give the same
        finite D_KL.
        """
        pytest.importorskip("anesthetic")
        from anesthetic import NestedSamples  # type: ignore[import-untyped]

        from tidal.inference._importance import compute_parameter_importance

        # Named-column path (manual reconstruction — matches
        # to_anesthetic_samples when chain_root is absent).
        imp_named = compute_parameter_importance(nested_result, n_bootstrap=10)
        for dkl in imp_named.marginal_d_kl.values():
            assert math.isfinite(dkl), "named-column path regressed"

        # Integer-indexed path — simulate what read_chains returns in
        # v2.0+ and monkeypatch to_anesthetic_samples to return it.
        n = nested_result.n_samples
        nlive = int(nested_result.metadata["nlive"])
        sorted_idx = np.argsort(nested_result.log_likelihood)
        sorted_logl = nested_result.log_likelihood[sorted_idx]
        logl_birth = np.full(n, -np.inf)
        for i in range(nlive, n):
            logl_birth[sorted_idx[i]] = sorted_logl[i - nlive]

        ns_int = NestedSamples(
            data=nested_result.samples,
            logL=nested_result.log_likelihood,
            logL_birth=logl_birth,
            columns=[0, 1],  # integer-indexed, mimics read_chains
        )
        import tidal.inference._importance as imp_mod

        original = imp_mod.to_anesthetic_samples
        try:
            imp_mod.to_anesthetic_samples = lambda _r: ns_int  # type: ignore[assignment]
            imp_int = compute_parameter_importance(nested_result, n_bootstrap=10)
        finally:
            imp_mod.to_anesthetic_samples = original

        for name, dkl in imp_int.marginal_d_kl.items():
            assert math.isfinite(dkl), (
                f"integer-column path regressed for '{name}' (got NaN)"
            )

        # Both paths should agree on the per-parameter ranking to within
        # ~50% (bootstrap noise is large at n_bootstrap=10).
        for name in nested_result.param_names:
            assert imp_named.marginal_d_kl[name] == pytest.approx(
                imp_int.marginal_d_kl[name],
                rel=0.5,
                abs=0.1,
            )


# ===================================================================
# Marginal D_KL prior transforms (#420)
# ===================================================================


class TestMarginalDKLPriorTransforms:
    """Regression tests for #420: marginal D_KL must transform each column
    into the space where its prior is uniform for EVERY supported prior
    kind, not just log_uniform.  The null test is exact: samples drawn
    from the prior itself have true marginal D_KL = 0, so any sizable
    reported value is estimator artifact (pre-fix: ~2.5 nats for
    arctan_uniform).
    """

    #: Histogram-KL bias for N samples in K bins is ~(K-1)/(2N) ≈ 1e-4
    #: at N=200k, K=40; 0.02 leaves two orders of magnitude of headroom.
    NULL_TOL = 0.02

    @pytest.mark.parametrize(
        ("kind", "lo", "hi"),
        [
            ("uniform", -1.0, 1.0),
            ("log_uniform", 1e-3, 1e3),
            ("arctan_uniform", -89.0, 89.0),
            ("normal", 0.0, 1.0),
        ],
    )
    def test_null_prior_equals_posterior(self, kind: str, lo: float, hi: float) -> None:
        """Samples from the prior report marginal D_KL ≈ 0 for every kind."""
        from tidal.inference._importance import (
            _hist_kl_vs_uniform,
            _uniformizing_transform,
        )
        from tidal.inference._prior import Prior

        rng = np.random.default_rng(42)
        n = 200_000
        samples = Prior(name="p", distribution=kind, low=lo, high=hi).sample(rng, n)
        transform = _uniformizing_transform(kind, lo, hi)
        assert transform is not None, f"no transform for supported kind {kind}"
        fn, a, b = transform
        kl, coverage = _hist_kl_vs_uniform(fn(samples), np.full(n, 1.0 / n), a, b)
        assert coverage == pytest.approx(1.0), "prior samples must lie in range"
        assert kl < self.NULL_TOL, f"{kind}: null D_KL {kl:.4f} not ~0"

    def test_positive_control_uniform(self) -> None:
        """A genuinely constrained posterior reports the analytic D_KL."""
        from tidal.inference._importance import (
            _hist_kl_vs_uniform,
            _uniformizing_transform,
        )

        rng = np.random.default_rng(3)
        lo, hi = -1.0, 1.0
        sigma = (hi - lo) / 20
        post = rng.normal(0.0, sigma, 200_000)
        post = post[(post > lo) & (post < hi)]
        transform = _uniformizing_transform("uniform", lo, hi)
        assert transform is not None
        fn, a, b = transform
        kl, _coverage = _hist_kl_vs_uniform(
            fn(post), np.full(len(post), 1.0 / len(post)), a, b
        )
        analytic = math.log((hi - lo) / (sigma * math.sqrt(2 * math.pi * math.e)))
        assert kl > 0.5
        assert kl == pytest.approx(analytic, abs=0.15)

    def test_arctan_transform_ignores_recorded_bounds(self) -> None:
        """Prior.sample ignores low/high for arctan_uniform (fixed
        eps-truncated theta range), so the uniformizing transform must
        too — using the recorded (-89, 89) as the histogram range is
        the original #420 bug.
        """
        from tidal.inference._importance import _uniformizing_transform
        from tidal.inference._prior import _ARCTAN_EPS

        t1 = _uniformizing_transform("arctan_uniform", -89.0, 89.0)
        t2 = _uniformizing_transform("arctan_uniform", -30.0, 30.0)
        assert t1 is not None
        assert t2 is not None
        assert (
            t1[1:]
            == t2[1:]
            == (
                -(math.pi / 2 - _ARCTAN_EPS),
                math.pi / 2 - _ARCTAN_EPS,
            )
        )

    # -- radial_angular ------------------------------------------------

    @staticmethod
    def _ra_record(q: np.ndarray) -> dict:
        return {
            "kind": "radial_angular",
            "names": ["c1", "c2", "c3", "c4"],
            "r_lo": 1e-3,
            "r_hi": 1e3,
            "face_idx": 1,
            "sub_tile": [1, 2, 1],
            "M": 2,
            "Q": np.asarray(q).tolist(),
        }

    @pytest.mark.parametrize("rotate", [False, True])
    def test_vectorized_sampler_matches_transform(self, rotate: bool) -> None:
        """_sample_radial_angular reproduces RadialAngularPrior.transform."""
        from tidal.inference._importance import _sample_radial_angular
        from tidal.inference._prior import RadialAngularPrior
        from tidal.inference._sphere import random_rotation

        q = random_rotation(4, seed=5) if rotate else np.eye(4)
        record = self._ra_record(q)
        prior = RadialAngularPrior(
            names=("c1", "c2", "c3", "c4"),
            r_lo=1e-3,
            r_hi=1e3,
            face_idx=1,
            sub_tile=(1, 2, 1),
            M=2,
            Q=q,
        )
        u_test = np.random.default_rng(7).random((100, 4))

        class _FixedRng:
            def random(self, shape: tuple[int, ...]) -> np.ndarray:
                assert shape == u_test.shape
                return u_test

        vec = _sample_radial_angular(record, 100, _FixedRng())  # type: ignore[arg-type]
        ref = np.array([prior.transform(u) for u in u_test])
        np.testing.assert_allclose(vec, ref, rtol=0, atol=1e-12)

    def test_radial_angular_null(self) -> None:
        """Posterior drawn from the radial_angular prior → D_KL ≈ 0 per
        coupling (pre-fix: these columns fell to a self-referential
        uniform fallback because the joint record has no 'name' key).
        """
        from tidal.inference._importance import (
            _RA_REFERENCE_SAMPLES,
            _RA_REFERENCE_SEED,
            _marginal_dkl_empirical,
            _sample_radial_angular,
        )

        record = self._ra_record(np.eye(4))
        reference = _sample_radial_angular(
            record, _RA_REFERENCE_SAMPLES, np.random.default_rng(_RA_REFERENCE_SEED)
        )
        n_post = 20_000
        posterior = _sample_radial_angular(record, n_post, np.random.default_rng(11))
        weights = np.full(n_post, 1.0 / n_post)
        for j in range(4):
            kl = _marginal_dkl_empirical(
                posterior[:, j], weights, reference[:, j], record["r_lo"]
            )
            assert kl < 0.05, f"c{j + 1}: null D_KL {kl:.4f} not ~0"

    def test_radial_angular_positive_control(self) -> None:
        """Constraining r to the top decade produces a clearly positive
        D_KL (guards against the metric collapsing to useless-zero,
        e.g. from linear binning over the 6-decade log-uniform range).
        """
        from tidal.inference._importance import (
            _RA_REFERENCE_SAMPLES,
            _RA_REFERENCE_SEED,
            _marginal_dkl_empirical,
            _sample_radial_angular,
        )

        record = self._ra_record(np.eye(4))
        reference = _sample_radial_angular(
            record, _RA_REFERENCE_SAMPLES, np.random.default_rng(_RA_REFERENCE_SEED)
        )
        n_post = 20_000
        u = np.random.default_rng(13).random((n_post, 4))
        u[:, 0] = 5 / 6 + u[:, 0] / 6  # top decade of the 6-decade r range

        class _FixedRng:
            def random(self, shape: tuple[int, ...]) -> np.ndarray:
                assert shape == u.shape
                return u

        posterior = _sample_radial_angular(record, n_post, _FixedRng())  # type: ignore[arg-type]
        weights = np.full(n_post, 1.0 / n_post)
        kl = _marginal_dkl_empirical(
            posterior[:, 0], weights, reference[:, 0], record["r_lo"]
        )
        assert kl > 0.5

    # -- prior_map handling --------------------------------------------

    def test_zero_low_bound_not_dropped(self, caplog: pytest.LogCaptureFixture) -> None:
        """Object-path priors with low == 0.0 must land in the prior map
        (pre-fix: `getattr(p, 'low', None) or ...` dropped falsy 0.0
        into the empirical fallback).  The fallback-with-metadata path
        now warns, so a clean log is the discriminator (the KL value
        itself barely differs: the empirical range ≈ the true bounds).
        """
        pytest.importorskip("anesthetic")
        import logging

        from tidal.inference._importance import compute_parameter_importance
        from tidal.inference._prior import Prior
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(9)
        n, nlive = 3000, 500
        prior = Prior(name="x", distribution="uniform", low=0.0, high=5.0)
        x = prior.sample(rng, n)
        log_l = np.zeros(n)
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=log_l,
            log_prior=np.zeros(n),
            param_names=["x"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            # Prior OBJECT (not dict) in metadata: exercises the getattr
            # path where the falsy-zero `or`-chain bug lived.
            metadata={"sampler": "test", "nlive": nlive, "priors": [prior]},
        )
        with caplog.at_level(logging.WARNING, logger="tidal.inference"):
            imp = compute_parameter_importance(result, n_bootstrap=10)
        assert "no usable prior record" not in caplog.text, (
            "low=0.0 prior fell through to the empirical fallback"
        )
        assert math.isfinite(imp.marginal_d_kl["x"])

    def test_integration_with_priors_metadata(self) -> None:
        """End-to-end through compute_parameter_importance with recorded
        priors metadata: an arctan_uniform column of prior-distributed
        samples reports ≈ 0 (pre-fix: ~2.4 nats), and radial_angular
        columns are finite (pre-fix path made them fallback values).

        Tolerances are looser than the helper-level nulls because the
        anesthetic reconstruction assigns geometric nested-sampling
        weights even for constant logL, capping the effective sample
        size at ~2*nlive; the histogram-KL bias is ~n_bins/(2*N_eff).
        """
        pytest.importorskip("anesthetic")
        from tidal.inference._importance import compute_parameter_importance
        from tidal.inference._prior import Prior
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(21)
        n, nlive = 6000, 1000
        arctan = Prior(
            name="alpha", distribution="arctan_uniform", low=-89.0, high=89.0
        )
        a_col = arctan.sample(rng, n)
        from tidal.inference._importance import _sample_radial_angular

        record = self._ra_record(np.eye(4))
        ra_cols = _sample_radial_angular(record, n, rng)
        samples = np.column_stack([a_col, ra_cols])
        log_l = np.zeros(n)  # posterior == prior
        result = InferenceResult(
            samples=samples,
            log_likelihood=log_l,
            log_prior=np.zeros(n),
            param_names=["alpha", "c1", "c2", "c3", "c4"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={
                "sampler": "test",
                "nlive": nlive,
                "priors": [
                    {
                        "kind": "scalar",
                        "name": "alpha",
                        "distribution": "arctan_uniform",
                        "low": -89.0,
                        "high": 89.0,
                    },
                    record,
                ],
            },
        )
        imp = compute_parameter_importance(result, n_bootstrap=10)
        assert imp.marginal_d_kl["alpha"] < 0.1, (
            f"arctan null reported {imp.marginal_d_kl['alpha']:.3f} nats "
            "(pre-#420 bug reported ~2.4)"
        )
        for name in ("c1", "c2", "c3", "c4"):
            assert math.isfinite(imp.marginal_d_kl[name])
            assert imp.marginal_d_kl[name] < 0.2
        # Consistency block present and healthy for a null posterior.
        assert imp.consistency
        assert imp.consistency["superadditivity_ok"]
        assert imp.consistency["saturated_params"] == []

    # -- prevention devices (post-merge review hardening) ---------------

    def test_every_valid_kind_has_uniformizing_transform(self) -> None:
        """#420-recurrence insurance: adding a distribution kind to
        VALID_DISTRIBUTIONS without teaching _uniformizing_transform about
        it must fail the suite, not silently fall back to the
        sample-range reference.
        """
        from tidal.inference._importance import _uniformizing_transform
        from tidal.inference._prior import VALID_DISTRIBUTIONS

        # Representative valid hyperparameters per kind.
        args = {
            "uniform": (-1.0, 1.0),
            "log_uniform": (1e-3, 1e3),
            "arctan_uniform": (0.0, 0.0),
            "normal": (0.0, 1.0),
        }
        assert set(args) == set(VALID_DISTRIBUTIONS), (
            "a prior kind was added or removed: extend this mapping AND "
            "_uniformizing_transform (see tidal/inference/_prior.py "
            "VALID_DISTRIBUTIONS)"
        )
        for kind in VALID_DISTRIBUTIONS:
            lo, hi = args[kind]
            assert _uniformizing_transform(kind, lo, hi) is not None, (
                f"no uniformizing transform for supported prior kind "
                f"'{kind}' — this reintroduces #420"
            )

    def test_superadditivity_check_fires_on_mis_referenced_prior(self) -> None:
        """Behavior test for the self-check itself: a deliberately wrong
        prior record (uniform over (-89, 89) claimed for arctan-shaped
        samples — the literal pre-#420 configuration) must produce
        superadditivity_ok == False at healthy n_eff.  Weakening the check
        (e.g. tolerance = inf) now fails a test instead of only deleting it.
        """
        pytest.importorskip("anesthetic")
        from tidal.inference._importance import compute_parameter_importance
        from tidal.inference._prior import Prior
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(31)
        n, nlive = 6000, 1000
        with pytest.warns(UserWarning, match="ignores its bounds"):
            arctan = Prior("a", "arctan_uniform", -89.0, 89.0)
        cols = [arctan.sample(rng, n) for _ in range(3)]
        samples = np.column_stack(cols)
        result = InferenceResult(
            samples=samples,
            log_likelihood=np.zeros(n),
            log_prior=np.zeros(n),
            param_names=["a1", "a2", "a3"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={
                "sampler": "test",
                "nlive": nlive,
                # WRONG on purpose: claims uniform(-89, 89) for arctan
                # samples — each marginal then reads ~2.4 nats while the
                # joint stays ~0, the exact pre-#420 signature.
                "priors": [
                    {
                        "kind": "scalar",
                        "name": nm,
                        "distribution": "uniform",
                        "low": -89.0,
                        "high": 89.0,
                    }
                    for nm in ("a1", "a2", "a3")
                ],
            },
        )
        imp = compute_parameter_importance(result, n_bootstrap=10)
        cons = imp.consistency
        assert cons["superadditivity_applicable"] is True
        assert cons["superadditivity_ok"] is False
        assert cons["sum_marginals"] > cons["joint_d_kl"] + 1.0

    def test_range_clipped_prior_is_flagged(self) -> None:
        """A recorded prior range narrower than the samples (mis-assigned
        record, the #434 class) must be flagged, not silently renormalized
        over the retained subset.
        """
        pytest.importorskip("anesthetic")
        from tidal.inference._importance import (
            compute_parameter_importance,
            format_importance_table,
        )
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(33)
        n, nlive = 6000, 1000
        x = rng.uniform(-2.0, 2.0, n)  # true support (-2, 2)
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=np.zeros(n),
            log_prior=np.zeros(n),
            param_names=["x"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={
                "sampler": "test",
                "nlive": nlive,
                # Recorded range covers only half the samples.
                "priors": [
                    {
                        "kind": "scalar",
                        "name": "x",
                        "distribution": "uniform",
                        "low": 0.0,
                        "high": 2.0,
                    }
                ],
            },
        )
        imp = compute_parameter_importance(result, n_bootstrap=10)
        clipped = imp.consistency["range_clipped"]
        assert "x" in clipped
        assert clipped["x"] == pytest.approx(0.5, abs=0.1)
        assert "outside the recorded prior range" in format_importance_table(imp)

    def test_zero_coverage_returns_nan_not_zero(self) -> None:
        """All posterior mass outside the recorded range must yield NaN,
        never a silent 0.0 (review finding H1).
        """
        from tidal.inference._importance import _hist_kl_vs_uniform

        w = np.full(100, 0.01)
        kl, coverage = _hist_kl_vs_uniform(np.full(100, 5.0), w, 0.0, 1.0)
        assert coverage == 0.0
        assert math.isnan(kl)

    def test_empirical_kl_degenerate_inputs_nan(self) -> None:
        """r_lo <= 0 (unvalidated dict records) and zero weight must yield
        NaN from the empirical estimator, never 0.0.
        """
        from tidal.inference._importance import _marginal_dkl_empirical

        col = np.linspace(-1, 1, 100)
        prior = np.linspace(-2, 2, 100)
        w = np.full(100, 0.01)
        assert math.isnan(_marginal_dkl_empirical(col, w, prior, 0.0))
        assert math.isnan(_marginal_dkl_empirical(col, w, prior, -1.0))
        assert math.isnan(_marginal_dkl_empirical(col, np.zeros(100), prior, 1e-3))

    def test_degenerate_weights_raise(self) -> None:
        """Non-finite or zero-sum weights must raise an actionable error
        instead of NaN-propagating into every marginal (review M7).
        """
        pytest.importorskip("anesthetic")
        import tidal.inference._importance as imp_mod

        class _BadNS:
            def get_weights(self):
                return np.zeros(100)

            def stats(self, nsamples: int):
                import pandas as pd

                return pd.DataFrame(
                    {
                        "D_KL": np.ones(nsamples),
                        "d_G": np.ones(nsamples),
                        "logZ": np.zeros(nsamples),
                    }
                )

        class _StubResult:
            method = "nested"
            param_names = ["x"]
            metadata = {}

        original = imp_mod.to_anesthetic_samples
        try:
            imp_mod.to_anesthetic_samples = lambda _r: _BadNS()  # type: ignore[assignment]
            with pytest.raises(ValueError, match="degenerate"):
                imp_mod.compute_parameter_importance(_StubResult(), n_bootstrap=5)  # type: ignore[arg-type]
        finally:
            imp_mod.to_anesthetic_samples = original

    def test_malformed_prior_record_degrades_not_aborts(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A schema-drifted radial_angular record (missing keys) must warn
        and fall back for its own parameters, not abort the computation
        (review M1: the abort was swallowed by save() at debug level,
        silently producing runs with no importance.json).
        """
        pytest.importorskip("anesthetic")
        import logging

        from tidal.inference._importance import compute_parameter_importance
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(37)
        n, nlive = 3000, 500
        x = rng.uniform(-1.0, 1.0, n)
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=np.zeros(n),
            log_prior=np.zeros(n),
            param_names=["c1"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={
                "sampler": "test",
                "nlive": nlive,
                # radial_angular record with every required key missing
                # except the discriminator.
                "priors": [{"kind": "radial_angular"}],
            },
        )
        with caplog.at_level(logging.WARNING, logger="tidal.inference"):
            imp = compute_parameter_importance(result, n_bootstrap=10)
        assert math.isfinite(imp.marginal_d_kl["c1"])  # fallback, not crash
        assert "malformed prior record" in caplog.text
        assert imp.consistency["fallback_params"] == ["c1"]
        assert imp.consistency["product_prior"] is False

    def test_importance_json_roundtrip_strict(self, tmp_path: Path) -> None:
        """The on-disk contract: importance.json must carry schema_version
        and the consistency block, parse under strict JSON (no bare NaN),
        and survive from_directory (review L6/L7).
        """
        pytest.importorskip("anesthetic")
        import json

        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(41)
        n, nlive = 3000, 500
        x = rng.uniform(-1.0, 1.0, n)
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=-0.5 * (x / 0.2) ** 2,
            log_prior=np.zeros(n),
            param_names=["x"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={
                "sampler": "test",
                "nlive": nlive,
                "priors": [
                    {
                        "kind": "scalar",
                        "name": "x",
                        "distribution": "uniform",
                        "low": -1.0,
                        "high": 1.0,
                    }
                ],
            },
        )
        out = tmp_path / "run"
        result.save(out)

        def _reject_constant(_s: str) -> float:
            msg = f"non-strict JSON constant: {_s}"
            raise ValueError(msg)

        raw = (out / "importance.json").read_text()
        data = json.loads(raw, parse_constant=_reject_constant)  # strict
        assert data["schema_version"] == 2
        assert "consistency" in data
        assert "n_eff" in data["consistency"]
        loaded = InferenceResult.from_directory(out)
        pi = loaded.metadata["parameter_importance"]
        assert pi["schema_version"] == 2
        assert "consistency" in pi

    def test_plot_importance_floor_aware(self, tmp_path: Path) -> None:
        """The bar chart must visually distinguish floor-dominated bars
        (gray + hatched, floor line drawn) — it is the only artifact of
        `tidal sample --importance` and previously carried zero caveats
        (#433 review gap).
        """
        pytest.importorskip("matplotlib")
        from tidal.inference._importance import ParameterImportanceResult
        from tidal.inference._visualize import plot_importance

        floored = ParameterImportanceResult(
            param_names=["a", "b"],
            d_kl=1.0,
            d_kl_err=0.1,
            d_g=1.0,
            d_g_err=0.1,
            marginal_d_kl={"a": 0.8, "b": 0.2},
            log_evidence=0.0,
            log_evidence_err=0.1,
            consistency={
                "n_eff": 21.0,
                "noise_floor": {"a": 0.93, "b": 0.93},
                "floor_dominated_params": ["a", "b"],
            },
        )
        out = tmp_path / "floored.png"
        plot_importance(floored, out)
        assert out.exists()
        assert out.stat().st_size > 0

        healthy = ParameterImportanceResult(
            param_names=["a"],
            d_kl=1.0,
            d_kl_err=0.1,
            d_g=1.0,
            d_g_err=0.1,
            marginal_d_kl={"a": 0.8},
            log_evidence=0.0,
            log_evidence_err=0.1,
            consistency={
                "n_eff": 5000.0,
                "noise_floor": {"a": 0.004},
                "floor_dominated_params": [],
            },
        )
        out2 = tmp_path / "healthy.png"
        plot_importance(healthy, out2)
        assert out2.exists()
        assert out2.stat().st_size > 0
        # The floor-annotated figure carries strictly more ink (hatching,
        # extra line, two-line title) — a cheap structural discriminator
        # that fails if the annotations are dropped.
        assert out.stat().st_size != out2.stat().st_size

    def test_posdep_probe_marker_and_single_warning(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """When the GH #421 guard makes the stability probe unavailable,
        the likelihood path must record an explicit stability_profile
        marker and warn exactly once per process — not swallow the guard
        at debug level and silently drop the stability columns
        (review H4).
        """
        import logging

        import tidal.inference._likelihood as lk

        original = lk._posdep_probe_warned
        lk._posdep_probe_warned = False
        try:
            with caplog.at_level(logging.WARNING, logger="tidal.inference"):
                meta1 = lk._posdep_probe_unavailable_meta(
                    NotImplementedError("posdep kinetic")
                )
                meta2 = lk._posdep_probe_unavailable_meta(
                    NotImplementedError("posdep kinetic")
                )
        finally:
            lk._posdep_probe_warned = original
        assert meta1 == {"stability_profile": "unavailable-posdep-kinetic"}
        assert meta2 == meta1
        assert caplog.text.count("WITHOUT the stability gate") == 1

    # -- consistency self-checks ---------------------------------------

    def test_superadditivity_warning_trips_on_buggy_style_marginals(
        self,
    ) -> None:
        """format_importance_table surfaces the impossible sum > joint
        signature (the way #420 announced itself in recorded outputs).
        """
        from tidal.inference._importance import (
            ParameterImportanceResult,
            format_importance_table,
        )

        buggy = ParameterImportanceResult(
            param_names=["a", "b"],
            d_kl=4.64,
            d_kl_err=0.05,
            d_g=2.0,
            d_g_err=0.1,
            marginal_d_kl={"a": 3.18, "b": 4.37},
            log_evidence=-10.0,
            log_evidence_err=0.5,
            consistency={
                "sum_marginals": 7.55,
                "joint_d_kl": 4.64,
                "superadditivity_ok": False,
                "product_prior": True,
                "saturated_params": ["b"],
                "note": "",
            },
        )
        table = format_importance_table(buggy)
        assert "WARNING" in table
        assert "exceeds joint D_KL" in table
        assert "resolution ceiling" in table

        healthy = ParameterImportanceResult(
            param_names=["a"],
            d_kl=2.0,
            d_kl_err=0.05,
            d_g=1.0,
            d_g_err=0.1,
            marginal_d_kl={"a": 1.5},
            log_evidence=-10.0,
            log_evidence_err=0.5,
            consistency={
                "sum_marginals": 1.5,
                "joint_d_kl": 2.0,
                "superadditivity_ok": True,
                "product_prior": True,
                "saturated_params": [],
                "note": "",
            },
        )
        assert "WARNING" not in format_importance_table(healthy)

    # -- noise-floor awareness (#433) ----------------------------------

    def test_low_n_eff_marks_floor_dominated(
        self, caplog: pytest.LogCaptureFixture
    ) -> None:
        """A chain whose weights concentrate on ~20 samples has a
        per-parameter noise floor of ~(40-1)/(2*20) ≈ 1 nat — every
        floor-level marginal must be flagged, not ranked (the T9 rescue
        chain failure mode).
        """
        pytest.importorskip("anesthetic")
        import logging

        from tidal.inference._importance import compute_parameter_importance
        from tidal.inference._prior import Prior
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(5)
        n, nlive = 2000, 100
        prior = Prior(name="x", distribution="uniform", low=-1.0, high=1.0)
        x = prior.sample(rng, n)
        # Steep logL ramp → nested weights concentrate on the top handful
        # of samples → tiny Kish n_eff, mimicking a floor-contaminated
        # rescue chain.
        log_l = 50.0 * np.linspace(0.0, 1.0, n) ** 8
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=log_l,
            log_prior=np.zeros(n),
            param_names=["x"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={"sampler": "test", "nlive": nlive, "priors": [prior]},
        )
        with caplog.at_level(logging.WARNING, logger="tidal.inference"):
            imp = compute_parameter_importance(result, n_bootstrap=10)

        cons = imp.consistency
        assert cons["n_eff"] < 200, "fixture failed to concentrate weights"
        assert "noise_floor" in cons
        assert cons["noise_floor"]["x"] == pytest.approx(
            39.0 / (2 * cons["n_eff"]), rel=1e-9
        )
        # The uniform-prior marginal of a prior-distributed column is pure
        # estimator bias here, so it must be flagged floor-dominated.
        assert cons["floor_dominated_params"] == ["x"]
        assert "noise floor" in caplog.text

        from tidal.inference._importance import format_importance_table

        table = format_importance_table(imp)
        assert "(<= floor)" in table
        assert "do not rank" in table

    def test_healthy_n_eff_no_floor_warning(self) -> None:
        """Equal-weight, large-N chains must not trip the floor flag."""
        pytest.importorskip("anesthetic")
        from tidal.inference._importance import (
            compute_parameter_importance,
            format_importance_table,
        )
        from tidal.inference._prior import Prior
        from tidal.inference._results import InferenceResult

        rng = np.random.default_rng(6)
        n, nlive = 6000, 1000
        prior = Prior(name="x", distribution="uniform", low=-1.0, high=1.0)
        # Genuinely constrained posterior (narrow normal inside the prior)
        # with flat-ish logL ordering so weights stay broad.
        x = np.clip(rng.normal(0.0, 0.05, n), -0.999, 0.999)
        log_l = -0.5 * (x / 0.05) ** 2
        result = InferenceResult(
            samples=x[:, None],
            log_likelihood=log_l,
            log_prior=np.zeros(n),
            param_names=["x"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={"sampler": "test", "nlive": nlive, "priors": [prior]},
        )
        imp = compute_parameter_importance(result, n_bootstrap=10)
        assert imp.consistency["floor_dominated_params"] == []
        assert "(<= floor)" not in format_importance_table(imp)


class TestEffectivePriorSupport:
    """GH #451/#425: one source of truth for the range a prior samples.

    Consumers must never re-derive the support from the recorded
    ``low``/``high`` — ``arctan_uniform`` ignores them, and every
    independent re-derivation of that fact so far has been wrong (the
    #420 estimator read them as a histogram range; the corner plot read
    them as degrees).
    """

    def test_arctan_support_is_the_sampled_range_not_the_bounds(self) -> None:
        from tidal.inference._prior import _ARCTAN_EPS, effective_support

        lo, hi = effective_support("arctan_uniform", -89.0, 89.0)
        expected = math.tan(math.pi / 2 - _ARCTAN_EPS)
        assert lo == pytest.approx(-expected)
        assert hi == pytest.approx(expected)
        # The two wrong answers this replaced.
        assert hi == pytest.approx(19.98, abs=0.01)
        assert hi != pytest.approx(89.0)
        assert hi != pytest.approx(math.tan(math.radians(89.0)), abs=1.0)

    def test_arctan_support_matches_what_sample_actually_draws(self) -> None:
        """The contract is empirical: draws must lie inside the support."""
        from tidal.inference._prior import Prior

        with pytest.warns(UserWarning, match="ignores its bounds"):
            prior = Prior(name="a", distribution="arctan_uniform", low=-89.0, high=89.0)
        lo, hi = prior.effective_support
        draws = prior.sample(np.random.default_rng(0), 20000)
        assert draws.min() >= lo
        assert draws.max() <= hi
        # And the support is tight, not a loose envelope.
        assert draws.max() > 0.9 * hi

    def test_bounded_kinds_report_their_bounds(self) -> None:
        from tidal.inference._prior import effective_support

        assert effective_support("uniform", -2.0, 3.0) == (-2.0, 3.0)
        assert effective_support("log_uniform", 1e-3, 1e3) == (1e-3, 1e3)

    def test_normal_support_is_infinite_not_silently_truncated(self) -> None:
        from tidal.inference._prior import effective_support

        lo, hi = effective_support("normal", 0.0, 1.0)
        assert lo == -math.inf
        assert hi == math.inf

    def test_unknown_kind_raises(self) -> None:
        from tidal.inference._prior import effective_support

        with pytest.raises(ValueError, match="Unknown distribution"):
            effective_support("cauchy", 0.0, 1.0)

    def test_every_valid_kind_has_a_support(self) -> None:
        """Adding a kind without teaching effective_support fails here."""
        from tidal.inference._prior import VALID_DISTRIBUTIONS, effective_support

        for kind in VALID_DISTRIBUTIONS:
            lo, hi = effective_support(kind, 0.5, 2.0)
            assert lo < hi


class TestPriorMapUsesEffectiveSupport:
    """The corner-plot prior map must report sampled ranges (GH #451)."""

    @staticmethod
    def _result(priors: list[dict[str, object]]) -> object:
        from tidal.inference._results import InferenceResult

        n = 32
        return InferenceResult(
            samples=np.zeros((n, 1)),
            log_likelihood=np.zeros(n),
            log_prior=np.zeros(n),
            param_names=["a"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata={"priors": priors},
        )

    def test_archived_record_without_effective_keys_is_corrected_on_read(
        self,
    ) -> None:
        """Old chains record only the unused bounds; do not trust them."""
        from tidal.inference._visualize import (
            _extract_prior_map,  # pyright: ignore[reportPrivateUsage]
        )

        result = self._result(
            [
                {
                    "kind": "scalar",
                    "name": "a",
                    "distribution": "arctan_uniform",
                    "low": -89.0,
                    "high": 89.0,
                }
            ]
        )
        dist, _lo, hi = _extract_prior_map(result)["a"]
        assert dist == "arctan_uniform"
        assert hi == pytest.approx(19.98, abs=0.01)
        # Neither the recorded bound nor the degrees misreading.
        assert hi != pytest.approx(89.0)
        assert hi < 30.0

    def test_recorded_effective_keys_are_used_verbatim(self) -> None:
        from tidal.inference._visualize import (
            _extract_prior_map,  # pyright: ignore[reportPrivateUsage]
        )

        result = self._result(
            [
                {
                    "kind": "scalar",
                    "name": "a",
                    "distribution": "arctan_uniform",
                    "low": -89.0,
                    "high": 89.0,
                    "effective_low": -19.98,
                    "effective_high": 19.98,
                }
            ]
        )
        _dist, lo, hi = _extract_prior_map(result)["a"]
        assert (lo, hi) == (-19.98, 19.98)

    def test_unbounded_support_is_omitted_not_truncated(self) -> None:
        from tidal.inference._visualize import (
            _extract_prior_map,  # pyright: ignore[reportPrivateUsage]
        )

        result = self._result(
            [
                {
                    "kind": "scalar",
                    "name": "a",
                    "distribution": "normal",
                    "low": 0.0,
                    "high": 1.0,
                }
            ]
        )
        assert "a" not in _extract_prior_map(result)

    def test_full_prior_axis_limits_use_the_sampled_range(self) -> None:
        """The GH #451 defect itself: panels drawn to +-57.3, not +-19.98."""
        matplotlib = pytest.importorskip("matplotlib")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from tidal.inference._visualize import (
            _set_full_prior_axis_limits,  # pyright: ignore[reportPrivateUsage]
        )

        fig, ax = plt.subplots()
        axes = np.array([[ax]], dtype=object)
        _set_full_prior_axis_limits(
            axes, ["a"], {"a": ("arctan_uniform", -19.98, 19.98)}
        )
        _lo, hi = ax.get_xlim()
        plt.close(fig)
        assert hi == pytest.approx(19.98, abs=0.01)
        assert hi < 30.0, "regressed to the tan(degrees) misreading"

    def test_arctan_prior_overlay_is_drawn(self) -> None:
        """#309's null check was unavailable for arctan params."""
        matplotlib = pytest.importorskip("matplotlib")
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from tidal.inference._visualize import (
            _overlay_priors,  # pyright: ignore[reportPrivateUsage]
        )

        result = self._result(
            [
                {
                    "kind": "scalar",
                    "name": "a",
                    "distribution": "arctan_uniform",
                    "low": -89.0,
                    "high": 89.0,
                }
            ]
        )
        fig, ax = plt.subplots()
        axes = np.array([[ax]], dtype=object)
        _overlay_priors(axes, result)
        lines = [ln for ln in ax.get_lines() if ln.get_label() == "prior"]
        assert len(lines) == 1, "no prior density drawn for arctan_uniform"
        xs = np.asarray(lines[0].get_xdata(), dtype=float)
        ys = np.asarray(lines[0].get_ydata(), dtype=float)
        # Cauchy shape: peaked at 0, decaying by 1/(1+x^2).
        peak = float(np.max(ys))
        assert float(ys[int(np.argmin(np.abs(xs)))]) == pytest.approx(peak, rel=1e-6)
        assert float(ys[0]) < 0.01 * peak
        plt.close(fig)


class TestFloorDominatedAccessor:
    """GH #433: every ranking consumer asks the same question one way."""

    def test_reads_the_consistency_block(self) -> None:
        from tidal.inference._importance import floor_dominated_params

        block = {
            "marginal_d_kl": {"a": 0.9, "b": 0.1},
            "consistency": {"floor_dominated_params": ["a"]},
        }
        assert floor_dominated_params(block) == frozenset({"a"})

    def test_pre_consistency_block_returns_empty(self) -> None:
        """A pre-v0.48.8 block has no floors computed — not 'no floors'."""
        from tidal.inference._importance import floor_dominated_params

        assert floor_dominated_params({"marginal_d_kl": {"a": 0.9}}) == frozenset()

    def test_malformed_consistency_does_not_raise(self) -> None:
        from tidal.inference._importance import floor_dominated_params

        assert floor_dominated_params({"consistency": "nonsense"}) == frozenset()
        assert (
            floor_dominated_params({"consistency": {"floor_dominated_params": None}})
            == frozenset()
        )

    def test_figure_ranking_drops_floor_marginals_but_keeps_cross(self) -> None:
        """The kl_carrier_corner selection must not rank on estimator noise."""
        import importlib.util
        import sys

        path = (
            Path(__file__).resolve().parents[1]
            / "scripts"
            / "figures"
            / "kl_carrier_corner.py"
        )
        spec = importlib.util.spec_from_file_location("_kl_carrier_corner", path)
        assert spec is not None
        assert spec.loader is not None
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        try:
            spec.loader.exec_module(mod)
        finally:
            sys.modules.pop(spec.name, None)

        importance = {
            "amp": {
                "marginal_d_kl": {"noise": 0.90, "real": 0.05},
                "consistency": {"floor_dominated_params": ["noise", "real"]},
            },
            "sup": {
                "marginal_d_kl": {"noise": 0.88, "real": 0.04},
                "consistency": {"floor_dominated_params": ["noise", "real"]},
            },
            # Prior-free, unaffected by #420/#433 — this is the real signal.
            "cross_amp_sup_kl": {"noise": 0.01, "real": 3.2},
        }
        assert mod._pick_top_k(importance, 1) == ["real"]


class TestNestedSaveWithoutPriorsWarns:
    """GH #434: a nested chain saved with no priors block must say so."""

    @staticmethod
    def _result(metadata: dict[str, object]) -> object:
        from tidal.inference._results import InferenceResult

        n = 40
        rng = np.random.default_rng(3)
        return InferenceResult(
            samples=rng.uniform(-1, 1, (n, 1)),
            log_likelihood=np.zeros(n),
            log_prior=np.zeros(n),
            param_names=["a"],
            method="nested",
            log_evidence=0.0,
            log_evidence_err=0.1,
            weights=np.ones(n) / n,
            metadata=metadata,
        )

    def test_missing_priors_warns_and_names_the_consequence(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        result = self._result({"sampler": "test", "nlive": 10})
        with caplog.at_level(logging.WARNING, logger="tidal.inference"):
            result.save(tmp_path)  # pyright: ignore[reportAttributeAccessIssue]
        assert "no priors metadata" in caplog.text
        assert "#434" in caplog.text

    def test_recorded_priors_do_not_warn(
        self,
        tmp_path: Path,
        caplog: pytest.LogCaptureFixture,
    ) -> None:
        from tidal.inference._prior import Prior

        prior = Prior(name="a", distribution="uniform", low=-1.0, high=1.0)
        result = self._result(
            {"sampler": "test", "nlive": 10, "priors": [prior]},
        )
        with caplog.at_level(logging.WARNING, logger="tidal.inference"):
            result.save(tmp_path)  # pyright: ignore[reportAttributeAccessIssue]
        assert "no priors metadata" not in caplog.text


# ===================================================================
# Analyze CLI (inference path)
# ===================================================================


class TestAnalyzeInference:
    """Test tidal analyze --inference --importance."""

    def test_parse_analyze_inference_flags(self) -> None:
        from tidal.cli import (
            _build_parser as build_parser,  # pyright: ignore[reportPrivateUsage]
        )

        parser = build_parser()
        args = parser.parse_args(
            [
                "analyze",
                "/tmp/test",
                "--inference",
                "--importance",
                "--n-bootstrap",
                "50",
            ],
        )
        assert args.inference is True
        assert args.importance is True
        assert args.n_bootstrap_importance == 50
