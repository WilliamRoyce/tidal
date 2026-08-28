"""Tests confirming the v3 architecture removed the Hwang-Noh perturbativity gate.

In v2 a sample with ``P_max > 0.5`` (and ``maximize`` / ``minimize`` likelihood
type) returned ``logL = -inf`` with ``run_status="non_perturbative"`` so it
never entered the live-points pool.  Per the v3 supervisor pivot (2026-05-08,
docs/V3_ARCHITECTURE.md), this gate is removed entirely: ``P_max`` is recorded
faithfully at all values; downstream analysis interprets whether the
linearized result is physically meaningful.

The internal ``compute_log_likelihood`` keeps a NaN/Inf safeguard returning
``-inf`` for genuinely undefined inputs (preserves external test contract);
the upstream ``_evaluate_likelihood`` routes that to the soft floor with
``run_status="metric_nan"`` (covered in ``test_likelihood_soft_floor_noise.py``
and the manual smoke at the end of this module).
"""

from __future__ import annotations

import math

from tidal.inference._likelihood import (
    LikelihoodConfig,
    compute_log_likelihood,
)


class TestHwangNohRemoved:
    def test_p_max_above_half_returns_finite_logl_maximize(self) -> None:
        """P_max=0.6 was -inf in v2; v3 returns the raw metric (= log A under
        baseline-formula, or just metric_value without one).
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        logl = compute_log_likelihood(0.6, lc)
        assert math.isfinite(logl)
        assert logl == 0.6  # maximize without baseline returns raw metric

    def test_p_max_above_half_returns_finite_logl_minimize(self) -> None:
        lc = LikelihoodConfig(metric="P_max", likelihood_type="minimize")
        logl = compute_log_likelihood(0.6, lc)
        assert math.isfinite(logl)
        assert logl == -0.6  # minimize negates the raw metric

    def test_p_max_above_unity_returns_finite_logl(self) -> None:
        """v3 specifically removed the upper cap.  P_max > 1 is faithful
        linearized-PDE output even though it violates probability
        conservation (interaction energies become negative to compensate).
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        for v in (1.5, 2.0, 5.0, 10.0):
            logl = compute_log_likelihood(v, lc)
            assert math.isfinite(logl), f"P_max={v} returned {logl}"
            assert logl == v

    def test_p_max_above_two_no_longer_capped(self) -> None:
        """v2 had a hard cap at P_max > 2.0 (the "non-perturbative regime"
        artifact protection).  v3 removed it entirely.
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        # In v2, compute_log_likelihood(3.0, lc) would have been -inf.  In v3:
        assert compute_log_likelihood(3.0, lc) == 3.0
        # Same at 100, 1e6, 1e20 — no cap.
        assert compute_log_likelihood(100.0, lc) == 100.0
        assert compute_log_likelihood(1e6, lc) == 1e6

    def test_baseline_formula_log_amplification_unchanged(self) -> None:
        """When baseline_formula is supplied, maximize returns log(P/baseline)
        — verify this still works for P > 0.5 (which v2 would have rejected
        before this branch fired).
        """
        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            baseline_formula="0.001",  # constant baseline ≈ P_GR for tiny B0
        )
        # Use eval_params={} since formula has no parameters
        logl = compute_log_likelihood(0.6, lc, eval_params={})
        # log(0.6 / 0.001) = log(600) ≈ 6.396
        assert math.isfinite(logl)
        assert abs(logl - math.log(600.0)) < 1e-9

    def test_nan_metric_still_returns_minus_inf_internally(self) -> None:
        """``compute_log_likelihood`` keeps the NaN safeguard for external
        callers (preserves the test_inference.py contract).  The upstream
        ``_evaluate_likelihood`` is responsible for routing this to the soft
        floor with ``run_status="metric_nan"`` (smoke-tested separately).
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert compute_log_likelihood(float("nan"), lc) == -math.inf
        assert compute_log_likelihood(float("inf"), lc) == -math.inf

    def test_no_gate_can_reintroduce_the_value_cap(self) -> None:
        """There is no configuration under which P_max > 0.5 is penalized.

        This used to compare permissive against ``--gated``: that flag
        controlled the upstream tachyonic gate, never the value-level
        Hwang-Noh check, and neither mode re-enabled the cap.  With the
        gate removed entirely (v0.49.6) there is only one mode left, and
        the property to pin is simply that a large P_max is scored
        faithfully.
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert compute_log_likelihood(0.7, lc) == 0.7
        assert compute_log_likelihood(5.0, lc) == 5.0
