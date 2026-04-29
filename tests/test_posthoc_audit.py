"""Tests for the post-hoc t-independence audit (#323)."""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING

import numpy as np
import pytest

from tidal.measurement._posthoc_audit import (
    AuditSummary,
    SampleAudit,
    _build_summary,
    _classify_sample,
    _eval_baseline_at_times,
    _write_audit_csv,
    select_audit_samples,
)

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------------
# _classify_sample — the policy gate
# ---------------------------------------------------------------------


class TestClassifySample:
    """Pin the per-sample classification logic.

    These tests are the authoritative specification of the
    Hwang-Noh + flag-for-review policy described in
    ``tidal/measurement/_posthoc_audit.py``'s module docstring.
    """

    def _classify(
        self,
        *,
        p1: float,
        p2: float,
        a1: float = 1.0,
        a2: float = 1.0,
        borderline: bool = False,
        p_max_lim: float = 0.5,
    ) -> tuple[str, str]:
        return _classify_sample(
            p1,
            p2,
            a1,
            a2,
            perturbativity_p_max=p_max_lim,
            a_log_ratio_threshold=math.log(3.0),
            p_growth_ratio_threshold=5.0,
            borderline_flag=borderline,
        )

    def test_pass_for_small_bounded_p(self) -> None:
        status, _ = self._classify(p1=0.01, p2=0.02, a1=1.5, a2=1.5)
        assert status == "pass"

    def test_hard_reject_for_p_above_hwang_noh(self) -> None:
        status, reason = self._classify(p1=0.001, p2=0.6)
        assert status == "hard_reject"
        assert "Hwang-Noh" in reason

    def test_hard_reject_dominates_borderline(self) -> None:
        """Hwang-Noh hard-reject takes priority over flag_review."""
        status, _ = self._classify(p1=0.001, p2=0.7, a1=1.0, a2=10.0, borderline=True)
        assert status == "hard_reject"

    def test_flag_review_on_a_ratio(self) -> None:
        """A_med ratio > 3 → flag_review (exponential-growth proxy)."""
        status, reason = self._classify(p1=0.01, p2=0.02, a1=1.0, a2=10.0)
        assert status == "flag_review"
        assert "log(A_med_2/A_med_1)" in reason

    def test_flag_review_on_p_growth_ratio(self) -> None:
        """P(2t)/P(t) > 5 catches slow growth even when both halves <0.5."""
        status, reason = self._classify(p1=0.001, p2=0.05, a1=1.0, a2=1.0)
        assert status == "flag_review"
        assert "P_max_2/P_max_1" in reason

    def test_flag_review_on_borderline_alone(self) -> None:
        """Probe borderline propagates without other signals firing."""
        status, reason = self._classify(p1=0.01, p2=0.02, borderline=True)
        assert status == "flag_review"
        assert "probe-borderline" in reason

    def test_reasons_are_concatenated_when_multiple(self) -> None:
        status, reason = self._classify(
            p1=0.001,
            p2=0.05,
            a1=1.0,
            a2=10.0,
            borderline=True,
        )
        assert status == "flag_review"
        assert reason.count(";") >= 2

    def test_perturbativity_threshold_is_per_profile(self) -> None:
        """A theory-tighter cutoff (e.g. 0.3 for parity-odd) flips the gate."""
        status, _ = self._classify(p1=0.001, p2=0.4, p_max_lim=0.3)
        assert status == "hard_reject"
        # Same data, default 0.5 cutoff → only flag_review (P-ratio > 5).
        status2, _ = self._classify(p1=0.001, p2=0.4, p_max_lim=0.5)
        assert status2 == "flag_review"


# ---------------------------------------------------------------------
# _eval_baseline_at_times — closed-form Gertsenshtein
# ---------------------------------------------------------------------


class TestEvalBaselineAtTimes:
    def test_evaluates_sin2_per_snapshot(self) -> None:
        params = {"kappa": 1.0, "B0": 0.01}
        times = np.array([0.0, 5.0, 10.0, 20.0])
        baseline = _eval_baseline_at_times(
            "sin(kappa*B0*t_end/2)**2",
            params,
            times,
        )
        # P_GR(t=0) = 0; P_GR(20) = sin²(0.1) ≈ 0.00997...
        assert baseline[0] == pytest.approx(0.0)
        assert baseline[-1] == pytest.approx(math.sin(0.1) ** 2)

    def test_unknown_param_yields_nan(self) -> None:
        baseline = _eval_baseline_at_times(
            "sin(unknown*t_end)**2",
            {},
            np.array([1.0]),
        )
        assert np.isnan(baseline[0])


# ---------------------------------------------------------------------
# select_audit_samples — strategy priority + bias-free coverage
# ---------------------------------------------------------------------


class _FakeInferenceResult:
    """Minimal stub for select_audit_samples. Mirrors InferenceResult fields."""

    def __init__(
        self,
        n_samples: int,
        n_params: int,
        *,
        log_likelihood: np.ndarray | None = None,
        log_prior: np.ndarray | None = None,
        weights: np.ndarray | None = None,
        borderline: np.ndarray | None = None,
        rng_seed: int = 1,
    ) -> None:
        rng = np.random.default_rng(rng_seed)
        self.n_samples = n_samples
        self.n_params = n_params
        self.samples = rng.uniform(-1, 1, size=(n_samples, n_params))
        self.log_likelihood = (
            log_likelihood
            if log_likelihood is not None
            else rng.normal(0, 1, size=n_samples)
        )
        self.log_prior = log_prior if log_prior is not None else np.zeros(n_samples)
        self.weights = weights
        self.metrics: dict[str, np.ndarray] = {}
        if borderline is not None:
            self.metrics["borderline_stability"] = borderline


class TestSelectAuditSamples:
    def test_borderline_takes_priority(self) -> None:
        """Borderline samples are selected even when not in top-N."""
        n = 100
        borderline = np.zeros(n, dtype=int)
        borderline[[5, 17, 42]] = 1  # Three borderline points
        log_lik = np.zeros(n)
        log_lik[80:] = 10.0  # Top 20 samples are at indices 80–99
        result = _FakeInferenceResult(
            n,
            3,
            log_likelihood=log_lik,
            borderline=borderline,
        )

        selected = select_audit_samples(result, n_top=10, n_stratified=0)
        assert 5 in selected
        assert selected[5] == "borderline"
        assert 17 in selected
        assert selected[17] == "borderline"
        assert 42 in selected
        assert selected[42] == "borderline"

    def test_top_n_picks_highest_log_posterior(self) -> None:
        n = 50
        log_lik = np.linspace(0, 10, n)
        result = _FakeInferenceResult(n, 2, log_likelihood=log_lik)
        selected = select_audit_samples(result, n_top=5, n_stratified=0)
        # Top 5 are indices 45..49.
        for idx in range(45, 50):
            assert idx in selected

    def test_map_is_always_included(self) -> None:
        n = 30
        log_lik = np.zeros(n)
        log_lik[12] = 100.0  # Clear MAP at index 12
        result = _FakeInferenceResult(n, 2, log_likelihood=log_lik)
        selected = select_audit_samples(
            result,
            n_top=0,
            n_stratified=0,
            include_map=True,
        )
        assert 12 in selected

    def test_priority_order_borderline_over_top(self) -> None:
        n = 100
        log_lik = np.zeros(n)
        log_lik[42] = 1000.0  # MAP at 42
        borderline = np.zeros(n, dtype=int)
        borderline[42] = 1  # Also borderline
        result = _FakeInferenceResult(
            n,
            2,
            log_likelihood=log_lik,
            borderline=borderline,
        )
        selected = select_audit_samples(result, n_top=10, n_stratified=0)
        # Borderline wins the priority race even though MAP would also pick it.
        assert selected[42] == "borderline"

    def test_stratified_uses_clustering_when_sklearn_present(self) -> None:
        """Stratified strategy at minimum picks one sample per cluster.

        Skips silently if scikit-learn is unavailable.
        """
        pytest.importorskip("sklearn")
        n = 100
        result = _FakeInferenceResult(n, 3)
        selected = select_audit_samples(
            result,
            n_top=0,
            n_stratified=20,
            n_clusters=4,
            include_map=False,
        )
        strats = set(selected.values())
        assert "stratified" in strats


# ---------------------------------------------------------------------
# _write_audit_csv — schema contract
# ---------------------------------------------------------------------


class TestWriteAuditCsv:
    def test_csv_header_pinned(self, tmp_path: Path) -> None:
        audits = [
            SampleAudit(
                index=0,
                parameters={"a": 0.1, "b": 0.2},
                p_max_first=0.01,
                p_max_second=0.02,
                p_max_total=0.02,
                a_med_first=1.0,
                a_med_second=1.05,
                a_max_first=1.0,
                a_max_second=1.05,
                a_med_log_ratio=0.05,
                gamma_eff_posthoc=0.005,
                status="pass",
                reason="ok",
                selection_strategy="top",
            ),
        ]
        out = tmp_path / "audit.csv"
        _write_audit_csv(audits, out, ["a", "b"])

        text = out.read_text(encoding="utf-8").splitlines()
        header = text[0].split(",")
        assert header[0] == "sample_index"
        assert "a" in header
        assert "b" in header
        for required in (
            "p_max_total",
            "a_med_log_ratio",
            "gamma_eff_posthoc",
            "audit_status",
            "audit_reason",
            "selection_strategy",
        ):
            assert required in header


# ---------------------------------------------------------------------
# _build_summary — aggregate report
# ---------------------------------------------------------------------


class TestBuildSummary:
    def test_summary_counts_and_contamination(self) -> None:
        audits = [
            _make_audit(0, "pass", 0.001),
            _make_audit(1, "flag_review", 0.05),
            _make_audit(2, "hard_reject", 0.30),
            _make_audit(3, "pass", 0.002),
            _make_audit(4, "sim_failed", float("nan")),
        ]

        class _Stub:
            n_samples = 100

        summary = _build_summary(
            _Stub(),  # type: ignore[arg-type]
            audits,
            profile_name="unit-ic-all-k-0.3",
            t_base=10.0,
            perturbativity_p_max=0.5,
            selection_strategies=("top", "borderline"),
        )
        assert summary.n_samples_audited == 5
        assert summary.n_pass == 2
        assert summary.n_flag_review == 1
        assert summary.n_hard_reject == 1
        assert summary.n_sim_failed == 1
        assert summary.contamination_rate == pytest.approx(2 / 5)
        assert summary.worst_sample_index == 2
        assert summary.worst_gamma_eff_posthoc == pytest.approx(0.30)
        assert summary.profile_name == "unit-ic-all-k-0.3"

    def test_summary_handles_no_finite_gammas(self) -> None:
        audits = [_make_audit(0, "sim_failed", float("nan"))]

        class _Stub:
            n_samples = 1

        summary = _build_summary(
            _Stub(),  # type: ignore[arg-type]
            audits,
            profile_name="unit-ic-all-k-0.3",
            t_base=10.0,
            perturbativity_p_max=0.5,
            selection_strategies=("top",),
        )
        assert summary.worst_sample_index is None

    def test_summary_serializes_to_json(self) -> None:
        """audit_summary.json round-trip — must be JSON-clean."""
        import dataclasses

        summary = AuditSummary(
            n_samples_total=10,
            n_samples_audited=5,
            n_pass=3,
            n_flag_review=1,
            n_hard_reject=1,
            n_sim_failed=0,
            contamination_rate=0.4,
            worst_sample_index=7,
            worst_gamma_eff_posthoc=0.12,
            profile_name="unit-ic-all-k-0.3",
            t_base=10.0,
            perturbativity_p_max=0.5,
            selection_strategies_used=("top", "borderline"),
        )
        # Should not raise on json.dumps.
        s = json.dumps(dataclasses.asdict(summary))
        loaded = json.loads(s)
        assert loaded["n_samples_total"] == 10
        assert loaded["selection_strategies_used"] == ["top", "borderline"]


# Helper for building synthetic audits.
def _make_audit(index: int, status: str, gamma: float) -> SampleAudit:
    return SampleAudit(
        index=index,
        parameters={},
        p_max_first=0.0,
        p_max_second=0.0,
        p_max_total=0.0,
        a_med_first=1.0,
        a_med_second=1.0,
        a_max_first=1.0,
        a_max_second=1.0,
        a_med_log_ratio=0.0,
        gamma_eff_posthoc=gamma,
        status=status,  # type: ignore[arg-type]
        reason="",
        selection_strategy="top",
    )
