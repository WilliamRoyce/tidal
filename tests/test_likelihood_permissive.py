"""Tests for the v3 permissive likelihood mode.

The v3 architecture (post-2026-05-08 supervisor pivot, see
``docs/V3_ARCHITECTURE.md``) defaults to ``permissive=True``: the pre-flight
tachyonic probe runs as a metadata measurement only — its verdict no longer
gates samples.  ``--gated`` (i.e. ``permissive=False``) preserves the v2 /
canonical-probe hard-rejection behaviour for reproducibility.

These tests verify the LikelihoodConfig surface and parser plumbing; the
end-to-end behaviour with a real simulation lives in the integration suite.
"""

from __future__ import annotations

import pytest

from tidal.inference._likelihood import (
    LikelihoodConfig,
    parse_likelihood,
)


class TestLikelihoodConfigDefaults:
    def test_permissive_defaults_true(self) -> None:
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert lc.permissive is True

    def test_soft_floor_noise_defaults_to_one(self) -> None:
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert lc.soft_floor_noise_sigma == 1.0

    def test_explicit_gated_construction(self) -> None:
        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            permissive=False,
        )
        assert lc.permissive is False

    def test_explicit_zero_noise_construction(self) -> None:
        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            soft_floor_noise_sigma=0.0,
        )
        assert lc.soft_floor_noise_sigma == 0.0


class TestParseLikelihoodPlumbing:
    def test_default_kwargs_propagate(self) -> None:
        """parse_likelihood without kwargs gives v3 defaults (permissive,
        sigma=1.0).
        """
        lc = parse_likelihood("P_max:maximize")
        assert lc.permissive is True
        assert lc.soft_floor_noise_sigma == 1.0

    def test_gated_via_parse_likelihood(self) -> None:
        lc = parse_likelihood("P_max:maximize", permissive=False)
        assert lc.permissive is False

    def test_soft_floor_noise_via_parse_likelihood(self) -> None:
        lc = parse_likelihood("P_max:maximize", soft_floor_noise_sigma=0.0)
        assert lc.soft_floor_noise_sigma == 0.0

    def test_kwargs_propagate_through_minimize(self) -> None:
        lc = parse_likelihood(
            "P_max:minimize",
            permissive=False,
            soft_floor_noise_sigma=2.5,
        )
        assert lc.likelihood_type == "minimize"
        assert lc.permissive is False
        assert lc.soft_floor_noise_sigma == 2.5

    def test_kwargs_propagate_through_gaussian(self) -> None:
        lc = parse_likelihood(
            "P_max:gaussian:0.5:0.1",
            permissive=False,
            soft_floor_noise_sigma=3.0,
        )
        assert lc.likelihood_type == "gaussian"
        assert lc.target == 0.5
        assert lc.sigma == 0.1
        assert lc.permissive is False
        assert lc.soft_floor_noise_sigma == 3.0

    def test_kwargs_propagate_through_threshold(self) -> None:
        lc = parse_likelihood(
            "P_max:threshold:0.01",
            permissive=False,
        )
        assert lc.likelihood_type == "threshold"
        assert lc.min_value == 0.01
        assert lc.permissive is False

    def test_kwargs_propagate_through_extremize(self) -> None:
        lc = parse_likelihood(
            "P_max:extremize",
            baseline_formula="0.5",
            permissive=False,
            soft_floor_noise_sigma=0.5,
        )
        assert lc.likelihood_type == "extremize"
        assert lc.permissive is False
        assert lc.soft_floor_noise_sigma == 0.5


class TestRunStatusTags:
    """Document the distinct ``run_status`` tags surfaced by v3.

    These are smoke checks — the actual tagging behaviour exercises the
    full simulation pipeline and lives in the integration tests.  Here we
    just enumerate the tags so the v3 architecture's commitment to
    distinct post-chain filtering is documented in code.
    """

    @pytest.mark.parametrize(
        "tag",
        [
            "success",
            "tachyonic_gated",
            "simulation_diverged",
            "simulation_failed",
            "metric_nan",
            "metric_missing",
            "logl_minus_inf",
            "below_noise_floor",
            "exception",
        ],
    )
    def test_run_status_tag_is_a_string(self, tag: str) -> None:
        """Each tag is a documented post-chain filter target.  See
        ``docs/V3_ARCHITECTURE.md`` Soft-penalty table.
        """
        assert isinstance(tag, str)
        assert tag  # non-empty
