"""Tests for the ``LikelihoodConfig`` surface and ``parse_likelihood`` plumbing.

Was ``test_likelihood_permissive.py``.  The v3 architecture
(post-2026-05-08 supervisor pivot, see ``docs/V3_ARCHITECTURE.md``) made
the pre-flight tachyonic probe a metadata-only measurement, keeping a
``permissive=False`` / ``--gated`` escape hatch that restored v2
hard-rejection.  That escape hatch was removed in v0.50.0: rejection on
tachyonic growth is abandoned policy, because growth cannot be
classified as physics or artifact without theory-level analysis
(PSALTer, GH #360).  Nothing now acts on the probe verdict, so there is
no gate to configure and the tests that covered one are gone.

What remains here is the soft-floor configuration surface, which is
unrelated to gating and was only ever colocated with it.
"""

from __future__ import annotations

from tidal.inference._likelihood import (
    LikelihoodConfig,
    parse_likelihood,
)


class TestLikelihoodConfigDefaults:
    def test_soft_floor_noise_defaults_to_one(self) -> None:
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert lc.soft_floor_noise_sigma == 1.0

    def test_explicit_zero_noise_construction(self) -> None:
        lc = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            soft_floor_noise_sigma=0.0,
        )
        assert lc.soft_floor_noise_sigma == 0.0

    def test_no_gate_to_configure(self) -> None:
        """The probe verdict is not a configurable rejection any more.

        Pins the v0.50.0 removal: a ``permissive`` field reappearing would
        mean the abandoned policy had been made reachable again.
        """
        lc = LikelihoodConfig(metric="P_max", likelihood_type="maximize")
        assert not hasattr(lc, "permissive")


class TestParseLikelihoodPlumbing:
    def test_default_kwargs_propagate(self) -> None:
        lc = parse_likelihood("P_max:maximize")
        assert lc.soft_floor_noise_sigma == 1.0

    def test_soft_floor_noise_via_parse_likelihood(self) -> None:
        lc = parse_likelihood("P_max:maximize", soft_floor_noise_sigma=0.0)
        assert lc.soft_floor_noise_sigma == 0.0

    def test_kwargs_propagate_through_minimize(self) -> None:
        lc = parse_likelihood("P_max:minimize", soft_floor_noise_sigma=2.5)
        assert lc.likelihood_type == "minimize"
        assert lc.soft_floor_noise_sigma == 2.5

    def test_kwargs_propagate_through_gaussian(self) -> None:
        lc = parse_likelihood("P_max:gaussian:0.5:0.1", soft_floor_noise_sigma=3.0)
        assert lc.likelihood_type == "gaussian"
        assert lc.target == 0.5
        assert lc.sigma == 0.1
        assert lc.soft_floor_noise_sigma == 3.0

    def test_kwargs_propagate_through_threshold(self) -> None:
        lc = parse_likelihood("P_max:threshold:0.01")
        assert lc.likelihood_type == "threshold"
        assert lc.min_value == 0.01

    def test_kwargs_propagate_through_extremize(self) -> None:
        lc = parse_likelihood(
            "P_max:extremize",
            baseline_formula="0.5",
            soft_floor_noise_sigma=0.5,
        )
        assert lc.likelihood_type == "extremize"
        assert lc.soft_floor_noise_sigma == 0.5
