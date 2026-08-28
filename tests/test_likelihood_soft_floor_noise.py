"""Tests for the v3 soft-penalty floor with optional Gaussian exploration noise.

Per the v3 architecture (post-2026-05-08 supervisor pivot, see
``docs/V3_ARCHITECTURE.md``), numerical-failure samples (sim divergence /
NaN / exception) return ``logL = SOFT_FLOOR_LOGL + Normal(0, sigma_explore)``
instead of the v2 ``-inf``.  The Gaussian noise gives the sampler a finite
gradient in the failure region so it doesn't see a flat plateau.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tidal.inference._likelihood import (
    SOFT_FLOOR_LOGL,
    _floor_rng,
    _soft_floor_logl,
    parse_likelihood,
)


class TestSoftFloorLogL:
    def test_sigma_zero_returns_bare_floor_exactly(self) -> None:
        """``--soft-floor-noise=0`` ablation path: deterministic -100."""
        v = _soft_floor_logl(0.0)
        assert v == SOFT_FLOOR_LOGL
        assert SOFT_FLOOR_LOGL == -100.0

    def test_sigma_negative_returns_bare_floor(self) -> None:
        """Defensive: negative sigma should also disable noise."""
        assert _soft_floor_logl(-1.0) == SOFT_FLOOR_LOGL

    def test_default_sigma_produces_distribution_around_floor(self) -> None:
        """Default sigma=1.0 noise samples a Normal around -100."""
        rng = np.random.default_rng(seed=42)
        samples = [_soft_floor_logl(1.0, rng=rng) for _ in range(2000)]
        mean = sum(samples) / len(samples)
        assert abs(mean - SOFT_FLOOR_LOGL) < 0.1, (
            f"sample mean {mean} drifted too far from {SOFT_FLOOR_LOGL}"
        )
        # Empirical sigma should be ~1.0
        std = float(np.std(np.array(samples)))
        assert 0.85 < std < 1.15, f"sample std {std} not ~1.0"

    def test_repeated_calls_with_default_rng_produce_different_values(self) -> None:
        """With no generator supplied, distinct calls give distinct values.

        This is the bare-helper contract only.  Production no longer relies
        on it: ``_evaluate_likelihood`` passes a ``theta``-derived generator
        so the floor is reproducible per parameter point (issue #408).  The
        gradient the sampler needs comes from the noise varying *across*
        samples, not from it varying between repeated evaluations of the
        same one.
        """
        floors = {_soft_floor_logl(1.0) for _ in range(50)}
        # Allow up to one duplicate by chance; >=49/50 unique is the bar
        assert len(floors) >= 48


class TestFloorRng:
    """The theta-derived generator behind reproducible soft-floor noise."""

    def test_same_seed_and_theta_reproduce(self) -> None:
        theta = [1.0, -2.5, 3.25]
        a = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(42, theta))
        b = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(42, theta))
        assert a == b

    def test_different_theta_gives_different_noise(self) -> None:
        a = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(42, [1.0, 2.0]))
        b = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(42, [1.0, 2.0001]))
        assert a != b

    def test_different_seed_gives_different_noise(self) -> None:
        theta = [1.0, 2.0]
        a = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(42, theta))
        b = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(43, theta))
        assert a != b

    def test_independent_of_evaluation_order(self) -> None:
        """The guard against regressing to counter-based seeding.

        A counter restarts at zero in every multiprocessing worker, so the
        k-th task on each worker would draw identical noise.  Deriving from
        theta means the value cannot depend on how many evaluations came
        before it.
        """
        theta = [0.5, 0.25]
        first = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(7, theta))
        for other in ([9.0], [1.0, 1.0], [3.3, 4.4, 5.5]):
            _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(7, other))
        later = _soft_floor_logl(1.0, SOFT_FLOOR_LOGL, _floor_rng(7, theta))
        assert first == later

    def test_config_carries_the_seed(self) -> None:
        assert parse_likelihood("P_max:maximize").noise_seed == 0
        assert parse_likelihood("P_max:maximize", noise_seed=7).noise_seed == 7
        assert parse_likelihood("P_max:minimize", noise_seed=7).noise_seed == 7

    def test_floor_well_below_typical_real_loglikelihood(self) -> None:
        """Real samples (typical logL > -50) should always out-compete the
        soft floor in PolyChord's posterior weight.
        """
        # exp(-100) << exp(-50): floor samples are effectively rejected
        # but visible to the sampler as a navigable surface.
        assert SOFT_FLOOR_LOGL < -50.0
        # And exp(SOFT_FLOOR_LOGL) is small enough to be effectively zero
        # in posterior weighting.
        assert math.exp(SOFT_FLOOR_LOGL) < 1e-40

    def test_explicit_rng_reproducibility(self) -> None:
        """Passing the same RNG seed twice should give identical floors."""
        rng1 = np.random.default_rng(seed=123)
        rng2 = np.random.default_rng(seed=123)
        v1 = _soft_floor_logl(1.0, rng=rng1)
        v2 = _soft_floor_logl(1.0, rng=rng2)
        assert v1 == v2

    @pytest.mark.parametrize("sigma", [0.5, 1.0, 2.0, 5.0])
    def test_noise_scales_with_sigma(self, sigma: float) -> None:
        """The empirical std of the noise should match the sigma parameter."""
        rng = np.random.default_rng(seed=999)
        samples = [_soft_floor_logl(sigma, rng=rng) for _ in range(2000)]
        std = float(np.std(np.array(samples)))
        # 10% tolerance on a 2000-sample empirical std
        assert abs(std - sigma) / sigma < 0.1, (
            f"empirical std {std} doesn't match sigma {sigma}"
        )


class TestEvaluateLikelihoodDeterminism:
    """End-to-end: the floor is reproducible through ``_evaluate_likelihood``.

    The helper-level tests above would still pass if a call site were wired
    up wrongly (say, passing the wrong variable into ``_floor_rng``), so
    this drives the real function.  A nonexistent spec path sends it down
    the ``run_status="exception"`` branch, which is one of the four floor
    returns.  See issue #408.
    """

    @staticmethod
    def _evaluate(theta: list[float], *, noise_seed: int = 42) -> tuple[float, str]:
        from argparse import Namespace
        from pathlib import Path

        from tidal.inference._likelihood import (
            LikelihoodConfig,
            _evaluate_likelihood,
        )

        config = LikelihoodConfig(
            metric="P_max",
            likelihood_type="maximize",
            soft_floor_noise_sigma=1.0,
            noise_seed=noise_seed,
        )
        logl, meta = _evaluate_likelihood(
            theta=theta,
            base_args=Namespace(),
            spec_path=Path("/nonexistent/spec.json"),
            param_names=["a", "b"],
            measurements={"peak_conversion"},
            source=None,
            target=None,
            threshold=0.0,
            likelihood_config=config,
            temp_dir=None,
            keep_sims=False,
            call_index=0,
        )
        return logl, str(meta["run_status"])

    def test_same_theta_gives_identical_logl(self) -> None:
        first, status = self._evaluate([1.0, 2.0])
        second, _ = self._evaluate([1.0, 2.0])
        assert status == "exception"
        assert first == second

    def test_different_theta_gives_different_logl(self) -> None:
        first, _ = self._evaluate([1.0, 2.0])
        other, _ = self._evaluate([1.0, 2.5])
        assert first != other

    def test_different_noise_seed_gives_different_logl(self) -> None:
        first, _ = self._evaluate([1.0, 2.0], noise_seed=42)
        other, _ = self._evaluate([1.0, 2.0], noise_seed=43)
        assert first != other

    def test_floor_value_is_near_the_configured_floor(self) -> None:
        logl, _ = self._evaluate([1.0, 2.0])
        assert abs(logl - SOFT_FLOOR_LOGL) < 10.0


class TestFinerTagsDoNotMoveTheChain:
    """GH #480: the new status tags must be diagnostics, nothing more.

    ``SimulationDivergedError`` and ``KineticEvaluationError`` used to fall
    into the bare ``except Exception`` and be tagged ``exception``.  They
    now get their own tags — but they take the SAME soft floor, with the
    same ``theta``-derived generator, so the logL is unchanged and no
    recorded chain could shift because of the finer labelling.

    This is the claim the change rests on, so it is asserted rather than
    argued: same ``theta``, three different failure causes, one logL.
    """

    @staticmethod
    def _evaluate_raising(
        exc: BaseException,
        theta: list[float],
        *,
        noise_seed: int = 42,
    ) -> tuple[float, str]:
        from argparse import Namespace
        from pathlib import Path

        import tidal.symbolic._spec_cache as spec_cache
        from tidal.inference._likelihood import (
            LikelihoodConfig,
            _evaluate_likelihood,  # pyright: ignore[reportPrivateUsage]
        )

        def _boom(*_a: object, **_k: object) -> object:
            raise exc

        original = spec_cache.load_spec_cached
        spec_cache.load_spec_cached = _boom  # type: ignore[assignment]
        try:
            logl, meta = _evaluate_likelihood(
                theta=theta,
                base_args=Namespace(),
                spec_path=Path("/nonexistent/spec.json"),
                param_names=["a", "b"],
                measurements={"peak_conversion"},
                source=None,
                target=None,
                threshold=0.0,
                likelihood_config=LikelihoodConfig(
                    metric="P_max",
                    likelihood_type="maximize",
                    soft_floor_noise_sigma=1.0,
                    noise_seed=noise_seed,
                ),
                temp_dir=None,
                keep_sims=False,
                call_index=0,
            )
        finally:
            spec_cache.load_spec_cached = original  # type: ignore[assignment]
        return logl, str(meta["run_status"])

    def test_divergence_is_tagged_but_scored_identically(self) -> None:
        from tidal.solver import SimulationDivergedError

        theta = [1.0, 2.0]
        diverged_logl, diverged_status = self._evaluate_raising(
            SimulationDivergedError("fields blew up"), theta
        )
        generic_logl, generic_status = self._evaluate_raising(
            RuntimeError("something else"), theta
        )

        assert diverged_status == "simulation_diverged"
        assert generic_status == "exception"
        # The tag is finer; the number is not different.
        assert diverged_logl == generic_logl

    def test_kinetic_error_is_tagged_but_scored_identically(self) -> None:
        from tidal.solver import KineticEvaluationError

        theta = [1.0, 2.0]
        kinetic_logl, kinetic_status = self._evaluate_raising(
            KineticEvaluationError("xi unbound"), theta
        )
        generic_logl, _ = self._evaluate_raising(RuntimeError("other"), theta)

        assert kinetic_status == "kinetic_error"
        assert kinetic_logl == generic_logl

    def test_the_new_branches_still_land_on_the_floor(self) -> None:
        from tidal.solver import SimulationDivergedError

        logl, _ = self._evaluate_raising(SimulationDivergedError("boom"), [1.0, 2.0])
        assert abs(logl - SOFT_FLOOR_LOGL) < 10.0
