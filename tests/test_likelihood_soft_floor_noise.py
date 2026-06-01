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
    _soft_floor_logl,
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
        """Without an explicit RNG, distinct calls give distinct floor values
        — what gives the sampler a gradient.
        """
        floors = {_soft_floor_logl(1.0) for _ in range(50)}
        # Allow up to one duplicate by chance; >=49/50 unique is the bar
        assert len(floors) >= 48

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
