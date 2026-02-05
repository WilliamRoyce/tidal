"""Tests for parameter name aliases in initial conditions.

This module verifies that initial condition classes (GaussianPulse, RingPulse2D)
correctly handle both 'width' and 'sigma' parameter names as aliases,
and validates NaN/Inf rejection for numeric parameters.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest

from torsion_gertsenshtein.kgsim import (
    GaussianPulse,
    GridConfig,
    RingPulse2D,
    make_grid,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


# ==================== Width/Sigma Alias Tests ====================


class TestWidthSigmaAliases:
    """Tests for width/sigma parameter aliases in initial conditions."""

    @pytest.fixture
    def grid_1d(self) -> CartesianGrid:
        """Create 1D grid for Gaussian tests.

        Returns
        -------
        CartesianGrid
            64-cell grid with bounds (0, 50).
        """
        return make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    @pytest.fixture
    def grid_2d(self) -> CartesianGrid:
        """Create 2D grid for RingPulse tests.

        Returns
        -------
        CartesianGrid
            32x32 grid with bounds (-10, 10) x (-10, 10).
        """
        return make_grid(
            GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
        )

    def test_gaussian_width_sigma_equivalent(self, grid_1d: CartesianGrid) -> None:
        """Verify width and sigma parameters produce identical results for GaussianPulse."""
        ic_width = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0])
        state_width = ic_width.build(grid_1d)

        ic_sigma = GaussianPulse(amplitude=1.0, sigma=2.0, center=[25.0])
        state_sigma = ic_sigma.build(grid_1d)

        assert np.allclose(state_width[0].data, state_sigma[0].data)
        assert np.allclose(state_width[1].data, state_sigma[1].data)

        # Both should store as width internally
        assert ic_width.width == 2.0  # noqa: PLR2004
        assert ic_sigma.width == 2.0  # noqa: PLR2004

    def test_gaussian_both_params_raises(self) -> None:
        """Verify specifying both width and sigma raises ValueError."""
        with pytest.raises(ValueError, match="Cannot specify both 'width' and 'sigma'"):
            GaussianPulse(amplitude=1.0, width=2.0, sigma=3.0, center=[25.0])

    def test_ring_width_sigma_equivalent(self, grid_2d: CartesianGrid) -> None:
        """Verify width and sigma parameters produce identical results for RingPulse2D."""
        ic_sigma = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.5)
        state_sigma = ic_sigma.build(grid_2d)

        ic_width = RingPulse2D(amplitude=1.0, initial_radius=5.0, width=1.5)
        state_width = ic_width.build(grid_2d)

        assert np.allclose(state_sigma[0].data, state_width[0].data)
        assert np.allclose(state_sigma[1].data, state_width[1].data)

        # Both should store as width internally
        assert ic_sigma.width == 1.5  # noqa: PLR2004
        assert ic_width.width == 1.5  # noqa: PLR2004

    def test_ring_both_params_raises(self) -> None:
        """Verify specifying both sigma and width raises ValueError."""
        with pytest.raises(ValueError, match="Cannot specify both 'sigma' and 'width'"):
            RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.5, width=2.5)

    @pytest.mark.parametrize("use_sigma", [True, False], ids=["sigma", "width"])
    def test_consistent_naming_across_ics(
        self, grid_1d: CartesianGrid, grid_2d: CartesianGrid, use_sigma: bool  # noqa: FBT001
    ) -> None:
        """Verify both ICs can use either width or sigma consistently."""
        if use_sigma:
            ic_gaussian = GaussianPulse(amplitude=1.0, sigma=2.0)
            ic_ring = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=2.0)
        else:
            ic_gaussian = GaussianPulse(amplitude=1.0, width=2.0)
            ic_ring = RingPulse2D(amplitude=1.0, initial_radius=5.0, width=2.0)

        state_gaussian = ic_gaussian.build(grid_1d)
        state_ring = ic_ring.build(grid_2d)

        assert len(state_gaussian) == 2  # noqa: PLR2004
        assert len(state_ring) == 2  # noqa: PLR2004
        assert ic_gaussian.width == 2.0  # noqa: PLR2004
        assert ic_ring.width == 2.0  # noqa: PLR2004


# ==================== NaN/Inf Validation Tests ====================


class TestNaNInfValidation:
    """Tests for NaN and infinity validation in initial condition parameters."""

    def test_gaussian_nan_amplitude_raises(self) -> None:
        """Verify NaN amplitude raises ValueError."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("nan"), width=2.0)

    def test_gaussian_inf_amplitude_raises(self) -> None:
        """Verify infinite amplitude raises ValueError."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("inf"), width=2.0)

    def test_gaussian_nan_width_raises(self) -> None:
        """Verify NaN width raises ValueError."""
        with pytest.raises(ValueError, match="width must be finite"):
            GaussianPulse(amplitude=1.0, width=float("nan"))

    def test_gaussian_inf_width_raises(self) -> None:
        """Verify infinite width raises ValueError."""
        with pytest.raises(ValueError, match="width must be finite"):
            GaussianPulse(amplitude=1.0, width=float("inf"))

    def test_gaussian_nan_initial_velocity_raises(self) -> None:
        """Verify NaN initial_velocity raises ValueError."""
        with pytest.raises(ValueError, match="initial_velocity must be finite"):
            GaussianPulse(amplitude=1.0, width=2.0, initial_velocity=float("nan"))

    def test_ring_nan_amplitude_raises(self) -> None:
        """Verify NaN amplitude raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            RingPulse2D(amplitude=float("nan"), initial_radius=5.0, sigma=1.0)

    def test_ring_inf_sigma_raises(self) -> None:
        """Verify infinite sigma raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="sigma must be finite"):
            RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=float("inf"))

    def test_ring_nan_initial_radius_raises(self) -> None:
        """Verify NaN initial_radius raises ValueError for RingPulse2D."""
        with pytest.raises(ValueError, match="initial_radius must be finite"):
            RingPulse2D(amplitude=1.0, initial_radius=float("nan"), sigma=1.0)

    def test_negative_inf_also_rejected(self) -> None:
        """Verify negative infinity is also rejected."""
        with pytest.raises(ValueError, match="amplitude must be finite"):
            GaussianPulse(amplitude=float("-inf"), width=2.0)
