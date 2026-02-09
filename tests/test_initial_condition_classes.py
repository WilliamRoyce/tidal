"""Tests for InitialCondition class-based API.

This module tests the GaussianPulse and RingPulse2D initial condition classes,
including parameter validation, grid compatibility, and custom subclassing.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from typing_extensions import override

from tidal.kgsim import (
    GaussianPulse,
    GridConfig,
    InitialCondition,
    RingPulse2D,
    make_grid,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


# ==================== GaussianPulse Tests ====================


class TestGaussianPulse:
    """Tests for the GaussianPulse initial condition class."""

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
        """Create 2D grid for 2D Gaussian tests.

        Returns
        -------
        CartesianGrid
            32x32 grid with bounds (-10, 10) x (-10, 10).
        """
        return make_grid(
            GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
        )

    def test_basic(self, grid_1d: CartesianGrid) -> None:
        """Verify GaussianPulse creates valid initial condition."""
        ic = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0])
        state = ic.build(grid_1d)

        expected_field_count = 2
        assert len(state) == expected_field_count
        assert state[0].data.shape == (64,)
        assert state[1].data.shape == (64,)

        # Check peak is roughly at center
        peak_idx = np.argmax(state[0].data)
        assert 28 < peak_idx < 36

        # Check amplitude
        assert np.isclose(state[0].data.max(), 1.0, atol=0.05)

    def test_with_velocity(self, grid_1d: CartesianGrid) -> None:
        """Verify initial velocity scales pi field correctly."""
        ic = GaussianPulse(
            amplitude=1.0, width=2.0, center=[25.0], initial_velocity=0.5
        )
        state = ic.build(grid_1d)

        assert np.allclose(state[1].data, 0.5 * state[0].data)

    def test_validates_width(self) -> None:
        """Verify GaussianPulse rejects non-positive width."""
        with pytest.raises(ValueError, match="width must be positive"):
            GaussianPulse(amplitude=1.0, width=0.0)

        with pytest.raises(ValueError, match="width must be positive"):
            GaussianPulse(amplitude=1.0, width=-1.0)

    def test_2d_grid(self, grid_2d: CartesianGrid) -> None:
        """Verify GaussianPulse works in 2D."""
        ic = GaussianPulse(amplitude=1.0, width=2.0, center=[0.0, 0.0])
        state = ic.build(grid_2d)

        assert state[0].data.shape == (32, 32)
        assert state[1].data.shape == (32, 32)
        assert np.isclose(state[0].data.max(), 1.0, atol=0.05)

    def test_auto_center(self, grid_1d: CartesianGrid) -> None:
        """Verify GaussianPulse defaults to grid center when center=None."""
        ic = GaussianPulse(amplitude=1.0, width=2.0, center=None)
        state = ic.build(grid_1d)

        peak_idx = np.argmax(state[0].data)
        assert 28 < peak_idx < 36

    def test_multiple_ics_same_grid(self, grid_1d: CartesianGrid) -> None:
        """Verify multiple ICs can use the same grid."""
        ic1 = GaussianPulse(amplitude=1.0, width=2.0, center=[15.0])
        ic2 = GaussianPulse(amplitude=0.5, width=3.0, center=[35.0])

        state1 = ic1.build(grid_1d)
        state2 = ic2.build(grid_1d)

        assert not np.allclose(state1[0].data, state2[0].data)
        assert np.isclose(state1[0].data.max(), 1.0, atol=0.01)
        assert np.isclose(state2[0].data.max(), 0.5, atol=0.01)


# ==================== RingPulse2D Tests ====================


class TestRingPulse2D:
    """Tests for the RingPulse2D initial condition class."""

    @pytest.fixture
    def grid_2d(self) -> CartesianGrid:
        """Create 2D grid for ring pulse tests.

        Returns
        -------
        CartesianGrid
            64x64 grid with bounds (-20, 20) x (-20, 20).
        """
        return make_grid(
            GridConfig(dim=2, shape=(64, 64), bounds=((-20.0, 20.0), (-20.0, 20.0)))
        )

    @pytest.fixture
    def grid_1d(self) -> CartesianGrid:
        """Create 1D grid for dimension validation tests.

        Returns
        -------
        CartesianGrid
            64-cell grid with bounds (0, 50).
        """
        return make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    def test_basic(self, grid_2d: CartesianGrid) -> None:
        """Verify RingPulse2D creates valid ring-shaped IC."""
        ic = RingPulse2D(amplitude=1.0, initial_radius=8.0, sigma=1.0)
        state = ic.build(grid_2d)

        assert state[0].data.shape == (64, 64)
        assert state[1].data.shape == (64, 64)
        assert np.isclose(state[0].data.max(), 1.0, atol=0.05)

        # Check center is lower than ring peak
        center_idx = 32
        center_value = state[0].data[center_idx, center_idx]
        assert center_value < 0.5 * state[0].data.max()

    def test_requires_2d_grid(self, grid_1d: CartesianGrid) -> None:
        """Verify RingPulse2D raises error for non-2D grid."""
        ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)

        with pytest.raises(ValueError, match="requires 2D grid"):
            ic.build(grid_1d)

    def test_validates_parameters(self) -> None:
        """Verify RingPulse2D parameter validation."""
        with pytest.raises(ValueError, match="sigma must be positive"):
            RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=0.0)

        with pytest.raises(ValueError, match="sigma must be positive"):
            RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=-1.0)

        with pytest.raises(ValueError, match="initial_radius must be non-negative"):
            RingPulse2D(amplitude=1.0, initial_radius=-5.0, sigma=1.0)


# ==================== Custom Subclass Tests ====================


class TestCustomSubclass:
    """Tests for custom InitialCondition subclasses."""

    def test_constant_ic(self) -> None:
        """Verify custom constant field IC works correctly."""

        class ConstantIC(InitialCondition):
            """Simple constant field IC for testing."""

            def __init__(self, value: float) -> None:
                self.value = value

            @override
            def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
                coords = self._get_coordinates(grid)
                return np.full(len(coords), self.value)

        grid = make_grid(GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),)))
        ic = ConstantIC(value=5.0)
        state = ic.build(grid)

        assert np.allclose(state[0].data, 5.0)
        assert np.allclose(state[1].data, 0.0)


# ==================== Distance Computation Tests ====================


class TestDistanceComputation:
    """Tests for distance computation helper method."""

    def test_distances_from_center(self) -> None:
        """Verify distance computation produces correct values."""
        grid = make_grid(
            GridConfig(dim=2, shape=(16, 16), bounds=((0.0, 10.0), (0.0, 10.0)))
        )

        ic = GaussianPulse(amplitude=1.0, width=2.0)

        # Test distance from center (accessing protected method for testing)
        distances = ic._compute_distances_from_center(grid, center=[5.0, 5.0])  # pyright: ignore[reportPrivateUsage]

        assert np.all(distances >= 0)
        assert distances.min() < 1.0  # Should have points near center
