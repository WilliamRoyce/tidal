"""Tests for InitialCondition class-based API."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from typing_extensions import override

from torsion_gertsenshtein.kgsim import (
    GaussianPulse,
    GridConfig,
    InitialCondition,
    RingPulse2D,
    make_grid,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


def test_gaussian_pulse_class_basic() -> None:
    """Test GaussianPulse class creates valid initial condition."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    ic = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0])
    state = ic.build(grid)

    # Check structure
    expected_field_count = 2
    assert len(state) == expected_field_count
    assert state[0].data.shape == (64,)
    assert state[1].data.shape == (64,)

    # Check peak is roughly at center
    peak_idx = np.argmax(state[0].data)
    peak_index_min = 28
    peak_index_max = 36
    assert (
        peak_index_min < peak_idx < peak_index_max
    )  # Should be near index 32 (center)

    # Check amplitude (with reasonable tolerance for grid discretization)
    assert np.isclose(state[0].data.max(), 1.0, atol=0.05)


def test_gaussian_pulse_with_velocity() -> None:
    """Test GaussianPulse with initial velocity."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    ic = GaussianPulse(amplitude=1.0, width=2.0, center=[25.0], initial_velocity=0.5)
    state = ic.build(grid)

    # Check pi field is scaled version of phi
    assert np.allclose(state[1].data, 0.5 * state[0].data)


def test_gaussian_pulse_validates_width() -> None:
    """Test GaussianPulse rejects non-positive width."""
    with pytest.raises(ValueError, match="width must be positive"):
        GaussianPulse(amplitude=1.0, width=0.0)

    with pytest.raises(ValueError, match="width must be positive"):
        GaussianPulse(amplitude=1.0, width=-1.0)


def test_gaussian_pulse_2d() -> None:
    """Test GaussianPulse works in 2D."""
    grid = make_grid(
        GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
    )

    ic = GaussianPulse(amplitude=1.0, width=2.0, center=[0.0, 0.0])
    state = ic.build(grid)

    # Check structure
    assert state[0].data.shape == (32, 32)
    assert state[1].data.shape == (32, 32)

    # Check peak near center (with reasonable tolerance for grid discretization)
    assert np.isclose(state[0].data.max(), 1.0, atol=0.05)


def test_ring_pulse_2d_basic() -> None:
    """Test RingPulse2D creates valid ring-shaped IC."""
    grid = make_grid(
        GridConfig(dim=2, shape=(64, 64), bounds=((-20.0, 20.0), (-20.0, 20.0)))
    )

    ic = RingPulse2D(amplitude=1.0, initial_radius=8.0, sigma=1.0)
    state = ic.build(grid)

    # Check structure
    assert state[0].data.shape == (64, 64)
    assert state[1].data.shape == (64, 64)

    # Check peak amplitude
    assert np.isclose(state[0].data.max(), 1.0, atol=0.05)

    # Check center (origin) is lower than peak (ring)
    center_idx = 32  # Middle of 64x64 grid
    center_value = state[0].data[center_idx, center_idx]
    assert center_value < 0.5 * state[0].data.max()


def test_ring_pulse_2d_requires_2d_grid() -> None:
    """Test RingPulse2D raises error for non-2D grid."""
    grid_1d = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)

    with pytest.raises(ValueError, match="requires 2D grid"):
        ic.build(grid_1d)


def test_ring_pulse_2d_validates_parameters() -> None:
    """Test RingPulse2D parameter validation."""
    with pytest.raises(ValueError, match="sigma must be positive"):
        RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=0.0)

    with pytest.raises(ValueError, match="sigma must be positive"):
        RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=-1.0)

    with pytest.raises(ValueError, match="initial_radius must be non-negative"):
        RingPulse2D(amplitude=1.0, initial_radius=-5.0, sigma=1.0)


def test_initial_condition_custom_subclass() -> None:
    """Test custom InitialCondition subclass."""

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
    assert np.allclose(state[1].data, 0.0)  # Default pi is zero


def test_initial_condition_distance_computation() -> None:
    """Test distance computation helper."""
    grid = make_grid(
        GridConfig(dim=2, shape=(16, 16), bounds=((0.0, 10.0), (0.0, 10.0)))
    )

    ic = GaussianPulse(amplitude=1.0, width=2.0)

    # Test distance from center (accessing protected method for testing purposes)
    distances = ic._compute_distances_from_center(grid, center=[5.0, 5.0])  # noqa: SLF001  # pyright: ignore[reportPrivateUsage]

    # Check distances are positive
    assert np.all(distances >= 0)

    # Check minimum distance is near center
    assert distances.min() < 1.0  # Should have points very close to center


def test_gaussian_pulse_auto_center() -> None:
    """Test GaussianPulse with automatic centering."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    # Create with no center specified (should use grid center)
    ic = GaussianPulse(amplitude=1.0, width=2.0, center=None)
    state = ic.build(grid)

    # Peak should be roughly in middle
    peak_idx = np.argmax(state[0].data)
    peak_index_min = 28
    peak_index_max = 36
    assert peak_index_min < peak_idx < peak_index_max


def test_multiple_ics_with_same_grid() -> None:
    """Test creating multiple ICs with same grid."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    ic1 = GaussianPulse(amplitude=1.0, width=2.0, center=[15.0])
    ic2 = GaussianPulse(amplitude=0.5, width=3.0, center=[35.0])

    state1 = ic1.build(grid)
    state2 = ic2.build(grid)

    # Check they're different
    assert not np.allclose(state1[0].data, state2[0].data)

    # Check amplitudes
    assert np.isclose(state1[0].data.max(), 1.0, atol=0.01)
    assert np.isclose(state2[0].data.max(), 0.5, atol=0.01)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
