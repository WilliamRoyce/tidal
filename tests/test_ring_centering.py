"""Tests to verify RingPulse2D uses grid center, not origin."""

import numpy as np
import pytest
from pde import CartesianGrid

from torsion_gertsenshtein.kgsim.initial_conditions import RingPulse2D


def test_ring_pulse_2d_centered_at_grid_midpoint() -> None:
    """Verify that RingPulse2D centers the ring at grid midpoint."""
    # Create grid centered at origin
    grid = CartesianGrid([(-10, 10), (-10, 10)], [64, 64])

    ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)
    state = ic.build(grid)
    phi = state[0].data

    # Grid center is at (0, 0) based on bounds
    # Find the maximum value location (should be at radius 5 from center)
    max_val = phi.max()

    # Check that the maximum is approximately equal to the amplitude
    # (it should be at the ring radius from the center)
    assert max_val == pytest.approx(1.0, abs=0.01)

    # Verify 4-fold symmetry about center
    # The field should be horizontally and vertically symmetric
    assert np.allclose(phi, phi[::-1, :]), "Not symmetric vertically"
    assert np.allclose(phi, phi[:, ::-1]), "Not symmetric horizontally"


def test_ring_pulse_2d_shifted_grid() -> None:
    """Verify RingPulse2D works correctly on a shifted grid."""
    # Grid centered at (20, 30) instead of origin
    grid = CartesianGrid([(10, 30), (20, 40)], [64, 64])

    ic = RingPulse2D(amplitude=2.0, initial_radius=3.0, sigma=0.5)
    state = ic.build(grid)
    phi = state[0].data

    # Maximum should still be at amplitude
    assert phi.max() == pytest.approx(2.0, abs=0.01)

    # Field should be centered at grid midpoint (20, 30)
    # Verify horizontal and vertical symmetry
    assert np.allclose(phi, phi[::-1, :]), "Not symmetric vertically"
    assert np.allclose(phi, phi[:, ::-1]), "Not symmetric horizontally"
