"""Tests to verify RingPulse2D uses grid center, not origin.

This module ensures that the RingPulse2D initial condition correctly
centers the ring at the grid midpoint regardless of grid bounds.
"""

from __future__ import annotations

import numpy as np
import pytest
from pde import CartesianGrid

from tidal.kgsim.initial_conditions import RingPulse2D

# ==================== Ring Centering Tests ====================


class TestRingCentering:
    """Tests for RingPulse2D centering behavior."""

    def test_centered_at_grid_midpoint(self) -> None:
        """Verify RingPulse2D centers the ring at grid midpoint.

        For a grid centered at origin, the ring should be centered at (0, 0)
        and exhibit 4-fold symmetry.
        """
        grid = CartesianGrid([(-10, 10), (-10, 10)], [64, 64])

        ic = RingPulse2D(amplitude=1.0, initial_radius=5.0, sigma=1.0)
        state = ic.build(grid)
        phi = state[0].data

        # Maximum should be approximately equal to the amplitude
        max_val = phi.max()
        assert max_val == pytest.approx(1.0, abs=0.01)

        # Verify 4-fold symmetry about center
        assert np.allclose(phi, phi[::-1, :]), "Not symmetric vertically"
        assert np.allclose(phi, phi[:, ::-1]), "Not symmetric horizontally"

    def test_shifted_grid(self) -> None:
        """Verify RingPulse2D works correctly on a shifted grid.

        Even when the grid is not centered at the origin, the ring should
        be centered at the grid midpoint and maintain symmetry.
        """
        # Grid centered at (20, 30) instead of origin
        grid = CartesianGrid([(10, 30), (20, 40)], [64, 64])

        ic = RingPulse2D(amplitude=2.0, initial_radius=3.0, sigma=0.5)
        state = ic.build(grid)
        phi = state[0].data

        # Maximum should still be at amplitude
        assert phi.max() == pytest.approx(2.0, abs=0.01)

        # Field should be centered at grid midpoint (20, 30)
        assert np.allclose(phi, phi[::-1, :]), "Not symmetric vertically"
        assert np.allclose(phi, phi[:, ::-1]), "Not symmetric horizontally"
