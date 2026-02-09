"""Tests for symmetry preservation in coupled Klein-Gordon simulations.

This module verifies that symmetric initial conditions in coupled field
systems remain symmetric during evolution, validating the numerical
accuracy of the coupled field solver.
"""

from __future__ import annotations

import numpy as np

from tidal.kgsim import GridConfig, SimulationConfig, make_grid, run
from tidal.kgsim.config import MultiFieldParams
from tidal.kgsim.equations import make_coupled_kg_pde
from tidal.kgsim.initial_conditions import multi_gaussian

# ==================== Coupled Field Symmetry Tests ====================


class TestCoupledSymmetry:
    """Tests for symmetry preservation in coupled Klein-Gordon simulations."""

    def test_two_field_symmetry(self) -> None:
        """Verify symmetric two-field coupled KG simulation preserves phi0 == phi1.

        When two coupled fields start with identical initial conditions and
        identical masses, they should remain identical throughout evolution.
        """
        # Create a grid with appropriate bounds for the test
        grid = make_grid(
            GridConfig(dim=1, shape=(32,), bounds=((0.0, 20.0),), periodic=True)
        )

        masses = [0.5, 0.5]
        g = 0.1
        coupling = [[0.0, g], [g, 0.0]]
        pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

        # Identical initial data for both fields -> symmetry should be preserved
        state = multi_gaussian(
            grid, amplitudes=[1.0, 1.0], widths=[5.0, 5.0], velocities=[0.0, 0.0]
        )

        simulation_config = SimulationConfig(
            t_end=2.0,
            dt=0.01,
            solver="explicit",
            method="RK4",
            backend="numpy",
            progress=False,
        )

        run(pde=pde, state=state, config=simulation_config)

        phi0 = np.asarray(state[0].data).ravel()
        phi1 = np.asarray(state[2].data).ravel()

        # Fields must match to numerical tolerance
        np.testing.assert_allclose(phi0, phi1, rtol=1e-6, atol=1e-8)
