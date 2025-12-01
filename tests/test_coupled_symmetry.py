from __future__ import annotations

import numpy as np

from torsion_gertsenshtein.kgsim import GridConfig, SimulationConfig, make_grid, run
from torsion_gertsenshtein.kgsim.config import MultiFieldParams
from torsion_gertsenshtein.kgsim.equations import make_coupled_kg_pde
from torsion_gertsenshtein.kgsim.initial_conditions import multi_gaussian


def test_two_field_symmetry() -> None:
    """Run a short symmetric two-field coupled KG simulation and assert phi0 == phi1."""
    # small, fast simulation
    grid = make_grid(
        GridConfig(dim=1, shape=(32,), bounds=((0.0, 20.0),), periodic=True)
    )

    masses = [0.5, 0.5]
    g = 0.1
    coupling = [[0.0, g], [g, 0.0]]
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # identical initial data for both fields -> symmetry should be preserved
    state = multi_gaussian(
        grid, amplitudes=[1.0, 1.0], widths=[5.0, 5.0], velocities=[0.0, 0.0]
    )

    simulation_config = SimulationConfig(
        t_end=2.0,
        dt=0.01,
        solver="explicit",  # fast and deterministic for tests
        method="RK4",
        backend="numpy",
        progress=False,
    )

    # run in-place; final fields live in `state`
    run(pde=pde, state=state, config=simulation_config)

    phi0 = np.asarray(state[0].data).ravel()
    phi1 = np.asarray(state[2].data).ravel()

    # fields must match to numerical tolerance
    np.testing.assert_allclose(phi0, phi1, rtol=1e-6, atol=1e-8)
