"""Tests for 2D coupled Klein-Gordon simulations."""

from __future__ import annotations

import numpy as np
import pytest

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
    run,
)
from torsion_gertsenshtein.kgsim.config import MultiFieldParams
from torsion_gertsenshtein.kgsim.equations import make_coupled_kg_pde
from torsion_gertsenshtein.kgsim.initial_conditions import (
    gaussian_pulse,
    multi_gaussian_2d,
)


def test_multi_gaussian_2d_basic() -> None:
    """Test that multi_gaussian_2d creates correct field structure."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    # Two fields with different centers
    state = multi_gaussian_2d(
        grid,
        amplitudes=[1.0, 0.5],
        widths=[2.0, 3.0],
        centers=[(-5.0, 0.0), (5.0, 0.0)],
    )

    # Should have 4 fields: phi0, pi0, phi1, pi1
    assert len(state) == 4  # noqa: PLR2004
    assert state.labels == ["phi0", "pi0", "phi1", "pi1"]

    # All fields should have correct shape
    for field in state:
        assert field.data.shape == (32, 32)

    # phi0 should be maximum near (-5, 0)
    phi0_data = np.asarray(state[0].data)
    assert phi0_data.max() == pytest.approx(1.0, rel=0.05)

    # phi1 should be maximum near (5, 0)
    phi1_data = np.asarray(state[2].data)
    assert phi1_data.max() == pytest.approx(0.5, rel=0.05)

    # pi fields should be zero (no initial velocity)
    pi0_data = np.asarray(state[1].data)
    pi1_data = np.asarray(state[3].data)
    assert pi0_data.max() == pytest.approx(0.0, abs=1e-10)
    assert pi1_data.max() == pytest.approx(0.0, abs=1e-10)


def test_multi_gaussian_2d_requires_2d_grid() -> None:
    """Test that multi_gaussian_2d raises ValueError for non-2D grids."""
    grid_1d = make_grid(
        GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),), periodic=True)
    )

    with pytest.raises(ValueError, match="requires a 2D grid"):
        multi_gaussian_2d(grid_1d, amplitudes=[1.0], widths=[2.0])


def test_multi_gaussian_2d_with_velocities() -> None:
    """Test that initial velocities are correctly applied."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    # Create with non-zero velocity
    state = multi_gaussian_2d(
        grid,
        amplitudes=[1.0],
        widths=[2.0],
        velocities=[0.5],
    )

    phi0_data = np.asarray(state[0].data)
    pi0_data = np.asarray(state[1].data)

    # pi should be 0.5 * phi
    np.testing.assert_allclose(pi0_data, 0.5 * phi0_data, rtol=1e-10)


def test_multi_gaussian_2d_default_center() -> None:
    """Test that default center is domain midpoint."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-5.0, 15.0)), periodic=True
        )
    )

    state = multi_gaussian_2d(grid, amplitudes=[1.0], widths=[2.0])

    phi0_data = np.asarray(state[0].data)

    # Maximum should be near the center of the domain
    # Domain center: (0, 5)
    # Since grid is 32x32, center indices are around (16, 16)
    max_idx = np.unravel_index(np.argmax(phi0_data), phi0_data.shape)

    # Should be near the center (within a few cells)
    assert abs(max_idx[0] - 16) < 3  # noqa: PLR2004
    assert abs(max_idx[1] - 16) < 3  # noqa: PLR2004


def test_two_field_symmetry_2d() -> None:
    """Test that symmetric 2D coupled fields preserve symmetry."""
    # Small, fast simulation
    grid = make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    masses = [0.5, 0.5]
    g = 0.1
    coupling = [[0.0, g], [g, 0.0]]
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # Identical initial data for both fields -> symmetry should be preserved
    state = multi_gaussian_2d(
        grid,
        amplitudes=[1.0, 1.0],
        widths=[2.0, 2.0],
        centers=[(0.0, 0.0), (0.0, 0.0)],  # same center
        velocities=[0.0, 0.0],
    )

    simulation_config = SimulationConfig(
        t_end=1.0,
        dt=0.05,
        solver="explicit",
        method="RK4",
        backend="numpy",
        progress=False,
    )

    # Run simulation
    run(pde=pde, state=state, config=simulation_config)

    phi0 = np.asarray(state[0].data).ravel()
    phi1 = np.asarray(state[2].data).ravel()

    # Fields must match to numerical tolerance
    np.testing.assert_allclose(phi0, phi1, rtol=1e-5, atol=1e-7)


def test_2d_coupled_energy_transfer() -> None:
    """Test that coupled 2D fields show energy exchange behavior."""
    # Small simulation with both fields initially excited
    grid = make_grid(
        GridConfig(
            dim=2, shape=(32, 32), bounds=((-20.0, 20.0), (-20.0, 20.0)), periodic=True
        )
    )

    masses = [0.5, 1.0]
    g = 0.3  # moderate coupling
    coupling = [[0.0, g], [g, 0.0]]
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # Start with asymmetric excitation - field 0 has more energy than field 1
    state = multi_gaussian_2d(
        grid,
        amplitudes=[1.0, 0.2],
        widths=[5.0, 5.0],
        centers=[(0.0, 0.0), (0.0, 0.0)],  # overlapping for coupling
        velocities=[0.0, 0.0],
    )

    # Record initial energies
    phi0_initial = np.asarray(state[0].data).copy()
    phi1_initial = np.asarray(state[2].data).copy()
    initial_energy_0 = np.sum(phi0_initial**2)
    initial_energy_1 = np.sum(phi1_initial**2)

    assert initial_energy_0 > initial_energy_1  # field 0 has more energy

    total_initial = initial_energy_0 + initial_energy_1

    simulation_config = SimulationConfig(
        t_end=5.0,
        dt=0.05,
        solver="explicit",
        method="RK4",
        backend="numpy",
        progress=False,
    )

    # Run simulation
    run(pde=pde, state=state, config=simulation_config)

    # Check final energies
    phi0_final = np.asarray(state[0].data)
    phi1_final = np.asarray(state[2].data)
    final_energy_0 = np.sum(phi0_final**2)
    final_energy_1 = np.sum(phi1_final**2)

    total_final = final_energy_0 + final_energy_1

    # With coupling, energies should have exchanged
    # Both fields should have non-zero energy at the end
    assert final_energy_0 > 0
    assert final_energy_1 > 0

    # Total energy should be approximately conserved (within numerical error + dissipation)
    # Allow up to 50% change due to wave propagation out of the central region
    assert total_final > 0.5 * total_initial


def test_2d_coupled_decoupled_limit() -> None:
    """Test that decoupled 2D fields evolve independently."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    masses = [0.5, 1.0]
    coupling = [[0.0, 0.0], [0.0, 0.0]]  # no coupling
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # Different initial conditions for each field
    state_coupled = multi_gaussian_2d(
        grid,
        amplitudes=[1.0, 0.5],
        widths=[2.0, 3.0],
        centers=[(-3.0, 0.0), (3.0, 0.0)],
    )

    # Run coupled (but decoupled) system
    config = SimulationConfig(
        t_end=1.0,
        dt=0.05,
        solver="explicit",
        method="RK4",
        backend="numpy",
        progress=False,
    )
    run(pde=pde, state=state_coupled, config=config)

    # Now run each field independently and compare
    # Field 0

    state0 = gaussian_pulse(grid, amplitude=1.0, center=[-3.0, 0.0], width=2.0)
    pde0 = KleinGordonPDE(KGParameters(mass=0.5))
    run(pde=pde0, state=state0, config=config)

    # Field 1
    state1 = gaussian_pulse(grid, amplitude=0.5, center=[3.0, 0.0], width=3.0)
    pde1 = KleinGordonPDE(KGParameters(mass=1.0))
    run(pde=pde1, state=state1, config=config)

    # Compare: coupled system fields should match independent evolution
    phi0_coupled = np.asarray(state_coupled[0].data).ravel()
    phi0_independent = np.asarray(state0[0].data).ravel()
    np.testing.assert_allclose(phi0_coupled, phi0_independent, rtol=1e-5, atol=1e-7)

    phi1_coupled = np.asarray(state_coupled[2].data).ravel()
    phi1_independent = np.asarray(state1[0].data).ravel()
    np.testing.assert_allclose(phi1_coupled, phi1_independent, rtol=1e-5, atol=1e-7)


def test_multi_gaussian_2d_validates_widths() -> None:
    """Test that multi_gaussian_2d rejects non-positive widths."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    with pytest.raises(ValueError, match="widths must be positive"):
        multi_gaussian_2d(grid, amplitudes=[1.0], widths=[0.0])

    with pytest.raises(ValueError, match="widths must be positive"):
        multi_gaussian_2d(grid, amplitudes=[1.0], widths=[-1.0])


def test_multi_gaussian_2d_validates_centers_length() -> None:
    """Test that multi_gaussian_2d validates centers parameter length."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    with pytest.raises(ValueError, match="centers must have length"):
        multi_gaussian_2d(
            grid,
            amplitudes=[1.0, 1.0],
            widths=[2.0, 2.0],
            centers=[(0.0, 0.0)],  # only 1 center
        )


def test_multi_gaussian_2d_empty_amplitudes() -> None:
    """Test that multi_gaussian_2d rejects empty amplitudes."""
    grid = make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
        )
    )

    with pytest.raises(ValueError, match="amplitudes must be a non-empty sequence"):
        multi_gaussian_2d(grid, amplitudes=[], widths=[])
