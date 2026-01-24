"""Tests for run_with_snapshots and MemoryStorage integration."""

from __future__ import annotations

import numpy as np
import pde as pde_module
import pytest
from pde import MemoryStorage

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
    run_with_snapshots,
)


def test_run_with_snapshots_basic() -> None:
    """Test that run_with_snapshots returns storage with correct snapshots."""
    # Simple 1D simulation
    grid_config = GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),), periodic=True)
    grid = make_grid(grid_config)
    state = gaussian_pulse(grid, amplitude=1.0, width=2.0, initial_velocity=0.0)
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=10.0, dt=0.1, solver="scipy", backend="numpy", progress=False
    )

    # Run with snapshots
    result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=1.0
    )

    # Check storage is MemoryStorage instance
    assert isinstance(storage, MemoryStorage)

    # Check we got snapshots
    assert len(storage) > 0
    min_expected_snapshots = 10
    assert (
        len(storage) >= min_expected_snapshots
    )  # Should have ~10-11 snapshots for t_end=10, interval=1

    # Check times are correct
    assert storage.times[0] == pytest.approx(0.0, abs=0.1)
    assert storage.times[-1] == pytest.approx(10.0, abs=0.1)

    # Check final state matches storage last entry
    assert np.allclose(result[0].data, storage[-1][0].data)  # type: ignore[index]
    assert np.allclose(result[1].data, storage[-1][1].data)  # type: ignore[index]


def test_run_with_snapshots_vs_manual_pattern() -> None:
    """Test that run_with_snapshots gives equivalent results to manual pattern."""
    grid_config = GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),), periodic=True)
    grid = make_grid(grid_config)
    state1 = gaussian_pulse(grid, amplitude=1.0, width=2.0, initial_velocity=0.0)
    state2 = gaussian_pulse(grid, amplitude=1.0, width=2.0, initial_velocity=0.0)
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=5.0, dt=0.1, solver="scipy", backend="numpy", progress=False
    )

    # Method 1: run_with_snapshots
    result1, storage = run_with_snapshots(
        pde=pde, state=state1, config=config, snapshot_interval=0.5
    )

    # Method 2: manual pattern
    snapshots: list[tuple[float, np.ndarray]] = []
    snapshots.append((0.0, state2[0].data.copy()))

    def record_phi(
        state_coll: pde_module.FieldCollection, t: float
    ) -> dict[str, object]:
        """Record phi field at each timestep."""
        snapshots.append((t, state_coll[0].data.copy()))  # type: ignore[index]
        return {}

    result2 = run(
        pde=pde,
        state=state2,
        config=config,
        snapshot_interval=0.5,
        extra_observer=record_phi,
    )

    # Compare final states
    assert np.allclose(result1[0].data, result2[0].data)
    assert np.allclose(result1[1].data, result2[1].data)

    # Compare snapshot counts (should be similar, might differ by 1)
    assert abs(len(storage) - len(snapshots)) <= 1

    # Compare snapshot data at matching times (skip first snapshot which may differ)
    # Storage includes t=0 automatically, manual pattern may not always match exactly
    # Compare the last few snapshots which should definitely match
    for i in range(-3, 0):  # Compare last 3 snapshots
        assert np.allclose(storage[i][0].data, snapshots[i][1], atol=1e-5)  # type: ignore[index]


def test_run_with_snapshots_different_intervals() -> None:
    """Test snapshot intervals are respected."""
    grid_config = GridConfig(dim=1, shape=(32,), bounds=((0.0, 20.0),), periodic=True)
    grid = make_grid(grid_config)
    state = gaussian_pulse(grid, amplitude=1.0, width=2.0)
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=10.0, dt=0.05, solver="scipy", backend="numpy", progress=False
    )

    # Test with interval=2.0
    _, storage_coarse = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=2.0
    )

    # Test with interval=0.5
    _, storage_fine = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=0.5
    )

    # Fine should have more snapshots than coarse
    assert len(storage_fine) > len(storage_coarse)

    # Rough check on snapshot counts
    min_coarse_snapshots = 5
    min_fine_snapshots = 20
    assert len(storage_coarse) >= min_coarse_snapshots  # ~10/2 = 5 snapshots
    assert len(storage_fine) >= min_fine_snapshots  # ~10/0.5 = 20 snapshots


def test_run_with_snapshots_preserves_field_structure() -> None:
    """Test that storage preserves FieldCollection structure."""
    grid_config = GridConfig(dim=1, shape=(32,), bounds=((0.0, 20.0),), periodic=True)
    grid = make_grid(grid_config)
    state = gaussian_pulse(grid, amplitude=1.0, width=2.0)
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=2.0, dt=0.1, solver="scipy", backend="numpy", progress=False
    )

    _result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=0.5
    )

    # Check each snapshot has correct structure
    expected_field_count = 2
    for i in range(len(storage)):
        snapshot = storage[i]
        # Should be FieldCollection with 2 fields
        assert isinstance(snapshot, pde_module.FieldCollection)
        assert len(snapshot) == expected_field_count
        # Check shapes match grid
        assert snapshot[0].data.shape == (32,)
        assert snapshot[1].data.shape == (32,)


def test_run_with_snapshots_2d() -> None:
    """Test run_with_snapshots works with 2D grids."""
    grid_config = GridConfig(
        dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)), periodic=True
    )
    grid = make_grid(grid_config)
    state = gaussian_pulse(grid, amplitude=1.0, width=2.0, center=[0.0, 0.0])
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    config = SimulationConfig(
        t_end=5.0, dt=0.1, solver="scipy", backend="numpy", progress=False
    )

    result, storage = run_with_snapshots(
        pde=pde, state=state, config=config, snapshot_interval=1.0
    )

    # Check storage
    min_snapshots_2d = 5
    assert len(storage) >= min_snapshots_2d
    # Check 2D shapes preserved
    snapshot_0 = storage[0]
    assert isinstance(snapshot_0, pde_module.FieldCollection)
    assert snapshot_0[0].data.shape == (32, 32)  # type: ignore[index]
    assert result[0].data.shape == (32, 32)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
