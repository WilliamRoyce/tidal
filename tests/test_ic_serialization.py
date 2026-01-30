"""Tests for serialization and thread safety of initial conditions."""

from __future__ import annotations

import copy
import pickle  # noqa: S403
import threading
from typing import TYPE_CHECKING

import numpy as np

from torsion_gertsenshtein.kgsim import (
    GaussianPulse,
    GridConfig,
    RingPulse2D,
    make_grid,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


def test_pickle_roundtrip() -> None:
    """Test that InitialCondition instances can be pickled and unpickled."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    # Create and use IC
    ic = GaussianPulse(amplitude=1.5, width=2.5, center=[25.0], initial_velocity=0.1)
    state_before = ic.build(grid)

    # Pickle and unpickle
    pickled = pickle.dumps(ic)
    ic_restored = pickle.loads(pickled)  # noqa: S301

    # Verify parameters preserved
    assert ic_restored.amplitude == 1.5  # noqa: PLR2004
    assert ic_restored.width == 2.5  # noqa: PLR2004
    assert ic_restored.center == [25.0]
    assert ic_restored.initial_velocity == 0.1  # noqa: PLR2004

    # Verify can build after unpickle
    state_after = ic_restored.build(grid)

    # Results should be identical
    assert np.allclose(state_before[0].data, state_after[0].data)
    assert np.allclose(state_before[1].data, state_after[1].data)


def test_deepcopy() -> None:
    """Test that InitialCondition instances can be deepcopied."""
    grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

    # Create and use IC
    ic = GaussianPulse(amplitude=2.0, width=3.0, center=[25.0])
    state_original = ic.build(grid)

    # Deep copy
    ic_copy = copy.deepcopy(ic)

    # Modify original parameters
    ic.amplitude = 999.0

    # Copy should be independent
    assert ic_copy.amplitude == 2.0  # noqa: PLR2004

    # Copy should still work
    state_copy = ic_copy.build(grid)
    assert np.allclose(state_original[0].data, state_copy[0].data)


def test_pickle_ringpulse2d() -> None:
    """Test pickle roundtrip for RingPulse2D."""
    grid = make_grid(
        GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
    )

    ic = RingPulse2D(amplitude=1.5, initial_radius=5.0, sigma=1.5)
    state_before = ic.build(grid)

    # Pickle roundtrip
    ic_restored = pickle.loads(pickle.dumps(ic))  # noqa: S301

    # Verify parameters and output
    assert ic_restored.amplitude == 1.5  # noqa: PLR2004
    assert ic_restored.initial_radius == 5.0  # noqa: PLR2004
    assert ic_restored.width == 1.5  # noqa: PLR2004

    state_after = ic_restored.build(grid)
    assert np.allclose(state_before[0].data, state_after[0].data)


def test_thread_safety_concurrent_builds() -> None:
    """Test that concurrent builds from different threads produce correct results."""
    ic = GaussianPulse(amplitude=1.0, width=2.0)
    grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0.0, 50.0),)))

    results: list[np.ndarray | None] = [None] * 10
    errors: list[Exception] = []

    def build_ic(index: int) -> None:
        """Build IC in a thread."""
        try:
            state = ic.build(grid)
            results[index] = state[0].data.copy()
        except Exception as e:  # noqa: BLE001
            errors.append(e)

    # Launch multiple threads building on same grid
    threads = [threading.Thread(target=build_ic, args=(i,)) for i in range(10)]

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

    # Check for errors
    assert not errors, f"Thread safety errors: {errors}"

    # All results should be identical
    assert all(r is not None for r in results), "Some threads failed to produce results"
    reference = results[0]
    assert reference is not None  # Type narrowing for type checker
    for i, result in enumerate(results[1:], start=1):
        assert result is not None  # Type narrowing for type checker
        assert np.allclose(
            result,
            reference,
        ), f"Thread {i} produced different result"


def test_backward_compatibility_custom_subclass() -> None:
    """Test that custom subclasses work without special initialization."""
    from torsion_gertsenshtein.kgsim import InitialCondition  # noqa: PLC0415

    class CustomIC(InitialCondition):
        """Custom IC with simple implementation."""

        def __init__(self, value: float = 1.0) -> None:
            """Initialize with a constant value."""
            self.value = value

        def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
            """Return constant field."""
            return np.full(grid.shape, self.value).ravel()

    grid = make_grid(GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),)))
    ic = CustomIC(value=42.0)

    state = ic.build(grid)
    assert state[0].data.shape == (32,)
    assert np.allclose(state[0].data, 42.0)
