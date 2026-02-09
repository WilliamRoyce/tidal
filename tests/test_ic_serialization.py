"""Tests for serialization and thread safety of initial conditions.

This module verifies that InitialCondition instances can be pickled,
deepcopied, and used safely from multiple threads. Also tests
backward compatibility for custom subclasses.
"""

from __future__ import annotations

import copy
import pickle  # noqa: S403  # pickle is intentional for testing serialization
import threading
from typing import TYPE_CHECKING

import numpy as np
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


# ==================== Serialization Tests ====================


class TestPickleSerialization:
    """Tests for pickle serialization of initial conditions."""

    def test_gaussian_roundtrip(self) -> None:
        """Verify GaussianPulse can be pickled and unpickled correctly."""
        grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

        ic = GaussianPulse(
            amplitude=1.5, width=2.5, center=[25.0], initial_velocity=0.1
        )
        state_before = ic.build(grid)

        # Pickle and unpickle
        pickled = pickle.dumps(ic)
        ic_restored = pickle.loads(pickled)  # noqa: S301

        # Verify parameters preserved
        assert ic_restored.amplitude == 1.5
        assert ic_restored.width == 2.5
        assert ic_restored.center == [25.0]
        assert ic_restored.initial_velocity == 0.1

        # Verify can build after unpickle
        state_after = ic_restored.build(grid)

        assert np.allclose(state_before[0].data, state_after[0].data)
        assert np.allclose(state_before[1].data, state_after[1].data)

    def test_ringpulse_roundtrip(self) -> None:
        """Verify RingPulse2D can be pickled and unpickled correctly."""
        grid = make_grid(
            GridConfig(dim=2, shape=(32, 32), bounds=((-10.0, 10.0), (-10.0, 10.0)))
        )

        ic = RingPulse2D(amplitude=1.5, initial_radius=5.0, sigma=1.5)
        state_before = ic.build(grid)

        ic_restored = pickle.loads(pickle.dumps(ic))  # noqa: S301

        assert ic_restored.amplitude == 1.5
        assert ic_restored.initial_radius == 5.0
        assert ic_restored.width == 1.5

        state_after = ic_restored.build(grid)
        assert np.allclose(state_before[0].data, state_after[0].data)


class TestDeepCopy:
    """Tests for deep copy behavior of initial conditions."""

    def test_deepcopy_independence(self) -> None:
        """Verify deepcopied IC is independent from original."""
        grid = make_grid(GridConfig(dim=1, shape=(64,), bounds=((0.0, 50.0),)))

        ic = GaussianPulse(amplitude=2.0, width=3.0, center=[25.0])
        state_original = ic.build(grid)

        ic_copy = copy.deepcopy(ic)

        # Modify original parameters
        ic.amplitude = 999.0

        # Copy should be independent
        assert ic_copy.amplitude == 2.0

        # Copy should still work and match original output
        state_copy = ic_copy.build(grid)
        assert np.allclose(state_original[0].data, state_copy[0].data)


# ==================== Thread Safety Tests ====================


class TestThreadSafety:
    """Tests for thread-safe usage of initial conditions."""

    def test_concurrent_builds(self) -> None:
        """Verify concurrent builds from different threads produce correct results."""
        ic = GaussianPulse(amplitude=1.0, width=2.0)
        grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0.0, 50.0),)))

        num_threads = 10
        results: list[np.ndarray | None] = [None] * num_threads
        errors: list[Exception] = []

        def build_ic(index: int) -> None:
            """Build IC in a thread."""
            try:
                state = ic.build(grid)
                results[index] = state[0].data.copy()
            except Exception as e:  # noqa: BLE001
                errors.append(e)

        threads = [
            threading.Thread(target=build_ic, args=(i,)) for i in range(num_threads)
        ]

        for thread in threads:
            thread.start()

        for thread in threads:
            thread.join()

        assert not errors, f"Thread safety errors: {errors}"
        assert all(r is not None for r in results), "Some threads failed"

        reference = results[0]
        assert reference is not None
        for i, result in enumerate(results[1:], start=1):
            assert result is not None
            assert np.allclose(result, reference), (
                f"Thread {i} produced different result"
            )


# ==================== Backward Compatibility Tests ====================


class TestBackwardCompatibility:
    """Tests for backward compatibility with custom subclasses."""

    def test_custom_subclass_works(self) -> None:
        """Verify custom InitialCondition subclasses work without special initialization."""

        class CustomIC(InitialCondition):
            """Custom IC with simple constant value implementation."""

            def __init__(self, value: float = 1.0) -> None:
                """Initialize with a constant value."""
                self.value = value

            @override
            def _compute_phi(self, grid: CartesianGrid) -> np.ndarray:
                """Return constant field."""
                return np.full(grid.shape, self.value).ravel()

        grid = make_grid(GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),)))
        ic = CustomIC(value=42.0)

        state = ic.build(grid)
        assert state[0].data.shape == (32,)
        assert np.allclose(state[0].data, 42.0)
