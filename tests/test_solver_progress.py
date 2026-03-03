"""Tests for SimulationProgress (tqdm-based progress bar)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.progress import SimulationProgress
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem


def _make_kg_spec() -> EquationSystem:
    """Klein-Gordon wave equation: d²φ/dt² = ∇²φ."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian",
                            "field": "phi_0",
                        },
                    ],
                },
            }
        ],
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                    "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                },
            ],
        },
    }
    return EquationSystem.from_dict(data)


class TestSimulationProgress:
    """Unit tests for the SimulationProgress class."""

    def test_basic_construction(self) -> None:
        """Progress bar can be created and closed without errors."""
        p = SimulationProgress(0.0, 10.0, solver_name="test", disable=True)
        p.update(5.0)
        p.finish()

    def test_disabled_is_noop(self) -> None:
        """When disable=True, update and finish succeed silently."""
        p = SimulationProgress(0.0, 10.0, disable=True)
        p.update(5.0)
        p.finish()

    def test_monotonic_tracking(self) -> None:
        """Progress never goes backwards (protects against IDA Jacobian FD)."""
        p = SimulationProgress(0.0, 10.0, disable=True)
        p.update(5.0)
        assert p._max_t == 5.0
        # Backwards t should be ignored
        p.update(4.0)
        assert p._max_t == 5.0
        # Forward t should advance
        p.update(7.0)
        assert p._max_t == 7.0
        p.finish()

    def test_finish_sets_complete(self) -> None:
        """finish() sets progress to 100%."""
        p = SimulationProgress(0.0, 10.0, disable=True)
        p.update(5.0)
        p.finish()
        assert p._pbar.n == p._pbar.total

    def test_nonzero_t_start(self) -> None:
        """Progress works correctly with non-zero t_start (e.g. resume)."""
        p = SimulationProgress(5.0, 15.0, disable=True)
        assert p._pbar.total == 10.0
        p.update(10.0)
        # Internal tqdm position should be 5.0 (halfway)
        assert p._pbar.n == pytest.approx(5.0)
        p.finish()
        assert p._pbar.n == pytest.approx(10.0)

    def test_solver_name_stored(self) -> None:
        """Solver name is stored on the progress object."""
        p = SimulationProgress(0.0, 10.0, solver_name="IDA", disable=True)
        # Disabled tqdm doesn't expose desc, but the progress object should
        # construct without error using the solver name
        p.update(5.0)
        p.finish()

    def test_default_desc_without_solver_name(self) -> None:
        """Without solver_name, construction succeeds."""
        p = SimulationProgress(0.0, 10.0, disable=True)
        p.update(5.0)
        p.finish()

    def test_set_phase(self) -> None:
        """set_phase updates the progress bar description."""
        p = SimulationProgress(0.0, 10.0, solver_name="CVODE", disable=True)
        p.set_phase("CVODE init")
        p.set_phase("CVODE")
        p.update(5.0)
        p.finish()

    def test_set_phase_disabled(self) -> None:
        """set_phase is a no-op when disabled."""
        p = SimulationProgress(0.0, 10.0, disable=True)
        p.set_phase("init")
        p.finish()


class TestProgressIntegration:
    """Integration tests: progress bar with actual solver calls."""

    def test_cvode_with_progress(self) -> None:
        """CVODE stepwise solver works with progress bar."""
        from tidal.solver.cvode import solve_cvode

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        p = SimulationProgress(0.0, 1.0, solver_name="CVODE", disable=True)
        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0),
            num_snapshots=5, progress=p,
        )
        assert result["success"]

    def test_leapfrog_with_progress(self) -> None:
        """Leapfrog solver works with progress bar."""
        from tidal.solver.leapfrog import solve_leapfrog

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        p = SimulationProgress(0.0, 1.0, solver_name="leapfrog", disable=True)
        result = solve_leapfrog(
            spec, grid, y0, t_span=(0.0, 1.0), dt=0.01,
            progress=p,
        )
        assert result["success"]

    def test_scipy_with_progress(self) -> None:
        """Scipy solver works with progress bar."""
        from tidal.solver.scipy_solver import solve_scipy

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        p = SimulationProgress(0.0, 1.0, solver_name="scipy", disable=True)
        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 1.0),
            num_snapshots=5, progress=p,
        )
        assert result["success"]

    def test_none_progress_unchanged_behavior(self) -> None:
        """progress=None (default) doesn't change solver behavior."""
        from tidal.solver.cvode import solve_cvode

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        # progress=None is the default — should work exactly as before
        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), num_snapshots=5,
        )
        assert result["success"]

    def test_cvode_stepwise_matches_solve(self) -> None:
        """CVODE stepwise results should match standard solve results."""
        from tidal.solver.cvode import solve_cvode

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Use a Gaussian pulse so there's actual dynamics
        y0 = np.zeros(layout.total_size)
        x = grid.coord_arrays()[0]
        y0[:grid.num_points] = np.exp(-(x - np.pi) ** 2)

        # Without progress (standard .solve path)
        result_std = solve_cvode(
            spec, grid, y0.copy(), t_span=(0.0, 1.0), num_snapshots=11,
        )

        # With progress (stepwise path)
        p = SimulationProgress(0.0, 1.0, disable=True)
        result_step = solve_cvode(
            spec, grid, y0.copy(), t_span=(0.0, 1.0), num_snapshots=11,
            progress=p,
        )

        assert result_std["success"]
        assert result_step["success"]
        np.testing.assert_allclose(
            result_std["y"][-1], result_step["y"][-1], rtol=1e-6,
        )
