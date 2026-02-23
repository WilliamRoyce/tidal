"""Tests for tidal.solver.scipy_solver — scipy solve_ivp adaptive ODE solver."""

from __future__ import annotations

from typing import Any

import numpy as np

from tidal.solver.grid import GridInfo
from tidal.solver.scipy_solver import solve_scipy
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
            "field_rates": {
                "phi_0": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_phi_0"},
                ]
            },
            "kinetic_matrix": {
                "entries": [{"i": 0, "j": 0, "value": 1.0}],
                "dimension": 1,
            },
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


class TestScipyBasic:
    def test_zero_state_stays_zero(self) -> None:
        """Zero state should remain zero."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            num_snapshots=5,
        )
        assert result["success"]
        np.testing.assert_allclose(result["y"][-1], 0, atol=1e-14)

    def test_dop853_default(self) -> None:
        """Default DOP853 method should complete successfully."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 0.5), bc="periodic",
            num_snapshots=5,
        )
        assert result["success"]
        assert len(result["t"]) == 5

    def test_standing_wave_conserved(self) -> None:
        """Standing wave sin(x)cos(t) should be preserved to tolerance."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        x = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        y0 = np.zeros(layout.total_size)
        y0[:64] = np.sin(x)

        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 2.0), bc="periodic",
            rtol=1e-8, atol=1e-10, num_snapshots=11,
            max_step=0.1,
        )
        assert result["success"]

        # Spatial discretization error (O(dx²)) dominates — use loose tolerance
        expected_phi = np.sin(x) * np.cos(2.0)
        actual_phi = result["y"][-1][:64]
        np.testing.assert_allclose(actual_phi, expected_phi, rtol=5e-3, atol=1e-3)

    def test_parametric_methods(self) -> None:
        """Multiple scipy methods should all complete."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        for method in ("RK45", "DOP853"):
            result = solve_scipy(
                spec, grid, y0, t_span=(0.0, 0.2), bc="periodic",
                method=method, num_snapshots=3, max_step=0.1,
            )
            assert result["success"], f"Method {method} failed"

    def test_result_dict_keys(self) -> None:
        """Result must have t, y, success, message keys."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 0.1), bc="periodic",
            num_snapshots=3,
        )
        assert set(result.keys()) == {"t", "y", "success", "message"}

    def test_snapshot_callback_called(self) -> None:
        """Snapshot callback should be invoked for each output time."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        callback_times: list[float] = []

        def cb(t: float, _y: np.ndarray) -> None:
            callback_times.append(t)

        result = solve_scipy(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            num_snapshots=11, snapshot_callback=cb,
        )
        assert result["success"]
        assert len(callback_times) == 11
