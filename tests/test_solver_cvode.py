"""Tests for tidal.solver.cvode — SUNDIALS CVODE adaptive ODE solver."""

from __future__ import annotations

from typing import Any

import numpy as np

from tidal.solver.cvode import solve_cvode
from tidal.solver.grid import GridInfo
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


class TestCVODEBasic:
    def test_zero_state_stays_zero(self) -> None:
        """Zero state should remain zero (trivial equilibrium)."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            num_snapshots=5,
        )
        assert result["success"]
        np.testing.assert_allclose(result["y"][-1], 0, atol=1e-14)

    def test_standing_wave_conserved(self) -> None:
        """Standing wave sin(x)cos(t) should be preserved to tolerance."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        x = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        y0 = np.zeros(layout.total_size)
        # phi_0 = sin(x), pi_phi_0 = 0 (cos(0) = 1 for velocity)
        y0[:64] = np.sin(x)

        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 2.0), bc="periodic",
            rtol=1e-8, atol=1e-10, num_snapshots=11,
        )
        assert result["success"]

        # At t=2, analytical: phi = sin(x)*cos(2), pi = -sin(x)*sin(2)
        # Spatial discretization error (O(dx²)) dominates — use loose tolerance
        expected_phi = np.sin(x) * np.cos(2.0)
        actual_phi = result["y"][-1][:64]
        np.testing.assert_allclose(actual_phi, expected_phi, rtol=5e-3, atol=1e-3)

    def test_result_dict_keys(self) -> None:
        """Result must have t, y, success, message keys."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 0.1), bc="periodic",
            num_snapshots=3,
        )
        assert set(result.keys()) == {"t", "y", "success", "message"}

    def test_snapshot_count(self) -> None:
        """Output time points should match num_snapshots."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            num_snapshots=21,
        )
        assert result["success"]
        assert len(result["t"]) == 21

    def test_bdf_and_adams_methods(self) -> None:
        """Both BDF and Adams methods should complete successfully."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        for method in ("BDF", "Adams"):
            result = solve_cvode(
                spec, grid, y0, t_span=(0.0, 0.5), bc="periodic",
                method=method, num_snapshots=5,
            )
            assert result["success"], f"Method {method} failed"

    def test_tighter_tolerance_more_accurate(self) -> None:
        """Lower rtol should produce smaller error."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        x = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        y0 = np.zeros(layout.total_size)
        y0[:64] = np.sin(x)

        expected_phi = np.sin(x) * np.cos(1.0)

        # Loose tolerance
        r1 = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            rtol=1e-4, atol=1e-6, num_snapshots=5,
        )
        err1 = float(np.max(np.abs(r1["y"][-1][:64] - expected_phi)))

        # Tight tolerance
        r2 = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            rtol=1e-10, atol=1e-12, num_snapshots=5,
        )
        err2 = float(np.max(np.abs(r2["y"][-1][:64] - expected_phi)))

        assert err2 < err1, f"Tight tolerance error {err2} not less than loose {err1}"

    def test_snapshot_callback_called(self) -> None:
        """Snapshot callback should be invoked for each output time."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        callback_times: list[float] = []

        def cb(t: float, _y: np.ndarray) -> None:
            callback_times.append(t)

        result = solve_cvode(
            spec, grid, y0, t_span=(0.0, 1.0), bc="periodic",
            num_snapshots=11, snapshot_callback=cb,
        )
        assert result["success"]
        assert len(callback_times) == 11
