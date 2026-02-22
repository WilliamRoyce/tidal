"""Tests for tidal.solver.leapfrog — Störmer-Verlet symplectic integrator."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.leapfrog import solve_leapfrog
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem


def _make_kg_spec() -> EquationSystem:
    """Klein-Gordon wave equation: d²φ/dt² = ∇²φ - m²φ."""
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


class TestLeapfrogBasic:
    def test_zero_state_stays_zero(self) -> None:
        """Zero state should remain zero (trivial equilibrium)."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        result = solve_leapfrog(
            spec, grid, y0, t_span=(0.0, 1.0), dt=0.01,
        )
        assert result["success"]
        np.testing.assert_allclose(result["y"][-1], 0, atol=1e-14)

    def test_standing_wave(self) -> None:
        """Standing wave sin(x)cos(t) should be preserved symplectically."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)

        # IC: phi = sin(x), pi = 0 (standing wave)
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        result = solve_leapfrog(
            spec, grid, y0,
            t_span=(0.0, 1.0),
            dt=0.01,
            snapshot_interval=1.0,
        )
        assert result["success"]

        # At t=1, analytic: phi = sin(x)*cos(1) ≈ 0.5403*sin(x)
        phi_final = result["y"][-1][0:n]
        expected = np.sin(x) * np.cos(1.0)
        np.testing.assert_allclose(phi_final, expected, atol=0.01)


class TestLeapfrogEnergy:
    def test_energy_conservation(self) -> None:
        """Leapfrog should preserve the discrete Hamiltonian."""
        from tidal.solver.operators import laplacian

        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)
        dx = grid.dx[0]

        # IC: phi = sin(x), pi = 0
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        # Collect the discrete Hamiltonian using the SAME laplacian as the
        # solver. H = ½ ∫ [π² - φ·∇²φ] dx (integration by parts on the
        # discrete 3-point stencil, NOT the central-difference gradient).
        energies: list[float] = []

        def _energy_callback(_t: float, y: np.ndarray) -> None:
            phi = y[0:n]
            pi = y[n:2 * n]
            lap_phi = laplacian(phi, grid)
            energy = 0.5 * dx * float(np.sum(pi**2 - phi * lap_phi))
            energies.append(energy)

        result = solve_leapfrog(
            spec, grid, y0,
            t_span=(0.0, 10.0),
            dt=0.01,
            snapshot_interval=0.1,
            snapshot_callback=_energy_callback,
        )
        assert result["success"]
        assert len(energies) > 10

        # Symplectic: shadow Hamiltonian preserved to O(dt²) ≈ 1e-4.
        e0 = energies[0]
        max_drift = max(abs(e - e0) / e0 for e in energies)
        assert max_drift < 1e-4, f"Energy drift {max_drift:.2e} exceeds tolerance"


class TestLeapfrogValidation:
    def test_warns_on_constraint_system(self) -> None:
        """Leapfrog should warn (not reject) systems with constraints."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [
                {"name": "A_0", "index": 0},
                {"name": "A_1", "index": 1},
            ],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"},
                        ],
                    },
                },
                {
                    "field": "A_1",
                    "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian", "field": "A_1"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        with pytest.warns(UserWarning, match="frozen"):
            result = solve_leapfrog(spec, grid, y0, t_span=(0.0, 0.1), dt=0.01)
        assert result["success"]

    def test_rejects_first_order_system(self) -> None:
        """Leapfrog should reject systems with first-order equations."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [{"name": "T_0", "index": 0}],
            "equations": [
                {
                    "field": "T_0",
                    "lhs": {"expression": "d_t(T_0)", "order": {"time": 1}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian", "field": "T_0"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        with pytest.raises(ValueError, match="first-order"):
            solve_leapfrog(spec, grid, y0, t_span=(0.0, 1.0), dt=0.01)

    def test_snapshot_callback(self) -> None:
        """Snapshot callback should be called at correct intervals."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        y0 = np.zeros(layout.total_size)
        callback_times: list[float] = []

        def _cb(t: float, _y: np.ndarray) -> None:
            callback_times.append(t)

        solve_leapfrog(
            spec, grid, y0,
            t_span=(0.0, 1.0),
            dt=0.1,
            snapshot_interval=0.5,
            snapshot_callback=_cb,
        )
        # Should have t=0.0, t≈0.5, t≈1.0 (3 snapshots)
        assert len(callback_times) >= 3
