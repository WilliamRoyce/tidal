"""Tests for tidal.solver.ida — SUNDIALS/IDA DAE solver integration."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.ida import build_residual_fn, solve_ida
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
            },
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


class TestResidualFunction:
    def test_residual_shape(self) -> None:
        """Residual function writes correct-sized output."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)
        resfn(0.0, y, yp, res)
        assert res.shape == (layout.total_size,)

    def test_residual_zero_state(self) -> None:
        """Zero state → zero residual (trivial equilibrium)."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)
        resfn(0.0, y, yp, res)
        np.testing.assert_allclose(res, 0, atol=1e-15)

    def test_residual_hamilton_first(self) -> None:
        """Hamilton's 1st: K*dq/dt - (pi - S) = 0 when dq/dt = pi."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        n = grid.num_points
        x = grid.axes_coords(0)

        # State: phi = sin(x), pi = cos(x)
        y = np.zeros(layout.total_size)
        y[0:n] = np.sin(x)  # phi
        y[n : 2 * n] = np.cos(x)  # pi

        # Correct yp for Hamilton's 1st: dq/dt = pi (K=I, S=0)
        yp = np.zeros(layout.total_size)
        yp[0:n] = np.cos(x)  # d(phi)/dt = pi = cos(x)
        yp[n : 2 * n] = 0  # d(pi)/dt = laplacian(phi) = -sin(x) — but we leave it 0

        res = np.zeros(layout.total_size)
        resfn(0.0, y, yp, res)

        # Field slot residual: K*yp - (pi - S) = cos(x) - cos(x) = 0
        np.testing.assert_allclose(res[0:n], 0, atol=1e-14)

        # Momentum slot residual: yp[pi] - laplacian(phi) = 0 - (-sin(x))
        # O(dx²) discretization error on 16-point grid
        np.testing.assert_allclose(res[n : 2 * n], np.sin(x), atol=0.02)


class TestIDAIntegration:
    @pytest.mark.slow
    def test_kg_standing_wave(self) -> None:
        """Klein-Gordon standing wave: sin(x)*cos(t) should stay coherent."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        n = grid.num_points
        x = grid.axes_coords(0)

        # IC: phi = sin(x), pi = 0 (standing wave)
        y0 = np.zeros(layout.total_size)
        y0[0:n] = np.sin(x)

        result = solve_ida(
            spec,
            grid,
            y0,
            t_span=(0.0, 1.0),
            num_snapshots=11,
            rtol=1e-6,
            atol=1e-8,
        )
        assert result["success"], f"IDA failed: {result['message']}"

        # At t=1, analytic: phi = sin(x)*cos(1) ≈ 0.5403*sin(x)
        phi_final = result["y"][-1][0:n]
        expected = np.sin(x) * np.cos(1.0)
        # Allow generous tolerance for 64-point grid + IDA defaults
        np.testing.assert_allclose(phi_final, expected, atol=0.05)
