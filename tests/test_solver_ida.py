"""Tests for tidal.solver.ida — SUNDIALS/IDA DAE solver integration."""

from __future__ import annotations

import warnings
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


def _make_coupled_spec_no_kinetic() -> EquationSystem:
    """Two coupled scalars with field_rates but NO kinetic_matrix."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
        "fields": [
            {"name": "phi_0", "index": 0},
            {"name": "chi_0", "index": 1},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            },
            {
                "field": "chi_0",
                "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "chi_0"},
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "phi_0": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_phi_0"},
                ],
                "chi_0": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_chi_0"},
                ],
            },
            # No kinetic_matrix — triggers field_rates fallback
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


def _make_coupled_spec_no_canonical() -> EquationSystem:
    """Two coupled scalars with NO canonical section at all."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
        "fields": [
            {"name": "phi_0", "index": 0},
            {"name": "chi_0", "index": 1},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            },
            {
                "field": "chi_0",
                "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "chi_0"},
                    ],
                },
            },
        ],
    }
    return EquationSystem.from_dict(data)


class TestFallbackWarnings:
    """Verify that fallback paths produce warnings instead of silent failure."""

    def test_field_rates_fallback_warns(self) -> None:
        """Path B: field_rates without kinetic_matrix should warn."""
        spec = _make_coupled_spec_no_kinetic()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(8, 8),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            resfn(0.0, y, yp, res)

        field_rates_warnings = [
            x for x in w if "field_rates fallback" in str(x.message)
        ]
        assert len(field_rates_warnings) >= 1, (
            "Expected warning about field_rates fallback"
        )

    def test_identity_k_fallback_warns_multi_field(self) -> None:
        """Path C: no kinetic_matrix, no field_rates, multi-field should warn."""
        spec = _make_coupled_spec_no_canonical()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(8, 8),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            resfn(0.0, y, yp, res)

        identity_warnings = [x for x in w if "assuming K = I" in str(x.message)]
        assert len(identity_warnings) >= 1, (
            "Expected warning about identity K assumption for multi-field system"
        )

    def test_single_field_identity_k_no_warning(self) -> None:
        """Path C with single field: K = I is always correct, no warning."""
        # Use existing KG spec but strip kinetic_matrix
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
            # No canonical section at all
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            resfn(0.0, y, yp, res)

        identity_warnings = [x for x in w if "assuming K = I" in str(x.message)]
        assert len(identity_warnings) == 0, (
            "Single-field system should not warn about K = I"
        )

    def test_fallback_warning_fires_once_per_field(self) -> None:
        """Warnings should only fire once per field, not per residual call."""
        spec = _make_coupled_spec_no_kinetic()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(8, 8),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        resfn = build_residual_fn(spec, layout, grid)

        y = np.zeros(layout.total_size)
        yp = np.zeros(layout.total_size)
        res = np.zeros(layout.total_size)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            resfn(0.0, y, yp, res)
            resfn(0.0, y, yp, res)
            resfn(0.0, y, yp, res)

        field_rates_warnings = [
            x for x in w if "field_rates fallback" in str(x.message)
        ]
        # 2 fields, each warned once = 2 warnings total (not 6)
        assert len(field_rates_warnings) == 2
