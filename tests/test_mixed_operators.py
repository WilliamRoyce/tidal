"""Tests for mixed time-space operator support in energy measurement.

Covers: parsing, EOM-based acceleration, constraint differentiation,
and full mixed_* factor evaluation for Hamiltonian energy computation.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.measurement import SimulationData
from tidal.measurement._energy import (
    _compute_acceleration_from_eom,
    _differentiate_constraint,
    _differentiate_constraint_twice,
    _evaluate_hamiltonian_factor,
    _evaluate_mixed_factor,
    _merge_parameters,
    _parse_mixed_operator,
    _resolve_time_derivative,
)
from tidal.symbolic.json_loader import EquationSystem, is_known_operator

# ============================================================
# Helpers
# ============================================================


def _make_wave_spec() -> EquationSystem:
    """Build a minimal EquationSystem with wave + constraint equations via from_dict.

    Fields:
      - f_0: dynamical, d²_t f_0 = laplacian_x(f_0)
      - f_1: dynamical, d²_t f_1 = laplacian_x(f_1)
      - c_0: constraint, c_0 = 0.5*gradient_x(v_f_0) + 0.5*gradient_x(v_f_1)
      - c_1: constraint (trace), c_1 = identity(f_0) + identity(f_1)
    """
    data: dict = {
        "spacetime": {
            "dimension": 2,
            "metric": "minkowski",
            "signature": [-1, 1],
            "coordinates": ["t", "x"],
        },
        "fields": [
            {"name": "f_0", "type": "scalar", "index": 0},
            {"name": "f_1", "type": "scalar", "index": 1},
            {"name": "c_0", "type": "scalar", "index": 2},
            {"name": "c_1", "type": "scalar", "index": 3},
        ],
        "equations": [
            {
                "field": "f_0",
                "lhs": {"expression": "d2_t(f_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "f_0"},
                    ],
                },
            },
            {
                "field": "f_1",
                "lhs": {"expression": "d2_t(f_1)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "f_1"},
                    ],
                },
            },
            {
                "field": "c_0",
                "lhs": {"expression": "c_0", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 0.5,
                            "operator": "gradient_x",
                            "field": "v_f_0",
                        },
                        {
                            "coefficient": 0.5,
                            "operator": "gradient_x",
                            "field": "v_f_1",
                        },
                    ],
                },
            },
            {
                "field": "c_1",
                "lhs": {"expression": "c_1", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "identity", "field": "f_0"},
                        {"coefficient": 1.0, "operator": "identity", "field": "f_1"},
                    ],
                },
            },
        ],
        "metadata": {},
    }
    return EquationSystem.from_dict(data)


def _make_sim_data(
    spec: EquationSystem,
    n_grid: int = 32,
    n_snapshots: int = 5,
) -> SimulationData:
    """Build synthetic SimulationData from an EquationSystem.

    Fields are sine waves, velocities are cosines (mimicking wave evolution).
    """
    dx = 10.0 / n_grid
    x = np.linspace(dx / 2, 10.0 - dx / 2, n_grid)
    times = np.linspace(0.0, 1.0, n_snapshots)

    fields: dict[str, np.ndarray] = {}
    velocities: dict[str, np.ndarray] = {}

    for eq in spec.equations:
        # Field snapshot: sin(2*pi*x/L)
        base = np.sin(2 * np.pi * x / 10.0)
        fields[eq.field_name] = np.stack([base * (1.0 + 0.01 * t) for t in times])

        if eq.time_derivative_order >= 2:
            vel_base = np.cos(2 * np.pi * x / 10.0)
            velocities[eq.field_name] = np.stack(
                [vel_base * (1.0 + 0.01 * t) for t in times]
            )

    return SimulationData(
        times=times,
        fields=fields,
        velocities=velocities,
        grid_spacing=(dx,),
        grid_bounds=((0.0, 10.0),),
        periodic=(True,),
        spec=spec,
        parameters={},
    )


# ============================================================
# Operator validation tests
# ============================================================


class TestMixedOperatorValidation:
    """Test that mixed_* operators are recognized by json_loader."""

    @pytest.mark.parametrize(
        "name",
        [
            "mixed_1_0_0_1",
            "mixed_2_0_0_0",
            "mixed_1_1_0",
            "mixed_1_0_1",
            "mixed_3_0_0_0",
        ],
    )
    def test_valid_mixed_operators(self, name: str) -> None:
        assert is_known_operator(name)

    @pytest.mark.parametrize(
        "name",
        [
            "mixed_",
            "mixed_a_0",
            "mixed_1",
        ],
    )
    def test_invalid_mixed_operators(self, name: str) -> None:
        assert not is_known_operator(name)


# ============================================================
# Parsing tests
# ============================================================


class TestParseMixedOperator:
    def test_3d_mixed_1(self) -> None:
        assert _parse_mixed_operator("mixed_1_0_0_1") == (1, (0, 0, 1))

    def test_3d_mixed_2(self) -> None:
        assert _parse_mixed_operator("mixed_2_0_0_0") == (2, (0, 0, 0))

    def test_2d_mixed(self) -> None:
        assert _parse_mixed_operator("mixed_1_1_0") == (1, (1, 0))

    def test_higher_order(self) -> None:
        assert _parse_mixed_operator("mixed_3_0_1_0") == (3, (0, 1, 0))


# ============================================================
# Acceleration from EOM tests
# ============================================================


class TestAccelerationFromEOM:
    """Test _compute_acceleration_from_eom evaluates E-L equation RHS."""

    def test_dynamical_field_acceleration(self) -> None:
        """d²_t f_0 = laplacian_x(f_0) should equal the discrete laplacian."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        accel = _compute_acceleration_from_eom(data, "f_0", 0, params)
        assert accel is not None

        # Manual laplacian computation (3-point stencil, periodic)
        f = data.fields["f_0"][0]
        dx = data.grid_spacing[0]
        lap = (np.roll(f, -1) + np.roll(f, 1) - 2 * f) / dx**2

        np.testing.assert_allclose(accel, lap, atol=1e-12)

    def test_constraint_field_returns_none(self) -> None:
        """Constraint fields don't have d²_t equations → returns None."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        assert _compute_acceleration_from_eom(data, "c_0", 0, params) is None

    def test_unknown_field_returns_none(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        assert _compute_acceleration_from_eom(data, "nonexistent", 0, params) is None


# ============================================================
# Time derivative resolution tests
# ============================================================


class TestResolveTimeDerivative:
    def test_order_0_returns_field(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _resolve_time_derivative(data, "f_0", 0, 0, params)
        assert result is not None
        np.testing.assert_array_equal(result, data.fields["f_0"][0])

    def test_order_1_returns_velocity(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _resolve_time_derivative(data, "f_0", 0, 1, params)
        assert result is not None
        np.testing.assert_array_equal(result, data.velocities["f_0"][0])

    def test_order_2_returns_eom_rhs(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _resolve_time_derivative(data, "f_0", 0, 2, params)
        assert result is not None

        # Should equal laplacian_x(f_0)
        accel = _compute_acceleration_from_eom(data, "f_0", 0, params)
        np.testing.assert_array_equal(result, accel)

    def test_order_3_returns_none(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        assert _resolve_time_derivative(data, "f_0", 0, 3, params) is None


# ============================================================
# Constraint differentiation tests
# ============================================================


class TestConstraintDifferentiation:
    def test_differentiate_velocity_source_constraint(self) -> None:
        """d_t(c_0) where c_0 = 0.5*gradient_x(v_f_0) + 0.5*gradient_x(v_f_1).

        d_t(c_0) = 0.5*gradient_x(accel_f_0) + 0.5*gradient_x(accel_f_1).
        """
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _differentiate_constraint(data, "c_0", 0, params)
        assert result is not None

        # Compute expected: 0.5*grad_x(lap_x(f_0)) + 0.5*grad_x(lap_x(f_1))
        dx = data.grid_spacing[0]
        expected = np.zeros(data.fields["f_0"][0].shape)
        for fname in ["f_0", "f_1"]:
            f = data.fields[fname][0]
            lap = (np.roll(f, -1) + np.roll(f, 1) - 2 * f) / dx**2
            grad_lap = (np.roll(lap, -1) - np.roll(lap, 1)) / (2 * dx)
            expected += 0.5 * grad_lap

        np.testing.assert_allclose(result, expected, atol=1e-10)

    def test_differentiate_field_source_constraint_twice(self) -> None:
        """d²_t(c_1) where c_1 = f_0 + f_1.

        d²_t(c_1) = accel_f_0 + accel_f_1 = lap_x(f_0) + lap_x(f_1).
        """
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _differentiate_constraint_twice(data, "c_1", 0, params)
        assert result is not None

        # Sum of laplacians of f_0 and f_1
        dx = data.grid_spacing[0]
        expected = np.zeros(data.fields["f_0"][0].shape)
        for fname in ["f_0", "f_1"]:
            f = data.fields[fname][0]
            expected += (np.roll(f, -1) + np.roll(f, 1) - 2 * f) / dx**2

        np.testing.assert_allclose(result, expected, atol=1e-12)


# ============================================================
# Full mixed factor evaluation tests
# ============================================================


class TestEvaluateMixedFactor:
    def test_mixed_1_gradient_of_velocity(self) -> None:
        """mixed_1_1 on dynamic field f_0 → gradient_x(v_f_0)."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        # mixed_1_1 = d_t d_x → gradient_x(velocity)
        result = _evaluate_mixed_factor("f_0", "mixed_1_1", data, 0, params)
        assert result is not None

        # Spatial gradient of velocity
        v = data.velocities["f_0"][0]
        dx = data.grid_spacing[0]
        expected = (np.roll(v, -1) - np.roll(v, 1)) / (2 * dx)

        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_mixed_2_acceleration(self) -> None:
        """mixed_2_0 on dynamic field → acceleration from EOM."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _evaluate_mixed_factor("f_0", "mixed_2_0", data, 0, params)
        assert result is not None

        # Expected: laplacian_x(f_0) (from EOM)
        f = data.fields["f_0"][0]
        dx = data.grid_spacing[0]
        expected = (np.roll(f, -1) + np.roll(f, 1) - 2 * f) / dx**2

        np.testing.assert_allclose(result, expected, atol=1e-12)

    def test_mixed_1_on_constraint_with_velocity_sources(self) -> None:
        """mixed_1_1 on constraint c_0 → gradient_x(d_t(c_0))."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        # mixed_1_1 on c_0: d_t(gradient_x(c_0)) = gradient_x(d_t(c_0))
        # Commuting: this applies d_t first, then gradient_x
        result = _evaluate_mixed_factor("c_0", "mixed_1_1", data, 0, params)
        # c_0 is a constraint with velocity sources → d_t(c_0) needs acceleration
        # The result should be computable
        assert result is not None

    def test_mixed_2_on_trace_constraint(self) -> None:
        """mixed_2_0 on trace constraint c_1 = f_0 + f_1 → accel_f_0 + accel_f_1."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _evaluate_mixed_factor("c_1", "mixed_2_0", data, 0, params)
        assert result is not None

        # Sum of laplacians (from EOM accelerations)
        dx = data.grid_spacing[0]
        expected = np.zeros(data.fields["f_0"][0].shape)
        for fname in ["f_0", "f_1"]:
            f = data.fields[fname][0]
            expected += (np.roll(f, -1) + np.roll(f, 1) - 2 * f) / dx**2

        np.testing.assert_allclose(result, expected, atol=1e-12)


# ============================================================
# Integration via _evaluate_hamiltonian_factor dispatch
# ============================================================


class TestHamiltonianFactorMixedDispatch:
    """Verify mixed operators route through _evaluate_hamiltonian_factor."""

    def test_mixed_routes_to_evaluator(self) -> None:
        """mixed_1_1 through _evaluate_hamiltonian_factor should not raise."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)

        result = _evaluate_hamiltonian_factor("f_0", "mixed_1_1", data, 0)
        assert result is not None
        assert result.shape == data.fields["f_0"][0].shape

    def test_mixed_2_routes_to_evaluator(self) -> None:
        spec = _make_wave_spec()
        data = _make_sim_data(spec)

        result = _evaluate_hamiltonian_factor("f_0", "mixed_2_0", data, 0)
        assert result is not None

    def test_mixed_dimension_collapse(self) -> None:
        """mixed_1_0_0_1 on a 1+1D system collapses z-derivative to x."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        # mixed_1_0_0_1 has 3 spatial indices but grid is 1D.
        # Should collapse: total spatial order = 1, applied to axis 0.
        result = _evaluate_mixed_factor("f_0", "mixed_1_0_0_1", data, 0, params)
        assert result is not None

        # Should equal gradient_x(velocity), same as mixed_1_1
        expected = _evaluate_mixed_factor("f_0", "mixed_1_1", data, 0, params)
        assert expected is not None
        np.testing.assert_array_equal(result, expected)

    def test_mixed_2_0_0_0_dimension_collapse(self) -> None:
        """mixed_2_0_0_0 on 1+1D system works (no spatial derivatives)."""
        spec = _make_wave_spec()
        data = _make_sim_data(spec)
        params = _merge_parameters(data)

        result = _evaluate_mixed_factor("f_0", "mixed_2_0_0_0", data, 0, params)
        assert result is not None

        # Should equal pure acceleration (same as mixed_2_0)
        expected = _evaluate_mixed_factor("f_0", "mixed_2_0", data, 0, params)
        assert expected is not None
        np.testing.assert_array_equal(result, expected)
