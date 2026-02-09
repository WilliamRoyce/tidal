"""Tests for generic Nth-order derivative support (Phase 12).

Verifies that the pipeline handles arbitrary derivative orders beyond
the hardcoded {0, 1, 2, 4}. Key use cases:
- 3rd-order (KdV dispersion)
- 5th-order (Kawahara equation)
- Directional higher-order derivatives in multi-D

Tests cover:
- Operator name validation (is_known_operator)
- Dynamic operator resolution in pde_builder
- Numeric accuracy of Nth-order derivatives
- JSON round-trip with generic operator names
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from tidal.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
    is_known_operator,
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    _op_nth_derivative,
    _parse_multi_axis_spec,
    create_initial_state,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """Fine 1D periodic grid for derivative accuracy tests."""
    return CartesianGrid([(0, 10)], 128, periodic=True)


@pytest.fixture
def grid_1d_coarse() -> CartesianGrid:
    """Coarse 1D periodic grid for fast tests."""
    return CartesianGrid([(0, 10)], 32, periodic=True)


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """2D periodic grid."""
    return CartesianGrid([(0, 10), (0, 10)], [32, 32], periodic=True)


@pytest.fixture
def kdv_dispersion_spec() -> EquationSystem:
    """Linearized KdV dispersion: d_t u = -d^3_x u.

    First-order in time, third-order spatial derivative.
    """
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("u",),
        equations=(
            ComponentEquation(
                field_name="u",
                field_index=0,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(
                        coefficient=-1.0,
                        operator="derivative_3_x",
                        field="u",
                    ),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test", "type": "kdv_dispersion"},
    )


@pytest.fixture
def kawahara_spec() -> EquationSystem:
    """Linearized Kawahara: d_t u = -d^3_x u + d^5_x u.

    First-order in time, with 3rd and 5th-order spatial derivatives.
    """
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("u",),
        equations=(
            ComponentEquation(
                field_name="u",
                field_index=0,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(
                        coefficient=-1.0,
                        operator="derivative_3_x",
                        field="u",
                    ),
                    OperatorTerm(
                        coefficient=1.0,
                        operator="derivative_5_x",
                        field="u",
                    ),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test", "type": "kawahara"},
    )


# ---------------------------------------------------------------------------
# Operator name validation (json_loader)
# ---------------------------------------------------------------------------


class TestIsKnownOperator:
    """Tests for is_known_operator function."""

    def test_static_operators_accepted(self) -> None:
        """All static operators are recognized."""
        static_ops = [
            "identity",
            "laplacian",
            "laplacian_x",
            "gradient_y",
            "cross_derivative_xy",
            "first_derivative_t",
            "biharmonic",
        ]
        for op in static_ops:
            assert is_known_operator(op), f"{op} should be known"

    def test_generic_single_axis_accepted(self) -> None:
        """Generic single-axis derivative patterns are accepted."""
        assert is_known_operator("derivative_3_x")
        assert is_known_operator("derivative_5_y")
        assert is_known_operator("derivative_7_z")
        assert is_known_operator("derivative_10_x")

    def test_generic_multi_axis_accepted(self) -> None:
        """Generic multi-axis derivative patterns are accepted."""
        assert is_known_operator("derivative_2x_1y")
        assert is_known_operator("derivative_3x_2z")

    def test_extended_axes_accepted(self) -> None:
        """Extended axis letters (w, v, u) for higher dimensions are accepted."""
        assert is_known_operator("derivative_3_w")
        assert is_known_operator("derivative_2_v")
        assert is_known_operator("derivative_1_u")
        assert is_known_operator("derivative_2w_1v")

    def test_unknown_operators_rejected(self) -> None:
        """Unknown operator names are rejected."""
        assert not is_known_operator("unknown")
        assert not is_known_operator("gradient_q")  # invalid axis
        assert not is_known_operator("derivative_x_3")  # wrong order
        assert not is_known_operator("derivative_3_q")  # invalid axis
        assert not is_known_operator("")
        assert not is_known_operator("derivative__x")  # missing number

    def test_operator_term_accepts_generic(self) -> None:
        """OperatorTerm.from_data accepts generic derivative names."""
        data = {
            "coefficient": -1.0,
            "operator": "derivative_3_x",
            "field": "u",
        }
        term = OperatorTerm.from_dict(data)
        assert term.operator == "derivative_3_x"
        assert term.coefficient == -1.0

    def test_operator_term_rejects_invalid(self) -> None:
        """OperatorTerm.from_data rejects truly unknown operators."""
        data = {
            "coefficient": 1.0,
            "operator": "quantum_leap",
            "field": "u",
        }
        with pytest.raises(ValueError, match="Unknown operator"):
            OperatorTerm.from_dict(data)


# ---------------------------------------------------------------------------
# Nth-order derivative operator (pde_builder)
# ---------------------------------------------------------------------------


class TestNthDerivativeOperator:
    """Tests for _op_nth_derivative factory function."""

    def test_first_derivative_matches_gradient(self, grid_1d: CartesianGrid) -> None:
        """_op_nth_derivative(0, 1) should match gradient_x."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d, data=np.sin(k * x))
        bc = "periodic"

        # Using nth derivative factory
        result_nth = _op_nth_derivative(0, 1)(field, bc)
        # Using gradient
        result_grad = field.gradient(bc=bc)[0]
        assert isinstance(result_grad, ScalarField)

        assert_allclose(result_nth.data, result_grad.data, atol=1e-10)

    def test_second_derivative_analytical(self, grid_1d: CartesianGrid) -> None:
        """_op_nth_derivative(0, 2) on sin(kx) gives -k^2 sin(kx)."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d, data=np.sin(k * x))
        bc = "periodic"

        result_nth = _op_nth_derivative(0, 2)(field, bc)
        expected = -(k**2) * np.sin(k * x)

        # Gradient-of-gradient uses a different stencil than field.laplace,
        # so we compare against the analytical result with a small tolerance.
        assert_allclose(result_nth.data, expected, rtol=0.005)

    def test_third_derivative_sin(self, grid_1d: CartesianGrid) -> None:
        """d^3/dx^3 sin(kx) = -k^3 cos(kx)."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d, data=np.sin(k * x))
        bc = "periodic"

        result = _op_nth_derivative(0, 3)(field, bc)
        expected = -(k**3) * np.cos(k * x)

        assert_allclose(result.data, expected, rtol=0.02)

    def test_fifth_derivative_sin(self, grid_1d: CartesianGrid) -> None:
        """d^5/dx^5 sin(kx) = k^5 cos(kx)."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d, data=np.sin(k * x))
        bc = "periodic"

        result = _op_nth_derivative(0, 5)(field, bc)
        expected = (k**5) * np.cos(k * x)

        # Higher-order derivatives need more tolerance on discrete grids
        assert_allclose(result.data, expected, rtol=0.15)

    def test_third_derivative_y_2d(self, grid_2d: CartesianGrid) -> None:
        """d^3/dy^3 sin(ky*y) = -ky^3 cos(ky*y) in 2D."""
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        ky = 2 * np.pi / 10
        field = ScalarField(grid_2d, data=np.sin(ky * y))
        bc = "periodic"

        result = _op_nth_derivative(1, 3)(field, bc)
        expected = -(ky**3) * np.cos(ky * y)

        assert_allclose(result.data, expected, rtol=0.1)


# ---------------------------------------------------------------------------
# Dynamic operator resolution in PDE
# ---------------------------------------------------------------------------


class TestDynamicOperatorResolution:
    """Tests for dynamic operator name resolution in PDEFromSpec._get_operator."""

    def test_derivative_3_x_resolves(self, grid_1d_coarse: CartesianGrid) -> None:
        """derivative_3_x resolves to 3rd-order derivative in x."""
        field = ScalarField(grid_1d_coarse, data=1.0)
        result = PDEFromSpec._get_operator("derivative_3_x", field, "periodic")
        # Uniform field: all derivatives = 0
        assert_allclose(result.data, 0.0, atol=1e-10)

    def test_derivative_5_x_resolves(self, grid_1d_coarse: CartesianGrid) -> None:
        """derivative_5_x resolves to 5th-order derivative in x."""
        field = ScalarField(grid_1d_coarse, data=1.0)
        result = PDEFromSpec._get_operator("derivative_5_x", field, "periodic")
        assert_allclose(result.data, 0.0, atol=1e-10)

    def test_derivative_3_y_requires_2d(self, grid_1d_coarse: CartesianGrid) -> None:
        """derivative_3_y requires at least 2D grid."""
        field = ScalarField(grid_1d_coarse, data=1.0)
        with pytest.raises(ValueError, match="requires at least 2D"):
            PDEFromSpec._get_operator("derivative_3_y", field, "periodic")

    def test_unknown_operator_raises(self, grid_1d_coarse: CartesianGrid) -> None:
        """Truly unknown operator still raises ValueError."""
        field = ScalarField(grid_1d_coarse, data=1.0)
        with pytest.raises(ValueError, match="Unknown operator"):
            PDEFromSpec._get_operator("quantum_leap", field, "periodic")


# ---------------------------------------------------------------------------
# KdV dispersion relation
# ---------------------------------------------------------------------------


class TestKdVDispersion:
    """Tests for linearized KdV equation: d_t u = -d^3_x u."""

    def test_kdv_state_size(self, kdv_dispersion_spec: EquationSystem) -> None:
        """First-order KdV: state_size = 1 (no momentum)."""
        assert kdv_dispersion_spec.state_size == 1

    def test_kdv_evolution_rate_shape(
        self, kdv_dispersion_spec: EquationSystem, grid_1d_coarse: CartesianGrid
    ) -> None:
        """Evolution rate returns 1 field."""
        pde = PDEFromSpec(kdv_dispersion_spec)
        state = create_initial_state(grid_1d_coarse, kdv_dispersion_spec)
        rates = pde.evolution_rate(state)
        assert len(rates) == 1

    def test_kdv_uniform_field(
        self, kdv_dispersion_spec: EquationSystem, grid_1d_coarse: CartesianGrid
    ) -> None:
        """Uniform field: d^3_x(const) = 0, so d_t u = 0."""
        pde = PDEFromSpec(kdv_dispersion_spec)
        state = FieldCollection([ScalarField(grid_1d_coarse, data=5.0)])
        rates = pde.evolution_rate(state)
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

    def test_kdv_sinusoidal(
        self, kdv_dispersion_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """sin(kx): d_t u = -d^3_x sin(kx) = -(-k^3 cos(kx)) = k^3 cos(kx)."""
        pde = PDEFromSpec(kdv_dispersion_spec)
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        u_data = np.sin(k * x)
        state = FieldCollection([ScalarField(grid_1d, data=u_data)])

        rates = pde.evolution_rate(state)

        # d_t u = -1 * d^3_x sin(kx) = -1 * (-k^3 cos(kx)) = k^3 cos(kx)
        expected = (k**3) * np.cos(k * x)
        assert_allclose(rates[0].data, expected, rtol=0.02)


# ---------------------------------------------------------------------------
# Kawahara equation
# ---------------------------------------------------------------------------


class TestKawahara:
    """Tests for linearized Kawahara: d_t u = -d^3_x u + d^5_x u."""

    def test_kawahara_uniform_field(
        self, kawahara_spec: EquationSystem, grid_1d_coarse: CartesianGrid
    ) -> None:
        """Uniform field: all derivatives = 0."""
        pde = PDEFromSpec(kawahara_spec)
        state = FieldCollection([ScalarField(grid_1d_coarse, data=3.0)])
        rates = pde.evolution_rate(state)
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

    def test_kawahara_sinusoidal(
        self, kawahara_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """sin(kx): d_t u = k^3 cos(kx) + k^5 cos(kx) = (k^3 + k^5) cos(kx)."""
        pde = PDEFromSpec(kawahara_spec)
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        u_data = np.sin(k * x)
        state = FieldCollection([ScalarField(grid_1d, data=u_data)])

        rates = pde.evolution_rate(state)

        # -1 * d^3_x sin(kx) + 1 * d^5_x sin(kx)
        # = -1 * (-k^3 cos(kx)) + 1 * (k^5 cos(kx))
        # = k^3 cos(kx) + k^5 cos(kx)
        expected = (k**3 + k**5) * np.cos(k * x)
        assert_allclose(rates[0].data, expected, rtol=0.15)


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------


class TestBackwardCompat:
    """Existing operators still work unchanged."""

    def test_static_operators_unchanged(self, grid_1d_coarse: CartesianGrid) -> None:
        """Static operators still resolve from registry, not dynamic."""
        x = cast("np.ndarray", grid_1d_coarse.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d_coarse, data=np.sin(k * x))
        bc = "periodic"

        # These should all resolve from the static registry (test verifies backward compatibility)
        PDEFromSpec._get_operator("laplacian", field, bc)
        PDEFromSpec._get_operator("identity", field, bc)
        PDEFromSpec._get_operator("gradient_x", field, bc)
        PDEFromSpec._get_operator("laplacian_x", field, bc)
        PDEFromSpec._get_operator("biharmonic", field, bc)


# ---------------------------------------------------------------------------
# Multi-axis derivative support (Issue #79)
# ---------------------------------------------------------------------------


class TestMultiAxisDerivatives:
    """Tests for multi-axis derivative parsing and resolution."""

    def test_parse_multi_axis_spec_2x_1y(self) -> None:
        """Parse '2x_1y' into [(0, 2), (1, 1)]."""
        result = _parse_multi_axis_spec("2x_1y")
        assert result == [(0, 2), (1, 1)]

    def test_parse_multi_axis_spec_1x_1y_1z(self) -> None:
        """Parse '1x_1y_1z' into [(0, 1), (1, 1), (2, 1)]."""
        result = _parse_multi_axis_spec("1x_1y_1z")
        assert result == [(0, 1), (1, 1), (2, 1)]

    def test_parse_multi_axis_spec_3x_2z(self) -> None:
        """Parse '3x_2z' into [(0, 3), (2, 2)]."""
        result = _parse_multi_axis_spec("3x_2z")
        assert result == [(0, 3), (2, 2)]

    def test_parse_multi_axis_spec_extended_axis(self) -> None:
        """Extended axis letter w (4th spatial dim) is accepted."""
        result = _parse_multi_axis_spec("2w")
        assert result == [(3, 2)]

    def test_parse_multi_axis_spec_invalid_axis(self) -> None:
        """Invalid axis letter raises ValueError."""
        with pytest.raises(ValueError, match="Invalid multi-axis"):
            _parse_multi_axis_spec("2q")

    def test_parse_multi_axis_spec_duplicate_axis_raises(self) -> None:
        """Duplicate axis (e.g. '2x_1x_1y') raises ValueError."""
        with pytest.raises(ValueError, match="Duplicate axis 'x'"):
            _parse_multi_axis_spec("2x_1x_1y")

    def test_derivative_2x_1y_resolves_2d(self, grid_2d: CartesianGrid) -> None:
        """derivative_2x_1y resolves on a 2D grid (uniform field -> 0)."""
        field = ScalarField(grid_2d, data=1.0)
        result = PDEFromSpec._get_operator("derivative_2x_1y", field, "periodic")
        assert_allclose(result.data, 0.0, atol=1e-10)

    def test_derivative_2x_1y_requires_2d(self, grid_1d_coarse: CartesianGrid) -> None:
        """derivative_2x_1y on a 1D grid raises ValueError."""
        field = ScalarField(grid_1d_coarse, data=1.0)
        with pytest.raises(ValueError, match="requires at least 2D"):
            PDEFromSpec._get_operator("derivative_2x_1y", field, "periodic")

    def test_derivative_2x_1y_numerical(self, grid_2d: CartesianGrid) -> None:
        """derivative_2x_1y on sin(kx*x)*sin(ky*y).

        d^2/dx^2 d/dy [sin(kx*x)*sin(ky*y)]
        = d^2/dx^2 [sin(kx*x)*ky*cos(ky*y)]
        = -kx^2 * ky * sin(kx*x) * cos(ky*y)
        """
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        kx = 2 * np.pi / 10
        ky = 2 * np.pi / 10

        field = ScalarField(grid_2d, data=np.sin(kx * x) * np.sin(ky * y))
        result = PDEFromSpec._get_operator("derivative_2x_1y", field, "periodic")
        expected = -(kx**2) * ky * np.sin(kx * x) * np.cos(ky * y)
        assert_allclose(result.data, expected, rtol=0.1)


# ---------------------------------------------------------------------------
# Biharmonic operator tests (Plan item #4)
# ---------------------------------------------------------------------------


class TestBiharmonicOperator:
    """Tests for the biharmonic operator ∇⁴ = ∇²(∇²)."""

    def test_biharmonic_analytical_1d(self, grid_1d: CartesianGrid) -> None:
        """∇⁴ sin(kx) = k⁴ sin(kx) in 1D."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10
        field = ScalarField(grid_1d, data=np.sin(k * x))
        bc = "periodic"

        result = PDEFromSpec._get_operator("biharmonic", field, bc)
        expected = (k**4) * np.sin(k * x)
        assert_allclose(result.data, expected, rtol=0.01)

    def test_biharmonic_uniform_field(self, grid_1d: CartesianGrid) -> None:
        """∇⁴(constant) = 0."""
        field = ScalarField(grid_1d, data=5.0)
        result = PDEFromSpec._get_operator("biharmonic", field, "periodic")
        assert_allclose(result.data, 0.0, atol=1e-10)

    def test_biharmonic_2d(self, grid_2d: CartesianGrid) -> None:
        """∇⁴ sin(kx)*sin(ky) = (kx²+ky²)² sin(kx)*sin(ky) in 2D."""
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        kx = 2 * np.pi / 10
        ky = 2 * np.pi / 10
        field = ScalarField(grid_2d, data=np.sin(kx * x) * np.sin(ky * y))
        bc = "periodic"

        result = PDEFromSpec._get_operator("biharmonic", field, bc)
        # ∇²f = -(kx²+ky²)f, so ∇⁴f = (kx²+ky²)²f
        expected = (kx**2 + ky**2) ** 2 * np.sin(kx * x) * np.sin(ky * y)
        assert_allclose(result.data, expected, rtol=0.02)

    def test_biharmonic_evolution_beam_equation(self, grid_1d: CartesianGrid) -> None:
        """Evolve d²u/dt² = -∇⁴u (Euler-Bernoulli beam) for one timestep.

        For u = sin(kx), RHS = -k⁴ sin(kx).
        After one timestep, momentum should decrease (negative acceleration).
        """
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("u",),
            equations=(
                ComponentEquation(
                    field_name="u",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(-1.0, "biharmonic", "u"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={"description": "beam equation"},
        )
        pde = PDEFromSpec(spec)
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        k = 2 * np.pi / 10

        # Initial: u = sin(kx), pi = 0
        u_field = ScalarField(grid_1d, data=np.sin(k * x))
        pi_field = ScalarField(grid_1d, data=0.0)
        state = FieldCollection([u_field, pi_field])

        rates = pde.evolution_rate(state)
        # d_t u = pi = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)
        # d_t pi = -∇⁴u = -k⁴ sin(kx)
        expected_accel = -(k**4) * np.sin(k * x)
        assert_allclose(rates[1].data, expected_accel, rtol=0.01)
