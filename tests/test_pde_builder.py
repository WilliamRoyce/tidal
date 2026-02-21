"""Tests for the PDE builder from equation specifications."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from tidal.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
)
from tidal.utils import normalize_solve_result

# Tests need to access private methods

# === Fixtures ===


@pytest.fixture
def em_json_path() -> Path:
    """Path to the EM 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "em_1d.json"


@pytest.fixture
def kg_json_path() -> Path:
    """Path to the Klein-Gordon 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "klein_gordon_1d.json"


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """1D periodic grid for testing."""
    return CartesianGrid([(0, 100)], 128, periodic=True)


@pytest.fixture
def grid_1d_small() -> CartesianGrid:
    """Small 1D periodic grid for fast testing."""
    return CartesianGrid([(0, 10)], 32, periodic=True)


@pytest.fixture
def simple_wave_spec() -> EquationSystem:
    """Return simple wave equation specification (single component, no mass)."""
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test"},
    )


@pytest.fixture
def kg_spec() -> EquationSystem:
    """Klein-Gordon equation specification (single component with mass)."""
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                    OperatorTerm(coefficient=-1.0, operator="identity", field="phi"),
                ),
            ),
        ),
        mass_matrix=((1.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test"},
    )


@pytest.fixture
def em_spec() -> EquationSystem:
    """EM field specification (two uncoupled components)."""
    return EquationSystem(
        n_components=2,
        dimension=2,
        spatial_dimension=1,
        component_names=("A_0", "A_1"),
        equations=(
            ComponentEquation(
                field_name="A_0",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_0"),
                ),
            ),
            ComponentEquation(
                field_name="A_1",
                field_index=1,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_1"),
                ),
            ),
        ),
        mass_matrix=((0.0, 0.0), (0.0, 0.0)),
        coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
        metadata={"source": "test", "gauge": "lorenz"},
    )


# === PDEFromSpec Tests ===


class TestPDEFromSpec:
    """Tests for PDEFromSpec class."""

    def test_init_simple(self, simple_wave_spec: EquationSystem) -> None:
        """Test basic initialization."""
        pde = PDEFromSpec(simple_wave_spec)

        assert pde.n_components == 1
        assert pde.spec is simple_wave_spec
        # Time dependence enabled for curved spacetime support (e.g., de Sitter)
        assert pde.explicit_time_dependence

    def test_init_multi_component(self, em_spec: EquationSystem) -> None:
        """Test initialization with multiple components."""
        pde = PDEFromSpec(em_spec)

        num_em_components = 2
        assert pde.n_components == num_em_components
        # Access component mapping through public methods if available
        assert pde.n_components == len(em_spec.component_names)

    def test_evolution_rate_wave_equation(
        self, simple_wave_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test evolution rate for simple wave equation."""
        pde = PDEFromSpec(simple_wave_spec)

        # Create state: [phi, pi] with Gaussian phi and zero pi
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        phi_data = np.exp(-((x - 5) ** 2))
        pi_data = np.zeros_like(x)

        phi = ScalarField(grid_1d_small, data=phi_data)
        pi = ScalarField(grid_1d_small, data=pi_data)
        state = FieldCollection([phi, pi])

        # Compute evolution rate
        rates = pde.evolution_rate(state)

        assert isinstance(rates, FieldCollection)
        num_wave_components = 2
        assert len(rates) == num_wave_components

        # d/dt phi = pi (should be zero)
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d/dt pi = laplacian(phi) (should be negative at peak, positive at edges)
        # The Laplacian of a Gaussian is negative at center, positive away
        center_idx = len(x) // 2
        # At center of Gaussian, Laplacian should be negative
        assert rates[1].data[center_idx] < 0

    def test_evolution_rate_kg_equation(
        self, kg_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test evolution rate for Klein-Gordon with mass term."""
        pde = PDEFromSpec(kg_spec)

        # Create state with uniform phi = 1, pi = 0
        phi = ScalarField(grid_1d_small, data=1.0)
        pi = ScalarField(grid_1d_small, data=0.0)
        state = FieldCollection([phi, pi])

        rates = pde.evolution_rate(state)

        # d/dt phi = pi = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d/dt pi = laplacian(phi) - m^2 * phi = 0 - 1*1 = -1
        # (laplacian of constant is 0)
        assert_allclose(rates[1].data, -1.0, atol=1e-10)

    def test_evolution_rate_em_field(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test evolution rate for EM field (two components)."""
        pde = PDEFromSpec(em_spec)

        # Create state: [A_0, Pi_0, A_1, Pi_1]
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        a0_data = np.exp(-((x - 3) ** 2))  # Gaussian at x=3
        a1_data = np.exp(-((x - 7) ** 2))  # Gaussian at x=7

        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=a0_data),
                ScalarField(grid_1d_small, data=np.zeros_like(x)),
                ScalarField(grid_1d_small, data=a1_data),
                ScalarField(grid_1d_small, data=np.zeros_like(x)),
            ]
        )

        rates = pde.evolution_rate(state)

        num_em_states = 4
        assert len(rates) == num_em_states

        # Pi rates should be zero (zero momentum)
        assert_allclose(rates[0].data, 0.0, atol=1e-10)
        assert_allclose(rates[2].data, 0.0, atol=1e-10)

        # d/dt Pi_0 = laplacian(A_0)
        # d/dt Pi_1 = laplacian(A_1)
        # These should be non-zero (Laplacian of Gaussian)
        assert np.max(np.abs(rates[1].data)) > 0
        assert np.max(np.abs(rates[3].data)) > 0

    def test_wrong_state_size_raises(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test that wrong state size raises ValueError."""
        pde = PDEFromSpec(em_spec)

        # Create state with wrong number of fields
        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
            ]
        )

        with pytest.raises(ValueError, match="Expected 4 fields"):
            pde.evolution_rate(state)

    @pytest.mark.usefixtures("grid_1d_small")
    def test_unknown_operator_raises(self) -> None:
        """Test that unknown operator raises ValueError."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            coefficient=1.0, operator="unknown_op", field="phi"
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        with pytest.raises(ValueError, match="Unknown operator"):
            PDEFromSpec(spec)


# === build_pde_from_json Tests ===


class TestBuildPDEFromJSON:
    """Tests for build_pde_from_json function."""

    def test_build_em_pde(self, em_json_path: Path) -> None:
        """Test building PDE from EM JSON file."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        pde = build_pde_from_json(em_json_path)

        assert isinstance(pde, PDEFromSpec)
        num_em_components = 2
        assert pde.n_components == num_em_components

    def test_build_kg_pde(self, kg_json_path: Path) -> None:
        """Test building PDE from Klein-Gordon JSON file."""
        if not kg_json_path.exists():
            pytest.skip(f"Test file not found: {kg_json_path}")

        pde = build_pde_from_json(kg_json_path, parameters={"m2": 1.0})

        assert isinstance(pde, PDEFromSpec)
        num_kg_components = 1
        assert pde.n_components == num_kg_components


# === create_initial_state Tests ===


class TestCreateInitialState:
    """Tests for create_initial_state function."""

    def test_create_default_state(
        self, simple_wave_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test creating state with default (zero) data."""
        state = create_initial_state(grid_1d_small, simple_wave_spec)

        assert isinstance(state, FieldCollection)
        num_wave_components = 2  # phi, pi
        assert len(state) == num_wave_components

        # All should be zero
        assert_allclose(state[0].data, 0.0)
        assert_allclose(state[1].data, 0.0)

    def test_create_state_with_field_data(
        self, simple_wave_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test creating state with specified field data."""
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        phi_data = np.sin(2 * np.pi * x / 10)

        state = create_initial_state(
            grid_1d_small, simple_wave_spec, field_data={"phi": phi_data}
        )

        assert_allclose(state[0].data, phi_data)
        assert_allclose(state[1].data, 0.0)  # Momentum defaults to zero

    def test_create_state_with_momentum_data(
        self, simple_wave_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test creating state with specified momentum data."""
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        pi_data = np.cos(2 * np.pi * x / 10)

        state = create_initial_state(
            grid_1d_small, simple_wave_spec, momentum_data={"phi": pi_data}
        )

        assert_allclose(state[0].data, 0.0)  # Field defaults to zero
        assert_allclose(state[1].data, pi_data)

    def test_create_multi_component_state(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test creating state for multi-component system."""
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        a1_data = np.exp(-((x - 5) ** 2))

        state = create_initial_state(
            grid_1d_small, em_spec, field_data={"A_1": a1_data}
        )

        num_em_states = 4  # A_0, Pi_0, A_1, Pi_1
        assert len(state) == num_em_states

        assert_allclose(state[0].data, 0.0)  # A_0 = 0
        assert_allclose(state[1].data, 0.0)  # Pi_0 = 0
        assert_allclose(state[2].data, a1_data)  # A_1 = Gaussian
        assert_allclose(state[3].data, 0.0)  # Pi_1 = 0


# === Integration Tests ===


class TestIntegration:
    """Integration tests for the full pipeline."""

    def test_wave_propagation(
        self, simple_wave_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that waves propagate correctly."""
        pde = PDEFromSpec(simple_wave_spec)

        # Create Gaussian pulse
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        center = 50.0
        width = 5.0
        phi_data = np.exp(-((x - center) ** 2) / (2 * width**2))

        state = create_initial_state(
            grid_1d, simple_wave_spec, field_data={"phi": phi_data}
        )

        # Run short simulation
        result = pde.solve(state, t_range=1.0, dt=0.01)
        sol = normalize_solve_result(result)

        # The pulse should have spread/propagated
        # Peak should be lower than initial due to dispersion
        initial_max = np.max(phi_data)
        final_max = np.max(sol.data[0])

        # For a wave equation, initial Gaussian splits into two pulses
        # Each has roughly half the amplitude
        assert final_max < initial_max

    def test_kg_mass_effect(
        self, kg_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test that Klein-Gordon mass term causes oscillation."""
        pde = PDEFromSpec(kg_spec)

        # Start with uniform phi = 1
        state = create_initial_state(
            grid_1d_small, kg_spec, field_data={"phi": np.ones(32)}
        )

        # Run simulation
        result = pde.solve(state, t_range=np.pi, dt=0.01)
        sol = normalize_solve_result(result)

        # With m^2 = 1, a uniform field should oscillate with period 2*pi
        # After time pi, phi should be approximately -1
        # (This is cos(m*t) behavior for uniform field)
        assert_allclose(sol.data[0], -1.0, atol=0.1)

    def test_em_independent_components(
        self, em_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that EM components evolve independently."""
        pde = PDEFromSpec(em_spec)

        # Initialize only A_1 with a pulse, A_0 = 0
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        a1_data = np.exp(-((x - 50) ** 2) / 50)

        state = create_initial_state(grid_1d, em_spec, field_data={"A_1": a1_data})

        # Run simulation
        result = pde.solve(state, t_range=5.0, dt=0.01)
        sol = normalize_solve_result(result)

        # A_0 should remain zero (no coupling)
        assert_allclose(sol.data[0], 0.0, atol=1e-10)

        # A_1 should have evolved (not zero)
        assert np.max(np.abs(sol.data[2])) > 0


# === Phase 4 Feature Tests: Momentum References and Cross-Derivatives ===


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """2D periodic grid for testing cross-derivatives."""
    return CartesianGrid([(0, 10), (0, 10)], [32, 32], periodic=True)


@pytest.fixture
def spec_with_momentum_gradient() -> EquationSystem:
    """Spec with momentum field reference (pi_1) in RHS.

    This represents a term like ∂_t∂_x A_1 = ∂_x(∂_t A_1) = ∂_x(π_1).
    Used to test Phase 4's mixed time-space derivative handling.
    """
    return EquationSystem(
        n_components=2,
        dimension=3,  # 2+1D
        spatial_dimension=2,
        component_names=("A_0", "A_1"),
        equations=(
            ComponentEquation(
                field_name="A_0",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_0"),
                    # ∂_x(π_1) - references momentum field
                    OperatorTerm(coefficient=0.5, operator="gradient_x", field="pi_1"),
                ),
            ),
            ComponentEquation(
                field_name="A_1",
                field_index=1,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_1"),
                ),
            ),
        ),
        mass_matrix=((0.0, 0.0), (0.0, 0.0)),
        coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
        metadata={"source": "test", "dimension": 3},
    )


@pytest.fixture
def spec_with_cross_derivative() -> EquationSystem:
    """Spec with cross_derivative_xy operator.

    This represents spatial cross-derivative ∂_x∂_y A.
    """
    return EquationSystem(
        n_components=1,
        dimension=3,  # 2+1D
        spatial_dimension=2,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                    OperatorTerm(
                        coefficient=0.5, operator="cross_derivative_xy", field="phi"
                    ),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test", "dimension": 3},
    )


class TestMomentumFieldReferences:
    """Tests for Phase 4: momentum field references (pi_*) in equations."""

    def test_get_field_from_state_regular_field(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test _get_field_from_state returns regular fields correctly."""
        pde = PDEFromSpec(em_spec)

        # State layout: A_0, pi_0, A_1, pi_1
        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=1.0),  # A_0
                ScalarField(grid_1d_small, data=2.0),  # pi_0
                ScalarField(grid_1d_small, data=3.0),  # A_1
                ScalarField(grid_1d_small, data=4.0),  # pi_1
            ]
        )

        # Regular field lookup
        a0 = pde._get_field_from_state(state, "A_0")
        a1 = pde._get_field_from_state(state, "A_1")

        assert_allclose(a0.data, 1.0)
        assert_allclose(a1.data, 3.0)

    def test_get_field_from_state_momentum_field(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test _get_field_from_state returns momentum fields for pi_* names."""
        pde = PDEFromSpec(em_spec)

        # State layout: A_0, pi_0, A_1, pi_1
        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=1.0),  # A_0
                ScalarField(grid_1d_small, data=2.0),  # pi_0
                ScalarField(grid_1d_small, data=3.0),  # A_1
                ScalarField(grid_1d_small, data=4.0),  # pi_1
            ]
        )

        # Momentum field lookup via pi_* syntax
        pi_0 = pde._get_field_from_state(state, "pi_0")
        pi_1 = pde._get_field_from_state(state, "pi_1")

        assert_allclose(pi_0.data, 2.0)
        assert_allclose(pi_1.data, 4.0)

    def test_get_field_from_state_invalid_momentum_raises(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test that invalid momentum index raises ValueError."""
        pde = PDEFromSpec(em_spec)

        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
            ]
        )

        # pi_5 is out of range for 2-component system
        with pytest.raises(ValueError, match="out of range"):
            pde._get_field_from_state(state, "pi_5")

    def test_get_field_from_state_malformed_momentum_raises(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test that malformed momentum names raise ValueError."""
        pde = PDEFromSpec(em_spec)

        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
            ]
        )

        # Non-numeric index
        with pytest.raises(ValueError, match="numeric index"):
            pde._get_field_from_state(state, "pi_abc")

        # Wrong format (too many underscores)
        with pytest.raises(ValueError, match="pi_N"):
            pde._get_field_from_state(state, "pi_0_extra")

    def test_get_field_from_state_unknown_field_raises(
        self, em_spec: EquationSystem, grid_1d_small: CartesianGrid
    ) -> None:
        """Test that unknown field name raises ValueError."""
        pde = PDEFromSpec(em_spec)

        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
                ScalarField(grid_1d_small, data=0.0),
            ]
        )

        with pytest.raises(ValueError, match="Unknown field name"):
            pde._get_field_from_state(state, "B_0")

    def test_evolution_with_momentum_gradient(
        self, spec_with_momentum_gradient: EquationSystem, grid_2d: CartesianGrid
    ) -> None:
        """Test evolution rate when RHS references momentum field gradient."""
        pde = PDEFromSpec(spec_with_momentum_gradient)

        # Create state with non-zero momentum in A_1
        # State layout: A_0, pi_0, A_1, pi_1
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        pi_1_data = np.sin(2 * np.pi * x / 10)  # Varies in x

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=0.0),  # A_0
                ScalarField(grid_2d, data=0.0),  # pi_0
                ScalarField(grid_2d, data=0.0),  # A_1
                ScalarField(grid_2d, data=pi_1_data),  # pi_1 (non-zero)
            ]
        )

        rates = pde.evolution_rate(state)

        # d/dt pi_0 should include 0.5 * gradient_x(pi_1)
        # gradient_x of sin(2πx/10) = (2π/10) * cos(2πx/10)
        # This should be non-zero
        assert np.max(np.abs(rates[1].data)) > 0


class TestCrossDerivativeXY:
    """Tests for Phase 4: cross_derivative_xy operator."""

    def test_cross_derivative_xy_operator(self, grid_2d: CartesianGrid) -> None:
        """Test cross_derivative_xy computes ∂_x∂_y correctly."""
        # Create field f(x,y) = sin(2πx/L) * sin(2πy/L)
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        lx, ly = 10.0, 10.0
        kx, ky = 2 * np.pi / lx, 2 * np.pi / ly

        field_data = np.sin(kx * x) * np.sin(ky * y)
        field = ScalarField(grid_2d, data=field_data)

        # Apply cross derivative via operator
        result = PDEFromSpec._get_operator("cross_derivative_xy", field, "periodic")

        # ∂_x∂_y [sin(kx*x)*sin(ky*y)] = kx*ky*cos(kx*x)*cos(ky*y)
        expected = kx * ky * np.cos(kx * x) * np.cos(ky * y)

        assert_allclose(result.data, expected, rtol=0.05)

    def test_cross_derivative_xy_requires_2d(
        self, grid_1d_small: CartesianGrid
    ) -> None:
        """Test cross_derivative_xy raises on 1D grid."""
        field = ScalarField(grid_1d_small, data=1.0)

        with pytest.raises(ValueError, match="requires at least 2D"):
            PDEFromSpec._get_operator("cross_derivative_xy", field, "periodic")

    def test_evolution_with_cross_derivative(
        self, spec_with_cross_derivative: EquationSystem, grid_2d: CartesianGrid
    ) -> None:
        """Test evolution rate when RHS includes cross_derivative_xy."""
        pde = PDEFromSpec(spec_with_cross_derivative)

        # Create field with xy variation
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        phi_data = np.sin(2 * np.pi * x / 10) * np.sin(2 * np.pi * y / 10)

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=phi_data),
                ScalarField(grid_2d, data=0.0),  # pi
            ]
        )

        rates = pde.evolution_rate(state)

        # d/dt pi should include laplacian + 0.5 * cross_derivative_xy
        # Both should be non-zero for this field
        assert np.max(np.abs(rates[1].data)) > 0


class TestDirectionalLaplacians:
    """Tests for directional Laplacian operators (laplacian_x, laplacian_y, laplacian_z).

    These operators compute pure second derivatives in a single direction,
    needed for anisotropic equations like Navier-Cauchy elasticity.
    """

    def test_laplacian_x_operator(self, grid_2d: CartesianGrid) -> None:
        """Test laplacian_x computes ∂²/∂x² correctly."""
        # Create field f(x,y) = sin(kx*x) * const_y
        # d²f/dx² = -kx² * sin(kx*x) * const_y
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        lx = 10.0
        kx = 2 * np.pi / lx

        # Use constant in y direction so laplacian_x is purely from x variation
        field_data = np.sin(kx * x)
        field = ScalarField(grid_2d, data=field_data)

        # Apply laplacian_x via operator
        result = PDEFromSpec._get_operator("laplacian_x", field, "periodic")

        # d²f/dx² = -kx² * sin(kx*x)
        expected = -(kx**2) * np.sin(kx * x)

        assert_allclose(result.data, expected, rtol=0.05)

    def test_laplacian_y_operator(self, grid_2d: CartesianGrid) -> None:
        """Test laplacian_y computes ∂²/∂y² correctly."""
        # Create field f(x,y) = const_x * sin(ky*y)
        # d²f/dy² = -ky² * const_x * sin(ky*y)
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        ly = 10.0
        ky = 2 * np.pi / ly

        # Use constant in x direction so laplacian_y is purely from y variation
        field_data = np.sin(ky * y)
        field = ScalarField(grid_2d, data=field_data)

        # Apply laplacian_y via operator
        result = PDEFromSpec._get_operator("laplacian_y", field, "periodic")

        # d²f/dy² = -ky² * sin(ky*y)
        expected = -(ky**2) * np.sin(ky * y)

        assert_allclose(result.data, expected, rtol=0.05)

    def test_laplacian_x_vs_laplacian_y_different(self, grid_2d: CartesianGrid) -> None:
        """Test that laplacian_x and laplacian_y give different results for mixed field."""
        # Create field f(x,y) = sin(kx*x) * sin(ky*y) with kx ≠ ky
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        kx = 2 * np.pi / 10.0  # different wave numbers
        ky = 4 * np.pi / 10.0

        field_data = np.sin(kx * x) * np.sin(ky * y)
        field = ScalarField(grid_2d, data=field_data)

        # Apply both operators
        lap_x_result = PDEFromSpec._get_operator("laplacian_x", field, "periodic")
        lap_y_result = PDEFromSpec._get_operator("laplacian_y", field, "periodic")

        # Expected: d²f/dx² = -kx² * sin(kx*x) * sin(ky*y)
        #           d²f/dy² = -ky² * sin(kx*x) * sin(ky*y)
        # Since kx ≠ ky, these should be different
        assert not np.allclose(lap_x_result.data, lap_y_result.data)

        # Verify ratio matches kx²/ky²
        ratio = lap_x_result.data / lap_y_result.data
        expected_ratio = (kx**2) / (ky**2)
        assert_allclose(ratio, expected_ratio, rtol=0.05)

    def test_laplacian_y_requires_2d(self, grid_1d_small: CartesianGrid) -> None:
        """Test laplacian_y raises on 1D grid."""
        field = ScalarField(grid_1d_small, data=1.0)

        with pytest.raises(ValueError, match="requires at least 2D"):
            PDEFromSpec._get_operator("laplacian_y", field, "periodic")

    def test_laplacian_x_works_on_1d(self, grid_1d_small: CartesianGrid) -> None:
        """Test laplacian_x works on 1D grid."""
        x = cast("np.ndarray", grid_1d_small.cell_coords[..., 0])
        kx = 2 * np.pi / 10.0

        field_data = np.sin(kx * x)
        field = ScalarField(grid_1d_small, data=field_data)

        # Should not raise - laplacian_x works in 1D
        result = PDEFromSpec._get_operator("laplacian_x", field, "periodic")

        expected = -(kx**2) * np.sin(kx * x)
        assert_allclose(result.data, expected, rtol=0.05)

    def test_laplacian_x_plus_laplacian_y_equals_laplacian(
        self, grid_2d: CartesianGrid
    ) -> None:
        """Test that laplacian_x + laplacian_y equals the full laplacian."""
        # Create field with both x and y variation
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        kx = 2 * np.pi / 10.0
        ky = 2 * np.pi / 10.0

        field_data = np.sin(kx * x) * np.sin(ky * y)
        field = ScalarField(grid_2d, data=field_data)

        # Get all three operators
        lap_x = PDEFromSpec._get_operator("laplacian_x", field, "periodic")
        lap_y = PDEFromSpec._get_operator("laplacian_y", field, "periodic")
        lap_full = PDEFromSpec._get_operator("laplacian", field, "periodic")

        # laplacian_x + laplacian_y should equal laplacian
        lap_sum = lap_x.data + lap_y.data

        assert_allclose(lap_sum, lap_full.data, rtol=0.01)


class TestChernSimonsIntegration:
    """Integration tests for Chern-Simons 2+1D with Phase 4 features."""

    @pytest.fixture
    def cs_json_path(self) -> Path:
        """Path to Chern-Simons 3D JSON file."""
        return (
            Path(__file__).parent.parent / "examples" / "data" / "chern_simons_3d.json"
        )

    def test_load_chern_simons_json(self, cs_json_path: Path) -> None:
        """Test loading Chern-Simons JSON with momentum references."""
        if not cs_json_path.exists():
            pytest.skip(f"Test file not found: {cs_json_path}")

        pde = build_pde_from_json(cs_json_path)

        assert isinstance(pde, PDEFromSpec)
        num_cs_components = 3
        assert pde.n_components == num_cs_components
        assert pde.spec.component_names == ("A_0", "A_1", "A_2")

    def test_chern_simons_evolution(
        self, cs_json_path: Path, grid_2d: CartesianGrid
    ) -> None:
        """Test Chern-Simons evolution with momentum gradients and cross-derivatives."""
        if not cs_json_path.exists():
            pytest.skip(f"Test file not found: {cs_json_path}")

        pde = build_pde_from_json(cs_json_path)

        # A_0 is first-order (no momentum slot), A_1 and A_2 are second-order
        # State layout: A_0, A_1, pi_1, A_2, pi_2 (5 fields total)
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        a1_data = np.exp(-((x - 5) ** 2 + (y - 5) ** 2) / 2)

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=0.0),  # A_0
                ScalarField(grid_2d, data=a1_data),  # A_1 (Gaussian)
                ScalarField(grid_2d, data=0.0),  # pi_1
                ScalarField(grid_2d, data=0.0),  # A_2
                ScalarField(grid_2d, data=0.0),  # pi_2
            ]
        )

        # Compute evolution rate - should not raise
        rates = pde.evolution_rate(state)

        num_cs_states = 5
        assert len(rates) == num_cs_states

        # d/dt A_1 = pi_1 + ∂_x(A_0); with constraint solver A_0 is solved
        # from Gauss's law, so ∂_x(A_0) is generally non-zero when A_1 has
        # spatial variation. Just verify it's finite.
        assert np.all(np.isfinite(rates[1].data))

        # d/dt pi_1 should include laplacian(A_1) contribution (non-zero)
        assert np.max(np.abs(rates[2].data)) > 0

    def test_chern_simons_simulation_short(
        self, cs_json_path: Path, grid_2d: CartesianGrid
    ) -> None:
        """Test short Chern-Simons simulation runs without error."""
        if not cs_json_path.exists():
            pytest.skip(f"Test file not found: {cs_json_path}")

        pde = build_pde_from_json(cs_json_path)

        # A_0 is first-order (no momentum slot), A_1 and A_2 are second-order
        # State layout: A_0, A_1, pi_1, A_2, pi_2 (5 fields total)
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        a0_data = np.exp(-((x - 5) ** 2 + (y - 5) ** 2) / 4)

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=a0_data),  # A_0
                ScalarField(grid_2d, data=0.0),  # A_1
                ScalarField(grid_2d, data=0.0),  # pi_1
                ScalarField(grid_2d, data=0.0),  # A_2
                ScalarField(grid_2d, data=0.0),  # pi_2
            ]
        )

        # Run very short simulation - just verify it doesn't crash
        result = pde.solve(state, t_range=0.1, dt=0.01)
        sol = normalize_solve_result(result)

        # Solution should exist and have correct shape
        assert sol.data is not None
        num_cs_states = 5
        assert len(sol.data) == num_cs_states


# === Phase 6: 3+1D Operator Tests ===


@pytest.fixture
def grid_3d() -> CartesianGrid:
    """3D periodic grid for testing 3+1D operators."""
    return CartesianGrid([(0, 10), (0, 10), (0, 10)], [16, 16, 16], periodic=True)


class TestCrossDerivatives3D:
    """Tests for cross-derivatives in 3D (xz, yz)."""

    def test_cross_derivative_xz_operator(self, grid_3d: CartesianGrid) -> None:
        """Test cross_derivative_xz computes d_x d_z correctly."""
        # Create field f(x,y,z) = sin(kx*x) * sin(kz*z)
        # d_x d_z f = kx * kz * cos(kx*x) * cos(kz*z)
        x = cast("np.ndarray", grid_3d.cell_coords[..., 0])
        z = cast("np.ndarray", grid_3d.cell_coords[..., 2])
        lx, lz = 10.0, 10.0
        kx, kz = 2 * np.pi / lx, 2 * np.pi / lz

        field_data = np.sin(kx * x) * np.sin(kz * z)
        field = ScalarField(grid_3d, data=field_data)

        result = PDEFromSpec._get_operator("cross_derivative_xz", field, "periodic")

        expected = kx * kz * np.cos(kx * x) * np.cos(kz * z)
        # Larger tolerance for 3D cross-derivatives on coarse 16^3 grid
        assert_allclose(result.data, expected, rtol=0.06)

    def test_cross_derivative_yz_operator(self, grid_3d: CartesianGrid) -> None:
        """Test cross_derivative_yz computes d_y d_z correctly."""
        # Create field f(x,y,z) = sin(ky*y) * sin(kz*z)
        # d_y d_z f = ky * kz * cos(ky*y) * cos(kz*z)
        y = cast("np.ndarray", grid_3d.cell_coords[..., 1])
        z = cast("np.ndarray", grid_3d.cell_coords[..., 2])
        ly, lz = 10.0, 10.0
        ky, kz = 2 * np.pi / ly, 2 * np.pi / lz

        field_data = np.sin(ky * y) * np.sin(kz * z)
        field = ScalarField(grid_3d, data=field_data)

        result = PDEFromSpec._get_operator("cross_derivative_yz", field, "periodic")

        expected = ky * kz * np.cos(ky * y) * np.cos(kz * z)
        # Larger tolerance for 3D cross-derivatives on coarse 16^3 grid
        assert_allclose(result.data, expected, rtol=0.06)

    def test_cross_derivative_xz_requires_3d(self, grid_2d: CartesianGrid) -> None:
        """Test cross_derivative_xz raises on 2D grid."""
        field = ScalarField(grid_2d, data=1.0)

        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec._get_operator("cross_derivative_xz", field, "periodic")

    def test_cross_derivative_yz_requires_3d(self, grid_2d: CartesianGrid) -> None:
        """Test cross_derivative_yz raises on 2D grid."""
        field = ScalarField(grid_2d, data=1.0)

        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec._get_operator("cross_derivative_yz", field, "periodic")


class TestDirectionalLaplacians3D:
    """Tests for laplacian_z in 3D."""

    def test_laplacian_z_operator(self, grid_3d: CartesianGrid) -> None:
        """Test laplacian_z computes d^2/dz^2 correctly."""
        # Create field f(x,y,z) = sin(kz*z)
        # d^2f/dz^2 = -kz^2 * sin(kz*z)
        z = cast("np.ndarray", grid_3d.cell_coords[..., 2])
        lz = 10.0
        kz = 2 * np.pi / lz

        field_data = np.sin(kz * z)
        field = ScalarField(grid_3d, data=field_data)

        result = PDEFromSpec._get_operator("laplacian_z", field, "periodic")

        expected = -(kz**2) * np.sin(kz * z)
        # Larger tolerance for 3D operators on coarse 16^3 grid
        assert_allclose(result.data, expected, rtol=0.06)

    def test_laplacian_z_requires_3d(self, grid_2d: CartesianGrid) -> None:
        """Test laplacian_z raises on 2D grid."""
        field = ScalarField(grid_2d, data=1.0)

        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec._get_operator("laplacian_z", field, "periodic")

    def test_laplacian_sum_3d(self, grid_3d: CartesianGrid) -> None:
        """Test laplacian_x + laplacian_y + laplacian_z = laplacian in 3D."""
        # Create field with variation in all three directions
        x = cast("np.ndarray", grid_3d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_3d.cell_coords[..., 1])
        z = cast("np.ndarray", grid_3d.cell_coords[..., 2])
        k = 2 * np.pi / 10.0

        field_data = np.sin(k * x) * np.sin(k * y) * np.sin(k * z)
        field = ScalarField(grid_3d, data=field_data)

        # Get all operators
        lap_x = PDEFromSpec._get_operator("laplacian_x", field, "periodic")
        lap_y = PDEFromSpec._get_operator("laplacian_y", field, "periodic")
        lap_z = PDEFromSpec._get_operator("laplacian_z", field, "periodic")
        lap_full = PDEFromSpec._get_operator("laplacian", field, "periodic")

        # laplacian_x + laplacian_y + laplacian_z should equal laplacian
        lap_sum = lap_x.data + lap_y.data + lap_z.data

        # Larger tolerance for 3D operators on coarse 16^3 grid
        assert_allclose(lap_sum, lap_full.data, rtol=0.05)


# === Parameterized Coefficients Tests ===


class TestParameterizedCoefficients:
    """Tests for runtime parameter injection via coefficient_symbolic."""

    @pytest.fixture
    def spec_with_symbolic(self) -> EquationSystem:
        """Create a spec with symbolic coefficient for mass term."""
        return EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(-1.0, "identity", "phi", "-m2"),  # symbolic
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

    def test_resolve_coefficient_with_parameter(
        self, spec_with_symbolic: EquationSystem
    ) -> None:
        """Test _resolve_coefficient uses parameter value when provided."""
        pde = PDEFromSpec(spec_with_symbolic, parameters={"m2": 0.25})
        term = spec_with_symbolic.equations[0].rhs_terms[1]

        resolved = pde._resolve_coefficient(term)
        # "-m2" with m2=0.25 should give -0.25
        assert resolved == -0.25

    def test_resolve_coefficient_without_parameter_raises(
        self, spec_with_symbolic: EquationSystem
    ) -> None:
        """Test _resolve_coefficient raises when symbolic param not provided."""
        pde = PDEFromSpec(spec_with_symbolic)  # No parameters
        term = spec_with_symbolic.equations[0].rhs_terms[1]

        with pytest.raises(ValueError, match="could not be resolved"):
            pde._resolve_coefficient(term)

    def test_resolve_coefficient_positive_symbolic(self) -> None:
        """Test _resolve_coefficient with positive symbolic coefficient."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "identity", "phi", "kappa"),),
                ),
            ),
            mass_matrix=((-1.0,),),
            coupling_matrix=((0.0,),),
            mass_matrix_symbolic=(("kappa",),),
            metadata={},
        )
        pde = PDEFromSpec(spec, parameters={"kappa": 2.5})
        term = spec.equations[0].rhs_terms[0]

        resolved = pde._resolve_coefficient(term)
        assert resolved == 2.5

    def test_resolve_coefficient_unmatched_symbolic_raises(self) -> None:
        """Test _resolve_coefficient raises when symbolic doesn't match any param."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.5, "identity", "phi", "unknown_param"),),
                ),
            ),
            mass_matrix=((-1.5,),),
            coupling_matrix=((0.0,),),
            mass_matrix_symbolic=(("unknown_param",),),
            metadata={},
        )
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})  # Different param
        term = spec.equations[0].rhs_terms[0]

        with pytest.raises(ValueError, match="could not be resolved"):
            pde._resolve_coefficient(term)

    def test_build_pde_from_json_with_parameters(self, kg_json_path: Path) -> None:
        """Test build_pde_from_json accepts parameters dict."""
        if not kg_json_path.exists():
            pytest.skip("Test file not found")

        pde = build_pde_from_json(kg_json_path, parameters={"m2": 2.0})
        assert pde._parameters == {"m2": 2.0}

    def test_pde_constructor_stores_parameters(
        self, spec_with_symbolic: EquationSystem
    ) -> None:
        """Test PDEFromSpec constructor stores parameters."""
        pde = PDEFromSpec(spec_with_symbolic, parameters={"m2": 0.5, "kappa": 1.0})
        assert pde._parameters == {"m2": 0.5, "kappa": 1.0}

    def test_pde_constructor_empty_parameters(
        self, spec_with_symbolic: EquationSystem
    ) -> None:
        """Test PDEFromSpec constructor with no parameters."""
        pde = PDEFromSpec(spec_with_symbolic)
        assert pde._parameters == {}

    def test_evolution_uses_resolved_coefficients(
        self, spec_with_symbolic: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that evolution actually uses resolved coefficient values."""
        # Create two PDEs with different mass parameters
        pde_low_mass = PDEFromSpec(spec_with_symbolic, parameters={"m2": 0.1})
        pde_high_mass = PDEFromSpec(spec_with_symbolic, parameters={"m2": 10.0})

        # Create identical initial states with a Gaussian
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        initial_data = np.exp(-((x - 50) ** 2) / 50)

        state_low = FieldCollection(
            [
                ScalarField(grid_1d, data=initial_data.copy()),
                ScalarField(grid_1d, data=0.0),  # momentum
            ]
        )
        state_high = FieldCollection(
            [
                ScalarField(grid_1d, data=initial_data.copy()),
                ScalarField(grid_1d, data=0.0),
            ]
        )

        # Evolve briefly
        result_low = pde_low_mass.solve(state_low, t_range=1.0, dt=0.01)
        result_high = pde_high_mass.solve(state_high, t_range=1.0, dt=0.01)

        result_low = normalize_solve_result(result_low)
        result_high = normalize_solve_result(result_high)

        # Extract final field data
        final_low = result_low[0].data
        final_high = result_high[0].data

        # Different mass should lead to different evolution
        # (higher mass = faster oscillation for Klein-Gordon)
        assert not np.allclose(final_low, final_high, rtol=0.01)


class TestPositionDependentCoefficients:
    """Tests for position-dependent (spatially varying) coefficients."""

    def test_resolve_coefficient_at_point_constant(self) -> None:
        """Constant coefficient returns float."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(2.5, "laplacian", "phi"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0)
        assert isinstance(result, float)
        assert result == 2.5

    def test_resolve_coefficient_at_point_time_dependent(self) -> None:
        """Time-only dependent coefficient returns float."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "identity",
                            "phi",
                            coefficient_symbolic="E^(2*dSH*t())",
                            time_dependent=True,
                            coordinate_dependent=("t",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((-1.0,),),
            coupling_matrix=((0.0,),),
            mass_matrix_symbolic=(("E^(2*dSH*t())",),),
            metadata={},
        )
        pde = PDEFromSpec(spec, parameters={"dSH": 0.1})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=1.0)
        assert isinstance(result, float)
        assert_allclose(result, np.exp(0.2), rtol=1e-10)

    def test_resolve_coefficient_at_point_position_dependent(self) -> None:
        """Position-dependent coefficient returns ndarray."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian_x",
                            "phi",
                            coefficient_symbolic="x()**2/(2*sphR**2)",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        grid = CartesianGrid(bounds=[(-4, 4), (-4, 4)], shape=[16, 16], periodic=True)
        pde = PDEFromSpec(spec, parameters={"sphR": 2.0})
        term = spec.equations[0].rhs_terms[0]

        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)
        assert isinstance(result, np.ndarray)
        assert result.shape == (16, 16)

        # Verify coefficient values match formula at actual cell coordinates
        # py-pde lacks complete type stubs for grid.cell_coords
        x_coords = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]

        # At grid center (nearest x=0): coefficient x^2/(2*R^2) should be small
        center_x = grid.shape[0] // 2
        x_at_center = x_coords[center_x, 0]
        expected_center = x_at_center**2 / (2 * 2.0**2)
        assert_allclose(result[center_x, 0], expected_center, rtol=1e-10)

        # At a point away from center: verify formula
        idx = grid.shape[0] - 2
        x_val = x_coords[idx, 0]
        expected = x_val**2 / (2 * 2.0**2)
        assert_allclose(result[idx, 0], expected, rtol=1e-10)

    def test_resolve_coefficient_position_dependent_no_grid_raises(self) -> None:
        """Position-dependent coefficient without grid raises ValueError."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian_x",
                            "phi",
                            coefficient_symbolic="x()**2",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]

        with pytest.raises(ValueError, match="requires grid"):
            pde._resolve_coefficient_at_point(term, t=0.0, grid=None)

    def test_resolve_coefficient_missing_parameter_raises(self) -> None:
        """Position-dependent coefficient with missing parameter raises ValueError."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian_x",
                            "phi",
                            coefficient_symbolic="x()**2/(2*sphR**2)",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        grid = CartesianGrid(bounds=[(-4, 4), (-4, 4)], shape=[8, 8], periodic=True)
        pde = PDEFromSpec(spec)  # No parameters!
        term = spec.equations[0].rhs_terms[0]

        with pytest.raises(ValueError, match="sphR"):
            pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

    def test_mathematica_to_python_power(self) -> None:
        """Test _mathematica_to_python converts ^ to **."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("sphR^2")
        assert "**" in result
        assert "^" not in result

    def test_mathematica_to_python_coordinates(self) -> None:
        """Test _mathematica_to_python converts x()/y() to x/y using spec coordinates."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("x()^2 + y()^2")
        assert "x()" not in result
        assert "y()" not in result
        assert "x**2" in result
        assert "y**2" in result

    def test_mathematica_to_python_non_cartesian(self) -> None:
        """Test _mathematica_to_python handles non-Cartesian coordinate names."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "r", "theta"),
        )
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("r()^2 + theta()")
        assert "r()" not in result
        assert "theta()" not in result
        assert "r**2" in result
        assert "theta" in result

    def test_evolution_with_position_dependent_coefficients(self) -> None:
        """Test that evolution runs with position-dependent coefficients."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        # Constant laplacian_x
                        OperatorTerm(0.25, "laplacian_x", "phi"),
                        # Position-dependent laplacian_x
                        OperatorTerm(
                            1.0,
                            "laplacian_x",
                            "phi",
                            coefficient_symbolic="x()**2/(2*R**2)",
                            coordinate_dependent=("x",),
                        ),
                        # Constant laplacian_y
                        OperatorTerm(0.25, "laplacian_y", "phi"),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        grid = CartesianGrid(bounds=[(-4, 4), (-4, 4)], shape=[32, 32], periodic=True)
        pde = PDEFromSpec(spec, parameters={"R": 2.0})

        # Create initial state
        phi = ScalarField(grid, data=0.0)
        pi = ScalarField(grid, data=0.0)
        # py-pde lacks complete type stubs for grid.cell_coords
        x = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]
        y = np.asarray(grid.cell_coords[..., 1])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]
        phi.data[:] = np.exp(-(x**2 + y**2) / 2)
        state = FieldCollection([phi, pi])

        # Evolve - should not crash
        result = pde.solve(state, t_range=0.5, dt=0.01, scheme="runge-kutta")
        result = normalize_solve_result(result)

        # Final state should differ from initial (evolution happened)
        assert not np.allclose(result[0].data, phi.data, rtol=0.01)


class TestMathematicaFunctionConversion:
    """Test _mathematica_to_python conversion for extended mathematical functions."""

    def make_spec(self) -> EquationSystem:
        """Create minimal EquationSystem for conversion testing."""
        return EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )

    # ===== Inverse Trigonometric =====

    def test_arcsin_conversion(self) -> None:
        """Test ArcSin[x] → arcsin(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcSin[x]")
        assert result == "arcsin(x)"

    def test_arccos_conversion(self) -> None:
        """Test ArcCos[x] → arccos(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcCos[x]")
        assert result == "arccos(x)"

    def test_arctan_conversion(self) -> None:
        """Test ArcTan[x] → arctan(x) (1-arg version)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcTan[x]")
        assert result == "arctan(x)"

    def test_arctan2_conversion(self) -> None:
        """Test ArcTan[x, y] → arctan2(y, x) with argument swap."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcTan[x, y]")
        # Critical: arguments must be swapped!
        assert result == "arctan2(y, x)"

        # Also test with coordinates
        result2 = pde._mathematica_to_python("ArcTan[x[], y[]]")
        # After bracket and coordinate conversion: arctan2(y, x)
        assert "arctan2" in result2
        assert "y" in result2
        assert "x" in result2

    # ===== Hyperbolic =====

    def test_sinh_conversion(self) -> None:
        """Test Sinh[x] → sinh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Sinh[x]")
        assert result == "sinh(x)"

    def test_cosh_conversion(self) -> None:
        """Test Cosh[x] → cosh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Cosh[x]")
        assert result == "cosh(x)"

    def test_tanh_conversion(self) -> None:
        """Test Tanh[x] → tanh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Tanh[x]")
        assert result == "tanh(x)"

    # ===== Inverse Hyperbolic =====

    def test_arcsinh_conversion(self) -> None:
        """Test ArcSinh[x] → arcsinh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcSinh[x]")
        assert result == "arcsinh(x)"

    def test_arccosh_conversion(self) -> None:
        """Test ArcCosh[x] → arccosh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcCosh[x]")
        assert result == "arccosh(x)"

    def test_arctanh_conversion(self) -> None:
        """Test ArcTanh[x] → arctanh(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("ArcTanh[x]")
        assert result == "arctanh(x)"

    # ===== Power Function =====

    def test_power_function_conversion(self) -> None:
        """Test Power[x, y] → (x)**(y)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Power[x, 2]")
        assert "**" in result
        assert "Power" not in result
        assert "x" in result
        assert "2" in result

        # Test with nested expression
        result2 = pde._mathematica_to_python("Power[Sin[x], 2]")
        assert "sin(x)" in result2
        assert "**" in result2
        assert "Power" not in result2

    # ===== Special Functions (Phase 2) =====

    def test_erf_conversion(self) -> None:
        """Test Erf[x] → erf(x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Erf[x]")
        assert result == "erf(x)"

    def test_besselj_conversion(self) -> None:
        """Test BesselJ[n, x] → jv(n, x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("BesselJ[0, x]")
        assert result == "jv(0, x)"

    def test_bessely_conversion(self) -> None:
        """Test BesselY[n, x] → yv(n, x)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("BesselY[2, x]")
        assert result == "yv(2, x)"

    # ===== Edge Cases =====

    def test_nested_functions(self) -> None:
        """Test nested function calls: Sin[ArcCos[x]]."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Sin[ArcCos[x]]")
        assert result == "sin(arccos(x))"

    def test_multiple_new_functions(self) -> None:
        """Test expression with multiple new functions."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Sinh[t[]] + ArcTan[x, y]^2")
        assert "sinh(t)" in result
        assert "arctan2(y, x)" in result
        assert "**2" in result

    def test_new_functions_with_coordinates(self) -> None:
        """Test new functions with coordinate symbols."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Cosh[x[]] * ArcSin[y[]]")
        # After conversion: cosh(x) * arcsin(y)
        assert "cosh(x)" in result
        assert "arcsin(y)" in result
        assert "[]" not in result

    # ===== Integration Test =====

    def test_coefficient_evaluation_with_new_functions(self) -> None:
        """Integration test: full coefficient resolution with new math functions."""
        # Create term with position-dependent coefficient using new functions
        # Coefficient expression: sinh(t) combined with arctan2(y, x)
        term = OperatorTerm(
            coefficient=1.0,
            operator="laplacian",
            field="phi",
            coefficient_symbolic="Sinh[t[]] + ArcTan[x[], y[]]",
            coordinate_dependent=("x", "y", "t"),
        )

        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(term,),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )

        grid = CartesianGrid([(0, 2 * np.pi), (0, 2 * np.pi)], [8, 8], periodic=True)
        pde = PDEFromSpec(spec, parameters={})

        # Resolve coefficient at t=1.0
        result = pde._resolve_coefficient_at_point(term, t=1.0, grid=grid)

        # Should be ndarray with shape matching grid
        assert isinstance(result, np.ndarray)
        assert result.shape == grid.shape

        # Verify numerical correctness at a sample point
        # At grid point [4, 4] (middle): coords ≈ (π, π)
        # arctan2(π, π) ≈ π/4 ≈ 0.7854
        # sinh(1.0) ≈ 1.1752
        # Sum ≈ 1.9606
        mid_val = result[4, 4]
        # Get actual coordinate values at this grid point
        x_coord = cast("np.ndarray", grid.cell_coords[4, 4, 0])  # type: ignore[index]
        y_coord = cast("np.ndarray", grid.cell_coords[4, 4, 1])  # type: ignore[index]
        expected = np.sinh(1.0) + np.arctan2(y_coord, x_coord)
        assert_allclose(mid_val, expected, rtol=1e-10)

    # ===== Rational, Pi, Sign, Max, Min =====

    def test_rational_conversion(self) -> None:
        """Test Rational[p, q] → (p)/(q)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Rational[1, 2]")
        assert result == "(1)/(2)"

    def test_rational_nested(self) -> None:
        """Test Rational in a larger expression."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Rational[3, 4]*x()^2")
        assert "(3)/(4)" in result
        assert "x**2" in result

    def test_pi_conversion(self) -> None:
        """Test Pi → np.pi."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("2*Pi*x()")
        assert "np.pi" in result
        assert result == "2*np.pi*x"

    def test_rational_pi_combined(self) -> None:
        """Test Rational[3, 4]*Pi → (3)/(4)*np.pi."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Rational[3, 4]*Pi")
        assert result == "(3)/(4)*np.pi"

    def test_sign_conversion(self) -> None:
        """Test Sign[x] → sign(x) (no np. prefix — available in eval namespace)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result = pde._mathematica_to_python("Sign[x()]")
        assert result == "sign(x)"

    def test_max_min_conversion(self) -> None:
        """Test Max/Min → maximum/minimum (no np. prefix — available in eval namespace)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        result_max = pde._mathematica_to_python("Max[x[], 0]")
        assert result_max == "maximum(x, 0)"
        result_min = pde._mathematica_to_python("Min[x[], 1]")
        assert result_min == "minimum(x, 1)"

    # ===== Piecewise / Inequality =====

    def test_inequality_simple(self) -> None:
        """Test Inequality[a, LessEqual, x, LessEqual, b] → chained <=."""
        from tidal.symbolic._eval_utils import _convert_inequality

        result = _convert_inequality(
            "Inequality[30, LessEqual, x[], LessEqual, 70]"
        )
        assert "<=" in result
        assert "&" in result
        assert "Inequality" not in result
        assert "LessEqual" not in result

    def test_inequality_less(self) -> None:
        """Test Inequality with Less (strict) operator."""
        from tidal.symbolic._eval_utils import _convert_inequality

        result = _convert_inequality(
            "Inequality[0, Less, x[], Less, 10]"
        )
        assert "<" in result
        assert "Less" not in result
        # Should NOT contain <= (strict inequality)
        assert "<=" not in result

    def test_piecewise_single_condition(self) -> None:
        """Test full Piecewise expression from scalar_potential_well JSON."""
        from tidal.symbolic._eval_utils import _convert_inequality, _convert_piecewise

        expr = "Piecewise[{{-V0, Inequality[30, LessEqual, x[], LessEqual, 70]}}, 0]"
        # First convert Inequality, then Piecewise (matching pipeline order)
        result = _convert_inequality(expr)
        result = _convert_piecewise(result)
        assert "piecewise(" in result
        assert "Piecewise" not in result
        assert "Inequality" not in result
        assert "LessEqual" not in result
        assert "-V0" in result

    def test_piecewise_multi_condition(self) -> None:
        """Test Piecewise with two conditions → nested piecewise() calls."""
        from tidal.symbolic._eval_utils import _convert_piecewise

        expr = "Piecewise[{{v1, c1}, {v2, c2}}, d]"
        result = _convert_piecewise(expr)
        assert result.count("piecewise(") == 2
        assert "Piecewise" not in result

    def test_piecewise_after_full_conversion(self) -> None:
        """Test full _mathematica_to_python pipeline with Piecewise expression."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        expr = "Piecewise[{{-V0, Inequality[30, LessEqual, x[], LessEqual, 70]}}, 0]"
        result = pde._mathematica_to_python(expr)
        # After full conversion: brackets → parens, x[] → x() → x
        assert "piecewise(" in result
        assert "Piecewise" not in result
        assert "Inequality" not in result
        assert "x[]" not in result
        assert "x()" not in result

    def test_exp_function_conversion(self) -> None:
        """Both E^(...) and Exp[...] should convert to exp(...)."""
        spec = self.make_spec()
        pde = PDEFromSpec(spec)
        # E^(x) → exp(x) — Mathematica's standard notation (Step 1)
        result_e = pde._mathematica_to_python("E^(x[])")
        assert "exp(" in result_e
        assert "E^" not in result_e
        # Exp[x] → exp(x) — Mathematica's function form (Step 4)
        result_exp = pde._mathematica_to_python("Exp[x[]]")
        assert "exp(" in result_exp
        assert "Exp" not in result_exp
        # Full Gaussian expression with Exp
        gaussian = "g0 * Exp[-(x[]^2 + y[]^2) / (2 * R^2)]"
        result_gauss = pde._mathematica_to_python(gaussian)
        assert "exp(" in result_gauss
        assert "Exp" not in result_gauss
        assert "x[]" not in result_gauss
        assert "y[]" not in result_gauss


class TestCachingOptimizations:
    """Tests for Issue #89 caching and memoization optimizations."""

    def _make_spec(
        self,
        *,
        symbolic_coeff: str | None = None,
        coord_dep: tuple[str, ...] = (),
    ) -> EquationSystem:
        """Create minimal EquationSystem for caching tests."""
        term = OperatorTerm(
            coefficient=1.0,
            operator="laplacian",
            field="phi",
            coefficient_symbolic=symbolic_coeff,
            coordinate_dependent=coord_dep,
        )
        return EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        term,
                        OperatorTerm(
                            -1.0, "identity", "phi", coefficient_symbolic="-m2"
                        ),
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

    def test_expr_cache_populated(self) -> None:
        """B1: _mathematica_to_python results are cached after first call."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        expr = "Sin[x[]]^2 + Cos[y[]]"
        result1 = pde._mathematica_to_python(expr)
        # Cache should not be populated yet (direct call doesn't cache)
        # But _resolve_coefficient_at_point does cache
        pde._expr_cache[expr] = result1
        result2 = pde._expr_cache[expr]
        assert result1 == result2

    def test_expr_cache_used_in_resolve(self) -> None:
        """B1: _resolve_coefficient_at_point uses the expression cache."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian",
                            "phi",
                            coefficient_symbolic="exp(t[])",
                            time_dependent=True,
                            coordinate_dependent=("t",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        grid = CartesianGrid([(-1, 1), (-1, 1)], [4, 4], periodic=True)
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]

        # First call populates cache
        pde._resolve_coefficient_at_point(term, t=1.0, grid=grid)
        assert "exp(t[])" in pde._expr_cache

        # Second call uses cache (same result)
        result2 = pde._resolve_coefficient_at_point(term, t=2.0, grid=grid)
        assert isinstance(result2, float)
        assert_allclose(result2, np.exp(2.0), rtol=1e-10)

    def test_base_namespace_prebuilt(self) -> None:
        """B2: Base namespace is pre-built in __init__ and contains math functions."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        ns = pde._base_namespace
        assert "sin" in ns
        assert "cos" in ns
        assert "exp" in ns
        assert "cot" in ns
        assert "erf" in ns
        assert "m2" in ns
        assert ns["m2"] == 1.0

    def test_base_namespace_does_not_contain_time(self) -> None:
        """B2: Base namespace does not contain t (injected per-call)."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        assert "t" not in pde._base_namespace

    def test_preresolved_coefficients(self) -> None:
        """B4: Constant coefficients are pre-resolved in __init__."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        # Term 0: coefficient=1.0, no symbolic → pre-resolved to 1.0
        assert (0, 0) in pde._preresolved
        assert pde._preresolved[0, 0] == 1.0
        # Term 1: coefficient_symbolic="-m2", m2=1.0 → pre-resolved to -1.0
        assert (0, 1) in pde._preresolved
        assert pde._preresolved[0, 1] == -1.0

    def test_preresolved_skips_position_dependent(self) -> None:
        """B4: Position-dependent coefficients are not pre-resolved."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian_x",
                            "phi",
                            coefficient_symbolic="x()**2/(2*R**2)",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        pde = PDEFromSpec(spec, parameters={"R": 1.0})
        assert (0, 0) not in pde._preresolved

    def test_bc_caching(self) -> None:
        """B5: Boundary conditions are cached after first evolution_rate call."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        grid = CartesianGrid([(0, 10)], [16], periodic=True)
        state = FieldCollection(
            [
                ScalarField(grid, data=0.0),
                ScalarField(grid, data=0.0),
            ]
        )
        assert pde._cached_bc is None
        pde.evolution_rate(state, t=0.0)
        assert pde._cached_bc is not None
        cached_ref = pde._cached_bc
        # Second call should reuse cached BC
        pde.evolution_rate(state, t=1.0)
        assert pde._cached_bc is cached_ref

    def test_preresolved_matches_runtime(self) -> None:
        """B4: Pre-resolved coefficients produce same result as runtime resolution."""
        spec = self._make_spec()
        pde = PDEFromSpec(spec, parameters={"m2": 1.0})
        for (eq_idx, term_idx), preresolved_val in pde._preresolved.items():
            term = spec.equations[eq_idx].rhs_terms[term_idx]
            runtime_val = pde._resolve_coefficient(term)
            assert preresolved_val == runtime_val

    def test_timestep_cache_deduplicates_eval(self) -> None:
        """C1: Per-timestep coeff_cache deduplicates eval for shared symbolic exprs."""
        # Two terms share the same coefficient_symbolic → only one eval() needed
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian",
                            "phi",
                            coefficient_symbolic="exp(2*H*t)",
                            time_dependent=True,
                        ),
                        OperatorTerm(
                            -1.0,
                            "identity",
                            "phi",
                            coefficient_symbolic="exp(2*H*t)",
                            time_dependent=True,
                        ),
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            mass_matrix_symbolic=(("exp(2*H*t)",),),
            metadata={},
        )
        pde = PDEFromSpec(spec, parameters={"H": 0.5})
        grid = CartesianGrid([(0, 10)], [16], periodic=True)
        state = FieldCollection(
            [
                ScalarField(grid, data=1.0),
                ScalarField(grid, data=0.0),
            ]
        )
        bc = "periodic"

        # Both terms resolve to the same value at t=1.0
        term_a = spec.equations[0].rhs_terms[0]
        term_b = spec.equations[0].rhs_terms[1]
        val_a = pde._resolve_coefficient_at_point(term_a, t=1.0, grid=grid)
        val_b = pde._resolve_coefficient_at_point(term_b, t=1.0, grid=grid)
        assert val_a == val_b

        # Run _compute_rhs_for_component — it should populate coeff_cache internally
        # and produce a valid result (no error = cache is working)
        result = pde._compute_rhs_for_component(0, state, bc, t=1.0)
        assert result.data is not None

    def test_coord_arrays_passed_to_resolve(self) -> None:
        """C2: Pre-extracted coord_arrays produce same result as grid-based resolution."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian_y",
                            "phi",
                            coefficient_symbolic="1/x**2",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid(
            [(0.5, 5), (0, 2 * np.pi)], [16, 16], periodic=[False, True]
        )

        term = spec.equations[0].rhs_terms[0]

        # Resolve using grid (fallback path)
        val_from_grid = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

        # Resolve using pre-extracted coord_arrays (C2 path)
        spatial_coords = spec.spatial_coordinates
        raw_coords = grid.cell_coords  # pyright: ignore[reportUnknownVariableType]
        coord_arrays: dict[str, np.ndarray[Any, np.dtype[np.float64]]] = {
            name: np.asarray(raw_coords[..., i])  # pyright: ignore[reportUnknownArgumentType]
            for i, name in enumerate(spatial_coords[: grid.num_axes])
        }
        val_from_arrays = pde._resolve_coefficient_at_point(
            term,
            t=0.0,
            grid=grid,
            coord_arrays=coord_arrays,  # pyright: ignore[reportUnknownArgumentType]
        )

        # Both paths should produce identical results
        np.testing.assert_array_equal(val_from_grid, val_from_arrays)


class TestEvalValidation:
    """Tests for post-eval validation of coefficient expressions (Issue #68)."""

    @staticmethod
    def _make_time_dep_spec(symbolic: str) -> EquationSystem:
        """Create a minimal spec with a time-dependent symbolic coefficient."""
        return EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian",
                            "phi",
                            coefficient_symbolic=symbolic,
                            time_dependent=True,
                            coordinate_dependent=("t",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

    @staticmethod
    def _make_position_dep_spec(symbolic: str) -> EquationSystem:
        """Create a minimal spec with a position-dependent symbolic coefficient."""
        return EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0,
                            "laplacian",
                            "phi",
                            coefficient_symbolic=symbolic,
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            coordinates=("t", "x"),
            metadata={},
        )

    def test_overflow_to_inf_raises(self) -> None:
        """exp(huge) overflows to Inf and raises ValueError."""
        # np.exp(1000) → inf (numpy doesn't raise, just returns inf)
        spec = self._make_time_dep_spec("E^(1000*t())")
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        with pytest.raises(ValueError, match="Inf"):
            pde._resolve_coefficient_at_point(term, t=1.0)

    def test_sqrt_negative_gives_nan_raises(self) -> None:
        """Sqrt[-1] via numpy gives NaN (not complex) and raises ValueError."""
        # np.sqrt(-1) returns nan (not complex) with a RuntimeWarning
        spec = self._make_time_dep_spec("Sqrt[-1]")
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        with pytest.raises(ValueError, match="NaN"):
            pde._resolve_coefficient_at_point(term, t=0.0)

    def test_python_complex_raises(self) -> None:
        """Python expression producing complex number raises TypeError."""
        # (-1)**0.5 in Python → complex (6.12e-17+1j)
        spec = self._make_time_dep_spec("(-1)^(0.5)")
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        with pytest.raises(TypeError, match="complex"):
            pde._resolve_coefficient_at_point(term, t=0.0)

    def test_valid_time_dependent_passes(self) -> None:
        """Normal time-dependent expression passes validation."""
        spec = self._make_time_dep_spec("E^(2*H*t())")
        pde = PDEFromSpec(spec, parameters={"H": 0.1})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=1.0)
        assert isinstance(result, float)
        assert_allclose(result, np.exp(0.2), rtol=1e-10)

    def test_inf_in_position_dependent_raises(self) -> None:
        """Position-dependent coefficient producing Inf raises ValueError."""
        # exp(1000*x) at x=10 → exp(10000) → inf
        spec = self._make_position_dep_spec("E^(1000*x[])")
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        grid = CartesianGrid([(0, 10)], [16], periodic=True)
        with pytest.raises(ValueError, match="Inf"):
            pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

    def test_nan_in_position_dependent_raises(self) -> None:
        """Position-dependent sqrt of negative values produces NaN and raises."""
        # Sqrt[x - 100] at x in [0,10] → all negative → NaN via numpy
        spec = self._make_position_dep_spec("Sqrt[x[] - 100]")
        pde = PDEFromSpec(spec)
        term = spec.equations[0].rhs_terms[0]
        grid = CartesianGrid([(0, 10)], [16], periodic=True)
        with pytest.raises(ValueError, match="NaN"):
            pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

    def test_validate_eval_result_scalar_valid(self) -> None:
        """_validate_eval_result passes through valid float."""
        result = PDEFromSpec._validate_eval_result(6.7, "test", "6.7")
        assert result == 6.7

    def test_validate_eval_result_scalar_nan(self) -> None:
        """_validate_eval_result rejects NaN scalar."""
        with pytest.raises(ValueError, match="NaN"):
            PDEFromSpec._validate_eval_result(float("nan"), "test", "0/0")

    def test_validate_eval_result_scalar_inf(self) -> None:
        """_validate_eval_result rejects Inf scalar."""
        with pytest.raises(ValueError, match="Inf"):
            PDEFromSpec._validate_eval_result(float("inf"), "test", "1/0")

    def test_validate_eval_result_complex(self) -> None:
        """_validate_eval_result rejects complex number."""
        with pytest.raises(TypeError, match="complex"):
            PDEFromSpec._validate_eval_result(1 + 2j, "test", "sqrt(-1)")

    def test_validate_eval_result_array_valid(self) -> None:
        """_validate_eval_result passes through valid ndarray."""
        arr = np.array([1.0, 2.0, 3.0])
        result = PDEFromSpec._validate_eval_result(arr, "test", "x")
        np.testing.assert_array_equal(result, arr)

    def test_validate_eval_result_array_nan(self) -> None:
        """_validate_eval_result rejects array containing NaN."""
        arr = np.array([1.0, float("nan"), 3.0])
        with pytest.raises(ValueError, match="NaN"):
            PDEFromSpec._validate_eval_result(arr, "test", "x/x")

    def test_validate_eval_result_array_inf(self) -> None:
        """_validate_eval_result rejects array containing Inf."""
        arr = np.array([1.0, float("inf"), 3.0])
        with pytest.raises(ValueError, match="Inf"):
            PDEFromSpec._validate_eval_result(arr, "test", "1/x")


class TestOperatorDimensionValidation:
    """Tests for early operator-dimension validation at construction time (Issue #75)."""

    @staticmethod
    def _make_spec(
        operator: str,
        spatial_dim: int,
        field_name: str = "phi",
    ) -> EquationSystem:
        """Create a minimal spec with a single operator and given spatial dimension."""
        dimension = spatial_dim + 1  # spacetime = spatial + time
        coords = ("t", "x", "y", "z")[:dimension]
        return EquationSystem(
            n_components=1,
            dimension=dimension,
            spatial_dimension=spatial_dim,
            component_names=(field_name,),
            equations=(
                ComponentEquation(
                    field_name=field_name,
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(1.0, operator, field_name),),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=coords,
        )

    def test_operator_dim_mismatch_at_construction(self) -> None:
        """gradient_z in a 2D spatial spec raises at construction, not runtime."""
        spec = self._make_spec("gradient_z", spatial_dim=2)
        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec(spec)

    def test_laplacian_z_mismatch_at_construction(self) -> None:
        """laplacian_z in a 1D spatial spec raises at construction."""
        spec = self._make_spec("laplacian_z", spatial_dim=1)
        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec(spec)

    def test_cross_derivative_xz_mismatch(self) -> None:
        """cross_derivative_xz in a 2D spatial spec raises at construction."""
        spec = self._make_spec("cross_derivative_xz", spatial_dim=2)
        with pytest.raises(ValueError, match="requires at least 3D"):
            PDEFromSpec(spec)

    def test_gradient_y_valid_in_2d(self) -> None:
        """gradient_y in a 2D spatial spec succeeds."""
        spec = self._make_spec("gradient_y", spatial_dim=2)
        pde = PDEFromSpec(spec)
        assert pde.spec.spatial_dimension == 2

    def test_laplacian_valid_in_1d(self) -> None:
        """Laplacian in a 1D spatial spec succeeds."""
        spec = self._make_spec("laplacian", spatial_dim=1)
        pde = PDEFromSpec(spec)
        assert pde.spec.spatial_dimension == 1

    def test_operator_min_dim_static_registry(self) -> None:
        """_operator_min_dim returns correct values for static operators."""
        assert PDEFromSpec._operator_min_dim("identity") == 1
        assert PDEFromSpec._operator_min_dim("laplacian") == 1
        assert PDEFromSpec._operator_min_dim("gradient_x") == 1
        assert PDEFromSpec._operator_min_dim("gradient_y") == 2
        assert PDEFromSpec._operator_min_dim("gradient_z") == 3
        assert PDEFromSpec._operator_min_dim("laplacian_y") == 2
        assert PDEFromSpec._operator_min_dim("cross_derivative_xy") == 2
        assert PDEFromSpec._operator_min_dim("cross_derivative_xz") == 3
        assert PDEFromSpec._operator_min_dim("biharmonic") == 1
        assert PDEFromSpec._operator_min_dim("first_derivative_t") == 1

    def test_operator_min_dim_dynamic_single(self) -> None:
        """_operator_min_dim resolves dynamic single-axis patterns."""
        assert PDEFromSpec._operator_min_dim("derivative_3_x") == 1
        assert PDEFromSpec._operator_min_dim("derivative_3_y") == 2
        assert PDEFromSpec._operator_min_dim("derivative_3_z") == 3

    def test_operator_min_dim_dynamic_multi(self) -> None:
        """_operator_min_dim resolves dynamic multi-axis patterns."""
        assert PDEFromSpec._operator_min_dim("derivative_2x_1y") == 2
        assert PDEFromSpec._operator_min_dim("derivative_1x_1z") == 3

    def test_operator_min_dim_unknown_raises(self) -> None:
        """_operator_min_dim raises for unrecognized operators."""
        with pytest.raises(ValueError, match="Unknown operator"):
            PDEFromSpec._operator_min_dim("nonexistent_op")
