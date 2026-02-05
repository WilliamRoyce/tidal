"""Tests for the PDE builder from equation specifications."""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
)
from torsion_gertsenshtein.symbolic.pde_builder import (
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
)
from torsion_gertsenshtein.utils import normalize_solve_result

# ruff: noqa: SLF001  # Tests need to access private methods

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
        assert not pde.explicit_time_dependence

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

    def test_unknown_operator_raises(self, grid_1d_small: CartesianGrid) -> None:
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

        pde = PDEFromSpec(spec)
        state = FieldCollection(
            [
                ScalarField(grid_1d_small, data=1.0),
                ScalarField(grid_1d_small, data=0.0),
            ]
        )

        with pytest.raises(ValueError, match="Unknown operator"):
            pde.evolution_rate(state)


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

        pde = build_pde_from_json(kg_json_path)

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

        # Create state with Gaussian pulse in A_1
        # State layout: A_0, pi_0, A_1, pi_1, A_2, pi_2
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        a1_data = np.exp(-((x - 5) ** 2 + (y - 5) ** 2) / 2)

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=0.0),  # A_0
                ScalarField(grid_2d, data=0.0),  # pi_0
                ScalarField(grid_2d, data=a1_data),  # A_1 (Gaussian)
                ScalarField(grid_2d, data=0.0),  # pi_1
                ScalarField(grid_2d, data=0.0),  # A_2
                ScalarField(grid_2d, data=0.0),  # pi_2
            ]
        )

        # Compute evolution rate - should not raise
        rates = pde.evolution_rate(state)

        num_cs_states = 6
        assert len(rates) == num_cs_states

        # d/dt A_1 = pi_1 = 0
        assert_allclose(rates[2].data, 0.0, atol=1e-10)

        # d/dt pi_1 should include laplacian(A_1) contribution (non-zero)
        assert np.max(np.abs(rates[3].data)) > 0

    def test_chern_simons_simulation_short(
        self, cs_json_path: Path, grid_2d: CartesianGrid
    ) -> None:
        """Test short Chern-Simons simulation runs without error."""
        if not cs_json_path.exists():
            pytest.skip(f"Test file not found: {cs_json_path}")

        pde = build_pde_from_json(cs_json_path)

        # Create minimal initial state
        x = cast("np.ndarray", grid_2d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_2d.cell_coords[..., 1])
        a0_data = np.exp(-((x - 5) ** 2 + (y - 5) ** 2) / 4)

        state = FieldCollection(
            [
                ScalarField(grid_2d, data=a0_data),  # A_0
                ScalarField(grid_2d, data=0.0),  # pi_0
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
        num_cs_states = 6
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
