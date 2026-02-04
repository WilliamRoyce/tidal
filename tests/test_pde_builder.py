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
