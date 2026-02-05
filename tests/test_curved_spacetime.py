"""Tests for curved spacetime functionality.

These tests verify the Phase 2 implementation of curved spacetime support,
including:
- first_derivative_t operator
- time_dependent coefficient flag
- Hubble friction in de Sitter expansion
"""

from pathlib import Path

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)
from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec, build_pde_from_json


# Path to test data
DATA_DIR = Path(__file__).parent.parent / "examples" / "data"


class TestFirstDerivativeTimeOperator:
    """Tests for the first_derivative_t operator."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run de_sitter_kg.wls first")
        return load_equation_system(json_path)

    @pytest.fixture
    def grid_1d(self) -> CartesianGrid:
        """Create a simple 1D grid."""
        return CartesianGrid(bounds=[(0, 10)], shape=[32], periodic=True)

    def test_de_sitter_has_first_derivative_t(self, de_sitter_spec: EquationSystem) -> None:
        """Verify de Sitter spec contains first_derivative_t operator."""
        eq = de_sitter_spec.equations[0]
        operators = [term.operator for term in eq.rhs_terms]
        assert "first_derivative_t" in operators

    def test_first_derivative_t_term_structure(self, de_sitter_spec: EquationSystem) -> None:
        """Verify first_derivative_t term has correct structure for 1+1D."""
        eq = de_sitter_spec.equations[0]
        first_deriv_term = next(
            term for term in eq.rhs_terms if term.operator == "first_derivative_t"
        )

        # Should have symbolic coefficient for Hubble friction
        assert first_deriv_term.coefficient_symbolic is not None
        # For 1+1D (n=1 spatial dim), coefficient should be -H, not -2H
        assert first_deriv_term.coefficient_symbolic == "-H"

    def test_first_derivative_t_evolution(
        self, de_sitter_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that first_derivative_t correctly uses momentum field."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"H": 0.1, "m2": 1.0})

        # Create state with field=0, momentum=1
        field = ScalarField(grid_1d, data=0.0, label="phi")
        momentum = ScalarField(grid_1d, data=1.0, label="pi")
        state = FieldCollection([field, momentum])

        # Compute evolution rate
        rates = pde.evolution_rate(state, t=0.0)

        # d/dt momentum should include contribution from first_derivative_t
        # first_derivative_t(phi) returns pi, coefficient is -2*H = -0.2
        # So this term contributes -0.2 * 1.0 = -0.2 to d(momentum)/dt
        d_momentum_dt = rates[1]

        # The momentum rate should be negative (Hubble friction)
        # Note: there are other terms too (laplacian, mass), but with phi=0
        # only the first_derivative_t term contributes
        assert np.all(d_momentum_dt.data < 0), "Hubble friction should give negative rate"


class TestTimeDependentCoefficients:
    """Tests for time-dependent coefficient handling."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run de_sitter_kg.wls first")
        return load_equation_system(json_path)

    def test_time_dependent_flag_parsing(self) -> None:
        """Test that time_dependent flag is correctly parsed from JSON."""
        # Create a term with time_dependent=True
        term_data = {
            "coefficient": 1.0,
            "operator": "identity",
            "field": "phi_0",
            "coefficient_symbolic": "-m2*exp(2*H*t)",
            "time_dependent": True,
        }
        term = OperatorTerm.from_dict(term_data)
        assert term.time_dependent is True

    def test_time_dependent_flag_default(self) -> None:
        """Test that time_dependent defaults to False."""
        term_data = {
            "coefficient": 1.0,
            "operator": "laplacian",
            "field": "phi_0",
        }
        term = OperatorTerm.from_dict(term_data)
        assert term.time_dependent is False

    def test_resolve_coefficient_hubble(self, de_sitter_spec: EquationSystem) -> None:
        """Test coefficient resolution for Hubble friction term in 1+1D."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"H": 0.5, "m2": 1.0})

        # Find the first_derivative_t term
        eq = de_sitter_spec.equations[0]
        hubble_term = next(
            term for term in eq.rhs_terms if term.operator == "first_derivative_t"
        )

        # Resolve at t=0 (time doesn't matter for -H which is constant)
        # For 1+1D (n=1), the symbolic coefficient is "-H", should resolve to -0.5
        coeff = pde._resolve_coefficient_at_time(hubble_term, t=0.0)

        # Should be -H = -0.5 (not -2H which is for 2+1D)
        assert coeff == pytest.approx(-0.5, rel=0.01)  # Tightened tolerance from 10% to 1%


class TestDeSitterSimulation:
    """Integration tests for de Sitter simulation."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run de_sitter_kg.wls first")
        return load_equation_system(json_path)

    @pytest.fixture
    def grid_1d(self) -> CartesianGrid:
        """Create a 1D grid."""
        return CartesianGrid(bounds=[(0, 50)], shape=[128], periodic=True)

    def test_de_sitter_energy_decreases(
        self, de_sitter_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that energy decreases due to Hubble friction."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"H": 0.2, "m2": 1.0})

        # Create Gaussian initial condition
        x = grid_1d.cell_coords[..., 0]
        gaussian = np.exp(-((x - 25) ** 2) / (2 * 3**2))

        field = ScalarField(grid_1d, data=gaussian, label="phi")
        momentum = ScalarField(grid_1d, data=0.0, label="pi")
        state = FieldCollection([field, momentum])

        # Compute energy proxy
        def energy(s: FieldCollection) -> float:
            phi = s[0].data
            pi = s[1].data
            return float(np.sum(phi**2 + pi**2))

        initial_energy = energy(state)

        # Run for a short time
        result = pde.solve(state, t_range=5.0, dt=0.01)

        final_energy = energy(result)

        # Energy should decrease due to Hubble friction
        assert final_energy < initial_energy, "Hubble friction should decrease energy"

    def test_de_sitter_amplitude_decays(
        self, de_sitter_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Test that amplitude decays due to Hubble friction."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"H": 0.3, "m2": 0.5})

        # Create Gaussian initial condition
        x = grid_1d.cell_coords[..., 0]
        gaussian = np.exp(-((x - 25) ** 2) / (2 * 3**2))

        field = ScalarField(grid_1d, data=gaussian, label="phi")
        momentum = ScalarField(grid_1d, data=0.0, label="pi")
        state = FieldCollection([field, momentum])

        initial_max = np.max(np.abs(state[0].data))

        # Run simulation
        result = pde.solve(state, t_range=10.0, dt=0.01)

        final_max = np.max(np.abs(result[0].data))

        # Amplitude should decay
        assert final_max < initial_max, "Amplitude should decay due to Hubble friction"


class TestStaticConformalEquivalence:
    """Tests for static conformal factor (Phase 1 verification)."""

    @pytest.fixture
    def conformal_spec(self) -> EquationSystem:
        """Load static conformal KG specification."""
        json_path = DATA_DIR / "conformal_kg_static.json"
        if not json_path.exists():
            pytest.skip(
                "conformal_kg_static.json not found - run conformal_kg_static.wls first"
            )
        return load_equation_system(json_path)

    def test_conformal_has_effective_mass(self, conformal_spec: EquationSystem) -> None:
        """Verify static conformal spec has effective mass term."""
        eq = conformal_spec.equations[0]

        # Find identity (mass) term
        identity_term = next(
            (term for term in eq.rhs_terms if term.operator == "identity"), None
        )

        assert identity_term is not None, "Should have identity (mass) term"
        # Effective mass should be m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4
        assert identity_term.coefficient == pytest.approx(-4.0, rel=0.01)

    def test_conformal_no_first_derivative_t(
        self, conformal_spec: EquationSystem
    ) -> None:
        """Static conformal should NOT have first_derivative_t (no Hubble friction)."""
        eq = conformal_spec.equations[0]
        operators = [term.operator for term in eq.rhs_terms]

        assert "first_derivative_t" not in operators
