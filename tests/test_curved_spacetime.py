"""Tests for curved spacetime functionality.

These tests verify the Phase 2 implementation of curved spacetime support,
including:
- first_derivative_t operator
- time_dependent coefficient flag
- Hubble friction in de Sitter expansion
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, ScalarField

from tidal.symbolic.json_loader import (
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)
from tidal.symbolic.pde_builder import PDEFromSpec

# Path to test data
DATA_DIR = Path(__file__).parent.parent / "examples" / "data"


class TestFirstDerivativeTimeOperator:
    """Tests for the first_derivative_t operator."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run `tidal derive de_sitter.toml` first")
        return load_equation_system(json_path)

    @pytest.fixture
    def grid_2d(self) -> CartesianGrid:
        """Create a simple 2D grid for 2+1D de Sitter."""
        return CartesianGrid(bounds=[(0, 10), (0, 10)], shape=[32, 32], periodic=True)

    def test_de_sitter_has_first_derivative_t(
        self, de_sitter_spec: EquationSystem
    ) -> None:
        """Verify de Sitter spec contains first_derivative_t operator."""
        eq = de_sitter_spec.equations[0]
        operators = [term.operator for term in eq.rhs_terms]
        assert "first_derivative_t" in operators

    def test_first_derivative_t_term_structure(
        self, de_sitter_spec: EquationSystem
    ) -> None:
        """Verify first_derivative_t term has correct structure for 2+1D."""
        eq = de_sitter_spec.equations[0]
        first_deriv_term = next(
            term for term in eq.rhs_terms if term.operator == "first_derivative_t"
        )

        # Should have symbolic coefficient for Hubble friction
        assert first_deriv_term.coefficient_symbolic is not None
        # For 2+1D (n=2 spatial dims), coefficient is -(n-1)*H = -H
        assert first_deriv_term.coefficient_symbolic == "-dSH"

    def test_first_derivative_t_evolution(
        self, de_sitter_spec: EquationSystem, grid_2d: CartesianGrid
    ) -> None:
        """Test that first_derivative_t correctly uses momentum field."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"dSH": 0.1, "dSm2": 1.0})

        # Create state with field=0, momentum=1
        field = ScalarField(grid_2d, data=0.0, label="phi")
        momentum = ScalarField(grid_2d, data=1.0, label="pi")
        state = FieldCollection([field, momentum])

        # Compute evolution rate
        rates = pde.evolution_rate(state, t=0.0)

        # d/dt momentum should include contribution from first_derivative_t
        # first_derivative_t(phi) returns pi, coefficient is -dSH = -0.1
        # So this term contributes -0.1 * 1.0 = -0.1 to d(momentum)/dt
        d_momentum_dt = rates[1]

        # The momentum rate should be negative (Hubble friction)
        # Note: there are other terms too (laplacian, mass), but with phi=0
        # only the first_derivative_t term contributes
        assert np.all(d_momentum_dt.data < 0), (
            "Hubble friction should give negative rate"
        )


class TestTimeDependentCoefficients:
    """Tests for time-dependent coefficient handling."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run `tidal derive de_sitter.toml` first")
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
        """Test coefficient resolution for Hubble friction term in 2+1D."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"dSH": 0.5, "dSm2": 1.0})

        # Find the first_derivative_t term
        eq = de_sitter_spec.equations[0]
        hubble_term = next(
            term for term in eq.rhs_terms if term.operator == "first_derivative_t"
        )

        # Resolve at t=0 (time doesn't matter for -dSH which is constant)
        # For 2+1D (n=2), the symbolic coefficient is "-dSH", should resolve to -0.5
        coeff = pde._resolve_coefficient_at_point(hubble_term, t=0.0)

        # Should be -dSH = -0.5
        assert coeff == pytest.approx(-0.5, rel=0.01)


class TestDeSitterSimulation:
    """Integration tests for de Sitter simulation (2+1D)."""

    @pytest.fixture
    def de_sitter_spec(self) -> EquationSystem:
        """Load de Sitter KG specification."""
        json_path = DATA_DIR / "de_sitter_kg.json"
        if not json_path.exists():
            pytest.skip("de_sitter_kg.json not found - run `tidal derive de_sitter.toml` first")
        return load_equation_system(json_path)

    @pytest.fixture
    def grid_2d(self) -> CartesianGrid:
        """Create a 2D grid for 2+1D de Sitter."""
        return CartesianGrid(bounds=[(0, 50), (0, 50)], shape=[64, 64], periodic=True)

    def test_de_sitter_energy_decreases(
        self, de_sitter_spec: EquationSystem, grid_2d: CartesianGrid
    ) -> None:
        """Test that energy decreases due to Hubble friction."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"dSH": 0.2, "dSm2": 1.0})

        # Create 2D Gaussian initial condition
        # py-pde cell_coords lacks type stubs
        coords = np.asarray(grid_2d.cell_coords)  # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]
        x = cast("np.ndarray", coords[..., 0])
        y = cast("np.ndarray", coords[..., 1])
        gaussian = np.exp(-((x - 25) ** 2 + (y - 25) ** 2) / (2 * 3**2))

        field = ScalarField(grid_2d, data=gaussian, label="phi")
        momentum = ScalarField(grid_2d, data=0.0, label="pi")
        state = FieldCollection([field, momentum])

        # Compute energy proxy
        def energy(s: FieldCollection) -> float:
            phi = s[0].data
            pi = s[1].data
            return float(np.sum(phi**2 + pi**2))

        initial_energy = energy(state)

        # Run for a short time
        result = pde.solve(state, t_range=5.0, dt=0.01)
        assert isinstance(result, FieldCollection)

        final_energy = energy(result)

        # Energy should decrease due to Hubble friction
        assert final_energy < initial_energy, "Hubble friction should decrease energy"

    def test_de_sitter_amplitude_decays(
        self, de_sitter_spec: EquationSystem, grid_2d: CartesianGrid
    ) -> None:
        """Test that amplitude decays due to Hubble friction."""
        pde = PDEFromSpec(de_sitter_spec, parameters={"dSH": 0.3, "dSm2": 0.5})

        # Create 2D Gaussian initial condition
        # py-pde cell_coords lacks type stubs
        coords = np.asarray(grid_2d.cell_coords)  # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]
        x = cast("np.ndarray", coords[..., 0])
        y = cast("np.ndarray", coords[..., 1])
        gaussian = np.exp(-((x - 25) ** 2 + (y - 25) ** 2) / (2 * 3**2))

        field = ScalarField(grid_2d, data=gaussian, label="phi")
        momentum = ScalarField(grid_2d, data=0.0, label="pi")
        state = FieldCollection([field, momentum])

        initial_max = np.max(np.abs(state[0].data))

        # Run simulation
        result = pde.solve(state, t_range=10.0, dt=0.01)
        assert isinstance(result, FieldCollection)

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
                "conformal_kg_static.json not found - run `tidal derive conformal_static.toml` first"
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
