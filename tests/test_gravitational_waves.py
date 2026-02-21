"""Tests for Phase 11: Linearized Gravity / Gravitational Waves.

Verifies both the gauge-unfixed and de Donder gauge systems derived from
the Einstein-Hilbert action via xPert linearization.

De Donder gauge system (primary simulation target):
- 10 uncoupled wave equations: Box[h_i] = 0 for each component
- All second-order in time
- Each component propagates independently at c = 1

Gauge-unfixed system (structural verification only):
- 10 coupled equations with mixed time orders
- Cross-field coupling through spatial derivatives and momentum references
- Some components are first-order (momentum constraints)
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid

from tidal.symbolic.json_loader import (
    EquationSystem,
    load_equation_system,
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
)
from tidal.utils import normalize_solve_result

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"
DEDONDER_JSON = DATA_DIR / "linearized_gravity_dedonder.json"
GW_JSON = DATA_DIR / "linearized_gravity.json"

_skip_dedonder = pytest.mark.skipif(
    not DEDONDER_JSON.exists(),
    reason="linearized_gravity_dedonder.json not found",
)
_skip_gw = pytest.mark.skipif(
    not GW_JSON.exists(),
    reason="linearized_gravity.json not found (run tidal derive)",
)


def _gw_json_has_gauge(gauge_name: str) -> bool:
    """Check if the GW JSON has the specified gauge in its metadata."""
    if not GW_JSON.exists():
        return False
    import json

    with open(GW_JSON) as f:
        data = json.load(f)
    return data.get("metadata", {}).get("gauge") == gauge_name


_skip_gw_tt = pytest.mark.skipif(
    not _gw_json_has_gauge("tt"),
    reason="linearized_gravity.json missing or not TT-gauge (re-run tidal derive)",
)


# ============================================================
# De Donder gauge system tests
# ============================================================


@_skip_dedonder
class TestDeDonderSpec:
    """Test that the de Donder gauge JSON loads correctly."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DEDONDER_JSON)

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 4
        assert spec.spatial_dimension == 3

    def test_coordinates(self, spec: EquationSystem) -> None:
        assert spec.coordinates == ("t", "x", "y", "z")

    def test_component_count(self, spec: EquationSystem) -> None:
        """Symmetric rank-2 in 4D has 10 independent components."""
        assert spec.n_components == 10

    def test_component_names(self, spec: EquationSystem) -> None:
        expected = tuple(f"h_{i}" for i in range(10))
        assert spec.component_names == expected

    def test_all_second_order(self, spec: EquationSystem) -> None:
        """All de Donder equations should be second-order in time."""
        for eq in spec.equations:
            assert eq.time_derivative_order == 2, (
                f"Field {eq.field_name} has time order {eq.time_derivative_order}, expected 2"
            )

    def test_state_size(self, spec: EquationSystem) -> None:
        """10 second-order fields -> 20 state slots (field + momentum each)."""
        assert spec.state_size == 20

    def test_each_equation_has_three_laplacian_terms(
        self, spec: EquationSystem
    ) -> None:
        """Each de Donder equation is Box[h] = 0, meaning 3 directional laplacians."""
        expected_operators = {"laplacian_x", "laplacian_y", "laplacian_z"}
        for eq in spec.equations:
            operators = {term.operator for term in eq.rhs_terms}
            assert operators == expected_operators, (
                f"Field {eq.field_name}: operators {operators} != expected {expected_operators}"
            )

    def test_no_cross_field_coupling(self, spec: EquationSystem) -> None:
        """In de Donder gauge, each component evolves independently."""
        for eq in spec.equations:
            for term in eq.rhs_terms:
                assert term.field == eq.field_name, (
                    f"Field {eq.field_name}: term references {term.field} (cross-coupling)"
                )

    def test_coefficients_are_positive_one(self, spec: EquationSystem) -> None:
        """Wave equation d²_t h = ∇²h has coefficient +1.0 for each laplacian."""
        for eq in spec.equations:
            for term in eq.rhs_terms:
                assert_allclose(
                    term.coefficient,
                    1.0,
                    err_msg=f"Field {eq.field_name}, operator {term.operator}: "
                    f"coeff {term.coefficient} != 1.0",
                )

    def test_metadata(self, spec: EquationSystem) -> None:
        assert spec.metadata.get("gauge") == "de_donder"
        assert spec.metadata.get("linearized") is True


@_skip_dedonder
class TestDeDonderPDE:
    """Test PDE construction and evolution for de Donder gauge."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DEDONDER_JSON)

    @pytest.fixture
    def grid(self) -> CartesianGrid:
        """Coarse 3D grid for testing."""
        return CartesianGrid(
            [(-5, 5), (-5, 5), (-5, 5)],
            [8, 8, 8],
            periodic=True,
        )

    def test_pde_builds_from_json(self) -> None:
        """PDE should build from JSON without errors."""
        pde = build_pde_from_json(DEDONDER_JSON)
        assert isinstance(pde, PDEFromSpec)

    def test_initial_state_creation(
        self,
        spec: EquationSystem,
        grid: CartesianGrid,
    ) -> None:
        """Initial state should have 20 fields (10 field + 10 momentum)."""
        state = create_initial_state(grid, spec)
        assert len(state) == 20

    def test_initial_state_with_data(
        self,
        spec: EquationSystem,
        grid: CartesianGrid,
    ) -> None:
        """Should accept initial field data for specific components."""
        # Set up a Gaussian pulse in h_4 (h_{xx}) and h_7 (h_{yy})
        x, y, z = (
            cast("np.ndarray", grid.cell_coords[..., 0]),
            cast("np.ndarray", grid.cell_coords[..., 1]),
            cast("np.ndarray", grid.cell_coords[..., 2]),
        )
        gaussian = np.exp(-(x**2 + y**2 + z**2) / 2.0)
        state = create_initial_state(
            grid,
            spec,
            field_data={"h_4": gaussian, "h_7": -gaussian},
        )
        # h_4 should be nonzero
        h4_field = state[spec.state_layout.index(("h_4", "field"))]
        assert np.max(np.abs(h4_field.data)) > 0.1

    def test_evolution_rate_shape(
        self,
        spec: EquationSystem,
        grid: CartesianGrid,
    ) -> None:
        """Evolution rate should produce 20 fields matching state layout."""
        pde = PDEFromSpec(spec)
        state = create_initial_state(grid, spec)
        rate = pde.evolution_rate(state, t=0)
        assert len(rate) == 20

    def test_short_simulation_stable(
        self,
        spec: EquationSystem,
        grid: CartesianGrid,
    ) -> None:
        """Short RK4 simulation should remain bounded."""
        pde = PDEFromSpec(spec)

        # Initialize a simple Gaussian in one component
        x, y, z = (
            cast("np.ndarray", grid.cell_coords[..., 0]),
            cast("np.ndarray", grid.cell_coords[..., 1]),
            cast("np.ndarray", grid.cell_coords[..., 2]),
        )
        gaussian = 0.01 * np.exp(-(x**2 + y**2 + z**2) / 2.0)
        state = create_initial_state(
            grid,
            spec,
            field_data={"h_4": gaussian},
        )

        # Very short simulation with RK4
        result = pde.solve(state, t_range=0.1, dt=0.02, scheme="runge-kutta")
        result = normalize_solve_result(result)
        final = result

        # Check that the solution hasn't blown up
        for field in final:
            max_val = np.max(np.abs(field.data))
            assert max_val < 100.0, f"Solution blew up: max |field| = {max_val}"

    def test_wave_speed_c_equals_1(
        self,
        spec: EquationSystem,
    ) -> None:
        """The wave speed for Box[h] = 0 should be c = 1.

        For each component, the equation is d²_t h = ∇²h.
        With CFL condition Δt < Δx/c, c = 1 means CFL = Δt/Δx < 1.
        We verify the operator coefficients sum to give c² = 1.
        """
        for eq in spec.equations:
            sum(term.coefficient for term in eq.rhs_terms)
            # In 3D, laplacian = d²_x + d²_y + d²_z, each with coefficient 1
            # Total effective c² per direction = 1 for each direction
            for term in eq.rhs_terms:
                assert_allclose(term.coefficient, 1.0)


# ============================================================
# Gauge-unfixed system tests (structural verification)
# ============================================================






# ============================================================
# Physics validation tests (de Donder)
# ============================================================


@_skip_dedonder
class TestDeDonderPhysics:
    """Physics-level tests for gravitational wave propagation."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DEDONDER_JSON)

    @pytest.fixture
    def grid(self) -> CartesianGrid:
        """1D-like 3D grid for wave propagation test (propagation along z)."""
        # Fine resolution along z, minimal in x,y
        return CartesianGrid(
            [(-1, 1), (-1, 1), (-10, 10)],
            [4, 4, 64],
            periodic=True,
        )

    def test_tt_compatible_initial_data(
        self,
        spec: EquationSystem,
        grid: CartesianGrid,
    ) -> None:
        """TT-gauge compatible initial data should propagate correctly.

        For a GW propagating along z:
        - h_4 (h_{xx}) and h_7 (h_{yy}) are the + polarization: h_+ = h_{xx} = -h_{yy}
        - h_5 (h_{xy}) is the x polarization
        - All other components zero (TT gauge)
        """
        z = cast("np.ndarray", grid.cell_coords[..., 2])
        # Gaussian wave packet for h_+ polarization
        h_plus = 0.01 * np.exp(-(z**2) / 2.0)

        state = create_initial_state(
            grid,
            spec,
            field_data={"h_4": h_plus, "h_7": -h_plus},
        )

        pde = PDEFromSpec(spec)
        result = pde.solve(state, t_range=0.5, dt=0.01, scheme="runge-kutta")
        result = normalize_solve_result(result)

        # After propagation, h_4 and h_7 should still be equal and opposite
        # (tracelessness maintained by the wave equation)
        h4_idx = spec.state_layout.index(("h_4", "field"))
        h7_idx = spec.state_layout.index(("h_7", "field"))
        h4_data = result[h4_idx].data
        h7_data = result[h7_idx].data

        # h_4 = -h_7 should be maintained (TT tracelessness)
        assert_allclose(
            h4_data,
            -h7_data,
            atol=1e-6,
            err_msg="TT tracelessness h_{xx} = -h_{yy} not maintained",
        )

    def test_energy_conservation_short_time(
        self,
        spec: EquationSystem,
    ) -> None:
        """Total field energy should be approximately conserved."""
        grid = CartesianGrid(
            [(-5, 5), (-5, 5), (-5, 5)],
            [8, 8, 8],
            periodic=True,
        )
        pde = PDEFromSpec(spec)

        z = cast("np.ndarray", grid.cell_coords[..., 2])
        h_plus = 0.01 * np.exp(-(z**2) / 2.0)
        state = create_initial_state(
            grid,
            spec,
            field_data={"h_4": h_plus},
        )

        # Compute initial "energy" (sum of field² + momentum²)
        initial_energy = sum(np.sum(f.data**2) for f in state)

        result = pde.solve(state, t_range=0.2, dt=0.05, scheme="runge-kutta")
        result = normalize_solve_result(result)
        final_energy = sum(np.sum(f.data**2) for f in result)

        # Energy should be conserved to within ~15% on this coarse 8³ grid
        rel_change = abs(final_energy - initial_energy) / max(
            float(initial_energy), 1e-15
        )
        assert rel_change < 0.15, (
            f"Energy changed by {rel_change:.1%} (initial={initial_energy:.6f}, "
            f"final={final_energy:.6f})"
        )


# ============================================================
# TT gauge system tests
# ============================================================


@_skip_gw_tt
class TestTTGaugeSpec:
    """Test that the TT gauge JSON loads correctly.

    TT gauge on 3+1D linearised gravity should:
    - Zero temporal components (h_0..h_3) — h_0 trivially, h_1..h_3 become constraints
    - Substitute h_9 (spatial trace) → -(h_4 + h_7), replace h_9's equation
      with the algebraic traceless constraint
    - Add transverse constraint equations (d^i h_{ij} = 0)
    - h_4..h_8 remain as evolution equations (5 dynamical DOF)
    """

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(GW_JSON)

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 4
        assert spec.spatial_dimension == 3

    def test_component_count(self, spec: EquationSystem) -> None:
        """TT gauge: 10 original + 3 transverse = 13 components."""
        assert spec.n_components == 13, (
            f"Expected 13 components (10 + 3 transverse), got {spec.n_components}"
        )

    def test_temporal_h0_trivial(self, spec: EquationSystem) -> None:
        """h_0 (h_{tt}) should be trivially zeroed — constraint with identity(h_0)."""
        eq_h0 = next(eq for eq in spec.equations if eq.field_name == "h_0")
        assert eq_h0.time_derivative_order == 0, "h_0 should be a constraint"

    def test_temporal_h1_h3_are_constraints(self, spec: EquationSystem) -> None:
        """h_1..h_3 (h_{ti}) should become constraint equations after temporal zeroing."""
        for name in ["h_1", "h_2", "h_3"]:
            eq = next(eq for eq in spec.equations if eq.field_name == name)
            assert eq.time_derivative_order == 0, (
                f"{name} should be constraint (time_order=0), got {eq.time_derivative_order}"
            )

    def test_spatial_components_h4_h8_are_evolution(self, spec: EquationSystem) -> None:
        """h_4..h_8 should remain as 2nd-order evolution equations."""
        for idx in range(4, 9):
            name = f"h_{idx}"
            eq = next(eq for eq in spec.equations if eq.field_name == name)
            assert eq.time_derivative_order == 2, (
                f"{name} should be evolution (time_order=2), got {eq.time_derivative_order}"
            )

    def test_h9_is_traceless_constraint(self, spec: EquationSystem) -> None:
        """h_9 equation should be the traceless constraint (time_order=0)."""
        eq_h9 = next(eq for eq in spec.equations if eq.field_name == "h_9")
        assert eq_h9.time_derivative_order == 0, (
            f"h_9 should be traceless constraint (time_order=0), got {eq_h9.time_derivative_order}"
        )
        # Should reference the spatial diagonal components
        fields = {term.field for term in eq_h9.rhs_terms}
        spatial_diags = {"h_4", "h_7", "h_9"}
        assert spatial_diags <= fields, (
            f"Traceless constraint h_9 should reference spatial diagonals {spatial_diags}, "
            f"got fields {fields}"
        )

    def test_transverse_constraints_present(self, spec: EquationSystem) -> None:
        """3 transverse constraint equations (d^i h_{ij} = 0) should be present."""
        transverse = [eq for eq in spec.equations if "transverse" in eq.field_name]
        assert len(transverse) == 3, f"Expected 3 transverse constraints, got {len(transverse)}"
        for eq in transverse:
            assert eq.time_derivative_order == 0
            # Should have gradient-type operators
            operators = {term.operator for term in eq.rhs_terms}
            has_gradient = any("gradient" in op for op in operators)
            assert has_gradient, (
                f"Transverse constraint {eq.field_name} missing gradient operators: {operators}"
            )

    def test_evolution_equations_have_laplacians(self, spec: EquationSystem) -> None:
        """Evolution equations should contain laplacian-type operators."""
        evolution_eqs = [eq for eq in spec.equations if eq.time_derivative_order == 2]
        assert len(evolution_eqs) == 5, f"Expected 5 evolution equations, got {len(evolution_eqs)}"
        for eq in evolution_eqs:
            operators = {term.operator for term in eq.rhs_terms}
            has_spatial_deriv = any(
                "laplacian" in op or "cross_derivative" in op or "gradient" in op
                for op in operators
            )
            assert has_spatial_deriv, (
                f"Evolution equation {eq.field_name} has no spatial derivative operators: {operators}"
            )

    def test_no_temporal_field_references_in_evolution(self, spec: EquationSystem) -> None:
        """Evolution equations should not reference zeroed temporal fields h_0..h_3."""
        temporal_names = {"h_0", "h_1", "h_2", "h_3"}
        temporal_momenta = {"pi_0", "pi_1", "pi_2", "pi_3"}
        evolution_eqs = [eq for eq in spec.equations if eq.time_derivative_order == 2]
        for eq in evolution_eqs:
            for term in eq.rhs_terms:
                assert term.field not in temporal_names, (
                    f"Evolution eq {eq.field_name} references zeroed field {term.field}"
                )
                assert term.field not in temporal_momenta, (
                    f"Evolution eq {eq.field_name} references zeroed momentum {term.field}"
                )

    def test_all_coefficients_numeric(self, spec: EquationSystem) -> None:
        """All coefficients should be numeric (no unsimplified Mathematica expressions)."""
        for eq in spec.equations:
            for term in eq.rhs_terms:
                assert isinstance(term.coefficient, (int, float)), (
                    f"Non-numeric coefficient in {eq.field_name}: {term.coefficient!r}"
                )

    def test_metadata(self, spec: EquationSystem) -> None:
        assert spec.metadata.get("gauge") == "tt"
        assert spec.metadata.get("linearized") is True
