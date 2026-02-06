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

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    EquationSystem,
    load_equation_system,
)
from torsion_gertsenshtein.symbolic.pde_builder import (
    PDEFromSpec,
    create_initial_state,
)

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"
DEDONDER_JSON = DATA_DIR / "linearized_gravity_dedonder.json"
UNFIXED_JSON = DATA_DIR / "linearized_gravity.json"


# ============================================================
# De Donder gauge system tests
# ============================================================


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

    def test_each_equation_has_three_laplacian_terms(self, spec: EquationSystem) -> None:
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
                    term.coefficient, 1.0,
                    err_msg=f"Field {eq.field_name}, operator {term.operator}: "
                            f"coeff {term.coefficient} != 1.0",
                )

    def test_metadata(self, spec: EquationSystem) -> None:
        assert spec.metadata.get("gauge") == "de_donder"
        assert spec.metadata.get("linearized") is True


class TestDeDonderPDE:
    """Test PDE construction and evolution for de Donder gauge."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DEDONDER_JSON)

    @pytest.fixture
    def grid(self) -> CartesianGrid:
        """Coarse 3D grid for testing."""
        return CartesianGrid(
            [(-5, 5), (-5, 5), (-5, 5)], [8, 8, 8], periodic=True,
        )

    def test_pde_builds_from_json(self) -> None:
        """PDE should build from JSON without errors."""
        from torsion_gertsenshtein.symbolic.pde_builder import build_pde_from_json
        pde = build_pde_from_json(DEDONDER_JSON)
        assert isinstance(pde, PDEFromSpec)

    def test_initial_state_creation(
        self, spec: EquationSystem, grid: CartesianGrid,
    ) -> None:
        """Initial state should have 20 fields (10 field + 10 momentum)."""
        state = create_initial_state(grid, spec)
        assert len(state) == 20

    def test_initial_state_with_data(
        self, spec: EquationSystem, grid: CartesianGrid,
    ) -> None:
        """Should accept initial field data for specific components."""
        # Set up a Gaussian pulse in h_4 (h_{xx}) and h_7 (h_{yy})
        x, y, z = grid.cell_coords[..., 0], grid.cell_coords[..., 1], grid.cell_coords[..., 2]
        gaussian = np.exp(-(x**2 + y**2 + z**2) / 2.0)
        state = create_initial_state(
            grid, spec,
            field_data={"h_4": gaussian, "h_7": -gaussian},
        )
        # h_4 should be nonzero
        h4_field = state[spec.state_layout.index(("h_4", "field"))]
        assert np.max(np.abs(h4_field.data)) > 0.1

    def test_evolution_rate_shape(
        self, spec: EquationSystem, grid: CartesianGrid,
    ) -> None:
        """Evolution rate should produce 20 fields matching state layout."""
        pde = PDEFromSpec(spec)
        state = create_initial_state(grid, spec)
        rate = pde.evolution_rate(state, t=0)
        assert len(rate) == 20

    def test_short_simulation_stable(
        self, spec: EquationSystem, grid: CartesianGrid,
    ) -> None:
        """Short RK4 simulation should remain bounded."""
        pde = PDEFromSpec(spec)

        # Initialize a simple Gaussian in one component
        x, y, z = grid.cell_coords[..., 0], grid.cell_coords[..., 1], grid.cell_coords[..., 2]
        gaussian = 0.01 * np.exp(-(x**2 + y**2 + z**2) / 2.0)
        state = create_initial_state(
            grid, spec,
            field_data={"h_4": gaussian},
        )

        # Very short simulation with RK4
        result = pde.solve(state, t_range=0.1, dt=0.02, scheme="runge-kutta")
        final = result

        # Check that the solution hasn't blown up
        for field in final:
            max_val = np.max(np.abs(field.data))
            assert max_val < 100.0, f"Solution blew up: max |field| = {max_val}"

    def test_wave_speed_c_equals_1(
        self, spec: EquationSystem,
    ) -> None:
        """The wave speed for Box[h] = 0 should be c = 1.

        For each component, the equation is d²_t h = ∇²h.
        With CFL condition Δt < Δx/c, c = 1 means CFL = Δt/Δx < 1.
        We verify the operator coefficients sum to give c² = 1.
        """
        for eq in spec.equations:
            total_coeff = sum(term.coefficient for term in eq.rhs_terms)
            # In 3D, laplacian = d²_x + d²_y + d²_z, each with coefficient 1
            # Total effective c² per direction = 1 for each direction
            for term in eq.rhs_terms:
                assert_allclose(term.coefficient, 1.0)


# ============================================================
# Gauge-unfixed system tests (structural verification)
# ============================================================


class TestGaugeUnfixedSpec:
    """Structural tests for the gauge-unfixed linearized Einstein equations.

    The gauge-unfixed system preserves the full constraint structure of
    linearized GR: Hamiltonian constraint (h_0, elliptic), momentum
    constraints (h_1-h_3), and evolution equations (h_4-h_9, second-order).
    Cross-field time derivatives are correctly classified as RHS terms.
    """

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(UNFIXED_JSON)

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 4
        assert spec.spatial_dimension == 3

    def test_component_count(self, spec: EquationSystem) -> None:
        assert spec.n_components == 10

    def test_component_names(self, spec: EquationSystem) -> None:
        expected = tuple(f"h_{i}" for i in range(10))
        assert spec.component_names == expected

    def test_mixed_time_orders(self, spec: EquationSystem) -> None:
        """Gauge-unfixed system should have mixed time derivative orders.

        In the non-trace-reversed formulation (as produced by xPert + ContractMetric),
        the linearized Einstein equations split into:
        - Constraint equations (time_order 0): h_0 (Hamiltonian constraint),
          h_5 (h_xy), h_6 (h_xz), h_8 (h_yz) — off-diagonal spatial components
          where ∂²_t cancels after trace subtraction
        - Evolution equations (time_order 2): h_1-h_3, h_4, h_7, h_9 —
          time-space and diagonal spatial components
        """
        orders = {eq.field_name: eq.time_derivative_order for eq in spec.equations}
        unique_orders = set(orders.values())
        # Should have at least two different time orders (0 and 2)
        assert len(unique_orders) >= 2, (
            f"Expected mixed time orders, got only {unique_orders}"
        )
        # h_0 (Hamiltonian constraint) should be elliptic (order 0)
        assert orders["h_0"] == 0, (
            f"Expected h_0 (Hamiltonian constraint) to have time_order 0, "
            f"got {orders['h_0']}"
        )
        # Evolution equations should be second-order
        assert 2 in unique_orders, "Expected some second-order (evolution) equations"
        # Diagonal spatial components h_4 (h_xx), h_7 (h_yy), h_9 (h_zz) are evolution
        for i in (4, 7, 9):
            assert orders[f"h_{i}"] == 2, (
                f"Expected h_{i} (evolution equation) to have time_order 2, "
                f"got {orders[f'h_{i}']}"
            )

    def test_cross_field_coupling_exists(self, spec: EquationSystem) -> None:
        """Gauge-unfixed equations should have cross-field coupling."""
        has_cross_coupling = False
        for eq in spec.equations:
            for term in eq.rhs_terms:
                # Check for reference to different field (not pi_ momentum)
                if (
                    term.field != eq.field_name
                    and not term.field.startswith("pi_")
                ):
                    has_cross_coupling = True
                    break
            if has_cross_coupling:
                break

        assert has_cross_coupling, (
            "Expected cross-field coupling in gauge-unfixed equations"
        )

    def test_momentum_references_exist(self, spec: EquationSystem) -> None:
        """Should reference momentum fields (pi_*) for mixed time-space derivatives."""
        has_momentum_ref = False
        for eq in spec.equations:
            for term in eq.rhs_terms:
                if term.field.startswith("pi_"):
                    has_momentum_ref = True
                    break
            if has_momentum_ref:
                break

        assert has_momentum_ref, "Expected momentum field references (pi_*)"

    def test_diverse_operators(self, spec: EquationSystem) -> None:
        """Gauge-unfixed equations should use multiple operator types."""
        all_operators: set[str] = set()
        for eq in spec.equations:
            for term in eq.rhs_terms:
                all_operators.add(term.operator)

        # Should have directional laplacians, gradients, cross-derivatives
        assert "laplacian_x" in all_operators
        assert "laplacian_y" in all_operators
        assert "laplacian_z" in all_operators
        # Should have gradients (from mixed time-space terms)
        has_gradient = any(op.startswith("gradient_") for op in all_operators)
        assert has_gradient, f"Expected gradient operators, got {all_operators}"
        # Should have cross-derivatives (from cross terms in linearized Einstein)
        has_cross = any(op.startswith("cross_derivative_") for op in all_operators)
        assert has_cross, f"Expected cross_derivative operators, got {all_operators}"

    def test_first_derivative_t_exists(self, spec: EquationSystem) -> None:
        """Some terms should have first_derivative_t operator (Hubble-like friction)."""
        has_fdt = False
        for eq in spec.equations:
            for term in eq.rhs_terms:
                if term.operator == "first_derivative_t":
                    has_fdt = True
                    break
            if has_fdt:
                break
        assert has_fdt, "Expected first_derivative_t terms in gauge-unfixed equations"

    def test_state_size_reflects_mixed_orders(self, spec: EquationSystem) -> None:
        """State size should be less than 20 due to constraint/first-order components."""
        # h_0 is elliptic (order 0, 1 slot), some may be first-order (1 slot each),
        # evolution eqs are second-order (2 slots each)
        # So state_size < 2 * n_components
        assert spec.state_size < 2 * spec.n_components, (
            f"Expected state_size < {2 * spec.n_components} due to non-second-order components, "
            f"got {spec.state_size}"
        )

    def test_metadata(self, spec: EquationSystem) -> None:
        assert spec.metadata.get("gauge") == "none"
        assert spec.metadata.get("linearized") is True

    def test_pde_builds(self) -> None:
        """PDE should build from gauge-unfixed JSON without errors."""
        from torsion_gertsenshtein.symbolic.pde_builder import build_pde_from_json
        pde = build_pde_from_json(UNFIXED_JSON)
        assert isinstance(pde, PDEFromSpec)


# ============================================================
# Physics validation tests (de Donder)
# ============================================================


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
            [(-1, 1), (-1, 1), (-10, 10)], [4, 4, 64], periodic=True,
        )

    def test_tt_compatible_initial_data(
        self, spec: EquationSystem, grid: CartesianGrid,
    ) -> None:
        """TT-gauge compatible initial data should propagate correctly.

        For a GW propagating along z:
        - h_4 (h_{xx}) and h_7 (h_{yy}) are the + polarization: h_+ = h_{xx} = -h_{yy}
        - h_5 (h_{xy}) is the x polarization
        - All other components zero (TT gauge)
        """
        z = grid.cell_coords[..., 2]
        # Gaussian wave packet for h_+ polarization
        h_plus = 0.01 * np.exp(-z**2 / 2.0)

        state = create_initial_state(
            grid, spec,
            field_data={"h_4": h_plus, "h_7": -h_plus},
        )

        pde = PDEFromSpec(spec)
        result = pde.solve(state, t_range=0.5, dt=0.01, scheme="runge-kutta")

        # After propagation, h_4 and h_7 should still be equal and opposite
        # (tracelessness maintained by the wave equation)
        h4_idx = spec.state_layout.index(("h_4", "field"))
        h7_idx = spec.state_layout.index(("h_7", "field"))
        h4_data = result[h4_idx].data
        h7_data = result[h7_idx].data

        # h_4 = -h_7 should be maintained (TT tracelessness)
        assert_allclose(
            h4_data, -h7_data, atol=1e-6,
            err_msg="TT tracelessness h_{xx} = -h_{yy} not maintained",
        )

    def test_energy_conservation_short_time(
        self, spec: EquationSystem,
    ) -> None:
        """Total field energy should be approximately conserved."""
        grid = CartesianGrid(
            [(-5, 5), (-5, 5), (-5, 5)], [8, 8, 8], periodic=True,
        )
        pde = PDEFromSpec(spec)

        z = grid.cell_coords[..., 2]
        h_plus = 0.01 * np.exp(-z**2 / 2.0)
        state = create_initial_state(
            grid, spec,
            field_data={"h_4": h_plus},
        )

        # Compute initial "energy" (sum of field² + momentum²)
        initial_energy = sum(np.sum(f.data**2) for f in state)

        result = pde.solve(state, t_range=0.2, dt=0.05, scheme="runge-kutta")
        final_energy = sum(np.sum(f.data**2) for f in result)

        # Energy should be conserved to within ~15% on this coarse 8³ grid
        rel_change = abs(final_energy - initial_energy) / max(initial_energy, 1e-15)
        assert rel_change < 0.15, (
            f"Energy changed by {rel_change:.1%} (initial={initial_energy:.6f}, "
            f"final={final_energy:.6f})"
        )
