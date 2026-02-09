"""Tests for mixed time-derivative order support (Phase 9).

Tests that the pipeline correctly handles:
- First-order equations (heat/diffusion, time_order=1)
- Constraint equations (elliptic, time_order=0)
- Mixed systems (some second-order, some first-order)
- Virtual momentum resolution for cross-references
- Backward compatibility with all-second-order systems
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
)
from torsion_gertsenshtein.symbolic.pde_builder import (
    PDEFromSpec,
    create_initial_state,
)
from torsion_gertsenshtein.utils import normalize_solve_result

# Tests need to access private methods


# === Fixtures ===


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """1D periodic grid for testing."""
    return CartesianGrid([(0, 10)], 64, periodic=True)


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """2D periodic grid for testing."""
    return CartesianGrid([(0, 10), (0, 10)], [32, 32], periodic=True)


@pytest.fixture
def heat_spec() -> EquationSystem:
    """First-order heat equation: d_t(phi) = laplacian(phi).

    Single component, first-order in time. State has 1 slot (no momentum).
    """
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test", "type": "heat"},
    )


@pytest.fixture
def constraint_spec() -> EquationSystem:
    """Constraint equation: 0 = laplacian(phi) + source.

    Single component, zeroth-order in time. State has 1 slot, no evolution.
    """
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=0,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test", "type": "constraint"},
    )


@pytest.fixture
def mixed_spec() -> EquationSystem:
    """Mixed system: one second-order component + one first-order component.

    phi: d2_t(phi) = laplacian(phi)        (wave, second-order)
    chi: d_t(chi) = laplacian(chi) + phi   (diffusion + forcing, first-order)

    State layout: [phi, pi_phi, chi] (3 slots total).
    """
    return EquationSystem(
        n_components=2,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi", "chi"),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                ),
            ),
            ComponentEquation(
                field_name="chi",
                field_index=1,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="chi"),
                    OperatorTerm(coefficient=1.0, operator="identity", field="phi"),
                ),
            ),
        ),
        mass_matrix=((0.0, 0.0), (0.0, 0.0)),
        coupling_matrix=((0.0, 0.0), (-1.0, 0.0)),
        metadata={"source": "test", "type": "mixed"},
    )


@pytest.fixture
def mixed_with_virtual_momentum_spec() -> EquationSystem:
    """Mixed system where second-order equation references first-order momentum.

    A_0: d_t(A_0) = gradient_x(A_1)                   (first-order)
    A_1: d2_t(A_1) = laplacian(A_1) + gradient_x(pi_0) (second-order, refs pi_0)

    pi_0 = d_t(A_0) is a virtual momentum (computed from A_0's RHS).
    State layout: [A_0, A_1, pi_1] (3 slots total).
    """
    return EquationSystem(
        n_components=2,
        dimension=2,
        spatial_dimension=1,
        component_names=("A_0", "A_1"),
        equations=(
            ComponentEquation(
                field_name="A_0",
                field_index=0,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="gradient_x", field="A_1"),
                ),
            ),
            ComponentEquation(
                field_name="A_1",
                field_index=1,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_1"),
                    OperatorTerm(coefficient=0.5, operator="gradient_x", field="pi_0"),
                ),
            ),
        ),
        mass_matrix=((0.0, 0.0), (0.0, 0.0)),
        coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
        metadata={"source": "test"},
    )


# === EquationSystem Property Tests ===


class TestEquationSystemStateProperties:
    """Tests for state_size, state_layout, and time_orders properties."""

    def test_all_second_order(self) -> None:
        """All second-order: state_size = 2 * n_components (backward compat)."""
        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("A_0", "A_1"),
            equations=(
                ComponentEquation(
                    "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                ),
                ComponentEquation(
                    "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )

        assert spec.time_orders == (2, 2)
        assert spec.state_size == 4
        assert spec.state_layout == (
            ("A_0", "field"),
            ("A_0", "momentum"),
            ("A_1", "field"),
            ("A_1", "momentum"),
        )

    def test_single_first_order(self, heat_spec: EquationSystem) -> None:
        """Single first-order: state_size = 1 (no momentum)."""
        assert heat_spec.time_orders == (1,)
        assert heat_spec.state_size == 1
        assert heat_spec.state_layout == (("phi", "field"),)

    def test_single_constraint(self, constraint_spec: EquationSystem) -> None:
        """Constraint (order=0): state_size = 1 (no momentum)."""
        assert constraint_spec.time_orders == (0,)
        assert constraint_spec.state_size == 1
        assert constraint_spec.state_layout == (("phi", "field"),)

    def test_mixed_orders(self, mixed_spec: EquationSystem) -> None:
        """Mixed (2nd + 1st): state_size = 3 (field + momentum + field)."""
        assert mixed_spec.time_orders == (2, 1)
        assert mixed_spec.state_size == 3
        assert mixed_spec.state_layout == (
            ("phi", "field"),
            ("phi", "momentum"),
            ("chi", "field"),
        )


# === PDEFromSpec Slot Map Tests ===


class TestSlotMaps:
    """Tests for _field_slot_map and _momentum_slot_map."""

    def test_all_second_order_slots(self) -> None:
        """Second-order: alternating field/momentum slots."""
        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("A_0", "A_1"),
            equations=(
                ComponentEquation(
                    "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                ),
                ComponentEquation(
                    "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )
        pde = PDEFromSpec(spec)

        assert pde._field_slot_map == {"A_0": 0, "A_1": 2}
        assert pde._momentum_slot_map == {"A_0": 1, "A_1": 3}

    def test_first_order_no_momentum_slot(self, heat_spec: EquationSystem) -> None:
        """First-order: no momentum slot."""
        pde = PDEFromSpec(heat_spec)

        assert pde._field_slot_map == {"phi": 0}
        assert pde._momentum_slot_map == {}

    def test_mixed_slots(self, mixed_spec: EquationSystem) -> None:
        """Mixed: phi has momentum, chi does not."""
        pde = PDEFromSpec(mixed_spec)

        assert pde._field_slot_map == {"phi": 0, "chi": 2}
        assert pde._momentum_slot_map == {"phi": 1}


# === Evolution Tests ===


class TestFirstOrderEvolution:
    """Tests for first-order (parabolic) equations."""

    def test_heat_equation_state_size(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Heat equation: state has 1 field (no momentum)."""
        state = create_initial_state(grid_1d, heat_spec)

        assert len(state) == 1

    def test_heat_equation_rate_size(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Heat equation: evolution_rate returns 1 field."""
        pde = PDEFromSpec(heat_spec)
        state = FieldCollection([ScalarField(grid_1d, data=1.0)])

        rates = pde.evolution_rate(state)

        assert len(rates) == 1

    def test_heat_equation_rate_is_laplacian(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Heat equation: d_t(phi) = laplacian(phi)."""
        pde = PDEFromSpec(heat_spec)

        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        kx = 2 * np.pi / 10.0
        phi_data = np.sin(kx * x)
        state = FieldCollection([ScalarField(grid_1d, data=phi_data)])

        rates = pde.evolution_rate(state)

        # d_t(phi) = laplacian(phi) = -kx^2 sin(kx*x)
        expected = -(kx**2) * np.sin(kx * x)
        assert_allclose(rates[0].data, expected, rtol=0.05)

    def test_heat_equation_diffusion(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Heat equation: Gaussian pulse should spread over time."""
        pde = PDEFromSpec(heat_spec)

        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        phi_data = np.exp(-((x - 5) ** 2))
        initial_peak = float(np.max(phi_data))
        state = FieldCollection([ScalarField(grid_1d, data=phi_data)])

        result = pde.solve(state, t_range=1.0, dt=0.01)
        sol = normalize_solve_result(result)

        # Peak should decrease (diffusion spreads the pulse)
        final_peak = float(np.max(sol[0].data))
        assert final_peak < initial_peak

    def test_pi_reference_for_first_order_raises(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Referencing pi_0 for a first-order field (without virtual momentum) raises."""
        pde = PDEFromSpec(heat_spec)
        state = FieldCollection([ScalarField(grid_1d, data=1.0)])

        # Direct lookup of pi_0 without virtual momenta should raise
        with pytest.raises(ValueError, match="no momentum in state"):
            pde._get_field_from_state(state, "pi_0")


class TestConstraintEvolution:
    """Tests for constraint (order=0) equations."""

    def test_constraint_state_size(
        self, constraint_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Constraint: state has 1 field."""
        state = create_initial_state(grid_1d, constraint_spec)
        assert len(state) == 1

    def test_constraint_rate_is_zero(
        self, constraint_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Constraint: evolution rate is zero (no time evolution)."""
        pde = PDEFromSpec(constraint_spec)
        state = FieldCollection([ScalarField(grid_1d, data=1.0)])

        rates = pde.evolution_rate(state)
        assert len(rates) == 1
        assert_allclose(rates[0].data, 0.0)

    def test_constraint_eps_configurable(self) -> None:
        """constraint_eps parameter controls the Laplacian zero-check threshold."""
        # Equation with a very small Laplacian coefficient
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=0,
                    rhs_terms=(
                        OperatorTerm(
                            coefficient=1e-10,
                            operator="laplacian",
                            field="phi",
                        ),
                    ),
                    constraint_solver=ConstraintSolverConfig(enabled=True),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={"source": "test"},
        )
        # Default eps=1e-14 should accept coeff=1e-10
        pde_ok = PDEFromSpec(spec)
        assert pde_ok._constraint_eps == 1e-14

        # Tight eps=1e-8 should reject coeff=1e-10
        import pytest

        pde_reject = PDEFromSpec(spec, constraint_eps=1e-8)
        grid = CartesianGrid([(0, 10)], 16, periodic=True)
        state = FieldCollection([ScalarField(grid, data=1.0)])
        with pytest.raises(ValueError, match="effectively zero"):
            pde_reject.evolution_rate(state)


class TestMixedOrderEvolution:
    """Tests for mixed time-derivative order systems."""

    def test_mixed_state_size(
        self, mixed_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Mixed (2nd + 1st): state has 3 slots."""
        state = create_initial_state(grid_1d, mixed_spec)
        num_mixed_slots = 3
        assert len(state) == num_mixed_slots

    def test_mixed_rate_size(
        self, mixed_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Mixed system: evolution_rate returns 3 rates."""
        pde = PDEFromSpec(mixed_spec)
        state = FieldCollection(
            [
                ScalarField(grid_1d, data=1.0),  # phi
                ScalarField(grid_1d, data=0.0),  # pi_phi
                ScalarField(grid_1d, data=0.0),  # chi
            ]
        )

        rates = pde.evolution_rate(state)
        num_mixed_slots = 3
        assert len(rates) == num_mixed_slots

    def test_mixed_wave_component_evolves(
        self, mixed_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Mixed: second-order phi component has d/dt phi = pi_phi."""
        pde = PDEFromSpec(mixed_spec)

        pi_value = 2.0
        state = FieldCollection(
            [
                ScalarField(grid_1d, data=0.0),  # phi
                ScalarField(grid_1d, data=pi_value),  # pi_phi
                ScalarField(grid_1d, data=0.0),  # chi
            ]
        )

        rates = pde.evolution_rate(state)

        # d/dt phi = pi_phi = 2.0
        assert_allclose(rates[0].data, pi_value, atol=1e-10)

    def test_mixed_first_order_component_evolves(
        self, mixed_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Mixed: first-order chi has d_t(chi) = laplacian(chi) + phi."""
        pde = PDEFromSpec(mixed_spec)

        # phi = 3.0 uniform, chi = 0 (uniform -> laplacian = 0)
        state = FieldCollection(
            [
                ScalarField(grid_1d, data=3.0),  # phi
                ScalarField(grid_1d, data=0.0),  # pi_phi
                ScalarField(grid_1d, data=0.0),  # chi (uniform)
            ]
        )

        rates = pde.evolution_rate(state)

        # d/dt chi = laplacian(0) + 1.0 * phi = 0 + 3.0 = 3.0
        assert_allclose(rates[2].data, 3.0, atol=1e-10)


class TestVirtualMomentum:
    """Tests for virtual momentum resolution (pi_N for first-order components)."""

    def test_virtual_momentum_in_evolution(
        self,
        mixed_with_virtual_momentum_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Second-order equation accessing pi_0 of first-order component."""
        pde = PDEFromSpec(mixed_with_virtual_momentum_spec)

        # State: [A_0, A_1, pi_1] (3 fields)
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        kx = 2 * np.pi / 10.0
        a1_data = np.sin(kx * x)

        state = FieldCollection(
            [
                ScalarField(grid_1d, data=0.0),  # A_0
                ScalarField(grid_1d, data=a1_data),  # A_1
                ScalarField(grid_1d, data=0.0),  # pi_1
            ]
        )

        # Should not raise — pi_0 is resolved as virtual momentum from A_0's RHS
        rates = pde.evolution_rate(state)

        num_mixed_slots = 3
        assert len(rates) == num_mixed_slots

        # A_0 equation: d_t(A_0) = gradient_x(A_1) = kx * cos(kx*x)
        # This is pi_0 (virtual momentum)
        expected_pi0 = kx * np.cos(kx * x)
        assert_allclose(rates[0].data, expected_pi0, rtol=0.05)

        # A_1 equation: d_t(pi_1) = laplacian(A_1) + 0.5 * gradient_x(pi_0)
        # laplacian(sin(kx*x)) = -kx^2 * sin(kx*x)
        # gradient_x(pi_0) = gradient_x(kx * cos(kx*x)) = -kx^2 * sin(kx*x)
        # Total: -kx^2 sin(kx*x) + 0.5 * (-kx^2 sin(kx*x)) = -1.5 kx^2 sin(kx*x)
        expected_dpi1 = -1.5 * kx**2 * np.sin(kx * x)
        assert_allclose(rates[2].data, expected_dpi1, rtol=0.05)

    def test_virtual_momentum_via_get_field_from_state(
        self,
        mixed_with_virtual_momentum_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """_get_field_from_state resolves pi_0 from virtual_momenta dict."""
        pde = PDEFromSpec(mixed_with_virtual_momentum_spec)

        state = FieldCollection(
            [
                ScalarField(grid_1d, data=0.0),
                ScalarField(grid_1d, data=0.0),
                ScalarField(grid_1d, data=0.0),
            ]
        )

        # Without virtual momenta, pi_0 should raise
        with pytest.raises(ValueError, match="no momentum in state"):
            pde._get_field_from_state(state, "pi_0")

        # With virtual momenta, pi_0 should resolve
        virtual = {"A_0": ScalarField(grid_1d, data=42.0)}
        result = pde._get_field_from_state(state, "pi_0", virtual)
        assert_allclose(result.data, 42.0)


# === create_initial_state Tests ===


class TestCreateInitialStateMixed:
    """Tests for create_initial_state with mixed time orders."""

    def test_first_order_no_momentum_slot(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """First-order: create_initial_state produces 1 field (no momentum)."""
        state = create_initial_state(grid_1d, heat_spec)
        assert len(state) == 1
        assert_allclose(state[0].data, 0.0)

    def test_mixed_state_layout(
        self, mixed_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Mixed: field_data and momentum_data map to correct slots."""
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        phi_data = np.sin(2 * np.pi * x / 10)
        chi_data = np.cos(2 * np.pi * x / 10)

        state = create_initial_state(
            grid_1d,
            mixed_spec,
            field_data={"phi": phi_data, "chi": chi_data},
            momentum_data={"phi": np.ones(64) * 5.0},
        )

        # State contains 3 slots: phi, pi_phi, and chi
        num_mixed_slots = 3
        assert len(state) == num_mixed_slots
        assert_allclose(state[0].data, phi_data)
        assert_allclose(state[1].data, 5.0)
        assert_allclose(state[2].data, chi_data)

    def test_momentum_data_ignored_for_first_order(
        self, heat_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Momentum data for first-order components is ignored (no momentum slot)."""
        state = create_initial_state(
            grid_1d,
            heat_spec,
            momentum_data={"phi": np.ones(64) * 99.0},
        )

        # Only 1 field, no momentum slot
        assert len(state) == 1
        assert_allclose(state[0].data, 0.0)  # Field defaults to zero


# === Backward Compatibility Tests ===


class TestBackwardCompatibility:
    """Verify all-second-order systems still work identically."""

    def test_kg_state_size(self) -> None:
        """KG spec: state_size == 2 * n_components."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    "phi",
                    0,
                    2,
                    (
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(-1.0, "identity", "phi"),
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        assert spec.state_size == 2
        assert spec.state_layout == (("phi", "field"), ("phi", "momentum"))

    def test_em_state_size(self) -> None:
        """EM spec: state_size == 4 (2 components * 2)."""
        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("A_0", "A_1"),
            equations=(
                ComponentEquation(
                    "A_0", 0, 2, (OperatorTerm(1.0, "laplacian", "A_0"),)
                ),
                ComponentEquation(
                    "A_1", 1, 2, (OperatorTerm(1.0, "laplacian", "A_1"),)
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
        )

        assert spec.state_size == 4

    def test_create_initial_state_backward_compat(self, grid_1d: CartesianGrid) -> None:
        """create_initial_state produces [field, momentum] pairs for second-order."""
        spec = EquationSystem(
            n_components=1,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    "phi", 0, 2, (OperatorTerm(1.0, "laplacian", "phi"),)
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        state = create_initial_state(
            grid_1d,
            spec,
            field_data={"phi": np.ones(64)},
            momentum_data={"phi": np.ones(64) * 2.0},
        )

        assert len(state) == 2
        assert_allclose(state[0].data, 1.0)
        assert_allclose(state[1].data, 2.0)
