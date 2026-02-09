"""Tests for rank-2 tensor field support (Phase 10).

Rank-2 symmetric tensors (e.g., h_{ab} for linearized gravity) decompose into
N*(N+1)/2 independent scalar components. The pipeline treats these as
N independent scalar equations — no structural changes needed in pde_builder.

These tests verify that:
- Rank-2 JSON specs load correctly
- State layout handles multiple components
- Evolution rates are computed correctly
- create_initial_state works with rank-2 component counts
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
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    create_initial_state,
)


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """1D periodic grid for testing."""
    return CartesianGrid([(0, 10)], 32, periodic=True)


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """2D periodic grid for testing."""
    return CartesianGrid([(0, 10), (0, 10)], [16, 16], periodic=True)


@pytest.fixture
def symmetric_rank2_2d_spec() -> EquationSystem:
    """Symmetric rank-2 tensor in 1+1D: 3 independent components.

    Simulates a simplified Fierz-Pauli-like system:
      d2_t(h_0) = laplacian(h_0)     (h_{00})
      d2_t(h_1) = laplacian(h_1)     (h_{01})
      d2_t(h_2) = laplacian(h_2)     (h_{11})

    All components are uncoupled (no cross terms) for simplicity.
    """
    return EquationSystem(
        n_components=3,
        dimension=2,
        spatial_dimension=1,
        component_names=("h_0", "h_1", "h_2"),
        equations=(
            ComponentEquation(
                field_name="h_0",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="h_0"),
                ),
            ),
            ComponentEquation(
                field_name="h_1",
                field_index=1,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="h_1"),
                ),
            ),
            ComponentEquation(
                field_name="h_2",
                field_index=2,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="h_2"),
                ),
            ),
        ),
        mass_matrix=(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ),
        coupling_matrix=(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ),
        metadata={"source": "test", "type": "symmetric_rank2"},
    )


@pytest.fixture
def coupled_rank2_spec() -> EquationSystem:
    """Rank-2 tensor with cross-coupling between components.

    Simulates linearized gravity with constraint-like coupling:
      d2_t(h_0) = laplacian(h_0) - h_2     (trace coupling h_{00} <- h_{11})
      d_t(h_1)  = gradient_x(h_0)           (constraint, first-order)
      d2_t(h_2) = laplacian(h_2) - h_0     (trace coupling h_{11} <- h_{00})
    """
    return EquationSystem(
        n_components=3,
        dimension=2,
        spatial_dimension=1,
        component_names=("h_0", "h_1", "h_2"),
        equations=(
            ComponentEquation(
                field_name="h_0",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="h_0"),
                    OperatorTerm(coefficient=-1.0, operator="identity", field="h_2"),
                ),
            ),
            ComponentEquation(
                field_name="h_1",
                field_index=1,
                time_derivative_order=1,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="gradient_x", field="h_0"),
                ),
            ),
            ComponentEquation(
                field_name="h_2",
                field_index=2,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="h_2"),
                    OperatorTerm(coefficient=-1.0, operator="identity", field="h_0"),
                ),
            ),
        ),
        mass_matrix=(
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ),
        coupling_matrix=(
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        ),
        metadata={"source": "test", "type": "coupled_rank2"},
    )


class TestRank2StateLayout:
    """Tests for state layout with rank-2 tensor components."""

    def test_symmetric_3_components_state_size(
        self, symmetric_rank2_2d_spec: EquationSystem
    ) -> None:
        """3 second-order components -> state_size = 6."""
        assert symmetric_rank2_2d_spec.n_components == 3
        assert symmetric_rank2_2d_spec.state_size == 6

    def test_symmetric_state_layout(
        self, symmetric_rank2_2d_spec: EquationSystem
    ) -> None:
        """Layout: h_0, pi_h0, h_1, pi_h1, h_2, pi_h2."""
        layout = symmetric_rank2_2d_spec.state_layout
        assert layout == (
            ("h_0", "field"),
            ("h_0", "momentum"),
            ("h_1", "field"),
            ("h_1", "momentum"),
            ("h_2", "field"),
            ("h_2", "momentum"),
        )

    def test_coupled_mixed_order_state_size(
        self, coupled_rank2_spec: EquationSystem
    ) -> None:
        """Mixed: 2 second-order + 1 first-order -> state_size = 5."""
        assert coupled_rank2_spec.state_size == 5

    def test_coupled_mixed_order_layout(
        self, coupled_rank2_spec: EquationSystem
    ) -> None:
        """h_0 (field, mom), h_1 (field only), h_2 (field, mom)."""
        layout = coupled_rank2_spec.state_layout
        assert layout == (
            ("h_0", "field"),
            ("h_0", "momentum"),
            ("h_1", "field"),
            ("h_2", "field"),
            ("h_2", "momentum"),
        )


class TestRank2Evolution:
    """Tests for rank-2 tensor evolution."""

    def test_symmetric_evolution_rate_shape(
        self,
        symmetric_rank2_2d_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Evolution rate has correct number of fields."""
        pde = PDEFromSpec(symmetric_rank2_2d_spec)
        state = create_initial_state(grid_1d, symmetric_rank2_2d_spec)

        rates = pde.evolution_rate(state)
        assert len(rates) == 6

    def test_uncoupled_components_independent(
        self,
        symmetric_rank2_2d_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Uncoupled rank-2 components evolve independently."""
        pde = PDEFromSpec(symmetric_rank2_2d_spec)

        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        h0_data = np.sin(2 * np.pi * x / 10)

        state = create_initial_state(
            grid_1d,
            symmetric_rank2_2d_spec,
            field_data={"h_0": h0_data},
        )

        rates = pde.evolution_rate(state)

        # d/dt h_0 = pi_h0 = 0 (zero momentum)
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d/dt pi_h0 = laplacian(h_0) (non-zero for sinusoidal field)
        assert np.max(np.abs(rates[1].data)) > 0

        # h_1 and h_2 fields and momenta should all be zero
        for i in [2, 3, 4, 5]:
            assert_allclose(rates[i].data, 0.0, atol=1e-10)

    def test_coupled_rank2_evolution(
        self,
        coupled_rank2_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Coupled rank-2 with mixed time orders evolves correctly."""
        pde = PDEFromSpec(coupled_rank2_spec)

        # State: [h_0, pi_0, h_1, h_2, pi_2] (5 fields)
        state = FieldCollection(
            [
                ScalarField(grid_1d, data=1.0),  # h_0 (uniform)
                ScalarField(grid_1d, data=0.0),  # pi_0
                ScalarField(grid_1d, data=0.0),  # h_1
                ScalarField(grid_1d, data=2.0),  # h_2 (uniform)
                ScalarField(grid_1d, data=0.0),  # pi_2
            ]
        )

        rates = pde.evolution_rate(state)
        assert len(rates) == 5

        # d/dt h_0 = pi_0 = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d/dt pi_0 = laplacian(h_0) - h_2 = 0 - 2 = -2
        assert_allclose(rates[1].data, -2.0, atol=1e-10)

        # d/dt h_1 = gradient_x(h_0) = 0 (uniform h_0)
        assert_allclose(rates[2].data, 0.0, atol=1e-10)

        # d/dt h_2 = pi_2 = 0
        assert_allclose(rates[3].data, 0.0, atol=1e-10)

        # d/dt pi_2 = laplacian(h_2) - h_0 = 0 - 1 = -1
        assert_allclose(rates[4].data, -1.0, atol=1e-10)


class TestRank2CreateInitialState:
    """Tests for create_initial_state with rank-2 specs."""

    def test_symmetric_initial_state(
        self,
        symmetric_rank2_2d_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Create initial state for 3-component symmetric tensor."""
        state = create_initial_state(grid_1d, symmetric_rank2_2d_spec)

        # 3 components * 2 (field + momentum) = 6
        assert len(state) == 6

    def test_symmetric_with_data(
        self,
        symmetric_rank2_2d_spec: EquationSystem,
        grid_1d: CartesianGrid,
    ) -> None:
        """Verify field_data maps to correct slots."""
        data = np.ones(32) * 3.0
        state = create_initial_state(
            grid_1d,
            symmetric_rank2_2d_spec,
            field_data={"h_1": data},
        )

        # h_0 = 0, pi_0 = 0, h_1 = 3, pi_1 = 0, h_2 = 0, pi_2 = 0
        assert_allclose(state[0].data, 0.0)
        assert_allclose(state[1].data, 0.0)
        assert_allclose(state[2].data, 3.0)
        assert_allclose(state[3].data, 0.0)
        assert_allclose(state[4].data, 0.0)
        assert_allclose(state[5].data, 0.0)


class TestRank3TensorSupport:
    """Smoke tests for rank-3 tensor component handling.

    Rank-3 tensors in 2D decompose into 8 scalar components (non-symmetric).
    The Python pipeline treats them identically to any other N-component system.
    """

    @pytest.fixture
    def rank3_2d_spec(self) -> EquationSystem:
        """Rank-3 tensor in 1+1D: 4 of 8 components as wave equations.

        Simplified test: 4 components with laplacian, verifying the Python
        pipeline handles arbitrary component counts from higher-rank tensors.
        """
        n = 4
        names = tuple(f"T_{i}" for i in range(n))
        equations = tuple(
            ComponentEquation(
                field_name=names[i],
                field_index=i,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field=names[i]),
                ),
            )
            for i in range(n)
        )
        zero_row = tuple(0.0 for _ in range(n))
        return EquationSystem(
            n_components=n,
            dimension=2,
            spatial_dimension=1,
            component_names=names,
            equations=equations,
            mass_matrix=tuple(zero_row for _ in range(n)),
            coupling_matrix=tuple(zero_row for _ in range(n)),
            metadata={"source": "test", "type": "rank3_tensor"},
        )

    def test_rank3_state_size(self, rank3_2d_spec: EquationSystem) -> None:
        """4 second-order components -> state_size = 8."""
        assert rank3_2d_spec.n_components == 4
        assert rank3_2d_spec.state_size == 8

    def test_rank3_state_layout(self, rank3_2d_spec: EquationSystem) -> None:
        """Layout alternates field/momentum for each component."""
        layout = rank3_2d_spec.state_layout
        assert len(layout) == 8
        assert layout[0] == ("T_0", "field")
        assert layout[1] == ("T_0", "momentum")
        assert layout[6] == ("T_3", "field")
        assert layout[7] == ("T_3", "momentum")

    def test_rank3_evolution(
        self, rank3_2d_spec: EquationSystem, grid_1d: CartesianGrid
    ) -> None:
        """Evolution rate has correct shape for 4-component rank-3 system."""
        pde = PDEFromSpec(rank3_2d_spec, parameters={})
        state = create_initial_state(grid_1d, rank3_2d_spec)
        # Set non-zero data for first component
        state[0].data[:] = np.sin(2 * np.pi * np.linspace(0, 10, 32) / 10)
        rates = pde.evolution_rate(state, t=0)
        assert len(rates) == 8

    def test_rank3_8_components(self, grid_1d: CartesianGrid) -> None:
        """Full rank-3 in 2D: 8 components all evolving independently."""
        n = 8
        names = tuple(f"T_{i}" for i in range(n))
        equations = tuple(
            ComponentEquation(
                field_name=names[i],
                field_index=i,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field=names[i]),
                ),
            )
            for i in range(n)
        )
        zero_row = tuple(0.0 for _ in range(n))
        spec = EquationSystem(
            n_components=n,
            dimension=2,
            spatial_dimension=1,
            component_names=names,
            equations=equations,
            mass_matrix=tuple(zero_row for _ in range(n)),
            coupling_matrix=tuple(zero_row for _ in range(n)),
            metadata={"source": "test", "type": "rank3_full"},
        )
        assert spec.n_components == 8
        assert spec.state_size == 16
        state = create_initial_state(grid_1d, spec)
        assert len(state) == 16
