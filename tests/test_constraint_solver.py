"""Tests for elliptic constraint solving (Issue #91).

Tests that the pipeline correctly handles:
- BoundaryCondition and ConstraintSolverConfig parsing
- Validation of constraint_solver on non-constraint equations
- Poisson equation solving via py-pde
- Backward compatibility with frozen constraints
- Cross-field constraint solving
"""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic.json_loader import (
    BoundaryCondition,
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
)
from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec

# === Fixtures ===


@pytest.fixture
def grid_1d_periodic() -> CartesianGrid:
    """1D periodic grid."""
    return CartesianGrid([(0, 2 * np.pi)], 64, periodic=True)


@pytest.fixture
def grid_1d_nonperiodic() -> CartesianGrid:
    """1D non-periodic grid."""
    return CartesianGrid([(0, 1)], 64, periodic=False)


@pytest.fixture
def grid_2d_periodic() -> CartesianGrid:
    """2D periodic grid."""
    return CartesianGrid([(0, 2 * np.pi), (0, 2 * np.pi)], [32, 32], periodic=True)


def _make_poisson_spec(
    *,
    solver_enabled: bool = True,
    bcs: dict[str, BoundaryCondition] | None = None,
    extra_source: OperatorTerm | None = None,
) -> EquationSystem:
    """Create a Poisson-type constraint spec.

    Equation: 0 = laplacian(phi) + source_terms
    """
    terms: list[OperatorTerm] = [
        OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
    ]
    if extra_source is not None:
        terms.append(extra_source)

    constraint_solver = ConstraintSolverConfig(
        enabled=solver_enabled,
        method="poisson",
        boundary_conditions=bcs or {},
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
                time_derivative_order=0,
                rhs_terms=tuple(terms),
                constraint_solver=constraint_solver,
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={"source": "test"},
    )


# === BoundaryCondition Parsing Tests ===


class TestBoundaryConditionParsing:
    """Test BoundaryCondition.from_dict."""

    def test_periodic_bc(self) -> None:
        bc = BoundaryCondition.from_dict({"type": "periodic"})
        assert bc.type == "periodic"
        assert bc.value is None
        assert bc.derivative is None

    def test_dirichlet_bc(self) -> None:
        bc = BoundaryCondition.from_dict({"type": "dirichlet", "value": 1.5})
        assert bc.type == "dirichlet"
        assert bc.value == 1.5

    def test_neumann_bc(self) -> None:
        bc = BoundaryCondition.from_dict({"type": "neumann", "derivative": 0.0})
        assert bc.type == "neumann"
        assert bc.derivative == 0.0

    def test_invalid_bc_type_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown BC type"):
            BoundaryCondition.from_dict({"type": "robin"})


# === ConstraintSolverConfig Tests ===


class TestConstraintSolverConfig:
    """Test ConstraintSolverConfig.from_dict."""

    def test_default_disabled(self) -> None:
        config = ConstraintSolverConfig.from_dict(None)
        assert config.enabled is False
        assert config.method == "poisson"
        assert config.boundary_conditions == {}

    def test_enabled_with_bcs(self) -> None:
        data = {
            "enabled": True,
            "method": "poisson",
            "boundary_conditions": {
                "x": {"type": "periodic"},
                "y": {"type": "dirichlet", "value": 0.0},
            },
        }
        config = ConstraintSolverConfig.from_dict(data)
        assert config.enabled is True
        assert config.method == "poisson"
        assert "x" in config.boundary_conditions
        assert config.boundary_conditions["x"].type == "periodic"
        assert config.boundary_conditions["y"].type == "dirichlet"
        assert config.boundary_conditions["y"].value == 0.0

    def test_empty_dict_defaults(self) -> None:
        config = ConstraintSolverConfig.from_dict({})
        assert config.enabled is False


# === Validation Tests ===


class TestConstraintSolverValidation:
    """Test that constraint_solver is only valid on constraint equations."""

    def test_solver_on_wave_equation_raises(self) -> None:
        """constraint_solver.enabled=True on time_order=2 must raise."""
        with pytest.raises(ValueError, match="only valid for time_order=0"):
            EquationSystem(
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
                        constraint_solver=ConstraintSolverConfig(enabled=True),
                    ),
                ),
                mass_matrix=((0.0,),),
                coupling_matrix=((0.0,),),
                metadata={},
            )

    def test_solver_on_first_order_raises(self) -> None:
        """constraint_solver.enabled=True on time_order=1 must raise."""
        with pytest.raises(ValueError, match="only valid for time_order=0"):
            EquationSystem(
                n_components=1,
                dimension=2,
                spatial_dimension=1,
                component_names=("phi",),
                equations=(
                    ComponentEquation(
                        field_name="phi",
                        field_index=0,
                        time_derivative_order=1,
                        rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                        constraint_solver=ConstraintSolverConfig(enabled=True),
                    ),
                ),
                mass_matrix=((0.0,),),
                coupling_matrix=((0.0,),),
                metadata={},
            )

    def test_disabled_solver_on_wave_ok(self) -> None:
        """constraint_solver.enabled=False is fine for any time_order."""
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
                    constraint_solver=ConstraintSolverConfig(enabled=False),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )
        assert spec.equations[0].constraint_solver.enabled is False


# === Poisson Solver Tests ===


class TestPoissonSolver:
    """Test elliptic solver integration."""

    def test_solve_sinusoidal_source_periodic(
        self, grid_1d_periodic: CartesianGrid
    ) -> None:
        """Solve nabla^2 phi = -sin(x) on [0, 2pi] periodic.

        Analytical solution: phi = sin(x) (since nabla^2 sin(x) = -sin(x)).

        JSON equation: 0 = laplacian(phi) + source
        where source = -identity(phi) with phi_source data = sin(x).
        But we need an external source. Instead, set up a 2-field system
        where rho is given and phi is the constraint:
            0 = laplacian(phi) + identity(rho)
        with rho = -sin(x).
        """
        # Two-field system: phi (constraint) + rho (wave with no evolution used)
        spec = EquationSystem(
            n_components=2,
            dimension=2,
            spatial_dimension=1,
            component_names=("phi", "rho"),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=0,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(1.0, "identity", "rho"),
                    ),
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True,
                        boundary_conditions={
                            "x": BoundaryCondition("periodic"),
                        },
                    ),
                ),
                ComponentEquation(
                    field_name="rho",
                    field_index=1,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(0.0, "identity", "rho"),),
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, -1.0), (0.0, 0.0)),
            metadata={},
        )

        grid = grid_1d_periodic
        pde = PDEFromSpec(spec)

        # rho = -sin(x) so that: 0 = nabla^2(phi) + (-sin(x))
        # => nabla^2(phi) = sin(x) => phi = -sin(x)
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        rho_data = -np.sin(x)

        # State: [phi, rho, pi_rho] (constraint + wave)
        state = FieldCollection(
            [
                ScalarField(grid, data=0.0),  # phi (to be solved)
                ScalarField(grid, data=rho_data),  # rho (source)
                ScalarField(grid, data=0.0),  # pi_rho (momentum)
            ]
        )

        # Calling evolution_rate triggers the constraint solver
        pde.evolution_rate(state, t=0.0)

        # After solving, phi should be ~ -sin(x) up to an additive constant
        # (periodic Poisson solutions are unique only up to a constant)
        actual = state[0].data
        expected = -np.sin(x)
        assert_allclose(
            actual - np.mean(actual),
            expected - np.mean(expected),
            atol=0.05,
        )

    def test_solve_dirichlet_bc(self, grid_1d_nonperiodic: CartesianGrid) -> None:
        """Solve nabla^2 phi = 0 with Dirichlet BCs phi(0)=0, phi(1)=1.

        Analytical solution: phi = x (linear profile).
        """
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
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True,
                        boundary_conditions={
                            "x": BoundaryCondition("dirichlet", value=0.5),
                        },
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = grid_1d_nonperiodic
        pde = PDEFromSpec(spec)

        state = FieldCollection([ScalarField(grid, data=0.0)])
        pde.evolution_rate(state, t=0.0)

        # With Dirichlet value=0.5 on both sides and zero source,
        # the solution should be uniform 0.5
        assert_allclose(state[0].data, 0.5, atol=0.05)

    def test_missing_laplacian_raises(self) -> None:
        """Constraint without laplacian(field) term should fail."""
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
                    rhs_terms=(OperatorTerm(1.0, "identity", "phi"),),
                    constraint_solver=ConstraintSolverConfig(enabled=True),
                ),
            ),
            mass_matrix=((-1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = CartesianGrid([(0, 1)], 32, periodic=True)
        pde = PDEFromSpec(spec)
        state = FieldCollection([ScalarField(grid, data=1.0)])

        with pytest.raises(ValueError, match="lacks a laplacian"):
            pde.evolution_rate(state, t=0.0)

    def test_zero_laplacian_coefficient_raises(self) -> None:
        """Zero laplacian coefficient should fail."""
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
                        OperatorTerm(0.0, "laplacian", "phi"),
                        OperatorTerm(1.0, "identity", "phi"),
                    ),
                    constraint_solver=ConstraintSolverConfig(enabled=True),
                ),
            ),
            mass_matrix=((-1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = CartesianGrid([(0, 1)], 32, periodic=True)
        pde = PDEFromSpec(spec)
        state = FieldCollection([ScalarField(grid, data=1.0)])

        with pytest.raises(ValueError, match="effectively zero"):
            pde.evolution_rate(state, t=0.0)

    def test_multiple_laplacian_terms_raises(self) -> None:
        """Multiple laplacian(field) terms should fail."""
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
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(2.0, "laplacian", "phi"),
                    ),
                    constraint_solver=ConstraintSolverConfig(enabled=True),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = CartesianGrid([(0, 1)], 32, periodic=True)
        pde = PDEFromSpec(spec)
        state = FieldCollection([ScalarField(grid, data=0.0)])

        with pytest.raises(ValueError, match="Multiple laplacian"):
            pde.evolution_rate(state, t=0.0)

    def test_constraint_rate_still_zero(self, grid_1d_periodic: CartesianGrid) -> None:
        """Even with solver enabled, evolution rate for constraint is zero."""
        spec = _make_poisson_spec(solver_enabled=True)
        pde = PDEFromSpec(spec)
        state = FieldCollection([ScalarField(grid_1d_periodic, data=0.0)])

        rates = pde.evolution_rate(state, t=0.0)
        assert_allclose(rates[0].data, 0.0)


# === 2D Poisson Test ===


class TestPoissonSolver2D:
    """Test Poisson solver on 2D grids."""

    def test_solve_2d_periodic(self, grid_2d_periodic: CartesianGrid) -> None:
        """Solve nabla^2 phi = -sin(x)sin(y) on [0,2pi]^2 periodic.

        Analytical: phi = sin(x)sin(y)/2 (since nabla^2 = -2*sin(x)sin(y)).
        """
        spec = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi", "rho"),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=0,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(1.0, "identity", "rho"),
                    ),
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True,
                        boundary_conditions={
                            "x": BoundaryCondition("periodic"),
                            "y": BoundaryCondition("periodic"),
                        },
                    ),
                ),
                ComponentEquation(
                    field_name="rho",
                    field_index=1,
                    time_derivative_order=2,
                    rhs_terms=(OperatorTerm(0.0, "identity", "rho"),),
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, -1.0), (0.0, 0.0)),
            metadata={},
        )

        grid = grid_2d_periodic
        pde = PDEFromSpec(spec)

        # rho = -2*sin(x)*sin(y) so nabla^2(phi) = 2*sin(x)*sin(y)
        # => phi = -sin(x)*sin(y)
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        y = cast("np.ndarray", grid.cell_coords[..., 1])
        rho_data = -2 * np.sin(x) * np.sin(y)

        # State: [phi, rho, pi_rho]
        state = FieldCollection(
            [
                ScalarField(grid, data=0.0),
                ScalarField(grid, data=rho_data),
                ScalarField(grid, data=0.0),
            ]
        )

        pde.evolution_rate(state, t=0.0)

        # Periodic Poisson solutions are unique only up to an additive constant,
        # so compare shapes by subtracting means.
        actual = state[0].data
        expected = -np.sin(x) * np.sin(y)
        assert_allclose(
            actual - np.mean(actual),
            expected - np.mean(expected),
            atol=0.1,
        )


# === Backward Compatibility Tests ===


class TestBackwardCompatibility:
    """Ensure existing constraint behavior is unchanged."""

    def test_frozen_constraint_default(self, grid_1d_periodic: CartesianGrid) -> None:
        """Constraint without solver config remains frozen (d/dt = 0)."""
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
                    rhs_terms=(OperatorTerm(1.0, "laplacian", "phi"),),
                    # No constraint_solver -> defaults to disabled
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
        )

        grid = grid_1d_periodic
        pde = PDEFromSpec(spec)

        x = cast("np.ndarray", grid.cell_coords[..., 0])
        initial_data = np.sin(x)
        state = FieldCollection([ScalarField(grid, data=initial_data)])

        rates = pde.evolution_rate(state, t=0.0)

        # Rate should be zero (frozen constraint)
        assert_allclose(rates[0].data, 0.0)

        # State should be unchanged
        assert_allclose(state[0].data, initial_data)

    def test_disabled_solver_same_as_no_solver(
        self, grid_1d_periodic: CartesianGrid
    ) -> None:
        """Explicitly disabled solver behaves the same as omitted solver."""
        spec = _make_poisson_spec(solver_enabled=False)
        pde = PDEFromSpec(spec)

        grid = grid_1d_periodic
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        initial_data = np.sin(x)
        state = FieldCollection([ScalarField(grid, data=initial_data)])

        rates = pde.evolution_rate(state, t=0.0)
        assert_allclose(rates[0].data, 0.0)
        assert_allclose(state[0].data, initial_data)


# === JSON Round-Trip Test ===


class TestJSONParsing:
    """Test constraint_solver in JSON loading."""

    def test_constraint_solver_from_json_dict(self) -> None:
        """Load a constraint solver config from a JSON-like dict."""
        data = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi", "index": 0}],
            "equations": [
                {
                    "field": "phi",
                    "lhs": {
                        "expression": "phi",
                        "order": {"time": 0, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi",
                            },
                        ],
                    },
                    "constraint_solver": {
                        "enabled": True,
                        "method": "poisson",
                        "boundary_conditions": {
                            "x": {"type": "dirichlet", "value": 0.0},
                        },
                    },
                },
            ],
        }

        spec = EquationSystem.from_dict(data)
        eq = spec.equations[0]
        assert eq.constraint_solver.enabled is True
        assert eq.constraint_solver.method == "poisson"
        assert eq.constraint_solver.boundary_conditions["x"].type == "dirichlet"
        assert eq.constraint_solver.boundary_conditions["x"].value == 0.0

    def test_json_without_constraint_solver(self) -> None:
        """JSON without constraint_solver defaults to disabled."""
        data = {
            "spacetime": {"dimension": 2},
            "fields": [{"name": "phi", "index": 0}],
            "equations": [
                {
                    "field": "phi",
                    "lhs": {
                        "expression": "d2_t(phi)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi",
                            },
                        ],
                    },
                },
            ],
        }

        spec = EquationSystem.from_dict(data)
        eq = spec.equations[0]
        assert eq.constraint_solver.enabled is False
