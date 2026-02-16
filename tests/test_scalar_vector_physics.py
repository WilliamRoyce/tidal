"""Physics tests for the scalar-vector coupling stress-test example.

This example stress-tests pipeline handling of:
- Mixed-rank cross-field coupling (scalar phi + vector A)
- Cross-field first_derivative_t and gradient operators
- Epsilon tensor + cross-field coupling in same system
- 4 symbolic constants simultaneously
- 4x4 mass/coupling matrices
- Constraint equation (A_0) with cross-field terms
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest
from pde import CartesianGrid

from tidal.symbolic.json_loader import (
    EquationSystem,
    load_equation_system,
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    create_initial_state,
)

if TYPE_CHECKING:
    from pde.grids.base import GridBase

_JSON_PATH = Path("examples/data/scalar_vector_coupling.json")
_DEFAULT_PARAMS = {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2}

pytestmark = pytest.mark.skipif(
    not _JSON_PATH.exists(),
    reason="scalar_vector_coupling.json not found (run tidal derive)",
)


@pytest.fixture
def spec() -> EquationSystem:
    return load_equation_system(_JSON_PATH)


@pytest.fixture
def pde(spec: EquationSystem) -> PDEFromSpec:
    return PDEFromSpec(spec, parameters=_DEFAULT_PARAMS)


@pytest.fixture
def grid() -> CartesianGrid:
    return CartesianGrid([(0, 10), (0, 10)], 16, periodic=True)


class TestJSONStructure:
    """Tests that the JSON spec loads correctly and has the right structure."""

    def test_json_loads_successfully(self, spec: EquationSystem) -> None:
        assert spec is not None
        assert spec.n_components == 4

    def test_field_names(self, spec: EquationSystem) -> None:
        names = list(spec.component_names)
        assert "phi_0" in names
        assert "A_0" in names
        assert "A_1" in names
        assert "A_2" in names

    def test_phi_is_dynamical(self, spec: EquationSystem) -> None:
        phi_eq = next(eq for eq in spec.equations if eq.field_name == "phi_0")
        assert phi_eq.time_derivative_order == 2

    def test_a0_is_constraint(self, spec: EquationSystem) -> None:
        a0_eq = next(eq for eq in spec.equations if eq.field_name == "A_0")
        assert a0_eq.time_derivative_order == 0

    def test_a1_a2_are_dynamical(self, spec: EquationSystem) -> None:
        for name in ("A_1", "A_2"):
            eq = next(eq for eq in spec.equations if eq.field_name == name)
            assert eq.time_derivative_order == 2, f"{name} should be dynamical"

    def test_spatial_dimension(self, spec: EquationSystem) -> None:
        assert spec.spatial_dimension == 2
        assert spec.dimension == 3


class TestCrossFieldOperators:
    """Tests that cross-field coupling operators are present and correct."""

    def test_phi_has_cross_field_first_derivative_t(self, spec: EquationSystem) -> None:
        """Phi equation has first_derivative_t(A_0) from div coupling."""
        phi_eq = next(eq for eq in spec.equations if eq.field_name == "phi_0")
        fdt_terms = [
            t for t in phi_eq.rhs_terms if t.operator == "first_derivative_t"
        ]
        assert len(fdt_terms) == 1
        assert fdt_terms[0].field == "A_0"

    def test_phi_has_cross_field_gradients(self, spec: EquationSystem) -> None:
        """Phi equation has gradient_x(A_1) and gradient_y(A_2)."""
        phi_eq = next(eq for eq in spec.equations if eq.field_name == "phi_0")
        grad_terms = [
            t
            for t in phi_eq.rhs_terms
            if t.operator.startswith("gradient_") and t.field != "phi_0"
        ]
        fields = {t.field for t in grad_terms}
        assert "A_1" in fields
        assert "A_2" in fields

    def test_a0_has_cross_field_first_derivative_t_phi(
        self, spec: EquationSystem
    ) -> None:
        """A_0 constraint has first_derivative_t(phi_0) from coupling."""
        a0_eq = next(eq for eq in spec.equations if eq.field_name == "A_0")
        fdt_terms = [
            t
            for t in a0_eq.rhs_terms
            if t.operator == "first_derivative_t" and t.field == "phi_0"
        ]
        assert len(fdt_terms) == 1

    def test_a1_has_cross_field_gradient_phi(self, spec: EquationSystem) -> None:
        """A_1 equation has gradient_x(phi_0) from coupling."""
        a1_eq = next(eq for eq in spec.equations if eq.field_name == "A_1")
        grad_phi = [
            t
            for t in a1_eq.rhs_terms
            if t.operator == "gradient_x" and t.field == "phi_0"
        ]
        assert len(grad_phi) == 1

    def test_a2_has_cross_field_gradient_phi(self, spec: EquationSystem) -> None:
        """A_2 equation has gradient_y(phi_0) from coupling."""
        a2_eq = next(eq for eq in spec.equations if eq.field_name == "A_2")
        grad_phi = [
            t
            for t in a2_eq.rhs_terms
            if t.operator == "gradient_y" and t.field == "phi_0"
        ]
        assert len(grad_phi) == 1


class TestSymbolicConstants:
    """Tests that all 4 symbolic constants are properly tracked."""

    def test_four_constants_discovered(self, spec: EquationSystem) -> None:
        from tidal.cli._inspect import discover_parameters

        params = discover_parameters(spec)
        assert "phim2" in params
        assert "Am2" in params
        assert "kCS" in params
        assert "gSV" in params

    def test_symbolic_coefficients_present(self, spec: EquationSystem) -> None:
        """Every equation has at least one symbolic coefficient."""
        for eq in spec.equations:
            symbolic = [
                t.coefficient_symbolic
                for t in eq.rhs_terms
                if t.coefficient_symbolic is not None
            ]
            assert len(symbolic) > 0, (
                f"{eq.field_name} has no symbolic coefficients"
            )


class TestMassMatrix:
    """Tests for the 4x4 mass/coupling matrices."""

    def test_mass_matrix_dimensions(self, spec: EquationSystem) -> None:
        assert len(spec.mass_matrix) == 4
        assert all(len(row) == 4 for row in spec.mass_matrix)

    def test_coupling_matrix_dimensions(self, spec: EquationSystem) -> None:
        assert len(spec.coupling_matrix) == 4
        assert all(len(row) == 4 for row in spec.coupling_matrix)

    def test_coupling_matrix_is_zero(self, spec: EquationSystem) -> None:
        """No identity cross-field terms => coupling matrix is all zero."""
        for row in spec.coupling_matrix:
            for val in row:
                assert val == 0.0

    def test_mass_matrix_symbolic_present(self, spec: EquationSystem) -> None:
        assert spec.mass_matrix_symbolic is not None
        assert len(spec.mass_matrix_symbolic) == 4


class TestPDEEvolution:
    """Tests that the PDE builds and evolves correctly."""

    def test_pde_builds(self, pde: PDEFromSpec) -> None:
        assert pde is not None

    def test_evolution_rate_finite(
        self, pde: PDEFromSpec, spec: EquationSystem, grid: GridBase
    ) -> None:
        state = create_initial_state(grid, spec)
        rates = pde.evolution_rate(state, t=0.0)
        assert np.all(np.isfinite(rates.data))

    def test_gaussian_phi_drives_vector_rates(
        self, pde: PDEFromSpec, spec: EquationSystem, grid: GridBase
    ) -> None:
        """A Gaussian in phi creates nonzero rates in A_1 and A_2
        via the gSV gradient coupling.
        """
        cc: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = grid.cell_coords  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        x: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 0]  # pyright: ignore[reportUnknownVariableType]
        y: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 1]  # pyright: ignore[reportUnknownVariableType]
        gauss = np.exp(-((x - 5) ** 2 + (y - 5) ** 2))  # pyright: ignore[reportUnknownArgumentType]
        state = create_initial_state(grid, spec, field_data={"phi_0": gauss})

        rates = pde.evolution_rate(state, t=0.0)

        # phi_0 rate (pi_phi) should be nonzero (laplacian of Gaussian)
        assert np.max(np.abs(rates.data[1])) > 0.01

        # pi_A_1 and pi_A_2 should be nonzero (from gSV * grad phi)
        assert np.max(np.abs(rates.data[4])) > 0.01
        assert np.max(np.abs(rates.data[6])) > 0.01

    def test_decoupling_limit(
        self, spec: EquationSystem, grid: GridBase
    ) -> None:
        """With gSV=0 and kCS=0, phi decouples from A.

        A Gaussian in phi should NOT drive rates in A components.
        """
        pde = PDEFromSpec(
            spec, parameters={"phim2": 1.0, "Am2": 0.5, "kCS": 0.0, "gSV": 0.0}
        )
        cc: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = grid.cell_coords  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        x: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 0]  # pyright: ignore[reportUnknownVariableType]
        y: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 1]  # pyright: ignore[reportUnknownVariableType]
        gauss = np.exp(-((x - 5) ** 2 + (y - 5) ** 2))  # pyright: ignore[reportUnknownArgumentType]
        state = create_initial_state(grid, spec, field_data={"phi_0": gauss})

        rates = pde.evolution_rate(state, t=0.0)

        # phi rates should still be nonzero (KG dynamics)
        assert np.max(np.abs(rates.data[1])) > 0.01

        # A rates should be zero (decoupled, zero IC)
        assert np.max(np.abs(rates.data[4])) < 1e-10
        assert np.max(np.abs(rates.data[6])) < 1e-10

    def test_short_evolution_stable(
        self, pde: PDEFromSpec, spec: EquationSystem, grid: GridBase
    ) -> None:
        """Short time evolution should remain finite."""
        cc: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = grid.cell_coords  # pyright: ignore[reportUnknownMemberType, reportUnknownVariableType]
        x: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 0]  # pyright: ignore[reportUnknownVariableType]
        y: np.ndarray[tuple[int, ...], np.dtype[np.float64]] = cc[..., 1]  # pyright: ignore[reportUnknownVariableType]
        gauss = 0.1 * np.exp(-((x - 5) ** 2 + (y - 5) ** 2))  # pyright: ignore[reportUnknownArgumentType]
        state = create_initial_state(grid, spec, field_data={"phi_0": gauss})

        raw_result = pde.solve(state, t_range=0.5, dt=0.01)
        assert raw_result is not None
        assert not isinstance(raw_result, tuple)
        result = raw_result
        assert np.all(np.isfinite(result.data))
        assert np.max(np.abs(result.data)) < 1e6
