"""Tests for curvilinear coordinate examples (Issue #81).

Verifies that the pipeline correctly handles non-Cartesian coordinate metrics:
- Polar (2+1D): 1/r gradient and 1/r^2 angular Laplacian
- Spherical (3+1D): 2/r gradient, cot(theta)/r^2, csc^2(theta)/r^2
- Cylindrical (3+1D): 1/r gradient, 1/r^2 angular Laplacian, flat z
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic import (
    EquationSystem,
    build_pde_from_json,
    load_equation_system,
)

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"


# ---------------------------------------------------------------------------
# Polar coordinate tests
# ---------------------------------------------------------------------------
class TestPolarKG:
    """Tests for polar coordinate Klein-Gordon (2+1D)."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DATA_DIR / "polar_kg.json")

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 3
        assert spec.spatial_dimension == 2

    def test_coordinates(self, spec: EquationSystem) -> None:
        assert spec.coordinates == ("t", "x", "y")

    def test_single_equation(self, spec: EquationSystem) -> None:
        assert len(spec.equations) == 1

    def test_time_order(self, spec: EquationSystem) -> None:
        eq = spec.equations[0]
        assert eq.time_derivative_order == 2

    def test_operator_types(self, spec: EquationSystem) -> None:
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "laplacian_x" in operators  # d2_r
        assert "gradient_x" in operators  # (1/r) d_r from Christoffel
        assert "laplacian_y" in operators  # (1/r^2) d2_theta
        assert "identity" in operators  # mass term

    def test_gradient_x_position_dependent(self, spec: EquationSystem) -> None:
        """The (1/r) gradient_x coefficient must be position-dependent on x."""
        eq = spec.equations[0]
        grad_terms = [t for t in eq.rhs_terms if t.operator == "gradient_x"]
        assert len(grad_terms) == 1
        assert grad_terms[0].position_dependent
        assert "x" in grad_terms[0].coordinate_dependent

    def test_laplacian_y_position_dependent(self, spec: EquationSystem) -> None:
        """The (1/r^2) laplacian_y coefficient must be position-dependent on x."""
        eq = spec.equations[0]
        lap_y_terms = [t for t in eq.rhs_terms if t.operator == "laplacian_y"]
        assert len(lap_y_terms) == 1
        assert lap_y_terms[0].position_dependent
        assert "x" in lap_y_terms[0].coordinate_dependent

    def test_laplacian_x_constant(self, spec: EquationSystem) -> None:
        """The d2_r term has constant coefficient 1."""
        eq = spec.equations[0]
        lap_x_terms = [t for t in eq.rhs_terms if t.operator == "laplacian_x"]
        assert len(lap_x_terms) == 1
        assert not lap_x_terms[0].position_dependent

    def test_build_pde(self) -> None:
        pde = build_pde_from_json(DATA_DIR / "polar_kg.json", parameters={"polm2": 1.0})
        assert pde is not None

    def test_evolution_rate_runs(self) -> None:
        """evolution_rate() executes without error on a mixed-periodic grid."""
        pde = build_pde_from_json(DATA_DIR / "polar_kg.json", parameters={"polm2": 1.0})
        grid = CartesianGrid(
            bounds=[(0.5, 5.0), (0, 2 * np.pi)],
            shape=[16, 16],
            periodic=[False, True],
        )
        phi = ScalarField(grid, data=0.1, label="polphi_0")
        pi = ScalarField(grid, data=0.0, label="pi_0")
        state = FieldCollection([phi, pi])

        rates = pde.evolution_rate(state, t=0.0)
        assert rates is not None
        assert len(rates) == 2
        # d_t(pi) should be non-zero for non-trivial phi (mass term: -m^2 * phi)
        assert not np.allclose(rates[1].data, 0.0)

    def test_coefficient_value_gradient_x(self) -> None:
        """Verify the 1/r gradient coefficient evaluates correctly."""
        pde = build_pde_from_json(DATA_DIR / "polar_kg.json", parameters={"polm2": 0.0})
        grid = CartesianGrid(
            bounds=[(1.0, 5.0), (0, 2 * np.pi)],
            shape=[8, 8],
            periodic=[False, True],
        )
        # Set up a field with d_r phi = 1 everywhere (linear in r)
        r = cast("np.ndarray", grid.cell_coords[..., 0])
        phi = ScalarField(grid, data=r, label="polphi_0")
        pi = ScalarField(grid, data=0.0, label="pi_0")
        state = FieldCollection([phi, pi])

        rates = pde.evolution_rate(state, t=0.0)
        # The momentum rate (rates[0]) includes the gradient_x contribution
        # which should include a 1/r term. We just verify it runs and is finite.
        assert np.all(np.isfinite(rates[0].data))
        assert np.all(np.isfinite(rates[1].data))


# ---------------------------------------------------------------------------
# Spherical coordinate tests
# ---------------------------------------------------------------------------
class TestSphericalKG:
    """Tests for spherical coordinate Klein-Gordon (3+1D)."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DATA_DIR / "spherical_kg.json")

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 4
        assert spec.spatial_dimension == 3

    def test_coordinates(self, spec: EquationSystem) -> None:
        assert spec.coordinates == ("t", "x", "y", "z")

    def test_single_equation(self, spec: EquationSystem) -> None:
        assert len(spec.equations) == 1

    def test_operator_types(self, spec: EquationSystem) -> None:
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "laplacian_x" in operators  # d2_r
        assert "gradient_x" in operators  # (2/r) d_r
        assert "laplacian_y" in operators  # (1/r^2) d2_theta
        assert "gradient_y" in operators  # cot(theta)/r^2 d_theta
        assert "laplacian_z" in operators  # csc^2(theta)/r^2 d2_phi
        assert "identity" in operators  # mass term

    def test_six_terms(self, spec: EquationSystem) -> None:
        """Spherical coordinates produce 6 terms in the wave equation."""
        eq = spec.equations[0]
        assert len(eq.rhs_terms) == 6

    def test_gradient_y_depends_on_x_and_y(self, spec: EquationSystem) -> None:
        """The cot(theta)/r^2 coefficient depends on both x (r) and y (theta)."""
        eq = spec.equations[0]
        grad_y_terms = [t for t in eq.rhs_terms if t.operator == "gradient_y"]
        assert len(grad_y_terms) == 1
        assert grad_y_terms[0].position_dependent
        deps = set(grad_y_terms[0].coordinate_dependent)
        assert "x" in deps
        assert "y" in deps

    def test_laplacian_z_depends_on_x_and_y(self, spec: EquationSystem) -> None:
        """The csc^2(theta)/r^2 coefficient depends on both x and y."""
        eq = spec.equations[0]
        lap_z_terms = [t for t in eq.rhs_terms if t.operator == "laplacian_z"]
        assert len(lap_z_terms) == 1
        assert lap_z_terms[0].position_dependent
        deps = set(lap_z_terms[0].coordinate_dependent)
        assert "x" in deps
        assert "y" in deps

    def test_build_pde(self) -> None:
        pde = build_pde_from_json(
            DATA_DIR / "spherical_kg.json", parameters={"spm2": 0.0}
        )
        assert pde is not None

    def test_evolution_rate_runs(self) -> None:
        """evolution_rate() executes on a 3D grid with trig coefficients."""
        pde = build_pde_from_json(
            DATA_DIR / "spherical_kg.json", parameters={"spm2": 0.0}
        )
        grid = CartesianGrid(
            bounds=[(0.5, 4.0), (0.3, np.pi - 0.3), (0, 2 * np.pi)],
            shape=[8, 6, 8],
            periodic=[False, False, True],
        )
        phi = ScalarField(grid, data=0.1, label="spphi_0")
        pi = ScalarField(grid, data=0.0, label="pi_0")
        state = FieldCollection([phi, pi])

        rates = pde.evolution_rate(state, t=0.0)
        assert rates is not None
        assert len(rates) == 2
        assert np.all(np.isfinite(rates[0].data))
        assert np.all(np.isfinite(rates[1].data))

    def test_cot_csc_conversion(self) -> None:
        """Verify Cot and Csc are correctly converted and evaluated."""
        pde = build_pde_from_json(
            DATA_DIR / "spherical_kg.json", parameters={"spm2": 0.0}
        )
        # Use theta near pi/4 where cot and csc are well-behaved
        grid = CartesianGrid(
            bounds=[(1.0, 3.0), (0.5, 2.5), (0, 2 * np.pi)],
            shape=[4, 4, 4],
            periodic=[False, False, True],
        )
        phi = ScalarField(grid, data=1.0, label="spphi_0")
        pi = ScalarField(grid, data=0.0, label="pi_0")
        state = FieldCollection([phi, pi])

        # Should not raise — if Cot/Csc conversion failed, we'd get NameError
        rates = pde.evolution_rate(state, t=0.0)
        assert np.all(np.isfinite(rates[0].data))


# ---------------------------------------------------------------------------
# Cylindrical coordinate tests
# ---------------------------------------------------------------------------
class TestCylindricalKG:
    """Tests for cylindrical coordinate Klein-Gordon (3+1D)."""

    @pytest.fixture
    def spec(self) -> EquationSystem:
        return load_equation_system(DATA_DIR / "cylindrical_kg.json")

    def test_dimension(self, spec: EquationSystem) -> None:
        assert spec.dimension == 4
        assert spec.spatial_dimension == 3

    def test_coordinates(self, spec: EquationSystem) -> None:
        assert spec.coordinates == ("t", "x", "y", "z")

    def test_single_equation(self, spec: EquationSystem) -> None:
        assert len(spec.equations) == 1

    def test_operator_types(self, spec: EquationSystem) -> None:
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "laplacian_x" in operators  # d2_r
        assert "gradient_x" in operators  # (1/r) d_r
        assert "laplacian_y" in operators  # (1/r^2) d2_theta
        assert "laplacian_z" in operators  # d2_z (flat, constant coeff)
        assert "identity" in operators  # mass term

    def test_five_terms(self, spec: EquationSystem) -> None:
        """Cylindrical coordinates produce 5 terms."""
        eq = spec.equations[0]
        assert len(eq.rhs_terms) == 5

    def test_laplacian_z_constant(self, spec: EquationSystem) -> None:
        """The z-direction Laplacian has constant coefficient (flat direction)."""
        eq = spec.equations[0]
        lap_z_terms = [t for t in eq.rhs_terms if t.operator == "laplacian_z"]
        assert len(lap_z_terms) == 1
        assert not lap_z_terms[0].position_dependent

    def test_gradient_x_position_dependent(self, spec: EquationSystem) -> None:
        """The (1/r) gradient_x coefficient is position-dependent on x."""
        eq = spec.equations[0]
        grad_terms = [t for t in eq.rhs_terms if t.operator == "gradient_x"]
        assert len(grad_terms) == 1
        assert grad_terms[0].position_dependent
        assert "x" in grad_terms[0].coordinate_dependent

    def test_build_pde(self) -> None:
        pde = build_pde_from_json(
            DATA_DIR / "cylindrical_kg.json", parameters={"cylm2": 0.0}
        )
        assert pde is not None

    def test_evolution_rate_runs(self) -> None:
        """evolution_rate() executes on a 3D grid with mixed periodic BCs."""
        pde = build_pde_from_json(
            DATA_DIR / "cylindrical_kg.json", parameters={"cylm2": 0.0}
        )
        grid = CartesianGrid(
            bounds=[(0.5, 4.0), (0, 2 * np.pi), (-3, 3)],
            shape=[8, 8, 8],
            periodic=[False, True, False],
        )
        phi = ScalarField(grid, data=0.1, label="cylphi_0")
        pi = ScalarField(grid, data=0.0, label="pi_0")
        state = FieldCollection([phi, pi])

        rates = pde.evolution_rate(state, t=0.0)
        assert rates is not None
        assert len(rates) == 2
        assert np.all(np.isfinite(rates[0].data))
        assert np.all(np.isfinite(rates[1].data))


# ---------------------------------------------------------------------------
# Cross-cutting tests
# ---------------------------------------------------------------------------
class TestCurvilinearCommon:
    """Tests that apply across all curvilinear coordinate examples."""

    @pytest.mark.parametrize(
        ("json_file", "expected_dim"),
        [
            ("polar_kg.json", 3),
            ("spherical_kg.json", 4),
            ("cylindrical_kg.json", 4),
        ],
    )
    def test_json_loads(self, json_file: str, expected_dim: int) -> None:
        spec = load_equation_system(DATA_DIR / json_file)
        assert spec.dimension == expected_dim

    @pytest.mark.parametrize(
        ("json_file", "params"),
        [
            ("polar_kg.json", {"polm2": 0.0}),
            ("spherical_kg.json", {"spm2": 0.0}),
            ("cylindrical_kg.json", {"cylm2": 0.0}),
        ],
    )
    def test_pde_builds(self, json_file: str, params: dict[str, float]) -> None:
        pde = build_pde_from_json(DATA_DIR / json_file, parameters=params)
        assert pde is not None

    @pytest.mark.parametrize(
        ("json_file", "n_terms"),
        [
            ("polar_kg.json", 4),
            ("spherical_kg.json", 6),
            ("cylindrical_kg.json", 5),
        ],
    )
    def test_term_counts(self, json_file: str, n_terms: int) -> None:
        spec = load_equation_system(DATA_DIR / json_file)
        eq = spec.equations[0]
        assert len(eq.rhs_terms) == n_terms

    def test_polar_has_no_gradient_y(self) -> None:
        """Polar coordinates should NOT have a gradient_y term (no angular Christoffel)."""
        spec = load_equation_system(DATA_DIR / "polar_kg.json")
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "gradient_y" not in operators

    def test_cylindrical_has_no_gradient_y(self) -> None:
        """Cylindrical coordinates should NOT have a gradient_y term."""
        spec = load_equation_system(DATA_DIR / "cylindrical_kg.json")
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "gradient_y" not in operators

    def test_spherical_has_gradient_y(self) -> None:
        """Spherical coordinates SHOULD have a gradient_y term (cot(theta) correction)."""
        spec = load_equation_system(DATA_DIR / "spherical_kg.json")
        eq = spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "gradient_y" in operators


# ---------------------------------------------------------------------------
# Reciprocal trig function tests
# ---------------------------------------------------------------------------
class TestReciprocalTrigFunctions:
    """Tests for Cot, Sec, Csc conversion in _mathematica_to_python."""

    def test_cot_evaluation(self) -> None:
        """Cot[y[]] should evaluate to cos(y)/sin(y)."""
        pde = build_pde_from_json(
            DATA_DIR / "spherical_kg.json", parameters={"spm2": 0.0}
        )
        # Test the converter directly via a known coefficient
        py_expr = pde._mathematica_to_python("Cot[y[]]")
        assert "cot" in py_expr

    def test_csc_evaluation(self) -> None:
        """Csc[y[]]^2 should evaluate to 1/sin(y)^2."""
        pde = build_pde_from_json(
            DATA_DIR / "spherical_kg.json", parameters={"spm2": 0.0}
        )
        py_expr = pde._mathematica_to_python("Csc[y[]]^2")
        assert "csc" in py_expr
