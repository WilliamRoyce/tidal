"""Tests for 3+1D validation (Phase 11).

Verifies that the pipeline handles 4D spacetime (t,x,y,z) correctly:
- JSON loads with dimension=4 and spatial_dimension=3
- State creation works on 3D grids
- Evolution rate computes correctly using laplacian_x/y/z
- Short simulation runs without errors
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system
from torsion_gertsenshtein.symbolic.pde_builder import (
    PDEFromSpec,
    create_initial_state,
)
from torsion_gertsenshtein.utils import normalize_solve_result


@pytest.fixture
def kg_3d_json_path() -> Path:
    """Path to the 3+1D Klein-Gordon JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "klein_gordon_3d.json"


@pytest.fixture
def grid_3d() -> CartesianGrid:
    """3D periodic grid for 3+1D testing."""
    return CartesianGrid(
        [(0, 10), (0, 10), (0, 10)],
        [16, 16, 16],
        periodic=True,
    )


@pytest.fixture
def grid_3d_small() -> CartesianGrid:
    """Small 3D periodic grid for fast testing."""
    return CartesianGrid(
        [(0, 10), (0, 10), (0, 10)],
        [8, 8, 8],
        periodic=True,
    )


class TestLoad3DJSON:
    """Tests for loading 3+1D equation specification."""

    def test_load_kg_3d(self, kg_3d_json_path: Path) -> None:
        """Test loading 3+1D Klein-Gordon JSON."""
        spec = load_equation_system(kg_3d_json_path)

        assert spec.dimension == 4  # noqa: PLR2004
        assert spec.spatial_dimension == 3  # noqa: PLR2004
        assert spec.n_components == 1
        assert spec.component_names == ("phi",)

    def test_kg_3d_equations(self, kg_3d_json_path: Path) -> None:
        """Test that equations have expected operators."""
        spec = load_equation_system(kg_3d_json_path)
        eq = spec.equations[0]

        assert eq.field_name == "phi"
        assert eq.time_derivative_order == 2  # noqa: PLR2004

        operators = {term.operator for term in eq.rhs_terms}
        assert "laplacian_x" in operators
        assert "laplacian_y" in operators
        assert "laplacian_z" in operators
        assert "identity" in operators

    def test_build_pde_3d(self, kg_3d_json_path: Path) -> None:
        """Test building PDE from 3+1D JSON."""
        pde = build_pde_from_json(kg_3d_json_path)

        assert isinstance(pde, PDEFromSpec)
        assert pde.n_components == 1
        assert pde.spec.spatial_dimension == 3  # noqa: PLR2004


class TestEvolution3D:
    """Tests for 3+1D evolution rate computation."""

    def test_evolution_rate_shape(
        self, kg_3d_json_path: Path, grid_3d_small: CartesianGrid
    ) -> None:
        """Evolution rate has correct number of fields on 3D grid."""
        pde = build_pde_from_json(kg_3d_json_path)
        state = create_initial_state(grid_3d_small, pde.spec)

        rates = pde.evolution_rate(state)

        assert len(rates) == 2  # noqa: PLR2004 — phi + momentum

    def test_uniform_field_mass_term(
        self, kg_3d_json_path: Path, grid_3d_small: CartesianGrid
    ) -> None:
        """Uniform field: laplacians = 0, so d_t pi = -m^2 phi."""
        pde = build_pde_from_json(kg_3d_json_path, parameters={"m2": 2.0})

        state = FieldCollection([
            ScalarField(grid_3d_small, data=1.0),  # phi = 1
            ScalarField(grid_3d_small, data=0.0),  # pi = 0
        ])

        rates = pde.evolution_rate(state)

        # d/dt phi = pi = 0
        assert_allclose(rates[0].data, 0.0, atol=1e-10)

        # d/dt pi = laplacian_x(1) + laplacian_y(1) + laplacian_z(1) - 2*1
        #         = 0 + 0 + 0 - 2 = -2
        assert_allclose(rates[1].data, -2.0, atol=1e-10)

    def test_sinusoidal_field_3d(
        self, kg_3d_json_path: Path, grid_3d: CartesianGrid
    ) -> None:
        """3D sinusoidal field: verify directional laplacians sum correctly."""
        pde = build_pde_from_json(kg_3d_json_path, parameters={"m2": 0.0})

        x = cast("np.ndarray", grid_3d.cell_coords[..., 0])
        y = cast("np.ndarray", grid_3d.cell_coords[..., 1])
        z = cast("np.ndarray", grid_3d.cell_coords[..., 2])
        kx = 2 * np.pi / 10.0
        ky = 2 * np.pi / 10.0
        kz = 2 * np.pi / 10.0

        phi_data = np.sin(kx * x) * np.sin(ky * y) * np.sin(kz * z)
        state = FieldCollection([
            ScalarField(grid_3d, data=phi_data),
            ScalarField(grid_3d, data=0.0),
        ])

        rates = pde.evolution_rate(state)

        # d/dt pi = laplacian_x + laplacian_y + laplacian_z (no mass)
        # = (-kx^2 - ky^2 - kz^2) * sin(kx*x)*sin(ky*y)*sin(kz*z)
        expected = -(kx**2 + ky**2 + kz**2) * phi_data
        # Larger tolerance for 3D operators on coarse 16^3 grid
        assert_allclose(rates[1].data, expected, rtol=0.06)

    def test_grid_dimension_mismatch_raises(
        self, kg_3d_json_path: Path
    ) -> None:
        """Passing 2D grid to 3+1D PDE raises ValueError."""
        pde = build_pde_from_json(kg_3d_json_path)
        grid_2d = CartesianGrid([(0, 10), (0, 10)], [8, 8], periodic=True)

        state = FieldCollection([
            ScalarField(grid_2d, data=0.0),
            ScalarField(grid_2d, data=0.0),
        ])

        with pytest.raises(ValueError, match="Grid dimension"):
            pde.evolution_rate(state)


class TestSimulation3D:
    """Integration test: short 3+1D simulation."""

    def test_gaussian_pulse_spreads_3d(
        self, kg_3d_json_path: Path, grid_3d_small: CartesianGrid
    ) -> None:
        """Gaussian pulse in 3D should spread over time."""
        pde = build_pde_from_json(kg_3d_json_path, parameters={"m2": 0.0})

        x = cast("np.ndarray", grid_3d_small.cell_coords[..., 0])
        y = cast("np.ndarray", grid_3d_small.cell_coords[..., 1])
        z = cast("np.ndarray", grid_3d_small.cell_coords[..., 2])

        center = 5.0
        width = 2.0
        phi_data = np.exp(
            -((x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2)
            / (2 * width**2)
        )
        initial_peak = float(np.max(phi_data))

        state = create_initial_state(
            grid_3d_small,
            pde.spec,
            field_data={"phi": phi_data},
        )

        # Run very short simulation
        result = pde.solve(state, t_range=0.5, dt=0.01)
        sol = normalize_solve_result(result)

        # Peak should decrease (pulse splits and propagates)
        final_peak = float(np.max(sol[0].data))
        assert final_peak < initial_peak
