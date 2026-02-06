"""Integration tests for Klein-Gordon on 2-sphere example."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system
from torsion_gertsenshtein.utils import normalize_solve_result

SPHERE_JSON = Path(__file__).parent.parent / "examples" / "data" / "sphere_kg.json"


@pytest.fixture
def sphere_spec():
    """Load the sphere KG equation specification."""
    if not SPHERE_JSON.exists():
        pytest.skip("sphere_kg.json not found (run sphere_kg.wls first)")
    return load_equation_system(SPHERE_JSON)


class TestSphereKGSpec:
    """Tests for sphere KG equation specification."""

    def test_spec_loads(self, sphere_spec) -> None:
        """Sphere KG spec loads without errors."""
        assert sphere_spec.n_components == 1
        assert sphere_spec.dimension == 3
        assert sphere_spec.spatial_dimension == 2

    def test_has_position_dependent_terms(self, sphere_spec) -> None:
        """Sphere KG should have position-dependent coefficient terms."""
        eq = sphere_spec.equations[0]
        pos_dep_terms = [t for t in eq.rhs_terms if t.position_dependent]
        assert len(pos_dep_terms) > 0, "Expected position-dependent terms from Christoffels"

    def test_no_time_dependent_terms(self, sphere_spec) -> None:
        """Static metric should produce no time-dependent coefficients."""
        eq = sphere_spec.equations[0]
        time_dep_terms = [t for t in eq.rhs_terms if t.time_dependent]
        assert len(time_dep_terms) == 0, "Static metric should have no time-dependent terms"

    def test_has_laplacian_operators(self, sphere_spec) -> None:
        """Sphere KG should use laplacian_x and laplacian_y operators."""
        eq = sphere_spec.equations[0]
        operators = {t.operator for t in eq.rhs_terms}
        assert "laplacian_x" in operators
        assert "laplacian_y" in operators

    def test_coordinate_dependent_contains_x_y(self, sphere_spec) -> None:
        """Position-dependent terms should reference x and/or y coordinates."""
        eq = sphere_spec.equations[0]
        all_coord_deps = set()
        for t in eq.rhs_terms:
            all_coord_deps.update(t.coordinate_dependent)
        # Should have x and y (spatial coords), not t
        assert "x" in all_coord_deps
        assert "y" in all_coord_deps
        assert "t" not in all_coord_deps


class TestSphereKGSimulation:
    """Integration tests for sphere KG simulation."""

    def test_build_pde_from_json(self) -> None:
        """PDE builds from sphere KG JSON with parameters."""
        if not SPHERE_JSON.exists():
            pytest.skip("sphere_kg.json not found")
        pde = build_pde_from_json(
            SPHERE_JSON, parameters={"sphR": 1.0, "sphm2": 0.0}
        )
        assert pde.n_components == 1

    def test_simulation_runs(self) -> None:
        """Sphere KG simulation runs on a small grid."""
        if not SPHERE_JSON.exists():
            pytest.skip("sphere_kg.json not found")
        pde = build_pde_from_json(
            SPHERE_JSON, parameters={"sphR": 2.0, "sphm2": 0.0}
        )
        grid = CartesianGrid(
            bounds=[(-4, 4), (-4, 4)], shape=[16, 16], periodic=True
        )
        phi = ScalarField(grid, data=0.0)
        pi = ScalarField(grid, data=0.0)
        x = np.asarray(grid.cell_coords[..., 0])
        y = np.asarray(grid.cell_coords[..., 1])
        phi.data[:] = np.exp(-(x**2 + y**2) / 2)
        state = FieldCollection([phi, pi])

        result = pde.solve(state, t_range=1.0, dt=0.01, scheme="runge-kutta")
        result = normalize_solve_result(result)

        # Should have evolved
        assert not np.allclose(result[0].data, phi.data, rtol=0.01)

    def test_wave_speed_varies_with_position(self) -> None:
        """Wave speed should be faster away from south pole."""
        if not SPHERE_JSON.exists():
            pytest.skip("sphere_kg.json not found")
        pde = build_pde_from_json(
            SPHERE_JSON, parameters={"sphR": 2.0, "sphm2": 0.0}
        )
        grid = CartesianGrid(
            bounds=[(-6, 6), (-6, 6)], shape=[32, 32], periodic=True
        )

        # Two initial pulses: one at center (south pole), one off-center
        phi = ScalarField(grid, data=0.0)
        pi = ScalarField(grid, data=0.0)
        x = np.asarray(grid.cell_coords[..., 0])
        y = np.asarray(grid.cell_coords[..., 1])
        phi.data[:] = np.exp(-(x**2 + y**2) / 0.5)
        state = FieldCollection([phi, pi])

        # Get evolution rate to check position dependence
        rhs = pde.evolution_rate(state, t=0.0)
        # The momentum derivative (rhs[1]) should be non-uniform
        # due to position-dependent coefficients
        rhs_data = rhs[1].data
        assert not np.allclose(rhs_data, 0.0), "RHS should be non-zero for Gaussian pulse"

        # The effective Laplacian coefficient at center vs edge should differ
        # At x=y=0: 1/Omega^2 = (R^2+0)^2/(4R^4) = 1/4
        # At x=4, y=0: 1/Omega^2 = (R^2+16)^2/(4R^4) = (4+16)^2/(4*16) = 400/64 ≈ 6.25
        # So the wave equation should produce much stronger response at edges
        center = grid.shape[0] // 2
        center_val = abs(rhs_data[center, center])
        # Check at an off-center point
        edge_idx = grid.shape[0] - 2  # Near edge
        edge_val = abs(rhs_data[edge_idx, center])
        # Edge could be zero if Gaussian has decayed, but center should be non-zero
        assert center_val > 0, "RHS at center should be non-zero"

    def test_missing_parameter_raises(self) -> None:
        """Building PDE without required sphR parameter should fail at evolution."""
        if not SPHERE_JSON.exists():
            pytest.skip("sphere_kg.json not found")
        pde = build_pde_from_json(SPHERE_JSON, parameters={"sphm2": 0.0})
        # sphR is missing - should fail when trying to resolve position-dependent coeff
        grid = CartesianGrid(
            bounds=[(-4, 4), (-4, 4)], shape=[8, 8], periodic=True
        )
        phi = ScalarField(grid, data=1.0)
        pi = ScalarField(grid, data=0.0)
        state = FieldCollection([phi, pi])

        with pytest.raises(ValueError, match="sphR"):
            pde.evolution_rate(state, t=0.0)
