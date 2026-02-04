"""End-to-end integration tests for the Lagrangian-to-PDE pipeline.

These tests verify the complete pipeline from JSON equation specification
to numerical simulation, using both Klein-Gordon and EM field examples.
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.kgsim.config import KGParameters
from torsion_gertsenshtein.kgsim.equations import KleinGordonPDE
from torsion_gertsenshtein.symbolic import (
    build_pde_from_json,
    load_equation_system,
)
from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec, create_initial_state
from torsion_gertsenshtein.utils import normalize_solve_result
from torsion_gertsenshtein.vectorfield import (
    ComponentFieldParams,
    ComponentGaussianPulse,
    create_gaussian_pulse_state,
)

# === Fixtures ===


@pytest.fixture
def em_json_path() -> Path:
    """Path to the EM 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "em_1d.json"


@pytest.fixture
def kg_json_path() -> Path:
    """Path to the Klein-Gordon 1D JSON file."""
    return Path(__file__).parent.parent / "examples" / "data" / "klein_gordon_1d.json"


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """1D periodic grid for testing."""
    return CartesianGrid([(0, 100)], 128, periodic=True)


@pytest.fixture
def grid_1d_fine() -> CartesianGrid:
    """Fine 1D periodic grid for accurate testing."""
    return CartesianGrid([(0, 100)], 256, periodic=True)


# === Pipeline Integration Tests ===


class TestEMPipeline:
    """Integration tests for the EM field pipeline."""

    def test_load_em_spec(self, em_json_path: Path) -> None:
        """Test that EM specification loads correctly."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        spec = load_equation_system(em_json_path)

        num_em_components = 2
        assert spec.n_components == num_em_components
        assert spec.component_names == ("A_0", "A_1")
        assert spec.metadata.get("gauge") == "lorenz"

    def test_build_em_pde(self, em_json_path: Path) -> None:
        """Test building PDE from EM specification."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        pde = build_pde_from_json(em_json_path)

        assert isinstance(pde, PDEFromSpec)
        num_em_components = 2
        assert pde.n_components == num_em_components

    def test_em_wave_propagation(
        self, em_json_path: Path, grid_1d: CartesianGrid
    ) -> None:
        """Test that EM waves propagate at c=1."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        spec = load_equation_system(em_json_path)
        pde = build_pde_from_json(em_json_path)

        # Create Gaussian pulse in A_1
        pulse = ComponentGaussianPulse(
            center=(50.0,),
            width=5.0,
            amplitude=1.0,
            active_components={"A_1": 1.0},
        )
        initial = pulse.create(grid_1d, spec)

        # Run short simulation
        t_end = 10.0
        result = pde.solve(initial, t_range=t_end, dt=0.01)
        sol = normalize_solve_result(result)

        # Initial pulse should have spread
        # In 1D, a Gaussian splits into two pulses at +/- c*t = +/- 10
        initial_max = np.max(initial.data[2])
        final_max = np.max(sol.data[2])

        # Peak amplitude should decrease as pulse splits
        assert final_max < initial_max

    def test_em_components_independent(
        self, em_json_path: Path, grid_1d: CartesianGrid
    ) -> None:
        """Test that A_0 and A_1 evolve independently (no coupling)."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        spec = load_equation_system(em_json_path)
        pde = build_pde_from_json(em_json_path)

        # Initialize only A_1, leave A_0 = 0
        pulse = ComponentGaussianPulse(
            center=(50.0,),
            width=5.0,
            amplitude=1.0,
            active_components={"A_1": 1.0},
        )
        initial = pulse.create(grid_1d, spec)

        # Verify A_0 is zero initially
        assert_allclose(initial.data[0], 0.0, atol=1e-10)

        # Run simulation
        result = pde.solve(initial, t_range=5.0, dt=0.01)
        sol = normalize_solve_result(result)

        # A_0 should remain zero (no coupling in the equations)
        assert_allclose(sol.data[0], 0.0, atol=1e-10)

        # A_1 should have evolved
        assert not np.allclose(sol.data[2], initial.data[2])


class TestKGPipelineValidation:
    """Validation tests comparing pipeline KG with existing implementation."""

    def test_load_kg_spec(self, kg_json_path: Path) -> None:
        """Test that KG specification loads correctly."""
        if not kg_json_path.exists():
            pytest.skip(f"Test file not found: {kg_json_path}")

        spec = load_equation_system(kg_json_path)

        assert spec.n_components == 1
        assert spec.component_names == ("phi_0",)

    def test_kg_mass_term(self, kg_json_path: Path, grid_1d: CartesianGrid) -> None:
        """Test that KG mass term causes oscillation."""
        if not kg_json_path.exists():
            pytest.skip(f"Test file not found: {kg_json_path}")

        spec = load_equation_system(kg_json_path)
        pde = build_pde_from_json(kg_json_path)

        # Initialize with uniform phi = 1
        initial = create_initial_state(
            grid_1d,
            spec,
            field_data={"phi_0": np.ones(128)},
        )

        # Run for half period (pi for m=1)
        result = pde.solve(initial, t_range=np.pi, dt=0.01)
        sol = normalize_solve_result(result)

        # Uniform field should oscillate: phi(t) = cos(m*t)
        # At t=pi, phi should be approximately -1
        assert_allclose(sol.data[0], -1.0, atol=0.1)

    def test_kg_pipeline_vs_existing(
        self, kg_json_path: Path, grid_1d_fine: CartesianGrid
    ) -> None:
        """Compare pipeline-generated KG PDE with existing KleinGordonPDE."""
        if not kg_json_path.exists():
            pytest.skip(f"Test file not found: {kg_json_path}")

        # Create PDE from pipeline
        spec = load_equation_system(kg_json_path)
        pipeline_pde = build_pde_from_json(kg_json_path)

        # Create existing KleinGordonPDE
        # Note: The JSON has m^2 = 1, so m = 1
        existing_pde = KleinGordonPDE(KGParameters(mass=1.0))

        # Create same initial conditions for both
        x = cast("np.ndarray", grid_1d_fine.cell_coords[..., 0])
        phi_data = np.exp(-((x - 50) ** 2) / 50)

        # Pipeline state: [phi, pi]
        pipeline_state = create_initial_state(
            grid_1d_fine, spec, field_data={"phi_0": phi_data}
        )

        # Existing state: [phi, pi] with field labels
        existing_state = FieldCollection(
            [
                ScalarField(grid_1d_fine, data=phi_data, label="phi"),
                ScalarField(grid_1d_fine, data=0.0, label="pi"),
            ]
        )

        # Run both for same duration
        t_end = 5.0
        dt = 0.01

        pipeline_result = pipeline_pde.solve(pipeline_state, t_range=t_end, dt=dt)
        existing_result = existing_pde.solve(existing_state, t_range=t_end, dt=dt)
        pipeline_result = normalize_solve_result(pipeline_result)
        existing_result = normalize_solve_result(existing_result)

        # Results should match closely
        # Allow some tolerance for numerical differences
        assert_allclose(
            pipeline_result.data[0],  # phi from pipeline
            existing_result.data[0],  # phi from existing
            rtol=0.01,  # 1% relative tolerance
            atol=0.01,  # 0.01 absolute tolerance
        )


class TestFullPipelineWorkflow:
    """Tests for the complete pipeline workflow."""

    def test_em_full_workflow(self, em_json_path: Path, grid_1d: CartesianGrid) -> None:
        """Test complete EM workflow from JSON to simulation."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        # 1. Load spec
        spec = load_equation_system(em_json_path)

        # 2. Create params
        params = ComponentFieldParams.from_equation_system(spec)
        assert params.is_massless

        # 3. Build PDE
        pde = build_pde_from_json(em_json_path)

        # 4. Create initial state
        state = create_gaussian_pulse_state(
            grid_1d,
            spec,
            center=(50.0,),
            width=5.0,
            component="A_1",
        )

        # 5. Verify state structure
        num_em_states = 4  # A_0, Pi_0, A_1, Pi_1
        assert len(state) == num_em_states
        assert_allclose(state[0].data, 0.0)  # A_0 = 0
        assert np.max(state[2].data) > 0  # A_1 has pulse

        # 6. Run simulation
        result = pde.solve(state, t_range=5.0, dt=0.01)
        sol = normalize_solve_result(result)

        # 7. Verify result
        num_em_states = 4
        assert len(sol) == num_em_states
        # Wave should have propagated
        assert not np.allclose(sol[2].data, state[2].data)

    def test_spec_component_params_consistency(self, em_json_path: Path) -> None:
        """Test that spec and ComponentFieldParams are consistent."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        spec = load_equation_system(em_json_path)
        params = ComponentFieldParams.from_equation_system(spec)

        assert params.n_components == spec.n_components
        assert params.component_names == spec.component_names
        assert params.spatial_dimension == spec.spatial_dimension


class TestEdgeCases:
    """Edge case and error handling tests."""

    def test_missing_json_file(self, tmp_path: Path) -> None:
        """Test handling of missing JSON file."""
        missing = tmp_path / "nonexistent.json"

        with pytest.raises(FileNotFoundError):
            build_pde_from_json(missing)

    def test_wrong_state_size(self, em_json_path: Path, grid_1d: CartesianGrid) -> None:
        """Test that wrong state size raises error."""
        if not em_json_path.exists():
            pytest.skip(f"Test file not found: {em_json_path}")

        pde = build_pde_from_json(em_json_path)

        # Create state with wrong number of fields (should be 4)
        wrong_state = FieldCollection(
            [
                ScalarField(grid_1d, data=0.0),
                ScalarField(grid_1d, data=0.0),
            ]
        )

        with pytest.raises(ValueError, match="Expected 4 fields"):
            pde.evolution_rate(wrong_state)
