"""Parametrized tests ensuring all example JSONs load and evolve correctly.

Every JSON file in examples/data/ is validated for:
1. Schema compliance (load_equation_system succeeds)
2. PDE construction (PDEFromSpec builds without error)
3. One-step evolution (evolution_rate produces finite output)

This catches regressions that would otherwise be invisible for untested examples.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from pde import CartesianGrid

from tidal.symbolic.json_loader import load_equation_system
from tidal.symbolic.pde_builder import PDEFromSpec, create_initial_state

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "data"

# Default parameters for JSONs with symbolic coefficients.
# Keys are JSON filenames; values are parameter dicts.
_DEFAULT_PARAMS: dict[str, dict[str, float]] = {
    "klein_gordon_1d.json": {"m2": 1.0},
    "klein_gordon_3d.json": {"m2": 1.0},
    "coupled_scalars.json": {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5},
    "chern_simons_3d.json": {"kappa": 0.5},
    "proca_1d.json": {"procaMassSquared": 1.0},
    "de_sitter_kg.json": {"dSm2": 1.0, "dSH": 0.1},
    "polar_kg.json": {"polm2": 1.0},
    "cylindrical_kg.json": {"cylm2": 1.0},
    "sphere_kg.json": {"sphm2": 1.0, "sphR": 1.0},
    "spherical_kg.json": {"spm2": 1.0},
    "conformal_kg_static.json": {},
    "massive_3form.json": {"m2": 1.0},
    "scalar_vector_coupling.json": {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2},
    "massive_gravity_3d.json": {"m2": 1.0},
    "coupled_proca_3d.json": {"mA2": 1.0, "mB2": 2.0, "gcoup": 0.5},
}

# Grid sizes per spatial dimension (small for speed).
_GRID_SIZE: dict[int, list[tuple[float, float]]] = {
    1: [(0.1, 10)],
    2: [(0.1, 10), (0.1, 10)],
    3: [(0.1, 5), (0.1, 5), (0.1, 5)],
}


def _collect_example_jsons() -> list[Path]:
    """Collect all JSON files in the examples/data directory."""
    if not EXAMPLES_DIR.exists():
        return []
    return sorted(EXAMPLES_DIR.glob("*.json"))


@pytest.mark.parametrize(
    "json_path",
    _collect_example_jsons(),
    ids=lambda p: p.stem,
)
class TestExampleJSON:
    """Validate every example JSON loads, builds, and evolves."""

    def test_load_schema(self, json_path: Path) -> None:
        """JSON loads without schema errors."""
        spec = load_equation_system(json_path)
        assert spec.n_components >= 1
        assert spec.dimension >= 2  # at least 1+1D
        assert len(spec.equations) == spec.n_components

    def test_build_and_evolve(self, json_path: Path) -> None:
        """PDE builds and produces finite evolution rates."""
        spec = load_equation_system(json_path)
        params = _DEFAULT_PARAMS.get(json_path.name, {})

        pde = PDEFromSpec(spec, parameters=params)

        spatial_dim = spec.spatial_dimension
        bounds = _GRID_SIZE.get(spatial_dim)
        if bounds is None:
            pytest.skip(f"No grid template for {spatial_dim}D")

        # Use non-periodic grid if any equation has a non-periodic constraint
        # solver (e.g., electrostatics with Dirichlet BCs).
        needs_dirichlet = any(
            eq.constraint_solver.enabled
            and any(
                bc.type == "dirichlet"
                for bc in eq.constraint_solver.boundary_conditions.values()
            )
            for eq in spec.equations
        )
        periodic = not needs_dirichlet

        n_cells = 16  # small for speed
        grid = CartesianGrid(bounds, n_cells, periodic=periodic)
        state = create_initial_state(grid, spec)

        rates = pde.evolution_rate(state, t=0.0)

        # All rate fields must be finite (no NaN/Inf)
        for i, field in enumerate(rates):
            assert np.all(np.isfinite(field.data)), (
                f"Non-finite values in rate field {i} for {json_path.name}"
            )
