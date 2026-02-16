"""Shared pytest fixtures for the test suite."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from tidal.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
)

if TYPE_CHECKING:
    from pde import CartesianGrid


# ==================== Grid Fixtures ====================


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """Create a standard 1D periodic grid for testing.

    Returns
    -------
    CartesianGrid
        1D grid with 32 cells from 0.0 to 10.0, periodic boundary conditions.
    """
    return make_grid(
        GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),), periodic=True)
    )


@pytest.fixture
def grid_1d_large() -> CartesianGrid:
    """Create a large 1D periodic grid for testing.

    Returns
    -------
    CartesianGrid
        1D grid with 512 cells from 0.0 to 200.0, periodic boundary conditions.
    """
    return make_grid(
        GridConfig(dim=1, shape=(512,), bounds=((0.0, 200.0),), periodic=True)
    )


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """Create a standard 2D periodic grid for testing.

    Returns
    -------
    CartesianGrid
        2D grid with 16x16 cells from (0.0, 0.0) to (10.0, 5.0), periodic.
    """
    return make_grid(
        GridConfig(
            dim=2, shape=(16, 16), bounds=((0.0, 10.0), (0.0, 5.0)), periodic=True
        )
    )


@pytest.fixture
def grid_2d_large() -> CartesianGrid:
    """Create a large 2D periodic grid for testing.

    Returns
    -------
    CartesianGrid
        2D grid with 128x128 cells from (0.0, 0.0) to (200.0, 200.0), periodic.
    """
    return make_grid(
        GridConfig(
            dim=2, shape=(128, 128), bounds=((0.0, 200.0), (0.0, 200.0)), periodic=True
        )
    )


# ==================== Configuration Fixtures ====================


@pytest.fixture
def grid_config_1d() -> GridConfig:
    """Create a standard 1D grid configuration.

    Returns
    -------
    GridConfig
        Configuration for 1D grid with 32 cells from 0.0 to 10.0, periodic.
    """
    return GridConfig(dim=1, shape=(32,), bounds=((0.0, 10.0),), periodic=True)


@pytest.fixture
def grid_config_2d() -> GridConfig:
    """Create a standard 2D grid configuration.

    Returns
    -------
    GridConfig
        Configuration for 2D grid with 16x16 cells, periodic.
    """
    return GridConfig(
        dim=2, shape=(16, 16), bounds=((0.0, 10.0), (0.0, 5.0)), periodic=True
    )


@pytest.fixture
def kg_params() -> KGParameters:
    """Create standard Klein-Gordon parameters.

    Returns
    -------
    KGParameters
        Klein-Gordon parameters with mass=0.5.
    """
    return KGParameters(mass=0.5)


@pytest.fixture
def kg_pde(kg_params: KGParameters) -> KleinGordonPDE:
    """Create a Klein-Gordon PDE instance.

    Parameters
    ----------
    kg_params : KGParameters
        Klein-Gordon parameters fixture.

    Returns
    -------
    KleinGordonPDE
        Klein-Gordon PDE with mass=0.5.
    """
    return KleinGordonPDE(params=kg_params)


@pytest.fixture
def simulation_config_fast() -> SimulationConfig:
    """Create a simulation configuration for fast tests.

    Returns
    -------
    SimulationConfig
        Configuration for quick simulation (t_end=5.0, scipy solver, no progress).
    """
    return SimulationConfig(
        t_end=5.0,
        dt=None,  # Adaptive
        solver="scipy",
        backend="numpy",
        progress=False,
    )


@pytest.fixture
def simulation_config_explicit() -> SimulationConfig:
    """Create a simulation configuration for explicit time stepping.

    Returns
    -------
    SimulationConfig
        Configuration with explicit solver, fixed dt=0.1.
    """
    return SimulationConfig(
        t_end=5.0,
        dt=0.1,
        solver="explicit",
        backend="numpy",
        progress=False,
    )


# ==================== CLI JSON Spec Fixtures ====================

_EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "data"


def _cli_json_fixture(name: str) -> Path:
    p = _EXAMPLES_DIR / name
    if not p.exists():
        pytest.skip(f"{name} not found")
    return p


@pytest.fixture
def klein_gordon_1d_json() -> Path:
    """Path to klein_gordon_1d.json, skip if absent."""
    return _cli_json_fixture("klein_gordon_1d.json")


@pytest.fixture
def klein_gordon_3d_json() -> Path:
    """Path to klein_gordon_3d.json, skip if absent."""
    return _cli_json_fixture("klein_gordon_3d.json")


@pytest.fixture
def polar_kg_json() -> Path:
    """Path to polar_kg.json, skip if absent."""
    return _cli_json_fixture("polar_kg.json")


@pytest.fixture
def chern_simons_json() -> Path:
    """Path to chern_simons_3d.json, skip if absent."""
    return _cli_json_fixture("chern_simons_3d.json")


@pytest.fixture
def coupled_scalars_json() -> Path:
    """Path to coupled_scalars.json, skip if absent."""
    return _cli_json_fixture("coupled_scalars.json")


@pytest.fixture
def massive_3form_json() -> Path:
    """Path to massive_3form.json, skip if absent."""
    return _cli_json_fixture("massive_3form.json")


# ==================== Inline JSON Spec Data ====================
#
# Self-contained JSON specs for core tests that must always run,
# regardless of whether examples/data/ has been populated via ``tidal derive``.

_KG_1D_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "lagrangian_expr": "-1/2 (d phi)^2 - m2/2 phi^2",
        "derived_from": "Euler-Lagrange",
        "gauge": "none",
        "linearized": False,
        "parameters": {"m2": 1.0},
    },
    "spacetime": {
        "dimension": 2,
        "signature": [-1, 1],
        "coordinates": ["t", "x"],
    },
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        }
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2"]]},
}

_COUPLED_SCALARS_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "lagrangian_expr": "-1/2 (d phi)^2 - mPhi2/2 phi^2 - 1/2 (d chi)^2 - mChi2/2 chi^2 - gCpl phi chi",
        "derived_from": "Euler-Lagrange",
        "gauge": "none",
        "linearized": False,
        "parameters": {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5},
    },
    "spacetime": {
        "dimension": 2,
        "signature": [-1, 1],
        "coordinates": ["t", "x"],
    },
    "fields": [
        {"name": "phi_0", "index": 0, "is_dynamical": True},
        {"name": "chi_0", "index": 1, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-mPhi2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        },
        {
            "field": "chi_0",
            "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-mChi2",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                ],
            },
        },
    ],
    "coupling": {
        "mass_matrix_symbolic": [["-mPhi2", None], [None, "-mChi2"]],
        "coupling_matrix_symbolic": [[None, "-gCpl"], ["-gCpl", None]],
    },
}

_CONSTRAINT_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "lagrangian_expr": "1/2 (nabla phi)^2 + rho * phi",
        "derived_from": "Euler-Lagrange",
        "gauge": "none",
        "linearized": True,
        "construction": "analytical",
    },
    "spacetime": {
        "dimension": 3,
        "signature": [-1, 1, 1],
        "coordinates": ["t", "x", "y"],
    },
    "fields": [
        {"name": "phi", "index": 0, "is_dynamical": True},
        {"name": "rho", "index": 1, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": "phi",
            "lhs": {"expression": "phi", "order": {"time": 0, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
                    {"coefficient": 1.0, "operator": "identity", "field": "rho"},
                ],
            },
            "constraint_solver": {
                "enabled": True,
                "method": "poisson",
                "boundary_conditions": {
                    "x": {"type": "dirichlet", "value": 0.0},
                    "y": {"type": "dirichlet", "value": 0.0},
                },
            },
        },
        {
            "field": "rho",
            "lhs": {"expression": "d2_t(rho)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 0.0, "operator": "identity", "field": "rho"},
                ],
            },
        },
    ],
    "coupling": {},
}


def _write_inline_json(
    tmp_path_factory: pytest.TempPathFactory, spec: dict[str, object], name: str
) -> Path:
    """Write an inline JSON spec dict to a temp file and return its path."""
    d = tmp_path_factory.getbasetemp() / "inline_specs"
    d.mkdir(exist_ok=True)
    p = d / name
    if not p.exists():
        p.write_text(json.dumps(spec, indent=2))
    return p


@pytest.fixture(scope="session")
def inline_kg_1d_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Inline 1D Klein-Gordon JSON spec (always available)."""
    return _write_inline_json(tmp_path_factory, _KG_1D_SPEC, "klein_gordon_1d.json")


@pytest.fixture(scope="session")
def inline_coupled_scalars_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Inline coupled-scalars JSON spec (always available)."""
    return _write_inline_json(
        tmp_path_factory, _COUPLED_SCALARS_SPEC, "coupled_scalars.json"
    )


@pytest.fixture(scope="session")
def inline_constraint_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Inline constraint (Poisson) JSON spec (always available)."""
    return _write_inline_json(
        tmp_path_factory, _CONSTRAINT_SPEC, "electrostatics_2d.json"
    )


# ==================== Directory Fixtures for tidal measure ====================


@pytest.fixture
def coupled_scalars_dir(inline_coupled_scalars_json: Path, tmp_path: Path) -> Path:
    """Run a short coupled_scalars simulation and save to directory for measurement tests.

    Returns the path to the generated snapshot directory.
    """
    from tidal.cli import main

    output = tmp_path / "coupled_scalars_out"
    ret = main([
        "simulate", str(inline_coupled_scalars_json),
        "--param", "mPhi2=1.0",
        "--param", "mChi2=4.0",
        "--param", "gCpl=0.5",
        "--t-end", "5.0",
        "--grid-shape", "32",
        "--output", str(output),
    ])
    assert ret == 0, "coupled_scalars simulation failed"
    assert output.is_dir(), "Snapshot directory was not created"
    return output
