"""Shared pytest fixtures for the test suite."""

from __future__ import annotations

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
def electrostatics_json() -> Path:
    """Path to electrostatics_2d.json, skip if absent."""
    return _cli_json_fixture("electrostatics_2d.json")


@pytest.fixture
def massive_3form_json() -> Path:
    """Path to massive_3form.json, skip if absent."""
    return _cli_json_fixture("massive_3form.json")


# ==================== NPZ Fixtures for tidal measure ====================


@pytest.fixture
def coupled_scalars_npz(coupled_scalars_json: Path, tmp_path: Path) -> Path:
    """Run a short coupled_scalars simulation and save NPZ for measurement tests.

    Returns the path to the generated NPZ file.
    """
    from tidal.cli import main

    output = tmp_path / "coupled_scalars.npz"
    ret = main([
        "simulate", str(coupled_scalars_json),
        "--param", "mPhi2=1.0",
        "--param", "mChi2=4.0",
        "--param", "gCpl=0.5",
        "--t-end", "5.0",
        "--grid-shape", "32",
        "--output", str(output),
    ])
    assert ret == 0, "coupled_scalars simulation failed"
    assert output.exists(), "NPZ file was not created"
    return output
