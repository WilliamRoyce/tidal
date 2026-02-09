"""Shared pytest fixtures for the test suite."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

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


# ==================== Common Test Values ====================


@pytest.fixture
def gaussian_widths_valid() -> list[float]:
    """Return valid Gaussian widths for testing.

    Returns
    -------
    list[float]
        List of positive width values.
    """
    return [1.0, 2.0, 5.0]


@pytest.fixture
def gaussian_widths_invalid() -> list[Sequence[float]]:
    """Invalid Gaussian widths for testing validation.

    Returns
    -------
    list[Sequence[float]]
        List of invalid width values (zero, negative).
    """
    return [[0.0], [-1.0], [-5.0]]


@pytest.fixture
def masses_coupled() -> list[float]:
    """Return standard masses for coupled field tests.

    Returns
    -------
    list[float]
        Two-field masses [0.25, 1.0].
    """
    return [0.25, 1.0]


@pytest.fixture
def coupling_matrix_symmetric() -> list[list[float]]:
    """Symmetric coupling matrix for coupled field tests.

    Returns
    -------
    list[list[float]]
        2x2 symmetric coupling matrix with off-diagonal elements = 0.2.
    """
    return [[0.0, 0.2], [0.2, 0.0]]


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
