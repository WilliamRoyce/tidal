"""Shared pytest fixtures for the test suite."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from tidal.solver.operators import set_fd_order, set_spectral

if TYPE_CHECKING:
    from collections.abc import Generator

# ==================== Module-state cleanup ====================


@pytest.fixture(autouse=True)
def _reset_operator_state() -> Generator[None, None, None]:  # pyright: ignore[reportUnusedFunction]
    """Reset FD order and spectral mode after every test.

    CLI tests call set_fd_order(4) (the CLI default) and/or
    set_spectral(True), which persist as module-level state and
    pollute subsequent tests that assume order 2 / FD mode.
    """
    yield
    set_fd_order(2)
    set_spectral(False)


# ==================== CLI JSON Spec Fixtures ====================

_EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "data"


def _cli_json_fixture(name: str) -> Path:
    p = _EXAMPLES_DIR / name
    if not p.exists():
        pytest.skip(f"{name} not found")
    return p


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
        },
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
    "canonical": {
        "hamiltonian_terms": [
            {
                "coefficient": 1.0,
                "factor_a": {"field": "phi_0", "operator": "identity"},
                "factor_b": {"field": "phi_0", "operator": "identity"},
                "coefficient_symbolic": "mPhi2/2",
            },
            {
                "coefficient": 1.0,
                "factor_a": {"field": "chi_0", "operator": "identity"},
                "factor_b": {"field": "chi_0", "operator": "identity"},
                "coefficient_symbolic": "mChi2/2",
            },
            {
                "coefficient": 1.0,
                "factor_a": {"field": "chi_0", "operator": "identity"},
                "factor_b": {"field": "phi_0", "operator": "identity"},
                "coefficient_symbolic": "gCpl",
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "phi_0", "operator": "gradient_x"},
                "factor_b": {"field": "phi_0", "operator": "gradient_x"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "chi_0", "operator": "gradient_x"},
                "factor_b": {"field": "chi_0", "operator": "gradient_x"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                "factor_b": {"field": "phi_0", "operator": "time_derivative"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "chi_0", "operator": "time_derivative"},
                "factor_b": {"field": "chi_0", "operator": "time_derivative"},
            },
        ],
    },
}

_EM_3D_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "lagrangian_expr": "-1/4 F[-a,-b] eta[a,c] eta[b,d] F[-c,-d]",
        "derived_from": "Euler-Lagrange",
        "gauge": "none",
        "linearized": False,
    },
    "spacetime": {
        "dimension": 3,
        "signature": [-1, 1, 1],
        "coordinates": ["t", "x", "y"],
    },
    "fields": [
        {"name": "A_0", "index": 0, "is_dynamical": True},
        {"name": "A_1", "index": 1, "is_dynamical": True},
        {"name": "A_2", "index": 2, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": "A_0",
            "lhs": {"expression": "A_0", "order": {"time": 0, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": -1.0, "operator": "laplacian_x", "field": "A_0"},
                    {"coefficient": 1.0, "operator": "gradient_x", "field": "v_A_1"},
                    {"coefficient": -1.0, "operator": "laplacian_y", "field": "A_0"},
                    {"coefficient": 1.0, "operator": "gradient_y", "field": "v_A_2"},
                ],
            },
            "constraint_solver": {
                "enabled": True,
                "method": "auto",
                "boundary_conditions": {
                    "x": {"type": "periodic"},
                    "y": {"type": "periodic"},
                },
            },
        },
        {
            "field": "A_1",
            "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_1"},
                    {
                        "coefficient": -1.0,
                        "operator": "cross_derivative_xy",
                        "field": "A_2",
                    },
                ],
            },
        },
        {
            "field": "A_2",
            "lhs": {"expression": "d2_t(A_2)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "cross_derivative_xy",
                        "field": "A_1",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_2"},
                ],
            },
        },
    ],
    "coupling": {},
    "canonical": {
        "hamiltonian_terms": [
            {
                "coefficient": -0.5,
                "factor_a": {"field": "A_0", "operator": "gradient_y"},
                "factor_b": {"field": "A_0", "operator": "gradient_y"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "A_1", "operator": "gradient_y"},
                "factor_b": {"field": "A_1", "operator": "gradient_y"},
            },
            {
                "coefficient": -0.5,
                "factor_a": {"field": "A_0", "operator": "gradient_x"},
                "factor_b": {"field": "A_0", "operator": "gradient_x"},
            },
            {
                "coefficient": -1.0,
                "factor_a": {"field": "A_1", "operator": "gradient_y"},
                "factor_b": {"field": "A_2", "operator": "gradient_x"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "A_2", "operator": "gradient_x"},
                "factor_b": {"field": "A_2", "operator": "gradient_x"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "A_1", "operator": "time_derivative"},
                "factor_b": {"field": "A_1", "operator": "time_derivative"},
            },
            {
                "coefficient": 0.5,
                "factor_a": {"field": "A_2", "operator": "time_derivative"},
                "factor_b": {"field": "A_2", "operator": "time_derivative"},
            },
        ],
        "hamiltonian_symbolic": "-Derivative[0, 0, 1][tidalA0][t[], x[], y[]]^2/2 + Derivative[0, 0, 1][tidalA1][t[], x[], y[]]^2/2 - Derivative[0, 1, 0][tidalA0][t[], x[], y[]]^2/2 - Derivative[0, 0, 1][tidalA1][t[], x[], y[]]*Derivative[0, 1, 0][tidalA2][t[], x[], y[]] + Derivative[0, 1, 0][tidalA2][t[], x[], y[]]^2/2 + Derivative[1, 0, 0][tidalA1][t[], x[], y[]]^2/2 + Derivative[1, 0, 0][tidalA2][t[], x[], y[]]^2/2",
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
    tmp_path_factory: pytest.TempPathFactory,
    spec: dict[str, object],
    name: str,
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
        tmp_path_factory,
        _COUPLED_SCALARS_SPEC,
        "coupled_scalars.json",
    )


@pytest.fixture(scope="session")
def inline_em_1d_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Inline 2+1D EM JSON spec — A_0 constraint + A_1/A_2 dynamical (always available)."""
    return _write_inline_json(tmp_path_factory, _EM_3D_SPEC, "em_3d.json")


@pytest.fixture(scope="session")
def inline_constraint_json(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Inline constraint (Poisson) JSON spec (always available)."""
    return _write_inline_json(
        tmp_path_factory,
        _CONSTRAINT_SPEC,
        "electrostatics_2d.json",
    )


# ==================== Directory Fixtures for tidal measure ====================


@pytest.fixture(scope="session")
def coupled_scalars_dir(
    inline_coupled_scalars_json: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Run a short coupled_scalars simulation and save to directory for measurement tests.

    Returns the path to the generated snapshot directory.
    """
    from tidal.cli import main

    output = tmp_path_factory.mktemp("measure") / "coupled_scalars_out"
    ret = main(
        [
            "simulate",
            str(inline_coupled_scalars_json),
            "--param",
            "mPhi2=1.0",
            "--param",
            "mChi2=4.0",
            "--param",
            "gCpl=0.5",
            "--t-end",
            "5.0",
            "--grid-shape",
            "32",
            "--output",
            str(output),
        ],
    )
    assert ret == 0, "coupled_scalars simulation failed"
    assert output.is_dir(), "Snapshot directory was not created"
    return output
