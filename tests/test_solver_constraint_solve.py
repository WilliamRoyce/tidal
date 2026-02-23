"""Tests for tidal.solver.constraint_solve — constraint pre-solve module."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.constraint_solve import (
    _classify_terms,
    _ConstraintTerms,
    _fft_solve_single,
    _matrix_solve,
    _probe_operator_matrix,
    _select_method,
    pre_solve_constraints,
)
from tidal.solver.grid import GridInfo
from tidal.solver.operators import AxisBCSpec, apply_operator, laplacian
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import (
    ConstraintSolverConfig,
    EquationSystem,
)

# ---------------------------------------------------------------------------
# Helper: inline spec builders
# ---------------------------------------------------------------------------


def _make_em_2d_spec() -> EquationSystem:
    """Build EM 2+1D spec inline (Maxwell A_0 constraint + A_1, A_2 wave eqs)."""
    data: dict[str, Any] = {
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
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "pi_1",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "pi_2",
                        },
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
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "pi_0",
                        },
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
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_2"},
                        {
                            "coefficient": -1.0,
                            "operator": "cross_derivative_xy",
                            "field": "A_1",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "pi_0",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "A_1": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_1"},
                    {"coefficient": 1.0, "operator": "gradient_x", "field": "A_0"},
                ],
                "A_2": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_2"},
                    {"coefficient": 1.0, "operator": "gradient_y", "field": "A_0"},
                ],
            },
            "kinetic_matrix": {
                "entries": [
                    {"i": 0, "j": 0, "value": 1.0},
                    {"i": 1, "j": 1, "value": 1.0},
                ],
                "dimension": 2,
            },
            "spatial_momenta": {
                "A_1": [
                    {"coefficient": -1.0, "operator": "gradient_x", "field": "A_0"}
                ],
                "A_2": [
                    {"coefficient": -1.0, "operator": "gradient_y", "field": "A_0"}
                ],
            },
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


def _make_chern_simons_spec() -> EquationSystem:
    """Build Chern-Simons 2+1D spec inline (A_0 constraint + A_1, A_2 with kappa coupling)."""
    data: dict[str, Any] = {
        "metadata": {"parameters": {"kappa": 0.5}},
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
                        {
                            "coefficient": 1.0,
                            "operator": "gradient_x",
                            "field": "A_2",
                            "coefficient_symbolic": "kappa",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "pi_1",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "A_1",
                            "coefficient_symbolic": "-kappa",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "pi_2",
                        },
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
                        {
                            "coefficient": 1.0,
                            "operator": "first_derivative_t",
                            "field": "A_2",
                            "coefficient_symbolic": "kappa",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "pi_0",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "A_0",
                            "coefficient_symbolic": "-kappa",
                        },
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
                            "coefficient": 1.0,
                            "operator": "gradient_x",
                            "field": "A_0",
                            "coefficient_symbolic": "kappa",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_2"},
                        {
                            "coefficient": -1.0,
                            "operator": "first_derivative_t",
                            "field": "A_1",
                            "coefficient_symbolic": "-kappa",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "cross_derivative_xy",
                            "field": "A_1",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "pi_0",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "A_1": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_1"},
                    {
                        "coefficient": -0.5,
                        "operator": "identity",
                        "field": "A_2",
                        "coefficient_symbolic": "-1/2*kappa",
                    },
                    {"coefficient": 1.0, "operator": "gradient_x", "field": "A_0"},
                ],
                "A_2": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_2"},
                    {
                        "coefficient": 0.5,
                        "operator": "identity",
                        "field": "A_1",
                        "coefficient_symbolic": "kappa/2",
                    },
                    {"coefficient": 1.0, "operator": "gradient_y", "field": "A_0"},
                ],
            },
            "kinetic_matrix": {
                "entries": [
                    {"i": 0, "j": 0, "value": 1.0},
                    {"i": 1, "j": 1, "value": 1.0},
                ],
                "dimension": 2,
            },
            "spatial_momenta": {
                "A_1": [
                    {"coefficient": -1.0, "operator": "gradient_x", "field": "A_0"}
                ],
                "A_2": [
                    {"coefficient": -1.0, "operator": "gradient_y", "field": "A_0"}
                ],
            },
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


def _make_kg_spec() -> EquationSystem:
    """Klein-Gordon: no constraints, purely dynamical."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            }
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "phi_0": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_phi_0"}
                ],
            },
            "kinetic_matrix": {
                "entries": [{"i": 0, "j": 0, "value": 1.0}],
                "dimension": 1,
            },
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


def _make_helmholtz_spec() -> EquationSystem:
    """Helmholtz constraint: (laplacian - m²)phi_0 = source.

    Non-singular (m² ≠ 0), so no zero-mode issue.
    """
    data: dict[str, Any] = {
        "spacetime": {
            "dimension": 3,
            "signature": [-1, 1, 1],
            "coordinates": ["t", "x", "y"],
        },
        "fields": [
            {"name": "phi_0", "index": 0},
            {"name": "rho_0", "index": 1},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "phi_0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                        {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                        {"coefficient": 1.0, "operator": "identity", "field": "rho_0"},
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
                "field": "rho_0",
                "lhs": {"expression": "d2_t(rho_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 0.0, "operator": "identity", "field": "rho_0"},
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "rho_0": [
                    {"coefficient": 1.0, "operator": "identity", "field": "pi_1"}
                ],
            },
            "kinetic_matrix": {
                "entries": [{"i": 0, "j": 0, "value": 1.0}],
                "dimension": 1,
            },
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# Unit tests: term classification
# ---------------------------------------------------------------------------


class TestClassifyTerms:
    def test_self_vs_source(self) -> None:
        """Self-terms reference constraint field; source-terms reference others."""
        from tidal.solver.coefficients import CoefficientEvaluator

        spec = _make_em_2d_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        coeff_eval = CoefficientEvaluator(spec, grid)

        eq = spec.equations[0]
        terms = _classify_terms(
            0, eq.rhs_terms, "A_0", coeff_eval, 0.0, eq.constraint_solver
        )

        # Self-terms: laplacian_x(A_0) + laplacian_y(A_0)
        assert len(terms.self_terms) == 2
        for coeff, op in terms.self_terms:
            assert op in {"laplacian_x", "laplacian_y"}
            assert coeff == 1.0

        # Source terms: gradient_x(pi_1) + gradient_y(pi_2)
        assert len(terms.source_terms) == 2
        for _coeff, op, field in terms.source_terms:
            assert op in {"gradient_x", "gradient_y"}
            assert field in {"pi_1", "pi_2"}

    def test_position_dependent_detection(self) -> None:
        """has_position_dependent_self should be False for constant coeffs."""
        from tidal.solver.coefficients import CoefficientEvaluator

        spec = _make_em_2d_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        coeff_eval = CoefficientEvaluator(spec, grid)

        eq = spec.equations[0]
        terms = _classify_terms(
            0, eq.rhs_terms, "A_0", coeff_eval, 0.0, eq.constraint_solver
        )
        assert not terms.has_position_dependent_self


# ---------------------------------------------------------------------------
# Unit tests: method selection
# ---------------------------------------------------------------------------


class TestSelectMethod:
    def _make_terms(
        self, *, periodic: bool = True, pos_dep: bool = False
    ) -> tuple[_ConstraintTerms, GridInfo]:
        grid = GridInfo(
            bounds=((0, 10), (0, 10)),
            shape=(8, 8),
            periodic=(periodic, periodic),
        )
        coeff: float | np.ndarray = 1.0
        if pos_dep:
            coeff = np.ones(grid.shape)

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="A_0",
            self_terms=[(coeff, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=pos_dep,
            config=config,
        )
        return terms, grid

    def test_fft_for_periodic_constant(self) -> None:
        terms, grid = self._make_terms(periodic=True, pos_dep=False)
        assert _select_method(terms, grid, None) == "fft"

    def test_matrix_for_nonperiodic(self) -> None:
        terms, grid = self._make_terms(periodic=False, pos_dep=False)
        assert _select_method(terms, grid, None) == "matrix"

    def test_matrix_for_position_dependent(self) -> None:
        terms, grid = self._make_terms(periodic=True, pos_dep=True)
        assert _select_method(terms, grid, None) == "matrix"

    def test_user_override(self) -> None:
        """Explicit method in config overrides auto-selection."""
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        config = ConstraintSolverConfig(enabled=True, method="matrix")
        terms = _ConstraintTerms(
            field_name="A_0",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )
        assert _select_method(terms, grid, None) == "matrix"

    def test_bc_override_nonperiodic(self) -> None:
        """Explicit non-periodic bc forces matrix even if grid is periodic."""
        terms, grid = self._make_terms(periodic=True, pos_dep=False)
        assert _select_method(terms, grid, ("neumann", "neumann")) == "matrix"

    def test_fft_with_axis_bc_spec_periodic(self) -> None:
        """AxisBCSpec periodic objects should route to FFT, not matrix."""
        terms, grid = self._make_terms(periodic=True, pos_dep=False)
        bc = (AxisBCSpec(periodic=True), AxisBCSpec(periodic=True))
        assert _select_method(terms, grid, bc) == "fft"


# ---------------------------------------------------------------------------
# Unit tests: FFT solver
# ---------------------------------------------------------------------------


class TestFFTSolver:
    def test_poisson_1d_residual(self) -> None:
        """Solve laplacian(u) = f on periodic 1D: verify FD residual is zero.

        Uses modified wavenumbers (exact for the FD stencil), so the solution
        satisfies the discrete Laplacian exactly, not the continuous one.
        """
        n = 64
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))
        x = grid.axes_coords(0)

        # Source: f = sin(2x), zero-mean
        f = np.sin(2 * x)

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        u_solved = _fft_solve_single(terms, grid, -f)

        # Verify the DISCRETE residual: laplacian_FD(u) - f ≈ 0
        residual = laplacian(u_solved, grid, "periodic") - f
        assert float(np.max(np.abs(residual))) < 1e-12

    def test_poisson_2d_residual(self) -> None:
        """Solve laplacian(u) = f on periodic 2D: verify FD residual is zero."""
        nx, ny = 32, 32
        lx, ly = 2 * np.pi, 2 * np.pi
        grid = GridInfo(
            bounds=((0, lx), (0, ly)), shape=(nx, ny), periodic=(True, True)
        )
        x, y = grid.coord_arrays()

        f = np.sin(x) * np.cos(2 * y)  # zero-mean

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        u_solved = _fft_solve_single(terms, grid, -f)

        residual = laplacian(u_solved, grid, "periodic") - f
        assert float(np.max(np.abs(residual))) < 1e-10

    def test_helmholtz_residual(self) -> None:
        """Solve (laplacian - m²)u = f: verify residual is zero."""
        n = 64
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))
        x = grid.axes_coords(0)

        m2 = 1.0
        f = np.sin(3 * x)

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian"), (-m2, "identity")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        u_solved = _fft_solve_single(terms, grid, -f)

        # Verify: laplacian(u) - m²*u - f ≈ 0
        residual = laplacian(u_solved, grid, "periodic") - m2 * u_solved - f
        assert float(np.max(np.abs(residual))) < 1e-12

    def test_zero_mode_handling(self) -> None:
        """Poisson with periodic BCs: zero-mean source → zero-mean solution."""
        n = 32
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))
        x = grid.axes_coords(0)

        f = np.cos(2 * x)  # zero mean

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        u = _fft_solve_single(terms, grid, -f)
        assert abs(np.mean(u)) < 1e-12

    def test_incompatible_source_raises(self) -> None:
        """Poisson with nonzero-mean source must raise ValueError."""
        n = 32
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))

        source = np.ones(grid.shape)  # nonzero mean → incompatible

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        with pytest.raises(ValueError, match="incompatible"):
            _fft_solve_single(terms, grid, source)


# ---------------------------------------------------------------------------
# Unit tests: operator probing (Tier 2)
# ---------------------------------------------------------------------------


class TestProbeMatrix:
    def test_laplacian_1d_periodic(self) -> None:
        """Probe-built Laplacian matrix matches expected tridiagonal structure."""
        n = 8
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))
        dx = grid.dx[0]

        mat = _probe_operator_matrix([(1.0, "laplacian")], grid, "periodic")

        # Dense version for inspection
        md: Any = mat.toarray()  # pyright: ignore[reportUnknownVariableType]
        # Diagonal should be -2/dx²
        np.testing.assert_allclose(np.diag(md), -2.0 / dx**2, atol=1e-12)  # pyright: ignore[reportUnknownArgumentType]
        # Off-diagonals should be 1/dx²
        for i in range(n):
            assert abs(md[i, (i + 1) % n] - 1.0 / dx**2) < 1e-12  # pyright: ignore[reportUnknownArgumentType]
            assert abs(md[i, (i - 1) % n] - 1.0 / dx**2) < 1e-12  # pyright: ignore[reportUnknownArgumentType]

    def test_probe_matches_fft_for_poisson(self) -> None:
        """Cross-validate: probing and FFT give same solution up to a constant.

        The Laplacian with periodic BCs has a one-dimensional null space
        (constants).  FFT pins the zero mode to 0 (zero-mean solution),
        while spsolve picks an arbitrary particular solution.  After
        removing the mean, both should agree to machine precision.
        """
        n = 32
        domain = 2 * np.pi
        grid = GridInfo(bounds=((0, domain),), shape=(n,), periodic=(True,))
        x = grid.axes_coords(0)

        f = np.sin(3 * x)

        config = ConstraintSolverConfig(enabled=True, method="auto")
        terms = _ConstraintTerms(
            field_name="u",
            self_terms=[(1.0, "laplacian")],
            source_terms=[],
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

        u_fft = _fft_solve_single(terms, grid, -f)

        mat = _probe_operator_matrix(terms.self_terms, grid, "periodic")
        u_mat = _matrix_solve(mat, -f, grid.shape)

        # Remove mean (gauge freedom for singular Poisson operator)
        u_fft_zm = u_fft - np.mean(u_fft)
        u_mat_zm = u_mat - np.mean(u_mat)
        np.testing.assert_allclose(u_mat_zm, u_fft_zm, atol=1e-10)

    def test_position_dependent_mass(self) -> None:
        """Matrix solve with spatially-varying mass: verify residual is small."""
        nx, ny = 16, 16
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(nx, ny),
            periodic=(True, True),
        )
        x_arr, y_arr = grid.coord_arrays()

        # m²(x) = 1 + 0.5*cos(x) — position-dependent
        m2_arr = 1.0 + 0.5 * np.cos(x_arr)

        self_terms: list[tuple[float | np.ndarray, str]] = [
            (1.0, "laplacian"),
            (-m2_arr, "identity"),  # NDArray coefficient
        ]

        # Build source: f = sin(x)*cos(y)
        source_rhs = np.sin(x_arr) * np.cos(y_arr)

        mat = _probe_operator_matrix(self_terms, grid, "periodic")
        u = _matrix_solve(mat, source_rhs, grid.shape)

        # Verify: apply operator to solution, check residual
        result = laplacian(u, grid, "periodic") - m2_arr * u
        residual = result + source_rhs  # should be ≈ 0
        assert float(np.max(np.abs(residual))) < 1e-10

    def test_neumann_bc(self) -> None:
        """Matrix solver with Neumann BCs produces valid solution."""
        n = 32
        grid = GridInfo(bounds=((0, 10),), shape=(n,), periodic=(False,))

        # Helmholtz: (laplacian - 1)u = f
        self_terms: list[tuple[float | np.ndarray, str]] = [
            (1.0, "laplacian"),
            (-1.0, "identity"),
        ]
        x = grid.axes_coords(0)
        source_rhs = np.sin(np.pi * x / 10)

        mat = _probe_operator_matrix(self_terms, grid, "neumann")
        u = _matrix_solve(mat, source_rhs, grid.shape)

        # Verify residual
        result = apply_operator("laplacian", u, grid, "neumann") - u
        residual = result + source_rhs
        assert float(np.max(np.abs(residual))) < 1e-10


# ---------------------------------------------------------------------------
# Integration tests: pre_solve_constraints
# ---------------------------------------------------------------------------


class TestPreSolveConstraints:
    def test_noop_for_kg(self) -> None:
        """KG spec (no constraints) returns y0 unchanged."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.random.default_rng(42).standard_normal(layout.total_size)

        y0_out = pre_solve_constraints(spec, grid, y0)
        np.testing.assert_array_equal(y0_out, y0)

    def test_noop_when_disabled(self) -> None:
        """Constraint with enabled=False is not solved."""
        data: dict[str, Any] = {
            "spacetime": {
                "dimension": 3,
                "signature": [-1, 1, 1],
                "coordinates": ["t", "x", "y"],
            },
            "fields": [
                {"name": "A_0", "index": 0},
                {"name": "A_1", "index": 1},
            ],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "A_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "A_0",
                            },
                        ],
                    },
                    "constraint_solver": {"enabled": False},
                },
                {
                    "field": "A_1",
                    "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "A_1",
                            },
                        ],
                    },
                },
            ],
            "canonical": {
                "hamiltonian_terms": [],
                "field_rates": {
                    "A_1": [
                        {"coefficient": 1.0, "operator": "identity", "field": "pi_1"}
                    ],
                },
                "kinetic_matrix": {
                    "entries": [{"i": 0, "j": 0, "value": 1.0}],
                    "dimension": 1,
                },
                "spatial_momenta": {},
                "hamiltonian_symbolic": "test",
            },
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.random.default_rng(42).standard_normal(layout.total_size)

        y0_out = pre_solve_constraints(spec, grid, y0)
        np.testing.assert_array_equal(y0_out, y0)

    def test_em_constraint_presolve(self) -> None:
        """EM constraint: with zero pi, A_0 should be zero."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)

        # All zeros: pi_1 = pi_2 = 0 → source = 0 → A_0 = 0
        y0 = np.zeros(layout.total_size)
        y0_out = pre_solve_constraints(spec, grid, y0, bc="periodic")

        # A_0 should remain zero
        a0_slot = layout.field_slot_map["A_0"]
        n = grid.num_points
        a0 = y0_out[a0_slot * n : (a0_slot + 1) * n]
        np.testing.assert_allclose(a0, 0.0, atol=1e-12)

    def test_em_constraint_nonzero_source(self) -> None:
        """EM constraint with nonzero pi: A_0 should be nonzero and satisfy residual."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        # Set pi_1 = cos(x)*cos(y) (zero-mean for compatibility)
        pi1_slot = layout.momentum_slot_map["A_1"]
        y0[pi1_slot * n : (pi1_slot + 1) * n] = (np.cos(x_arr) * np.cos(y_arr)).ravel()

        y0_out = pre_solve_constraints(spec, grid, y0, bc="periodic")

        # A_0 should be nonzero
        a0_slot = layout.field_slot_map["A_0"]
        a0 = y0_out[a0_slot * n : (a0_slot + 1) * n].reshape(grid.shape)
        assert float(np.max(np.abs(a0))) > 1e-6

        # Verify residual: laplacian(A_0) + source ≈ 0
        pi1 = y0_out[pi1_slot * n : (pi1_slot + 1) * n].reshape(grid.shape)
        lap_a0 = apply_operator("laplacian", a0, grid, "periodic")
        source = -apply_operator("gradient_x", pi1, grid, "periodic")
        residual = lap_a0 + source
        assert float(np.max(np.abs(residual))) < 1e-10

    def test_chern_simons_presolve(self) -> None:
        """Chern-Simons: Gaussian IC on A_1 → A_0 nonzero after pre-solve."""
        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        # Set A_1 = Gaussian centered at (25, 25)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        y0_out = pre_solve_constraints(
            spec, grid, y0, bc="periodic", parameters={"kappa": 0.5}
        )

        # A_0 should be nonzero (kappa source from A_1)
        a0_slot = layout.field_slot_map["A_0"]
        a0 = y0_out[a0_slot * n : (a0_slot + 1) * n].reshape(grid.shape)
        assert float(np.max(np.abs(a0))) > 1e-6

        # A_1 should be unchanged (Gaussian)
        a1_out = y0_out[a1_slot * n : (a1_slot + 1) * n].reshape(grid.shape)
        np.testing.assert_allclose(a1_out, gaussian, atol=1e-14)

    def test_helmholtz_presolve(self) -> None:
        """Helmholtz constraint (laplacian - m²)phi_0 = source."""
        spec = _make_helmholtz_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, _ = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        # Set rho_0 = cos(x) (zero-mean)
        rho_slot = layout.field_slot_map["rho_0"]
        y0[rho_slot * n : (rho_slot + 1) * n] = np.cos(x_arr).ravel()

        y0_out = pre_solve_constraints(spec, grid, y0, bc="periodic")

        phi_slot = layout.field_slot_map["phi_0"]
        phi = y0_out[phi_slot * n : (phi_slot + 1) * n].reshape(grid.shape)

        # Verify residual: laplacian(phi) - phi + rho ≈ 0
        rho = y0_out[rho_slot * n : (rho_slot + 1) * n].reshape(grid.shape)
        residual = laplacian(phi, grid, "periodic") - phi + rho
        assert float(np.max(np.abs(residual))) < 1e-10

    def test_does_not_mutate_input(self) -> None:
        """pre_solve_constraints should not modify the input y0."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.ones(layout.total_size)
        y0_copy = y0.copy()

        pre_solve_constraints(spec, grid, y0, bc="periodic")
        np.testing.assert_array_equal(y0, y0_copy)


# ---------------------------------------------------------------------------
# Integration test: IDA with pre-solve
# ---------------------------------------------------------------------------


def _has_sundials() -> bool:
    """Check if sksundae is importable."""
    try:
        import sksundae  # noqa: F401  # pyright: ignore[reportMissingTypeStubs, reportUnusedImport]
    except ImportError:
        return False
    else:
        return True


class TestIDAWithPreSolve:
    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_helmholtz_ida_succeeds(self) -> None:
        """IDA solve with Helmholtz constraint pre-solve succeeds.

        Uses a Helmholtz (laplacian - m²) constraint, which is non-singular
        and thus has a well-conditioned Jacobian for IDA.
        """
        from tidal.solver.ida import solve_ida

        spec = _make_helmholtz_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, _ = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        # Set rho_0 = cos(x) (source for constraint)
        rho_slot = layout.field_slot_map["rho_0"]
        y0[rho_slot * n : (rho_slot + 1) * n] = np.cos(x_arr).ravel()

        result = solve_ida(
            spec,
            grid,
            y0,
            (0.0, 1.0),
            bc="periodic",
            num_snapshots=5,
        )

        assert result["success"], f"IDA failed: {result['message']}"
        assert len(result["t"]) == 5

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_chern_simons_ida_succeeds(self) -> None:
        """Full IDA solve on Chern-Simons: pre-solve + gauge fixing.

        CS has a Poisson constraint (singular Laplacian + periodic BCs).
        The constraint pre-solve finds consistent y0, and the IDA residual
        pins mean(A_0)=0 to make the Jacobian non-singular.  This is the
        standard gauge-fixing approach used in FEniCS/Firedrake.
        """
        from tidal.solver.ida import solve_ida

        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        result = solve_ida(
            spec,
            grid,
            y0,
            (0.0, 1.0),
            bc="periodic",
            parameters={"kappa": 0.5},
            num_snapshots=5,
        )

        assert result["success"], f"IDA failed: {result['message']}"
        assert len(result["t"]) == 5

        # Gauge fixing: A_0[0] pinned to zero throughout evolution
        a0_slot = layout.field_slot_map["A_0"]
        a0_final = result["y"][-1][a0_slot * n : (a0_slot + 1) * n]
        assert abs(a0_final[0]) < 1e-10


class TestGaugeRegularizationWarnings:
    """Verify that gauge regularization emits user-visible warnings."""

    def test_fft_singular_mode_warns(self) -> None:
        """FFT pre-solve emits UserWarning when zero-mode is regularized."""
        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        with pytest.warns(UserWarning, match="singular mode"):
            pre_solve_constraints(
                spec, grid, y0, bc="periodic", parameters={"kappa": 0.5}
            )

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_gauge_pin_warns(self) -> None:
        """IDA solver emits UserWarning for gauge regularization (point-pinning)."""
        from tidal.solver.ida import solve_ida

        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        with pytest.warns(UserWarning, match="pinning A_0"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.1),
                bc="periodic",
                parameters={"kappa": 0.5},
                num_snapshots=3,
            )

    def test_fft_with_axis_bc_spec_periodic(self) -> None:
        """AxisBCSpec periodic objects should still trigger FFT pre-solve."""
        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        # Pass AxisBCSpec instead of string — should still find FFT path
        bc = (AxisBCSpec(periodic=True), AxisBCSpec(periodic=True))
        with pytest.warns(UserWarning, match="singular mode"):
            pre_solve_constraints(spec, grid, y0, bc=bc, parameters={"kappa": 0.5})

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_gauge_fix_with_axis_bc_spec(self) -> None:
        """IDA gauge-fix detection works with AxisBCSpec periodic BCs."""
        from tidal.solver.ida import solve_ida

        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        gaussian = np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / (2 * 5**2))
        a1_slot = layout.field_slot_map["A_1"]
        y0[a1_slot * n : (a1_slot + 1) * n] = gaussian.ravel()

        bc = (AxisBCSpec(periodic=True), AxisBCSpec(periodic=True))
        with pytest.warns(UserWarning, match="pinning A_0"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.1),
                bc=bc,
                parameters={"kappa": 0.5},
                num_snapshots=3,
            )

    def test_helmholtz_no_gauge_warning(self) -> None:
        """Helmholtz constraint (non-singular) emits no gauge regularization warning."""
        spec = _make_helmholtz_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, _ = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        rho_slot = layout.field_slot_map["rho_0"]
        y0[rho_slot * n : (rho_slot + 1) * n] = np.cos(x_arr).ravel()

        import warnings as _warnings

        with _warnings.catch_warnings():
            _warnings.simplefilter("error", UserWarning)
            # Should NOT warn — Helmholtz has mass term, not singular
            pre_solve_constraints(spec, grid, y0, bc="periodic")


# ---------------------------------------------------------------------------
# No-self-term constraint tests (IDA freezes field at zero)
# ---------------------------------------------------------------------------


def _make_no_self_term_spec() -> EquationSystem:
    """Spec with a constraint equation that has no self-referencing terms.

    Models the situation from linearized gravity: constraint field h_0 has
    an equation involving only gradients of h_1 (a dynamical field),
    not h_0 itself.  h_1 is a dynamical wave field.
    """
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [
            {"name": "h_0", "index": 0},
            {"name": "h_1", "index": 1},
        ],
        "equations": [
            {
                "field": "h_0",
                "lhs": {"expression": "h_0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        # Only references h_1, NOT h_0
                        {
                            "coefficient": -0.5,
                            "operator": "gradient_x",
                            "field": "h_1",
                        },
                    ],
                },
            },
            {
                "field": "h_1",
                "lhs": {
                    "expression": "d2_t(h_1)",
                    "order": {"time": 2},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian",
                            "field": "h_1",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "h_1": [
                    {
                        "coefficient": 1.0,
                        "operator": "identity",
                        "field": "pi_1",
                    },
                ],
            },
            "kinetic_matrix": {
                "entries": [{"i": 0, "j": 0, "value": 1.0}],
                "dimension": 1,
            },
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


def _make_multi_no_self_term_spec() -> EquationSystem:
    """Spec with TWO no-self-term constraint equations.

    Models the gravitational-waves-like scenario: two constraint fields
    (c_0, c_1) each constrain a dynamical field (h_0) via gradients, but
    neither constraint field appears in its own equation.

    c_0's constraint: gradient_x(h_0) = 0
    c_1's constraint: gradient_y(h_0) = 0

    Both must be jointly satisfied: h_0's IC must have zero gradients
    in both x and y directions (i.e. h_0 must be constant or zero).
    """
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
        "fields": [
            {"name": "c_0", "index": 0},
            {"name": "c_1", "index": 1},
            {"name": "h_0", "index": 2},
        ],
        "equations": [
            {
                "field": "c_0",
                "lhs": {"expression": "c_0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "h_0",
                        },
                    ],
                },
            },
            {
                "field": "c_1",
                "lhs": {"expression": "c_1", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "h_0",
                        },
                    ],
                },
            },
            {
                "field": "h_0",
                "lhs": {
                    "expression": "d2_t(h_0)",
                    "order": {"time": 2},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian",
                            "field": "h_0",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
            "field_rates": {
                "h_0": [
                    {
                        "coefficient": 1.0,
                        "operator": "identity",
                        "field": "pi_2",
                    },
                ],
            },
            "kinetic_matrix": {
                "entries": [{"i": 0, "j": 0, "value": 1.0}],
                "dimension": 1,
            },
            "spatial_momenta": {},
            "hamiltonian_symbolic": "test",
        },
    }
    return EquationSystem.from_dict(data)


class TestNoSelfTermConstraints:
    """Test IDA handling of constraints with no self-referencing terms."""

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_no_self_term_freezes_field(self) -> None:
        """IDA freezes no-self-term constraint field at zero."""
        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Give h_1 a Gaussian initial condition
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.coord_arrays()[0]
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.warns(UserWarning, match="no self-referencing"):
            result = solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
            )

        assert result["success"], result["message"]

        # h_0 should be frozen at zero throughout
        h0_slot = layout.field_slot_map["h_0"]
        for yi in result["y"]:
            h0_data = yi[h0_slot * n : (h0_slot + 1) * n]
            assert np.allclose(h0_data, 0.0, atol=1e-12), (
                f"h_0 should be frozen at zero, got max={np.abs(h0_data).max()}"
            )

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_no_self_term_warns(self) -> None:
        """IDA emits UserWarning for no-self-term constraints."""
        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        with pytest.warns(UserWarning, match="Freezing 'h_0' at zero"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.1),
                bc="periodic",
                num_snapshots=3,
            )

    def test_sparsity_no_self_term_diagonal(self) -> None:
        """Sparsity pattern has diagonal block for no-self-term constraints."""
        from tidal.solver.sparsity import build_jacobian_sparsity

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")

        # h_0 is slot 0 — should have diagonal entries (identity coupling)
        h0_slot = layout.field_slot_map["h_0"]
        h0_start = h0_slot * n
        h0_end = (h0_slot + 1) * n

        # Extract h_0's row block
        h0_block: Any = pattern[h0_start:h0_end, h0_start:h0_end]  # pyright: ignore[reportUnknownVariableType]

        # Should have diagonal entries (at minimum)
        for i in range(n):
            assert h0_block[i, i] != 0, (
                f"Missing diagonal entry at ({i},{i}) for no-self-term field"
            )

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_ic_consistency_warns_when_violated(self) -> None:
        """IDA emits UserWarning when IC violates a no-self-term constraint."""
        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_1 to a Gaussian — h_0's equation (gradient_x of h_1)
        # will have nonzero residual → IC consistency warning
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.coord_arrays()[0]
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.warns(UserWarning, match="does not satisfy subsidiary"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
            )

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_ic_consistency_no_warn_when_satisfied(self) -> None:
        """No IC consistency warning when initial data satisfies constraint."""
        import warnings as _warnings

        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)
        # All zeros → constraint RHS = gradient_x(h_1) = 0 → satisfied

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.1),
                bc="periodic",
                num_snapshots=3,
            )

        ic_warns = [
            w for w in caught if "does not satisfy subsidiary" in str(w.message)
        ]
        assert len(ic_warns) == 0, "Should not warn when IC satisfies constraint"

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_multi_constraint_joint_ic_check(self) -> None:
        """Multiple no-self-term constraints: single warning lists all violations.

        Two constraint fields (c_0, c_1) each impose gradient conditions
        on h_0.  A Gaussian IC on h_0 violates both — the warning should
        mention both c_0 and c_1 in a single message.
        """
        from tidal.solver.ida import solve_ida

        spec = _make_multi_no_self_term_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_0 to a 2D Gaussian — violates both gradient_x and gradient_y
        h0_slot = layout.field_slot_map["h_0"]
        x_arr, y_arr = grid.coord_arrays()
        y0[h0_slot * n : (h0_slot + 1) * n] = np.exp(
            -((x_arr - 5) ** 2 + (y_arr - 5) ** 2) / 2
        ).ravel()

        import warnings as _warnings

        with _warnings.catch_warnings(record=True) as caught:
            _warnings.simplefilter("always")
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
            )

        # Should get exactly ONE summary warning mentioning both constraints
        ic_warns = [
            w for w in caught if "does not satisfy subsidiary" in str(w.message)
        ]
        assert len(ic_warns) == 1, f"Expected 1 summary warning, got {len(ic_warns)}"
        msg = str(ic_warns[0].message)
        assert "c_0" in msg, "Warning should mention c_0"
        assert "c_1" in msg, "Warning should mention c_1"
