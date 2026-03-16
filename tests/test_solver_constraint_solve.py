"""Tests for tidal.solver.constraint_solve — constraint pre-solve module."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.constraint_solve import (
    _classify_terms,
    _ConstraintTerms,
    _fft_solve_single,
    _find_connected_components,
    _matrix_solve,
    _probe_operator_matrix,
    _select_method,
    ensure_consistent_ic,
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
                            "field": "v_A_1",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_y", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_y",
                            "field": "v_A_2",
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
        "canonical": {
            "hamiltonian_terms": [],
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
                            "field": "v_A_1",
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
                            "field": "v_A_2",
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
                            "coefficient": -0.25,
                            "operator": "identity",
                            "field": "A_1",
                            "coefficient_symbolic": "-1/4*kappa^2",
                        },
                        {"coefficient": 1.0, "operator": "identity", "field": "A_1"},
                        {
                            "coefficient": -0.5,
                            "operator": "gradient_y",
                            "field": "A_0",
                            "coefficient_symbolic": "-1/2*kappa",
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
                            "coefficient": -0.25,
                            "operator": "identity",
                            "field": "A_2",
                            "coefficient_symbolic": "-1/4*kappa^2",
                        },
                        {
                            "coefficient": -0.5,
                            "operator": "identity",
                            "field": "v_A_2",
                            "coefficient_symbolic": "-1/2*kappa",
                        },
                        {
                            "coefficient": 0.5,
                            "operator": "gradient_x",
                            "field": "A_0",
                            "coefficient_symbolic": "kappa/2",
                        },
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
        "canonical": {
            "hamiltonian_terms": [],
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

        # Source terms: gradient_x(pi_A_1) + gradient_y(pi_A_2)
        assert len(terms.source_terms) == 2
        for _coeff, op, field in terms.source_terms:
            assert op in {"gradient_x", "gradient_y"}
            assert field in {"v_A_1", "v_A_2"}

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
# Unit tests: _find_connected_components
# ---------------------------------------------------------------------------


class TestFindConnectedComponents:
    """Unit tests for connected component decomposition of constraint groups."""

    @staticmethod
    def _make_group(
        name: str,
        source_fields: list[str] | None = None,
    ) -> _ConstraintTerms:
        """Create a minimal _ConstraintTerms for component testing."""
        config = ConstraintSolverConfig(enabled=True, method="auto")
        source_terms = [(1.0, "identity", f) for f in (source_fields or [])]
        return _ConstraintTerms(
            field_name=name,
            self_terms=[(1.0, "identity")],
            source_terms=source_terms,
            eq_idx=0,
            has_position_dependent_self=False,
            config=config,
        )

    def test_single_independent(self) -> None:
        """One constraint with no cross-refs → 1 component of size 1."""
        g = self._make_group("A", source_fields=["x", "y"])
        components = _find_connected_components([g])
        assert len(components) == 1
        assert len(components[0]) == 1
        assert components[0][0].field_name == "A"

    def test_all_independent(self) -> None:
        """Three constraints with no mutual refs → 3 components."""
        groups = [
            self._make_group("A", source_fields=["x"]),
            self._make_group("B", source_fields=["y"]),
            self._make_group("C", source_fields=["z"]),
        ]
        components = _find_connected_components(groups)
        assert len(components) == 3

    def test_all_coupled(self) -> None:
        """Three constraints, A→B, B→C → 1 component of size 3."""
        groups = [
            self._make_group("A", source_fields=["B"]),
            self._make_group("B", source_fields=["C"]),
            self._make_group("C", source_fields=["x"]),
        ]
        components = _find_connected_components(groups)
        assert len(components) == 1
        assert len(components[0]) == 3

    def test_two_components(self) -> None:
        """A↔B coupled, C independent → 2 components."""
        groups = [
            self._make_group("A", source_fields=["B", "x"]),
            self._make_group("B", source_fields=["A", "y"]),
            self._make_group("C", source_fields=["z"]),
        ]
        components = _find_connected_components(groups)
        assert len(components) == 2
        names = [{g.field_name for g in c} for c in components]
        assert {"A", "B"} in names
        assert {"C"} in names

    def test_massive_gravity_pattern(self) -> None:
        """Mimics massive_gravity: h_0 independent, h_1↔h_2 coupled.

        h_0 references only non-constraint fields (h_3, h_5).
        h_1 references h_2 (constraint). h_2 references h_1 (constraint).
        """
        groups = [
            self._make_group("h_0", source_fields=["h_3", "h_5"]),
            self._make_group("h_1", source_fields=["h_2", "v_h_0"]),
            self._make_group("h_2", source_fields=["h_1", "v_h_0"]),
        ]
        components = _find_connected_components(groups)
        assert len(components) == 2
        names = [{g.field_name for g in c} for c in components]
        assert {"h_0"} in names
        assert {"h_1", "h_2"} in names

    def test_name_map_resolves(self) -> None:
        """Cross-constraint ref via name_map is detected."""
        groups = [
            self._make_group("A", source_fields=["alias_B"]),
            self._make_group("B", source_fields=["x"]),
        ]
        name_map = {"alias_B": "B"}
        components = _find_connected_components(groups, name_map)
        assert len(components) == 1  # A→B via alias

    def test_self_loop(self) -> None:
        """Self-referencing constraint does not create phantom edges."""
        g = self._make_group("A", source_fields=["A", "x"])
        components = _find_connected_components([g])
        assert len(components) == 1
        assert len(components[0]) == 1

    def test_empty_input(self) -> None:
        """No groups → empty result."""
        assert _find_connected_components([]) == []

    def test_shared_non_constraint_field(self) -> None:
        """Non-constraint field shared by multiple constraints → no coupling."""
        groups = [
            self._make_group("A", source_fields=["x", "y"]),
            self._make_group("B", source_fields=["x", "z"]),
        ]
        components = _find_connected_components(groups)
        assert len(components) == 2


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

        # All zeros: pi_A_1 = pi_A_2 = 0 → source = 0 → A_0 = 0
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
        # Set pi_A_1 = cos(x)*cos(y) (zero-mean for compatibility)
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

        # allow_inconsistent_ic=True so constraint violation is a warning, not error
        import warnings as _warnings

        with _warnings.catch_warnings(record=True):
            _warnings.simplefilter("always")
            result = solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
                allow_inconsistent_ic=True,
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
    def test_ida_ic_consistency_errors_when_violated(self) -> None:
        """IDA raises ValueError when IC violates a no-self-term constraint (strict mode)."""
        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_1 to a Gaussian — h_0's equation (gradient_x of h_1)
        # will have nonzero residual → IC consistency error
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.coord_arrays()[0]
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.raises(ValueError, match="does not satisfy constraint"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
            )

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_ic_consistency_warns_when_lenient(self) -> None:
        """IDA emits UserWarning with allow_inconsistent_ic=True."""
        from tidal.solver.ida import solve_ida

        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.coord_arrays()[0]
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.warns(UserWarning, match="does not satisfy constraint"):
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
                allow_inconsistent_ic=True,
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
            w for w in caught if "does not satisfy constraint" in str(w.message)
        ]
        assert len(ic_warns) == 0, "Should not warn when IC satisfies constraint"

    @pytest.mark.skipif(not _has_sundials(), reason="sksundae not available")
    def test_ida_multi_constraint_joint_ic_check(self) -> None:
        """Multiple no-self-term constraints: ValueError lists all violations.

        Two constraint fields (c_0, c_1) each impose gradient conditions
        on h_0.  A Gaussian IC on h_0 violates both — the error should
        mention both c_0 and c_1 in the message.
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

        with pytest.raises(ValueError, match="does not satisfy constraint") as exc_info:
            solve_ida(
                spec,
                grid,
                y0,
                (0.0, 0.5),
                bc="periodic",
                num_snapshots=3,
            )

        msg = str(exc_info.value)
        assert "c_0" in msg, "Error should mention c_0"
        assert "c_1" in msg, "Error should mention c_1"


# ---------------------------------------------------------------------------
# Helper: GW-like spec for subsidiary constraint cascade tests
# ---------------------------------------------------------------------------


def _make_gw_like_spec() -> EquationSystem:
    """Simplified GW-like spec with subsidiary constraint cascade.

    Models the TT gauge constraint structure:
    - h_a: dynamical wave equation (d2_t = laplacian)
    - h_b: subsidiary constraint: 0 = -gradient_x(h_a) - gradient_x(h_b)
      → when h_a = Gaussian, h_b = -Gaussian
    - h_c: traceless constraint: 0 = identity(h_a) + identity(h_b) + identity(h_c)
      → after h_b determined, h_c = -(h_a + h_b) = 0
    """
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [
            {"name": "h_a", "index": 0},
            {"name": "h_b", "index": 1},
            {"name": "h_c", "index": 2},
        ],
        "equations": [
            # h_a: dynamical wave equation
            {
                "field": "h_a",
                "lhs": {
                    "expression": "d2_t(h_a)",
                    "order": {"time": 2},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian",
                            "field": "h_a",
                        },
                    ],
                },
            },
            # h_b: subsidiary transverse constraint
            # 0 = -gradient_x(h_a) - gradient_x(h_b)
            {
                "field": "h_b",
                "lhs": {"expression": "h_b", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "h_a",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "h_b",
                        },
                    ],
                },
            },
            # h_c: traceless constraint
            # 0 = identity(h_a) + identity(h_b) + identity(h_c)
            {
                "field": "h_c",
                "lhs": {"expression": "h_c", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "h_a",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "h_b",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "h_c",
                        },
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
        },
    }
    return EquationSystem.from_dict(data)


def _make_velocity_constraint_spec() -> EquationSystem:
    """Spec with a constraint referencing only velocities.

    When velocities are zero (default), constraint is trivially satisfied.
    """
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [
            {"name": "h_0", "index": 0},
            {"name": "h_1", "index": 1},
        ],
        "equations": [
            # h_0: velocity constraint (0 = gradient_x(v_h_1))
            {
                "field": "h_0",
                "lhs": {"expression": "h_0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -0.5,
                            "operator": "gradient_x",
                            "field": "v_h_1",
                        },
                    ],
                },
            },
            # h_1: dynamical wave equation
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
        },
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# Tests for ensure_consistent_ic — unified constraint IC solver
# ---------------------------------------------------------------------------


class TestEnsureConsistentIC:
    """Tests for the unified constraint IC solver."""

    def test_no_constraints_passthrough(self) -> None:
        """No constraint equations → y0 returned unchanged."""
        spec = _make_kg_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        n = grid.num_points
        gaussian = np.exp(-((x - 5) ** 2) / 2)
        slot = layout.field_slot_map["phi_0"]
        y0[slot * n : (slot + 1) * n] = gaussian.ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")
        np.testing.assert_array_equal(result, y0)

    def test_standard_constraint_em(self) -> None:
        """EM-like standard constraint (enabled=True) → same as pre_solve."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        y0 = np.zeros(layout.total_size)

        # Set v_A_1 to a Gaussian (source for A_0 constraint)
        va1_slot = layout.velocity_slot_map["A_1"]
        x_arr, _ = grid.coord_arrays()
        y0[va1_slot * n : (va1_slot + 1) * n] = np.exp(-((x_arr - 5) ** 2) / 2).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # A_0 should be non-zero (solved from source)
        a0_slot = layout.field_slot_map["A_0"]
        a0_data = result[a0_slot * n : (a0_slot + 1) * n].reshape(grid.shape)
        assert float(np.max(np.abs(a0_data))) > 1e-5, (
            "A_0 should be non-zero after constraint solve"
        )

    def test_subsidiary_gradient_solve(self) -> None:
        """Subsidiary constraint: gradient_x(h_a) + gradient_x(h_b) = 0.

        When h_a = Gaussian, h_b should be determined as -Gaussian.
        """
        spec = _make_gw_like_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_a to Gaussian
        ha_slot = layout.field_slot_map["h_a"]
        x = grid.axes_coords(0)
        gaussian = np.exp(-((x - 5) ** 2) / 2)
        y0[ha_slot * n : (ha_slot + 1) * n] = gaussian.ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # h_b should be approximately -h_a (up to FFT zero-mode)
        hb_slot = layout.field_slot_map["h_b"]
        hb_data = result[hb_slot * n : (hb_slot + 1) * n]
        ha_data = result[ha_slot * n : (ha_slot + 1) * n]

        # Check that gradient_x(h_a) + gradient_x(h_b) ≈ 0
        from tidal.solver.operators import gradient

        grad_ha = gradient(ha_data.reshape(grid.shape), 0, grid, "periodic")
        grad_hb = gradient(hb_data.reshape(grid.shape), 0, grid, "periodic")
        residual = grad_ha + grad_hb
        assert float(np.max(np.abs(residual))) < 1e-10, (
            f"Subsidiary constraint residual too large: {np.max(np.abs(residual)):.2e}"
        )

    def test_cascade_two_step(self) -> None:
        """Cascade: h_b determined from h_a, then h_c from h_a + h_b.

        h_b's transverse constraint: gradient_x(h_a) + gradient_x(h_b) = 0
        h_c's traceless constraint: identity(h_a) + identity(h_b) + identity(h_c) = 0

        First iteration determines h_b, second determines h_c.
        """
        spec = _make_gw_like_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_a to Gaussian
        ha_slot = layout.field_slot_map["h_a"]
        x = grid.axes_coords(0)
        gaussian = np.exp(-((x - 5) ** 2) / 2)
        y0[ha_slot * n : (ha_slot + 1) * n] = gaussian.ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # After cascade, h_c = -(h_a + h_b)
        hc_slot = layout.field_slot_map["h_c"]
        ha_data = result[ha_slot * n : (ha_slot + 1) * n]
        hb_slot = layout.field_slot_map["h_b"]
        hb_data = result[hb_slot * n : (hb_slot + 1) * n]
        hc_data = result[hc_slot * n : (hc_slot + 1) * n]

        # h_c should satisfy: h_a + h_b + h_c = 0
        trace_residual = ha_data + hb_data + hc_data
        assert float(np.max(np.abs(trace_residual))) < 1e-10, (
            f"Traceless constraint residual: {np.max(np.abs(trace_residual)):.2e}"
        )

    def test_all_satisfied_noop(self) -> None:
        """All constraints already satisfied → y0 unchanged."""
        spec = _make_gw_like_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)
        # All zeros → all constraints trivially satisfied

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")
        np.testing.assert_array_equal(result, y0)

    def test_violated_strict_raises(self) -> None:
        """Constraint violated with no free fields → ValueError in strict mode."""
        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # h_1 is determined (Gaussian) but constraint 0 = gradient_x(h_1) violated
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.axes_coords(0)
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.raises(ValueError, match="does not satisfy constraint"):
            ensure_consistent_ic(spec, grid, y0, bc="periodic", strict=True)

    def test_violated_lenient_warns(self) -> None:
        """Constraint violated with strict=False → UserWarning, no error."""
        spec = _make_no_self_term_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.axes_coords(0)
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        with pytest.warns(UserWarning, match="does not satisfy constraint"):
            result = ensure_consistent_ic(
                spec,
                grid,
                y0,
                bc="periodic",
                strict=False,
            )

        # Should return y0 unchanged (no field could be solved)
        np.testing.assert_array_equal(result, y0)

    def test_velocity_refs_auto_satisfied(self) -> None:
        """Velocity-referencing constraint: auto-satisfied when velocities are zero."""
        spec = _make_velocity_constraint_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_1 to Gaussian (only field IC, not velocity)
        h1_slot = layout.field_slot_map["h_1"]
        x = grid.axes_coords(0)
        y0[h1_slot * n : (h1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        # Should not raise — velocity constraint is trivially satisfied
        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")
        # h_1 should be unchanged
        h1_out = result[h1_slot * n : (h1_slot + 1) * n]
        h1_in = y0[h1_slot * n : (h1_slot + 1) * n]
        np.testing.assert_array_equal(h1_out, h1_in)

    def test_em_backward_compatible(self) -> None:
        """EM spec: ensure_consistent_ic produces same result as pre_solve_constraints."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        y0 = np.zeros(layout.total_size)

        # Set v_A_1 = Gaussian as source
        va1_slot = layout.velocity_slot_map["A_1"]
        x_arr, _ = grid.coord_arrays()
        source = np.exp(-((x_arr - 5) ** 2) / 2)
        y0[va1_slot * n : (va1_slot + 1) * n] = source.ravel()

        # pre_solve_constraints (existing function)
        y0_old = pre_solve_constraints(spec, grid, y0.copy(), bc="periodic")

        # ensure_consistent_ic (new unified function)
        y0_new = ensure_consistent_ic(spec, grid, y0.copy(), bc="periodic")

        # Both should produce the same A_0 solution
        a0_slot = layout.field_slot_map["A_0"]
        a0_old = y0_old[a0_slot * n : (a0_slot + 1) * n]
        a0_new = y0_new[a0_slot * n : (a0_slot + 1) * n]
        np.testing.assert_allclose(a0_new, a0_old, atol=1e-12)

    def test_cs_backward_compatible(self) -> None:
        """Chern-Simons: ensure_consistent_ic matches pre_solve_constraints."""
        spec = _make_chern_simons_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        y0 = np.zeros(layout.total_size)

        # Set v_A_1 and A_2 as sources (CS has both velocity and field sources)
        va1_slot = layout.velocity_slot_map["A_1"]
        x_arr, _ = grid.coord_arrays()
        y0[va1_slot * n : (va1_slot + 1) * n] = np.exp(-((x_arr - 5) ** 2) / 2).ravel()

        y0_old = pre_solve_constraints(
            spec,
            grid,
            y0.copy(),
            bc="periodic",
            parameters={"kappa": 0.5},
        )
        y0_new = ensure_consistent_ic(
            spec,
            grid,
            y0.copy(),
            bc="periodic",
            parameters={"kappa": 0.5},
        )

        a0_slot = layout.field_slot_map["A_0"]
        a0_old = y0_old[a0_slot * n : (a0_slot + 1) * n]
        a0_new = y0_new[a0_slot * n : (a0_slot + 1) * n]
        np.testing.assert_allclose(a0_new, a0_old, atol=1e-12)

    def test_helmholtz_constraint(self) -> None:
        """Helmholtz (laplacian - m²): non-singular constraint solve.

        Unlike Poisson (singular at k=0), Helmholtz has no zero-mode issue.
        Ensures the non-singular path works in ensure_consistent_ic.
        """
        spec = _make_helmholtz_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        # Set rho_0 = cos(x)*cos(y) as source
        rho_slot = layout.field_slot_map["rho_0"]
        y0[rho_slot * n : (rho_slot + 1) * n] = (np.cos(x_arr) * np.cos(y_arr)).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # phi_0 should be non-zero (solved from rho_0)
        phi_slot = layout.field_slot_map["phi_0"]
        phi_data = result[phi_slot * n : (phi_slot + 1) * n].reshape(grid.shape)
        assert float(np.max(np.abs(phi_data))) > 1e-5

        # Verify discrete residual: (laplacian - 1)phi_0 + rho_0 ≈ 0
        from tidal.solver.operators import laplacian as lap_op

        lap_phi = lap_op(phi_data, grid, "periodic")
        rho_data = y0[rho_slot * n : (rho_slot + 1) * n].reshape(grid.shape)
        residual = lap_phi - phi_data + rho_data
        assert float(np.max(np.abs(residual))) < 1e-10, (
            f"Helmholtz residual: {np.max(np.abs(residual)):.2e}"
        )

    def test_neumann_bc_matrix_path(self) -> None:
        """Non-periodic BCs force the matrix solver path.

        Helmholtz with Neumann BCs — cannot use FFT, must use matrix probe.
        """
        spec = _make_helmholtz_spec()
        grid = GridInfo(
            bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
            shape=(16, 16),
            periodic=(False, False),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points
        x_arr, y_arr = grid.coord_arrays()

        y0 = np.zeros(layout.total_size)
        rho_slot = layout.field_slot_map["rho_0"]
        # Source compatible with Neumann BCs (zero normal derivative at boundary)
        y0[rho_slot * n : (rho_slot + 1) * n] = (np.cos(x_arr) * np.cos(y_arr)).ravel()

        bc = ("neumann", "neumann")
        result = ensure_consistent_ic(spec, grid, y0, bc=bc)

        # phi_0 should be solved (non-zero)
        phi_slot = layout.field_slot_map["phi_0"]
        phi_data = result[phi_slot * n : (phi_slot + 1) * n].reshape(grid.shape)
        assert float(np.max(np.abs(phi_data))) > 1e-5

    def test_identity_only_self_term(self) -> None:
        """Identity-only self-term: direct algebraic solve.

        Constraint: 0 = identity(h_a) + identity(h_b) + identity(h_c)
        When h_a and h_b are known, h_c = -(h_a + h_b) — solved directly
        via identity operator (trivial inversion).
        """
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "h_a", "index": 0},
                {"name": "h_b", "index": 1},
                {"name": "h_c", "index": 2},
            ],
            "equations": [
                {
                    "field": "h_a",
                    "lhs": {"expression": "d2_t(h_a)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_a",
                            },
                        ],
                    },
                },
                {
                    "field": "h_b",
                    "lhs": {"expression": "d2_t(h_b)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_b",
                            },
                        ],
                    },
                },
                {
                    "field": "h_c",
                    "lhs": {"expression": "h_c", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "h_a",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "h_b",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "h_c",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        gaussian = np.exp(-((x - 5) ** 2) / 2)

        # Set h_a and h_b to different Gaussians
        ha_slot = layout.field_slot_map["h_a"]
        hb_slot = layout.field_slot_map["h_b"]
        y0[ha_slot * n : (ha_slot + 1) * n] = gaussian.ravel()
        y0[hb_slot * n : (hb_slot + 1) * n] = (0.5 * gaussian).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # h_c should be -(h_a + h_b)
        hc_slot = layout.field_slot_map["h_c"]
        ha_data = result[ha_slot * n : (ha_slot + 1) * n]
        hb_data = result[hb_slot * n : (hb_slot + 1) * n]
        hc_data = result[hc_slot * n : (hc_slot + 1) * n]
        expected = -(ha_data + hb_data)
        np.testing.assert_allclose(hc_data, expected, atol=1e-12)

    def test_disabled_self_term_constraint(self) -> None:
        """Disabled constraint with self-terms: Phase 2 still solves it.

        Even when constraint_solver.enabled=False, Phase 2 of
        ensure_consistent_ic can solve for the constraint field if it
        has self-terms and is the only free field.
        """
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "psi_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "phi_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "psi_0",
                            },
                        ],
                    },
                    "constraint_solver": {"enabled": False},
                },
                {
                    "field": "psi_0",
                    "lhs": {"expression": "d2_t(psi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "psi_0",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        psi_slot = layout.field_slot_map["psi_0"]
        y0[psi_slot * n : (psi_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # phi_0 should be -psi_0 (from identity(phi_0) + identity(psi_0) = 0)
        phi_slot = layout.field_slot_map["phi_0"]
        phi_data = result[phi_slot * n : (phi_slot + 1) * n]
        psi_data = result[psi_slot * n : (psi_slot + 1) * n]
        np.testing.assert_allclose(phi_data, -psi_data, atol=1e-12)

    def test_underdetermined_graceful_skip(self) -> None:
        """Underdetermined constraint (2+ free fields): skipped gracefully.

        When a constraint references two free (zero) fields and there's no
        way to determine either one, the solver should skip it without error.
        """
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "c_0", "index": 0},
                {"name": "h_a", "index": 1},
                {"name": "h_b", "index": 2},
            ],
            "equations": [
                # Both h_a and h_b are free → underdetermined
                {
                    "field": "c_0",
                    "lhs": {"expression": "c_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "h_a",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "h_b",
                            },
                        ],
                    },
                },
                {
                    "field": "h_a",
                    "lhs": {"expression": "d2_t(h_a)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_a",
                            },
                        ],
                    },
                },
                {
                    "field": "h_b",
                    "lhs": {"expression": "d2_t(h_b)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_b",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)

        # Both h_a and h_b are zero → constraint is trivially satisfied
        # but underdetermined (can't solve for either individually)
        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")
        np.testing.assert_array_equal(result, y0)

    def test_phase2_stall_emits_warning(self) -> None:
        """Phase 2 stall emits actionable warning about undetermined fields."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "c_0", "index": 0},
                {"name": "h_a", "index": 1},
                {"name": "h_b", "index": 2},
            ],
            "equations": [
                {
                    "field": "c_0",
                    "lhs": {"expression": "c_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "h_a",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "h_b",
                            },
                        ],
                    },
                },
                {
                    "field": "h_a",
                    "lhs": {"expression": "d2_t(h_a)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_a",
                            },
                        ],
                    },
                },
                {
                    "field": "h_b",
                    "lhs": {"expression": "d2_t(h_b)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_b",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        y0 = np.zeros(StateLayout.from_spec(spec, grid.num_points).total_size)

        with pytest.warns(UserWarning, match="Constraint propagation stalled"):
            ensure_consistent_ic(spec, grid, y0, bc="periodic")

    def test_multi_independent_subsidiary(self) -> None:
        """Two independent subsidiary constraints, each solving for a different field.

        Constraint 1: gradient_x(h_a) + gradient_x(h_b) = 0  → solve h_b
        Constraint 2: gradient_x(h_a) + gradient_x(h_c) = 0  → solve h_c

        Both determined independently from h_a in a single iteration.
        """
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "h_a", "index": 0},
                {"name": "h_b", "index": 1},
                {"name": "h_c", "index": 2},
                {"name": "c_0", "index": 3},
                {"name": "c_1", "index": 4},
            ],
            "equations": [
                {
                    "field": "h_a",
                    "lhs": {"expression": "d2_t(h_a)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_a",
                            },
                        ],
                    },
                },
                {
                    "field": "h_b",
                    "lhs": {"expression": "d2_t(h_b)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_b",
                            },
                        ],
                    },
                },
                {
                    "field": "h_c",
                    "lhs": {"expression": "d2_t(h_c)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "h_c",
                            },
                        ],
                    },
                },
                # Subsidiary constraint: gradient_x(h_a) + gradient_x(h_b) = 0
                {
                    "field": "c_0",
                    "lhs": {"expression": "c_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "h_a",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "h_b",
                            },
                        ],
                    },
                },
                # Subsidiary constraint: gradient_x(h_a) + gradient_x(h_c) = 0
                {
                    "field": "c_1",
                    "lhs": {"expression": "c_1", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "h_a",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "h_c",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        gaussian = np.exp(-((x - 5) ** 2) / 2)
        ha_slot = layout.field_slot_map["h_a"]
        y0[ha_slot * n : (ha_slot + 1) * n] = gaussian.ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # Both h_b and h_c should satisfy gradient constraints
        from tidal.solver.operators import gradient

        ha_data = result[ha_slot * n : (ha_slot + 1) * n].reshape(grid.shape)
        hb_slot = layout.field_slot_map["h_b"]
        hb_data = result[hb_slot * n : (hb_slot + 1) * n].reshape(grid.shape)
        hc_slot = layout.field_slot_map["h_c"]
        hc_data = result[hc_slot * n : (hc_slot + 1) * n].reshape(grid.shape)

        grad_ha = gradient(ha_data, 0, grid, "periodic")
        grad_hb = gradient(hb_data, 0, grid, "periodic")
        grad_hc = gradient(hc_data, 0, grid, "periodic")

        # gradient(h_a) + gradient(h_b) ≈ 0
        res_b = grad_ha + grad_hb
        assert float(np.max(np.abs(res_b))) < 1e-10, (
            f"c_0 residual: {np.max(np.abs(res_b)):.2e}"
        )
        # gradient(h_a) + gradient(h_c) ≈ 0
        res_c = grad_ha + grad_hc
        assert float(np.max(np.abs(res_c))) < 1e-10, (
            f"c_1 residual: {np.max(np.abs(res_c)):.2e}"
        )

    def test_coupled_constraint_enabled(self) -> None:
        """Two coupled constraints (enabled=True) with mutual references.

        Models coupled Proca: A_0 and B_0 each reference each other.
        """
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "A_0", "index": 0},
                {"name": "B_0", "index": 1},
                {"name": "A_1", "index": 2},
                {"name": "B_1", "index": 3},
            ],
            "equations": [
                # A_0 constraint: laplacian(A_0) - A_0 + 0.5*B_0 = source_A
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
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "A_0",
                            },
                            {
                                "coefficient": 0.5,
                                "operator": "identity",
                                "field": "B_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "v_A_1",
                            },
                        ],
                    },
                    "constraint_solver": {
                        "enabled": True,
                        "method": "auto",
                        "boundary_conditions": {"x": {"type": "periodic"}},
                    },
                },
                # B_0 constraint: laplacian(B_0) - B_0 + 0.5*A_0 = source_B
                {
                    "field": "B_0",
                    "lhs": {"expression": "B_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "B_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "B_0",
                            },
                            {
                                "coefficient": 0.5,
                                "operator": "identity",
                                "field": "A_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "gradient_x",
                                "field": "v_B_1",
                            },
                        ],
                    },
                    "constraint_solver": {
                        "enabled": True,
                        "method": "auto",
                        "boundary_conditions": {"x": {"type": "periodic"}},
                    },
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
                {
                    "field": "B_1",
                    "lhs": {"expression": "d2_t(B_1)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "B_1",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x = grid.axes_coords(0)
        # Set v_A_1 source (for A_0)
        va1_slot = layout.velocity_slot_map["A_1"]
        y0[va1_slot * n : (va1_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # Both A_0 and B_0 should be solved (non-zero due to coupling)
        a0_slot = layout.field_slot_map["A_0"]
        a0_data = result[a0_slot * n : (a0_slot + 1) * n]
        assert float(np.max(np.abs(a0_data))) > 1e-6, "A_0 should be non-zero"

        b0_slot = layout.field_slot_map["B_0"]
        b0_data = result[b0_slot * n : (b0_slot + 1) * n]
        # B_0 is coupled to A_0 — should also be non-zero
        assert float(np.max(np.abs(b0_data))) > 1e-6, (
            "B_0 should be non-zero from coupling to A_0"
        )

    def test_does_not_mutate_input(self) -> None:
        """ensure_consistent_ic does not modify the input y0 array."""
        spec = _make_gw_like_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        ha_slot = layout.field_slot_map["h_a"]
        x = grid.axes_coords(0)
        y0[ha_slot * n : (ha_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        y0_orig = y0.copy()
        _ = ensure_consistent_ic(spec, grid, y0, bc="periodic")
        np.testing.assert_array_equal(y0, y0_orig)

    def test_multiple_violations_all_reported(self) -> None:
        """When multiple constraints are violated, all are listed in the error."""
        spec = _make_multi_no_self_term_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        h0_slot = layout.field_slot_map["h_0"]
        x_arr, y_arr = grid.coord_arrays()
        y0[h0_slot * n : (h0_slot + 1) * n] = np.exp(
            -((x_arr - 5) ** 2 + (y_arr - 5) ** 2) / 2
        ).ravel()

        with pytest.raises(ValueError, match="does not satisfy constraint") as exc_info:
            ensure_consistent_ic(spec, grid, y0, bc="periodic", strict=True)

        # Both c_0 and c_1 should be mentioned
        msg = str(exc_info.value)
        assert "c_0" in msg, "Error should mention c_0"
        assert "c_1" in msg, "Error should mention c_1"


# ---------------------------------------------------------------------------
# Integration test: real GW plane-wave JSON
# ---------------------------------------------------------------------------


class TestEnsureConsistentICRealJSON:
    """Integration tests using real derived JSON specifications."""

    def test_gw_plane_wave_cascade(self) -> None:
        """GW 1D (2 fields after full constraint elimination): set h_7 and verify.

        After derive-time constraint elimination:
          - Traceless: h_9 eliminated, h_4 = -h_7 (degenerate pair)
          - Transverse: h_6, h_8 eliminated (gradient_z = 0 → field = 0)
        Only h_5 (cross) and h_7 (plus) survive — both free wave equations.
        """
        import json
        from pathlib import Path

        json_path = Path("examples/data/gw_plane_wave_1d.json")
        if not json_path.exists():
            pytest.skip("GW plane-wave JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        # Should have exactly 2 fields after full elimination
        field_names = {eq.field_name for eq in spec.equations}
        assert field_names == {"h_5", "h_7"}, (
            f"Expected 2-field system (h_5, h_7) but got {field_names}"
        )

        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set h_7 (plus polarization) to a Gaussian
        h7_slot = layout.field_slot_map["h_7"]
        x = grid.axes_coords(0)
        gaussian = np.exp(-((x - 5) ** 2) / 2)
        y0[h7_slot * n : (h7_slot + 1) * n] = gaussian.ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # h_7 should retain the Gaussian IC (no constraints to modify it)
        h7_data = result[h7_slot * n : (h7_slot + 1) * n]
        np.testing.assert_allclose(h7_data, gaussian.ravel(), atol=1e-12)

    def test_gw_plane_wave_h5_polarization(self) -> None:
        """GW: setting h_5 (cross polarization) instead of h_7.

        h_5 is the other physical DOF (cross polarization). Setting h_5
        should not affect h_7 — they are independent wave equations.
        """
        import json
        from pathlib import Path

        json_path = Path("examples/data/gw_plane_wave_1d.json")
        if not json_path.exists():
            pytest.skip("GW plane-wave JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        h5_slot = layout.field_slot_map["h_5"]
        x = grid.axes_coords(0)
        y0[h5_slot * n : (h5_slot + 1) * n] = np.exp(-((x - 5) ** 2) / 2).ravel()

        # Should not raise — no constraints in 2-field system
        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # h_7 should remain zero (independent wave equation)
        h7_slot = layout.field_slot_map["h_7"]
        h7_data = result[h7_slot * n : (h7_slot + 1) * n]
        assert np.allclose(h7_data, 0.0, atol=1e-12), (
            "h_7 should remain zero when only h_5 is set"
        )

    def test_em_json_constraint(self) -> None:
        """EM JSON: standard Gauss's law constraint solve from real spec."""
        import json
        from pathlib import Path

        json_path = Path("examples/data/em_3d.json")
        if not json_path.exists():
            pytest.skip("EM JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # Set v_A_1 as source
        va1_slot = layout.velocity_slot_map["A_1"]
        x_arr, _ = grid.coord_arrays()
        y0[va1_slot * n : (va1_slot + 1) * n] = np.exp(-((x_arr - 5) ** 2) / 2).ravel()

        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # A_0 should be solved (non-zero)
        a0_slot = layout.field_slot_map["A_0"]
        a0_data = result[a0_slot * n : (a0_slot + 1) * n]
        assert float(np.max(np.abs(a0_data))) > 1e-5

    def test_massive_gravity_connected_components(self) -> None:
        """Massive gravity: h_0 skipped (no self-terms), h_1↔h_2 coupled.

        With the correct linearized Einstein tensor (from √|g|·L
        perturbation), the tt equation (h_0) doesn't involve h_0 — it's
        a compatibility constraint on spatial components h_3/h_4/h_5.
        The solver skips h_0 (no self-terms) and solves h_1↔h_2 as a
        coupled screened-Poisson system.

        IC: traceless spatial perturbation h_3 = -h_5 = cos(k*x)*cos(k*y)
        with k_x = k_y satisfies the h_0 constraint exactly (the
        -m2*(h_3+h_5) term vanishes and the mixed Laplacians cancel).
        """
        import json
        from pathlib import Path

        json_path = Path("examples/data/massive_gravity_3d.json")
        if not json_path.exists():
            pytest.skip("Massive gravity JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        grid = GridInfo(
            bounds=((0, 10), (0, 10)),
            shape=(16, 16),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x_arr, y_arr = grid.coord_arrays()
        k = 2 * np.pi / 10  # k_x = k_y = 2π/L
        wave = 0.5 * np.cos(k * x_arr) * np.cos(k * y_arr)

        # Traceless spatial IC: h_3 (h_xx) = -h_5 (h_yy)
        h3_slot = layout.field_slot_map["h_3"]
        h5_slot = layout.field_slot_map["h_5"]
        y0[h3_slot * n : (h3_slot + 1) * n] = wave.ravel()
        y0[h5_slot * n : (h5_slot + 1) * n] = (-wave).ravel()

        result = ensure_consistent_ic(
            spec,
            grid,
            y0,
            bc="periodic",
            parameters={"m2": 1.0},
        )

        # h_0 remains at zero (not constrained by its own equation)
        h0_slot = layout.field_slot_map["h_0"]
        h0_data = result[h0_slot * n : (h0_slot + 1) * n]
        assert float(np.max(np.abs(h0_data))) < 1e-10, (
            "h_0 should remain zero (equation has no self-terms)"
        )

    def test_massive_gravity_inconsistent_ic(self) -> None:
        """Massive gravity: IC that violates the h_0 compatibility constraint.

        Setting only h_3 (with h_4 = h_5 = 0) violates the constraint
        0 = -m2*(h_3+h_5) + 0.5*∇²_x(h_5) + 0.5*∇²_y(h_3) - ∂²_xy(h_4).
        With strict=True this raises ValueError; with strict=False it
        warns but proceeds.
        """
        import json
        from pathlib import Path

        json_path = Path("examples/data/massive_gravity_3d.json")
        if not json_path.exists():
            pytest.skip("Massive gravity JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        grid = GridInfo(
            bounds=((0, 50), (0, 50)),
            shape=(32, 32),
            periodic=(True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        x_arr, y_arr = grid.coord_arrays()
        h3_slot = layout.field_slot_map["h_3"]
        y0[h3_slot * n : (h3_slot + 1) * n] = (
            np.exp(-((x_arr - 25) ** 2 + (y_arr - 25) ** 2) / 50)
            * np.sin(2 * np.pi * x_arr / 50)
        ).ravel()

        # strict=True should raise for the violated h_0 constraint
        with pytest.raises(ValueError, match="does not satisfy"):
            ensure_consistent_ic(
                spec,
                grid,
                y0,
                bc="periodic",
                parameters={"m2": 1.0},
                strict=True,
            )

        # strict=False should warn but succeed
        with pytest.warns(UserWarning, match="does not satisfy"):
            result = ensure_consistent_ic(
                spec,
                grid,
                y0,
                bc="periodic",
                parameters={"m2": 1.0},
                strict=False,
            )

        # h_0 should remain at zero (equation doesn't involve h_0)
        h0_slot = layout.field_slot_map["h_0"]
        h0_data = result[h0_slot * n : (h0_slot + 1) * n]
        assert float(np.max(np.abs(h0_data))) < 1e-10

    def test_gw_3d_with_complete_tt_ic(self) -> None:
        """GW 3+1D: complete TT gauge IC (h_4 + h_7 = -h_4) passes all constraints.

        Setting both h_4 (h_xx) and h_7 (h_yy) = -h_4 for z-propagating
        + polarization satisfies the traceless constraint and allows the
        cascade to determine h_9 = 0. Transverse constraints are trivially
        satisfied since h_4/h_7 depend only on z.
        """
        import json
        from pathlib import Path

        json_path = Path("examples/data/linearized_gravity.json")
        if not json_path.exists():
            pytest.skip("Linearized gravity JSON not available")

        with json_path.open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)

        grid = GridInfo(
            bounds=((0, 4), (0, 4), (0, 10)),
            shape=(4, 4, 16),
            periodic=(True, True, True),
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        n = grid.num_points

        y0 = np.zeros(layout.total_size)
        # z-dependent wave packet for h_4 (h_xx)
        _, _, z_arr = grid.coord_arrays()
        wave = np.exp(-((z_arr - 5) ** 2) / 2) * np.cos(0.63 * (z_arr - 5))

        h4_slot = layout.field_slot_map["h_4"]
        y0[h4_slot * n : (h4_slot + 1) * n] = wave.ravel()

        # TT gauge: h_7 (h_yy) = -h_4 (traceless for z-propagation)
        h7_slot = layout.field_slot_map["h_7"]
        y0[h7_slot * n : (h7_slot + 1) * n] = (-wave).ravel()

        # Should succeed: h_9 = -(h_4 + h_7) = 0, transverse trivially OK
        result = ensure_consistent_ic(spec, grid, y0, bc="periodic")

        # h_4 and h_7 should be preserved exactly (user-provided IC)
        h4_data = result[h4_slot * n : (h4_slot + 1) * n]
        np.testing.assert_allclose(h4_data, wave.ravel(), atol=1e-14)
        h7_data = result[h7_slot * n : (h7_slot + 1) * n]
        np.testing.assert_allclose(h7_data, (-wave).ravel(), atol=1e-14)

        # h_9 should be solved to ≈ 0 (traceless: h_4 + h_7 + h_9 = 0)
        h9_slot = layout.field_slot_map["h_9"]
        h9_data = result[h9_slot * n : (h9_slot + 1) * n]
        assert np.allclose(h9_data, 0.0, atol=1e-10), (
            f"h_9 should be ≈ 0 for TT gauge, got max={np.max(np.abs(h9_data)):.2e}"
        )
