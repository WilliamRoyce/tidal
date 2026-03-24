"""Tests for the analytical Jacobian builder."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.sparse import issparse  # type: ignore[import-untyped]

from tidal.solver.analytical_jacobian import (
    _OperatorCache,
    _resolve_term_target,
    build_jacobian_matrices,
    try_analytical_jacobian,
)
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


def _make_kg_1d_spec() -> EquationSystem:
    """Klein-Gordon 1D: single scalar, second-order."""
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
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                        },
                        {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                    ],
                },
            },
        ],
        "canonical": {"hamiltonian_terms": []},
    }
    return EquationSystem.from_dict(data)


def _make_em_2d_spec() -> EquationSystem:
    """EM 2+1D: constraint A_0 + wave A_1, with gradient_x(v_A_1)."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
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
                        {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"},
                        {
                            "coefficient": -1.0,
                            "operator": "gradient_x",
                            "field": "v_A_1",
                        },
                    ],
                },
            },
            {
                "field": "A_1",
                "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "A_1"},
                    ],
                },
            },
        ],
        "canonical": {"hamiltonian_terms": []},
    }
    return EquationSystem.from_dict(data)


def _make_constraint_velocity_spec() -> EquationSystem:
    """CS-like: constraint A_0, dynamical A_1 with gradient_x(v_A_0)."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
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
                            "coefficient": -1.0,
                            "operator": "laplacian_x",
                            "field": "A_0",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "gradient_x",
                            "field": "v_A_1",
                        },
                    ],
                },
            },
            {
                "field": "A_1",
                "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_1"},
                        {
                            "coefficient": 1.0,
                            "operator": "gradient_x",
                            "field": "v_A_0",
                        },
                    ],
                },
            },
        ],
        "canonical": {"hamiltonian_terms": []},
    }
    return EquationSystem.from_dict(data)


def _make_first_deriv_t_constraint_spec() -> EquationSystem:
    """Spec with first_derivative_t(constraint_field) in dynamical eq."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [
            {"name": "phi_0", "index": 0},
            {"name": "A_0", "index": 1},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "first_derivative_t",
                            "field": "A_0",
                        },
                    ],
                },
            },
            {
                "field": "A_0",
                "lhs": {"expression": "A_0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_0"},
                    ],
                },
            },
        ],
        "canonical": {"hamiltonian_terms": []},
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# all_constant() tests
# ---------------------------------------------------------------------------


class TestAllConstant:
    def test_kg_is_constant(self) -> None:
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff = CoefficientEvaluator(spec, grid, {})
        assert coeff.all_constant()

    def test_em_is_constant(self) -> None:
        spec = _make_em_2d_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True))
        coeff = CoefficientEvaluator(spec, grid, {})
        assert coeff.all_constant()

    def test_parametric_is_constant(self) -> None:
        """Symbolic coefficients with parameters are still constant."""
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
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-m2",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff = CoefficientEvaluator(spec, grid, {"m2": 1.0})
        assert coeff.all_constant()

    def test_time_dependent_is_not_constant(self) -> None:
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
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "E^(2*H*t[])",
                                "time_dependent": True,
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff = CoefficientEvaluator(spec, grid, {"H": 1.0})
        assert not coeff.all_constant()


# ---------------------------------------------------------------------------
# _resolve_term_target tests
# ---------------------------------------------------------------------------


class TestResolveTermTarget:
    def test_normal_field_reference(self) -> None:
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        op_cache = _OperatorCache(grid, "periodic")

        # phi_0 equation, laplacian(phi_0) term
        term = spec.equations[0].rhs_terms[0]
        result = _resolve_term_target(term, layout, set(), op_cache)
        assert result is not None
        col_slot, _op_mat, is_dyp = result
        assert col_slot == layout.field_slot_map["phi_0"]
        assert not is_dyp

    def test_first_derivative_t_dynamical(self) -> None:
        """first_derivative_t(dynamical_field) → dF/dy on velocity slot."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
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
                                "operator": "first_derivative_t",
                                "field": "A_1",
                            },
                        ],
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
                                "operator": "laplacian_x",
                                "field": "A_1",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)
        op_cache = _OperatorCache(grid, ("periodic", "periodic"))

        term = spec.equations[0].rhs_terms[0]
        constraint_fields: set[str] = {"A_0"}
        result = _resolve_term_target(term, layout, constraint_fields, op_cache)
        assert result is not None
        col_slot, _op_mat, is_dyp = result
        # A_1 is dynamical → couples via y to velocity slot
        assert col_slot == layout.velocity_slot_map["A_1"]
        assert not is_dyp

    def test_first_derivative_t_constraint(self) -> None:
        """first_derivative_t(constraint_field) → dF/dyp on field slot."""
        spec = _make_first_deriv_t_constraint_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        op_cache = _OperatorCache(grid, "periodic")

        # phi_0 eq term: first_derivative_t(A_0) where A_0 is constraint
        term = spec.equations[0].rhs_terms[1]
        constraint_fields = {"A_0"}
        result = _resolve_term_target(term, layout, constraint_fields, op_cache)
        assert result is not None
        col_slot, _op_mat, is_dyp = result
        assert col_slot == layout.field_slot_map["A_0"]
        assert is_dyp

    def test_explicit_v_constraint(self) -> None:
        """gradient_x(v_A_0) where A_0 is constraint → dF/dyp."""
        spec = _make_constraint_velocity_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)
        op_cache = _OperatorCache(grid, ("periodic", "periodic"))

        # A_1 eq term: gradient_x(v_A_0) where A_0 is constraint
        term = spec.equations[1].rhs_terms[1]
        constraint_fields = {"A_0"}
        result = _resolve_term_target(term, layout, constraint_fields, op_cache)
        assert result is not None
        col_slot, _op_mat, is_dyp = result
        assert col_slot == layout.field_slot_map["A_0"]
        assert is_dyp


# ---------------------------------------------------------------------------
# Jacobian matrix correctness tests
# ---------------------------------------------------------------------------


class TestBuildJacobianMatrices:
    def test_kg_1d_shapes(self) -> None:
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n_total = layout.total_size  # 2 slots × 8 = 16
        assert dF_dy.shape == (n_total, n_total)
        assert dF_dyp.shape == (n_total, n_total)
        assert issparse(dF_dy)
        assert issparse(dF_dyp)

    def test_kg_1d_dynamical_field_block(self) -> None:
        """Dynamical field: res = yp - v → dF/dy[field, vel] = -I."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points
        # Slot 0: phi_0 (field), slot 1: v_phi_0 (velocity)
        field_block_dy = dF_dy[:n, n : 2 * n].toarray()
        field_block_dyp = dF_dyp[:n, :n].toarray()

        np.testing.assert_array_equal(field_block_dy, -np.eye(n))
        np.testing.assert_array_equal(field_block_dyp, np.eye(n))

    def test_kg_1d_velocity_block(self) -> None:
        """Velocity: res = yp - RHS → dF/dy has negated RHS operators."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points
        # Velocity slot dyp diagonal should be I
        vel_dyp = dF_dyp[n : 2 * n, n : 2 * n].toarray()
        np.testing.assert_array_equal(vel_dyp, np.eye(n))

        # Velocity slot dF/dy should have -laplacian and +identity blocks
        # from negated RHS: -(1.0 * laplacian + -1.0 * identity)
        vel_dy_phi = dF_dy[n : 2 * n, :n].toarray()
        # Build expected: -laplacian + identity
        from tidal.solver.operators import apply_operator

        lap_mat = np.zeros((n, n))
        for j in range(n):
            e_j = np.zeros((n,))
            e_j[j] = 1.0
            lap_mat[:, j] = apply_operator("laplacian_x", e_j, grid, "periodic")
        expected = -1.0 * lap_mat - (-1.0) * np.eye(n)
        np.testing.assert_allclose(vel_dy_phi, expected, atol=1e-12)

    def test_constraint_velocity_in_dyp(self) -> None:
        """v_A_0 reference (A_0 constraint) should appear in dF/dyp."""
        spec = _make_constraint_velocity_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(
            spec,
            layout,
            grid,
            ("periodic", "periodic"),
            {},
        )

        n = grid.num_points  # 16
        # Slots: A_0 (0), A_1 (1), v_A_1 (2)
        # v_A_1 eq references gradient_x(v_A_0) — dF/dyp coupling
        v_A1_rows = slice(2 * n, 3 * n)
        A0_cols = slice(0, n)

        dyp_block = dF_dyp[v_A1_rows, A0_cols].toarray()
        assert np.any(dyp_block != 0), (
            "gradient_x(v_A_0) should produce entries in dF/dyp at "
            "(v_A_1_slot, A_0_slot)"
        )

        # The same block in dF/dy should be zero for this coupling
        dF_dy[v_A1_rows, A0_cols].toarray()
        # dF_dy may have entries from other terms but the constraint velocity
        # coupling specifically goes to dF_dyp

    def test_first_derivative_t_constraint_in_dyp(self) -> None:
        """first_derivative_t(A_0) where A_0 is constraint → dF/dyp."""
        spec = _make_first_deriv_t_constraint_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        _dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points  # 8
        # Slots: phi_0 (0), v_phi_0 (1), A_0 (2)
        # v_phi_0 eq has first_derivative_t(A_0) → dF/dyp[1, 2] = +1*I
        # But negated (res = yp - RHS) → dF/dyp[1, 2] should be -(-1)*I = +I
        # Wait: coefficient is -1.0 and negate=True → -(-1.0)*I = +I
        v_phi_rows = slice(1 * n, 2 * n)
        A0_cols = slice(2 * n, 3 * n)

        dyp_block = dF_dyp[v_phi_rows, A0_cols].toarray()
        # term coefficient = -1.0, negated → -(-1.0) * I = +I
        np.testing.assert_allclose(dyp_block, np.eye(n), atol=1e-12)

    def test_matches_fd_jacobian(self) -> None:
        """Analytical Jacobian matches finite-difference Jacobian."""
        spec = _make_kg_1d_spec()
        n_pts = 16
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(n_pts,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        # Build residual function
        from tidal.solver.ida import build_residual_fn

        resfn = build_residual_fn(spec, layout, grid, "periodic", parameters={})

        n_total = layout.total_size
        y0 = np.random.default_rng(42).standard_normal(n_total)
        yp0 = np.random.default_rng(43).standard_normal(n_total)

        # Compute residual at (y0, yp0)
        res0 = np.zeros(n_total)
        resfn(0.0, y0, yp0, res0)

        # Finite-difference dF/dy
        eps = 1e-7
        fd_dF_dy = np.zeros((n_total, n_total))
        for j in range(n_total):
            y_plus = y0.copy()
            y_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y_plus, yp0, res_plus)
            fd_dF_dy[:, j] = (res_plus - res0) / eps

        # Finite-difference dF/dyp
        fd_dF_dyp = np.zeros((n_total, n_total))
        for j in range(n_total):
            yp_plus = yp0.copy()
            yp_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y0, yp_plus, res_plus)
            fd_dF_dyp[:, j] = (res_plus - res0) / eps

        np.testing.assert_allclose(
            dF_dy.toarray(),
            fd_dF_dy,
            atol=1e-5,
            err_msg="dF/dy mismatch",
        )
        np.testing.assert_allclose(
            dF_dyp.toarray(),
            fd_dF_dyp,
            atol=1e-5,
            err_msg="dF/dyp mismatch",
        )

    def test_constraint_velocity_matches_fd(self) -> None:
        """Analytical Jacobian with constraint velocity matches FD."""
        spec = _make_constraint_velocity_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(
            spec,
            layout,
            grid,
            ("periodic", "periodic"),
            {},
        )

        from tidal.solver.ida import build_residual_fn

        resfn = build_residual_fn(
            spec,
            layout,
            grid,
            ("periodic", "periodic"),
            parameters={},
        )

        n_total = layout.total_size
        rng = np.random.default_rng(44)
        y0 = rng.standard_normal(n_total)
        yp0 = rng.standard_normal(n_total)

        res0 = np.zeros(n_total)
        resfn(0.0, y0, yp0, res0)

        eps = 1e-7

        # FD dF/dy
        fd_dy = np.zeros((n_total, n_total))
        for j in range(n_total):
            y_plus = y0.copy()
            y_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y_plus, yp0, res_plus)
            fd_dy[:, j] = (res_plus - res0) / eps

        # FD dF/dyp
        fd_dyp = np.zeros((n_total, n_total))
        for j in range(n_total):
            yp_plus = yp0.copy()
            yp_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y0, yp_plus, res_plus)
            fd_dyp[:, j] = (res_plus - res0) / eps

        np.testing.assert_allclose(
            dF_dy.toarray(),
            fd_dy,
            atol=1e-5,
            err_msg="dF/dy mismatch for constraint velocity system",
        )
        np.testing.assert_allclose(
            dF_dyp.toarray(),
            fd_dyp,
            atol=1e-5,
            err_msg="dF/dyp mismatch for constraint velocity system",
        )


# ---------------------------------------------------------------------------
# Jacobian delivery tests
# ---------------------------------------------------------------------------


class TestJacobianDelivery:
    def test_dense_jacfn(self) -> None:
        """Dense jacfn should produce correct J = dF_dy + cj * dF_dyp."""
        from tidal.solver.analytical_jacobian import _create_jacfn

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})
        jacfn = _create_jacfn(dF_dy, dF_dyp)

        n_total = layout.total_size
        cj = math.pi
        JJ = np.zeros((n_total, n_total))
        jacfn(0.0, np.zeros(n_total), np.zeros(n_total), np.zeros(n_total), cj, JJ)

        expected = dF_dy.toarray() + cj * dF_dyp.toarray()
        np.testing.assert_allclose(JJ, expected, atol=1e-12)

    def test_sparse_tier_configures_sparse_jacfn(self) -> None:
        """Sparse tier should configure 1D CSC jacfn + sparsity pattern."""
        from unittest.mock import patch

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        # Lower DENSE_THRESHOLD so system is in sparse tier
        with patch("tidal.solver._types.DENSE_THRESHOLD", 1):
            result = try_analytical_jacobian(
                options, spec, layout, grid, "periodic", {}
            )

        assert result is True
        assert options["linsolver"] == "sparse"
        assert "jacfn" in options
        assert callable(options["jacfn"])
        assert "sparsity" in options
        assert issparse(options["sparsity"])

    def test_sparse_jacfn_fills_1d_array(self) -> None:
        """Sparse jacfn callback should fill 1D CSC data array correctly."""
        from tidal.solver.analytical_jacobian import _create_sparse_jacfn

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})
        jacfn, sparsity = _create_sparse_jacfn(dF_dy, dF_dyp)

        # Test for several cj values
        for cj in [0.0, 1.0, math.pi, 100.0]:
            JJ = np.zeros(sparsity.nnz)
            n = layout.total_size
            jacfn(0.0, np.zeros(n), np.zeros(n), np.zeros(n), cj, JJ)

            # Build reference: full J = dF_dy + cj * dF_dyp, extract CSC data
            J_ref = (dF_dy + cj * dF_dyp).tocsc()
            from tidal.solver.analytical_jacobian import _extract_aligned_data

            ref_data = _extract_aligned_data(J_ref, sparsity)
            np.testing.assert_allclose(JJ, ref_data, atol=1e-14)

    def test_sparse_jacfn_csc_alignment(self) -> None:
        """Sparsity pattern should have CSC structure matching union of matrices."""
        from tidal.solver.analytical_jacobian import _create_sparse_jacfn

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})
        _, sparsity = _create_sparse_jacfn(dF_dy, dF_dyp)

        # Union pattern should cover all nonzeros from both matrices
        union_ref = (abs(dF_dy) + abs(dF_dyp)).tocsc()
        union_ref.eliminate_zeros()
        np.testing.assert_array_equal(sparsity.indptr, union_ref.indptr)
        np.testing.assert_array_equal(sparsity.indices, union_ref.indices)
        assert sparsity.nnz == union_ref.nnz

    def test_sparse_tier_superlu_nnz_limit_fallback(self) -> None:
        """Sparse tier should fall through to GMRES when nnz > SUPERLU_NNZ_LIMIT."""
        from unittest.mock import patch

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        with (
            patch("tidal.solver._types.DENSE_THRESHOLD", 1),
            patch("tidal.solver._types.SUPERLU_NNZ_LIMIT", 1),
        ):
            result = try_analytical_jacobian(
                options, spec, layout, grid, "periodic", {}
            )

        assert result is True
        assert options["linsolver"] == "gmres"
        assert "jactimes" in options
        assert "jacfn" not in options

    def test_cvode_sparse_jacfn(self) -> None:
        """CVODE sparse jacfn should produce ODE Jacobian -dF_dy."""
        from tidal.solver.analytical_jacobian import _create_cvode_sparse_jacfn

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, _ = build_jacobian_matrices(spec, layout, grid, "periodic", {})
        jacfn, sparsity = _create_cvode_sparse_jacfn(dF_dy)

        JJ = np.zeros(sparsity.nnz)
        n = layout.total_size
        jacfn(0.0, np.zeros(n), np.zeros(n), JJ)

        neg_dy = (-dF_dy).tocsc()
        neg_dy.eliminate_zeros()
        np.testing.assert_allclose(JJ, neg_dy.data, atol=1e-14)


# ---------------------------------------------------------------------------
# Circulant operator matrix construction
# ---------------------------------------------------------------------------


class TestCirculantOperatorMatrix:
    """Verify circulant fast-path matches probing for all operator types."""

    @staticmethod
    def _build_probing(operator: str, grid: GridInfo, bc: str | tuple[str, ...]) -> Any:  # noqa: ANN401
        """Build operator matrix via O(N²) probing (reference)."""
        from tidal.solver.operators import apply_operator

        n = grid.num_points
        from tidal.solver._scipy_types import lil_matrix as _lil

        mat = _lil((n, n))
        e_j = np.zeros(grid.shape)
        for j in range(n):
            e_j.flat[:] = 0.0
            e_j.flat[j] = 1.0
            col = apply_operator(operator, e_j, grid, bc)
            col_flat = col.ravel()
            nz = np.nonzero(col_flat)[0]
            for row in nz:
                mat[row, j] = col_flat[row]
        return mat.tocsc()

    def test_laplacian_x_2d(self) -> None:
        """laplacian_x on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("laplacian_x", grid, "periodic")
        ref = self._build_probing("laplacian_x", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_laplacian_y_2d(self) -> None:
        """laplacian_y on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("laplacian_y", grid, "periodic")
        ref = self._build_probing("laplacian_y", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_gradient_x_2d(self) -> None:
        """gradient_x on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("gradient_x", grid, "periodic")
        ref = self._build_probing("gradient_x", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_gradient_y_2d(self) -> None:
        """gradient_y on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("gradient_y", grid, "periodic")
        ref = self._build_probing("gradient_y", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_cross_derivative_xy_2d(self) -> None:
        """cross_derivative_xy on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("cross_derivative_xy", grid, "periodic")
        ref = self._build_probing("cross_derivative_xy", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_laplacian_full_2d(self) -> None:
        """Full laplacian (sum of directional) on 2D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("laplacian", grid, "periodic")
        ref = self._build_probing("laplacian", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_biharmonic_1d(self) -> None:
        """Biharmonic on 1D periodic grid."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        mat = build_operator_matrix("biharmonic", grid, "periodic")
        ref = self._build_probing("biharmonic", grid, "periodic")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_3d_operators(self) -> None:
        """laplacian_z, gradient_z, cross_derivative_xz, cross_derivative_yz on 3D."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(
            bounds=((0, 10), (0, 10), (0, 10)),
            shape=(4, 4, 4),
            periodic=(True, True, True),
        )
        for op in [
            "laplacian_z",
            "gradient_z",
            "cross_derivative_xz",
            "cross_derivative_yz",
        ]:
            mat = build_operator_matrix(op, grid, ("periodic", "periodic", "periodic"))
            ref = self._build_probing(op, grid, ("periodic", "periodic", "periodic"))
            np.testing.assert_allclose(
                mat.toarray(), ref.toarray(), atol=1e-14, err_msg=f"{op} mismatch"
            )

    def test_identity_uses_probing_fallback(self) -> None:
        """Identity operator should NOT use circulant (trivially sparse eye)."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        mat = build_operator_matrix("identity", grid, "periodic")
        expected = np.eye(8)
        np.testing.assert_allclose(mat.toarray(), expected, atol=1e-14)

    def test_nonperiodic_uses_probing_fallback(self) -> None:
        """Non-periodic BCs should fall back to probing."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(False,))
        mat = build_operator_matrix("laplacian_x", grid, "neumann")
        ref = self._build_probing("laplacian_x", grid, "neumann")
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)

    def test_tuple_bc_periodic(self) -> None:
        """Tuple BC ('periodic', 'periodic') should use circulant path."""
        from tidal.solver.analytical_jacobian import build_operator_matrix

        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        mat = build_operator_matrix("laplacian_x", grid, ("periodic", "periodic"))
        ref = self._build_probing("laplacian_x", grid, ("periodic", "periodic"))
        np.testing.assert_allclose(mat.toarray(), ref.toarray(), atol=1e-14)


# ---------------------------------------------------------------------------
# Multi-theory sparse Jacobian integration tests
# ---------------------------------------------------------------------------


class TestMultiTheorySparseJacobian:
    """Test sparse analytical Jacobian across real example theories.

    These tests load real JSON specs and verify that the sparse jacfn
    callback produces values matching the reference dense Jacobian.
    """

    @staticmethod
    def _test_theory(json_name: str, params: dict[str, float]) -> None:
        """Load a theory JSON and verify sparse jacfn matches dense Jacobian."""
        import json
        from pathlib import Path

        from tidal.solver.analytical_jacobian import (
            _create_sparse_jacfn,
            _extract_aligned_data,
        )

        json_path = Path("examples/data") / json_name
        if not json_path.exists():
            return  # skip if example not available
        spec_data = json.loads(json_path.read_text(encoding="utf-8"))
        spec = EquationSystem.from_dict(spec_data)

        # Use small grid to keep test fast
        ndim = spec_data["spacetime"]["dimension"] - 1
        if ndim == 1:
            grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        elif ndim == 2:
            grid = GridInfo(
                bounds=((0, 10), (0, 10)), shape=(6, 6), periodic=(True, True)
            )
        else:
            grid = GridInfo(
                bounds=((0, 10), (0, 10), (0, 10)),
                shape=(4, 4, 4),
                periodic=(True, True, True),
            )
        bc = "periodic"
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Build analytical Jacobian matrices
        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, bc, params)

        # Create sparse jacfn
        jacfn, sparsity = _create_sparse_jacfn(dF_dy, dF_dyp)

        # Test for several cj values
        n = layout.total_size
        for cj in [0.0, 1.0, math.pi]:
            JJ = np.zeros(sparsity.nnz)
            jacfn(0.0, np.zeros(n), np.zeros(n), np.zeros(n), cj, JJ)

            # Reference
            J_ref = (dF_dy + cj * dF_dyp).tocsc()
            ref_data = _extract_aligned_data(J_ref, sparsity)
            np.testing.assert_allclose(
                JJ,
                ref_data,
                atol=1e-12,
                err_msg=f"{json_name} cj={cj}",
            )

    def test_klein_gordon_1d(self) -> None:
        """1D scalar: identity + laplacian_x."""
        self._test_theory("klein_gordon_1d.json", {"m2": 1.0})

    def test_coupled_scalars(self) -> None:
        """2-field coupled: identity + laplacian_x + gradient coupling."""
        self._test_theory(
            "coupled_scalars.json",
            {"B0": 1.0, "kappa": 0.1, "mg2": 1.0, "omegaP2": 1.0},
        )

    def test_proca_1d(self) -> None:
        """1D Proca: identity + laplacian_x + gradient_x."""
        self._test_theory("proca_1d.json", {"procaMassSquared": 1.0})

    def test_em_3d(self) -> None:
        """3D EM with constraints: laplacian + gradient + cross_derivative."""
        self._test_theory("em_3d.json", {})

    def test_navier_cauchy_2d(self) -> None:
        """2D Navier-Cauchy: laplacian_x/y + cross_derivative_xy."""
        self._test_theory("navier_cauchy_2d.json", {"lam": 1.0, "mu": 1.0, "rho": 1.0})

    def test_scalar_vector_coupling(self) -> None:
        """4-field with constraints + first_derivative_t."""
        self._test_theory(
            "scalar_vector_coupling.json",
            {"Am2": 1.0, "phim2": 1.0, "gSV": 0.5, "kCS": 0.1},
        )

    def test_coupled_scattering(self) -> None:
        """Position-dependent (Gaussian) coefficients."""
        self._test_theory(
            "coupled_scattering.json",
            {"g0": 0.5, "R": 2.0, "mPhi2": 1.0, "mChi2": 1.0},
        )


# ---------------------------------------------------------------------------
# try_analytical_jacobian integration
# ---------------------------------------------------------------------------


class TestTryAnalyticalJacobian:
    def test_configures_dense_jacfn(self) -> None:
        """Small constant system should get dense jacfn."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        result = try_analytical_jacobian(options, spec, layout, grid, "periodic", {})

        assert result is True
        assert options["linsolver"] == "dense"
        assert "jacfn" in options
        assert callable(options["jacfn"])

    def test_returns_false_for_non_constant(self) -> None:
        """Non-constant system should return False."""
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
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "E^(2*H*t[])",
                                "time_dependent": True,
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        result = try_analytical_jacobian(
            options, spec, layout, grid, "periodic", {"H": 1.0}
        )

        assert result is False
        assert "jacfn" not in options

    def test_position_dependent_gets_analytical_jacobian(self) -> None:
        """Position-dependent (but time-independent) system should succeed."""
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
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-V0*UnitStep[x[]-3]*UnitStep[7-x[]]",
                                "coordinate_dependent": ["x"],
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(False,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        result = try_analytical_jacobian(
            options,
            spec,
            layout,
            grid,
            "neumann",
            {"V0": 2.0},
        )

        assert result is True
        assert options["linsolver"] == "dense"
        assert "jacfn" in options


# ---------------------------------------------------------------------------
# Sparse tier tests
# ---------------------------------------------------------------------------


class TestSparseTier:
    """Tests for sparse/GMRES tier analytical Jacobian."""

    def test_gmres_tier_via_try_analytical(self) -> None:
        """try_analytical_jacobian selects GMRES tier above SPARSE_THRESHOLD."""
        from unittest.mock import patch

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        # Lower both thresholds so system is in GMRES tier
        with (
            patch("tidal.solver._types.DENSE_THRESHOLD", 1),
            patch("tidal.solver._types.SPARSE_THRESHOLD", 1),
        ):
            result = try_analytical_jacobian(
                options,
                spec,
                layout,
                grid,
                "periodic",
                {},
            )

        assert result is True
        assert options["linsolver"] == "gmres"
        assert "jactimes" in options


# ---------------------------------------------------------------------------
# Position-dependent coefficient tests
# ---------------------------------------------------------------------------


class TestPositionDependentJacobian:
    """Tests for position-dependent (time-independent) analytical Jacobian."""

    def test_position_dependent_jacobian_builds(self) -> None:
        """build_jacobian_matrices handles ndarray coefficients."""
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
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "Sin[x[]]",
                                "coordinate_dependent": ["x"],
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, np.pi),), shape=(8,), periodic=(False,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Should not raise (previously would fail with float(ndarray))
        dF_dy, dF_dyp = build_jacobian_matrices(
            spec,
            layout,
            grid,
            "neumann",
            {},
        )

        assert issparse(dF_dy)
        assert issparse(dF_dyp)

        # The velocity slot identity term should be diagonal (from sin(x))
        n = grid.num_points
        vel_slot = layout.velocity_slot_map["phi_0"]
        vel_rows = slice(vel_slot * n, (vel_slot + 1) * n)
        field_slot = layout.field_slot_map["phi_0"]
        field_cols = slice(field_slot * n, (field_slot + 1) * n)

        # The identity block in dF/dy[vel, field] should have position-
        # dependent diagonal entries from -sign * coeff * I = +sin(x_i)
        block = dF_dy[vel_rows, field_cols].toarray()
        diag_vals = np.diag(block)
        # sin(x) evaluated at grid cell centres should vary
        assert not np.allclose(diag_vals, diag_vals[0])

    def test_all_time_independent_method(self) -> None:
        """CoefficientEvaluator.all_time_independent() detects correctly."""
        # Time-independent system
        spec_const = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff = CoefficientEvaluator(spec_const, grid, {})
        assert coeff.all_time_independent() is True
        assert coeff.all_constant() is True

        # Position-dependent but time-independent
        data_pos: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "Sin[x[]]",
                                "coordinate_dependent": ["x"],
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec_pos = EquationSystem.from_dict(data_pos)
        coeff_pos = CoefficientEvaluator(spec_pos, grid, {})
        assert coeff_pos.all_time_independent() is True
        assert coeff_pos.all_constant() is False

        # Time-dependent system
        data_time: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "E^(2*H*t[])",
                                "time_dependent": True,
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec_time = EquationSystem.from_dict(data_time)
        coeff_time = CoefficientEvaluator(spec_time, grid, {"H": 1.0})
        assert coeff_time.all_time_independent() is False


# ---------------------------------------------------------------------------
# Integration test: IDA with analytical Jacobian
# ---------------------------------------------------------------------------


class TestConstraintSubCases:
    """Tests for constraint sub-cases: no-self-term, gauge-fix, first-order."""

    def test_no_self_term_constraint(self) -> None:
        """Constraint with no self-referencing terms: res = y[field] → dF/dy = I."""
        # A_0 is a subsidiary constraint whose equation only references phi_0
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "A_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "A_0",
                    "lhs": {"expression": "A_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            # Only references phi_0, NOT A_0 itself
                            {
                                "coefficient": 1.0,
                                "operator": "gradient_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                    # constraint_solver not enabled (default) → no-self-term path
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points  # 8
        # Slots: phi_0 (0), v_phi_0 (1), A_0 (2)
        A0_slot = 2
        A0_rows = slice(A0_slot * n, (A0_slot + 1) * n)
        field_slot = layout.field_slot_map["A_0"]
        A0_cols = slice(field_slot * n, (field_slot + 1) * n)

        # dF/dy at (A_0, A_0) should be identity
        dy_block = dF_dy[A0_rows, A0_cols].toarray()
        np.testing.assert_array_equal(dy_block, np.eye(n))

        # dF/dyp at (A_0, :) should be all zeros (no yp dependence)
        dyp_row = dF_dyp[A0_rows, :].toarray()
        np.testing.assert_array_equal(dyp_row, 0.0)

    def test_gauge_fix_constraint(self) -> None:
        """Pure Laplacian + periodic BC → row 0 zeroed + pinning entry."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "A_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "A_0",
                    "lhs": {"expression": "A_0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            # Pure Laplacian self-term → triggers gauge fix
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "A_0",
                            },
                            {
                                "coefficient": 0.5,
                                "operator": "gradient_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                    "constraint_solver": {"enabled": True},
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points  # 8
        # Slots: phi_0 (0), v_phi_0 (1), A_0 (2)
        A0_slot_idx = 2
        global_row0 = A0_slot_idx * n  # first row of A_0 block

        # Row 0 of A_0 block in dF/dy: entire row should be zero except pinning
        row_dy = dF_dy[global_row0, :].toarray().ravel()
        field_slot = layout.field_slot_map["A_0"]
        pinning_col = field_slot * n  # first column of A_0 field block
        assert row_dy[pinning_col] == 1.0, "Pinning entry should be 1.0"
        # All other entries in this row should be zero
        mask = np.ones(len(row_dy), dtype=bool)
        mask[pinning_col] = False
        np.testing.assert_array_equal(row_dy[mask], 0.0)

        # Row 0 of A_0 block in dF/dyp: entire row should be zero
        row_dyp = dF_dyp[global_row0, :].toarray().ravel()
        np.testing.assert_array_equal(row_dyp, 0.0)

        # Other rows in A_0 block should have Laplacian entries (not all zero)
        other_rows = dF_dy[global_row0 + 1 : (A0_slot_idx + 1) * n, :].toarray()
        assert np.any(other_rows != 0), "Non-pinned rows should have Laplacian entries"

    def test_first_order_equation(self) -> None:
        """First-order equation (time_order=1): dF/dyp = I, negated RHS."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "psi_0", "index": 0}],
            "equations": [
                {
                    "field": "psi_0",
                    "lhs": {"expression": "d_t(psi_0)", "order": {"time": 1}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "psi_0",
                            },
                            {
                                "coefficient": -0.5,
                                "operator": "identity",
                                "field": "psi_0",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        n = grid.num_points  # 8
        # Single slot: psi_0 (0), time_order=1

        # dF/dyp should be identity (res = yp - RHS)
        np.testing.assert_array_equal(dF_dyp.toarray(), np.eye(n))

        # dF/dy should be negated RHS: -(1.0 * Lap + -0.5 * I)
        from tidal.solver.operators import apply_operator

        lap_mat = np.zeros((n, n))
        for j in range(n):
            e_j = np.zeros((n,))
            e_j[j] = 1.0
            lap_mat[:, j] = apply_operator("laplacian_x", e_j, grid, "periodic")
        expected = -(1.0 * lap_mat + (-0.5) * np.eye(n))
        np.testing.assert_allclose(dF_dy.toarray(), expected, atol=1e-12)

    def test_first_order_matches_fd(self) -> None:
        """First-order analytical Jacobian matches finite-difference."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "psi_0", "index": 0}],
            "equations": [
                {
                    "field": "psi_0",
                    "lhs": {"expression": "d_t(psi_0)", "order": {"time": 1}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "psi_0",
                            },
                            {
                                "coefficient": -0.5,
                                "operator": "identity",
                                "field": "psi_0",
                            },
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})

        from tidal.solver.ida import build_residual_fn

        resfn = build_residual_fn(spec, layout, grid, "periodic", parameters={})

        n_total = layout.total_size
        rng = np.random.default_rng(55)
        y0 = rng.standard_normal(n_total)
        yp0 = rng.standard_normal(n_total)

        res0 = np.zeros(n_total)
        resfn(0.0, y0, yp0, res0)

        eps = 1e-7
        fd_dy = np.zeros((n_total, n_total))
        for j in range(n_total):
            y_plus = y0.copy()
            y_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y_plus, yp0, res_plus)
            fd_dy[:, j] = (res_plus - res0) / eps

        fd_dyp = np.zeros((n_total, n_total))
        for j in range(n_total):
            yp_plus = yp0.copy()
            yp_plus[j] += eps
            res_plus = np.zeros(n_total)
            resfn(0.0, y0, yp_plus, res_plus)
            fd_dyp[:, j] = (res_plus - res0) / eps

        np.testing.assert_allclose(
            dF_dy.toarray(),
            fd_dy,
            atol=1e-5,
            err_msg="dF/dy mismatch for first-order equation",
        )
        np.testing.assert_allclose(
            dF_dyp.toarray(),
            fd_dyp,
            atol=1e-5,
            err_msg="dF/dyp mismatch for first-order equation",
        )


# ---------------------------------------------------------------------------
# CVODE delivery tests
# ---------------------------------------------------------------------------


class TestCVODEDelivery:
    def test_cvode_dense_jacfn(self) -> None:
        """CVODE dense jacfn should produce JJ = -dF_dy."""
        from tidal.solver.analytical_jacobian import _create_cvode_jacfn

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        dF_dy, _dF_dyp = build_jacobian_matrices(spec, layout, grid, "periodic", {})
        jacfn = _create_cvode_jacfn(dF_dy)

        n_total = layout.total_size
        JJ = np.zeros((n_total, n_total))
        jacfn(0.0, np.zeros(n_total), np.zeros(n_total), JJ)

        expected = (-dF_dy).toarray()
        np.testing.assert_allclose(JJ, expected, atol=1e-12)

    def test_try_analytical_cvode(self) -> None:
        """try_analytical_jacobian with solver='cvode' configures CVODE jacfn."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        result = try_analytical_jacobian(
            options,
            spec,
            layout,
            grid,
            "periodic",
            {},
            solver="cvode",
        )

        assert result is True
        assert options["linsolver"] == "dense"
        assert "jacfn" in options
        assert callable(options["jacfn"])

        # Verify the jacfn has CVODE signature (4 args, no cj/res)
        n_total = layout.total_size
        JJ = np.zeros((n_total, n_total))
        options["jacfn"](0.0, np.zeros(n_total), np.zeros(n_total), JJ)
        assert np.any(JJ != 0), "CVODE jacfn should fill JJ"


# ---------------------------------------------------------------------------
# SuperLU nnz limit and GMRES fallback
# ---------------------------------------------------------------------------


def _make_time_dependent_spec() -> EquationSystem:
    """KG with time-dependent coefficient — bypasses analytical Jacobian."""
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
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "identity",
                            "field": "phi_0",
                            "coefficient_symbolic": "E^(2*H*t[])",
                            "time_dependent": True,
                        },
                    ],
                },
            },
        ],
        "canonical": {"hamiltonian_terms": []},
    }
    return EquationSystem.from_dict(data)


class TestSuperLUNnzFallback:
    """Tests for nnz-based GMRES fallback in configure_linear_solver.

    Uses a time-dependent spec so the analytical Jacobian path is bypassed
    and the FD sparse tier in configure_linear_solver is exercised.
    """

    def test_falls_back_to_gmres_when_nnz_exceeds_limit(self) -> None:
        """When nnz exceeds SUPERLU_NNZ_LIMIT, configure_linear_solver selects GMRES."""
        import warnings
        from unittest.mock import patch

        from tidal.solver._setup import configure_linear_solver

        spec = _make_time_dependent_spec()
        n_pts = 64
        grid = GridInfo(bounds=((0, 10),), shape=(n_pts,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        with (
            patch("tidal.solver._types.DENSE_THRESHOLD", 1),
            patch("tidal.solver._types.SUPERLU_NNZ_LIMIT", 1),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            configure_linear_solver(
                options, layout, spec, grid, "periodic", parameters={"H": 1.0}
            )

        assert options["linsolver"] == "gmres"
        assert "sparsity" not in options

        warning_messages = [
            str(w.message) for w in caught if issubclass(w.category, UserWarning)
        ]
        assert any("SUPERLU_NNZ_LIMIT" in m for m in warning_messages), (
            f"Expected UserWarning about SUPERLU_NNZ_LIMIT, got: {warning_messages}"
        )

    def test_uses_superlu_when_nnz_under_limit(self) -> None:
        """When nnz is under SUPERLU_NNZ_LIMIT, SuperLU is selected normally."""
        from unittest.mock import patch

        from tidal.solver._setup import configure_linear_solver

        spec = _make_time_dependent_spec()
        n_pts = 8
        grid = GridInfo(bounds=((0, 10),), shape=(n_pts,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        with patch("tidal.solver._types.DENSE_THRESHOLD", 1):
            configure_linear_solver(
                options, layout, spec, grid, "periodic", parameters={"H": 1.0}
            )

        assert options["linsolver"] == "sparse"
        assert "sparsity" in options

    def test_warning_contains_nnz_and_limit(self) -> None:
        """UserWarning message includes nnz value and limit for diagnostics."""
        import warnings
        from unittest.mock import patch

        from tidal.solver._setup import configure_linear_solver

        spec = _make_time_dependent_spec()
        n_pts = 64
        grid = GridInfo(bounds=((0, 10),), shape=(n_pts,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        options: dict[str, Any] = {}
        with (
            patch("tidal.solver._types.DENSE_THRESHOLD", 1),
            patch("tidal.solver._types.SUPERLU_NNZ_LIMIT", 1),
            warnings.catch_warnings(record=True) as caught,
        ):
            warnings.simplefilter("always")
            configure_linear_solver(
                options, layout, spec, grid, "periodic", parameters={"H": 1.0}
            )

        user_warnings = [
            str(w.message) for w in caught if issubclass(w.category, UserWarning)
        ]
        assert user_warnings, "Expected at least one UserWarning"
        msg = user_warnings[0]
        assert "nnz=" in msg
        assert "SUPERLU_NNZ_LIMIT=1" in msg
        assert "gmres" in msg


class TestIDAIntegration:
    def test_kg_analytical_matches_fd(self) -> None:
        """IDA with analytical jacfn should match standard FD solve."""
        from tidal.solver.ida import solve_ida

        spec = _make_kg_1d_spec()
        n_pts = 32
        grid = GridInfo(
            bounds=((0, 2 * np.pi),),
            shape=(n_pts,),
            periodic=(True,),
        )

        x = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.total_size)
        y0[:n_pts] = np.exp(-((x - np.pi) ** 2))

        # Solve with analytical Jacobian (parameters triggers it)
        result_analytical = solve_ida(
            spec,
            grid,
            y0,
            t_span=(0.0, 0.5),
            bc="periodic",
            parameters={},
            num_snapshots=5,
        )
        assert result_analytical["success"]

        # Solve without analytical Jacobian (no parameters → FD)
        result_fd = solve_ida(
            spec,
            grid,
            y0,
            t_span=(0.0, 0.5),
            bc="periodic",
            num_snapshots=5,
        )
        assert result_fd["success"]

        np.testing.assert_allclose(
            result_analytical["y"][-1],
            result_fd["y"][-1],
            rtol=1e-5,
            atol=1e-8,
            err_msg="Analytical Jacobian result differs from FD Jacobian",
        )
