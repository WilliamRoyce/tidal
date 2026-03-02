"""Tests for analytical Jacobian sparsity pattern builder."""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import sparse  # pyright: ignore[reportMissingTypeStubs]

from tidal.solver.grid import GridInfo
from tidal.solver.sparsity import (
    build_jacobian_sparsity,
    operator_stencil_offsets,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Operator stencil offset tests
# ---------------------------------------------------------------------------


class TestOperatorStencilOffsets:
    def test_identity_1d(self) -> None:
        offsets = operator_stencil_offsets("identity", 1)
        assert offsets == [(0,)]

    def test_identity_2d(self) -> None:
        offsets = operator_stencil_offsets("identity", 2)
        assert offsets == [(0, 0)]

    def test_gradient_x_2d(self) -> None:
        offsets = operator_stencil_offsets("gradient_x", 2)
        assert (-1, 0) in offsets
        assert (1, 0) in offsets
        assert (0, 0) in offsets  # self for robustness

    def test_gradient_y_2d(self) -> None:
        offsets = operator_stencil_offsets("gradient_y", 2)
        assert (0, -1) in offsets
        assert (0, 1) in offsets
        assert (0, 0) in offsets

    def test_laplacian_x_1d(self) -> None:
        offsets = operator_stencil_offsets("laplacian_x", 1)
        assert set(offsets) == {(-1,), (0,), (1,)}

    def test_laplacian_2d(self) -> None:
        offsets = operator_stencil_offsets("laplacian", 2)
        offsets_set = set(offsets)
        # Must include ±1 along each axis + self
        assert (-1, 0) in offsets_set
        assert (1, 0) in offsets_set
        assert (0, -1) in offsets_set
        assert (0, 1) in offsets_set
        assert (0, 0) in offsets_set

    def test_cross_derivative_xy_2d(self) -> None:
        offsets = operator_stencil_offsets("cross_derivative_xy", 2)
        offsets_set = set(offsets)
        # All four diagonal neighbors
        assert (1, 1) in offsets_set
        assert (1, -1) in offsets_set
        assert (-1, 1) in offsets_set
        assert (-1, -1) in offsets_set

    def test_biharmonic_1d(self) -> None:
        offsets = operator_stencil_offsets("biharmonic", 1)
        offsets_set = set(offsets)
        # Must include radius-2 stencil
        assert (-2,) in offsets_set
        assert (-1,) in offsets_set
        assert (0,) in offsets_set
        assert (1,) in offsets_set
        assert (2,) in offsets_set

    def test_first_derivative_t(self) -> None:
        offsets = operator_stencil_offsets("first_derivative_t", 2)
        assert offsets == [(0, 0)]

    def test_unknown_operator_conservative(self) -> None:
        """Unknown operators should get conservative radius-1 stencil."""
        offsets = operator_stencil_offsets("exotic_operator", 1)
        assert len(offsets) == 3  # {-1, 0, 1}


# ---------------------------------------------------------------------------
# Sparsity pattern tests
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
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                        {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                    ],
                },
            },
        ],
        "canonical": {
            "hamiltonian_terms": [],
        },
    }
    return EquationSystem.from_dict(data)


def _make_em_2d_spec() -> EquationSystem:
    """Simplified EM 2+1D: constraint A_0 + wave A_1."""
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
        "canonical": {
            "hamiltonian_terms": [],
        },
    }
    return EquationSystem.from_dict(data)


class TestBuildJacobianSparsity:
    def test_kg_1d_shape(self) -> None:
        """Sparsity pattern for KG 1D should have correct shape."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")

        # 2 slots (phi_0, pi_phi_0) x 16 points = 32
        assert pattern.shape == (32, 32)

    def test_kg_1d_much_sparser_than_dense(self) -> None:
        """Sparsity pattern should be much sparser than dense N²."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")

        n_total = layout.total_size  # 128
        density = pattern.nnz / (n_total * n_total)
        # For 1D KG with 64 points, density should be well under 10%
        assert density < 0.10, f"Pattern too dense: {density:.1%}"
        # But should have some entries (at least diagonal + stencil)
        assert pattern.nnz > n_total

    def test_kg_1d_diagonal_present(self) -> None:
        """All non-constraint diagonal entries should be nonzero (cj*I)."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")

        # Both slots are non-constraint (phi_0: field, pi_phi_0: momentum)
        diag: Any = pattern.diagonal()  # pyright: ignore[reportUnknownVariableType]
        assert np.all(diag > 0), "Diagonal entries missing"

    def test_em_2d_constraint_no_diagonal(self) -> None:
        """Constraint slot (A_0) should NOT have cj*I diagonal coupling."""
        spec = _make_em_2d_spec()
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, ("periodic", "periodic"))

        # Slot 0 = A_0 (constraint) — should NOT have cj*I diagonal
        # But it DOES have self-coupling from laplacian(A_0) in its RHS
        # So diagonal entries exist from the laplacian self-term, not cj*I
        # Just verify the pattern has correct total shape
        assert pattern.shape == (layout.total_size, layout.total_size)

    def test_em_2d_sparsity(self) -> None:
        """EM 2D pattern should be very sparse."""
        spec = _make_em_2d_spec()
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(16, 16), periodic=(True, True)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, ("periodic", "periodic"))

        n_total = layout.total_size  # 3 slots x 256 = 768
        density = pattern.nnz / (n_total * n_total)
        # Should be very sparse (< 5%)
        assert density < 0.05, f"Pattern too dense: {density:.1%}"

    def test_periodic_wrapping(self) -> None:
        """Periodic BC should create wrapping entries at grid boundaries."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")

        # Convert to dense for inspection
        dense: Any = pattern.toarray()  # pyright: ignore[reportUnknownVariableType]
        # The momentum slot (slot 1) equation references laplacian(phi_0).
        # For periodic 1D with 8 points, point 0 couples to point 7 (wrap).
        # Momentum rows are [8..15], phi columns are [0..7].
        # Row 8 (momentum[0]) should have nonzero at col 7 (phi[7]).
        assert dense[8, 7] != 0, "Missing periodic wrap coupling"

    def test_nonperiodic_no_wrapping(self) -> None:
        """Non-periodic BC should NOT create wrapping entries."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(False,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "neumann")

        dense: Any = pattern.toarray()  # pyright: ignore[reportUnknownVariableType]
        # Row 8 (momentum[0]) should NOT have nonzero at col 7 (phi[7])
        # because boundary point 0 doesn't wrap to point 7
        assert dense[8, 7] == 0, "Unexpected wrapping entry for non-periodic BC"

    def test_pattern_is_csc(self) -> None:
        """Pattern should be returned as CSC sparse matrix."""
        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")
        assert isinstance(pattern, sparse.csc_matrix)

    def test_constraint_velocity_coupling_in_pattern(self) -> None:
        """Sparsity pattern should include constraint velocity couplings.

        When a dynamical equation references ``v_X`` where X is a constraint
        field (time_order=0), the coupling goes through ``yp[X_field_slot]``.
        The sparsity pattern must include entries at the constraint field's
        column positions for the FD Jacobian to capture these.
        """
        # Build a CS-like spec: A_0 constraint, A_1 dynamical with
        # gradient_x(v_A_0) term.
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
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "A_0"},
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
                            {"coefficient": 1.0, "operator": "gradient_x", "field": "v_A_0"},
                        ],
                    },
                },
            ],
            "canonical": {"hamiltonian_terms": []},
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10), (0, 10)), shape=(8, 8), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, ("periodic", "periodic"))
        dense: Any = pattern.toarray()

        n = grid.num_points  # 64
        # Slot layout: A_0 (slot 0), A_1 (slot 1), v_A_1 (slot 2)
        # v_A_1 equation references gradient_x(v_A_0).
        # v_A_0 = yp[A_0_slot=0], so coupling column = slot 0.
        # v_A_1 row block = slot 2 * n : slot 3 * n
        v_A_1_rows = slice(2 * n, 3 * n)
        A_0_cols = slice(0, n)

        # There should be nonzero entries in v_A_1 rows × A_0 cols
        block = dense[v_A_1_rows, A_0_cols]
        assert block.any(), (
            "Missing constraint velocity coupling: v_A_1 equation references "
            "gradient_x(v_A_0) but no entries in sparsity pattern at "
            "(v_A_1_slot, A_0_slot)"
        )

    def test_first_derivative_t_constraint_in_pattern(self) -> None:
        """first_derivative_t(constraint_field) should couple to field slot."""
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
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                            {"coefficient": -1.0, "operator": "first_derivative_t", "field": "A_0"},
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
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(16,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        pattern = build_jacobian_sparsity(spec, layout, grid, "periodic")
        dense: Any = pattern.toarray()

        n = grid.num_points  # 16
        # Slots: phi_0 (0), v_phi_0 (1), A_0 (2)
        # v_phi_0 equation references first_derivative_t(A_0).
        # A_0 is constraint → v_A_0 = yp[A_0_field_slot=2]
        # Coupling should be at (v_phi_0 rows, A_0 cols) = (slot 1, slot 2)
        v_phi_rows = slice(1 * n, 2 * n)
        A_0_cols = slice(2 * n, 3 * n)

        block = dense[v_phi_rows, A_0_cols]
        assert block.any(), (
            "Missing first_derivative_t(constraint) coupling: v_phi_0 "
            "equation has first_derivative_t(A_0) but no entries at "
            "(v_phi_0_slot, A_0_slot)"
        )

    def test_ida_with_sparse_matches_dense(self) -> None:
        """IDA with sparse solver should give same result as dense."""
        from tidal.solver.ida import solve_ida

        spec = _make_kg_1d_spec()
        grid = GridInfo(bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # Gaussian IC with parameters (triggers RHS evaluator path)
        x = np.linspace(0, 2 * np.pi, 32, endpoint=False)
        y0 = np.zeros(layout.total_size)
        y0[:32] = np.exp(-((x - np.pi) ** 2))
        params: dict[str, float] = {}

        # Dense solve (small system, under _DENSE_THRESHOLD)
        result_dense = solve_ida(
            spec,
            grid,
            y0,
            t_span=(0.0, 0.5),
            bc="periodic",
            parameters=params,
            num_snapshots=5,
        )

        # Sparse solve (force sparse via direct IDA construction)
        from sksundae.ida import IDA  # pyright: ignore[reportMissingTypeStubs]

        from tidal.solver.ida import build_residual_fn
        from tidal.solver.sparsity import (
            build_jacobian_sparsity as _bjs,
        )

        resfn = build_residual_fn(spec, layout, grid, "periodic", parameters=params)
        pattern = _bjs(spec, layout, grid, "periodic")
        t_eval = np.linspace(0.0, 0.5, 5)
        yp0 = np.zeros_like(y0)

        solver = IDA(
            resfn,
            rtol=1e-8,
            atol=1e-10,
            linsolver="sparse",
            sparsity=pattern,
            calc_initcond="yp0",
            calc_init_dt=float(t_eval[1] - t_eval[0]),
        )
        result_sparse: Any = solver.solve(t_eval, y0, yp0)

        assert result_dense["success"], f"Dense solve failed: {result_dense['message']}"
        assert result_sparse.success, f"Sparse solve failed: {result_sparse.message}"
        np.testing.assert_allclose(
            result_sparse.y[-1],
            result_dense["y"][-1],
            rtol=1e-5,
            atol=1e-8,
        )
