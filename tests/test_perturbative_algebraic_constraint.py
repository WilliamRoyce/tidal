"""Toy-theory end-to-end test: algebraic constraint promoted by a correction.

This test exercises the full Phase 6A chain (Gap B + Gap C) on the
smallest spec that matters for R̃² PGT-class theories: a two-field
system where one field is dynamical and the other is algebraic at
b₅=0 but has a parameter-dependent kinetic coefficient that formally
promotes it to second order.

Structure
---------
- φ : dynamical (time_order=2), Klein-Gordon plus laplacian.
- h : LHS ``d2_t(h)`` with ``kinetic_coefficient_symbolic = "b5"``. At
  b₅=0 Gap B (`EquationSystem.base_spec(["b5"])`) detects the
  collapsing kinetic and demotes h to the algebraic constraint
  ``1·h - g·φ = 0``.
- Correction on φ: ``b5·h`` (order_in_eps=1) — after Schur
  substitution of h via the Pass 0 recovery, this is a source on the
  φ evolution.

Expectations
------------
1. Pass 0: ``h⁰(x,t) = g·φ⁰(x,t)`` exactly (Schur recovery, no
   approximation).
2. Pass 1: ``φ¹`` is non-zero and drives a non-zero ``h¹`` via the
   same recovery operator — ``h¹ = g·φ¹`` to machine precision.

Together these imply that the v6 loader correctly treats the JSON
spec (which *looks* second-order in h to a naive reader) as the
ghost-free algebraic-constraint system at b₅=0 and layers the b₅
correction on top iteratively.
"""

from __future__ import annotations

import copy

import numpy as np

from tidal.solver.grid import GridInfo
from tidal.solver.perturbative_driver import PerturbativeSolver
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

_G = 0.7  # coupling strength φ → h
_M2 = 1.0  # φ mass²

_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"m2": _M2, "b5": 0.01, "g": _G},
        "perturbation": {"small_parameters": ["b5"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "phi", "index": 0, "is_dynamical": True},
        {"name": "h", "index": 1, "is_dynamical": True},
    ],
    "equations": [
        # phi: ∂²_t phi = ∇²phi - m² phi + b5·g·h  (correction source)
        {
            "field": "phi",
            "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi"},
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi",
                        "coefficient_symbolic": "-m2",
                    },
                    # Correction: b5 · g · h (order-1 source that references
                    # the constraint field).  The _build_source_matrix_k
                    # path Schur-substitutes h via the recovery operator.
                    {
                        "coefficient": 1.0,
                        "operator": "identity",
                        "field": "h",
                        "coefficient_symbolic": "b5*g",
                        "order_in_eps": 1,
                    },
                ],
            },
        },
        # h: LHS d²_t(h) with kinetic b5.  At b5=0 Gap B demotes to
        # algebraic constraint h = g·phi.  Base RHS terms:
        #   identity(h) self-term (Schur constraint form) with
        #     coefficient +1
        #   identity(phi) coupling with coefficient -g (so base eq is
        #     1·h − g·phi = 0, i.e. h = g·phi).
        {
            "field": "h",
            "lhs": {
                "expression": "d2_t(h)",
                "order": {"time": 2, "space": 0},
                "kinetic_coefficient_symbolic": "b5",
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": 1.0,
                        "operator": "identity",
                        "field": "h",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi",
                        "coefficient_symbolic": "-g",
                    },
                ],
            },
        },
    ],
    "coupling": {},
}


def _make_spec() -> EquationSystem:
    return EquationSystem.from_dict(copy.deepcopy(_SPEC))


def _setup() -> tuple[
    EquationSystem,
    PerturbativeSolver,
    GridInfo,
    np.ndarray,
    StateLayout,
]:
    """Build spec + solver + matching IC + the Pass 0 base layout.

    The returned ``layout`` is derived from ``solver.base_spec`` (with
    h demoted), which is what both Pass 0 and Pass 1 outputs live in.
    """
    spec = _make_spec()
    solver = PerturbativeSolver(spec)
    n_grid = 32
    length = 2 * np.pi
    grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
    # Pass 0 / Pass 1 outputs live in the base_spec layout (h demoted
    # to algebraic constraint).
    layout = StateLayout.from_spec(solver.base_spec, grid.num_points)

    # Single-Fourier-mode IC on phi: sin(x).  phi velocity = 0.
    x = np.linspace(0.0, length, n_grid, endpoint=False)
    y0 = np.zeros(layout.num_slots * grid.num_points)
    phi_slot = layout.field_slot_map["phi"]
    y0[phi_slot * n_grid : (phi_slot + 1) * n_grid] = np.sin(x)
    return spec, solver, grid, y0, layout


class TestGapBLHSDemotion:
    """Gap B: JSON has d²_t(h) LHS, driver should treat h as algebraic."""

    def test_base_spec_demotes_h_when_b5_vanishes(self) -> None:
        spec = _make_spec()
        # Driver's base_spec:
        base = spec.base_spec(["b5"])
        h_eq = next(eq for eq in base.equations if eq.field_name == "h")
        assert h_eq.time_derivative_order == 0, (
            "Gap B should have demoted h to algebraic (time_order=0) "
            f"but got {h_eq.time_derivative_order}"
        )
        assert h_eq.kinetic_coefficient_symbolic is None

    def test_perturbative_solver_uses_demoted_base(self) -> None:
        spec = _make_spec()
        solver = PerturbativeSolver(spec)
        h_eq = next(eq for eq in solver.base_spec.equations if eq.field_name == "h")
        assert h_eq.time_derivative_order == 0


class TestPass0SchurRecovery:
    """At Pass 0, h⁰(t) = g · φ⁰(t) exactly — Schur base recovery."""

    def test_pass0_constraint_matches_g_times_phi(self) -> None:
        _spec, solver, grid, y0, layout = _setup()
        res = solver.solve(
            y0,
            grid,
            (0.0, 1.0),
            order=0,  # Pass 0 only
            parameters={"m2": _M2, "b5": 0.01, "g": _G},
            num_snapshots=6,
        )
        n_grid = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]

        # Extract the physical-space trajectories and compare pointwise.
        y = res.total["y"]
        for ti in range(y.shape[0]):
            phi = y[ti, phi_slot * n_grid : (phi_slot + 1) * n_grid]
            h = y[ti, h_slot * n_grid : (h_slot + 1) * n_grid]
            np.testing.assert_allclose(h, _G * phi, atol=1e-10)

    def test_pass0_inconsistent_ic_is_projected_onto_constraint(self) -> None:
        """R5.4 / #283: an IC that violates ``h = g·phi`` at t=0 must
        either be rejected or projected onto the constraint manifold.
        Silent violation (driver runs and h drifts from g·phi) would
        be a bug.

        The constraint IC solver is expected to project: it replaces
        the user-provided h with ``g·phi`` at t=0 so the evolution
        starts on the constraint manifold. This test documents that
        behaviour as the observed contract.
        """
        _spec, solver, grid, y0, layout = _setup()
        n_grid = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]
        # Violate the constraint at t=0: h₀ = 0.3·sin(x), but g·phi₀ = 0.7·sin(x).
        y0 = y0.copy()
        y0[h_slot * n_grid : (h_slot + 1) * n_grid] = 0.3 * np.sin(
            np.linspace(0.0, 2 * np.pi, n_grid, endpoint=False)
        )

        res = solver.solve(
            y0,
            grid,
            (0.0, 1.0),
            order=0,
            parameters={"m2": _M2, "b5": 0.01, "g": _G},
            num_snapshots=6,
        )
        y = res.total["y"]
        # Assert constraint holds at every snapshot — inclusive of t=0,
        # which documents that projection (not silent drift) is what
        # the driver guarantees.
        max_violation = 0.0
        for ti in range(y.shape[0]):
            phi = y[ti, phi_slot * n_grid : (phi_slot + 1) * n_grid]
            h = y[ti, h_slot * n_grid : (h_slot + 1) * n_grid]
            v = float(np.max(np.abs(h - _G * phi)))
            max_violation = max(max_violation, v)
        assert max_violation < 1e-10, (
            f"Inconsistent IC: constraint h = g·phi violated by "
            f"max |h − g·phi| = {max_violation:.3e}. The driver must "
            f"either project onto the constraint manifold OR raise "
            f"loudly; silent violation is a bug. See #283."
        )


class TestPass1ConstraintAndDynamicalAgree:
    """Pass 1 h¹(t) = g · φ¹(t) via Schur recovery (Gap C)."""

    def test_pass1_h_matches_g_times_phi(self) -> None:
        _spec, solver, grid, y0, layout = _setup()
        res = solver.solve(
            y0,
            grid,
            (0.0, 1.0),
            order=1,
            parameters={"m2": _M2, "b5": 0.01, "g": _G},
            num_snapshots=6,
        )

        # Inspect the Pass 1 contribution only (orders[1]["y"]).
        pass1 = res.orders[1]
        n_grid = grid.num_points
        assert pass1["y"].shape[1] == layout.num_slots * n_grid
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]

        # At later times φ¹ must be non-zero (the correction drives it).
        phi1_last = pass1["y"][-1, phi_slot * n_grid : (phi_slot + 1) * n_grid]
        assert np.max(np.abs(phi1_last)) > 1e-8, (
            "Pass 1 φ¹ identically zero — correction source not wired"
        )

        # h¹ must track g·φ¹ to machine precision (Schur recovery).
        for ti in range(pass1["y"].shape[0]):
            phi1 = pass1["y"][ti, phi_slot * n_grid : (phi_slot + 1) * n_grid]
            h1 = pass1["y"][ti, h_slot * n_grid : (h_slot + 1) * n_grid]
            np.testing.assert_allclose(h1, _G * phi1, atol=1e-10)

    def test_combined_total_includes_both_passes(self) -> None:
        """Total = Pass 0 + Pass 1 should differ from Pass 0 alone."""
        _spec, solver, grid, y0, _layout = _setup()
        res_full = solver.solve(
            y0,
            grid,
            (0.0, 1.0),
            order=1,
            parameters={"m2": _M2, "b5": 0.01, "g": _G},
            num_snapshots=5,
        )
        res_base = solver.solve(
            y0,
            grid,
            (0.0, 1.0),
            order=0,
            parameters={"m2": _M2, "b5": 0.01, "g": _G},
            num_snapshots=5,
        )
        diff = np.max(np.abs(res_full.total["y"] - res_base.total["y"]))
        assert diff > 1e-8, (
            "Pass 1 correction is identically zero; "
            "the driver is not actually applying the b5 source"
        )
