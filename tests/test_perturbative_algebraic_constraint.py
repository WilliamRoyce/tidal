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
        behavior as the observed contract.
        """
        _spec, solver, grid, y0, layout = _setup()
        n_grid = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]
        # Violate the constraint at t=0: h₀ = 0.3·sin(x), but g·phi₀ = 0.7·sin(x).
        y0 = y0.copy()
        y0[h_slot * n_grid : (h_slot + 1) * n_grid] = 0.3 * np.sin(
            np.linspace(0.0, 2 * np.pi, n_grid, endpoint=False),
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
            f"max |h - g*phi| = {max_violation:.3e}. The driver must "
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

        # R8 / #290: the recovered h¹ now includes the LHS-feedback
        # contribution. Analytical light-mode eigenvector:
        #   h_total = g·phi·(1 - ε·ω²) + O(ε²)
        # so pass1_h = g·pass1_phi + ε·g·(-ω²)·phi⁰, with ε = b5.
        omega2 = _M2 + 1.0  # k=1, so ω² = m² + 1 = 2
        b5 = 0.01
        pass0_y = res.orders[0]["y"]
        for ti in range(pass1["y"].shape[0]):
            phi1 = pass1["y"][ti, phi_slot * n_grid : (phi_slot + 1) * n_grid]
            h1 = pass1["y"][ti, h_slot * n_grid : (h_slot + 1) * n_grid]
            phi0 = pass0_y[ti, phi_slot * n_grid : (phi_slot + 1) * n_grid]
            expected_lhs_feedback = -b5 * omega2 * _G * phi0
            actual_lhs_feedback = h1 - _G * phi1
            np.testing.assert_allclose(
                actual_lhs_feedback,
                expected_lhs_feedback,
                atol=1e-10,
            )

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


class TestAugmentedRecoveryMatchesMpmath:
    """R8.6 / #290: validate augmented Schur recovery against an mpmath
    50-dps eigenvector of the 2x2 dispersion on the synthetic spec.

    The full EOM (before expansion):

        b5 * d²_t(h) = h - g*phi
        d²_t(phi) = (laplacian - m²)*phi + b5*g*h

    At k=1 with plane-wave ansatz, the light-mode eigenvector has
    ``h/phi = g/(1 + b5*omega²)``. Expanding in b5 gives the
    analytical O(ε¹) correction ``h = g*phi*(1 - b5*omega²) + O(b5²)``.

    We compute this ratio at arbitrary precision via mpmath and compare
    against the driver's Pass 0 + Pass 1 output.
    """

    def test_light_mode_ratio_matches_analytical_at_eps_1e_3(self) -> None:
        _spec, solver, grid, y0, layout = _setup()
        params = {"m2": _M2, "b5": 1e-3, "g": _G}
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters=params,
            num_snapshots=11,
        )
        n = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]
        phi_tot = res.total["y"][:, phi_slot * n : (phi_slot + 1) * n]
        h_tot = res.total["y"][:, h_slot * n : (h_slot + 1) * n]

        import mpmath  # type: ignore[import-untyped]

        mpmath.mp.dps = 50
        b5 = mpmath.mpf(params["b5"])
        m2 = mpmath.mpf(params["m2"])
        g = mpmath.mpf(params["g"])
        k2 = mpmath.mpf(1)
        # Light-mode omega² from the 2x2 dispersion det: solve
        #   b5*omega^4 - (1 + a*b5)*omega^2 + (a - b5*g^2) = 0  with a = m² + k².
        # Take the smaller root (light mode).
        a = m2 + k2  # type: ignore[reportUnknownVariableType]
        # Light-mode dispersion: ω² = a + ε·g² + O(ε²) from
        # (a - ω²)·(1 + ε·ω²) + ε·g² = 0 (from det of 2x2 block).
        # Exact eigenvector ratio: h/phi = g / (1 + ε·ω²).
        # Driver Pass 0 + Pass 1 truncates at O(ε¹):
        #   h/phi ≈ g·(1 − ε·ω²_0) = g·(1 − ε·a) + O(ε²)
        ratio_truncated = g * (1 - b5 * a)  # type: ignore[reportUnknownVariableType]
        ratio_analytical = float(ratio_truncated)  # type: ignore[reportUnknownArgumentType]

        # Driver output: sample ratio h_total / phi_total at the
        # snapshot / grid point with largest |phi| (away from nodes).
        max_idx = np.argmax(np.abs(phi_tot.ravel()))
        ti_sample, xi_sample = np.unravel_index(max_idx, phi_tot.shape)
        ratio_driver = h_tot[ti_sample, xi_sample] / phi_tot[ti_sample, xi_sample]
        rel_err = abs(ratio_driver - ratio_analytical) / abs(ratio_analytical)
        assert rel_err < 1e-6, (
            f"driver h/phi ratio {ratio_driver:.8f}, mpmath-truncated "
            f"{ratio_analytical:.8f}, rel_err {rel_err:.3e} exceeds 1e-6. "
            f"R8 augmented recovery may be off. See #290."
        )

    def test_light_mode_ratio_matches_exact_at_eps_1e_3(self) -> None:
        """At ε=1e-3, the EXACT light-mode eigenvector g/(1+ε·ω²)
        differs from the O(ε¹) truncation by O(ε²) ≈ 4e-6. The driver's
        output (truncated at O(ε)) should match the exact mpmath
        eigenvector to within that O(ε²) tolerance.
        """
        _spec, solver, grid, y0, layout = _setup()
        params = {"m2": _M2, "b5": 1e-3, "g": _G}
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters=params,
            num_snapshots=11,
        )
        n = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]
        phi_tot = res.total["y"][:, phi_slot * n : (phi_slot + 1) * n]
        h_tot = res.total["y"][:, h_slot * n : (h_slot + 1) * n]

        import mpmath  # type: ignore[import-untyped]

        mpmath.mp.dps = 50
        b5 = mpmath.mpf(params["b5"])
        g = mpmath.mpf(params["g"])
        a = mpmath.mpf(params["m2"]) + mpmath.mpf(1)  # type: ignore[reportUnknownVariableType]
        # Light-mode ω² from dispersion: (a - ω²)(1 + ε·ω²) + ε·g² = 0
        # Solve exactly via quadratic: ε·ω⁴ + (1 - ε·a)·ω² - (a + ε·g²) = 0
        A = b5
        B = 1 - b5 * a  # type: ignore[reportUnknownVariableType]
        C = -(a + b5 * g * g)  # type: ignore[reportUnknownVariableType]
        disc = mpmath.sqrt(B * B - 4 * A * C)  # type: ignore[reportUnknownVariableType]
        # Light mode = smaller positive root.
        roots = [(-B + disc) / (2 * A), (-B - disc) / (2 * A)]  # type: ignore[reportUnknownVariableType]
        omega2_light_mp = min((r for r in roots if r > 0), key=abs)  # type: ignore[reportUnknownVariableType]
        ratio_exact = g / (1 + b5 * omega2_light_mp)  # type: ignore[reportUnknownVariableType]
        ratio_analytical = float(ratio_exact)  # type: ignore[reportUnknownArgumentType]

        max_idx = np.argmax(np.abs(phi_tot.ravel()))
        ti_sample, xi_sample = np.unravel_index(max_idx, phi_tot.shape)
        ratio_driver = h_tot[ti_sample, xi_sample] / phi_tot[ti_sample, xi_sample]
        rel_err = abs(ratio_driver - ratio_analytical) / abs(ratio_analytical)
        # O(ε²) truncation error ~ (ε·ω²)² = 4e-6 at ε=1e-3. Allow 1e-5.
        assert rel_err < 1e-5, (
            f"driver h/phi {ratio_driver:.10f}, mpmath exact "
            f"{ratio_analytical:.10f}, rel_err {rel_err:.3e}. Expected "
            f"O(ε²) agreement (~4e-6). See #290."
        )

    def test_light_mode_ratio_matches_truncation_at_eps_1e_2(self) -> None:
        """Same check at ε=1e-2. Tolerance looser (1e-4) because O(ε²)
        truncation error scales as (b5·ω²)² which is ~4e-4 at this ε.
        """
        _spec, solver, grid, y0, layout = _setup()
        params = {"m2": _M2, "b5": 1e-2, "g": _G}
        res = solver.solve(
            y0,
            grid,
            (0.0, 0.5),
            order=1,
            parameters=params,
            num_snapshots=11,
        )
        n = grid.num_points
        phi_slot = layout.field_slot_map["phi"]
        h_slot = layout.field_slot_map["h"]
        phi_tot = res.total["y"][:, phi_slot * n : (phi_slot + 1) * n]
        h_tot = res.total["y"][:, h_slot * n : (h_slot + 1) * n]

        import mpmath  # type: ignore[import-untyped]

        mpmath.mp.dps = 50
        b5 = mpmath.mpf(params["b5"])
        m2 = mpmath.mpf(params["m2"])
        g = mpmath.mpf(params["g"])
        k2 = mpmath.mpf(1)
        a = m2 + k2  # type: ignore[reportUnknownVariableType]
        # Pass-0 + Pass-1 corresponds to h/phi ≈ g·(1 − b5·ω²_0) at
        # O(ε¹), with ω²_0 = a (leading-order). Higher-order ω²
        # corrections kick in at O(ε²) which are below our tolerance.
        ratio_mp_truncated = g * (1 - b5 * a)  # type: ignore[reportUnknownVariableType]
        ratio_analytical = float(ratio_mp_truncated)  # type: ignore[reportUnknownArgumentType]

        max_idx = np.argmax(np.abs(phi_tot.ravel()))
        ti_sample, xi_sample = np.unravel_index(max_idx, phi_tot.shape)
        ratio_driver = h_tot[ti_sample, xi_sample] / phi_tot[ti_sample, xi_sample]
        rel_err = abs(ratio_driver - ratio_analytical) / abs(ratio_analytical)
        assert rel_err < 1e-4, (
            f"driver h/phi ratio {ratio_driver:.8f}, mpmath truncated "
            f"{ratio_analytical:.8f}, rel_err {rel_err:.3e} exceeds 1e-4."
        )
