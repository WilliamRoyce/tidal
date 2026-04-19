"""Stage 4 end-to-end tests for the Pass 1 Duhamel modal solver.

Validates :func:`tidal.solver.modal.solve_modal_pass1` against closed-form
analytical solutions of driven linear oscillators: the base Pass 0
system (unperturbed Klein-Gordon) is solved by the existing modal
eigendecomposition, and a small correction term — here an additional
mass shift ``ε·φ`` — is solved by the Duhamel kernel. The sum of the
two passes is compared to the exact full solution at the same ε.

Reference: the iterative expansion
    q(t) = q⁽⁰⁾(t) + ε·q⁽¹⁾(t) + O(ε²)
with q⁽⁰⁾ driving the Pass 1 source. See Parker & Simon 1993
(literature/gr-qc_9211002/) and the v6 implementation plan.
"""

from __future__ import annotations

import copy
from typing import Any, cast

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.modal import solve_modal, solve_modal_pass1
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Shared KG base spec (same shape as test_solver_modal._KG_1D_SPEC)
# ---------------------------------------------------------------------------

_KG_BASE: dict[str, object] = {
    "metadata": {"source": "inline-test", "parameters": {"m2": 1.0, "eps": 0.1}},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    # Base mass: -m2 · φ  (order 0)
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    # Base laplacian: +∇²φ  (order 0)
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                    # Correction: -eps · φ  (order 1 — mass shift)
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        }
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2 - eps"]]},
}


def _make_spec(data: dict[str, Any]) -> EquationSystem:
    return EquationSystem.from_dict(copy.deepcopy(data))


class TestSolveModalPass1AnalyticalMatch:
    """Driven Klein-Gordon: analytical answer via iterative expansion."""

    @pytest.fixture
    def setup(self) -> dict[str, Any]:
        spec = _make_spec(_KG_BASE)
        n = 64
        length = 2 * np.pi
        grid = GridInfo(shape=(n,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        m2 = 1.0
        eps = 0.05
        k_mode = 1.0
        x = np.linspace(0.0, length, n, endpoint=False)

        # IC: phi = sin(x), v = 0 — single Fourier mode
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n] = np.sin(k_mode * x)
        return {
            "spec": spec,
            "grid": grid,
            "layout": layout,
            "y0": y0,
            "m2": m2,
            "eps": eps,
            "k_mode": k_mode,
            "x": x,
            "n": n,
        }

    def test_pass1_combined_matches_exact_to_eps2(self, setup: dict[str, Any]) -> None:
        """q⁽⁰⁾ + ε·q⁽¹⁾ matches exact full-ε solution to O(ε²)."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]
        m2 = setup["m2"]
        eps = setup["eps"]
        setup["k_mode"]
        setup["x"]
        n = setup["n"]
        t_end = 2.0
        num_snap = 21

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)

        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters={"m2": m2, "eps": eps},
                num_snapshots=num_snap,
                return_eigendata=True,
            ),
        )
        t_eval = pass0["t"]
        eigendata = pass0["eigendata"]

        pass1 = solve_modal_pass1(
            eigendata,
            correction_spec,
            grid,
            t_eval,
            parameters={"m2": m2, "eps": eps},
        )

        # Combined approximation: q_total = q⁽⁰⁾ + q⁽¹⁾
        # (ε is already baked into the correction coefficients by the
        # coefficient evaluator, so no extra ε multiplier is needed.)
        q_total = pass0["y"] + pass1["y"]

        # Exact: evolve the FULL spec (both order=0 and order=1 terms) to
        # ω_full = sqrt(k² + m² + ε). Use modal solver as ground truth.
        full = cast(
            "dict[str, Any]",
            solve_modal(
                spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters={"m2": m2, "eps": eps},
                num_snapshots=num_snap,
            ),
        )

        # The iterative approximation error is O(ε²). With ε=0.05 and
        # ω·t ~ (sqrt(1+1+0.05))·2 ≈ 2.9, the expected truncation error
        # is roughly ε² · (ω·t) · max|q| ≈ 2.5e-3 · 3 ≈ 7e-3. Use 1e-2.
        final_err = float(np.max(np.abs(q_total[-1, :n] - full["y"][-1, :n])))
        assert final_err < 1e-2, (
            f"Pass 0 + Pass 1 truncation error {final_err:.2e} exceeds O(ε²) bound"
        )

        # Sanity: the pure Pass 0 error is the omitted O(ε) mass-shift term,
        # much larger than the Pass 1-corrected error. Verify that the
        # correction substantially improves accuracy.
        pass0_only_err = float(np.max(np.abs(pass0["y"][-1, :n] - full["y"][-1, :n])))
        assert final_err < 0.25 * pass0_only_err, (
            f"Pass 1 failed to improve accuracy: "
            f"pass0_err={pass0_only_err:.3e}, combined_err={final_err:.3e}"
        )

    def test_pass1_zero_source_returns_zero(self, setup: dict[str, Any]) -> None:
        """With no order-1 terms, Pass 1 must return the zero solution."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]

        base_spec = spec.filter_by_order(0)
        empty_spec = spec.filter_by_order(99)  # no such order → all-empty rhs
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 1.0),
                parameters={"m2": 1.0, "eps": 0.05},
                num_snapshots=5,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            empty_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.05},
        )
        assert np.allclose(pass1["y"], 0.0, atol=1e-14)

    def test_pass1_initial_condition_is_zero(self, setup: dict[str, Any]) -> None:
        """Pass 1 by construction satisfies q⁽¹⁾(0) = 0."""
        spec = setup["spec"]
        grid = setup["grid"]
        y0 = setup["y0"]

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 1.0),
                parameters={"m2": 1.0, "eps": 0.05},
                num_snapshots=3,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            correction_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.05},
        )
        # At t = t_eval[0], the Pass 1 correction is zero (IC).
        assert np.max(np.abs(pass1["y"][0])) < 1e-14


class TestPass1NearDegeneracy:
    """Phase 6C.3: Pass 1 kernel stability in the near-degenerate regime.

    The closed-form Duhamel kernel G(λ, μ; t) has a removable
    singularity at μ = λ (degenerate eigenvalues).  Near-degenerate
    modes arise naturally in coupled systems where the base operator
    eigenvalues cluster — the correction source then couples modes
    whose frequencies nearly match.  The modal solver handles this
    via a Taylor branch inside _duhamel_kernel; Phase 4 already
    covers the scalar kernel sweep against mpmath to 1e-13.

    This test exercises the stability of the FULL pipeline
    (_build_source_matrix_k → _evolve_duhamel_per_mode → Pass 1
    output) when the source frequency matches the base eigenfrequency.
    The test uses a scalar KG system where the single-mode base
    frequency is ω₀ = √(m² + k²) and the correction re-introduces
    the same identity source at the same mode — effectively a
    resonant (μ = λ) case.  The test asserts the output remains
    finite and the Taylor branch did not produce NaN/Inf.
    """

    def test_resonant_source_stays_finite(self) -> None:
        spec = _make_spec(_KG_BASE)
        n_grid = 32
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        # IC: phi = sin(x) — excites a single Fourier mode at k=1.
        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n_grid] = np.sin(x)

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)

        # Large t to stress the long-time behaviour.
        t_end = 5.0
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters={"m2": 1.0, "eps": 0.1},
                num_snapshots=11,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            correction_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.1},
        )

        assert np.all(np.isfinite(pass1["y"])), (
            "Pass 1 resonant output contains NaN/Inf — Taylor branch of "
            "the Duhamel kernel failed on a near-degenerate source"
        )
        # The resonant correction grows ∝ t (secular), so amplitude
        # should be bounded by (eps · ω² · t) × initial_amplitude.
        max_amp = float(np.max(np.abs(pass1["y"])))
        assert max_amp < 10.0, (
            f"Resonant Pass 1 amplitude {max_amp:.3e} implausibly large; "
            f"secular growth formula predicts ~eps*omega^2*t ≈ 0.1*2*5 = 1"
        )

    def test_very_large_t_remains_finite(self) -> None:
        """At very large t, the secular growth eventually violates
        validity but the numerics should still produce finite output,
        letting the validity monitor flag the regime breakdown.
        """
        spec = _make_spec(_KG_BASE)
        n_grid = 16
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n_grid] = np.sin(x)

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 100.0),
                parameters={"m2": 1.0, "eps": 0.01},
                num_snapshots=3,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            correction_spec,
            grid,
            pass0["t"],
            parameters={"m2": 1.0, "eps": 0.01},
        )
        assert np.all(np.isfinite(pass1["y"]))


# ---------------------------------------------------------------------------
# Two-field spec whose base theory has two independent Pass 0 blocks
# (phi and chi decoupled at ε=0) and whose correction links them.
# The Pass 1 Duhamel path does not yet support cross-block source
# coupling; the solver must refuse cleanly with NotImplementedError.
# ---------------------------------------------------------------------------

_CROSS_BLOCK_SPEC: dict[str, object] = {
    "metadata": {
        "source": "inline-test",
        "parameters": {"mPhi2": 1.0, "mChi2": 2.0, "eps": 0.05},
        "perturbation": {"small_parameters": ["eps"], "order": 1},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "phi_0", "index": 0, "is_dynamical": True},
        {"name": "chi_0", "index": 1, "is_dynamical": True},
    ],
    "equations": [
        # phi_0: Klein-Gordon with mass mPhi2 + order-1 coupling to chi_0.
        {
            "field": "phi_0",
            "lhs": {
                "expression": "d2_t(phi_0)",
                "order": {"time": 2, "space": 0},
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-mPhi2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                    # Order-1 cross-coupling: eps·chi_0 on phi_0's RHS.
                    # This links the previously-independent phi and chi
                    # sectors via the correction.
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        },
        {
            "field": "chi_0",
            "lhs": {
                "expression": "d2_t(chi_0)",
                "order": {"time": 2, "space": 0},
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-mChi2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                ],
            },
        },
    ],
    "coupling": {},
}


class TestCrossBlockCouplingRaises:
    """Phase 6C.4: Duhamel evolution refuses cross-block source coupling."""

    def test_cross_block_correction_raises(self) -> None:
        """A correction linking phi (block A) to chi (block B) must raise.

        The Duhamel kernel operates per eigenblock; Pass 0
        independently eigendecomposed phi and chi. A correction term
        whose source references chi from inside phi's equation would
        require a joint eigenbasis, which the current implementation
        does not build. It should refuse with NotImplementedError and
        a clear diagnostic message so users (or future maintainers)
        know this is a documented limitation, not a silent failure.
        """
        spec = _make_spec(_CROSS_BLOCK_SPEC)
        n_grid = 16
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n_grid] = np.sin(x)

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)

        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 0.5),
                parameters={"mPhi2": 1.0, "mChi2": 2.0, "eps": 0.05},
                num_snapshots=3,
                return_eigendata=True,
            ),
        )
        with pytest.raises(NotImplementedError, match="[Cc]ross-block"):
            solve_modal_pass1(
                pass0["eigendata"],
                correction_spec,
                grid,
                pass0["t"],
                parameters={"mPhi2": 1.0, "mChi2": 2.0, "eps": 0.05},
            )

    def test_cross_block_guard_tolerates_schur_tail_noise(self) -> None:
        """v6 R1.4 / #275: the guard uses a scale-relative threshold.

        Absolute 1e-14 is below double-precision roundoff of compound
        operations like ``coeff * mult[m] * recovery_matrix[m, c_idx, :]``.
        Schur-recovering constraint fields on ill-conditioned S_cc can
        legitimately produce O(1e-12) tail noise in off-block columns;
        that should NOT raise. Only coupling magnitudes meaningfully
        above roundoff (relative to max|M_src|) are genuine.

        Direct unit test: craft an M_src matrix where the off-block
        tail noise is below ``max|M_src| * 1e-10`` and verify the
        guard stays silent.
        """
        from tidal.solver.modal import (
            _evolve_duhamel_per_mode,
        )

        # Base spec: two independent KG fields → two eigenblocks.
        base_data = copy.deepcopy(_CROSS_BLOCK_SPEC)
        # Strip the order-1 coupling so the base spec is actually 2
        # independent blocks.
        phi_eq = base_data["equations"][0]  # type: ignore[index]
        phi_eq["rhs"]["terms"] = [  # type: ignore[index]
            t for t in phi_eq["rhs"]["terms"] if t.get("order_in_eps", 0) == 0
        ]
        base_spec = _make_spec(base_data)
        n_grid = 16
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(base_spec, grid.num_points)

        x = np.linspace(0.0, length, n_grid, endpoint=False)
        # Non-zero IC on BOTH fields so both blocks survive the zero-IC
        # pruning in solve_modal (otherwise chi_0's block is dropped).
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[:n_grid] = np.sin(x)  # phi_0
        y0[2 * n_grid : 3 * n_grid] = 0.5 * np.sin(x)  # chi_0

        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, 0.5),
                parameters={"mPhi2": 1.0, "mChi2": 2.0},
                num_snapshots=3,
                return_eigendata=True,
            ),
        )

        # Hand-craft M_src_k with legitimate in-block source + tiny
        # off-block tail noise from a hypothetical Schur recovery.
        eigendata = pass0["eigendata"]
        n_slots = layout.num_slots
        rfft_last = grid.shape[-1] // 2 + 1
        n_modes = int(np.prod([*grid.shape[:-1], rfft_last]))
        m_src = np.zeros((n_modes, n_slots, n_slots), dtype=np.complex128)
        blocks = eigendata["blocks"]
        assert len(blocks) >= 2, (
            f"Need two independent blocks for this test; got {len(blocks)}"
        )
        b0_idx = int(blocks[0]["slot_indices"][0])
        b1_idx = int(blocks[1]["slot_indices"][0])
        m_src[:, b0_idx, b0_idx] = 0.1  # in-block source, magnitude 0.1
        m_src[:, b0_idx, b1_idx] = 1e-13  # off-block tail noise
        # max_scale = 0.1, atol = max(1e-14, 0.1 * 1e-10) = 1e-11.
        # cross = 1e-13 < 1e-11 → guard must NOT raise.

        t_eval = pass0["t"]
        _ = _evolve_duhamel_per_mode(eigendata, m_src, t_eval, layout, grid)

        # Regression: bumping the tail to well above the relative
        # threshold re-triggers the guard.
        m_src[:, b0_idx, b1_idx] = 1e-3  # ratio 1e-2 > 1e-10 → must raise
        with pytest.raises(NotImplementedError, match="[Cc]ross-block"):
            _evolve_duhamel_per_mode(eigendata, m_src, t_eval, layout, grid)
