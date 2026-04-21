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

# ruff: noqa: RUF002, RUF012, RUF043 — math symbols in docstrings; ClassVar in test fixtures.

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


class TestPass1Resonance:
    """Phase 6C.3 (renamed R5.3 / #282): Pass 1 kernel stability when the
    correction source resonates with the base eigenfrequency.

    Pre-R5.3 this class was mis-labelled ``TestPass1NearDegeneracy``,
    but it actually tests *resonance* (correction sources a mode that
    matches the base frequency), not *near-degeneracy* (two base
    eigenvalues close to each other). They are unrelated in general.
    The true Taylor-branch crossover is tested by
    :class:`TestPass1NearDegeneracy` below (R5.3).

    The test uses a scalar KG system where the single-mode base
    frequency is ω₀ = √(m² + k²) and the correction re-introduces
    the same identity source at the same mode — effectively a
    resonant (μ = λ) case. The Duhamel kernel ``G(λ, λ, t) = t·exp(λt)``
    is trivially exact in both branches; this exercises the full
    pipeline path but does NOT stress the Taylor-branch math.
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
            t
            for t in phi_eq["rhs"]["terms"]  # type: ignore[reportUnknownVariableType]
            if t.get("order_in_eps", 0) == 0
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
        # #293: _evolve_duhamel_per_mode now takes a dict keyed by
        # operator time-derivative order. Identity sources have order 0.
        m_src_by_order = {0: m_src}
        _ = _evolve_duhamel_per_mode(eigendata, m_src_by_order, t_eval, layout, grid)

        # Regression: bumping the tail to well above the relative
        # threshold re-triggers the guard.
        m_src[:, b0_idx, b1_idx] = 1e-3  # ratio 1e-2 > 1e-10 → must raise
        with pytest.raises(NotImplementedError, match="[Cc]ross-block"):
            _evolve_duhamel_per_mode(eigendata, m_src_by_order, t_eval, layout, grid)


class TestPass1NearDegeneracy:
    """R5.3 / #282: Pass 1 Duhamel kernel in the Taylor-branch regime.

    ``_duhamel_kernel(λ, μ; t)`` has two numerical branches:

    - Direct:  ``G = exp(λt) · (exp(z) - 1) / z``   with ``z = (μ-λ)·t``
    - Taylor:  series expansion around ``z = 0`` to avoid catastrophic
      cancellation when ``|z|`` is small.

    The crossover window is ``|z| ∈ [1e-8, 1e-5]``. The existing
    scalar-level sweep in ``test_modal_duhamel_degeneracy.py`` validates
    the kernel math to 1e-13 in that window. This class adds the
    missing **pipeline-level** check: construct crafted eigendata with
    a block whose two eigenvalues are close (Δλ = 2e-11), sweep
    ``t`` so ``|z|`` traverses the crossover, feed the inputs through
    ``_evolve_duhamel_per_mode`` (the public Duhamel evolver), and
    compare against a scalar mpmath reference.

    Pre-R5.3 ``TestPass1NearDegeneracy`` tested resonance (μ = λ
    exactly), not near-degeneracy. Renamed to :class:`TestPass1Resonance`.
    """

    @staticmethod
    def _craft_eigendata(
        lam_i: complex, lam_j: complex, alpha: np.ndarray, grid: GridInfo
    ) -> tuple[dict[str, Any], StateLayout]:
        """Build a minimal 2-slot eigendata + layout for direct use in
        _evolve_duhamel_per_mode. The block has diagonal D = diag(λ_i,
        λ_j), V = V⁻¹ = I, so the test drives the kernel math in
        isolation from the spec loader / companion-matrix machinery.
        """
        from tidal.solver.state import SlotInfo

        bs = 2
        rfft_last = grid.shape[-1] // 2 + 1
        n_modes = int(np.prod([*grid.shape[:-1], rfft_last]))
        v_mat = np.broadcast_to(
            np.eye(bs, dtype=np.complex128), (n_modes, bs, bs)
        ).copy()
        v_inv = v_mat.copy()
        d_diag = np.broadcast_to(
            np.array([lam_i, lam_j], dtype=np.complex128), (n_modes, bs)
        ).copy()
        alpha_full = np.broadcast_to(alpha, (n_modes, bs)).astype(np.complex128).copy()
        block = {
            "slot_indices": (0, 1),
            "V": v_mat,
            "V_inv": v_inv,
            "D_diag": d_diag,
            "alpha": alpha_full,
        }
        slot_a = SlotInfo(name="a", field_name="a", kind="field", time_order=2)
        slot_b = SlotInfo(name="b", field_name="b", kind="field", time_order=2)
        layout = StateLayout(
            slots=(slot_a, slot_b),
            num_points=grid.num_points,
            field_slot_map={"a": 0, "b": 1},
            velocity_slot_map={},
            dynamical_fields=("a", "b"),
        )
        eigendata = {"blocks": [block], "state_layout": layout}
        return eigendata, layout

    @pytest.mark.parametrize(
        ("lam_diff", "t_end"),
        [
            (2e-11, 1.0),  # |z| = 2e-11 → Taylor branch
            (2e-7, 50.0),  # |z| = 1e-5  → Taylor/direct crossover
            (2e-5, 500.0),  # |z| = 1e-2  → direct branch
            (2e-3, 500.0),  # |z| = 1.0   → direct branch
        ],
    )
    def test_pipeline_matches_mpmath_across_crossover(
        self, lam_diff: float, t_end: float
    ) -> None:
        """Sweep |z| = |Δλ|·t across the Taylor-branch crossover and
        assert ``_evolve_duhamel_per_mode`` agrees with an mpmath
        scalar reference to 1e-10.

        The math: for a single block with diagonal eigenvalues
        ``diag(λ_i, λ_j)``, IC amplitudes ``α = (1, 0)``, and source
        matrix ``β = [[0, 1], [0, 0]]`` (source into slot 0 from
        slot 1), the Pass 1 amplitude at slot 0 is
        ``z_0(t) = β_{01} · α_1 · G(λ_i, λ_j; t)``. With α_1 = 0 the
        result is zero — useless. Flip: α = (0, 1), β = [[0, 1],[0, 0]]
        → z_0(t) = G(λ_i, λ_j; t). That's the scalar we compare
        against mpmath's direct evaluation of the kernel.
        """
        import mpmath  # type: ignore[import-untyped]

        from tidal.solver.modal import (
            _evolve_duhamel_per_mode,
        )

        # Closely-spaced real-part-zero eigenvalues (oscillatory).
        lam_i = 1.0j
        lam_j = lam_i + lam_diff  # purely real offset, keeps |Δλ| = lam_diff

        # α = (0, 1) so only slot 1 is excited at t=0.
        alpha0 = np.array([0.0, 1.0], dtype=np.complex128)
        # Need at least shape=(2,) so rfft has ≥1 mode.
        grid = GridInfo(shape=(2,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
        eigendata, layout = self._craft_eigendata(lam_i, lam_j, alpha0, grid)

        # M_src such that β = V⁻¹·M·V = [[0,1],[0,0]] (slot 0 <- slot 1).
        # With V = I, β = M. Set M[0,1] = 1.
        rfft_last = grid.shape[-1] // 2 + 1
        n_modes = int(np.prod([*grid.shape[:-1], rfft_last]))
        m_src = np.zeros((n_modes, 2, 2), dtype=np.complex128)
        m_src[:, 0, 1] = 1.0

        t_eval = np.array([0.0, t_end], dtype=np.float64)
        # #293: _evolve_duhamel_per_mode takes a dict keyed by operator
        # time-derivative order. This test crafts a bare identity-form
        # M_src representing a source at order 0.
        _, _y_phys, y_hat_snap = _evolve_duhamel_per_mode(
            eigendata, {0: m_src}, t_eval, layout, grid
        )

        # Pipeline output at snapshot 1, slot 0, mode 0:
        # β_{01} · α_1 · G(λ_i, λ_j; t_end) → compare with scalar mpmath.
        pipeline_z = complex(y_hat_snap[1, 0, 0])

        # mpmath reference: G(λ_i, λ_j; t) = (exp(μt) - exp(λt)) / (μ - λ)
        # (equivalent form with t factored out: exp(λt)·t·expm1(z)/z,
        # matching the direct branch in _duhamel_kernel after the
        # off-by-t fix committed in Stage 4).
        mp = mpmath.mp
        mp.dps = 50
        z = (lam_j - lam_i) * t_end
        # Primary reference: direct formula, arbitrary-precision.
        num = mpmath.exp(lam_j * t_end) - mpmath.exp(lam_i * t_end)  # type: ignore[reportUnknownVariableType]
        den = mpmath.mpc(lam_j - lam_i)  # type: ignore[reportUnknownVariableType]
        G_ref_mp = num / den  # type: ignore[reportUnknownVariableType]
        G_ref_complex = complex(float(G_ref_mp.real), float(G_ref_mp.imag))  # type: ignore[reportUnknownArgumentType]
        _ = z  # keep z in scope for the diagnostic message below

        rel_err = abs(pipeline_z - G_ref_complex) / (abs(G_ref_complex) + 1e-30)
        assert rel_err < 1e-10, (
            f"lam_diff={lam_diff:.1e}, t_end={t_end:.1e}, "
            f"|z|={abs(z):.3e}: pipeline={pipeline_z!r} vs "
            f"mpmath={G_ref_complex!r}, rel_err={rel_err:.3e}. "
            f"Taylor-branch crossover broken (#282)."
        )


class TestPass1TimeDerivativeTargetingDynamical:
    """#293: correction operators with ``time_order > 0`` targeting a
    dynamical field must carry the ``λⁿ`` eigenvalue factor when
    projected into the Pass 1 source. The pre-fix
    ``_build_source_matrix_k`` used only the spatial multiplier and
    silently dropped ``time_order``, producing Pass 1 output that was
    wrong by a factor of ``1 / λⁿ`` for any ``d^n_t(dyn_field)`` term.

    Not triggered by any shipped theory (every O(ε¹) b5-coupled
    correction targets a constraint field), but a latent bug for any
    future theory.
    """

    _SPEC_WITH_DT_CORRECTION: dict[str, object] = {
        "metadata": {
            "source": "inline-test-293",
            "parameters": {"m2": 1.0, "gamma": 0.05},
            "perturbation": {"small_parameters": ["gamma"], "order": 1},
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
                "lhs": {
                    "expression": "d2_t(phi_0)",
                    "order": {"time": 2, "space": 0},
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        # Base Klein-Gordon: d²_t(phi) = ∇²phi − m²phi.
                        {
                            "coefficient": -1.0,
                            "operator": "identity",
                            "field": "phi_0",
                            "coefficient_symbolic": "-m2",
                        },
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                        },
                        # Order-1 damping correction: −γ · ∂_t(phi). The
                        # operator ``first_derivative_t`` has
                        # ``_OperatorDecomp.time_order == 1``, so its
                        # contribution must scale by ``λ`` in the
                        # Pass 1 source matrix.
                        {
                            "coefficient": -1.0,
                            "operator": "first_derivative_t",
                            "field": "phi_0",
                            "coefficient_symbolic": "-gamma",
                            "order_in_eps": 1,
                        },
                    ],
                },
            }
        ],
        "coupling": {},
    }

    def test_first_derivative_t_on_dynamical_target_scales_by_lambda(
        self,
    ) -> None:
        """Pass 0 + Pass 1 of a damped KG matches the analytical
        first-order-in-γ eigenvalue shift.

        Base: d²_t(phi) + m²phi − ∇²phi = 0. Eigenvalues per mode:
        ``λ = ±i·ω_0`` with ``ω_0² = m² + k²``.

        Correction: add ``-γ · ∂_t(phi)``. Full dispersion at O(γ¹):
        ``λ = ±i·ω_0 − γ/2 + O(γ²)`` (damped oscillator). So the
        Pass 0 + Pass 1 trajectory at small γ should decay like
        ``exp(-γ·t/2) · cos(ω_0·t)`` up to O(γ²).

        Without the λ factor (pre-fix), Pass 1 source would be the
        identity operator scaled by γ, yielding ``λ = ±i·ω_0 · (1 -
        γ/ω_0)`` — a DIFFERENT correction at O(γ).
        """
        import copy as _copy

        spec = _make_spec(_copy.deepcopy(self._SPEC_WITH_DT_CORRECTION))
        n_grid = 32
        length = 2 * np.pi
        grid = GridInfo(shape=(n_grid,), bounds=((0.0, length),), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = np.linspace(0.0, length, n_grid, endpoint=False)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        # IC: phi(x, 0) = sin(x), v_phi = 0 → k=1 mode only.
        y0[:n_grid] = np.sin(x)

        base_spec = spec.filter_by_order(0)
        correction_spec = spec.filter_by_order(1)

        m2 = 1.0
        gamma = 1e-3  # small: stay well within O(γ²) truncation error
        t_end = 1.0
        params = {"m2": m2, "gamma": gamma}

        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                base_spec,
                grid,
                y0,
                t_span=(0.0, t_end),
                parameters=params,
                num_snapshots=11,
                return_eigendata=True,
            ),
        )
        pass1 = solve_modal_pass1(
            pass0["eigendata"],
            correction_spec,
            grid,
            pass0["t"],
            parameters=params,
        )
        _ = pass0["y"] + pass1["y"]

        omega0 = float(np.sqrt(m2 + 1.0))  # k=1
        t_vals = pass0["t"]
        # Analytical: damped oscillator d²Φ + ω₀²Φ = -γ·dΦ. With
        # initial conditions Φ(0)=1, Φ̇(0)=0, the exact solution is
        #   Φ(t) = exp(-γt/2) · [cos(Ωt) + (γ/2Ω)·sin(Ωt)]
        # with Ω² = ω₀² - γ²/4 ≈ ω₀². Expanding to O(γ¹):
        #   Φ(t) ≈ cos(ω₀t) + γ·[-t/2·cos(ω₀t) + 1/(2ω₀)·sin(ω₀t)] + O(γ²)
        # Pass 0 alone yields sin(x)·cos(ω₀t); the O(γ¹) correction is
        # the square bracket (times sin(x)). Both terms are required:
        # the first is the amplitude damping, the second is the phase-
        # lag/velocity-IC accommodation.
        cos_term = np.cos(omega0 * t_vals)
        sin_term = np.sin(omega0 * t_vals)
        pass1_phi_factor = gamma * (
            -0.5 * t_vals * cos_term + (1.0 / (2.0 * omega0)) * sin_term
        )
        pass1_expected = np.outer(pass1_phi_factor, np.sin(x))

        # Extract phi slot (index 0) from the flattened state.
        pass1_phi = pass1["y"][:, :n_grid]

        max_analytic = float(np.max(np.abs(pass1_expected)))
        max_diff = float(np.max(np.abs(pass1_phi - pass1_expected)))
        rel_err = max_diff / (max_analytic + 1e-30)

        # O(γ²) truncation error is ~gamma² = 1e-6 relative; allow 1e-3
        # to tolerate the Duhamel snapshot discretisation.
        assert rel_err < 1e-3, (
            f"Pass 1 first_derivative_t correction does not match the "
            f"damped-oscillator expansion. max_diff={max_diff:.3e}, "
            f"max_analytic={max_analytic:.3e}, rel_err={rel_err:.3e}. "
            f"Pre-fix (before #293) the Pass 1 source was missing the "
            f"λ factor so this test would fail by O(1)."
        )

        # Sanity: verify the Pass 1 output is non-trivial. Without the
        # fix applied, Pass 1 from a d_t source on an identity-only base
        # would still be non-zero, so this is a weak check but at least
        # catches a completely broken pipeline.
        assert np.max(np.abs(pass1_phi)) > 1e-6, (
            "Pass 1 output is essentially zero; source wiring is broken."
        )
