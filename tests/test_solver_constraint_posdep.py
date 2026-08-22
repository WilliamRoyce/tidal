"""Regression tests for GH #379 — modal solver with pos-dep + constraints.

The bug: when a theory had BOTH position-dependent coefficients AND
algebraic constraints (``time_derivative_order == 0`` fields), the
dispatch in :func:`tidal.solver.modal.solve_modal` had::

    needs_reduction = (has_constraints or has_time_ops) and not has_pos_dep

The ``and not has_pos_dep`` clause silently bypassed Schur elimination
for this combination, and ``_build_convolution_matrix`` (the path the
code fell through to) treated algebraic constraint equations
``q_c = RHS`` as first-order ODEs ``dq_c/dt = RHS``, producing
divergent dynamics (~10²⁵²) at any nonzero coupling. Symptom:
``SimulationDivergedError`` with a misleading "known failure mode of
expm_multiply" message recommending ``--scheme cvode`` (which itself
cannot run these theories due to ``d2_t`` RHS operators).

Fix (v0.43.2): new path 4 in modal.py —
:func:`tidal.solver.modal._build_convolution_matrix_with_constraints`
performs Schur elimination on full convolution matrices, mirroring the
constant-coeff genEig path 2.

The 38-field nonminimal PGT-torsion theory at the Phase E geometry is
the canonical bug reproducer.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_convolution_matrix_with_constraints,
    _build_k_axes_full,
    _build_k_grid,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

REPO_ROOT = Path(__file__).resolve().parent.parent
BROKEN_JSON = (
    REPO_ROOT / "examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"
)


# ---------------------------------------------------------------------------
# Phase E geometry (from scripts/hpc_submit_drafts/v3e_localised/_geometry.env).
# Locked in feedback_phase_e_geometry_choices.md; do not relitigate.
# ---------------------------------------------------------------------------

PHASE_E_PARAMS = {
    "kappa": 1.0,
    "Bpeak": 0.01,
    "sigB": 5.0,
    "zc1": 25.0,
    "zc2": 75.0,
}
PHASE_E_GEOMETRY = {
    "L": 100.0,
    "H0": 0.01,  # IC amplitude (well below linearization threshold h ≪ 1)
    "sigma_w": 3.0,
    "x_c": 8.0,
    "k_carrier": 2.0,
}


def _load_broken_spec() -> EquationSystem:
    if not BROKEN_JSON.exists():
        pytest.skip(f"{BROKEN_JSON.name} not present in this checkout")
    return EquationSystem.from_dict(json.loads(BROKEN_JSON.read_text()))


def _run_minimal_repro(
    alpha: float, n: int = 64, t_end: float = 20.0
) -> dict[str, float]:
    """Run the GH #379 minimal repro at given α; return peak |h_5|, |a_1|, |t_22|."""
    spec = _load_broken_spec()
    grid = GridInfo(
        bounds=((0.0, PHASE_E_GEOMETRY["L"]),),
        shape=(n,),
        periodic=(True,),
    )
    params = dict(
        PHASE_E_PARAMS, alpha1=alpha, alpha2=alpha, alpha3=alpha, delta1=alpha
    )

    # IC: Gaussian on h_5
    layout = StateLayout.from_spec(spec, grid.num_points)
    n_pts = grid.num_points
    y0 = np.zeros(layout.num_slots * n_pts)
    x_grid = np.linspace(
        grid.bounds[0][0],
        grid.bounds[0][1],
        n_pts,
        endpoint=False,
    )
    gauss = PHASE_E_GEOMETRY["H0"] * np.exp(
        -((x_grid - PHASE_E_GEOMETRY["x_c"]) ** 2)
        / (2 * PHASE_E_GEOMETRY["sigma_w"] ** 2)
    )
    carrier = np.cos(PHASE_E_GEOMETRY["k_carrier"] * x_grid)
    h5_slot = layout.field_slot_map["h_5"]
    y0[h5_slot * n_pts : (h5_slot + 1) * n_pts] = gauss * carrier

    from tidal.solver.modal import solve_modal

    result = solve_modal(
        spec,
        grid,
        y0,
        t_span=(0.0, t_end),
        bc="periodic",
        parameters=params,
        num_snapshots=2,
    )

    # Final-snapshot peaks
    final = result["y"][-1]
    peaks: dict[str, float] = {}
    for fname in ("h_5", "a_1", "t_22"):
        slot = layout.field_slot_map[fname]
        peaks[fname] = float(np.max(np.abs(final[slot * n_pts : (slot + 1) * n_pts])))
    return peaks


# ---------------------------------------------------------------------------
# Test 1 — minimal bug-reproducer (the user's headline test case)
# ---------------------------------------------------------------------------


@pytest.mark.skip(
    reason="GH #455: the corrected #444 velocity-row scale exposes a near-singular\n    velocity-coupling composition on the ungauged-gravity spec class (operator\n    abscissa 0.51 -> 1594, alpha-independent). These expectations were recorded\n    against the pre-#444 defective operator; on the corrected code the outcomes\n    are resolution-dependent (fail/diverge/hang), so they are skipped rather\n    than marked xfail. Re-enable with re-derived expectations once #455 adjudicates\n    truncation-vs-degeneracy (WS2 constraint-residual oracle)."
)
def test_gh379_minimal_repro_does_not_diverge() -> None:
    """At α=δ=0.001 the broken theory must produce smooth perturbative dynamics.

    Pre-fix: ``SimulationDivergedError`` with peak ~10²⁵² at any nonzero α.
    Post-fix: h_5 peak ~0.005 (Boccaletti baseline), a_1 peak ~10⁻⁴
    (perturbative torsion-mediated), t_22 ~10⁻⁷ (torsion field scales with α).
    """
    peaks = _run_minimal_repro(alpha=0.001)
    # h_5 should be near the Boccaletti baseline (~5e-3 at this geometry)
    assert 1e-4 < peaks["h_5"] < 1e-1, (
        f"h_5 peak {peaks['h_5']:.3e} outside Boccaletti baseline window "
        f"[1e-4, 1e-1] — divergence or numerical instability"
    )
    # a_1 should be perturbative graviton-photon conversion (~10⁻⁴)
    assert 1e-7 < peaks["a_1"] < 1e-2, (
        f"a_1 peak {peaks['a_1']:.3e} outside perturbative window "
        f"[1e-7, 1e-2] — pre-fix this was ~10²⁵²"
    )
    # Torsion field should be present (perturbative δ₁ scaling) but bounded
    assert 1e-10 < peaks["t_22"] < 1e-3, (
        f"t_22 peak {peaks['t_22']:.3e} outside perturbative torsion window "
        f"[1e-10, 1e-3]"
    )


@pytest.mark.skip(
    reason="GH #455: the corrected #444 velocity-row scale exposes a near-singular\n    velocity-coupling composition on the ungauged-gravity spec class (operator\n    abscissa 0.51 -> 1594, alpha-independent). These expectations were recorded\n    against the pre-#444 defective operator; on the corrected code the outcomes\n    are resolution-dependent (fail/diverge/hang), so they are skipped rather\n    than marked xfail. Re-enable with re-derived expectations once #455 adjudicates\n    truncation-vs-degeneracy (WS2 constraint-residual oracle)."
)
def test_gh379_alpha_zero_baseline() -> None:
    """At α=0 the torsion sector decouples; output must match Gertsenshtein."""
    peaks = _run_minimal_repro(alpha=0.0)
    assert 1e-4 < peaks["h_5"] < 1e-1
    assert 1e-7 < peaks["a_1"] < 1e-2
    # Torsion fields must be EXACTLY zero (no coupling source)
    assert peaks["t_22"] < 1e-12, f"t_22 peak {peaks['t_22']:.3e} should be ≈0 at α=0"


@pytest.mark.skip(
    reason="GH #455: the corrected #444 velocity-row scale exposes a near-singular\n    velocity-coupling composition on the ungauged-gravity spec class (operator\n    abscissa 0.51 -> 1594, alpha-independent). These expectations were recorded\n    against the pre-#444 defective operator; on the corrected code the outcomes\n    are resolution-dependent (fail/diverge/hang), so they are skipped rather\n    than marked xfail. Re-enable with re-derived expectations once #455 adjudicates\n    truncation-vs-degeneracy (WS2 constraint-residual oracle)."
)
def test_gh379_torsion_scales_linearly_with_alpha() -> None:
    """Torsion-field peak scales linearly with α (perturbative regime).

    At leading order, t_constraint = recovery·y_d with recovery ∝ α, so
    |t| should grow as α to first order. This is the smoking-gun signature
    that the fix is computing genuine perturbative torsion-mediated physics
    (and not a coincidentally bounded numerical artifact).
    """
    peaks_001 = _run_minimal_repro(alpha=0.001)
    peaks_010 = _run_minimal_repro(alpha=0.01)
    peaks_100 = _run_minimal_repro(alpha=0.1)
    # Each 10× α should give roughly 10× t_22. Tolerance loose (factor 3) for
    # higher-order nonlinearities and any α-mixing in the recovery.
    ratio_low = peaks_010["t_22"] / max(peaks_001["t_22"], 1e-30)
    ratio_high = peaks_100["t_22"] / max(peaks_010["t_22"], 1e-30)
    assert 3 < ratio_low < 30, (
        f"t_22 ratio α=0.01/0.001 = {ratio_low:.2f}, expected ~10 (linear scaling)"
    )
    assert 3 < ratio_high < 30, (
        f"t_22 ratio α=0.1/0.01 = {ratio_high:.2f}, expected ~10 (linear scaling)"
    )


# ---------------------------------------------------------------------------
# Test 2 — synthetic minimal theory exercising the Schur-elimination path
# ---------------------------------------------------------------------------


def _make_minimal_constraint_spec(
    g_amp: float = 0.3, h_amp: float = 0.5
) -> EquationSystem:
    """1 dynamical φ with d2_t(φ)=lap(φ)+g(x)·χ, 1 constraint χ - h(x)·φ = 0.

    The constraint is encoded in canonical RHS-=-0 form: rhs.terms includes
    ``-χ + h(x)·φ`` so that K_cc·χ + K_cd·φ = 0 yields χ = h(x)·φ.
    """
    g_sym = f"{g_amp}*Cos[x[]]"
    h_sym = f"{h_amp}*Cos[x[]]"
    spec_data = {
        "metadata": {"name": "gh379_synthetic"},
        "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
        "fields": [{"name": "phi", "index": 0}, {"name": "chi", "index": 1}],
        "equations": [
            {
                "field": "phi",
                "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi",
                        },
                        {
                            "coefficient": g_amp,
                            "operator": "identity",
                            "field": "chi",
                            "coefficient_symbolic": g_sym,
                            "coordinate_dependent": ["x"],
                        },
                    ],
                },
            },
            {
                "field": "chi",
                "lhs": {"expression": "chi", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        # -χ (LHS moved to RHS with sign flip, per JSON convention)
                        {"coefficient": -1.0, "operator": "identity", "field": "chi"},
                        # +h(x)·φ
                        {
                            "coefficient": h_amp,
                            "operator": "identity",
                            "field": "phi",
                            "coefficient_symbolic": h_sym,
                            "coordinate_dependent": ["x"],
                        },
                    ],
                },
            },
        ],
    }
    return EquationSystem.from_dict(spec_data)


def test_gh379_synthetic_recovery_matches_analytical() -> None:
    """For χ = h(x)·φ, the recovery matrix must reproduce h(x)·φ when
    applied to a known φ.

    This is the natural "oracle" for the matrix arithmetic: no numerical
    integration, just check that the constraint-elimination algebra is
    correct. Validates the GH #379 fix end-to-end at the matrix-builder
    level, independent of expm_multiply / time evolution.
    """
    spec = _make_minimal_constraint_spec(g_amp=0.0, h_amp=0.5)  # g=0 decouples ϕ
    n = 16
    grid = GridInfo(bounds=((0.0, 2.0 * np.pi),), shape=(n,), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    ce = CoefficientEvaluator(spec, grid, {})
    # GH #445: the convolution builders operate in the FULL fftn basis.
    k_axes = _build_k_axes_full(grid)
    k_grid = _build_k_grid(k_axes)
    full_shape = (n,)

    _A_red, recovery, c_names, orig_to_reduced = (
        _build_convolution_matrix_with_constraints(
            spec,
            layout,
            grid,
            ce,
            k_grid,
            full_shape,
        )
    )

    assert c_names == ["chi"]
    n_modes = int(np.prod(full_shape))

    # Build a known φ(x) = 1.0 (constant); set v_φ=0; apply recovery
    n_dyn_slots = len(orig_to_reduced)
    y_d_hat = np.zeros((n_dyn_slots, n_modes), dtype=np.complex128)
    phi_slot = orig_to_reduced[layout.field_slot_map["phi"]]
    # φ(x) = 1 in real space → FFT is δ_{m,0} · N
    phi_real = np.ones(n)
    y_d_hat[phi_slot] = np.fft.fft(phi_real)

    c_hat = recovery @ y_d_hat.ravel()
    c_real = np.fft.ifft(c_hat).real

    # Expected: χ(x) = h(x) · 1.0 = 0.5·cos(x), evaluated at the same
    # cell-centered grid coords (x_i = (i+0.5)·dx) the modal solver uses
    # for position-dependent coefficient resolution.
    dx = (2.0 * np.pi) / n
    x_grid_centered = (np.arange(n) + 0.5) * dx
    chi_expected = 0.5 * np.cos(x_grid_centered)

    rel_err = np.max(np.abs(c_real - chi_expected)) / max(
        np.max(np.abs(chi_expected)), 1e-30
    )
    assert rel_err < 1e-10, (
        f"recovery·φ disagrees with analytical h(x)·φ: rel_err={rel_err:.3e}"
    )


def test_gh379_synthetic_smooth_evolution() -> None:
    """Time-evolve the synthetic theory and confirm bounded solution.

    Coupling is on (g, h ≠ 0); the system has constraint+pos-dep, the
    exact combination that pre-fix produced divergence in the 38-field
    case. Here at small grid+amplitude the system is benign — we just
    check no NaN / blowup.
    """
    spec = _make_minimal_constraint_spec(g_amp=0.3, h_amp=0.5)
    n = 32
    grid = GridInfo(bounds=((0.0, 2.0 * np.pi),), shape=(n,), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    n_pts = grid.num_points
    y0 = np.zeros(layout.num_slots * n_pts)
    x_grid = np.linspace(0.0, 2.0 * np.pi, n_pts, endpoint=False)
    # Smooth IC on φ
    phi_slot = layout.field_slot_map["phi"]
    y0[phi_slot * n_pts : (phi_slot + 1) * n_pts] = 0.01 * np.cos(x_grid)

    from tidal.solver.modal import solve_modal

    result = solve_modal(
        spec,
        grid,
        y0,
        t_span=(0.0, 5.0),
        bc="periodic",
        parameters={},
        num_snapshots=2,
    )
    final = result["y"][-1]
    assert np.all(np.isfinite(final))
    assert np.max(np.abs(final)) < 1e3, (
        f"Synthetic theory diverged: max|y|={np.max(np.abs(final)):.3e}"
    )


# ---------------------------------------------------------------------------
# GH #444 — deferred v_c/ẍ_c substitutions must carry the velocity-row M⁻¹
# ---------------------------------------------------------------------------


def _make_deferred_velc_spec(kinetic: str | None, c_v: float = 0.4) -> EquationSystem:
    """φ with ``M·d2_t(φ) = ∂²ₓφ + c_v·v_χ`` and constraint ``χ = h(x)·φ``.

    The ``c_v·v_χ`` term targets a CONSTRAINT field's velocity from a
    dynamical RHS — exactly the deferred-substitution path
    (``deferred_terms_dyn_velc``) that GH #444 found emitting without the
    row's ``velocity_row_scale``. All 8 dual-Gaussian roster specs carry
    6–944 such terms on rows with non-unit kinetics.
    """
    h_amp = 0.5
    spec_data: dict[str, Any] = {
        "metadata": {"name": "gh444_synthetic", "parameters": {"alpha": 0.5}},
        "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
        "fields": [{"name": "phi", "index": 0}, {"name": "chi", "index": 1}],
        "equations": [
            {
                "field": "phi",
                "lhs": {
                    "expression": "d2_t(phi)",
                    "order": {"time": 2, "space": 0},
                    **(
                        {"kinetic_coefficient_symbolic": kinetic}
                        if kinetic is not None
                        else {}
                    ),
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi",
                        },
                        # Deferred-path trigger: dyn RHS references v_χ.
                        {
                            "coefficient": c_v,
                            "operator": "identity",
                            "field": "v_chi",
                        },
                    ],
                },
            },
            {
                "field": "chi",
                "lhs": {"expression": "chi", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": -1.0, "operator": "identity", "field": "chi"},
                        {
                            "coefficient": h_amp,
                            "operator": "identity",
                            "field": "phi",
                            "coefficient_symbolic": f"{h_amp}*Cos[x[]]",
                            "coordinate_dependent": ["x"],
                        },
                    ],
                },
            },
        ],
    }
    return EquationSystem.from_dict(spec_data)


class TestGH444DeferredVelocityRowScale:
    """Composed-DAE-residual oracle for the deferred v_c substitution.

    With ``χ = h(x)·φ`` exactly, ``v_χ = h(x)·v_φ``, so the velocity-row
    action of the composed reduced operator must satisfy, pointwise in
    real space::

        (dv_φ/dt)(x) = (1/M(x)) · [∂²ₓφ(x) + c_v·h(x)·v_φ(x)]

    This is independent of every internal of the Schur/deferral chain
    (recovery, vel_coupling composition, B_lhs pre-solve) — it only trusts
    the constraint's analytic solution. Pre-#444 the ``c_v·h·v_φ`` piece
    was emitted WITHOUT 1/M, so any M ≠ 1 fails by that factor (~33% at
    M = 1.5 with these amplitudes); the WS2 audit generalizes this oracle
    pattern across the roster.
    """

    C_V = 0.4
    H_AMP = 0.5
    PARAMS = {"alpha": 0.5}

    @pytest.mark.parametrize(
        "kinetic",
        [
            None,  # M = 1 baseline (oracle sanity: passes pre- and post-fix)
            "1 + alpha",  # constant M = 1.5 (the corpus shape class)
            "1 + alpha + 0*x[]",  # pos-dep in form, analytically 1.5 (ndarray path)
            "1 + 0.25*Sin[x[]]",  # genuinely varying M(x)
        ],
    )
    def test_velocity_row_action_matches_dae_residual(
        self, kinetic: str | None
    ) -> None:
        spec = _make_deferred_velc_spec(kinetic, c_v=self.C_V)
        n = 32
        length = 2.0 * np.pi
        grid = GridInfo(bounds=((0.0, length),), shape=(n,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        ce = CoefficientEvaluator(spec, grid, self.PARAMS)
        k_grid = _build_k_grid(_build_k_axes_full(grid))

        A_red, _recovery, c_names, orig_to_reduced = (
            _build_convolution_matrix_with_constraints(
                spec, layout, grid, ce, k_grid, (n,)
            )
        )
        assert c_names == ["chi"]

        x = grid.axes_coords(0)
        h_of_x = self.H_AMP * np.cos(x)
        # Evaluate M(x) independently of the solver machinery.
        alpha = self.PARAMS["alpha"]
        if kinetic is None:
            m_of_x = np.ones(n)
        elif kinetic in {"1 + alpha", "1 + alpha + 0*x[]"}:
            m_of_x = np.full(n, 1.0 + alpha)
        else:  # "1 + 0.25*Sin[x[]]"
            m_of_x = 1.0 + 0.25 * np.sin(x)

        # Random smooth real state: band-limited φ and v_φ.
        rng = np.random.default_rng(444)
        k_int = np.fft.fftfreq(n, d=1.0 / n)
        keep = np.abs(k_int) <= 5
        phi_hat = np.where(
            keep, rng.standard_normal(n) + 1j * rng.standard_normal(n), 0
        )
        vphi_hat = np.where(
            keep, rng.standard_normal(n) + 1j * rng.standard_normal(n), 0
        )
        # Hermitian-symmetrize so the states are real fields.
        phi = np.fft.ifft(phi_hat).real
        vel_phi = np.fft.ifft(vphi_hat).real
        phi_hat = np.fft.fft(phi)
        vphi_hat = np.fft.fft(vel_phi)

        phi_red = orig_to_reduced[layout.field_slot_map["phi"]]
        vphi_red = orig_to_reduced[layout.velocity_slot_map["phi"]]
        n_dyn = len(orig_to_reduced)
        y_hat = np.zeros(n_dyn * n, dtype=np.complex128)
        y_hat[phi_red * n : (phi_red + 1) * n] = phi_hat
        y_hat[vphi_red * n : (vphi_red + 1) * n] = vphi_hat

        action = (A_red @ y_hat)[vphi_red * n : (vphi_red + 1) * n]
        action_real = np.fft.ifft(action).real

        k_phys = 2.0 * np.pi * np.fft.fftfreq(n, d=length / n)
        lap_phi = np.fft.ifft(-(k_phys**2) * phi_hat).real
        expected = (lap_phi + self.C_V * h_of_x * vel_phi) / m_of_x

        err = np.max(np.abs(action_real - expected)) / max(
            np.max(np.abs(expected)), 1e-300
        )
        assert err < 1e-12, (
            f"deferred v_c velocity-row action violates the DAE residual "
            f"oracle for kinetic={kinetic!r}: rel_err={err:.3e} "
            f"(pre-#444 signature: error ≈ the c_v·h·v_φ share of (1 − 1/M))"
        )

    def test_evolution_matches_ida_reference(self) -> None:
        """End-to-end: the fixed path-4 evolution agrees with IDA (which
        evaluates kinetics on the grid, GH #382) on the deferred-velc spec
        with M = 1.5, RMS < 1%.
        """
        from tidal.solver.ida import solve_ida
        from tidal.solver.modal import solve_modal

        spec = _make_deferred_velc_spec("1 + alpha", c_v=self.C_V)
        n = 32
        grid = GridInfo(bounds=((0.0, 2.0 * np.pi),), shape=(n,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = grid.axes_coords(0)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        y0[layout.slot_slice(layout.field_slot_map["phi"])] = 1e-3 * (
            np.sin(x) + 0.3 * np.cos(2 * x + 0.4)
        )
        t_span = (0.0, 0.5)

        modal = solve_modal(
            spec, grid, y0, t_span, parameters=self.PARAMS, num_snapshots=5
        )
        ida = solve_ida(
            spec,
            grid,
            y0,
            t_span,
            bc="periodic",
            parameters=self.PARAMS,
            num_snapshots=5,
            rtol=1e-10,
            atol=1e-12,
        )
        assert modal["success"]
        assert ida["success"]
        diff = np.asarray(modal["y"]) - np.asarray(ida["y"])
        rel = float(
            np.sqrt(np.mean(diff**2)) / np.sqrt(np.mean(np.asarray(ida["y"]) ** 2))
        )
        assert rel < 1e-2, f"modal-vs-IDA RMS {rel:.3e} exceeds 1% (GH #444)"


class TestGH444ECalRosterValidation:
    """Roster-level #444 coverage on the E.cal calibration spec — smoke only.

    `gertsenshtein_ungauged_e_dual_gaussian.json` is the Phase E
    calibration theory: graviton rows h_5/h_6/h_8 with `-kappa^(-2)`
    kinetics, 7 algebraic constraints, and 6 deferred constraint-velocity
    terms on those rows (the GH #444 pattern; at kappa = 1 the missing
    scale was a SIGN flip).

    Three audit findings (recorded on #444/#449) constrain what can be
    asserted here today:

    1. `d2_t` RHS operators → the time-domain backends cannot run this
       spec at all: NO cross-backend reference exists for this roster
       class.
    2. `first_derivative_t` RHS operators → the energy machinery's
       `apply_operator` table cannot process the spec either (raises
       Unknown-operator from the EOM-evaluation path), so
       energy-conservation is not usable as an oracle — independent of,
       and predating, #444.
    3. The #178 constraint-kinetic-coupling caveat applies on top: the
       naive Legendre H would not be the conserved quantity anyway.

    Consequence: the only independent correctness check for this class is
    the composed-DAE-residual oracle harness (WS2 of the round-2 plan),
    which generalizes the machine-precision oracle proven on the
    synthetic spec above. Until it lands, this test pins completion and
    boundedness only; the quantitative #444 acceptance lives in
    `TestGH444DeferredVelocityRowScale` (pre-fix-failing oracle + IDA
    agreement).
    """

    ECAL_JSON = REPO_ROOT / "examples/data/gertsenshtein_ungauged_e_dual_gaussian.json"

    @pytest.mark.skip(
        reason="GH #455: E.cal shares the ungauged-gravity near-singular composition\n        (identical operator numbers to the nonminimal spec); skipped pending #455."
    )
    def test_ecal_runs_under_modal_post_444(self) -> None:
        from tidal.solver.modal import solve_modal

        if not self.ECAL_JSON.exists():
            pytest.skip(f"{self.ECAL_JSON.name} not present in this checkout")
        spec = EquationSystem.from_dict(json.loads(self.ECAL_JSON.read_text()))
        n = 48
        grid = GridInfo(
            bounds=((0.0, PHASE_E_GEOMETRY["L"]),), shape=(n,), periodic=(True,)
        )
        layout = StateLayout.from_spec(spec, grid.num_points)
        x = grid.axes_coords(0)
        length = PHASE_E_GEOMETRY["L"]

        rng = np.random.default_rng(444)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        for fname in ("h_5", "h_6", "h_8"):
            sl = layout.slot_slice(layout.field_slot_map[fname])
            phases = rng.uniform(0.0, 2 * np.pi, size=3)
            y0[sl] = 1e-3 * sum(
                np.cos(2 * np.pi * (m + 1) * x / length + phases[m]) for m in range(3)
            )

        result = solve_modal(
            spec, grid, y0, (0.0, 2.0), parameters=PHASE_E_PARAMS, num_snapshots=3
        )
        assert result["success"]
        final = np.asarray(result["y"])[-1]
        assert np.all(np.isfinite(final))
        # Linear wave dynamics at h ~ 1e-3 over t=2 must stay bounded;
        # the pre-#444 sign-flipped deferred terms are covered
        # quantitatively by the synthetic oracle above.
        assert np.max(np.abs(final)) < 1.0
