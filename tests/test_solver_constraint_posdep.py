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

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_convolution_matrix_with_constraints,
    _build_k_axes,
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


def test_gh379_alpha_zero_baseline() -> None:
    """At α=0 the torsion sector decouples; output must match Gertsenshtein."""
    peaks = _run_minimal_repro(alpha=0.0)
    assert 1e-4 < peaks["h_5"] < 1e-1
    assert 1e-7 < peaks["a_1"] < 1e-2
    # Torsion fields must be EXACTLY zero (no coupling source)
    assert peaks["t_22"] < 1e-12, f"t_22 peak {peaks['t_22']:.3e} should be ≈0 at α=0"


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
        "coupling": {"matrix": [[0, 0], [0, 0]]},
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
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape = (n // 2 + 1,)

    _A_red, recovery, c_names, orig_to_reduced = (
        _build_convolution_matrix_with_constraints(
            spec,
            layout,
            grid,
            ce,
            k_grid,
            rfft_shape,
        )
    )

    assert c_names == ["chi"]
    n_modes = int(np.prod(rfft_shape))

    # Build a known φ(x) = 1.0 (constant); set v_φ=0; apply recovery
    n_dyn_slots = len(orig_to_reduced)
    y_d_hat = np.zeros((n_dyn_slots, n_modes), dtype=np.complex128)
    phi_slot = orig_to_reduced[layout.field_slot_map["phi"]]
    # φ(x) = 1 in real space → FFT is δ_{m,0} · N
    phi_real = np.ones(n)
    y_d_hat[phi_slot] = np.fft.rfft(phi_real)

    c_hat = recovery @ y_d_hat.ravel()
    c_real = np.fft.irfft(c_hat, n=n)

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
