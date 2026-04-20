"""R̃² PGT end-to-end regression test (v6 Phase 6C.1 + post-audit R5.1).

Validates that the full v6 pipeline works on the flagship R̃² PGT
theory (``torsion_gertsenshtein.json``):

- Re-derived JSON loads through the post-Ostrogradsky loader
  (Phase 6B) without migration errors.
- ``base_spec(["b5"])`` correctly demotes ``h_4``/``h_7``/``h_9``
  from ``time_order=4`` (kinetic ``2*b5``) to algebraic constraints.
- Pass 0 on a non-gauge photon IC (``a_2 = 0.01 sin(x)``) produces
  finite state evolution AND excites the Gertsenshtein cross-
  coupling: Schur-recovered ``h_4`` / ``h_7`` / ``h_9`` pick up
  B₀-proportional amplitude from the dynamical photon sector.
- All Pass 1 correction terms that cannot be routed at O(ε¹) are
  dropped with a structured diagnostic (R1.1 / #272); every drop
  must be of category ``"row-missing"`` — NEVER of category
  ``"column-missing"`` (that would be a bug).
- In the shipped JSON, every b5-coupled RHS term lives on a
  constraint-field row (h_0/h_4/h_7/h_9/t_*). Therefore at O(ε¹)
  all corrections legitimately collapse to row-missing drops and
  Pass 1 on dynamical fields is **identically zero**. This is a
  theory-level property of R̃² PGT: its b5 coupling does not
  propagate to linear photon / graviton waves until O(ε²), which
  the current driver does not compute (#273). The test enforces
  this invariant so regressions where Pass 1 "magically" picks
  up spurious dynamical contributions are caught.
- Pass 0 trajectory is insensitive to the value of b5 (the
  cleanest check that ``base_spec`` isolates the b5=0 dynamics).

Marked ``@pytest.mark.slow`` because the full-grid simulate takes
several seconds; not in the default fast lane.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.perturbative_driver import PerturbativeSolver
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import load_equation_system

_JSON_PATH = (
    Path(__file__).parent.parent / "examples" / "data" / "torsion_gertsenshtein.json"
)


pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def rtilde_spec() -> object:
    """Load the re-derived torsion_gertsenshtein spec once per module."""
    return load_equation_system(_JSON_PATH)


def test_metadata_has_perturbation(rtilde_spec: Any) -> None:
    """Re-derivation emitted the perturbation metadata block."""
    pert = rtilde_spec.metadata.get("perturbation")
    assert pert == {"small_parameters": ["b5"], "order": 1}


def test_order_in_eps_tags_present(rtilde_spec: Any) -> None:
    """RHS terms containing b5 are tagged order_in_eps=1."""
    n_order_1 = 0
    n_order_0 = 0
    for eq in rtilde_spec.equations:
        for term in eq.rhs_terms:
            if term.order_in_eps == 1:
                n_order_1 += 1
            elif term.order_in_eps == 0:
                n_order_0 += 1
    assert n_order_1 > 0, "At least some RHS terms must be order-1"
    assert n_order_0 > 0, "Base (order-0) terms must also exist"


def test_base_spec_demotes_h_4_h_7_h_9(rtilde_spec: Any) -> None:
    """Gap B: h_4/h_7/h_9 demoted from time_order=4 to algebraic."""
    base = rtilde_spec.base_spec()
    promoted = {"h_4", "h_7", "h_9"}
    pre = {
        eq.field_name: eq.time_derivative_order
        for eq in rtilde_spec.equations
        if eq.field_name in promoted
    }
    post = {
        eq.field_name: eq.time_derivative_order
        for eq in base.equations
        if eq.field_name in promoted
    }
    assert pre == {"h_4": 4, "h_7": 4, "h_9": 4}
    assert post == {"h_4": 0, "h_7": 0, "h_9": 0}


def _photon_ic(solver: PerturbativeSolver, grid: GridInfo) -> np.ndarray:
    """Non-gauge photon IC on a_2 that couples to gravitons via B₀.

    a_2 is a dynamical photon polarisation; the b5=0 base theory
    has a Gertsenshtein-style coupling between a_2 and the demoted
    h_{0,4,7,9} constraint fields via the background B₀ (identity
    terms with coefficient ``±B0/2``, ``±B0²/8``). Exciting a_2
    with a k=1 Fourier mode probes Schur recovery end-to-end.
    """
    layout = StateLayout.from_spec(solver.base_spec, grid.num_points)
    n = grid.num_points
    x = np.linspace(0.0, 2 * np.pi, n, endpoint=False)
    y0 = np.zeros(layout.num_slots * n)
    a2_slot = layout.field_slot_map["a_2"]
    y0[a2_slot * n : (a2_slot + 1) * n] = 0.01 * np.sin(x)
    return y0


def _rtilde_params(b5: float) -> dict[str, float]:
    """Baseline parameter set for R̃² PGT at non-zero b5."""
    return {
        "b5": b5,
        "kappa": 1.0,
        "B0": 0.1,
        "alpha1": 0.1,
        "alpha2": 0.1,
        "alpha3": 0.1,
    }


def test_pass0_exhibits_gertsenshtein_mixing(rtilde_spec: Any) -> None:
    """Pass 0 on a photon IC produces B₀-mediated h-constraint excitation.

    This is the physical Gertsenshtein effect at b5=0: a_2 (photon)
    couples to h_{0,4,7,9} via identity terms proportional to B₀ / B₀²
    in the constraint equations. Schur recovery populates the h
    constraint slots. If this check fails, the Pass 0 Schur path is
    broken regardless of b5.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(rtilde_spec)
    layout = StateLayout.from_spec(solver.base_spec, grid.num_points)
    y0 = _photon_ic(solver, grid)
    res = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.5),
        order=0,  # Pass 0 only — we are testing the base Schur path
        parameters=_rtilde_params(b5=1e-3),
        num_snapshots=5,
    )
    pass0 = res.orders[0]
    a2_slot = layout.field_slot_map["a_2"]
    a2 = pass0["y"][-1, a2_slot * n : (a2_slot + 1) * n]
    assert np.all(np.isfinite(pass0["y"]))
    assert np.max(np.abs(a2)) > 1e-3, (
        "Photon IC should evolve dynamically (a_2 is time_order=2 in base)."
    )

    # Schur-recovered h_4 must pick up B₀-proportional amplitude from
    # the photon via the identity(a_2) term in h_4's constraint eq.
    h4_slot = layout.field_slot_map["h_4"]
    h4 = pass0["y"][-1, h4_slot * n : (h4_slot + 1) * n]
    # Expected magnitude: a_2 amplitude (1e-2) × B₀ coupling (~0.1).
    assert 1e-5 < np.max(np.abs(h4)) < 1e-1, (
        f"Schur-recovered h_4 has implausible magnitude "
        f"{np.max(np.abs(h4)):.3e}. Expected ~B₀·a_2 scale "
        f"(1e-3 .. 1e-2). Gap C / Schur base recovery may be broken."
    )


def test_pass1_has_no_column_missing_drops(rtilde_spec: Any) -> None:
    """R1.1 / #272 + R8 / #290: no column-missing or unclassified drops.

    After R8 augmented Schur recovery (#290), the 159 O(ε¹)
    corrections on demoted constraint rows are routed to
    ``_compute_constraint_source_hat`` and applied via the augmented
    recovery — no longer dropped in the Pass 1 dynamical source
    path. The remaining drop category is ``"row-routed-to-augmented"``
    which is informational, not a bug.

    A ``"column-missing"`` drop means a correction term's TARGET
    field is neither dynamical nor a Schur-recoverable constraint —
    that's always a bug.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(rtilde_spec)
    y0 = _photon_ic(solver, grid)
    res = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.5),
        order=1,
        parameters=_rtilde_params(b5=1e-3),
        num_snapshots=3,
    )
    drops = res.validity.get("correction_drops") or []
    column_missing = [d for d in drops if d["reason"].startswith("column-missing")]
    assert not column_missing, (
        f"Found {len(column_missing)} column-missing drops — these "
        f"indicate a correction term references an unresolved field "
        f"(bug, #272). Records: {column_missing[:3]}"
    )
    # All remaining drops should be the informational R8-routed kind.
    unknown = [
        d
        for d in drops
        if not (
            d["reason"].startswith("row-routed-to-augmented")
            or d["reason"].startswith("column-missing")
        )
    ]
    assert not unknown, (
        f"Found {len(unknown)} drops of unrecognised category; expected "
        f"only row-routed-to-augmented (R8 / #290). First few: {unknown[:3]}"
    )


def test_pass1_h_constraints_nonzero(rtilde_spec: Any) -> None:
    """R8 / #290: augmented Schur recovery populates h_{4,7,9} Pass 1.

    With R8, the 159 O(ε¹) constraint-row corrections (plus the
    LHS-feedback from the removed b5·d^n_t(h_c) kinetics) are
    applied via ``_compute_constraint_source_hat`` + S_cc_inv. Pass 1
    on the demoted constraint fields must now be non-trivial.

    Pre-R8 Pass 1 was identically zero on these fields — the bug
    that motivated the R8 remediation.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(rtilde_spec)
    y0 = _photon_ic(solver, grid)
    res = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.5),
        order=1,
        parameters=_rtilde_params(b5=1e-3),
        num_snapshots=3,
    )
    pass1 = res.orders[1]
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(solver.base_spec, grid.num_points)
    for fname in ("h_4", "h_7", "h_9"):
        slot = layout.field_slot_map[fname]
        vals = pass1["y"][-1, slot * n : (slot + 1) * n]
        assert np.all(np.isfinite(vals))
        max_abs = float(np.max(np.abs(vals)))
        assert max_abs > 1e-9, (
            f"Pass 1 {fname}: max|y| = {max_abs:.3e}, expected "
            f"> 1e-9 after R8 augmented recovery (#290). This was "
            f"identically zero pre-R8 — regression check."
        )


def test_pass1_linear_in_b5_when_excited(rtilde_spec: Any) -> None:
    """R8 / #290: Pass 1 amplitude scales linearly in b5.

    With the augmented Schur recovery landed, Pass 1 on demoted
    constraint fields (h_{4,7,9}) is non-zero and proportional to
    b5 (the small parameter). Doubling b5 doubles Pass 1 amplitude
    to within the O(ε²) truncation error.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(rtilde_spec)
    y0 = _photon_ic(solver, grid)
    res_a = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.3),
        order=1,
        parameters=_rtilde_params(b5=1e-3),
        num_snapshots=3,
    )
    res_b = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.3),
        order=1,
        parameters=_rtilde_params(b5=2e-3),
        num_snapshots=3,
    )
    p1a = res_a.orders[1]["y"]
    p1b = res_b.orders[1]["y"]
    scale_a = float(np.max(np.abs(p1a)))
    scale_b = float(np.max(np.abs(p1b)))
    assert scale_a > 1e-9, (
        f"Pass 1 at b5=1e-3 is essentially zero (max={scale_a:.3e}). "
        f"R8 augmented recovery should produce non-trivial output."
    )
    ratio = scale_b / scale_a
    assert abs(ratio - 2.0) < 0.1, (
        f"Pass 1 amplitude ratio at 2×b5 / 1×b5 = {ratio:.3f}, "
        f"expected 2.0 ± 0.1 for a linear O(ε) correction. "
        f"Indicates non-linear or miscompiled b5 dependence."
    )


def test_pass0_insensitive_to_b5(rtilde_spec: Any) -> None:
    """Changing b5 by a factor of 2 leaves Pass 0 unchanged to rounding.

    Pass 0 only sees order-0 RHS terms (no b5 factors); Pass 1
    absorbs all b5 dependence. If base_spec leaked b5 into Pass 0,
    doubling b5 would shift the Pass 0 trajectory.
    """
    n = 32
    grid = GridInfo(shape=(n,), bounds=((0.0, 2 * np.pi),), periodic=(True,))
    solver = PerturbativeSolver(rtilde_spec)
    y0 = _photon_ic(solver, grid)
    res_small = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.3),
        order=0,
        parameters=_rtilde_params(b5=1e-3),
        num_snapshots=3,
    )
    res_big = solver.solve(
        y0,
        grid,
        t_span=(0.0, 0.3),
        order=0,
        parameters=_rtilde_params(b5=2e-3),
        num_snapshots=3,
    )
    max_diff = float(np.max(np.abs(res_small.total["y"] - res_big.total["y"])))
    assert max_diff < 1e-10, (
        f"Pass 0 trajectory shifts by {max_diff:.3e} when b5 doubles — "
        f"suggests base_spec is leaking b5 into the unperturbed system."
    )
