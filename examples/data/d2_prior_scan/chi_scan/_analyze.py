r"""Reproduce the χ-scan numbers cited in ``docs/tex/pgt_stability_priors.tex``.

Mirrors the structure of the δ₁ scan one directory up
(``examples/data/d2_prior_scan/_analyze.py``) but uses the canonical
probe directly — there is no ``tidal sweep`` step here because the
χ-stability question lives entirely in the pre-flight probe (the
γ_eff = log‖expm(M·t_test) y₀‖/t_test reading), not in any actual
simulation.  The probe is microseconds to milliseconds per call, so
the full four-table scan (≈ 200 calls) runs in seconds and re-runs on
demand to verify the TeX-cited numbers haven't drifted.

Outputs (written next to this file):

* ``chi_single_param.csv``       — the 5-anchor × 12-χ scan at δ₁=0.
  Confirms Finding~1 (\\cref{sec:d21-chi-scan}, lines~489–497):
  γ_eff = 0.101 for every combination; χ alone is benign.
* ``chi_delta1_joint.csv``       — the 5×5 (δ₁,χ) grid at the
  working point.  Confirms the rectangular joint-instability
  boundary (\\cref{tab:d21-joint-scan-working}, lines~505–516).
* ``chi_threshold_by_anchor.csv``— per-anchor χ-threshold scan at
  δ₁=0.005 (\\cref{sec:d21-priors}, lines~531-541).  P5 is the
  most-sensitive anchor and sets χ★_5 = 0.018.
* ``chi_corner_verification.csv``— P5 corner check at δ₁=0.025
  (\\cref{tab:d21-p5-corner-verification}, lines~552–559).

Replaces the wide-format ``chi_scan_table_{a,b,c}.csv`` from
commit 84c4da9: same physics, finer χ resolution (the prior scan
couldn't pin the 0.018 vs 0.020 distinction at sigfig precision; the
finer scan in ``chi_threshold_by_anchor.csv`` resolves it).

Run: ``uv run python examples/data/d2_prior_scan/chi_scan/_analyze.py``
"""

from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

# The probe's matrix-exponential path emits ``RuntimeWarning: overflow``
# from inside scipy when the constraint-Schur matrix is singular.
# That's the *expected* signal — the probe correctly reports
# ``max_excess = inf`` in that case — so suppress the noise here for
# clean stdout when re-running the scan.
warnings.filterwarnings("ignore", category=RuntimeWarning)

from tidal.measurement._stability import check_conversion_stability
from tidal.solver.grid import GridInfo
from tidal.symbolic._spec_cache import load_spec_cached

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SPEC_PATH = REPO / "examples/data/torsion_gertsenshtein_general_nonminimal.json"

# Canonical campaign settings — must match the working invocation
# documented in ``docs/tex/pgt_stability_priors.tex`` and the
# D2.0/D2.1 submission scripts.
GRID = GridInfo(shape=(256,), bounds=((0.0, 100.0),), periodic=(True,))
IC_K = 2 * math.pi / 100.0  # = 0.06283185307179587
KAPPA = 1.0
B0 = 0.01

# The five (β,ξ)-anchor points pinned in
# ``pgt_stability_priors.tex`` (\\cref{tab:d20-scan-points},
# lines~220-225).  IDs match the TeX paragraph labels.
ANCHORS: dict[str, dict[str, float]] = {
    "P1": {"beta1": 0.000, "beta2": -0.600, "beta3": +0.500, "xi": 1.000},
    "P2": {"beta1": +0.206, "beta2": -2.623, "beta3": -0.600, "xi": 0.011},
    "P3": {"beta1": -1.180, "beta2": -1.784, "beta3": -0.258, "xi": 6.030},
    "P4": {"beta1": -0.100, "beta2": -0.450, "beta3": +0.000, "xi": 1.000},
    "P5": {"beta1": +0.000, "beta2": -1.500, "beta3": +0.000, "xi": 1.000},
}

# χ values used by the single-anchor scan.  Covers the prior interior
# at fine resolution (0.005 step in [0, 0.025]) and the wide range up
# to 0.5 to confirm the δ₁=0 invariance generalizes beyond the
# immediate prior.  Includes the cited prior-edge values.
CHI_VALUES: tuple[float, ...] = (
    0.0,
    0.005,
    0.010,
    0.015,
    0.018,
    0.020,
    0.025,
    0.050,
    0.100,
    0.200,
    0.300,
    0.500,
)

# Working-point joint-scan grids (\\cref{tab:d21-joint-scan-working}).
JOINT_DELTA1: tuple[float, ...] = (0.0, 0.001, 0.005, 0.012, 0.020)
JOINT_CHI: tuple[float, ...] = (0.0, 0.005, 0.010, 0.020, 0.050)

# Per-anchor threshold scan (\\cref{sec:d21-priors}, lines~531-541).
# Covers the cited rejection points {0.018, 0.040} and bridging
# values for bracketing at sigfig precision.
THRESHOLD_DELTA1 = 0.005
THRESHOLD_CHI: tuple[float, ...] = (
    0.005,
    0.010,
    0.015,
    0.018,
    0.020,
    0.025,
    0.030,
    0.035,
    0.040,
    0.045,
    0.050,
)

# Prior-corner verification (\\cref{tab:d21-p5-corner-verification}).
# Includes the finer-scan values 0.023, 0.024 from commit 383e20a so
# the pinned threshold at ≈ 0.0235 reproduces.
CORNER_DELTA1 = 0.025
CORNER_CHI: tuple[float, ...] = (
    0.0,
    0.005,
    0.010,
    0.015,
    0.018,
    0.020,
    0.023,
    0.024,
    0.025,
)


def _parameters(
    *,
    beta1: float,
    beta2: float,
    beta3: float,
    xi: float,
    chi: float,
    delta1: float,
) -> dict[str, float]:
    r"""Full parameter dict for the general spec at the campaign defaults.

    The general spec carries ``eta``, ``zeta1``, ``zeta2``, ``zeta3``
    as free symbols; the D2.0/D2.1 invocations all fix these at zero
    (\\cref{sec:d20-priors}, "Operational gotcha").
    """
    return {
        "kappa": KAPPA,
        "B0": B0,
        "beta1": beta1,
        "beta2": beta2,
        "beta3": beta3,
        "xi": xi,
        "chi": chi,
        "delta1": delta1,
        "eta": 0.0,
        "zeta1": 0.0,
        "zeta2": 0.0,
        "zeta3": 0.0,
    }


def _probe(spec: object, params: dict[str, float]) -> tuple[float, bool]:
    """Run the canonical probe and return ``(max_excess, stable)``."""
    result = check_conversion_stability(
        spec,  # type: ignore[arg-type]
        GRID,
        params,
        source="h_5",
        target="a_1",
        ic_wavevector=IC_K,
    )
    return float(result.max_excess), bool(result.stable)


def _format_excess(x: float) -> str:
    if not math.isfinite(x):
        return "inf"
    return f"{x:.6g}"


def run_single_param(spec: object) -> Path:
    """Single-parameter χ scan at δ₁=0 across all 5 anchors."""
    out = HERE / "chi_single_param.csv"
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "anchor_id",
                "beta1",
                "beta2",
                "beta3",
                "xi",
                "chi",
                "max_excess",
                "stable",
            ],
        )
        for anchor_id, anchor in ANCHORS.items():
            for chi in CHI_VALUES:
                params = _parameters(**anchor, chi=chi, delta1=0.0)
                max_excess, stable = _probe(spec, params)
                w.writerow(
                    [
                        anchor_id,
                        f"{anchor['beta1']:.4f}",
                        f"{anchor['beta2']:.4f}",
                        f"{anchor['beta3']:.4f}",
                        f"{anchor['xi']:.4f}",
                        f"{chi:.4f}",
                        _format_excess(max_excess),
                        int(stable),
                    ],
                )
    return out


def run_joint_scan(spec: object) -> Path:
    """Joint (δ₁, χ) scan at the working point P1."""
    out = HERE / "chi_delta1_joint.csv"
    anchor = ANCHORS["P1"]
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chi", "delta1", "max_excess", "stable"])
        for delta1 in JOINT_DELTA1:
            for chi in JOINT_CHI:
                params = _parameters(**anchor, chi=chi, delta1=delta1)
                max_excess, stable = _probe(spec, params)
                w.writerow(
                    [
                        f"{chi:.4f}",
                        f"{delta1:.4f}",
                        _format_excess(max_excess),
                        int(stable),
                    ],
                )
    return out


def run_threshold_by_anchor(spec: object) -> Path:
    r"""Cross-anchor χ-threshold scan at δ₁=0.005.

    P2 is excluded — the small-ξ anchor is already singular at any
    non-zero δ₁ in the χ=0 column of \\cref{tab:d20-general-spec-scan}.
    """
    out = HERE / "chi_threshold_by_anchor.csv"
    targets = ("P1", "P3", "P4", "P5")
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["anchor_id", "delta1", "chi", "max_excess", "stable"])
        for anchor_id in targets:
            anchor = ANCHORS[anchor_id]
            for chi in THRESHOLD_CHI:
                params = _parameters(**anchor, chi=chi, delta1=THRESHOLD_DELTA1)
                max_excess, stable = _probe(spec, params)
                w.writerow(
                    [
                        anchor_id,
                        f"{THRESHOLD_DELTA1:.4f}",
                        f"{chi:.4f}",
                        _format_excess(max_excess),
                        int(stable),
                    ],
                )
    return out


def run_corner_verification(spec: object) -> Path:
    """Prior-corner verification at P5, δ₁=0.025."""
    out = HERE / "chi_corner_verification.csv"
    anchor = ANCHORS["P5"]
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["anchor_id", "delta1", "chi", "max_excess", "stable"])
        for chi in CORNER_CHI:
            params = _parameters(**anchor, chi=chi, delta1=CORNER_DELTA1)
            max_excess, stable = _probe(spec, params)
            w.writerow(
                [
                    "P5",
                    f"{CORNER_DELTA1:.4f}",
                    f"{chi:.4f}",
                    _format_excess(max_excess),
                    int(stable),
                ],
            )
    return out


def _summarize(out_paths: dict[str, Path]) -> None:
    """Print a compact summary tying each output back to its TeX claim."""
    print("\n=== TeX-cited numbers — reproducibility check ===\n")

    # 1. Single-parameter δ₁=0 invariance
    rows = list(csv.DictReader(out_paths["single"].open(encoding="utf-8")))
    excesses = sorted({r["max_excess"] for r in rows})
    invariant = all(e != "inf" and abs(float(e) - 0.1008) < 0.0001 for e in excesses)
    print(
        f"[1] Single-parameter (5 anchors × {len(CHI_VALUES)} χ) at δ₁=0:\n"
        f"    distinct max_excess values across {len(rows)} rows: {excesses}\n"
        f"    TeX claim: γ_eff = 0.1008 for every combination "
        f"({'PASS' if invariant else 'CHECK — values drift'})",
    )

    # 2. Joint working-point boundary
    rows = list(csv.DictReader(out_paths["joint"].open(encoding="utf-8")))
    inf_cells = [(r["delta1"], r["chi"]) for r in rows if r["max_excess"] == "inf"]
    finite_cells = [r for r in rows if r["max_excess"] != "inf"]
    print(
        f"\n[2] Joint (δ₁,χ) at P1 working point:\n"
        f"    finite-γ cells: {len(finite_cells)} / {len(rows)} "
        f"(expected: 21 / 25 — all χ ∈ {{0,0.005,0.01,0.02}} stable + "
        f"χ=0.05/δ₁=0 cell)\n"
        f"    inf cells (δ₁,χ): {inf_cells}",
    )

    # 3. Per-anchor χ-threshold at δ₁=0.005
    rows = list(csv.DictReader(out_paths["threshold"].open(encoding="utf-8")))
    print("\n[3] Per-anchor χ-threshold at δ₁=0.005 (smallest χ where stable=0):")
    for anchor_id in ("P1", "P3", "P4", "P5"):
        anchor_rows = [r for r in rows if r["anchor_id"] == anchor_id]
        first_unstable = next(
            (r for r in anchor_rows if int(r["stable"]) == 0),
            None,
        )
        if first_unstable is None:
            print(
                f"    {anchor_id}: stable through χ={anchor_rows[-1]['chi']} "
                f"(no rejection in scan range)",
            )
        else:
            print(
                f"    {anchor_id}: first rejection at χ={first_unstable['chi']} "
                f"(γ_eff = {first_unstable['max_excess']})",
            )
    print(
        "    TeX claim: P1≈0.040, P3>0.050, P4≈0.040, P5≈0.018 "
        "(P5 binding → χ★_5 = 0.018)",
    )

    # 4. Corner verification at P5/δ₁=0.025
    rows = list(csv.DictReader(out_paths["corner"].open(encoding="utf-8")))
    print("\n[4] P5 corner verification at δ₁=0.025:")
    for r in rows:
        print(
            f"    χ={r['chi']}: γ_eff={r['max_excess']}, stable={r['stable']}",
        )
    finer_rows = [r for r in rows if r["chi"] in {"0.0230", "0.0240"}]
    if len(finer_rows) >= 2:
        print(
            f"    Finer scan pin: χ=0.023 stable={finer_rows[0]['stable']}, "
            f"χ=0.024 stable={finer_rows[1]['stable']} "
            f"(TeX cites threshold ≈ 0.0235)",
        )


def main() -> None:
    spec = load_spec_cached(SPEC_PATH)
    print(f"Loaded spec: {SPEC_PATH.name}")

    out_paths = {
        "single": run_single_param(spec),
        "joint": run_joint_scan(spec),
        "threshold": run_threshold_by_anchor(spec),
        "corner": run_corner_verification(spec),
    }
    for label, path in out_paths.items():
        n_rows = sum(1 for _ in path.open(encoding="utf-8")) - 1  # minus header
        print(f"  wrote {path.name} ({n_rows} rows) [{label}]")

    _summarize(out_paths)


if __name__ == "__main__":
    main()
