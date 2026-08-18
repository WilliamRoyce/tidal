r"""Shapiro ζ₁/ζ₂/ζ₃ prior sensitivity scan for the D2.2 campaign.

Follows the per-model recipe in ``docs/tex/pgt_stability_priors.tex``
§"Per-model sensitivity scan recipe", using the D2.1 χ-scan
(``examples/data/d2_prior_scan/chi_scan/_analyze.py``) as a template.

The five anchor points re-use the D2.0/D2.1 MAP estimates plus spread
points covering the stable β-region:
  P1  D2.0 Bahamonde MAP  β₁=−0.10, β₂=−1.55, β₃=+0.96, ξ=1.85
  P2  D2.1 Barker MAP     β₁=+0.73, β₂=−1.91, β₃=−0.67, ξ=0.0105
  P3  D1 working + ξ=1   β₁=0.00, β₂=−0.60, β₃=+0.50, ξ=1.000
  P4  Boundary spread     β₁=−0.10, β₂=−0.45, β₃=0.00, ξ=1.000
  P5  Interior spread     β₁=+0.00, β₂=−1.50, β₃=0.00, ξ=1.000

Scans performed:
  zeta{i}_single_param.csv    — ζᵢ ∈ logspace(−4,+1,41) both signs, δ₁=0
                                 across all five anchors; threshold ζᵢ* where
                                 γ_eff > 0.3 first appears.
  zeta{i}_delta1_joint.csv    — 5×5 (δ₁, ζᵢ) grid at P1 to identify
                                 rectangular vs multiplicative structure.
  zeta{i}_zeta{j}_pair.csv    — 5×5 (ζᵢ, ζⱼ) at P1/δ₁=0 for (i,j) ∈
                                 {(1,2),(1,3),(2,3)}.
  corner_verification.csv     — 9 prior-corner points spanning the proposed
                                 prior edges; acceptance labels.

Run:  ``uv run python examples/data/d2_prior_scan/zeta_scan/_analyze.py``
"""

from __future__ import annotations

import csv
import math
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

from tidal.measurement._stability import check_conversion_stability
from tidal.solver.grid import GridInfo
from tidal.symbolic._spec_cache import load_spec_cached

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SPEC_PATH = REPO / "examples/data/torsion_gertsenshtein_general_nonminimal.json"

GRID = GridInfo(shape=(256,), bounds=((0.0, 100.0),), periodic=(True,))
IC_K = 2 * math.pi / 100.0  # 0.06283185307179587
KAPPA = 1.0
B0 = 0.01

# Five anchor points: D2.0 MAP, D2.1 MAP, D1 working + ξ=1, plus two spread
# across the stable β-triangle (reusing chi_scan P4/P5 which are already
# verified stable).
ANCHORS: dict[str, dict[str, float]] = {
    "P1": {"beta1": -0.100, "beta2": -1.550, "beta3": +0.960, "xi": 1.850},
    "P2": {"beta1": +0.730, "beta2": -1.910, "beta3": -0.670, "xi": 0.0105},
    "P3": {"beta1": 0.000, "beta2": -0.600, "beta3": +0.500, "xi": 1.000},
    "P4": {"beta1": -0.100, "beta2": -0.450, "beta3": 0.000, "xi": 1.000},
    "P5": {"beta1": 0.000, "beta2": -1.500, "beta3": 0.000, "xi": 1.000},
}

# Single-parameter scan values: logspace magnitudes, both signs + zero.
_MAGS = np.logspace(-4, 1, 41)
SINGLE_PARAM_VALUES: tuple[float, ...] = tuple(
    sorted({0.0, *list(_MAGS), *list(-_MAGS)})
)

# Joint (δ₁, ζᵢ) grid values — set {0, ±0.005, ±0.02}.
JOINT_DELTA1: tuple[float, ...] = (-0.020, -0.005, 0.0, 0.005, 0.020)

# Stability threshold for max_excess — probe uses conservative mode.
GAMMA_THRESHOLD = 0.3


def _parameters(
    *,
    beta1: float,
    beta2: float,
    beta3: float,
    xi: float,
    delta1: float = 0.0,
    zeta1: float = 0.0,
    zeta2: float = 0.0,
    zeta3: float = 0.0,
) -> dict[str, float]:
    return {
        "kappa": KAPPA,
        "B0": B0,
        "beta1": beta1,
        "beta2": beta2,
        "beta3": beta3,
        "xi": xi,
        "chi": 0.0,
        "eta": 0.0,
        "delta1": delta1,
        "zeta1": zeta1,
        "zeta2": zeta2,
        "zeta3": zeta3,
    }


def _probe(spec: object, params: dict[str, float]) -> tuple[float, bool]:
    result = check_conversion_stability(
        spec,  # type: ignore[arg-type]
        GRID,
        params,
        source="h_5",
        target="a_1",
        ic_wavevector=IC_K,
    )
    return float(result.max_excess), bool(result.stable)


def _fmt(x: float) -> str:
    if not math.isfinite(x):
        return "inf"
    return f"{x:.6g}"


# ---------------------------------------------------------------------------
# Scan 1: single-parameter  ζᵢ at δ₁=0, all anchors
# ---------------------------------------------------------------------------


def run_single_param(spec: object) -> dict[str, Path]:
    """ζᵢ-only scan at δ₁=0 across all anchors. Returns path per ζᵢ."""
    paths: dict[str, Path] = {}
    for zi_name in ("zeta1", "zeta2", "zeta3"):
        out = HERE / f"{zi_name}_single_param.csv"
        with out.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "anchor_id",
                    "beta1",
                    "beta2",
                    "beta3",
                    "xi",
                    zi_name,
                    "max_excess",
                    "stable",
                ]
            )
            for anchor_id, anchor in ANCHORS.items():
                for zval in SINGLE_PARAM_VALUES:
                    kw = {zi_name: zval}
                    params = _parameters(**anchor, **kw)
                    exc, stable = _probe(spec, params)
                    w.writerow(
                        [
                            anchor_id,
                            f"{anchor['beta1']:.4f}",
                            f"{anchor['beta2']:.4f}",
                            f"{anchor['beta3']:.4f}",
                            f"{anchor['xi']:.5f}",
                            f"{zval:.8g}",
                            _fmt(exc),
                            int(stable),
                        ]
                    )
        paths[zi_name] = out
    return paths


def _find_threshold(csv_path: Path) -> float | None:
    """Return the smallest |ζ| where stable=0, or None if all stable."""
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    unstable = [
        abs(float(r[csv_path.stem.replace("_single_param", "")]))
        for r in rows
        if int(r["stable"]) == 0
    ]
    return min(unstable) if unstable else None


# ---------------------------------------------------------------------------
# Scan 2: joint (δ₁, ζᵢ) at P1
# ---------------------------------------------------------------------------


def _joint_zeta_values(threshold: float | None) -> tuple[float, ...]:
    """5 ζ values centered on 0: {0, ±half, ±full}."""
    if threshold is None:
        # benign alone up to max scan range — use a conservative round number
        half = 0.5
        full = 1.0
    else:
        half = threshold / 2.0
        full = threshold
    return (-full, -half, 0.0, half, full)


def run_joint_scans(
    spec: object,
    thresholds: dict[str, float | None],
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    anchor = ANCHORS["P1"]
    for zi_name, thr in thresholds.items():
        zvals = _joint_zeta_values(thr)
        out = HERE / f"{zi_name}_delta1_joint.csv"
        with out.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["delta1", zi_name, "max_excess", "stable"])
            for d1 in JOINT_DELTA1:
                for zval in zvals:
                    kw = {zi_name: zval}
                    params = _parameters(**anchor, delta1=d1, **kw)
                    exc, stable = _probe(spec, params)
                    w.writerow([f"{d1:.4f}", f"{zval:.8g}", _fmt(exc), int(stable)])
        paths[zi_name] = out
    return paths


# ---------------------------------------------------------------------------
# Scan 3: pairwise (ζᵢ, ζⱼ) at P1, δ₁=0
# ---------------------------------------------------------------------------


def run_pairwise_scans(
    spec: object,
    thresholds: dict[str, float | None],
) -> dict[tuple[str, str], Path]:
    paths: dict[tuple[str, str], Path] = {}
    anchor = ANCHORS["P1"]
    pairs = [("zeta1", "zeta2"), ("zeta1", "zeta3"), ("zeta2", "zeta3")]
    for zi_name, zj_name in pairs:
        zi_vals = _joint_zeta_values(thresholds[zi_name])
        zj_vals = _joint_zeta_values(thresholds[zj_name])
        short = f"{zi_name}_{zj_name}_pair"
        out = HERE / f"{short}.csv"
        with out.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow([zi_name, zj_name, "max_excess", "stable"])
            for ziv in zi_vals:
                for zjv in zj_vals:
                    params = _parameters(
                        **anchor, delta1=0.0, **{zi_name: ziv, zj_name: zjv}
                    )
                    exc, stable = _probe(spec, params)
                    w.writerow([f"{ziv:.8g}", f"{zjv:.8g}", _fmt(exc), int(stable)])
        paths[zi_name, zj_name] = out
    return paths


# ---------------------------------------------------------------------------
# Scan 4: corner verification
# ---------------------------------------------------------------------------


def run_corner_verification(
    spec: object,
    prior_bounds: dict[str, float],
) -> Path:
    """9 points spanning proposed prior corners with acceptance labels."""
    d1_max = 0.025  # D2.0 prior edge
    z1 = prior_bounds["zeta1"]
    z2 = prior_bounds["zeta2"]
    z3 = prior_bounds["zeta3"]

    corners = [
        # (anchor_id, delta1, zeta1, zeta2, zeta3, label)
        ("P1", +d1_max, +z1, 0.0, 0.0, "P1 d1+ z1+"),
        ("P1", +d1_max, -z1, 0.0, 0.0, "P1 d1+ z1-"),
        ("P1", +d1_max, 0.0, +z2, 0.0, "P1 d1+ z2+"),
        ("P1", +d1_max, 0.0, 0.0, +z3, "P1 d1+ z3+"),
        ("P5", +d1_max, +z1, 0.0, 0.0, "P5 d1+ z1+ (interior)"),
        ("P5", +d1_max, 0.0, +z2, 0.0, "P5 d1+ z2+ (interior)"),
        ("P5", +d1_max, 0.0, 0.0, +z3, "P5 d1+ z3+ (interior)"),
        ("P2", +d1_max, +z1, 0.0, 0.0, "P2 d1+ z1+ (small-xi)"),
        ("P1", +d1_max, +z1, +z2, +z3, "P1 d1+ all-zeta+"),
    ]

    out = HERE / "corner_verification.csv"
    with out.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "label",
                "anchor_id",
                "delta1",
                "zeta1",
                "zeta2",
                "zeta3",
                "max_excess",
                "stable",
            ]
        )
        for anchor_id, d1, z1v, z2v, z3v, label in corners:
            anchor = ANCHORS[anchor_id]
            params = _parameters(
                **anchor,
                delta1=d1,
                zeta1=z1v,
                zeta2=z2v,
                zeta3=z3v,
            )
            exc, stable = _probe(spec, params)
            w.writerow(
                [
                    label,
                    anchor_id,
                    f"{d1:.4f}",
                    f"{z1v:.6g}",
                    f"{z2v:.6g}",
                    f"{z3v:.6g}",
                    _fmt(exc),
                    int(stable),
                ]
            )
    return out


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------


def _classify_joint(csv_path: Path, zi_name: str) -> str:
    """Classify joint structure as 'rectangular', 'multiplicative', or 'benign'."""
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    inf_cells = [
        (float(r["delta1"]), float(r[zi_name]))
        for r in rows
        if r["max_excess"] == "inf"
    ]
    if not inf_cells:
        return "benign (no rejection in joint grid)"
    # Rectangular: rejection appears when BOTH delta1 != 0 AND |zeta| > thr
    # Multiplicative: rejection appears when product delta1*zeta > thr
    # Check: is every inf cell at delta1 != 0?
    d1_zero_infs = [c for c in inf_cells if abs(c[0]) < 1e-9]
    if d1_zero_infs:
        return f"NOT rectangular — instability at δ₁=0 at {d1_zero_infs}"
    # Check: for non-zero d1 rows, at what zeta does it first go inf?
    thresholds_by_d1: dict[float, float] = {}
    for d1 in {c[0] for c in inf_cells}:
        min_zabs = min(abs(c[1]) for c in inf_cells if c[0] == d1)
        thresholds_by_d1[d1] = min_zabs
    if len({round(v, 4) for v in thresholds_by_d1.values()}) <= 2:
        return f"rectangular (rejection at δ₁≠0 + |ζ|≥threshold, threshold ~{min(thresholds_by_d1.values()):.4f})"
    return f"complex (threshold varies with δ₁): {thresholds_by_d1}"


def _classify_pairwise(csv_path: Path, zi_name: str, zj_name: str) -> str:
    rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
    inf_cells = [
        (float(r[zi_name]), float(r[zj_name])) for r in rows if r["max_excess"] == "inf"
    ]
    if not inf_cells:
        return "benign (no rejection in pairwise grid)"
    both_zero = [c for c in inf_cells if abs(c[0]) < 1e-9 and abs(c[1]) < 1e-9]
    one_zero = [c for c in inf_cells if abs(c[0]) < 1e-9 or abs(c[1]) < 1e-9]
    return (
        f"has rejection ({len(inf_cells)} cells); "
        f"at-(0,0): {len(both_zero)}, on-axis: {len(one_zero)}"
    )


def _summarize(
    joint_paths: dict[str, Path],
    pair_paths: dict[tuple[str, str], Path],
    corner_path: Path,
    thresholds: dict[str, float | None],
    prior_bounds: dict[str, float],
) -> None:
    print("\n=== ζ-scan reproducibility summary ===\n")
    print("Single-parameter thresholds (γ_eff > 0.3 criterion):")
    for zi, thr in thresholds.items():
        pb = prior_bounds[zi]
        if thr is None:
            print(
                f"  {zi}: benign across full logspace range (no threshold found) → prior = ±{pb:.4g}"
            )
        else:
            print(
                f"  {zi}*_empirical = {thr:.4g} → prior = ±{pb:.4g} (×0.5 safety margin)"
            )
    print()
    print("Joint (δ₁, ζᵢ) structure at P1:")
    for zi, path in joint_paths.items():
        verdict = _classify_joint(path, zi)
        print(f"  {zi}: {verdict}")
    print()
    print("Pairwise (ζᵢ, ζⱼ) structure at P1, δ₁=0:")
    for (zi, zj), path in pair_paths.items():
        verdict = _classify_pairwise(path, zi, zj)
        print(f"  ({zi}, {zj}): {verdict}")
    print()
    print("Corner verification:")
    rows = list(csv.DictReader(corner_path.open(encoding="utf-8")))
    for r in rows:
        status = "PASS" if int(r["stable"]) == 1 else "FAIL"
        print(f"  [{status}] {r['label']}: γ_eff={r['max_excess']}")


def _write_summary_md(
    joint_paths: dict[str, Path],
    pair_paths: dict[tuple[str, str], Path],
    corner_path: Path,
    thresholds: dict[str, float | None],
    prior_bounds: dict[str, float],
) -> Path:
    out = HERE / "SUMMARY.md"
    lines: list[str] = []
    lines.extend(
        (
            "# ζ₁/ζ₂/ζ₃ Prior Sensitivity Scan — Summary\n",
            "Scan date: generated by `examples/data/d2_prior_scan/zeta_scan/_analyze.py`  \nRecipe: `docs/tex/pgt_stability_priors.tex` sec. Per-model sensitivity scan recipe  \nTemplate: `examples/data/d2_prior_scan/chi_scan/_analyze.py`\n",
            "## Per-ζᵢ Thresholds\n",
            "| Parameter | ζᵢ∗_empirical | ζᵢ∗_prior (×0.5) | Prior range |",
            "|-----------|------------|--------------|-------------|",
        )
    )
    for zi in ("zeta1", "zeta2", "zeta3"):
        thr = thresholds[zi]
        pb = prior_bounds[zi]
        sym = f"ζ{'₁₂₃'[int(zi[-1]) - 1]}"
        if thr is None:
            lines.append(
                f"| {sym} | benign (>10) | {pb:.4g} (default) | [−{pb:.4g}, +{pb:.4g}] |"
            )
        else:
            lines.append(f"| {sym} | {thr:.4g} | {pb:.4g} | [−{pb:.4g}, +{pb:.4g}] |")
    lines.extend(("", "## Joint Structure (δ₁, ζᵢ) at P1\n"))
    for zi, path in joint_paths.items():
        verdict = _classify_joint(path, zi)
        lines.append(f"- **{zi}**: {verdict}")
    lines.extend(("", "## Pairwise Structure (ζᵢ, ζⱼ) at P1, δ₁=0\n"))
    for (zi, zj), path in pair_paths.items():
        verdict = _classify_pairwise(path, zi, zj)
        lines.append(f"- **({zi}, {zj})**: {verdict}")
    lines.extend(("", "## Corner Verification\n"))
    rows = list(csv.DictReader(corner_path.open(encoding="utf-8")))
    lines.extend(
        (
            "| Label | δ₁ | ζ₁ | ζ₂ | ζ₃ | γ_eff | Accept |",
            "|-------|----|----|----|----|-------|--------|",
        )
    )
    all_pass = True
    for r in rows:
        status = "✓" if int(r["stable"]) == 1 else "✗"
        if int(r["stable"]) == 0:
            all_pass = False
        lines.append(
            f"| {r['label']} | {r['delta1']} | {r['zeta1']} | {r['zeta2']} | {r['zeta3']}"
            f" | {r['max_excess']} | {status} |"
        )
    lines.append("")
    verdict = (
        "All 9 prior-corner points are stable."
        if all_pass
        else "WARNING: some corners are unstable — tighten priors."
    )
    lines.extend((f"**{verdict}**\n", "## Recommended Priors and Constraints\n", "```"))
    for zi in ("zeta1", "zeta2", "zeta3"):
        pb = prior_bounds[zi]
        lines.append(f'--prior "{zi}=uniform:-{pb:.4g}:{pb:.4g}"')
    lines.extend(
        (
            "```",
            "\nNo additional `--constraint` clauses are expected (rectangular joint structure, verified by joint + pairwise scans).",
        )
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def main() -> None:
    spec = load_spec_cached(SPEC_PATH)
    print(f"Loaded spec: {SPEC_PATH.name}")
    print(
        f"Single-param scan: {len(SINGLE_PARAM_VALUES)} values × 5 anchors × 3 ζ = "
        f"{len(SINGLE_PARAM_VALUES) * 5 * 3} probe calls"
    )

    # --- Scan 1: single-parameter ---
    print("\n[1/4] Running single-parameter scans ...")
    single_paths = run_single_param(spec)
    for path in single_paths.values():
        n = sum(1 for _ in path.open()) - 1
        print(f"  wrote {path.name} ({n} rows)")

    # Derive thresholds
    thresholds: dict[str, float | None] = {}
    for zi_name, path in single_paths.items():
        thresholds[zi_name] = _find_threshold(path)
    prior_bounds: dict[str, float] = {
        zi: (thr / 2.0 if thr is not None else 1.0) for zi, thr in thresholds.items()
    }

    # --- Scan 2: joint (δ₁, ζᵢ) ---
    print("\n[2/4] Running joint (δ₁, ζᵢ) scans at P1 ...")
    joint_paths = run_joint_scans(spec, thresholds)
    for path in joint_paths.values():
        n = sum(1 for _ in path.open()) - 1
        print(f"  wrote {path.name} ({n} rows)")

    # --- Scan 3: pairwise (ζᵢ, ζⱼ) ---
    print("\n[3/4] Running pairwise (ζᵢ, ζⱼ) scans at P1, δ₁=0 ...")
    pair_paths = run_pairwise_scans(spec, thresholds)
    for path in pair_paths.values():
        n = sum(1 for _ in path.open()) - 1
        print(f"  wrote {path.name} ({n} rows)")

    # --- Scan 4: corner verification ---
    print("\n[4/4] Running corner verification ...")
    corner_path = run_corner_verification(spec, prior_bounds)
    n = sum(1 for _ in corner_path.open()) - 1
    print(f"  wrote {corner_path.name} ({n} rows)")

    # --- Summary ---
    _summarize(
        joint_paths,
        pair_paths,
        corner_path,
        thresholds,
        prior_bounds,
    )

    # --- SUMMARY.md ---
    md_path = _write_summary_md(
        joint_paths,
        pair_paths,
        corner_path,
        thresholds,
        prior_bounds,
    )
    print(f"\n  wrote {md_path.name}")


if __name__ == "__main__":
    main()
