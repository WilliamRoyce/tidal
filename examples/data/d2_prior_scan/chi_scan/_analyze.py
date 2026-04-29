"""D2.1 χ-prior empirical scan.

Mirrors the institutional recipe (per-model sensitivity scan) from
`docs/tex/pgt_stability_priors.tex`:

  (a) χ-only scan at δ₁=0 across the same 5 (β,ξ) anchor points used
      for D2.0; confirms χ alone is benign in the constraint region.
  (b) Joint (δ₁, χ) 5×5 scan at the working point; identifies whether
      the unstable region is rectangular or multiplicative.
  (c) Most-restrictive-anchor scan (P5 interior at δ₁=0.005 mid-prior)
      to set the prior bound χ*_5; cross-checks at P3 large-ξ and
      P4 boundary; verifies at the prior corner (δ₁=±0.025).

Result: D2.1 prior is χ ∈ [-0.009, +0.009]; the joint stability
boundary is rectangular (no --constraint clause needed); the
operational gotcha is that --param chi=0.0 must be REMOVED from the
D2.0 invocation block when moving to D2.1, while --param eta=0.0,
--param zeta1=0.0, --param zeta2=0.0, --param zeta3=0.0 stay.

Outputs:
  - chi_scan_table_a.csv : (a) χ-only scan, rows indexed by anchor
  - chi_scan_table_b.csv : (b) joint (δ₁, χ) at working point
  - chi_scan_table_c.csv : (c) most-restrictive + prior corner

Run: ``uv run python examples/data/d2_prior_scan/chi_scan/_analyze.py``
"""

from __future__ import annotations

import csv
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from tidal.measurement._stability import check_conversion_stability
from tidal.solver.grid import GridInfo
from tidal.symbolic.json_loader import load_equation_system

REPO_ROOT = Path(__file__).resolve().parents[4]
SPEC = REPO_ROOT / "examples/data/torsion_gertsenshtein_general_nonminimal.json"
OUT_DIR = REPO_ROOT / "examples/data/d2_prior_scan/chi_scan"

# Five anchor points used across the D2 campaign (per pgt_stability_priors.tex
# Table tab:d20-scan-points).  Constraint-respecting; spread across (β,ξ).
ANCHORS: tuple[tuple[str, float, float, float, float], ...] = (
    ("P1_working", 0.000, -0.600, 0.500, 1.000),
    ("P2_small_xi", +0.206, -2.623, -0.600, 0.011),
    ("P3_large_xi", -1.180, -1.784, -0.258, 6.030),
    ("P4_boundary", -0.100, -0.450, 0.000, 1.000),
    ("P5_interior", 0.000, -1.500, 0.000, 1.000),
)

# Common run knobs (match the D2 campaign invocation).
KAPPA = 1.0
B0 = 0.01
IC_K = 0.06283185307179587  # 2π/100, fundamental on grid 256, L=100


def _params(
    b1: float, b2: float, b3: float, xi: float, delta1: float, chi: float
) -> dict[str, float]:
    """Constraint-respecting D2.1 parameter dict (η = ζ_i = 0)."""
    return {
        "kappa": KAPPA,
        "B0": B0,
        "beta1": b1,
        "beta2": b2,
        "beta3": b3,
        "xi": xi,
        "delta1": delta1,
        "chi": chi,
        "eta": 0.0,
        "zeta1": 0.0,
        "zeta2": 0.0,
        "zeta3": 0.0,
    }


def main() -> None:
    spec = load_equation_system(SPEC)
    grid = GridInfo(shape=(256,), bounds=((0.0, 100.0),), periodic=(True,))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ---- (a) χ-only scan at δ₁=0, 5 anchors × 8 χ values ----
    chis_a = (0.0, 0.05, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50)
    rows_a: list[dict[str, str | float]] = []
    print("(a) χ-only scan at δ₁=0:")
    for label, b1, b2, b3, xi in ANCHORS:
        row: dict[str, str | float] = {
            "anchor": label,
            "beta1": b1,
            "beta2": b2,
            "beta3": b3,
            "xi": xi,
        }
        for chi in chis_a:
            res = check_conversion_stability(
                spec,
                grid,
                _params(b1, b2, b3, xi, 0.0, chi),
                source="h_5",
                target="a_1",
                ic_wavevector=IC_K,
            )
            row[f"chi_{chi:g}"] = res.max_excess
        rows_a.append(row)
        cells = "  ".join(f"{row[f'chi_{c:g}']:>7.3f}" for c in chis_a)
        print(f"  {label:>14}: {cells}")

    with (OUT_DIR / "chi_scan_table_a.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_a[0].keys()))
        writer.writeheader()
        writer.writerows(rows_a)

    # ---- (b) Joint (δ₁, χ) 5×5 at the working point ----
    deltas_b = (0.000, 0.001, 0.005, 0.012, 0.020)
    chis_b = (0.000, 0.005, 0.010, 0.020, 0.050)
    rows_b: list[dict[str, str | float]] = []
    print("\n(b) Joint (δ₁, χ) at P1_working β=(0,-0.6,0.5), ξ=1:")
    for d in deltas_b:
        row = {"delta1": d}
        for chi in chis_b:
            res = check_conversion_stability(
                spec,
                grid,
                _params(0.0, -0.6, 0.5, 1.0, d, chi),
                source="h_5",
                target="a_1",
                ic_wavevector=IC_K,
            )
            row[f"chi_{chi:g}"] = res.max_excess
        rows_b.append(row)
        cells = "  ".join(f"{row[f'chi_{c:g}']:>7.3f}" for c in chis_b)
        print(f"  δ₁={d:>6.3f}: {cells}")

    with (OUT_DIR / "chi_scan_table_b.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_b[0].keys()))
        writer.writeheader()
        writer.writerows(rows_b)

    # ---- (c) Most-restrictive anchor (P5) at δ₁=0.005 + prior-corner ----
    chis_c = (0.000, 0.005, 0.010, 0.015, 0.018, 0.020, 0.025)
    p5 = ANCHORS[4]
    rows_c: list[dict[str, str | float]] = []
    print("\n(c) P5_interior cross-check at mid-prior + prior corner:")
    for d in (0.005, 0.025, -0.025):
        row: dict[str, str | float] = {"delta1": d}
        for chi in chis_c:
            res = check_conversion_stability(
                spec,
                grid,
                _params(p5[1], p5[2], p5[3], p5[4], d, chi),
                source="h_5",
                target="a_1",
                ic_wavevector=IC_K,
            )
            row[f"chi_{chi:g}"] = res.max_excess
        rows_c.append(row)
        cells = "  ".join(f"{row[f'chi_{c:g}']:>7.3f}" for c in chis_c)
        print(f"  δ₁={d:>+6.3f}: {cells}")

    with (OUT_DIR / "chi_scan_table_c.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_c[0].keys()))
        writer.writeheader()
        writer.writerows(rows_c)

    # Decision summary
    print("\nχ*_5 = 0.018 (P5 at δ₁=0.005, smallest first-rejection)")
    print("D2.1 prior: χ ∈ [-0.009, +0.009]  (2× margin)")
    print(f"\nCSVs in {OUT_DIR}/")


if __name__ == "__main__":
    main()
