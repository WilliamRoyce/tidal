#!/usr/bin/env python3
"""Add amplification columns to sweep results.csv.

Computes:
  C₀     = P_max / B₀²           (conversion coefficient)
  C_EM   = P_EM / B₀²            (standard GR baseline)
  A      = C₀ / C_EM = P_max / P_EM  (amplification factor)
  log10_A = log10(A)              (for divergent colormap, centered at 0)
  P_valid = True if P_max < 0.1   (conservative linear regime check)

Usage:
  python postprocess_sweep.py <results.csv> [B0]
  python postprocess_sweep.py examples/data/nonminimal_heatmap_d1_a2/results.csv 0.0001
"""

import csv
import math
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python postprocess_sweep.py <results.csv> [B0]")
        sys.exit(1)

    path = Path(sys.argv[1])
    B0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0001

    # Gertsenshtein baseline: P_EM = sin²(κ B₀ t_end / 2), κ=1, t_end=50
    kappa = 1.0
    t_end = 50.0
    P_EM = math.sin(kappa * B0 * t_end / 2) ** 2
    C_EM = P_EM / B0**2

    print(f"B0={B0}, P_EM={P_EM:.4e}, C_EM={C_EM:.2f}")

    with Path(path).open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    n_valid = 0
    n_amplified = 0
    max_A = 0.0

    for r in rows:
        pmax_str = r.get("P_max", "")
        try:
            pmax = float(pmax_str)
        except (ValueError, TypeError):
            pmax = float("nan")

        if P_EM > 0 and pmax == pmax:  # not NaN
            A = pmax / P_EM
            valid = 0 < pmax < 0.1
            r["C0"] = f"{pmax / B0**2:.4e}"
            r["A"] = f"{A:.6e}" if valid else ""
            r["log10_A"] = f"{math.log10(A):.4f}" if A > 0 and valid else ""
            r["P_valid"] = "true" if valid else "false"
            if valid and max_A < A:
                max_A = A
            if valid and A > 1:
                n_amplified += 1
            if 0 < pmax < 0.1:
                n_valid += 1
        else:
            r["C0"] = ""
            r["A"] = ""
            r["log10_A"] = ""
            r["P_valid"] = "false"

    fieldnames = list(rows[0].keys())
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    print(f"Processed {len(rows)} rows:")
    print(f"  Valid (P<0.1): {n_valid}")
    print(f"  Amplified (A>1): {n_amplified}")
    print(f"  Max A: {max_A:.1f}")
    print(f"  Written: {path}")


if __name__ == "__main__":
    main()
