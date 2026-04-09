#!/usr/bin/env python3
"""Analyze sweep results across all phases.

Classifies points as:
  - stable: A < A_SUSPECT (default 100)
  - suspect: A_SUSPECT <= A < A_TACHYONIC
  - tachyonic: A >= A_TACHYONIC (default 1e4) or run_status == "tachyonic"

Usage:
  python analyze_sweeps.py [SWEEP_DIR] [--phase 1|2|3]
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Thresholds
A_SUSPECT = 100  # above this, likely artifact
A_TACHYONIC = 1e4  # above this, certainly tachyonic
P_GR = math.sin(1.0 * 0.01 * 50.0 / 2) ** 2  # 0.06121

BASE_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_sweeps")  # noqa: S108


def load_results(path: Path) -> list[dict]:
    csv_path = path / "results.csv"
    if not csv_path.exists():
        return []
    with csv_path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def classify(rows: list[dict]) -> dict:
    """Classify rows into stable/suspect/tachyonic."""
    stable, suspect, tachyonic, error = [], [], [], []
    for r in rows:
        if r["run_status"] == "tachyonic":
            tachyonic.append(r)
            continue
        if r["run_status"] != "success":
            error.append(r)
            continue
        try:
            A = float(r["P_max"]) / P_GR
        except (ValueError, KeyError):
            error.append(r)
            continue
        if A >= A_TACHYONIC:
            tachyonic.append(r)
        elif A >= A_SUSPECT:
            suspect.append(r)
        else:
            stable.append(r)
    return {
        "stable": stable,
        "suspect": suspect,
        "tachyonic": tachyonic,
        "error": error,
    }


# ---- Phase 1 Analysis ----
print("=" * 80)
print("PHASE 1: 1D Parameter Slices")
print("=" * 80)
params = ["beta1", "beta2", "beta3", "xi", "delta1", "chi", "zeta1", "zeta2", "zeta3"]
phase1_dir = BASE_DIR / "phase1"

for param in params:
    rows = load_results(phase1_dir / param)
    if not rows:
        continue
    c = classify(rows)
    Avals = [float(r["P_max"]) / P_GR for r in c["stable"]]
    if not Avals:
        print(
            f"{param:8s}: {len(c['stable'])}s {len(c['tachyonic'])}t — no stable data"
        )
        continue
    print(
        f"{param:8s}: {len(c['stable']):3d}s {len(c['suspect']):2d}? {len(c['tachyonic']):2d}t "
        f"| A=[{min(Avals):.4f}, {max(Avals):.4f}]"
    )

# ---- Phase 2 Analysis ----
print()
print("=" * 80)
print("PHASE 2: Pairwise 2D Interactions")
print("=" * 80)
phase2_dir = BASE_DIR / "phase2"
if phase2_dir.exists():
    for subdir in sorted(phase2_dir.iterdir()):
        if not subdir.is_dir():
            continue
        rows = load_results(subdir)
        if not rows:
            continue
        c = classify(rows)
        Avals = [float(r["P_max"]) / P_GR for r in c["stable"]]
        name = subdir.name
        if Avals:
            Amin, Amax = min(Avals), max(Avals)
            varies = abs(Amax - Amin) > 0.01
            tag = "VARIES" if varies else "null"
            print(
                f"{name:20s}: {len(c['stable']):3d}s {len(c['suspect']):2d}? "
                f"{len(c['tachyonic']):2d}t | A=[{Amin:.4f}, {Amax:.4f}] {tag}"
            )
        else:
            print(
                f"{name:20s}: {len(c['stable']):3d}s {len(c['tachyonic']):2d}t — no stable data"
            )

# ---- Phase 3 Analysis ----
print()
print("=" * 80)
print("PHASE 3: Full 9D LHS Survey")
print("=" * 80)
phase3_dir = BASE_DIR / "phase3" / "lhs1000"
rows = load_results(phase3_dir)
if rows:
    c = classify(rows)
    Avals = [float(r["P_max"]) / P_GR for r in c["stable"]]
    print(f"Total: {len(rows)} runs")
    print(
        f"  Stable: {len(c['stable'])}, Suspect: {len(c['suspect'])}, "
        f"Tachyonic: {len(c['tachyonic'])}, Error: {len(c['error'])}"
    )
    if Avals:
        print(f"  A range (stable): {min(Avals):.4f} — {max(Avals):.4f}")
        print(f"  A > 1: {sum(1 for a in Avals if a > 1)}")
        print(f"  A > 2: {sum(1 for a in Avals if a > 2)}")
        print(f"  A > 10: {sum(1 for a in Avals if a > 10)}")
        print(f"  A < 0.5: {sum(1 for a in Avals if a < 0.5)} (suppression)")
else:
    print("No Phase 3 data yet.")

# ---- Generate Phase 2 scatter plots ----
if phase2_dir.exists():
    pairs = sorted(
        [d for d in phase2_dir.iterdir() if d.is_dir() and (d / "results.csv").exists()]
    )
    varying_pairs = []
    for subdir in pairs:
        rows = load_results(subdir)
        c = classify(rows)
        Avals = [float(r["P_max"]) / P_GR for r in c["stable"]]
        if Avals and (max(Avals) - min(Avals)) > 0.01:
            varying_pairs.append(subdir)

    if varying_pairs:
        ncols = min(3, len(varying_pairs))
        nrows = (len(varying_pairs) + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
        if nrows == 1 and ncols == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = axes[np.newaxis, :]
        elif ncols == 1:
            axes = axes[:, np.newaxis]

        for idx, subdir in enumerate(varying_pairs):
            ax = axes[idx // ncols][idx % ncols]
            p1, p2 = subdir.name.split("_", 1)
            rows = load_results(subdir)
            c = classify(rows)

            # Stable points colored by log10(A)
            if c["stable"]:
                x = [float(r[p1]) for r in c["stable"]]
                y = [float(r[p2]) for r in c["stable"]]
                A = [float(r["P_max"]) / P_GR for r in c["stable"]]
                log_a = [math.log10(max(a, 1e-6)) for a in A]
                sc = ax.scatter(
                    x, y, c=log_a, cmap="RdYlBu_r", s=12, alpha=0.8, vmin=-1, vmax=2
                )
                plt.colorbar(sc, ax=ax, label="log10(A)")

            # Tachyonic points
            if c["tachyonic"]:
                xt = [float(r.get(p1, 0)) for r in c["tachyonic"]]
                yt = [float(r.get(p2, 0)) for r in c["tachyonic"]]
                ax.scatter(
                    xt, yt, c="red", s=15, marker="x", alpha=0.7, label="tachyonic"
                )
                ax.legend(fontsize=7)

            ax.set_xlabel(p1)
            ax.set_ylabel(p2)
            ax.set_title(f"{p1} x {p2}")

        # Hide empty axes
        for idx in range(len(varying_pairs), nrows * ncols):
            axes[idx // ncols][idx % ncols].set_visible(False)

        plt.suptitle("Phase 2: Pairwise — A = P/P_GR (stable points)", fontsize=12)
        plt.tight_layout()
        outpath = phase2_dir / "phase2_scatter.png"
        plt.savefig(outpath, dpi=150)
        print(f"\nScatter plot saved: {outpath}")
