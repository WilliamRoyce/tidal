#!/usr/bin/env python3
"""Analyze Phase 1 results: 1D parameter slices.

Generates:
  1. Summary table with P_max range and amplification for each parameter
  2. Individual P_max vs parameter plots for varying parameters
  3. t_end independence verification for top amplification points

Usage:
  python analyze_phase1.py [PHASE1_DIR]
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt

PHASE1_DIR = (
    Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_sweeps/phase1")  # noqa: S108
)
PARAMS = ["beta1", "beta2", "beta3", "xi", "delta1", "chi", "zeta1", "zeta2", "zeta3"]
P_GR = math.sin(1.0 * 0.01 * 50.0 / 2) ** 2  # sin^2(kappa * B0 * t_end / 2)

print(f"P_GR(t=50) = {P_GR:.6f}")
print()


def load_sweep(param: str) -> list[dict]:
    p = PHASE1_DIR / param / "results.csv"
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


# Summary table
print("Phase 1 Summary: 1D Parameter Slices (50 arctan LHS each)")
print("=" * 80)
print(
    f"{'Param':8s} {'Success':>7s} {'Tach':>5s} {'P_min':>10s} {'P_max':>10s} "
    f"{'A_min':>8s} {'A_max':>8s} {'Effect':>10s}"
)
print("-" * 80)

varying_params = []
for param in PARAMS:
    rows = load_sweep(param)
    if not rows:
        print(f"{param:8s}: NOT FOUND")
        continue
    success = [r for r in rows if r["run_status"] == "success"]
    tach = [r for r in rows if r["run_status"] == "tachyonic"]
    pvals = [float(r["P_max"]) for r in success if r.get("P_max")]

    if not pvals:
        print(f"{param:8s}: {len(success):7d} {len(tach):5d}  NO P DATA")
        continue

    pmin, pmax = min(pvals), max(pvals)
    A_min = pmin / P_GR
    A_max = pmax / P_GR
    varies = abs(pmax - pmin) > 1e-8
    effect = "VARIES" if varies else "null"
    if varies:
        varying_params.append(param)

    print(
        f"{param:8s} {len(success):7d} {len(tach):5d} {pmin:10.6f} {pmax:10.6f} "
        f"{A_min:8.4f} {A_max:8.4f} {effect:>10s}"
    )

print()
print(
    f"Parameters with effect: {', '.join(varying_params) if varying_params else 'NONE'}"
)
print(
    f"Parameters with no effect: {', '.join(p for p in PARAMS if p not in varying_params)}"
)

# Plot varying parameters
if varying_params:
    fig, axes = plt.subplots(
        1, len(varying_params), figsize=(5 * len(varying_params), 4)
    )
    if len(varying_params) == 1:
        axes = [axes]

    for ax, param in zip(axes, varying_params, strict=False):
        rows = load_sweep(param)
        success = [r for r in rows if r["run_status"] == "success"]
        tach = [r for r in rows if r["run_status"] == "tachyonic"]

        x_s = [float(r[param]) for r in success]
        y_s = [float(r["P_max"]) / P_GR for r in success]

        ax.scatter(x_s, y_s, c="steelblue", s=15, alpha=0.7, label="stable")
        if tach:
            x_t = [float(r[param]) for r in tach]
            ax.scatter(
                x_t,
                [0] * len(x_t),
                c="red",
                s=15,
                marker="x",
                alpha=0.7,
                label="tachyonic",
            )

        ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5, label="GR baseline")
        ax.set_xlabel(param)
        ax.set_ylabel("A = P / P_GR")
        ax.set_title(f"Phase 1: {param}")
        ax.legend(fontsize=8)
        ax.set_yscale("log")

    plt.tight_layout()
    outpath = PHASE1_DIR / "phase1_varying_params.png"
    plt.savefig(outpath, dpi=150)
    print(f"\nPlot saved: {outpath}")

# Also plot ALL parameters in a compact grid
fig2, axes2 = plt.subplots(3, 3, figsize=(12, 10))
for idx, param in enumerate(PARAMS):
    ax = axes2[idx // 3][idx % 3]
    rows = load_sweep(param)
    if not rows:
        ax.set_title(f"{param}: no data")
        continue
    success = [r for r in rows if r["run_status"] == "success"]
    x = [float(r[param]) for r in success]
    y = [float(r["P_max"]) / P_GR for r in success]

    ax.scatter(x, y, c="steelblue", s=8, alpha=0.7)
    ax.axhline(y=1.0, color="gray", linestyle="--", alpha=0.5)
    ax.set_xlabel(param, fontsize=9)
    ax.set_ylabel("A", fontsize=9)
    ax.set_title(param, fontsize=10, fontweight="bold")

    varies = max(y) - min(y) > 0.01 if y else False
    if varies:
        ax.set_yscale("log")
        ax.patch.set_facecolor("#fff3e0")
    else:
        ax.set_ylim(0.9, 1.1)

plt.suptitle("Phase 1: All Parameters — A = P/P_GR vs parameter value", fontsize=12)
plt.tight_layout()
outpath2 = PHASE1_DIR / "phase1_all_params.png"
plt.savefig(outpath2, dpi=150)
print(f"Plot saved: {outpath2}")
