#!/usr/bin/env python3
"""Generate amplification heatmap with stability boundary hatching.

Produces a 2D heatmap of log₁₀(A) = log₁₀(C₀/C_EM) where:
  A > 1 (red): torsion amplifies the Gertsenshtein effect
  A = 1 (white): standard GR baseline
  A < 1 (blue): torsion suppresses the Gertsenshtein effect
  Hatched grey: unstable/diverged parameters (instability region)

Usage:
  python plot_amplification_heatmap.py <results_dir> [output.png]
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import contextlib

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch


def main() -> None:
    results_dir = Path(sys.argv[1])
    output = sys.argv[2] if len(sys.argv) > 2 else str(results_dir / "heatmap_A.png")

    csv_path = results_dir / "results.csv"
    with Path(csv_path).open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    # Detect swept parameters (exactly 2 expected)
    # Read sweep.json for parameter names
    import json

    with Path(results_dir / "sweep.json").open(encoding="utf-8") as f:
        sweep_meta = json.load(f)

    swept = sweep_meta.get("swept_parameters", sweep_meta.get("swept_params", {}))
    param_names = list(swept.keys())
    if len(param_names) != 2:
        print(f"Expected 2 swept params, got {len(param_names)}: {param_names}")
        sys.exit(1)

    p1_name, p2_name = param_names

    # Extract B0 from fixed params
    fixed = {}
    for r in rows[:1]:
        for k, v in r.items():
            if k not in {
                p1_name,
                p2_name,
                "P_max",
                "P_final",
                "P_max_time",
                "run_status",
                "error_message",
                "solver_exit_code",
                "wall_time_s",
                "grid_shape",
                "t_end",
                "dt",
                "scheme",
                "bc",
                "C0",
                "A",
                "log10_A",
                "P_valid",
            }:
                with contextlib.suppress(ValueError, TypeError):
                    fixed[k] = float(v)
    B0 = fixed.get("B0", 0.0001)
    kappa = fixed.get("kappa", 1.0)
    t_end = fixed.get("t_end", 50.0)
    P_EM = math.sin(kappa * B0 * t_end / 2) ** 2

    # Build grid
    p1_vals = sorted({float(r[p1_name]) for r in rows})
    p2_vals = sorted({float(r[p2_name]) for r in rows})
    p1_idx = {v: i for i, v in enumerate(p1_vals)}
    p2_idx = {v: i for i, v in enumerate(p2_vals)}
    n1, n2 = len(p1_vals), len(p2_vals)

    log10_A_grid = np.full((n2, n1), np.nan)
    diverged_grid = np.zeros((n2, n1), dtype=bool)

    for r in rows:
        i = p2_idx[float(r[p2_name])]
        j = p1_idx[float(r[p1_name])]

        status = r.get("run_status", "success")
        pmax_str = r.get("P_max", "")

        if status != "success" or not pmax_str:
            diverged_grid[i, j] = True
            continue

        pmax = float(pmax_str)
        if pmax <= 0 or pmax >= 0.1:
            # Invalid: either zero or outside linear regime
            diverged_grid[i, j] = True
            continue

        A = pmax / P_EM
        if A > 0:
            log10_A_grid[i, j] = math.log10(A)

    valid = log10_A_grid[np.isfinite(log10_A_grid)]
    n_valid = len(valid)
    n_diverged = int(diverged_grid.sum())
    print(f"Grid: {n1}×{n2} = {n1 * n2} cells")
    print(f"Valid: {n_valid}, Diverged/invalid: {n_diverged}")
    if n_valid > 0:
        print(f"log10(A) range: {valid.min():.1f} to {valid.max():.1f}")
        print(f"A range: {10 ** valid.min():.2e} to {10 ** valid.max():.2e}")

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10, 7))

    # Clamp color range for readability
    vmin = max(valid.min(), -8) if n_valid > 0 else -5
    vmax = min(valid.max(), 5) if n_valid > 0 else 5
    # Ensure center is between vmin and vmax
    if vmin >= 0:
        vmin = -0.5
    if vmax <= 0:
        vmax = 0.5
    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    # Main heatmap (valid data)
    im = ax.pcolormesh(
        p1_vals,
        p2_vals,
        log10_A_grid,
        shading="nearest",
        cmap="RdBu_r",
        norm=norm,
    )

    # Hatched overlay for diverged/invalid regions
    # Create a masked array: 1 where diverged, NaN where valid
    hatch_grid = np.where(diverged_grid, 1.0, np.nan)
    ax.pcolormesh(
        p1_vals,
        p2_vals,
        hatch_grid,
        shading="nearest",
        cmap="Greys",
        vmin=0,
        vmax=2,
        alpha=0.3,
    )
    # Add cross-hatching via contourf with hatch pattern
    # Use a filled contour at the diverged boundary
    X, Y = np.meshgrid(p1_vals, p2_vals)
    ax.contourf(
        X,
        Y,
        diverged_grid.astype(float),
        levels=[0.5, 1.5],
        colors="none",
        hatches=["///"],
    )

    # Labels
    ax.set_xlabel(f"{p1_name}", fontsize=12)
    ax.set_ylabel(f"{p2_name}", fontsize=12)

    # Build title from fixed params
    fixed_str = ", ".join(
        f"{k}={v:g}"
        for k, v in sorted(fixed.items())
        if k not in {p1_name, p2_name, "B0", "kappa"}
        and k in {"alpha1", "alpha2", "alpha3", "xi"}
    )
    ax.set_title(
        f"log₁₀(A): Amplification of Gertsenshtein Effect\n"
        f"({fixed_str}, κ={kappa:g}, B₀={B0:g})",
        fontsize=13,
    )

    # Colorbar
    cb = fig.colorbar(im, ax=ax, label="log₁₀(A) = log₁₀(C₀/C_EM)", extend="both")
    cb.ax.axhline(0, color="k", linewidth=0.8)

    # Legend for hatching
    handles = [
        Patch(
            facecolor="white", edgecolor="grey", hatch="///", label="Unstable / invalid"
        ),
    ]
    ax.legend(handles=handles, loc="lower left", fontsize=9, framealpha=0.9)

    plt.tight_layout()
    plt.savefig(output, dpi=150)
    print(f"Saved: {output}")


if __name__ == "__main__":
    main()
