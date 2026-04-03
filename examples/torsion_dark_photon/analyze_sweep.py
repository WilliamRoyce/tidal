#!/usr/bin/env python3
"""Analyze dark photon torsion sweep results using C₀ = P/B₀² metric.

Reads results.csv from a sweep output directory, computes the conversion
coefficient C₀ = P_final/B₀² and amplification A = √(C₀/C_EM), and
generates a matplotlib plot with the Gertsenshtein baseline.

Usage:
    python analyze_sweep.py <sweep_dir> [--output <png_path>]

The primary metric is the conversion coefficient C₀, not raw P_max.
See project_sweep_measurement_strategy.md (issue #206) for rationale.

Refs:
    Hwang & Noh (2310.04150) — linearized regime validity
    Domcke & Garcia-Cely (2301.02072) — Gertsenshtein thin-magnet formula
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_results(sweep_dir: Path) -> list[dict[str, float | str]]:
    """Load results.csv and parse numeric columns."""
    csv_path = sweep_dir / "results.csv"
    if not csv_path.exists():
        print(f"ERROR: {csv_path} not found", file=sys.stderr)
        sys.exit(1)
    rows: list[dict[str, float | str]] = []
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            parsed: dict[str, float | str] = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (ValueError, TypeError):
                    parsed[k] = v
            rows.append(parsed)
    return rows


def compute_c0(rows: list[dict]) -> list[dict]:
    """Add C_0 and A columns to each row."""
    for row in rows:
        b0 = float(row.get("B0", 0.01))
        kappa = float(row.get("kappa", 1.0))
        t_end = float(row.get("t_end", 50.0))

        # P_final is the physics measurement at fixed interaction length D = t_end
        p_final = row.get("P_final")
        try:
            p_final = float(p_final)  # type: ignore[arg-type]
        except (ValueError, TypeError):
            # Fall back to P_max if P_final not available or empty
            try:
                p_final = float(row.get("P_max", "nan"))  # type: ignore[arg-type]
            except (ValueError, TypeError):
                p_final = float("nan")

        # Conversion coefficient: C₀ = P / B₀²
        c0 = p_final / b0**2 if b0 > 0 else float("nan")

        # Gertsenshtein baseline: C_EM = sin²(κ B₀ t_end / 2) / B₀²
        c_em = math.sin(kappa * b0 * t_end / 2) ** 2 / b0**2 if b0 > 0 else 0.0

        # Amplification: A = √(C₀ / C_EM)
        a_factor = math.sqrt(c0 / c_em) if c_em > 0 and c0 > 0 else float("nan")

        row["C_0"] = c0
        row["C_EM"] = c_em
        row["A"] = a_factor
        row["P_final_used"] = p_final

    return rows


def filter_successful(rows: list[dict]) -> list[dict]:
    """Keep only successful, non-diverged runs."""
    good = []
    for row in rows:
        status = row.get("run_status", "success")
        c0 = row.get("C_0", float("nan"))
        # Filter: successful AND C0 within 100x baseline (exclude diverged)
        c_em = row.get("C_EM", 1.0)
        if status == "success" and isinstance(c0, float) and c0 < 100 * c_em:
            good.append(row)
    return good


def detect_swept_params(rows: list[dict]) -> list[str]:
    """Identify which parameters were swept (have >1 unique value)."""
    if not rows:
        return []
    # Check standard sweep parameters
    candidates = ["alpha", "xi", "deltam", "B0", "kappa"]
    swept = []
    for p in candidates:
        vals = {row.get(p) for row in rows if p in row}
        if len(vals) > 1:
            swept.append(p)
    return swept


def plot_1d(rows: list[dict], param: str, output: Path, c_em: float) -> None:
    """Generate 1D C₀ vs parameter plot."""
    x = np.array([float(r[param]) for r in rows])
    c0 = np.array([float(r["C_0"]) for r in rows])
    a_vals = np.array([float(r["A"]) for r in rows])

    sort_idx = np.argsort(x)
    x, c0, a_vals = x[sort_idx], c0[sort_idx], a_vals[sort_idx]

    fig, ax1 = plt.subplots(figsize=(8, 5))

    ax1.plot(x, c0, "o-", color="tab:blue", label="C₀ = P/B₀² (TIDAL)")
    ax1.axhline(
        c_em,
        color="tab:orange",
        linestyle="--",
        label=f"C_EM = {c_em:.1f} (Gertsenshtein)",
    )
    ax1.set_xlabel(param)
    ax1.set_ylabel("C₀ = P_final / B₀²")
    ax1.set_title(f"Conversion coefficient vs {param}")
    ax1.legend(loc="upper left")

    # Secondary y-axis: amplification A
    ax2 = ax1.twinx()
    ax2.plot(
        x,
        a_vals,
        "s:",
        color="tab:green",
        alpha=0.6,
        markersize=4,
        label="A = √(C₀/C_EM)",
    )
    ax2.axhline(1.0, color="tab:green", linestyle=":", alpha=0.3)
    ax2.set_ylabel("Amplification A = κ_eff / κ", color="tab:green")
    ax2.tick_params(axis="y", labelcolor="tab:green")
    ax2.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {output}")


def plot_2d(rows: list[dict], params: list[str], output: Path, c_em: float) -> None:
    """Generate 2D C₀ heatmap."""
    p1, p2 = params[0], params[1]
    x = sorted({float(r[p1]) for r in rows})
    y = sorted({float(r[p2]) for r in rows})

    grid = np.full((len(y), len(x)), np.nan)
    for r in rows:
        xi = x.index(float(r[p1]))
        yi = y.index(float(r[p2]))
        grid[yi, xi] = float(r["C_0"])

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(x, y, grid, shading="nearest", cmap="RdBu_r")
    ax.set_xlabel(p1)
    ax.set_ylabel(p2)
    ax.set_title(f"C₀ = P/B₀² over ({p1}, {p2})")
    cbar = fig.colorbar(im, ax=ax, label="C₀ = P_final / B₀²")
    # Mark the Gertsenshtein baseline on the colorbar
    cbar.ax.axhline(c_em, color="black", linestyle="--", linewidth=1.5)

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)
    print(f"  Plot saved: {output}")


def print_summary(rows: list[dict], all_rows: list[dict], c_em: float) -> None:
    """Print analysis summary."""
    n_total = len(all_rows)
    n_good = len(rows)
    n_diverged = n_total - n_good

    if not rows:
        print(f"  All {n_total} runs diverged or failed — no stable results.")
        return

    c0_vals = [float(r["C_0"]) for r in rows]
    a_vals = [float(r["A"]) for r in rows]

    print(f"  Runs: {n_good}/{n_total} stable ({n_diverged} diverged/failed)")
    print(f"  C_EM (Gertsenshtein baseline): {c_em:.2f}")
    print(f"  C₀ range: [{min(c0_vals):.2f}, {max(c0_vals):.2f}]")
    print(f"  A range:  [{min(a_vals):.4f}, {max(a_vals):.4f}]")
    amplification_threshold = 1.01
    if max(a_vals) > amplification_threshold:
        print(f"  >>> Torsion amplification detected: A_max = {max(a_vals):.4f}")
    elif n_good == n_total:
        print("  >>> No amplification: all C₀ ≈ C_EM (torsion is transparent)")


def main() -> None:  # noqa: D103
    parser = argparse.ArgumentParser(
        description="Analyze torsion sweep results using C₀ = P/B₀² metric"
    )
    parser.add_argument("sweep_dir", type=Path, help="Sweep output directory")
    parser.add_argument("--output", type=Path, default=None, help="Output PNG path")
    args = parser.parse_args()

    print(f"=== Conversion Coefficient Analysis: {args.sweep_dir.name} ===")

    # Load and process
    all_rows = load_results(args.sweep_dir)
    all_rows = compute_c0(all_rows)
    good_rows = filter_successful(all_rows)

    # Baseline
    c_em = float(good_rows[0]["C_EM"]) if good_rows else 612.0

    # Summary
    print_summary(good_rows, all_rows, c_em)

    # Detect dimensionality
    swept = detect_swept_params(all_rows)
    if not swept:
        print("  No swept parameters detected.")
        return

    # Determine output path
    output = args.output or args.sweep_dir / "analysis_C0.png"

    # Plot
    if len(swept) == 1:
        if good_rows:
            plot_1d(good_rows, swept[0], output, c_em)
        else:
            print("  No stable runs to plot.")
    elif len(swept) == 2:  # noqa: PLR2004
        if good_rows:
            plot_2d(good_rows, swept, output, c_em)
        else:
            print("  No stable runs to plot.")
    else:
        print(f"  {len(swept)} swept params ({swept}) — use scatter matrix for >2D")
        # Still print C₀ values for each run
        for r in good_rows[:10]:
            vals = ", ".join(f"{p}={r[p]:.3g}" for p in swept)
            print(f"    {vals}: C₀={r['C_0']:.2f}, A={r['A']:.4f}")


if __name__ == "__main__":
    main()
