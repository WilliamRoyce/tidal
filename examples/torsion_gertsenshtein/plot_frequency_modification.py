#!/usr/bin/env python3
"""Fig 7: Frequency modification — P(t) and A(t) time-series.

This is the one publication figure that cannot be generated via CLI because
it requires loading raw NPY simulation output arrays and computing P(t)
manually from field energies.

Usage:
  python plot_frequency_modification.py [TS_DIR] [OUTPUT]
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Constants ---
KAPPA = 1.0
B0 = 0.001
TS_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_tests")  # noqa: S108
OUTPUT = (
    Path(sys.argv[2]) if len(sys.argv) > 2 else Path("fig7_frequency_modification.png")
)


def load_timeseries(dirname: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load h_5 and a_1 time series from a simulation output directory."""
    d = TS_DIR / dirname
    return np.load(d / "times.npy"), np.load(d / "h_5.npy"), np.load(d / "a_1.npy")


def compute_p(dirname: str, dx: float = 100.0 / 256) -> tuple[np.ndarray, np.ndarray]:
    """Compute conversion probability P(t) = E_a1 / (E_h5 + E_a1)."""
    times, h5, a1 = load_timeseries(dirname)
    e_h5 = 0.5 * np.sum(h5**2, axis=tuple(range(1, h5.ndim))) * dx
    e_a1 = 0.5 * np.sum(a1**2, axis=tuple(range(1, a1.ndim))) * dx
    total = e_h5 + e_a1
    p_vals = np.where(total > 1e-20, e_a1 / total, 0.0)
    return times, p_vals


def main() -> None:
    fig, (ax_pt, ax_at) = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(
        "Physical Effect: Frequency Modification, Not Amplitude Amplification",
        fontsize=14,
        fontweight="bold",
    )

    # GR analytical baseline
    t_an = np.linspace(0, 200, 1000)
    p_gr_an = np.sin(KAPPA * B0 * t_an / 2) ** 2

    sims = [
        ("v2_ts_d1_0", r"$\delta_1 = 0$ (GR)", "#2E7D32"),
        ("v2_d1_11p5", r"$\delta_1 = 11.5$", "#1565C0"),
        ("v2_d1_64p5", r"$\delta_1 = 64.5$", "#D84315"),
    ]

    # Panel A: P(t)
    ax_pt.plot(
        t_an,
        p_gr_an,
        color="#9E9E9E",
        linestyle="--",
        linewidth=1.5,
        label="GR analytical",
        zorder=1,
    )

    for dirname, label, color in sims:
        try:
            times, p_vals = compute_p(dirname)
            ax_pt.plot(
                times, p_vals, color=color, linewidth=1.5, label=label, alpha=0.85,
            )
        except FileNotFoundError:
            print(f"  Warning: {dirname} not found, skipping")

    ax_pt.set_xlabel("t", fontsize=11)
    ax_pt.set_ylabel(r"P(t) = E$_{a_1}$ / (E$_{h_5}$ + E$_{a_1}$)", fontsize=11)
    ax_pt.set_title("Conversion Probability P(t)", fontsize=12, fontweight="bold")
    ax_pt.legend(fontsize=9, loc="upper left")
    ax_pt.set_xlim(0, 200)
    ax_pt.set_ylim(-0.005, None)

    # Annotate GR first max
    t_gr_max = math.pi / (KAPPA * B0)
    ax_pt.annotate(
        f"GR 1st max: t={t_gr_max:.0f}",
        xy=(200, p_gr_an[-1]),
        xytext=(150, 0.85),
        fontsize=8,
        color="#9E9E9E",
        arrowprops={"arrowstyle": "->", "color": "#9E9E9E"},
    )

    # Panel B: A(t) = P(t) / P_GR(t)
    for dirname, label, color in sims[1:]:
        try:
            times, p_vals = compute_p(dirname)
            p_gr_at_t = np.sin(KAPPA * B0 * times / 2) ** 2
            valid = p_gr_at_t > 1e-15
            a_ratio = np.full_like(p_vals, np.nan)
            a_ratio[valid] = p_vals[valid] / p_gr_at_t[valid]
            ax_at.plot(
                times[valid],
                a_ratio[valid],
                color=color,
                linewidth=1.5,
                label=label,
                alpha=0.85,
            )
        except FileNotFoundError:
            pass

    ax_at.axhline(
        y=1.0,
        color="#2E7D32",
        linestyle="--",
        linewidth=1.5,
        label="GR (A=1)",
        alpha=0.6,
    )
    ax_at.set_xlabel("t", fontsize=11)
    ax_at.set_ylabel(r"A(t) = P(t) / P$_{GR}$(t)", fontsize=11)
    ax_at.set_title("Amplification Ratio A(t)", fontsize=12, fontweight="bold")
    ax_at.legend(fontsize=9)
    ax_at.set_xlim(0, 200)
    ax_at.set_yscale("log")
    ax_at.set_ylim(0.1, 300)
    ax_at.text(
        0.97,
        0.97,
        "A is NOT constant:\ncoupling modifies\noscillation frequency,\nnot amplitude cap",
        transform=ax_at.transAxes,
        ha="right",
        va="top",
        fontsize=8,
        bbox={"boxstyle": "round,pad=0.4", "facecolor": "#FFF3E0", "alpha": 0.9},
    )

    plt.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUTPUT, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUTPUT}")


if __name__ == "__main__":
    main()
