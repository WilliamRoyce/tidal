#!/usr/bin/env python3
"""Comprehensive visualization of the PGT nonminimal sweep campaign.

Generates 7 publication-quality figures from Phase 1-4 sweep results:
  Fig 1: Phase 1 overview (3x3 grid, which parameters matter)
  Fig 2: Active parameters detail (delta1, zeta2, zeta3)
  Fig 3: Phase 2 pairwise interactions (coupling x coupling)
  Fig 4: Chi toxicity (divergence rates)
  Fig 5: Phase 3 tornado + distribution
  Fig 6: Phase 3 scatter matrix (active parameters)
  Fig 7: Frequency modification P(t) time series

Usage:
  python visualize_sweep_campaign.py [DATA_DIR] [OUTPUT_DIR]

Defaults:
  DATA_DIR   = /tmp/tidal_sweeps
  OUTPUT_DIR = /tmp/tidal_sweeps/figures
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import operator

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr

# --- Constants ---
KAPPA = 1.0
B0 = 0.001
T_END = 50.0
P_GR = math.sin(KAPPA * B0 * T_END / 2) ** 2  # 6.25e-4

PARAMS = ["beta1", "beta2", "beta3", "xi", "delta1", "chi", "zeta1", "zeta2", "zeta3"]
ACTIVE_PARAMS = ["delta1", "zeta2", "zeta3"]
NULL_PARAMS = ["beta1", "beta2", "beta3", "xi", "chi", "zeta1"]

# Greek labels for plots
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "xi": r"$\xi$",
    "delta1": r"$\delta_1$",
    "chi": r"$\chi$",
    "zeta1": r"$\zeta_1$",
    "zeta2": r"$\zeta_2$",
    "zeta3": r"$\zeta_3$",
}

DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_sweeps_v2")  # noqa: S108
OUT_DIR = Path(sys.argv[2]) if len(sys.argv) > 2 else DATA_DIR / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Time-series simulation directories
TS_DIR = Path("/tmp/tidal_tests")  # noqa: S108


# --- Data Loading Helpers ---


def load_phase1(param: str) -> list[dict]:
    p = DATA_DIR / "phase1" / param / "results.csv"
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_phase2(pair_name: str) -> list[dict]:
    p = DATA_DIR / "phase2" / pair_name / "results.csv"
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def load_phase3() -> list[dict]:
    p = DATA_DIR / "phase3" / "lhs1000" / "results.csv"
    if not p.exists():
        return []
    with p.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def get_success(rows: list[dict]) -> list[dict]:
    return [r for r in rows if r.get("run_status") == "success"]


def get_a_values(rows: list[dict]) -> np.ndarray:
    return np.array([float(r["P_max"]) / P_GR for r in rows])


def load_timeseries(dirname: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load h_5 and a_1 time series, return (times, h5, a1)."""
    d = TS_DIR / dirname
    return np.load(d / "times.npy"), np.load(d / "h_5.npy"), np.load(d / "a_1.npy")


# --- Figure 1: Phase 1 Overview ---


def fig1_phase1_overview() -> None:
    fig, axes = plt.subplots(3, 3, figsize=(14, 11))
    fig.suptitle(
        "Phase 1: Which Parameters Affect Conversion?", fontsize=14, fontweight="bold"
    )

    for idx, param in enumerate(PARAMS):
        ax = axes[idx // 3][idx % 3]
        rows = get_success(load_phase1(param))
        if not rows:
            ax.set_title(f"{PARAM_LABELS[param]}: no data", fontsize=11)
            continue

        x = np.array([float(r[param]) for r in rows])
        a_vals = get_a_values(rows)
        is_active = param in ACTIVE_PARAMS

        if is_active:
            ax.set_facecolor("#FFF8E1")
            ax.scatter(
                x, a_vals, c="#D84315", s=20, alpha=0.8, zorder=3, edgecolors="none"
            )
            ax.set_yscale("log")
            ax.set_ylim(0.5, max(a_vals) * 2)
            ax.axhspan(0.9, 1.1, color="#E8F5E9", alpha=0.4, zorder=1)
            # Annotate max
            i_max = np.argmax(a_vals)
            ax.annotate(
                f"A={a_vals[i_max]:.1f}",
                (x[i_max], a_vals[i_max]),
                textcoords="offset points",
                xytext=(5, 8),
                fontsize=8,
                color="#D84315",
                fontweight="bold",
            )
        else:
            ax.set_facecolor("#F5F5F5")
            ax.scatter(x, a_vals, c="#78909C", s=12, alpha=0.6, edgecolors="none")
            ax.set_ylim(0.95, 1.05)

        ax.axhline(y=1.0, color="#2E7D32", linestyle="--", alpha=0.6, linewidth=1)
        ax.set_xlabel(PARAM_LABELS[param], fontsize=11)
        ax.set_ylabel("A = P/P$_{GR}$", fontsize=9)

        tag = "ACTIVE" if is_active else "null"
        ax.set_title(
            f"{PARAM_LABELS[param]}  [{tag}]",
            fontsize=11,
            fontweight="bold" if is_active else "normal",
        )

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    path = OUT_DIR / "fig1_phase1_overview.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 2: Active Parameters Detail ---


def fig2_active_detail() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle(
        "Active Coupling Parameters: Conversion Modification",
        fontsize=14,
        fontweight="bold",
    )

    colors_map = plt.colormaps["RdYlBu_r"]

    for ax, param in zip(axes, ACTIVE_PARAMS, strict=False):
        rows = get_success(load_phase1(param))
        if not rows:
            continue

        x = np.array([float(r[param]) for r in rows])
        a_vals = get_a_values(rows)
        log_a = np.log10(np.clip(a_vals, 1e-3, None))

        # Normalize colors: center at log10(1) = 0
        vmax = max(abs(log_a.min()), abs(log_a.max()))
        norm = plt.Normalize(-vmax, vmax)

        sc = ax.scatter(
            x,
            a_vals,
            c=log_a,
            cmap=colors_map,
            norm=norm,
            s=40,
            alpha=0.85,
            edgecolors="k",
            linewidths=0.3,
        )
        ax.set_yscale("log")
        ax.axhline(
            y=1.0,
            color="#2E7D32",
            linestyle="--",
            alpha=0.6,
            linewidth=1.5,
            label="GR baseline",
        )

        # Annotate max and min
        i_max = np.argmax(a_vals)
        i_min = np.argmin(a_vals)
        ax.annotate(
            f"max A={a_vals[i_max]:.1f}\n({param}={x[i_max]:.1f})",
            (x[i_max], a_vals[i_max]),
            textcoords="offset points",
            xytext=(8, -15),
            fontsize=8,
            color="#B71C1C",
            fontweight="bold",
            arrowprops={"arrowstyle": "->", "color": "#B71C1C", "lw": 0.8},
        )
        ax.annotate(
            f"min A={a_vals[i_min]:.3f}",
            (x[i_min], a_vals[i_min]),
            textcoords="offset points",
            xytext=(8, 8),
            fontsize=7,
            color="#1565C0",
        )

        ax.set_xlabel(PARAM_LABELS[param], fontsize=12)
        ax.set_ylabel("A = P / P$_{GR}$", fontsize=11)
        ax.set_title(PARAM_LABELS[param], fontsize=13, fontweight="bold")
        ax.legend(fontsize=8, loc="upper left")
        plt.colorbar(sc, ax=ax, label="log$_{10}$(A)", shrink=0.8)

    plt.tight_layout(rect=(0, 0, 1, 0.93))
    path = OUT_DIR / "fig2_active_detail.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 3: Phase 2 Pairwise Interactions ---


def fig3_pairwise() -> None:
    pairs = [
        ("delta1", "zeta2"),
        ("delta1", "zeta3"),
        ("zeta2", "zeta3"),
        ("delta1", "xi"),
        ("zeta2", "beta2"),
        ("delta1", "beta2"),
    ]
    ncols = 3
    nrows = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 10))
    fig.suptitle(
        "Phase 2: Pairwise Parameter Interactions", fontsize=14, fontweight="bold"
    )

    for idx, (p1, p2) in enumerate(pairs):
        ax = axes[idx // ncols][idx % ncols]
        pair_name = f"{p1}_{p2}"
        rows = get_success(load_phase2(pair_name))
        if not rows:
            ax.set_title(f"{PARAM_LABELS[p1]} x {PARAM_LABELS[p2]}: no data")
            continue

        x = np.array([float(r[p1]) for r in rows])
        y = np.array([float(r[p2]) for r in rows])
        a_vals = get_a_values(rows)
        # Filter artifacts
        mask = a_vals < 1e4
        x, y, a_vals = x[mask], y[mask], a_vals[mask]

        log_a = np.log10(np.clip(a_vals, 1e-3, None))
        vmax = max(abs(log_a.min()), abs(log_a.max()), 1.0)
        norm = plt.Normalize(-vmax, vmax)

        sc = ax.scatter(
            x,
            y,
            c=log_a,
            cmap="RdYlBu_r",
            norm=norm,
            s=15,
            alpha=0.8,
            edgecolors="none",
        )
        plt.colorbar(sc, ax=ax, label="log$_{10}$(A)", shrink=0.8)

        ax.set_xlabel(PARAM_LABELS[p1], fontsize=10)
        ax.set_ylabel(PARAM_LABELS[p2], fontsize=10)
        ax.set_title(
            f"{PARAM_LABELS[p1]} x {PARAM_LABELS[p2]}", fontsize=11, fontweight="bold"
        )

    plt.tight_layout(rect=(0, 0, 1, 0.94))
    path = OUT_DIR / "fig3_pairwise.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 4: Chi Toxicity ---


def fig4_chi_toxicity() -> None:
    phase2_dir = DATA_DIR / "phase2"
    if not phase2_dir.exists():
        print("  Skipping fig4: no Phase 2 data")
        return

    pair_names = sorted(
        d.name
        for d in phase2_dir.iterdir()
        if d.is_dir() and (d / "results.csv").exists()
    )

    names = []
    success_rates = []
    diverge_rates = []
    has_chi = []

    for pair_name in pair_names:
        rows = load_phase2(pair_name)
        if not rows:
            continue
        n_total = len(rows)
        n_success = sum(1 for r in rows if r.get("run_status") == "success")
        n_diverged = sum(1 for r in rows if r.get("run_status") == "diverged")

        p1, p2 = pair_name.split("_", 1)
        label = f"{PARAM_LABELS.get(p1, p1)} x {PARAM_LABELS.get(p2, p2)}"
        names.append(label)
        success_rates.append(100 * n_success / n_total)
        diverge_rates.append(100 * n_diverged / n_total)
        has_chi.append("chi" in pair_name)

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle(
        r"Phase 2: $\chi$ Kinetic Mixing Causes System Instability",
        fontsize=14,
        fontweight="bold",
    )

    y_pos = np.arange(len(names))
    colors = ["#D32F2F" if c else "#1976D2" for c in has_chi]

    bars = ax.barh(
        y_pos, diverge_rates, color=colors, alpha=0.85, edgecolor="white", linewidth=0.5
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Divergence Rate (%)", fontsize=11)
    ax.set_xlim(0, 105)
    ax.axvline(x=50, color="grey", linestyle=":", alpha=0.4)

    # Add percentage labels
    for bar, rate in zip(bars, diverge_rates, strict=False):
        if rate > 5:
            ax.text(
                bar.get_width() - 3,
                bar.get_y() + bar.get_height() / 2,
                f"{rate:.0f}%",
                ha="right",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )
        elif rate > 0:
            ax.text(
                bar.get_width() + 1,
                bar.get_y() + bar.get_height() / 2,
                f"{rate:.1f}%",
                ha="left",
                va="center",
                fontsize=7,
                color="#333",
            )

    # Legend
    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor="#D32F2F", label=r"$\chi$-containing pair"),
            Patch(facecolor="#1976D2", label=r"Non-$\chi$ pair"),
        ],
        loc="lower right",
        fontsize=10,
    )

    ax.invert_yaxis()
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    path = OUT_DIR / "fig4_chi_toxicity.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 5: Phase 3 Tornado + Distribution ---


def fig5_tornado_distribution() -> None:
    rows = get_success(load_phase3())
    if not rows:
        print("  Skipping fig5: no Phase 3 data")
        return

    a_vals = get_a_values(rows)
    # Filter artifacts
    clean_mask = a_vals < 1e4
    a_clean = a_vals[clean_mask]
    clean_rows = [r for r, m in zip(rows, clean_mask, strict=False) if m]

    fig, (ax_tornado, ax_hist) = plt.subplots(
        1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1]}
    )
    fig.suptitle(
        "Phase 3: 8D Parameter Space Survey (1000 LHS Samples)",
        fontsize=14,
        fontweight="bold",
    )

    # --- Tornado chart ---
    sweep_params = [
        "beta1",
        "beta2",
        "beta3",
        "xi",
        "delta1",
        "zeta1",
        "zeta2",
        "zeta3",
    ]
    impacts = []
    for param in sweep_params:
        vals = np.array([float(r[param]) for r in clean_rows])
        rho, pval = spearmanr(np.abs(vals), a_clean)
        impacts.append((param, abs(rho), pval))

    impacts.sort(key=operator.itemgetter(1))
    param_names_sorted = [i[0] for i in impacts]
    rho_sorted = [i[1] for i in impacts]
    pval_sorted = [i[2] for i in impacts]

    y_pos = np.arange(len(impacts))
    colors = [
        "#D84315" if p < 0.001 else "#FFA726" if p < 0.05 else "#BDBDBD"
        for p in pval_sorted
    ]

    ax_tornado.barh(y_pos, rho_sorted, color=colors, edgecolor="white", linewidth=0.5)
    ax_tornado.set_yticks(y_pos)
    ax_tornado.set_yticklabels(
        [PARAM_LABELS[p] for p in param_names_sorted], fontsize=11
    )
    ax_tornado.set_xlabel("Spearman |$\\rho$| with A", fontsize=11)
    ax_tornado.set_title("Parameter Sensitivity", fontsize=12, fontweight="bold")
    ax_tornado.axvline(x=0.1, color="grey", linestyle=":", alpha=0.4)

    # Annotate significance
    for y, (_, rho, pval) in zip(y_pos, impacts, strict=False):
        sig = (
            "***"
            if pval < 0.001
            else "**"
            if pval < 0.01
            else "*"
            if pval < 0.05
            else ""
        )
        if sig:
            ax_tornado.text(
                rho + 0.005,
                y,
                sig,
                va="center",
                fontsize=10,
                color="#D84315",
                fontweight="bold",
            )

    from matplotlib.patches import Patch

    ax_tornado.legend(
        handles=[
            Patch(facecolor="#D84315", label="p < 0.001"),
            Patch(facecolor="#FFA726", label="p < 0.05"),
            Patch(facecolor="#BDBDBD", label="not significant"),
        ],
        fontsize=8,
        loc="lower right",
    )

    # --- Histogram ---
    log_a = np.log10(np.clip(a_clean, 1e-6, None))
    ax_hist.hist(
        log_a, bins=40, color="#1976D2", alpha=0.8, edgecolor="white", linewidth=0.5
    )
    ax_hist.axvline(
        x=0, color="#2E7D32", linestyle="--", linewidth=2, label="A = 1 (GR)"
    )
    ax_hist.set_xlabel("log$_{10}$(A)", fontsize=11)
    ax_hist.set_ylabel("Count", fontsize=11)
    ax_hist.set_title(
        f"Amplification Distribution (n={len(a_clean)})", fontsize=12, fontweight="bold"
    )
    ax_hist.legend(fontsize=10)

    # Add text annotations
    n_near_gr = np.sum((a_clean > 0.9) & (a_clean < 1.1))
    n_gt2 = np.sum(a_clean > 2)
    n_gt10 = np.sum(a_clean > 10)
    ax_hist.text(
        0.97,
        0.95,
        f"Near GR: {n_near_gr} ({100 * n_near_gr / len(a_clean):.0f}%)\n"
        f"A > 2: {n_gt2} ({100 * n_gt2 / len(a_clean):.0f}%)\n"
        f"A > 10: {n_gt10} ({100 * n_gt10 / len(a_clean):.0f}%)",
        transform=ax_hist.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "#E3F2FD", "alpha": 0.8},
    )

    plt.tight_layout(rect=(0, 0, 1, 0.93))
    path = OUT_DIR / "fig5_tornado_distribution.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 6: Phase 3 Scatter Matrix ---


def fig6_scatter_matrix() -> None:
    rows = get_success(load_phase3())
    if not rows:
        print("  Skipping fig6: no Phase 3 data")
        return

    a_vals = get_a_values(rows)
    clean_mask = a_vals < 1e4
    clean_rows = [r for r, m in zip(rows, clean_mask, strict=False) if m]
    a_clean = a_vals[clean_mask]
    log_a = np.log10(np.clip(a_clean, 1e-6, None))

    params = ACTIVE_PARAMS  # delta1, zeta2, zeta3
    n = len(params)
    fig, axes = plt.subplots(n, n, figsize=(12, 11))
    fig.suptitle(
        "Phase 3: Active Parameter Space (colored by log$_{10}$(A))",
        fontsize=14,
        fontweight="bold",
    )

    vmax = max(abs(log_a.min()), abs(log_a.max()), 1.0)
    norm = plt.Normalize(-vmax, vmax)
    cmap = plt.colormaps["RdYlBu_r"]

    param_data = {p: np.array([float(r[p]) for r in clean_rows]) for p in params}

    for i in range(n):
        for j in range(n):
            ax = axes[i][j]
            if i == j:
                # Diagonal: marginal (parameter vs A)
                ax.scatter(
                    param_data[params[i]],
                    a_clean,
                    c=log_a,
                    cmap=cmap,
                    norm=norm,
                    s=8,
                    alpha=0.7,
                    edgecolors="none",
                )
                ax.set_yscale("log")
                ax.axhline(
                    y=1.0, color="#2E7D32", linestyle="--", alpha=0.5, linewidth=1
                )
                ax.set_ylabel("A", fontsize=9)
            else:
                # Off-diagonal: pairwise scatter
                ax.scatter(
                    param_data[params[j]],
                    param_data[params[i]],
                    c=log_a,
                    cmap=cmap,
                    norm=norm,
                    s=8,
                    alpha=0.7,
                    edgecolors="none",
                )

            if i == n - 1:
                ax.set_xlabel(PARAM_LABELS[params[j]], fontsize=10)
            if j == 0 and i != 0:
                ax.set_ylabel(PARAM_LABELS[params[i]], fontsize=10)
            elif j == 0 and i == 0:
                ax.set_ylabel("A", fontsize=10)

    # Add colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes((0.94, 0.15, 0.02, 0.7))
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cbar_ax, label="log$_{10}$(A)")

    plt.tight_layout(rect=(0, 0, 0.92, 0.94))
    path = OUT_DIR / "fig6_scatter_matrix.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Figure 7: Frequency Modification ---


def fig7_frequency_modification() -> None:
    fig, (ax_pt, ax_at) = plt.subplots(1, 2, figsize=(14, 5.5))
    fig.suptitle(
        "Physical Effect: Frequency Modification, Not Amplitude Amplification",
        fontsize=14,
        fontweight="bold",
    )

    dx = 100.0 / 256

    # Compute P(t) from simulation data
    def compute_p_from_sim(dirname: str) -> tuple[np.ndarray, np.ndarray]:
        times, h5, a1 = load_timeseries(dirname)
        p_vals = []
        for i in range(len(times)):
            e_h5 = 0.5 * np.sum(h5[i] ** 2) * dx
            e_a1 = 0.5 * np.sum(a1[i] ** 2) * dx
            total = e_h5 + e_a1
            p_vals.append(e_a1 / total if total > 1e-20 else 0.0)
        return times, np.array(p_vals)

    # GR analytical
    t_an = np.linspace(0, 200, 1000)
    p_gr_an = np.sin(KAPPA * B0 * t_an / 2) ** 2

    # Load simulation data
    sims = [
        ("v2_ts_d1_0", r"$\delta_1 = 0$ (GR)", "#2E7D32", "-"),
        ("v2_d1_11p5", r"$\delta_1 = 11.5$", "#1565C0", "-"),
        ("v2_d1_64p5", r"$\delta_1 = 64.5$", "#D84315", "-"),
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

    for dirname, label, color, ls in sims:
        try:
            times, p_vals = compute_p_from_sim(dirname)
            ax_pt.plot(
                times,
                p_vals,
                color=color,
                linestyle=ls,
                linewidth=1.5,
                label=label,
                alpha=0.85,
            )
        except FileNotFoundError:
            print(f"    Warning: {dirname} not found, skipping")

    ax_pt.set_xlabel("t", fontsize=11)
    ax_pt.set_ylabel("P(t) = E$_{a_1}$ / (E$_{h_5}$ + E$_{a_1}$)", fontsize=11)
    ax_pt.set_title("Conversion Probability P(t)", fontsize=12, fontweight="bold")
    ax_pt.legend(fontsize=9, loc="upper left")
    ax_pt.set_xlim(0, 200)
    ax_pt.set_ylim(-0.005, None)

    # Annotate first max for delta1=64.5
    try:
        times64, p64 = compute_p_from_sim("v2_d1_64p5")
        from scipy.signal import argrelmax

        maxima = argrelmax(p64, order=3)[0]
        if len(maxima) > 0:
            t_max = times64[maxima[0]]
            p_max = p64[maxima[0]]
            ax_pt.annotate(
                f"1st max: t={t_max:.0f}",
                (t_max, p_max),
                textcoords="offset points",
                xytext=(20, -10),
                fontsize=8,
                color="#D84315",
                arrowprops={"arrowstyle": "->", "color": "#D84315"},
            )
        # GR first max
        t_gr_max = math.pi / (KAPPA * B0)
        ax_pt.annotate(
            f"GR 1st max: t={t_gr_max:.0f}",
            xy=(200, p_gr_an[-1]),
            xytext=(150, 0.85),
            fontsize=8,
            color="#9E9E9E",
            arrowprops={"arrowstyle": "->", "color": "#9E9E9E"},
        )
    except FileNotFoundError:
        pass

    # Panel B: A(t) = P(t) / P_GR(t)
    for dirname, label, color, _ls in sims[1:]:  # skip GR (A=1 by definition)
        try:
            times, p_vals = compute_p_from_sim(dirname)
            p_gr_at_t = np.sin(KAPPA * B0 * times / 2) ** 2
            # Avoid division by zero at t=0
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
    ax_at.set_ylabel("A(t) = P(t) / P$_{GR}$(t)", fontsize=11)
    ax_at.set_title("Amplification Ratio A(t)", fontsize=12, fontweight="bold")
    ax_at.legend(fontsize=9)
    ax_at.set_xlim(0, 200)
    ax_at.set_yscale("log")
    ax_at.set_ylim(0.1, 300)

    # Annotate key insight
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
    path = OUT_DIR / "fig7_frequency_modification.png"
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# --- Main ---

if __name__ == "__main__":
    print("Generating sweep campaign visualizations...")
    print(f"  Data:   {DATA_DIR}")
    print(f"  Output: {OUT_DIR}")
    print(f"  P_GR = {P_GR:.6f}")
    print()

    print("Figure 1: Phase 1 Overview")
    fig1_phase1_overview()

    print("Figure 2: Active Parameters Detail")
    fig2_active_detail()

    print("Figure 3: Pairwise Interactions")
    fig3_pairwise()

    print("Figure 4: Chi Toxicity")
    fig4_chi_toxicity()

    print("Figure 5: Tornado + Distribution")
    fig5_tornado_distribution()

    print("Figure 6: Scatter Matrix")
    fig6_scatter_matrix()

    print("Figure 7: Frequency Modification")
    fig7_frequency_modification()

    print(f"\nAll figures saved to: {OUT_DIR}")
