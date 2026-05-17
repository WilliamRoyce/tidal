r"""Figure C.6 — Analytical-Jacobian tier speedup over the finite-difference proxy.

Per-Newton-iteration wall time of the three analytical-Jacobian delivery tiers
(dense, sparse-CSC, GMRES jactimes) compared against the finite-difference
proxy cost (n_colors × one sparse mat-vec) across a resolution sweep for three
representative theories.

Layout: 1×3 row of subplots (one per theory), spanning \textwidth via
\begin{figure*}.  Within each panel four series are plotted on log–log axes:
dense (filled circles), sparse (filled squares), GMRES (filled triangles), and
FD proxy (open diamonds, grey).  Vertical dashed lines mark the production
auto-selection boundaries at n_total = 2 000 (dense→sparse) and
n_total = 200 000 (sparse→GMRES).  Shaded tier-label bands are drawn in the
first panel only to avoid clutter.

Data source:  benchmark_results/canonical/jac_speedup.json
Output:       manuscript/figures/figC6_jac_speedup.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex (fig:AnalyticalJacobianSpeedup)
"""

from __future__ import annotations

import argparse
import json
import math
import operator
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{stix}",
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.0,
        "axes.linewidth": 0.5,
    }
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "jac_speedup.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC6_jac_speedup.pdf"

# PRD figure* width (two columns).
_FIG_WIDTH = 6.75
_FIG_HEIGHT = 2.4

_THEORY_ORDER = ["coupled_scalars", "gertsenshtein", "navier_cauchy_2d"]
_TITLES = {
    "coupled_scalars": r"(a) Plasma--graviton (Raffelt--Stodolsky)",
    "gertsenshtein": r"(b) Einstein--Maxwell Gertsenshtein",
    "navier_cauchy_2d": r"(c) Navier--Cauchy elasticity ($D{=}2{+}1$)",
}

# Production auto-selection boundaries (mirror tidal/solver/_types.py).
_DENSE_THRESHOLD = 2_000
_SPARSE_THRESHOLD = 200_000

_DENSE_COLOR = "#1f77b4"  # blue  — filled circles
_SPARSE_COLOR = "#ff7f0e"  # orange — filled squares
_GMRES_COLOR = "#2ca02c"  # green  — filled triangles
_FD_COLOR = "#888888"  # grey   — open diamonds


def _plot(data: dict, out_path: Path) -> None:
    by_theory: dict[str, list[dict]] = {}
    for row in data["results"]:
        by_theory.setdefault(row["theory"], []).append(row)
    for rows in by_theory.values():
        rows.sort(key=operator.itemgetter("n_total"))

    fig, axes = plt.subplots(1, 3, figsize=(_FIG_WIDTH, _FIG_HEIGHT))

    for ax_idx, (ax, theory) in enumerate(zip(axes, _THEORY_ORDER, strict=True)):
        rows = by_theory.get(theory, [])
        if not rows:
            ax.text(
                0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes
            )
            continue

        n_totals = np.array([r["n_total"] for r in rows])

        nan = float("nan")
        dense_mean = np.array([r.get("dense_s_mean", nan) for r in rows])
        dense_std = np.array([r.get("dense_s_std", 0.0) for r in rows])
        sparse_mean = np.array([r.get("sparse_s_mean", nan) for r in rows])
        sparse_std = np.array([r.get("sparse_s_std", 0.0) for r in rows])
        gmres_mean = np.array([r.get("gmres_s_mean", nan) for r in rows])
        gmres_std = np.array([r.get("gmres_s_std", 0.0) for r in rows])
        fd_mean = np.array([r.get("fd_s_mean", nan) for r in rows])
        fd_std = np.array([r.get("fd_s_std", 0.0) for r in rows])

        # Dense — omit NaN points.
        mask_dense = ~np.array([math.isnan(v) for v in dense_mean])

        def _plot_series(ax, x, y, yerr, mask, marker, color, mfc, label) -> None:
            xm, ym, em = x[mask], y[mask], yerr[mask]
            if xm.size == 0:
                return
            ax.errorbar(
                xm,
                ym,
                yerr=em,
                marker=marker,
                color=color,
                mfc=mfc,
                ls="none",
                ms=4,
                elinewidth=0.7,
                capsize=2.0,
                capthick=0.7,
                label=label,
                zorder=3,
            )

        all_mask = np.ones(len(rows), dtype=bool)

        _plot_series(
            ax,
            n_totals,
            dense_mean,
            dense_std,
            mask_dense,
            "o",
            _DENSE_COLOR,
            _DENSE_COLOR,
            "Dense analytical",
        )
        _plot_series(
            ax,
            n_totals,
            sparse_mean,
            sparse_std,
            all_mask,
            "s",
            _SPARSE_COLOR,
            _SPARSE_COLOR,
            "Sparse analytical",
        )
        _plot_series(
            ax,
            n_totals,
            gmres_mean,
            gmres_std,
            all_mask,
            "^",
            _GMRES_COLOR,
            _GMRES_COLOR,
            "GMRES analytical",
        )
        _plot_series(
            ax, n_totals, fd_mean, fd_std, all_mask, "D", _FD_COLOR, "none", "FD proxy"
        )

        # Axis limits — set x before drawing vertical lines.
        ax.set_xscale("log")
        ax.set_yscale("log")

        x_lo = n_totals[0] * 0.7
        x_hi = n_totals[-1] * 1.4
        ax.set_xlim(x_lo, x_hi)

        # Vertical tier-boundary lines.
        for n_thresh in (_DENSE_THRESHOLD, _SPARSE_THRESHOLD):
            if x_lo < n_thresh < x_hi:
                ax.axvline(n_thresh, ls="--", lw=0.7, color="#555555", zorder=2)

        # Tier-label shaded bands only in the first panel.
        if ax_idx == 0:
            _ylo_ax, _yhi_ax = ax.get_ylim()
            # Recompute from data to get reasonable bounds.
            all_vals = np.concatenate(
                [
                    dense_mean[mask_dense],
                    sparse_mean,
                    gmres_mean,
                    fd_mean,
                ]
            )
            all_vals = all_vals[~np.isnan(all_vals)]
            if all_vals.size:
                all_vals.min() * 0.3
                all_vals.max() * 3.0
            else:
                _ylo_data, _yhi_data = 1e-9, 1e-3

            band_alpha = 0.07
            band_color = "#000000"
            # Dense region
            ax.axvspan(
                x_lo, _DENSE_THRESHOLD, alpha=band_alpha, color=band_color, zorder=1
            )
            # Sparse region
            ax.axvspan(
                _DENSE_THRESHOLD,
                min(_SPARSE_THRESHOLD, x_hi),
                alpha=band_alpha * 0.5,
                color=band_color,
                zorder=1,
            )
            # GMRES region
            if x_hi > _SPARSE_THRESHOLD:
                ax.axvspan(
                    _SPARSE_THRESHOLD,
                    x_hi,
                    alpha=band_alpha,
                    color=band_color,
                    zorder=1,
                )

            # Tier text annotations at the top of the panel.
            for n_mid, label in [
                (min(_DENSE_THRESHOLD * 0.5, x_hi * 0.9), "auto: dense"),
                (
                    math.sqrt(_DENSE_THRESHOLD * min(_SPARSE_THRESHOLD, x_hi)),
                    "auto: sparse",
                ),
                (min(_SPARSE_THRESHOLD * 2.0, x_hi * 0.9), "auto: GMRES"),
            ]:
                if x_lo < n_mid < x_hi:
                    ax.text(
                        n_mid,
                        1.0,
                        label,
                        transform=ax.get_xaxis_transform(),
                        ha="center",
                        va="bottom",
                        fontsize=6.5,
                        color="#333333",
                        rotation=0,
                    )

        ax.set_title(_TITLES[theory], fontsize=8, pad=4)
        ax.set_xlabel(r"$n_{\rm total}$")

        ax.tick_params(which="major", direction="in", top=False, right=False)
        ax.tick_params(which="minor", direction="in", top=False, right=False, length=2)
        ax.yaxis.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
        ax.xaxis.grid(visible=True, which="major", ls=":", alpha=0.35, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel(r"Per-call wall time (s)")

    # Single shared legend above the row.
    handles, labels_leg = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels_leg,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncols=4,
        fontsize=8,
        frameon=False,
        handlelength=1.5,
        handletextpad=0.4,
        columnspacing=1.0,
    )

    fig.tight_layout(pad=0.3, rect=(0, 0, 1, 0.93))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)
    print(f"wrote {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    with args.data.open() as fh:
        data = json.load(fh)
    _plot(data, args.out)


if __name__ == "__main__":
    main()
