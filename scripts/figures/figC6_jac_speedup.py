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
from _palette import IBM_PALETTE

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

_DENSE_COLOR = IBM_PALETTE["blue"]  # filled circles
_SPARSE_COLOR = IBM_PALETTE["orange"]  # filled squares
_GMRES_COLOR = IBM_PALETTE["purple"]  # filled triangles (tab:green → IBM purple)
_FD_COLOR = "#888888"  # grey — open diamonds


def _plot(data: dict, out_path: Path) -> None:
    by_theory: dict[str, list[dict]] = {}
    for row in data["results"]:
        by_theory.setdefault(row["theory"], []).append(row)
    for rows in by_theory.values():
        rows.sort(key=operator.itemgetter("n_total"))

    # Global x-range: union of all data points, extended so the
    # sparse->GMRES threshold (200K) lands inside every panel.
    all_n = np.array(
        [r["n_total"] for r in data["results"]],
        dtype=float,
    )
    if all_n.size == 0:
        all_n = np.array([1.0, _SPARSE_THRESHOLD * 2.0])
    global_x_lo = float(all_n.min()) * 0.7
    global_x_hi = max(float(all_n.max()), _SPARSE_THRESHOLD * 1.5) * 1.4

    # Global y-range: union of all four series across all theories.
    all_y_vals: list[float] = []
    for r in data["results"]:
        for key in ("dense_s_mean", "sparse_s_mean", "gmres_s_mean", "fd_s_mean"):
            v = r.get(key, float("nan"))
            if isinstance(v, (int, float)) and not math.isnan(v) and v > 0:
                all_y_vals.append(float(v))
    if all_y_vals:
        global_y_lo = min(all_y_vals) * 0.3
        global_y_hi = max(all_y_vals) * 3.0
    else:
        global_y_lo, global_y_hi = 1e-9, 1e-3

    fig, axes = plt.subplots(1, 3, figsize=(_FIG_WIDTH, _FIG_HEIGHT), sharey=True)

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

        # Shared axis limits + scales — set x before drawing vertical lines.
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(global_x_lo, global_x_hi)
        ax.set_ylim(global_y_lo, global_y_hi)

        # Vertical tier-boundary lines — guard still passes for both thresholds
        # because the global x-range was extended past _SPARSE_THRESHOLD * 1.5.
        for n_thresh in (_DENSE_THRESHOLD, _SPARSE_THRESHOLD):
            if global_x_lo < n_thresh < global_x_hi:
                ax.axvline(n_thresh, ls="--", lw=0.7, color="#555555", zorder=2)

        # Tier-zone labels: italicised, near the bottom of the centre panel
        # only.  Avoids the title-clash and shading-ambiguity issues of the
        # previous design.
        if ax_idx == 1:
            for n_mid, label in (
                (math.sqrt(global_x_lo * _DENSE_THRESHOLD), "dense"),
                (math.sqrt(_DENSE_THRESHOLD * _SPARSE_THRESHOLD), "sparse"),
                (math.sqrt(_SPARSE_THRESHOLD * global_x_hi), "GMRES"),
            ):
                if global_x_lo < n_mid < global_x_hi:
                    ax.text(
                        n_mid,
                        0.02,
                        label,
                        transform=ax.get_xaxis_transform(),
                        ha="center",
                        va="bottom",
                        fontsize=6.5,
                        style="italic",
                        color="#555555",
                    )

        ax.set_title(_TITLES[theory], fontsize=8, pad=4)
        ax.set_xlabel(r"$n_{\rm total}$")

        # Force log-scale minor ticks (2x, 3x, ..., 9x within each decade)
        # to render on both axes — matplotlib auto-suppresses these when
        # the major tick range spans more than ~6 decades, which our
        # 9-decade y-range triggers.
        ax.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=12))
        ax.xaxis.set_minor_locator(
            mpl.ticker.LogLocator(base=10, subs=tuple(range(2, 10)), numticks=120)
        )
        ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10, numticks=12))
        ax.yaxis.set_minor_locator(
            mpl.ticker.LogLocator(base=10, subs=tuple(range(2, 10)), numticks=120)
        )
        ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())

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
