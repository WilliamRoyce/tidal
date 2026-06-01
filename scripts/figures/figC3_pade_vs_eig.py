r"""Figure C.3 — Padé (scipy.linalg.expm) vs eigendecomposition per-mode
wall time across a resolution sweep, for three representative modal-eligible
theories.

Layout: 1×3 row of subplots (one per theory), spanning the wide-float
\\textwidth via \\begin{figure*}. Within each panel the Padé path is
solid blue with filled circles and the eigendecomposition path is dashed
orange with open squares; the y-axis is per-mode wall time on a log scale
and the x-axis is the grid size N (also log).

Data source:  benchmark_results/canonical/pade_vs_eig.json
Output:       manuscript/figures/figC3_pade_vs_eig.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex (tab:PadeBenchmark)
"""

from __future__ import annotations

import argparse
import json
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
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "pade_vs_eig.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC3_pade_vs_eig.pdf"

# Wide float — matches PRD figure* width (2× \columnwidth).
_FIG_WIDTH = 6.75
_FIG_HEIGHT = 2.3

_THEORY_ORDER = ["coupled_scalars", "gertsenshtein", "massive_gravity_3d"]
_TITLES = {
    "coupled_scalars": r"(a) Plasma--graviton (Raffelt--Stodolsky)",
    "gertsenshtein": r"(b) Einstein--Maxwell Gertsenshtein",
    "massive_gravity_3d": r"(c) Fierz--Pauli ($D{=}2{+}1$)",
}


def _n_to_float(label: str) -> float:
    """Convert a label like '256' or '32^2' to a numeric value for the x-axis."""
    if "^" in label:
        base, exp = label.split("^")
        return float(base) ** float(exp)
    return float(label)


def _plot(data: dict, out_path: Path) -> None:
    by_theory: dict[str, list[dict]] = {}
    for row in data["results"]:
        by_theory.setdefault(row["theory"], []).append(row)
    for rows in by_theory.values():
        rows.sort(key=lambda r: _n_to_float(r["n"]))

    fig, axes = plt.subplots(1, 3, figsize=(_FIG_WIDTH, _FIG_HEIGHT), sharey=True)

    pade_color = IBM_PALETTE["blue"]
    eig_color = IBM_PALETTE["orange"]

    for ax, theory in zip(axes, _THEORY_ORDER, strict=True):
        rows = by_theory.get(theory, [])
        if not rows:
            ax.text(0.5, 0.5, "no data", ha="center", va="center")
            continue
        ns = np.array([_n_to_float(r["n"]) for r in rows])
        labels = [r["n"] for r in rows]
        pade_mean = np.array([r.get("pade_s_mean", r["pade_s"]) for r in rows])
        pade_std = np.array([r.get("pade_s_std", 0.0) for r in rows])
        eig_mean = np.array([r.get("eig_s_mean", r["eig_s"]) for r in rows])
        eig_std = np.array([r.get("eig_s_std", 0.0) for r in rows])

        ax.errorbar(
            ns,
            pade_mean,
            yerr=pade_std,
            marker="o",
            color=pade_color,
            ls="none",
            ms=4,
            elinewidth=0.7,
            capsize=2.0,
            capthick=0.7,
            label=r"Pad\'e",
        )
        ax.errorbar(
            ns,
            eig_mean,
            yerr=eig_std,
            marker="s",
            mfc="none",
            color=eig_color,
            ls="none",
            ms=4,
            elinewidth=0.7,
            capsize=2.0,
            capthick=0.7,
            label="Eigendecomposition",
        )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(_TITLES[theory], fontsize=8, pad=4)
        ax.set_xlabel("$N$")
        # Label only the power-of-two N values; intermediate N values
        # (192, 384, 768, …) still show markers but their ticks are
        # unlabelled to avoid crowding.

        def is_pow2(v):
            return v > 0 and (round(v) & (round(v) - 1)) == 0

        def is_pow2_sq(v):
            return is_pow2(round(v**0.5))  # 2D theories

        tick_filter = is_pow2_sq if "2+1" in _TITLES[theory] else is_pow2
        major_ticks = [n for n in ns if tick_filter(round(n))]
        major_labels = [
            lbl for n, lbl in zip(ns, labels, strict=False) if tick_filter(round(n))
        ]
        ax.set_xticks(major_ticks)
        ax.set_xticklabels(
            [f"${round(_n_to_float(lbl))}$" for lbl in major_labels], fontsize=7
        )
        # Intermediate N values: keep them as minor ticks (no label).
        minor_ticks = [n for n in ns if not tick_filter(round(n))]
        ax.set_xticks(minor_ticks, minor=True)
        ax.tick_params(which="major", direction="in", top=False, right=False)
        ax.tick_params(which="minor", direction="in", top=False, right=False, length=2)
        ax.yaxis.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
        ax.xaxis.grid(visible=True, which="major", ls=":", alpha=0.35, lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)

    axes[0].set_ylabel(r"Per-mode wall time (s)")

    # Single shared legend above the row.
    handles, labels_legend = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels_legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 1.0),
        ncols=2,
        fontsize=8,
        frameon=False,
        handlelength=2.0,
        handletextpad=0.5,
        columnspacing=1.2,
    )

    fig.tight_layout(pad=0.2, rect=(0, 0, 1, 0.94))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    with args.data.open() as fh:
        data = json.load(fh)
    _plot(data, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
