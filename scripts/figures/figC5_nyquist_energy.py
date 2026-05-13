"""Figure C.5 — Relative energy drift |dE/E| under the Fourier modal solver,
with and without the Nyquist-zero IC projection, on two coupled-scalar
configurations.

Layout: single log-log panel, x = grid resolution N, y = |dE/E|. Two
theories distinguished by color; the projection state ("with proj" vs
"without proj") distinguished by linestyle (solid vs dashed) and marker
fill (filled vs open).

Data source:  benchmark_results/canonical/nyquist_energy.json
Output:       manuscript/figures/figC5_nyquist_energy.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex (tab:NyquistEnergy)
"""

from __future__ import annotations

import argparse
import json
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
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.0,
        "axes.linewidth": 0.5,
    }
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "nyquist_energy.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC5_nyquist_energy.pdf"

_FIG_WIDTH = 3.375
_FIG_HEIGHT = 2.65

_THEORY_STYLE = {
    "coupled_scalars": {
        "color": "#1f77b4",
        "marker": "o",
        "label": r"Plasma--graviton, $\omega_P^2{=}0$",
    },
    "coupled_scalars_massive": {
        "color": "#ff7f0e",
        "marker": "s",
        "label": r"Plasma--graviton, $\omega_P^2{=}1$",
    },
}


def _plot(data: dict, out_path: Path) -> None:
    by_theory: dict[str, list[dict]] = {}
    for row in data["results"]:
        by_theory.setdefault(row["theory"], []).append(row)
    for rows in by_theory.values():
        rows.sort(key=operator.itemgetter("n"))

    fp_eps = 2.0**-52
    fig, ax = plt.subplots(figsize=(_FIG_WIDTH, _FIG_HEIGHT))

    n_vals_all: list[int] = []
    for theory, style in _THEORY_STYLE.items():
        rows = by_theory.get(theory, [])
        if not rows:
            continue
        ns = np.array([r["n"] for r in rows])
        n_vals_all.extend(ns.tolist())
        fixed = np.array([r["abs_dE_over_E_fixed"] for r in rows])
        nofix = np.array([r["abs_dE_over_E_nofixed"] for r in rows])

        # "With proj" — filled marker (the good case)
        ax.plot(
            ns,
            fixed,
            marker=style["marker"],
            color=style["color"],
            ls="none",
            ms=4,
            label=f"{style['label']} (w/ proj.)",
        )
        # "Without proj" — open marker (the bad case)
        ax.plot(
            ns,
            nofix,
            marker=style["marker"],
            mfc="none",
            color=style["color"],
            ls="none",
            ms=4,
            label=f"{style['label']} (w/o proj.)",
        )

    # Machine-epsilon floor reference.
    n_vals_sorted = sorted(set(n_vals_all))
    ax.axhline(fp_eps, color="gray", ls=":", lw=0.6, alpha=0.6)
    ax.text(
        n_vals_sorted[-1] * 1.3,
        fp_eps * 2.5,
        r"$\varepsilon_{\rm mach}$",
        color="gray",
        fontsize=6.5,
        ha="right",
        va="bottom",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Grid resolution $N$")
    ax.set_ylabel(r"$|\Delta E / E|$")
    # Major ticks at powers of two; intermediate N values as minor ticks.
    pow2_ticks = [n for n in n_vals_sorted if (n & (n - 1)) == 0]
    minor_ticks = [n for n in n_vals_sorted if n not in pow2_ticks]
    ax.set_xticks(pow2_ticks)
    ax.set_xticklabels([str(n) for n in pow2_ticks])
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xlim(n_vals_sorted[0] * 0.7, n_vals_sorted[-1] * 1.4)
    ax.tick_params(which="minor", direction="in", top=False, right=False, length=2)
    ax.yaxis.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
    ax.xaxis.grid(visible=True, which="major", ls=":", alpha=0.35, lw=0.5)
    ax.legend(
        bbox_to_anchor=(0.5, 1.02),
        loc="lower center",
        ncols=1,
        fontsize=6.5,
        frameon=False,
        handlelength=2.0,
        handletextpad=0.4,
        borderpad=0.0,
        labelspacing=0.2,
    )
    ax.tick_params(which="both", direction="in", top=False, right=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.2)
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
