"""Figure C.2 — FD-stencil and spectral convergence on the smooth gradient test.

Plots L_2 error vs N on log-log axes for FD orders 2, 4, 6 and the FFT-based
pseudo-spectral path, on the test problem f(x) = sin(x) over a periodic
[0, 2pi] grid. Dashed guide lines show the theoretical algebraic rates
N^{-2}, N^{-4}, N^{-6}.

Data source:  benchmark_results/canonical/fd_convergence.json
Output:       manuscript/figures/figC2_fd_convergence.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex:236
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Match the manuscript: revtex4-2/PRD, \usepackage{stix}, single-column 3.375 in.
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
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "fd_convergence.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC2_fd_convergence.pdf"

# PRD column width = 3.375 in (85 mm). Height chosen for readability on log-log.
_FIG_WIDTH = 3.375
_FIG_HEIGHT = 2.65

SCHEME_STYLE = {
    "fd_o2": {"label": "FD order 2", "marker": "o", "color": "#1f77b4"},
    "fd_o4": {"label": "FD order 4", "marker": "s", "color": "#ff7f0e"},
    "fd_o6": {"label": "FD order 6", "marker": "^", "color": "#2ca02c"},
    "spectral": {"label": "Spectral", "marker": "D", "color": "#d62728"},
}


def _group_by_scheme(rows: list[dict]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    out: dict[str, list[tuple[int, float]]] = {}
    for row in rows:
        out.setdefault(row["scheme"], []).append((row["n"], row["l2_error"]))
    grouped: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for scheme, pairs in out.items():
        pairs.sort()
        ns = np.array([p[0] for p in pairs])
        errs = np.array([p[1] for p in pairs])
        grouped[scheme] = (ns, errs)
    return grouped


def _plot(data: dict, out_path: Path) -> None:
    grouped = _group_by_scheme(data["results"])

    fig, ax = plt.subplots(figsize=(_FIG_WIDTH, _FIG_HEIGHT))

    # Guide lines span the full data range (N=16 to N_max).
    n_values = sorted({row["n"] for row in data["results"]})
    n_guide = np.array([n_values[0], n_values[-1]])

    for order, color in [(2, "#1f77b4"), (4, "#ff7f0e"), (6, "#2ca02c")]:
        scheme = f"fd_o{order}"
        if scheme not in grouped:
            continue
        ns, errs = grouped[scheme]
        anchor_idx = len(ns) // 2
        anchor_n = ns[anchor_idx]
        anchor_err = errs[anchor_idx]
        guide = anchor_err * (anchor_n / n_guide) ** order
        ax.plot(n_guide, guide, ls="--", color=color, alpha=0.45, lw=0.8)
        # Label the slope at the geometric midpoint of the guide line.
        mid_n = float(np.sqrt(n_guide[0] * n_guide[-1]))
        mid_e = float(anchor_err * (anchor_n / mid_n) ** order)
        ax.text(
            mid_n,
            mid_e * 2.2,
            rf"$N^{{-{order}}}$",
            color=color,
            fontsize=7,
            ha="center",
            va="bottom",
        )

    for scheme, (ns, errs) in grouped.items():
        style = SCHEME_STYLE[scheme]
        mask = errs > 0  # hide round-off floor so log plot stays tidy
        ax.plot(
            ns[mask],
            errs[mask],
            marker=style["marker"],
            color=style["color"],
            label=style["label"],
            lw=1.0,
            ms=4,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Grid resolution $N$")
    ax.set_ylabel(r"$L_2$ error, $\partial_x \sin x$")
    ax.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
    ax.legend(frameon=False, loc="lower left", fontsize=8, handlelength=1.5)
    ax.tick_params(which="both", direction="in", top=False, right=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.3)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
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
