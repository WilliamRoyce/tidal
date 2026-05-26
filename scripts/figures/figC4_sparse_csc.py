"""Figure C.4 — Sparse-CSC vs dense `expm_multiply` wall time on the
localised-Gertsenshtein convolution matrix.

Single log-log panel. Dense and sparse paths are plotted as separate
series; the vertical gap between them at any N is the speedup.
Dense is capped at N=1024 (matrix exceeds 1 GB above that and
expm_multiply becomes infeasible); sparse continues to N=4096.

Data source:  benchmark_results/canonical/sparse_csc_vs_dense.json
Output:       manuscript/figures/figC4_sparse_csc.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex (tab:SparseCSCBenchmark)
"""

from __future__ import annotations

import argparse
import json
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
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "sparse_csc_vs_dense.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC4_sparse_csc.pdf"

_FIG_WIDTH = 3.375
_FIG_HEIGHT = 2.65


def _plot(data: dict, out_path: Path) -> None:
    dense = sorted(
        [r for r in data["results"] if r["path"] == "dense"],
        key=operator.itemgetter("n"),
    )
    sparse = sorted(
        [r for r in data["results"] if r["path"] == "sparse"],
        key=operator.itemgetter("n"),
    )

    fig, ax = plt.subplots(figsize=(_FIG_WIDTH, _FIG_HEIGHT))

    if dense:
        ns_d = np.array([r["n"] for r in dense])
        t_d = np.array([r.get("wall_s_mean", r["wall_s"]) for r in dense])
        s_d = np.array([r.get("wall_s_std", 0.0) for r in dense])
        ax.errorbar(
            ns_d,
            t_d,
            yerr=s_d,
            marker="o",
            color=IBM_PALETTE["blue"],
            ls="none",
            ms=4,
            elinewidth=0.7,
            capsize=2.0,
            capthick=0.7,
            label="Dense",
        )
    if sparse:
        ns_s = np.array([r["n"] for r in sparse])
        t_s = np.array([r.get("wall_s_mean", r["wall_s"]) for r in sparse])
        s_s = np.array([r.get("wall_s_std", 0.0) for r in sparse])
        ax.errorbar(
            ns_s,
            t_s,
            yerr=s_s,
            marker="s",
            mfc="none",
            color=IBM_PALETTE["orange"],
            ls="none",
            ms=4,
            elinewidth=0.7,
            capsize=2.0,
            capthick=0.7,
            label="Sparse-CSC",
        )

    all_ns = sorted({r["n"] for r in data["results"]})
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Grid resolution $N$")
    ax.set_ylabel(r"Wall time per call (s)")
    # Major ticks at powers of two; intermediate N values appear as minor
    # ticks (no label) so the axis stays readable with dense sampling.
    pow2_ticks = [n for n in all_ns if (n & (n - 1)) == 0]
    minor_ticks = [n for n in all_ns if n not in pow2_ticks]
    ax.set_xticks(pow2_ticks)
    ax.set_xticklabels([str(n) for n in pow2_ticks])
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xlim(all_ns[0] * 0.7, all_ns[-1] * 1.4)
    ax.tick_params(which="minor", direction="in", top=False, right=False, length=2)
    ax.yaxis.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
    ax.xaxis.grid(visible=True, which="major", ls=":", alpha=0.35, lw=0.5)
    ax.legend(
        loc="upper left",
        fontsize=8,
        frameon=False,
        handlelength=2.0,
        handletextpad=0.5,
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
