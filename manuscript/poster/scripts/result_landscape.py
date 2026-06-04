r"""Hero result: amplification-vs-suppression landscape for the poster.

Regenerates the dominant parity-even sector's amp/sup overlay corner in the
modern Cambridge poster palette (Cherry = amplification, Crest = suppression),
restricted to a few informative axes including the headline curvature-torsion
contraction. Reuses the report's anesthetic loader (`load_chains`) so the report
figures in manuscript/figures/ are untouched; only the colours/style change.

Chains: hpc_results/29682868/t7_{amp,sup}_v2 (parity-even closure).
Axis labels use schematic operator contractions (not bare parameter names).

Output: manuscript/poster/figures/result_landscape.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
REPORT_FIG_DIR = REPO_ROOT / "scripts" / "figures"
sys.path.insert(0, str(REPORT_FIG_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import CONTOUR_LEVELS, load_chains  # report anesthetic loader
from _poster_palette import CHERRY, CREST, apply_poster_style

AMP = REPO_ROOT / "hpc_results" / "29682868" / "t7_amp_v2"
SUP = REPO_ROOT / "hpc_results" / "29682868" / "t7_sup_v2"
OUT = Path(__file__).resolve().parents[1] / "figures" / "result_landscape.pdf"

# Full 18-param column list for the T7 (chi-closure) chains (column naming).
ALL_PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "xi",
    "delta1",
    "zeta1",
    "zeta2",
    "zeta3",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi5",
    "chi6",
    "chi7",
    "chi8",
    "chi9",
    "chi10",
]

# Restricted axes for poster legibility; lead with the headline R~ x grad-T
# full contraction. Labels are schematic operator contractions.
PLOT_PARAMS = ["chi1", "zeta1", "delta1", "beta3"]
PARAM_LABELS = {
    "chi1": r"$\nabla T\!\cdot\!\tilde{R}$",
    "zeta1": r"$\nabla T\!\cdot\!F$",
    "delta1": r"$\tilde{R}_{[\mu\nu]}F^{\mu\nu}$",
    "beta3": r"$T_a T^a$",
}
OVERLAY_ALPHA = 0.55


def main() -> None:
    apply_poster_style()
    OUT.parent.mkdir(parents=True, exist_ok=True)

    amp = load_chains(AMP, params=ALL_PARAMS, param_labels=PARAM_LABELS)
    sup = load_chains(SUP, params=ALL_PARAMS, param_labels=PARAM_LABELS)

    axes = amp.plot_2d(
        PLOT_PARAMS,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=CHERRY,
        alpha=OVERLAY_ALPHA,
    )
    sup.plot_2d(
        axes, kinds="kde", levels=CONTOUR_LEVELS, color=CREST, alpha=OVERLAY_ALPHA
    )

    fig = axes.iloc[0, 0].figure
    for ax in axes.values.flatten():
        if ax is not None:
            ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

    handles = [
        mpatches.Patch(color=CHERRY, alpha=OVERLAY_ALPHA, label="amplification"),
        mpatches.Patch(color=CREST, alpha=OVERLAY_ALPHA, label="suppression"),
    ]
    n = len(PLOT_PARAMS)
    axes.iloc[0, 0].legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(n - 0.5, 1.0),
        bbox_transform=axes.iloc[0, 0].transAxes,
        frameon=True,
        facecolor="white",
        edgecolor="#B5BDC8",
        framealpha=1.0,
        fontsize=20,
    )

    fig.set_size_inches(8.5, 8.0)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
