r"""Hero result: amplification-vs-suppression landscape for the poster.

Regenerates the dominant parity-even (chi-closure / T7) sector's amp-vs-sup
overlay corner in the modern Cambridge poster palette (Cherry = amplification,
Crest = suppression), restricted to the genuinely most-informative couplings.

The restricted axis set is the principled top-6 from `_top_k_union`
(scripts/figures/overlay_corner_pair.py) -- ranked by max(amp marginal, sup
marginal, cross amp/sup KL) across the propagating + non-propagating chains. This
is the SAME ranking the manuscript's restricted chi-closure figure uses; chi1
(the full curvature-torsion contraction) ranks #1 via its amp/sup cross-distinction.

Only the propagating amp/sup posteriors are drawn (the xi=0 control is a verbal
discussion point, kept off the poster). Report figures in manuscript/figures/ are
untouched -- we only reuse the loaders.

Output: manuscript/poster/figures/result_landscape.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
REPORT_FIG_DIR = REPO_ROOT / "scripts" / "figures"
sys.path.insert(0, str(REPORT_FIG_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import CONTOUR_LEVELS, load_chains
from _poster_palette import CHERRY, CREST, apply_poster_style
from overlay_corner_pair import _top_k_union

# Chi-closure (T7) propagating amp/sup chains + the NP amp importance used by the
# principled ranking.
AMP = REPO_ROOT / "hpc_results" / "29682868" / "t7_amp_v2"
SUP = REPO_ROOT / "hpc_results" / "29682868" / "t7_sup_v2"
NP_AMP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_amp_v1"
OUT = Path(__file__).resolve().parents[1] / "figures" / "result_landscape.pdf"

# Full 18-param column list for the T7 chains (column naming for the loader).
ALL_PARAMS = [
    "beta1", "beta2", "beta3", "xi", "delta1",
    "zeta1", "zeta2", "zeta3",
    "chi1", "chi2", "chi3", "chi4", "chi5",
    "chi6", "chi7", "chi8", "chi9", "chi10",
]
NP_PARAMS = [p for p in ALL_PARAMS if p != "xi"]

# Schematic operator-contraction labels (no bare parameter names), reused from the
# talk-deck restricted figure.
LABELS = {
    "beta1": r"$T_{abc}T^{abc}$",
    "beta2": r"$T_{abc}T^{bac}$",
    "beta3": r"$T_a T^a$",
    "delta1": r"$\tilde{R}_{[\mu\nu]}F^{\mu\nu}$",
    "zeta1": r"$\nabla T\!\cdot\!F$",
    "chi1": r"$\nabla T\!\cdot\!\tilde{R}$",
    "chi2": r"$\nabla T\!\cdot\!\tilde{R}\,(acbd)$",
    "chi5": r"$\mathrm{tr}\,\nabla T\!\cdot\!\tilde{R}$",
    "chi7": r"$\nabla T\!\cdot\!\tilde{R}\,(cdab)$",
}
OVERLAY_ALPHA = 0.55


def main() -> None:
    apply_poster_style()
    OUT.parent.mkdir(parents=True, exist_ok=True)

    # Principled top-6 axes (same ranking as the manuscript restricted figure).
    top = _top_k_union(
        [AMP / "parameter_importance.json", NP_AMP / "parameter_importance.json"],
        6,
        NP_PARAMS,
    )
    print(f"top-6 carrier params: {top}")
    labels = {n: LABELS.get(n, n) for n in ALL_PARAMS}

    amp = load_chains(AMP, params=ALL_PARAMS, param_labels=labels)
    sup = load_chains(SUP, params=ALL_PARAMS, param_labels=labels)

    axes = amp.plot_2d(top, kinds="kde", levels=CONTOUR_LEVELS,
                       color=CHERRY, alpha=OVERLAY_ALPHA)
    sup.plot_2d(axes, kinds="kde", levels=CONTOUR_LEVELS,
                color=CREST, alpha=OVERLAY_ALPHA)

    fig = axes.iloc[0, 0].figure
    for ax in axes.values.flatten():
        if ax is not None:
            ax.tick_params(bottom=False, left=False,
                           labelbottom=False, labelleft=False)

    handles = [
        mpatches.Patch(color=CHERRY, alpha=OVERLAY_ALPHA, label="amplification"),
        mpatches.Patch(color=CREST, alpha=OVERLAY_ALPHA, label="suppression"),
    ]
    n = len(top)
    axes.iloc[0, 0].legend(handles=handles, loc="upper right",
                           bbox_to_anchor=(n - 0.5, 1.0),
                           bbox_transform=axes.iloc[0, 0].transAxes,
                           frameon=True, facecolor="white",
                           edgecolor="#B5BDC8", framealpha=1.0, fontsize=22)

    fig.set_size_inches(9.5, 9.0)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
