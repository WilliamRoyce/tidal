r"""Hero result for the poster — the talk's chi-closure restricted corner.

Replicates the content of the talk figure
(scripts/figures/overlay_chi_closure_pair_restricted_for_talk.py): the chi-closure
sector restricted to the principled top-6 couplings, with FOUR overlaid posteriors
-- propagating amplification/suppression + the non-propagating (xi=0) control --
so the plot is not sparse.

Poster styling (differs from the talk/report):
- plain SCHEMATIC family labels (R~ grad-T, T^2 ...), not explicit contractions;
- sans fonts consistent with the poster (Lato; mathtext dejavusans; no usetex);
- light-blue (#D1F9F1) figure + panel background to blend into the poster box;
- Cambridge 4-colour scheme (prop = Cherry/Crest, control = WarmCherry/WarmCrest).

Report figures in manuscript/figures/ are untouched -- only the loaders/ranking
helper are reused.

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

from _corner_style import CONTOUR_LEVELS, load_chains
from _poster_palette import (
    CAM,
    CHERRY,
    CREST,
    WARM_CHERRY,
    WARM_CREST,
    apply_poster_style,
)
from overlay_corner_pair import _top_k_union

# Chi-closure (T7) propagating + non-propagating chains (same as the talk figure).
PROP_AMP = REPO_ROOT / "hpc_results" / "29682868" / "t7_amp_v2"
PROP_SUP = REPO_ROOT / "hpc_results" / "29682868" / "t7_sup_v2"
NP_AMP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_amp_v1"
NP_SUP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_sup_v1"
OUT = Path(__file__).resolve().parents[1] / "figures" / "result_landscape.pdf"

PROP_PARAMS = [
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
NP_PARAMS = [p for p in PROP_PARAMS if p != "xi"]

# Explicit operator-contraction labels (no parameter-name prefixes, per the
# no-parameter-name rule). The double-width Results box gives room for these.
SCHEMATIC = {
    "beta1": r"$T_{abc}T^{abc}$",
    "beta2": r"$T_{abc}T^{bac}$",
    "beta3": r"$T_a T^a$",
    "delta1": r"$\tilde{R}_{[\mu\nu]}F^{\mu\nu}$",
    "zeta1": r"$(\nabla T)\!\cdot\!F$",
    "zeta2": r"$(\nabla T)\!\cdot\!F\,'$",
    "zeta3": r"$(\nabla T)\!\cdot\!F\,''$",
    "chi1": r"$(\nabla T)\!\cdot\!\tilde{R}$",
    "chi2": r"$(\nabla T)\!\cdot\!\tilde{R}\,(acbd)$",
    "chi5": r"$\mathrm{tr}(\nabla T)\!\cdot\!\tilde{R}$",
    "chi7": r"$(\nabla T)\!\cdot\!\tilde{R}\,(cdab)$",
}
for _c in range(1, 11):
    SCHEMATIC.setdefault(f"chi{_c}", r"$(\nabla T)\!\cdot\!\tilde{R}$")

OVERLAY_ALPHA = 0.55


def main() -> None:
    apply_poster_style()
    OUT.parent.mkdir(parents=True, exist_ok=True)

    top = _top_k_union(
        [PROP_AMP / "parameter_importance.json", NP_AMP / "parameter_importance.json"],
        6,
        NP_PARAMS,
    )
    print(f"top-6 carrier params: {top}")

    prop_labels = {n: SCHEMATIC.get(n, n) for n in PROP_PARAMS}
    np_labels = {n: SCHEMATIC.get(n, n) for n in NP_PARAMS}

    prop_amp = load_chains(PROP_AMP, params=PROP_PARAMS, param_labels=prop_labels)
    prop_sup = load_chains(PROP_SUP, params=PROP_PARAMS, param_labels=prop_labels)
    np_amp = load_chains(NP_AMP, params=NP_PARAMS, param_labels=np_labels)
    np_sup = load_chains(NP_SUP, params=NP_PARAMS, param_labels=np_labels)

    axes = prop_amp.plot_2d(
        top, kinds="kde", levels=CONTOUR_LEVELS, color=CHERRY, alpha=OVERLAY_ALPHA
    )
    prop_sup.plot_2d(
        axes, kinds="kde", levels=CONTOUR_LEVELS, color=CREST, alpha=OVERLAY_ALPHA
    )
    np_amp.plot_2d(
        axes, kinds="kde", levels=CONTOUR_LEVELS, color=WARM_CHERRY, alpha=OVERLAY_ALPHA
    )
    np_sup.plot_2d(
        axes, kinds="kde", levels=CONTOUR_LEVELS, color=WARM_CREST, alpha=OVERLAY_ALPHA
    )

    fig = axes.iloc[0, 0].figure
    for ax in axes.values.flatten():
        if ax is not None:
            ax.set_facecolor(CAM["light_blue"])
            ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)

    handles = [
        mpatches.Patch(color=CHERRY, alpha=OVERLAY_ALPHA, label="amplification"),
        mpatches.Patch(color=CREST, alpha=OVERLAY_ALPHA, label="suppression"),
        mpatches.Patch(
            color=WARM_CHERRY, alpha=OVERLAY_ALPHA, label=r"amplification ($\xi=0$)"
        ),
        mpatches.Patch(
            color=WARM_CREST, alpha=OVERLAY_ALPHA, label=r"suppression ($\xi=0$)"
        ),
    ]
    n = len(top)
    axes.iloc[0, 0].legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(n - 0.4, 1.0),
        bbox_transform=axes.iloc[0, 0].transAxes,
        frameon=True,
        facecolor=CAM["light_blue"],
        edgecolor=CAM["slate2"],
        framealpha=1.0,
        fontsize=20,
    )

    fig.set_size_inches(13.0, 12.0)
    fig.savefig(OUT, bbox_inches="tight", facecolor=CAM["light_blue"])
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
