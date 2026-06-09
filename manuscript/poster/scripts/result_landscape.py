r"""Hero result for the poster — the talk's chi-closure restricted corner.

Replicates the content of the talk figure
(scripts/figures/overlay_chi_closure_pair_restricted_for_talk.py): the chi-closure
sector restricted to the principled top-6 couplings, with FOUR overlaid posteriors
-- propagating amplification/suppression + the non-propagating control -- so the
plot is not sparse.

Determinism note: anesthetic's 2-D contour plotting randomly subsamples the weighted
chains (triangular_sample_compression_2d -> np.random.choice, unseeded), so contours
otherwise vary every render. We seed np.random and pass a large `ncompress` so the
contours use (near-)all samples -> stable, smooth, true-density contours.

Poster styling (matches the talk figure):
- schematic family labels (R . grad-T, T^2);
- the talk's four DISTINCT Cambridge colours (Crest/Purple/Cherry/Indigo);
- (propagating)/(non-propagating) legend;
- light-blue (#D1F9F1) figure + panel background to blend into the poster box.

Report figures in manuscript/figures/ are untouched -- only the loaders/ranking
helper are reused.

Output: manuscript/poster/figures/result_landscape.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
REPORT_FIG_DIR = REPO_ROOT / "scripts" / "figures"
sys.path.insert(0, str(REPORT_FIG_DIR))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import CONTOUR_LEVELS, load_chains
from _poster_palette import CAM, CHERRY, CREST, INDIGO, PURPLE, apply_poster_style
from overlay_corner_pair import _top_k_union

# Chi-closure (T7) propagating + non-propagating chains (same as the talk figure).
PROP_AMP = REPO_ROOT / "hpc_results" / "29682868" / "t7_amp_v2"
PROP_SUP = REPO_ROOT / "hpc_results" / "29682868" / "t7_sup_v2"
NP_AMP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_amp_v1"
NP_SUP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_sup_v1"
OUT = Path(__file__).resolve().parents[1] / "figures" / "result_landscape.pdf"

# Reproducibility + faithful contours: fixed seed, and use ~all samples (the
# chains hold ~1300-1400 samples, so 1000 is safely below every population).
SEED = 1
NCOMPRESS = 1000

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

# Schematic family labels (calligraphic, suppressed indices) -- matches the talk
# figure. Every chi_i collapses to R . grad-T; beta_i to T^2. mathtext-safe (no \!).
SCHEMATIC = {
    "beta1": r"$\mathcal{T}^2$",
    "beta2": r"$\mathcal{T}^2$",
    "beta3": r"$\mathcal{T}^2$",
    "delta1": r"$\mathcal{R}\cdot F$",
    "zeta1": r"$\nabla\mathcal{T}\cdot F$",
    "zeta2": r"$\nabla\mathcal{T}\cdot F$",
    "zeta3": r"$\nabla\mathcal{T}\cdot F$",
}
for _c in range(1, 11):
    SCHEMATIC[f"chi{_c}"] = r"$\mathcal{R}\cdot\nabla\mathcal{T}$"

OVERLAY_ALPHA = 0.55

# Four DISTINCT Cambridge colours, matching the talk figure's assignment.
COL_PROP_AMP = CREST
COL_PROP_SUP = PURPLE
COL_NP_AMP = CHERRY
COL_NP_SUP = INDIGO


def main() -> None:
    np.random.seed(SEED)  # deterministic anesthetic subsampling
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

    kw = {
        "kinds": "kde",
        "levels": CONTOUR_LEVELS,
        "alpha": OVERLAY_ALPHA,
        "ncompress": NCOMPRESS,
    }
    axes = prop_amp.plot_2d(top, color=COL_PROP_AMP, **kw)
    prop_sup.plot_2d(axes, color=COL_PROP_SUP, **kw)
    np_amp.plot_2d(axes, color=COL_NP_AMP, **kw)
    np_sup.plot_2d(axes, color=COL_NP_SUP, **kw)

    fig = axes.iloc[0, 0].figure
    # Explicitly force the schematic axis-label size AFTER the corner is built, so
    # anesthetic/matplotlib cannot silently reset it to a default (rcParams alone
    # was unreliable here). Large for A0 viewing distance.
    LABEL_FONTSIZE = 28
    for ax in axes.values.flatten():
        if ax is not None:
            ax.set_facecolor(CAM["light_blue"])
            ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
            ax.xaxis.label.set_fontsize(LABEL_FONTSIZE)
            ax.yaxis.label.set_fontsize(LABEL_FONTSIZE)

    handles = [
        mpatches.Patch(
            color=COL_PROP_AMP, alpha=OVERLAY_ALPHA, label="amplification (propagating)"
        ),
        mpatches.Patch(
            color=COL_PROP_SUP, alpha=OVERLAY_ALPHA, label="suppression (propagating)"
        ),
        mpatches.Patch(
            color=COL_NP_AMP,
            alpha=OVERLAY_ALPHA,
            label="amplification (non-propagating)",
        ),
        mpatches.Patch(
            color=COL_NP_SUP, alpha=OVERLAY_ALPHA, label="suppression (non-propagating)"
        ),
    ]
    n = len(top)
    axes.iloc[0, 0].legend(
        handles=handles,
        loc="upper right",
        # Push the legend into the empty upper-right region (clear of the panels);
        # its right edge sits at the grid's right edge rather than overlapping the
        # diagonal/contour panels.
        bbox_to_anchor=(n + 0.5, 1.0),
        bbox_transform=axes.iloc[0, 0].transAxes,
        frameon=True,
        facecolor=CAM["light_blue"],
        edgecolor=CAM["slate2"],
        framealpha=1.0,
        fontsize=26,
    )

    fig.set_size_inches(13.0, 12.0)
    fig.savefig(OUT, bbox_inches="tight", facecolor=CAM["light_blue"])
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
