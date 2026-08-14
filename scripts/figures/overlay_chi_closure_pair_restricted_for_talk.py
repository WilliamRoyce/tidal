r"""Talk-deck variant of the manuscript's Fig. 9.

Re-runs the chi-closure restricted overlay corner (18D parity-even closure
vs. xi=0 non-propagating-torsion control, top-6 by D_KL) for the practice
talk:

  - reuses chain paths, param lists, top-K ranking from
    scripts/figures/overlay_corner_pair.py
  - overrides labels with literal operator-contraction notation (no
    `beta_i:` / `chi_i:` prefix), drawn verbatim from the manuscript
    enumeration in research/lagrangian_enumeration/general_quadratic_lagrangian.tex
  - bumps in-figure font sizes for presentation distance
  - swaps the four overlay colours from IBM to the 2024 Cambridge palette
    (Dark Crest / Dark Purple / Dark Cherry / Dark Indigo) so the corner
    accents match the slide accents
  - saves with transparent figure background so the corner sits cleanly
    on the cream Cambridge beamer background (axes panels keep their
    white background)

Output: docs/talks/practice_2026_04/figures/overlay_chi_closure_pair_restricted_for_talk.pdf
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import FIG_WIDTH
from _palette import CAMBRIDGE_PALETTE
from overlay_corner_pair import _top_k_union, render_overlay_pair

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("overlay_chi_closure_pair_restricted_for_talk")


# Operator-contraction labels — no parameter-symbol prefix.  Each label
# shows the actual index contraction the parameter multiplies, in the
# notation used by the manuscript source (general_quadratic_lagrangian.tex
# eq:r_dt_explicit / eq:torsion_squared / eq:curv_em).
# Schematic class-name labels (not specific contractions) so the labels
# stay legible at presentation distance and the audience reads "which
# sector dominates" rather than parsing index-permutation differences.
# Multiple panels may share the same label by design.
LABEL_OVERRIDES: dict[str, str] = {
    # T^2 torsion-mass class
    "beta1": r"$\mathcal{T}^2$",
    "beta2": r"$\mathcal{T}^2$",
    "beta3": r"$\mathcal{T}^2$",
    # Non-minimal Riemann-Cartan x EM
    "delta1": r"$\mathcal{R}\!\cdot\! F$",
    # Derivative torsion x EM
    "zeta1": r"$\nabla\mathcal{T}\!\cdot\! F$",
    "zeta2": r"$\nabla\mathcal{T}\!\cdot\! F$",
    "zeta3": r"$\nabla\mathcal{T}\!\cdot\! F$",
    # Curvature x derivative-torsion class (every chi_i collapses to the same
    # schematic class).
    "chi1": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi2": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi3": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi4": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi5": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi6": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi7": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi8": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi9": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
    "chi10": r"$\mathcal{R}\!\cdot\!\nabla\mathcal{T}$",
}


# Cambridge palette overlay assignment (sixth-pass: brighter, more
# distinguishable than the dark variants at 50% alpha).
#   propagating amp  -> Crest   (orange-coral, warm)
#   propagating sup  -> Purple
#   NP-control amp   -> Cherry  (warm-pair with Crest)
#   NP-control sup   -> Indigo
CAMBRIDGE_OVERLAY_COLOURS = {
    "prop_amp": CAMBRIDGE_PALETTE["crest"],
    "prop_sup": CAMBRIDGE_PALETTE["purple"],
    "np_amp": CAMBRIDGE_PALETTE["cherry"],
    "np_sup": CAMBRIDGE_PALETTE["indigo"],
}

# Legend labels — replace internal `xi=0 control` jargon with `non-propagating`.
LEGEND_LABELS = {
    "prop_amp": "amplification (propagating)",
    "prop_sup": "suppression (propagating)",
    "np_amp": "amplification (non-propagating)",
    "np_sup": "suppression (non-propagating)",
}


# Chain paths (mirror Pair A of overlay_corner_pair.py).
CHI_PROP_PARAMS = [
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
CHI_NP_PARAMS = [
    "beta1",
    "beta2",
    "beta3",
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
CHI_PROP_AMP = REPO_ROOT / "hpc_results" / "29682868" / "t7_amp_v2"
CHI_PROP_SUP = REPO_ROOT / "hpc_results" / "29682868" / "t7_sup_v2"
CHI_NP_AMP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_amp_v1"
CHI_NP_SUP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_sup_v1"

OUT_PATH = (
    REPO_ROOT
    / "docs"
    / "talks"
    / "practice_2026_04"
    / "figures"
    / "overlay_chi_closure_pair_restricted_for_talk.pdf"
)


def main() -> None:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Bump in-figure font sizes for slide / presentation distance.
    # render_overlay_pair() calls apply_style() internally, which sets the
    # manuscript-tuned RCPARAMS; override afterwards by wrapping the call.
    # Simplest path: update rcParams now, render, then matplotlib resets per
    # process exit.  rcParams set here propagate through apply_style() because
    # apply_style() merges into existing rcParams via mpl.rcParams.update.
    chi_top = _top_k_union(
        [
            CHI_PROP_AMP / "parameter_importance.json",
            CHI_NP_AMP / "parameter_importance.json",
        ],
        6,
        CHI_NP_PARAMS,
    )
    log.info("[talk] top-6 carrier params: %s", chi_top)

    render_overlay_pair(
        out_path=OUT_PATH,
        plot_params=chi_top,
        prop_amp_dir=CHI_PROP_AMP,
        prop_sup_dir=CHI_PROP_SUP,
        prop_param_names=CHI_PROP_PARAMS,
        np_amp_dir=CHI_NP_AMP,
        np_sup_dir=CHI_NP_SUP,
        np_param_names=CHI_NP_PARAMS,
        fig_width=FIG_WIDTH,
        param_label_overrides=LABEL_OVERRIDES,
        transparent=True,
        colours=CAMBRIDGE_OVERLAY_COLOURS,
        legend_labels=LEGEND_LABELS,
        strip_ylabels=False,
        strip_xlabels=True,
        ylabel_rotation=0,
        rcparams_overrides={
            "axes.labelsize": 18,
            "legend.fontsize": 18,
        },
        legend_anchor=(7.5, 1.0),
        # Deterministic, true-density contours: anesthetic randomly subsamples
        # the weighted chains (triangular_sample_compression_2d) so contours
        # otherwise drift run-to-run.  Seed + ncompress ≈ all samples (chains
        # hold ~1300–1400) -> pixel-identical, smooth, faithful contours.
        seed=1,
        ncompress=1000,
    )


if __name__ == "__main__":
    main()
