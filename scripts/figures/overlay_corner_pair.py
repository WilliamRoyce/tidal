r"""Overlay corner plot for a paired-sector comparison.

Plots one corner with four overlaid posteriors:
  - amplification (propagating)   — IBM magenta
  - suppression  (propagating)    — IBM yellow
  - amplification ($\xi=0$ control) — IBM purple
  - suppression  ($\xi=0$ control)  — IBM orange

Two pairs are rendered:
  - Pair A: $\chi$-closure (18D propagating) vs NP $\chi$-closure (17D control)
  - Pair B: YM-PGT non-minimal union (9D propagating) vs NP control (8D)

For each pair, both a full-corner overlay (all 17/18 axes for A, all 8/9 for B)
and a top-K restricted-corner overlay (K=6) are rendered.

The script is standalone: it does not modify the existing `kl_carrier_corner.py`
or `overlay_corner()` pipeline. Outputs:
  - manuscript/figures/overlay_chi_closure_pair_full.pdf
  - manuscript/figures/overlay_chi_closure_pair_restricted.pdf
  - manuscript/figures/overlay_ym_union_pair_full.pdf
  - manuscript/figures/overlay_ym_union_pair_restricted.pdf
"""

from __future__ import annotations

import json
import logging
import math
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import (
    AMP_COLOR,
    COLUMN_WIDTH,
    CONTOUR_LEVELS,
    FIG_WIDTH,
    MAX_NCOMPRESS,
    OVERLAY_ALPHA,
    SUP_COLOR,
    apply_style,
    load_chains,
)
from _palette import IBM_PALETTE

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("overlay_corner_pair")

# Four-color scheme: existing AMP/SUP for propagating; IBM purple/orange for NP control.
NP_AMP_COLOR = IBM_PALETTE["purple"]
NP_SUP_COLOR = IBM_PALETTE["orange"]

PARAM_LABEL: dict[str, str] = {
    "beta1": r"\beta_1",
    "beta2": r"\beta_2",
    "beta3": r"\beta_3",
    "xi": r"\xi",
    "chi": r"\chi",
    "delta1": r"\delta_1",
    "zeta1": r"\zeta_1",
    "zeta2": r"\zeta_2",
    "zeta3": r"\zeta_3",
    "chi1": r"\chi_1",
    "chi2": r"\chi_2",
    "chi3": r"\chi_3",
    "chi4": r"\chi_4",
    "chi5": r"\chi_5",
    "chi6": r"\chi_6",
    "chi7": r"\chi_7",
    "chi8": r"\chi_8",
    "chi9": r"\chi_9",
    "chi10": r"\chi_{10}",
}


def _label(name: str) -> str:
    if name in PARAM_LABEL:
        return f"${PARAM_LABEL[name]}$"
    return rf"$\mathtt{{{name}}}$"


def _score_from_imp(imp: dict, name: str) -> float:
    amp = imp.get("amp", {}).get("marginal_d_kl", {}) or {}
    sup = imp.get("sup", {}).get("marginal_d_kl", {}) or {}
    cross = imp.get("cross_amp_sup_kl", {}) or {}
    vals = [amp.get(name, 0.0), sup.get(name, 0.0), cross.get(name, 0.0)]
    vals = [v if v is not None and math.isfinite(v) else 0.0 for v in vals]
    return max(vals)


def _top_k(amp_imp_path: Path, k: int, allowed: list[str]) -> list[str]:
    """Pick top-K parameters by combined max(amp marginal, sup marginal, cross)."""
    with amp_imp_path.open() as fh:
        imp = json.load(fh)
    allowed_set = set(allowed)
    names = [n for n in imp.get("amp", {}).get("marginal_d_kl", {}) if n in allowed_set]
    return sorted(names, key=lambda n: -_score_from_imp(imp, n))[:k]


def _top_k_union(paths: list[Path], k: int, allowed: list[str]) -> list[str]:
    """Pick top-K parameters by max of scores across multiple chains.

    For each parameter in `allowed`, score = max over all chains of that
    chain's combined max(amp_marginal, sup_marginal, cross). Returns the
    top-K names ranked by this max-score.
    """
    imps = []
    for p in paths:
        with p.open() as fh:
            imps.append(json.load(fh))
    allowed_set = set(allowed)
    names = set()
    for imp in imps:
        names.update(
            n for n in imp.get("amp", {}).get("marginal_d_kl", {}) if n in allowed_set
        )

    def score(n: str) -> float:
        return max(_score_from_imp(imp, n) for imp in imps)

    return sorted(names, key=lambda n: -score(n))[:k]


def render_overlay_pair(
    *,
    out_path: Path,
    plot_params: list[str],
    prop_amp_dir: Path,
    prop_sup_dir: Path,
    prop_param_names: list[str],
    np_amp_dir: Path,
    np_sup_dir: Path,
    np_param_names: list[str],
    fig_width: float,
    param_label_overrides: dict[str, str] | None = None,
    transparent: bool = False,
    colours: dict[str, str] | None = None,
    legend_labels: dict[str, str] | None = None,
    strip_ylabels: bool = False,
    strip_xlabels: bool = False,
    ylabel_rotation: float | None = None,
    rcparams_overrides: dict | None = None,
    legend_anchor: tuple[float, float] | None = None,
    seed: int | None = 1,
    ncompress: int | None = None,  # deprecated; dynamic per-chain ncompress wins
) -> None:
    apply_style()
    # Apply caller's rcparams overrides AFTER apply_style() so they actually
    # stick — apply_style() unconditionally updates mpl.rcParams from RCPARAMS
    # which would otherwise overwrite any rcParams set by the caller before
    # the render call.
    if rcparams_overrides:
        import matplotlib as mpl

        mpl.rcParams.update(rcparams_overrides)

    # Deterministic, full-detail contours. Anesthetic's 2-D plotting
    # randomly subsamples the weighted chains
    # (triangular_sample_compression_2d -> np.random.choice, unseeded);
    # passing ncompress = len(chain) per call makes that subsample a
    # no-op so EVERY weighted sample contributes to the 2D KDE. Reseed
    # np.random before each plot_2d as a belt-and-braces guard against
    # any other unseeded np.random consumer inside anesthetic.
    import numpy as np

    seed_value = seed if seed is not None else 1
    _ = ncompress  # accepted for backwards compat with talk-script caller; ignored in favour of dynamic per-chain value

    # Resolve four overlay colours (prop_amp, prop_sup, np_amp, np_sup),
    # defaulting to module-level constants when not overridden.
    amp = (colours or {}).get("prop_amp", AMP_COLOR)
    sup = (colours or {}).get("prop_sup", SUP_COLOR)
    np_amp_ = (colours or {}).get("np_amp", NP_AMP_COLOR)
    np_sup_ = (colours or {}).get("np_sup", NP_SUP_COLOR)

    def _resolve_label(name: str) -> str:
        if param_label_overrides and name in param_label_overrides:
            return param_label_overrides[name]
        return _label(name)

    prop_amp = load_chains(
        prop_amp_dir,
        params=prop_param_names,
        param_labels={n: _resolve_label(n) for n in prop_param_names},
    )
    prop_sup = load_chains(
        prop_sup_dir,
        params=prop_param_names,
        param_labels={n: _resolve_label(n) for n in prop_param_names},
    )
    np_amp = load_chains(
        np_amp_dir,
        params=np_param_names,
        param_labels={n: _resolve_label(n) for n in np_param_names},
    )
    np_sup = load_chains(
        np_sup_dir,
        params=np_param_names,
        param_labels={n: _resolve_label(n) for n in np_param_names},
    )

    # Canvas spans the full `plot_params` list. NP chains may lack some
    # columns (e.g. xi). When every plot_params entry is in the NP list
    # (restricted-corner case) we plot NP directly on `axes`, matching
    # the original behavior. When some columns are propagating-only
    # (full-corner case with xi) we slice axes to the NP subset so the
    # propagating-only column shows propagating contours only.
    prop_only_cols = [p for p in plot_params if p not in np_param_names]
    np_plot_params = [p for p in plot_params if p in np_param_names]

    np.random.seed(seed_value)  # noqa: NPY002 -- must set the legacy global RNG anesthetic reads
    axes = prop_amp.plot_2d(
        plot_params,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=amp,
        alpha=OVERLAY_ALPHA,
        ncompress=min(MAX_NCOMPRESS, max(1, len(prop_amp) - 1)),
    )
    np.random.seed(seed_value)  # noqa: NPY002 -- see rationale above
    prop_sup.plot_2d(
        axes,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=sup,
        alpha=OVERLAY_ALPHA,
        ncompress=min(MAX_NCOMPRESS, max(1, len(prop_sup) - 1)),
    )
    if prop_only_cols:
        # Full-corner case (some columns propagating-only): slice axes
        # so NP only overlays the columns it has data for.
        np_axes = axes.loc[np_plot_params, np_plot_params]
        np.random.seed(seed_value)  # noqa: NPY002 -- see rationale above
        np_amp.plot_2d(
            np_axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=np_amp_,
            alpha=OVERLAY_ALPHA,
            ncompress=min(MAX_NCOMPRESS, max(1, len(np_amp) - 1)),
        )
        np.random.seed(seed_value)  # noqa: NPY002 -- see rationale above
        np_sup.plot_2d(
            np_axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=np_sup_,
            alpha=OVERLAY_ALPHA,
            ncompress=min(MAX_NCOMPRESS, max(1, len(np_sup) - 1)),
        )
    else:
        # Restricted case (all plot_params in NP list): plot NP directly
        # on the same `axes` object, matching the pre-round-15 behavior.
        np.random.seed(seed_value)  # noqa: NPY002 -- see rationale above
        np_amp.plot_2d(
            axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=np_amp_,
            alpha=OVERLAY_ALPHA,
            ncompress=min(MAX_NCOMPRESS, max(1, len(np_amp) - 1)),
        )
        np.random.seed(seed_value)  # noqa: NPY002 -- see rationale above
        np_sup.plot_2d(
            axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=np_sup_,
            alpha=OVERLAY_ALPHA,
            ncompress=min(MAX_NCOMPRESS, max(1, len(np_sup) - 1)),
        )

    fig = axes.iloc[0, 0].figure

    for ax in axes.values.flatten():
        if ax is None:
            continue
        ax.tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
        if strip_ylabels:
            ax.set_ylabel("")
        if strip_xlabels:
            ax.set_xlabel("")
        if ylabel_rotation is not None and ax.get_ylabel():
            ax.yaxis.label.set_rotation(ylabel_rotation)
            ax.yaxis.label.set_horizontalalignment("right")
            ax.yaxis.label.set_verticalalignment("center")
            ax.yaxis.labelpad = 10

    ll = legend_labels or {}
    lbl_prop_amp = ll.get("prop_amp", "amplification (propagating)")
    lbl_prop_sup = ll.get("prop_sup", "suppression (propagating)")
    lbl_np_amp = ll.get("np_amp", r"amplification ($\xi=0$ control)")
    lbl_np_sup = ll.get("np_sup", r"suppression ($\xi=0$ control)")
    legend_handles = [
        mpatches.Patch(color=amp, alpha=OVERLAY_ALPHA, label=lbl_prop_amp),
        mpatches.Patch(color=sup, alpha=OVERLAY_ALPHA, label=lbl_prop_sup),
        mpatches.Patch(color=np_amp_, alpha=OVERLAY_ALPHA, label=lbl_np_amp),
        mpatches.Patch(color=np_sup_, alpha=OVERLAY_ALPHA, label=lbl_np_sup),
    ]
    ax_anchor = axes.iloc[0, 0]
    n = len(plot_params)
    if legend_anchor is not None:
        loc, anchor = "upper right", legend_anchor
    elif n <= 3:
        loc, anchor = "upper left", (1.05, 1.0)
    else:
        loc, anchor = "upper right", (n - 0.5, 1.0)
    ax_anchor.legend(
        handles=legend_handles,
        loc=loc,
        bbox_to_anchor=anchor,
        bbox_transform=ax_anchor.transAxes,
        frameon=True,
        facecolor="none",
        edgecolor="#aaaaaa",
        framealpha=1.0,
    )

    fig.set_size_inches(fig_width, fig_width * 0.95)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if transparent:
        fig.patch.set_alpha(0.0)
    fig.savefig(out_path, format="pdf", bbox_inches="tight", transparent=transparent)
    plt.close(fig)
    log.info("[overlay_pair] wrote %s", out_path)


def main() -> None:
    figs_dir = REPO_ROOT / "manuscript" / "figures"

    # Pair A: chi-closure (18D propagating) vs NP-chi-closure (17D control)
    chi_prop_params = [
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
    chi_np_params = [
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
    chi_full_plot = (
        chi_prop_params  # 18 axes; xi column shows propagating-only contours
    )
    chi_prop_amp = REPO_ROOT / "hpc_results/29682868/t7_amp_v2"
    chi_prop_sup = REPO_ROOT / "hpc_results/29682868/t7_sup_v2"
    chi_np_amp = REPO_ROOT / "hpc_results/29705560/np_ceven_amp_v1"
    chi_np_sup = REPO_ROOT / "hpc_results/29705560/np_ceven_sup_v1"

    # Full overlay (17 axes)
    render_overlay_pair(
        out_path=figs_dir / "overlay_chi_closure_pair_full.pdf",
        plot_params=chi_full_plot,
        prop_amp_dir=chi_prop_amp,
        prop_sup_dir=chi_prop_sup,
        prop_param_names=chi_prop_params,
        np_amp_dir=chi_np_amp,
        np_sup_dir=chi_np_sup,
        np_param_names=chi_np_params,
        fig_width=FIG_WIDTH,
    )

    # Restricted top-K overlay (K=6) — rank by max across propagating amp
    # and NP amp parameter_importance.json, so the chosen subset includes
    # the most-structured couplings from both chains rather than only the
    # propagating chain's top.
    chi_top = _top_k_union(
        [
            chi_prop_amp / "parameter_importance.json",
            chi_np_amp / "parameter_importance.json",
        ],
        6,
        chi_np_params,
    )
    log.info("[chi_pair] top-6 carrier params: %s", chi_top)
    render_overlay_pair(
        out_path=figs_dir / "overlay_chi_closure_pair_restricted.pdf",
        plot_params=chi_top,
        prop_amp_dir=chi_prop_amp,
        prop_sup_dir=chi_prop_sup,
        prop_param_names=chi_prop_params,
        np_amp_dir=chi_np_amp,
        np_sup_dir=chi_np_sup,
        np_param_names=chi_np_params,
        fig_width=FIG_WIDTH,
    )

    # Pair B: YM-PGT non-minimal union (9D propagating) vs NP control (8D)
    union_prop_params = [
        "beta1",
        "beta2",
        "beta3",
        "xi",
        "delta1",
        "chi",
        "zeta1",
        "zeta2",
        "zeta3",
    ]
    union_np_params = [
        "beta1",
        "beta2",
        "beta3",
        "delta1",
        "chi",
        "zeta1",
        "zeta2",
        "zeta3",
    ]
    union_full_plot = (
        union_prop_params  # 9 axes; xi column shows propagating-only contours
    )
    union_prop_amp = REPO_ROOT / "hpc_results/29468763/d23_full_amp_v3"
    union_prop_sup = REPO_ROOT / "hpc_results/29471255/d23_full_sup_v3"
    union_np_amp = REPO_ROOT / "hpc_results/29700462/np_amp_v1"
    union_np_sup = REPO_ROOT / "hpc_results/29700462/np_sup_v1"

    render_overlay_pair(
        out_path=figs_dir / "overlay_ym_union_pair_full.pdf",
        plot_params=union_full_plot,
        prop_amp_dir=union_prop_amp,
        prop_sup_dir=union_prop_sup,
        prop_param_names=union_prop_params,
        np_amp_dir=union_np_amp,
        np_sup_dir=union_np_sup,
        np_param_names=union_np_params,
        fig_width=FIG_WIDTH,
    )

    union_top = _top_k_union(
        [
            union_prop_amp / "parameter_importance.json",
            union_np_amp / "parameter_importance.json",
        ],
        6,
        union_np_params,
    )
    log.info("[union_pair] top-6 carrier params: %s", union_top)
    render_overlay_pair(
        out_path=figs_dir / "overlay_ym_union_pair_restricted.pdf",
        plot_params=union_top,
        prop_amp_dir=union_prop_amp,
        prop_sup_dir=union_prop_sup,
        prop_param_names=union_prop_params,
        np_amp_dir=union_np_amp,
        np_sup_dir=union_np_sup,
        np_param_names=union_np_params,
        fig_width=COLUMN_WIDTH,
    )


if __name__ == "__main__":
    main()
