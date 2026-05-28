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
from matplotlib.ticker import MaxNLocator

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import (
    AMP_COLOR,
    COLUMN_WIDTH,
    CONTOUR_LEVELS,
    FIG_WIDTH,
    HIGH_DIM_THRESHOLD,
    OVERLAY_ALPHA,
    SUP_COLOR,
    apply_style,
    load_chains,
)
from _palette import IBM_PALETTE

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("overlay_corner_pair")

# Four-colour scheme: existing AMP/SUP for propagating; IBM purple/orange for NP control.
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


def _top_k(amp_imp_path: Path, k: int, allowed: list[str]) -> list[str]:
    """Pick top-K parameters by combined max(amp marginal, sup marginal, cross)."""
    with amp_imp_path.open() as fh:
        imp = json.load(fh)
    amp = imp.get("amp", {}).get("marginal_d_kl", {}) or {}
    sup = imp.get("sup", {}).get("marginal_d_kl", {}) or {}
    cross = imp.get("cross_amp_sup_kl", {}) or {}
    allowed_set = set(allowed)
    names = [n for n in amp if n in allowed_set] or [n for n in sup if n in allowed_set]

    def score(n: str) -> float:
        vals = [amp.get(n, 0.0), sup.get(n, 0.0), cross.get(n, 0.0)]
        vals = [v if v is not None and math.isfinite(v) else 0.0 for v in vals]
        return max(vals)

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
) -> None:
    apply_style()

    {n: _label(n) for n in plot_params}

    prop_amp = load_chains(
        prop_amp_dir,
        params=prop_param_names,
        param_labels={n: _label(n) for n in prop_param_names},
    )
    prop_sup = load_chains(
        prop_sup_dir,
        params=prop_param_names,
        param_labels={n: _label(n) for n in prop_param_names},
    )
    np_amp = load_chains(
        np_amp_dir,
        params=np_param_names,
        param_labels={n: _label(n) for n in np_param_names},
    )
    np_sup = load_chains(
        np_sup_dir,
        params=np_param_names,
        param_labels={n: _label(n) for n in np_param_names},
    )

    # All four are plotted on the same `plot_params` columns. The NP chains
    # do not contain xi so plot_params must be a subset of np_param_names.
    axes = prop_amp.plot_2d(
        plot_params,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=AMP_COLOR,
        alpha=OVERLAY_ALPHA,
    )
    for ns, color in [
        (prop_sup, SUP_COLOR),
        (np_amp, NP_AMP_COLOR),
        (np_sup, NP_SUP_COLOR),
    ]:
        ns.plot_2d(
            axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=color,
            alpha=OVERLAY_ALPHA,
        )

    fig = axes.iloc[0, 0].figure

    if len(plot_params) >= HIGH_DIM_THRESHOLD:
        tick_label_pt = 5 if len(plot_params) >= 20 else 6
        for ax in axes.values.flatten():
            if ax is None:
                continue
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.tick_params(axis="both", labelsize=tick_label_pt)
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_horizontalalignment("right")

    legend_handles = [
        mpatches.Patch(
            color=AMP_COLOR, alpha=OVERLAY_ALPHA, label="amplification (propagating)"
        ),
        mpatches.Patch(
            color=SUP_COLOR, alpha=OVERLAY_ALPHA, label="suppression (propagating)"
        ),
        mpatches.Patch(
            color=NP_AMP_COLOR,
            alpha=OVERLAY_ALPHA,
            label=r"amplification ($\xi=0$ control)",
        ),
        mpatches.Patch(
            color=NP_SUP_COLOR,
            alpha=OVERLAY_ALPHA,
            label=r"suppression ($\xi=0$ control)",
        ),
    ]
    ax_anchor = axes.iloc[0, 0]
    n = len(plot_params)
    if n <= 3:
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
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
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
        chi_np_params  # 17 axes; the 18th (xi) is excluded since NP has no xi
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

    # Restricted top-K overlay (K=6) — rank by propagating amp's parameter_importance
    chi_top = _top_k(chi_prop_amp / "parameter_importance.json", 6, chi_np_params)
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
    union_full_plot = union_np_params  # 8 axes
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

    union_top = _top_k(union_prop_amp / "parameter_importance.json", 6, union_np_params)
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
