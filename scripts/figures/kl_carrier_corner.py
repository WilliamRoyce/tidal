r"""Carrier sub-corner: anesthetic overlay corner over the top-$K$
parameters of each theory ranked by combined $D_\mathrm{KL}$.

Reads `parameter_importance.json` next to each amp chain (emitted by
`scripts/analysis/recompute_parameter_kl.py`), picks the top-$K$
parameters by combined ranking
(``max(amp_marginal, sup_marginal, cross_amp_sup)``), and renders the
existing-style overlay corner over just those columns.  The output PDF
filename is `manuscript/figures/kl_carrier_corner_<theory>.pdf` for each
processed theory.

Per-theory $K$ values:
  - $K = $ all params when the theory has $\leq 6$ parameters (the full
    corner is the carrier corner).
  - $K = 6$ otherwise.

Driven by the same `d_kl_chains.toml` config as the recompute walker
so it picks up new chains (e.g. localised geometry) without code edits.

Run from the repo root:
    uv run python scripts/figures/kl_carrier_corner.py
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from pathlib import Path

import matplotlib as mpl

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "figures"))
sys.path.insert(0, str(REPO_ROOT))

try:
    import tomllib  # type: ignore[import]
except ImportError:
    import tomli as tomllib  # type: ignore[import]

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

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("kl_carrier_corner")


# Same PARAM_LABEL family as scripts/analysis/extract_d_kl.py — kept in
# sync so the figure axes match the table parameter columns.
PARAM_LABEL: dict[str, str] = {
    "beta1": r"\beta_1",
    "beta2": r"\beta_2",
    "beta3": r"\beta_3",
    "alpha1": r"\alpha_1",
    "alpha2": r"\alpha_2",
    "alpha3": r"\alpha_3",
    "xi": r"\xi",
    "chi": r"\chi",
    "delta1": r"\delta_1",
    "zeta1": r"\zeta_1",
    "zeta2": r"\zeta_2",
    "zeta3": r"\zeta_3",
    "d14": r"d_{14}",
    "d15": r"d_{15}",
    "d17": r"d_{17}",
    "d19": r"d_{19}",
    "d20": r"d_{20}",
    "d21": r"d_{21}",
    "mA2": r"m_A^{2}",
    "deltam": r"\delta_m",
    "rho": r"\rho",
    "sigma": r"\sigma",
}


def _label(name: str) -> str:
    if name in PARAM_LABEL:
        return f"${PARAM_LABEL[name]}$"
    for prefix, sym in (
        ("tildechi", r"\tilde\chi"),
        ("tildezeta", r"\tilde\zeta"),
        ("chi", r"\chi"),
        ("xi", r"\xi"),
        ("zeta", r"\zeta"),
    ):
        if name.startswith(prefix) and name[len(prefix) :].isdigit():
            idx = int(name[len(prefix) :])
            return f"${sym}_{{{idx}}}$"
    return rf"$\mathtt{{{name}}}$"


def _pick_top_k(importance: dict, k: int) -> list[str]:
    """Rank parameters by combined max(amp marginal, sup marginal, cross)."""
    amp = importance.get("amp", {}).get("marginal_d_kl", {}) or {}
    sup = importance.get("sup", {}).get("marginal_d_kl", {}) or {}
    cross = importance.get("cross_amp_sup_kl", {}) or {}
    names = list(amp.keys() or sup.keys() or cross.keys())

    def score(name: str) -> float:
        return max(
            amp.get(name, 0.0) if math.isfinite(amp.get(name, 0.0) or 0.0) else 0.0,
            sup.get(name, 0.0) if math.isfinite(sup.get(name, 0.0) or 0.0) else 0.0,
            cross.get(name, 0.0) if math.isfinite(cross.get(name, 0.0) or 0.0) else 0.0,
        )

    return sorted(names, key=lambda n: -score(n))[:k]


def _process(entry: dict) -> None:
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt

    theory = entry["theory"]
    amp_path = REPO_ROOT / entry["amp_path"]
    sup_path = REPO_ROOT / entry["sup_path"]
    json_path = amp_path / "parameter_importance.json"
    if not json_path.is_file():
        log.warning("[%s] missing parameter_importance.json — skipping", theory)
        return
    with json_path.open() as f:
        imp = json.load(f)

    n_params = imp.get("n_params", 0)
    k = min(n_params, 6)
    plot_params = _pick_top_k(imp, k)
    if not plot_params:
        log.warning("[%s] no rankable parameters — skipping", theory)
        return

    # Load with the FULL ordered parameter list from the TOML (matches the
    # chain column count); plot only the carrier subset.
    full_params = list(entry["param_names"])
    full_labels = {name: _label(name) for name in full_params}

    # wide=true in the TOML means the figure is included at \textwidth (figure*)
    # in the manuscript; all others are included at \columnwidth.  The figure
    # must be generated at the matching width so that absolute-pt font sizes
    # land at the intended pt values in the final PDF rather than being
    # halved by a 7in→3.375in scale-down.
    fig_w = FIG_WIDTH if entry.get("wide") else COLUMN_WIDTH

    apply_style()
    amp_ns = load_chains(amp_path, params=full_params, param_labels=full_labels)
    sup_ns = load_chains(sup_path, params=full_params, param_labels=full_labels)

    log.info("[%s] rendering top-%d carrier sub-corner: %s", theory, k, plot_params)
    axes = amp_ns.plot_2d(
        plot_params,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=AMP_COLOR,
        alpha=OVERLAY_ALPHA,
    )
    sup_ns.plot_2d(
        axes,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=SUP_COLOR,
        alpha=OVERLAY_ALPHA,
    )
    fig = axes.iloc[0, 0].figure

    # For wide (high-D) carrier corners apply the same tick-density and
    # tick-label-size fixes as overlay_corner() does for full-chain plots.
    # Only applies when wide=true — column-width figures carrying ≤6 params
    # in the sub-corner get tick labels removed instead (see block below).
    if entry.get("wide") and n_params >= HIGH_DIM_THRESHOLD:
        from matplotlib.ticker import MaxNLocator

        tick_label_pt = 5 if n_params >= 20 else 6
        for ax in axes.values.flatten():
            if ax is None:
                continue
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.tick_params(axis="both", labelsize=tick_label_pt)
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_horizontalalignment("right")

    # Single-column sub-corners showing >5 params: remove numeric tick labels.
    if k > 5 and not entry.get("wide"):
        for ax in axes.values.flatten():
            if ax is None:
                continue
            ax.tick_params(labelbottom=False, labelleft=False)

    # Add legend anchored to the top-right of the corner grid.

    legend_handles = [
        mpatches.Patch(color=AMP_COLOR, alpha=OVERLAY_ALPHA, label="amplification"),
        mpatches.Patch(color=SUP_COLOR, alpha=OVERLAY_ALPHA, label="suppression"),
    ]
    ax_anchor = axes.iloc[0, 0]
    n_plot = len(plot_params)
    loc = "upper left" if n_plot <= 3 else "upper right"
    anchor = (1.05, 1.0) if n_plot <= 3 else (n_plot - 0.5, 1.0)
    ax_anchor.legend(
        handles=legend_handles,
        loc=loc,
        bbox_to_anchor=anchor,
        bbox_transform=ax_anchor.transAxes,
        frameon=True,
        facecolor="none",
        edgecolor="#aaaaaa",
        framealpha=1.0,
        fontsize=mpl.rcParams["legend.fontsize"],
    )

    out_pdf = REPO_ROOT / "manuscript" / "figures" / f"kl_carrier_corner_{theory}.pdf"
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.set_size_inches(fig_w, fig_w * 0.95)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)
    log.info("[%s] wrote %s", theory, out_pdf.relative_to(REPO_ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--config",
        type=Path,
        default=REPO_ROOT / "scripts" / "analysis" / "d_kl_chains.toml",
    )
    parser.add_argument(
        "--theory",
        type=str,
        default=None,
        help="If set, only render the named theory.",
    )
    args = parser.parse_args()

    with args.config.open("rb") as f:
        cfg = tomllib.load(f)

    for entry in cfg.get("chain", []):
        if args.theory and entry["theory"] != args.theory:
            continue
        try:
            _process(entry)
        except Exception:
            log.exception("[%s] uncaught render error", entry.get("theory", "?"))


if __name__ == "__main__":
    main()
