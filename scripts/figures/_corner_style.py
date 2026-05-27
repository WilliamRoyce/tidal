r"""Shared style and plotting helpers for the §4 results-section corner plots.

All per-theory scripts in ``scripts/figures/corner_*.py`` import from this
module. Changes here propagate to every corner plot, ensuring style
consistency across the §4 figure set.

The amplification-rewarding (``+log A``) posterior is drawn in IBM magenta and the
suppression-rewarding (``-log A``) posterior in IBM yellow, per the App J
convention of \\cref{NestedScore} (see
``manuscript/sections/appendices/inference_architecture.tex``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from anesthetic import read_chains
from matplotlib.ticker import MaxNLocator

if TYPE_CHECKING:
    from pathlib import Path

# Manuscript style — absolute pt sizes set the physical text size on the
# rendered PDF and DO NOT scale with figure inches.
# All sizes are chosen for figures rendered at their native width (COLUMN_WIDTH
# or FIG_WIDTH) so that they appear at these pt sizes in the final PDF.
# High-D corners apply a post-hoc tick-label shrink in overlay_corner().
RCPARAMS = {
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{stix}",
    "font.family": "serif",
    "font.size": 10,  # matches caption body
    "axes.labelsize": 12,  # parameter labels (β, χ, ξ) comfortably above body text
    "legend.fontsize": 10,  # matches body text — no longer tiny relative to figure
    "xtick.labelsize": 8,  # adequate for small corners; high-D post-processing shrinks further
    "ytick.labelsize": 8,  # matched to xtick
    "lines.linewidth": 1.0,
    "axes.linewidth": 0.5,
}

# Default figure width (PRD two-column, \textwidth). Per-script callers may
# override via the fig_width kwarg of overlay_corner(): small-N corners
# (≤ 9D) render at single-column width (3.375 in) so the absolute-pt labels
# read proportionally larger; high-N corners (T6 full, T7, T9) use the
# default two-column width.
FIG_WIDTH = 7.0
COLUMN_WIDTH = 3.375
HIGH_DIM_THRESHOLD = 12  # n_params at which the high-D tick fixes apply

# Posterior colours sourced from the central IBM colorblind palette.
# Magenta = amplification-rewarding likelihood (+log A);
# yellow  = suppression-rewarding likelihood (-log A); per App J
# \cref{NestedScore}.  See `_palette.py` for the full palette.
from _palette import AMP_COLOR, SUP_COLOR

# Filled contour levels: 2sigma then 1sigma (anesthetic ≥2.0.0-beta.10
# requires descending order, outermost-to-innermost).
CONTOUR_LEVELS = [0.954, 0.683]

# Alpha so the union of overlaid amp+sup posteriors reads cleanly.
OVERLAY_ALPHA = 0.5


def apply_style() -> None:
    """Apply the shared matplotlib rcParams."""
    mpl.rcParams.update(RCPARAMS)


def _resolve_chains_root(chains_dir: Path) -> str:
    """Resolve the anesthetic file-root for a chain directory.

    Accepts either the chain parent directory (with `_chains/tidal*` inside)
    or the `_chains` directory directly. Returns the path-string anesthetic
    expects as `read_chains` argument.
    """
    candidates = [
        chains_dir / "_chains" / "tidal",
        chains_dir / "tidal",
        chains_dir,
    ]
    for cand in candidates:
        if any(cand.parent.glob(f"{cand.name}*.txt")) or any(
            cand.parent.glob(f"{cand.name}*")
        ):
            return str(cand)
    msg = f"No chain files found in {chains_dir} or its _chains subdir"
    raise FileNotFoundError(msg)


def load_chains(
    chains_dir: Path,
    params: list[str] | None = None,
    param_labels: dict[str, str] | None = None,
):
    """Load a PolyChord chain (full or partial) via anesthetic.

    ``chains_dir`` may be the chain parent directory or the ``_chains``
    subdirectory. PolyChord output here lacks a ``.paramnames`` file, so
    ``params`` and ``param_labels`` are passed through to ``read_chains`` to
    name the columns and supply LaTeX labels for the corner-plot axes.
    Anesthetic handles truncated dead-chain files transparently, returning a
    ``NestedSamples`` whose ``logZ`` and ``plot_2d`` methods both work on
    partial output.
    """
    root = _resolve_chains_root(chains_dir)
    read_kwargs: dict = {}
    if params is not None:
        read_kwargs["columns"] = params
    if param_labels is not None:
        read_kwargs["labels"] = param_labels
    return read_chains(root, **read_kwargs)


def overlay_corner(
    *,
    amp_chains_dir: Path | None,
    sup_chains_dir: Path | None,
    params: list[str],
    param_labels: dict[str, str] | None = None,
    out_path: Path,
    title: str | None = None,
    height_ratio: float = 0.95,
    fig_width: float = FIG_WIDTH,
    prior_samples=None,
    legend_kw: dict | None = None,
) -> None:
    """Render an overlaid amp+sup corner plot to PDF.

    Either ``amp_chains_dir`` or ``sup_chains_dir`` may be ``None`` for
    single-posterior plots (e.g. partial campaigns with only one direction
    completed); the function then renders a single-colour plot.

    ``fig_width`` controls the rendered figure width in inches; small-N
    corner plots (≤ 9D) should pass ``fig_width=COLUMN_WIDTH`` so the
    absolute-pt axis labels read proportionally larger on the rendered
    figure. Default is two-column ``FIG_WIDTH``.

    ``prior_samples`` (optional) is an anesthetic ``Samples``/``NestedSamples``
    object holding samples from the prior; when provided, the prior is drawn
    as a low-alpha bottom layer (Legner ``fig:TorCprior`` template). Default
    is ``None`` — for our arctan-uniform / log-uniform priors the visual
    contribution is marginal, so the prior overlay is opt-in.

    ``legend_kw`` (optional) overrides any keyword arguments forwarded to
    ``ax.legend()``.  Use to reposition the legend for unusual grid shapes
    (e.g. 2-parameter plots where the default upper-right anchor overlaps the
    plotted panels).
    """
    from _palette import PRIOR_ALPHA, PRIOR_COLOR

    apply_style()
    sources = []
    if prior_samples is not None:
        sources.append(("prior", prior_samples, PRIOR_COLOR, PRIOR_ALPHA))
    if amp_chains_dir is not None:
        sources.append(
            (
                "amp",
                load_chains(amp_chains_dir, params=params, param_labels=param_labels),
                AMP_COLOR,
                OVERLAY_ALPHA,
            )
        )
    if sup_chains_dir is not None:
        sources.append(
            (
                "sup",
                load_chains(sup_chains_dir, params=params, param_labels=param_labels),
                SUP_COLOR,
                OVERLAY_ALPHA,
            )
        )
    if not sources:
        msg = "At least one of amp_chains_dir / sup_chains_dir is required"
        raise ValueError(msg)

    # Initialise axes from the first available source, then overlay any
    # additional sources on the same axes. Anesthetic's plot_2d returns
    # an AxesDataFrame (pandas-like), not a (fig, axes) tuple.  Each source
    # in `sources` is a (tag, samples, color, alpha) 4-tuple so the prior
    # overlay can use a faint alpha while amp/sup keep their primary alpha.
    _, first_ns, first_color, first_alpha = sources[0]
    # Corner plots show ONLY the lower triangle + diagonal marginals; the
    # upper triangle would just mirror the same 2D distributions and adds
    # no new information. The 'kde' shortcut puts 1D KDE on the diagonal
    # and 2D KDE in the lower triangle (no upper triangle).
    axes = first_ns.plot_2d(
        params,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=first_color,
        alpha=first_alpha,
    )
    for _, ns, color, alpha in sources[1:]:
        ns.plot_2d(
            axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=color,
            alpha=alpha,
        )

    fig = axes.iloc[0, 0].figure
    # Anesthetic propagates param_labels via the labels= kwarg to read_chains
    # (set in load_chains above) — no per-axis label-loop required here.

    # High-D tick fixes: for corners with ≥ HIGH_DIM_THRESHOLD parameters,
    # apply 45° diagonal rotation on x-tick labels and limit major-tick
    # density to 3 ticks per axis. Tick labels (the numeric values) are
    # shrunk per-panel because the panel width on a \textwidth figure
    # scales as 1/n_params; without the shrinkage, tick numerals overlap
    # by ~20 panels. The PARAMETER labels ($\beta_1$, $\chi_3$, etc.) are
    # NOT shrunk — they read at the RCPARAMS axes.labelsize default.
    #   12–19 params (T7 18D, NP.T7 17D): tick labels 6 pt
    #   ≥ 20 params  (T9 32D, T6 20D):    tick labels 5 pt
    if len(params) >= HIGH_DIM_THRESHOLD:
        tick_label_pt = 5 if len(params) >= 20 else 6
        for ax in axes.values.flatten():
            if ax is None:
                continue
            ax.xaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=3, prune="both"))
            ax.tick_params(axis="both", labelsize=tick_label_pt)
            for label in ax.get_xticklabels():
                label.set_rotation(45)
                label.set_horizontalalignment("right")

    # Legend hosted INSIDE the empty upper-right region of the corner-plot
    # grid, anchored to the topmost diagonal axis (axes.iloc[0,0]) in its
    # transAxes coordinates. The horizontal offset (n_params − 0.5) is in
    # units of axes.iloc[0,0]'s width (panel-widths), placing the legend's
    # upper-right corner at roughly the right edge of the corner-plot's
    # column span, at the top of the topmost row. The legend then extends
    # leftward and downward into the empty upper-triangle region. This is
    # robust to bbox_inches='tight' cropping because the anchor scales with
    # the axes grid itself, not the figure margins.
    legend_handles = []
    if any(tag == "prior" for tag, *_ in sources):
        legend_handles.append(
            mpatches.Patch(color=PRIOR_COLOR, alpha=PRIOR_ALPHA, label="prior")
        )
    if any(tag == "amp" for tag, *_ in sources):
        legend_handles.append(
            mpatches.Patch(color=AMP_COLOR, alpha=OVERLAY_ALPHA, label="amplification")
        )
    if any(tag == "sup" for tag, *_ in sources):
        legend_handles.append(
            mpatches.Patch(color=SUP_COLOR, alpha=OVERLAY_ALPHA, label="suppression")
        )
    if legend_handles:
        ax_anchor = axes.iloc[0, 0]
        n_params = len(params)
        # Default: anchor the legend inside the empty upper-right region.
        # For very small grids (≤ 3 params) the upper-right cell is only one
        # panel wide; the legend can bleed into the 2D joint panel below it.
        # Place it above the figure instead (y > 1 in transAxes → above the
        # top edge; bbox_inches='tight' then includes it cleanly).
        if n_params <= 3:
            loc, anchor = "lower center", (0.5, 1.08)
        else:
            loc, anchor = "upper right", (n_params - 0.5, 1.0)
        leg_kw: dict = {
            "handles": legend_handles,
            "loc": loc,
            "bbox_to_anchor": anchor,
            "bbox_transform": ax_anchor.transAxes,
            "frameon": False,
            "fontsize": mpl.rcParams["legend.fontsize"],
        }
        if legend_kw:
            leg_kw.update(legend_kw)
        ax_anchor.legend(**leg_kw)

    if title is not None:
        fig.suptitle(title, y=1.02)
    fig.set_size_inches(fig_width, fig_width * height_ratio)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
