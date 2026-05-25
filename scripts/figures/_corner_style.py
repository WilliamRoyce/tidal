r"""Shared style and plotting helpers for the §4 results-section corner plots.

All per-theory scripts in ``scripts/figures/corner_*.py`` import from this
module. Changes here propagate to every corner plot, ensuring style
consistency across the §4 figure set.

The amplification-rewarding (``+log A``) posterior is drawn in red and the
suppression-rewarding (``-log A``) posterior in blue, per the App J
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
# rendered PDF and DO NOT scale with figure inches. Parameter (axis) labels
# are slightly larger than caption body so they read as prominent on the
# rendered figure, even when the figure is rendered at single-column width.
RCPARAMS = {
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}\usepackage{stix}",
    "font.family": "serif",
    "font.size": 10,  # was 9 — matches caption body
    "axes.labelsize": 11,  # was 9 — parameter labels (β, χ, ξ) prominent over caption body
    "legend.fontsize": 9,  # was 8 — matches caption body
    "xtick.labelsize": 6,  # reduced from 8 — diagonal labels on dense corners must not overlap
    "ytick.labelsize": 6,  # matched to xtick
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

# Posterior colours. Red = amplification-rewarding likelihood (+log A);
# blue = suppression-rewarding likelihood (-log A); per App J
# \cref{NestedScore}.
AMP_COLOR = "#d62728"
SUP_COLOR = "#1f77b4"

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
) -> None:
    """Render an overlaid amp+sup corner plot to PDF.

    Either ``amp_chains_dir`` or ``sup_chains_dir`` may be ``None`` for
    single-posterior plots (e.g. partial campaigns with only one direction
    completed); the function then renders a single-colour plot.

    ``fig_width`` controls the rendered figure width in inches; small-N
    corner plots (≤ 9D) should pass ``fig_width=COLUMN_WIDTH`` so the
    absolute-pt axis labels read proportionally larger on the rendered
    figure. Default is two-column ``FIG_WIDTH``.
    """
    apply_style()
    sources = []
    if amp_chains_dir is not None:
        sources.append(
            (
                "amp",
                load_chains(amp_chains_dir, params=params, param_labels=param_labels),
                AMP_COLOR,
            )
        )
    if sup_chains_dir is not None:
        sources.append(
            (
                "sup",
                load_chains(sup_chains_dir, params=params, param_labels=param_labels),
                SUP_COLOR,
            )
        )
    if not sources:
        msg = "At least one of amp_chains_dir / sup_chains_dir is required"
        raise ValueError(msg)

    # Initialise axes from the first available posterior, then overlay any
    # additional posteriors on the same axes. Anesthetic's plot_2d returns
    # an AxesDataFrame (pandas-like), not a (fig, axes) tuple.
    _, first_ns, first_color = sources[0]
    # Corner plots show ONLY the lower triangle + diagonal marginals; the
    # upper triangle would just mirror the same 2D distributions and adds
    # no new information. The 'kde' shortcut puts 1D KDE on the diagonal
    # and 2D KDE in the lower triangle (no upper triangle).
    axes = first_ns.plot_2d(
        params,
        kinds="kde",
        levels=CONTOUR_LEVELS,
        color=first_color,
        alpha=OVERLAY_ALPHA,
    )
    for _, ns, color in sources[1:]:
        ns.plot_2d(
            axes,
            kinds="kde",
            levels=CONTOUR_LEVELS,
            color=color,
            alpha=OVERLAY_ALPHA,
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
    # NOT shrunk — they read at the RCPARAMS \axes.labelsize\ default and
    # must stay at that size for legibility.
    #   12–19 params (T7 18D, NP.T7 17D): tick labels 5 pt
    #   ≥ 20 params  (T9 32D, T6 20D):    tick labels 4 pt
    if len(params) >= HIGH_DIM_THRESHOLD:
        tick_label_pt = 4 if len(params) >= 20 else 5
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
    if any(tag == "amp" for tag, _, _ in sources):
        legend_handles.append(
            mpatches.Patch(color=AMP_COLOR, alpha=OVERLAY_ALPHA, label="amplification")
        )
    if any(tag == "sup" for tag, _, _ in sources):
        legend_handles.append(
            mpatches.Patch(color=SUP_COLOR, alpha=OVERLAY_ALPHA, label="suppression")
        )
    if legend_handles:
        ax_anchor = axes.iloc[0, 0]
        n_params = len(params)
        ax_anchor.legend(
            handles=legend_handles,
            loc="upper right",
            bbox_to_anchor=(n_params - 0.5, 1.0),
            bbox_transform=ax_anchor.transAxes,
            frameon=False,
            fontsize=mpl.rcParams["legend.fontsize"],
        )

    if title is not None:
        fig.suptitle(title, y=1.02)
    fig.set_size_inches(fig_width, fig_width * height_ratio)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
