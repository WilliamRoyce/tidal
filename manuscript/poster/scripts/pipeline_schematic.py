r"""Compact TIDAL pipeline schematic for the poster (Section 2).

Four steps in a 2x2 serpentine (fits the narrow left column better than a 4-wide
row): Lagrangian -> linearised EOM (xAct) -> simulation (spectral solver) ->
nested sampling. Drawn in the modern Cambridge palette to match the LaTeX theme.
A new, poster-specific figure (does NOT reuse the report's figA1_architecture).

  [Lagrangian]      -->  [Linearised EOM]
                                |
                                v
  [Nested sampling] <--  [Simulation]

Output: manuscript/poster/figures/pipeline_schematic.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _poster_palette import CAM

OUT = Path(__file__).resolve().parents[1] / "figures" / "pipeline_schematic.pdf"

# One bold header per box (no sub-caption); long headers wrap across lines so
# they fit comfortably inside the box. Keyed by 2x2 grid position.
TITLES = {
    "TL": r"Lagrangian  $\mathscr{L}$",
    "TR": "Linearised\nEOM",
    "BR": "Simulation",
    "BL": "Nested\nsampling",
}

# Box geometry (data units). Tight around the (fixed-pt) labels: the figure uses a
# constant data-unit -> inch scale (figsize = (xlim, ylim) below), so shrinking these
# genuinely hugs the text rather than rescaling everything (a fixed-width figsize
# would keep the box:text ratio constant).
BOX_W, BOX_H = 2.6, 0.85
GX, GY = 0.9, 0.85  # horizontal / vertical gaps between boxes
MARGIN = 0.25
# FancyBboxPatch padding -- the visible rounded box extends this far beyond its
# rect on every side, so arrows must start/end OUTSIDE it to avoid overlapping,
# and the axis limits must clear it so no box edge is clipped.
BOX_PAD = 0.08
ARROW_CLEAR = 0.14  # gap between an arrow tip/tail and the box edge


def main() -> None:
    col0 = MARGIN
    col1 = MARGIN + BOX_W + GX
    bot = MARGIN
    top = MARGIN + BOX_H + GY
    xlim = col1 + BOX_W + MARGIN
    ylim = top + BOX_H + MARGIN

    # Box lower-left corners by grid position.
    pos = {
        "TL": (col0, top),
        "TR": (col1, top),
        "BR": (col1, bot),
        "BL": (col0, bot),
    }

    # Constant data-unit -> inch scale so box size is set relative to the fixed-pt
    # text (see the BOX_W/BOX_H note): 1 data unit = 1 inch.
    fig, ax = plt.subplots(figsize=(xlim, ylim))
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, ylim)
    ax.axis("off")
    # Light-blue background so the figure blends into its light-blue poster box
    # (this script does not call apply_poster_style(), so set it explicitly).
    fig.patch.set_facecolor(CAM["light_blue"])
    ax.set_facecolor(CAM["light_blue"])

    for key, (x, y) in pos.items():
        box = FancyBboxPatch(
            (x, y),
            BOX_W,
            BOX_H,
            boxstyle=f"round,pad={BOX_PAD},rounding_size=0.18",
            facecolor=CAM["light_blue"],
            edgecolor=CAM["dark_blue"],
            lw=2.5,
        )
        ax.add_patch(box)
        ax.text(
            x + BOX_W / 2,
            y + BOX_H / 2,
            TITLES[key],
            ha="center",
            va="center",
            color=CAM["dark_blue"],
            fontsize=20,
            fontweight="bold",
            linespacing=1.05,
        )

    def _arrow(p0: tuple[float, float], p1: tuple[float, float]) -> None:
        ax.add_patch(
            FancyArrowPatch(
                p0,
                p1,
                arrowstyle="-|>",
                mutation_scale=28,
                color=CAM["warm_blue"],
                lw=3.0,
            )
        )

    cy_top = top + BOX_H / 2
    cy_bot = bot + BOX_H / 2
    cx_right = col1 + BOX_W / 2
    # A: TL -> TR (rightwards along the top row)
    _arrow(
        (col0 + BOX_W + BOX_PAD + ARROW_CLEAR, cy_top),
        (col1 - BOX_PAD - ARROW_CLEAR, cy_top),
    )
    # B: TR -> BR (downwards along the right column)
    _arrow(
        (cx_right, top - BOX_PAD - ARROW_CLEAR),
        (cx_right, bot + BOX_H + BOX_PAD + ARROW_CLEAR),
    )
    # C: BR -> BL (leftwards along the bottom row)
    _arrow(
        (col1 - BOX_PAD - ARROW_CLEAR, cy_bot),
        (col0 + BOX_W + BOX_PAD + ARROW_CLEAR, cy_bot),
    )

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", facecolor=CAM["light_blue"])
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
