r"""Compact TIDAL pipeline schematic for the poster (Section 2).

Three steps: Lagrangian -> linearised EOM (xAct) -> Bayesian coupling-space
scan. Drawn in the modern Cambridge palette to match the LaTeX theme. A new,
poster-specific figure (does NOT reuse the report's figA1_architecture).

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

# One bold header per step (no sub-caption); long headers wrap across lines so
# they fit comfortably inside the box.
TITLES = [
    r"Lagrangian  $\mathscr{L}$",
    "Linearised\nEOM",
    "Nested\nsampling",
]

# Box geometry (data units).
BOX_W, BOX_H, GAP, X0 = 3.6, 1.7, 1.0, 0.25
# FancyBboxPatch padding -- the visible rounded box extends this far beyond its
# rect on every side, so arrows must start/end OUTSIDE it to avoid overlapping,
# and the axis x-limit must clear it so no box edge is clipped.
BOX_PAD = 0.08
ARROW_CLEAR = 0.14  # gap between an arrow tip/tail and the box edge


def main() -> None:
    n = len(TITLES)
    # Total content width incl. a right margin, so the last box never clips.
    xlim = X0 + n * BOX_W + (n - 1) * GAP + X0
    fig, ax = plt.subplots(figsize=(12, 2.6))
    ax.set_xlim(0, xlim)
    ax.set_ylim(0, 2.6)
    ax.axis("off")
    # Light-blue background so the figure blends into its light-blue poster box
    # (this script does not call apply_poster_style(), so set it explicitly).
    fig.patch.set_facecolor(CAM["light_blue"])
    ax.set_facecolor(CAM["light_blue"])

    y0 = 0.45
    cy = y0 + BOX_H / 2  # vertical centre, shared by titles and arrows
    x = X0
    lefts = []
    for title in TITLES:
        box = FancyBboxPatch(
            (x, y0),
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
            cy,
            title,
            ha="center",
            va="center",
            color=CAM["dark_blue"],
            fontsize=26,
            fontweight="bold",
            linespacing=1.05,
        )
        lefts.append(x)
        x += BOX_W + GAP

    for i in range(n - 1):
        # Start just outside the right edge of box i; end just outside the left
        # edge of box i+1 (account for the rounded-box padding + clearance) so the
        # arrow sits cleanly in the gap and does not overlap either box.
        x0 = lefts[i] + BOX_W + BOX_PAD + ARROW_CLEAR
        x1 = lefts[i + 1] - BOX_PAD - ARROW_CLEAR
        ax.add_patch(
            FancyArrowPatch(
                (x0, cy),
                (x1, cy),
                arrowstyle="-|>",
                mutation_scale=28,
                color=CAM["warm_blue"],
                lw=3.0,
            )
        )

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight", facecolor=CAM["light_blue"])
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
