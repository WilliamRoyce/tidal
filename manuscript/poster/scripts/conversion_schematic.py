r"""Graphical abstract: graviton -> photon conversion in a magnetic field.

The 30-second-pitch hero image for the poster. A gravitational wave enters a
region of background magnetic field B0 and emerges (partly) as a photon. Drawn
in the modern Cambridge palette to match the LaTeX theme.

Output: manuscript/poster/figures/conversion_schematic.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _poster_palette import CAM, CHERRY, CREST, apply_poster_style

OUT = Path(__file__).resolve().parents[1] / "figures" / "conversion_schematic.pdf"


def _wave(ax, x0, x1, y, amp, n, color, lw=3.0) -> None:
    x = np.linspace(x0, x1, 400)
    ax.plot(
        x,
        y + amp * np.sin(2 * np.pi * n * (x - x0) / (x1 - x0)),
        color=color,
        lw=lw,
        solid_capstyle="round",
    )


def main() -> None:
    apply_poster_style()
    fig, ax = plt.subplots(figsize=(11, 4.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")

    # --- background magnetic-field region -----------------------------------
    region = Rectangle(
        (4.3, 0.4),
        3.4,
        4.2,
        facecolor=CAM["light_blue"],
        edgecolor=CAM["warm_blue"],
        lw=2.5,
        zorder=1,
    )
    ax.add_patch(region)
    for xb in np.linspace(4.7, 7.3, 5):
        ax.annotate(
            "",
            xy=(xb, 4.35),
            xytext=(xb, 0.65),
            arrowprops={"arrowstyle": "-|>", "color": CAM["warm_blue"], "lw": 2.0},
            zorder=2,
        )
    ax.text(
        6.0,
        4.95,
        r"$B_0$",
        color=CAM["dark_blue"],
        fontsize=30,
        ha="center",
        va="bottom",
        fontweight="bold",
    )

    # --- incoming gravitational wave (left) ---------------------------------
    _wave(ax, 0.4, 4.3, 2.5, 0.55, 5, CAM["dark_blue"])
    ax.text(
        2.3, 3.7, "gravitational wave", color=CAM["dark_blue"], fontsize=22, ha="center"
    )
    ax.text(0.5, 2.5, r"$h$", color=CAM["dark_blue"], fontsize=30, va="center")

    # --- outgoing photon (right) --------------------------------------------
    _wave(ax, 7.7, 11.6, 2.5, 0.55, 7, CHERRY)
    ax.text(9.7, 3.7, "photon", color=CHERRY, fontsize=22, ha="center")
    ax.text(11.55, 2.5, r"$\gamma$", color=CHERRY, fontsize=30, va="center", ha="right")

    # --- conversion arrow through the region --------------------------------
    ax.add_patch(
        FancyArrowPatch(
            (4.3, 2.5),
            (7.7, 2.5),
            arrowstyle="-|>",
            mutation_scale=28,
            color=CREST,
            lw=3.0,
            zorder=3,
        )
    )

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
