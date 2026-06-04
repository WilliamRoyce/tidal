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

STEPS = [
    (r"Lagrangian  $\mathscr{L}$", "any gravity theory\nwith torsion"),
    ("Linearised EOM", "symbolic derivation\n(xAct)"),
    ("Bayesian scan", "where is the\nconversion amplified?"),
]


def main() -> None:
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 3.2)
    ax.axis("off")

    box_w, box_h, gap = 3.2, 2.0, 1.0
    x = 0.3
    centres = []
    for title, sub in STEPS:
        box = FancyBboxPatch(
            (x, 0.6),
            box_w,
            box_h,
            boxstyle="round,pad=0.08,rounding_size=0.18",
            facecolor=CAM["light_blue"],
            edgecolor=CAM["dark_blue"],
            lw=2.5,
        )
        ax.add_patch(box)
        cx = x + box_w / 2
        ax.text(
            cx,
            1.95,
            title,
            ha="center",
            va="center",
            color=CAM["dark_blue"],
            fontsize=24,
            fontweight="bold",
        )
        ax.text(
            cx, 1.15, sub, ha="center", va="center", color=CAM["slate4"], fontsize=17
        )
        centres.append((x, cx))
        x += box_w + gap

    for i in range(len(STEPS) - 1):
        x0 = centres[i][0] + box_w
        x1 = centres[i + 1][0]
        ax.add_patch(
            FancyArrowPatch(
                (x0, 1.6),
                (x1, 1.6),
                arrowstyle="-|>",
                mutation_scale=26,
                color=CAM["warm_blue"],
                lw=3.0,
            )
        )

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
