r"""Central IBM colorblind-friendly palette for §4 figures.

All figure scripts in `scripts/figures/` should import their colours from this
module rather than hardcoding hex values, so the palette stays consistent
across the manuscript.

Palette source: https://www.color-hex.com/color-palette/1044488
(also cited by Legner–Handley–Barker 2507.09228 "SL would like to thank the
IBM colour blind friendly palette for the color scheme used in the figures").

Role assignment for our overlay-corner figures:
- AMP (amplification posterior, +log A direction)  → ORANGE
- SUP (suppression posterior, -log A direction)    → BLUE
- PRIOR (optional third overlay layer)             → MAGENTA at low alpha
- HIGHLIGHT (carrier-parameter accent in tables)   → YELLOW
- RESERVED (future fourth overlay, e.g. localised) → PURPLE

The amp / sup assignment replaces the previous matplotlib-default
`#d62728 / #1f77b4` red/blue pair with the colorblind-safe orange/blue
variant; same semantic role, accessibility improved.
"""

from __future__ import annotations

# Five-colour IBM colorblind palette (hex, RGB).
IBM_PALETTE = {
    "yellow": "#ffb000",
    "orange": "#fe6100",
    "magenta": "#dc267f",
    "purple": "#785ef0",
    "blue": "#648fff",
}

AMP_COLOR = IBM_PALETTE["orange"]
SUP_COLOR = IBM_PALETTE["blue"]
PRIOR_COLOR = IBM_PALETTE["magenta"]
HIGHLIGHT_COLOR = IBM_PALETTE["yellow"]
RESERVED_COLOR = IBM_PALETTE["purple"]

# Default alphas tuned for the overlay-corner visual: amp / sup carry the
# primary signal; prior (when drawn) is a faint reference cloud.
AMP_ALPHA = 0.5
SUP_ALPHA = 0.5
PRIOR_ALPHA = 0.15


def ibm_sequential_cmap(name: str = "IBMwarm"):
    """Sequential colormap built from the IBM palette: yellow → orange → magenta.

    Perceptually monotonic in luminance (yellow brightest, magenta darkest),
    colorblind-safe, and tonally consistent with the IBM categorical palette
    used elsewhere in the manuscript. Use as a drop-in for viridis/cividis
    on sequential log-scale heatmaps.
    """
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        name,
        [IBM_PALETTE["yellow"], IBM_PALETTE["orange"], IBM_PALETTE["magenta"]],
    )
