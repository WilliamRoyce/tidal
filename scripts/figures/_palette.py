r"""Central IBM colorblind-friendly palette for §4 figures.

All figure scripts in `scripts/figures/` should import their colors from this
module rather than hardcoding hex values, so the palette stays consistent
across the manuscript.

Palette source: https://www.color-hex.com/color-palette/1044488
(also cited by Legner–Handley–Barker 2507.09228 "SL would like to thank the
IBM color blind friendly palette for the color scheme used in the figures").

Role assignment for our overlay-corner figures:
- AMP (amplification posterior, +log A direction)  → MAGENTA
- SUP (suppression posterior, -log A direction)    → YELLOW
- PRIOR (optional third overlay layer)             → PURPLE (indigo) at low alpha
- HIGHLIGHT (carrier-parameter accent in tables)   → ORANGE
- RESERVED (future fourth overlay, e.g. localized) → BLUE

The yellow/magenta pairing follows Legner–Handley–Barker (2507.09228) and
gives stronger contrast than the previous orange/blue pair; the purple
(indigo) prior layer reads as a distinct underlay behind both posteriors.
"""

from __future__ import annotations

# Five-color IBM colorblind palette (hex, RGB).
IBM_PALETTE = {
    "yellow": "#ffb000",
    "orange": "#fe6100",
    "magenta": "#dc267f",
    "purple": "#785ef0",
    "blue": "#648fff",
}

AMP_COLOR = IBM_PALETTE["magenta"]
SUP_COLOR = IBM_PALETTE["yellow"]
PRIOR_COLOR = IBM_PALETTE["purple"]
HIGHLIGHT_COLOR = IBM_PALETTE["orange"]
RESERVED_COLOR = IBM_PALETTE["blue"]

# Default alphas tuned for the overlay-corner visual: amp / sup carry the
# primary signal; prior (when drawn) is a faint reference cloud.
AMP_ALPHA = 0.5
SUP_ALPHA = 0.5
PRIOR_ALPHA = 0.15


# 2024 Cambridge brand palette (Johnson Banks refresh).
# Source: cam.ac.uk/brand-resources/guidelines/colour.
# Used by the talk-deck figure variants so the corner-plot accents match
# the slide accents; manuscript figures continue to use IBM_PALETTE.
CAMBRIDGE_PALETTE = {
    # Core palette
    "blue": "#8EE8D8",
    "light_blue": "#D1F9F1",
    "warm_blue": "#00BDB6",
    "dark_blue": "#133844",
    "slate4": "#232830",
    # Crest family (orange / coral)
    "light_crest": "#FFE2C8",
    "warm_crest": "#FFC392",
    "crest": "#FD8153",
    "dark_crest": "#DD3025",
    # Cherry family (pink / red)
    "light_cherry": "#F2CAD8",
    "warm_cherry": "#E18AAC",
    "cherry": "#CD3572",
    "dark_cherry": "#911449",
    # Purple family
    "light_purple": "#F2ECF8",
    "warm_purple": "#D1B7EB",
    "purple": "#A368DF",
    "dark_purple": "#681FB1",
    # Indigo family
    "light_indigo": "#EBEDFB",
    "warm_indigo": "#B0B9F1",
    "indigo": "#5366E0",
    "dark_indigo": "#29347A",
}
