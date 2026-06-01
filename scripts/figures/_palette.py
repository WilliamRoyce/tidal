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
