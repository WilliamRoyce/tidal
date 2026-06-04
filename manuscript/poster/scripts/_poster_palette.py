r"""Modern-Cambridge poster palette + shared matplotlib style.

Mirrors the role assignment of ``scripts/figures/_palette.py`` but recolours to
the 2024 University of Cambridge identity (Johnson Banks refresh) so the poster
figures match the LaTeX theme (``theme/beamercolorthemecam.sty``).

Source of hex values: https://www.cam.ac.uk/brand-resources/guidelines/colour

Role assignment (kept identical to the report figures, recoloured):
- AMP  (amplification, +log A)  -> Cherry  (pink/magenta)
- SUP  (suppression,  -log A)  -> Crest   (coral/orange)
- structure / rules            -> Dark Blue / Warm Blue
- body text                    -> Slate 4
"""

from __future__ import annotations

# ---- Core Cambridge palette -------------------------------------------------
CAM = {
    "blue": "#8EE8D8",  # anchor / signature
    "light_blue": "#D1F9F1",
    "warm_blue": "#00BDB6",
    "dark_blue": "#133844",
    "slate4": "#232830",  # body text
    "slate2": "#B5BDC8",
}

# ---- Figure-coherent accent pair (secondary families) -----------------------
CHERRY = "#CD3572"  # amplification, propagating (pink/magenta)
CHERRY_DARK = "#911449"
WARM_CHERRY = "#E18AAC"  # amplification, non-propagating control (light cherry)
CREST = "#FD8153"  # suppression, propagating (coral/orange)
CREST_DARK = "#DD3025"
WARM_CREST = "#FFC392"  # suppression, non-propagating control (light crest)

# ---- Role aliases (match _palette.py interface) -----------------------------
AMP_COLOR = CHERRY
SUP_COLOR = CREST
PRIOR_COLOR = CAM["slate2"]
HIGHLIGHT_COLOR = CAM["warm_blue"]

AMP_ALPHA = 0.55
SUP_ALPHA = 0.55
PRIOR_ALPHA = 0.18

# ---- Shared matplotlib rcParams (sans-serif to match Lato; no usetex) -------
# Larger absolute sizes than the report figures: poster figures are viewed from
# a distance, and we do NOT use the report's stix/usetex serif preamble.
POSTER_RCPARAMS = {
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.sans-serif": ["Lato", "DejaVu Sans", "Arial"],
    "font.size": 16,
    "axes.labelsize": 18,
    "axes.titlesize": 20,
    "axes.edgecolor": CAM["slate4"],
    "axes.labelcolor": CAM["slate4"],
    "text.color": CAM["slate4"],
    "xtick.color": CAM["slate4"],
    "ytick.color": CAM["slate4"],
    "legend.fontsize": 16,
    "lines.linewidth": 2.0,
    "axes.linewidth": 1.0,
    # Math labels (T^2, R~ grad T) in a sans font consistent with Lato.
    "mathtext.fontset": "dejavusans",
    # Light-blue background so the figure blends into its light-blue poster box.
    "figure.facecolor": CAM["light_blue"],
    "savefig.facecolor": CAM["light_blue"],
    "axes.facecolor": CAM["light_blue"],
}


def apply_poster_style() -> None:
    """Apply the shared poster matplotlib rcParams."""
    import matplotlib as mpl

    mpl.rcParams.update(POSTER_RCPARAMS)
