r"""Figure K.2 — YM--PGT atlas, up/down overlay on the IBM magenta/yellow pair.

Reads survey output from ``hpc_results/29890993/atlas_t5_2/`` (campaign tag
``t5_2`` corresponds to the YM--PGT union sector of ``results.tex``
\cref{NonMinimalLagrangian} — the union of the Bahamonde, Barker, and
Shapiro families on top of the EM--RC base, with eight BSM couplings
$\beta_{1,2,3}$, $\xi$, $\chi$, $\zeta_{1,2,3}$) and writes
``manuscript/figures/figK2_atlas_ympgt.pdf`` via TIDAL's atlas plotter
in ``overlay_up_down`` mode.

The overlay collapses the up- and down-face posteriors of each
coupling axis into a single panel, halving the panel count from 2N to
N.  For N=8 the layout is 2 x 4: top row axes 1--4, bottom row axes
5--8.  Each panel overlays the positive-sign face (``\uparrow``, IBM
magenta ``#dc267f``) and the negative-sign face (``\downarrow``, IBM
yellow ``#ffb000``) at alpha 0.5 — the same convention as the
manuscript's results-section corner plots (see
``scripts/figures/_palette.py``).

Reproducible by anyone with the survey directory; if the path is
missing, the script exits with a clear pull instruction.

Data source:
  hpc_results/29890993/atlas_t5_2/

Output:
  manuscript/figures/figK2_atlas_ympgt.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SURVEY = REPO_ROOT / "hpc_results" / "29890993" / "atlas_t5_2"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figK2_atlas_ympgt.pdf"

# IBM colorblind palette — same hex pair used by the results-section
# corner plots' AMP_COLOR / SUP_COLOR in `scripts/figures/_palette.py`.
# Magenta carries the positive-axis (\uparrow) overlay, yellow the
# negative-axis (\downarrow), so the K.2 atlas reads in the same
# visual language as the rest of the manuscript's corner plots.  Both
# are drawn at alpha 0.5 (`OVERLAY_ALPHA` in `_corner_style.py`) so
# overlapping support mixes to a darker tone.
COLOR_UP = "#dc267f"
COLOR_DOWN = "#ffb000"
COLOR_ALPHA = 0.5


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--survey-dir",
        type=Path,
        default=DEFAULT_SURVEY,
        help=f"Survey directory (default: {DEFAULT_SURVEY})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUT,
        help=f"Output PDF path (default: {DEFAULT_OUT})",
    )
    args = parser.parse_args()

    if not args.survey_dir.exists():
        msg = (
            f"error: {args.survey_dir} not found — pull the YM--PGT atlas "
            "campaign (tag t5_2) first via\n"
            "  bash scripts/hpc_shuttle.sh pull 29890993\n"
            "or pass --survey-dir to point at an existing YM--PGT survey."
        )
        print(msg, file=sys.stderr)
        sys.exit(1)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Import here so a missing tidal install doesn't blow up the
    # --help message above.
    from tidal.inference._atlas import plot_atlas

    plot_atlas(
        args.survey_dir,
        args.output,
        overlay_up_down=True,
        color_up=COLOR_UP,
        color_down=COLOR_DOWN,
        color_alpha=COLOR_ALPHA,
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
