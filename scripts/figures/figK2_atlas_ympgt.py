r"""Figure K.2 — YM--PGT atlas as cubed-sphere methodology illustration.

Reads survey output from ``hpc_results/29890993/atlas_t5_2/`` (campaign tag
``t5_2`` corresponds to the YM--PGT union sector of ``results.tex``
\cref{NonMinimalLagrangian} — the union of the Bahamonde, Barker, and
Shapiro families on top of the EM--RC base, with eight BSM couplings
$\beta_{1,2,3}$, $\xi$, $\chi$, $\zeta_{1,2,3}$) and writes
``manuscript/figures/figK2_atlas_ympgt.pdf`` via TIDAL's atlas plotter
with the 4x4 block-pair layout (``layout_cols=4``).

The 4x4 layout: top two rows hold axes 1--4 (up over down, column-wise);
bottom two rows hold axes 5--8 similarly. Generalises the default 2 x N
strip for high-N theories so the figure fits one page.

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
LAYOUT_COLS = 4  # 4x4 block-pair layout for N=8


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

    plot_atlas(args.survey_dir, args.output, layout_cols=LAYOUT_COLS)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
