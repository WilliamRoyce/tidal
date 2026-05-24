r"""Corner plot for T4 Ricci--electromagnetism cross-term — D1 v3 publication chains.

Amplification (29189748) and suppression (29189761) posteriors overlaid in the
App J \\cref{NestedScore} sign convention. Four-dimensional prior
$\\{\\alpha_1, \\alpha_2, \\alpha_3, \\delta_1\\}$.

Data sources:
  hpc_results/29189748/d1_amp_v3/
  hpc_results/29189761/d1_sup_v3/
Output:
  manuscript/figures/corner_t4_ricci_em.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29189748" / "d1_amp_v3"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29189761" / "d1_sup_v3"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t4_ricci_em.pdf"

PARAMS = ["alpha1", "alpha2", "alpha3", "delta1"]
PARAM_LABELS = {
    "alpha1": r"$\alpha_1$",
    "alpha2": r"$\alpha_2$",
    "alpha3": r"$\alpha_3$",
    "delta1": r"$\delta_1$",
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--amp", type=Path, default=DEFAULT_AMP)
    parser.add_argument("--sup", type=Path, default=DEFAULT_SUP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    overlay_corner(
        amp_chains_dir=args.amp,
        sup_chains_dir=args.sup,
        params=PARAMS,
        param_labels=PARAM_LABELS,
        out_path=args.output,
        fig_width=COLUMN_WIDTH,
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
