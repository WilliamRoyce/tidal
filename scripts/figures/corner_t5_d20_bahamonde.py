r"""Corner plot for T5.0 Yang--Mills--PGT (Bahamonde 5D) — D2.0 v3 publication chains.

Amplification (29507332) and suppression (29471255) posteriors overlaid in the
App J \\cref{NestedScore} sign convention. Five-dimensional prior
$\\{\\beta_1, \\beta_2, \\beta_3, \\xi, \\delta_1\\}$.

Data sources:
  hpc_results/29507332/d20_bahamonde_amp_v3_pub/
  hpc_results/29471255/d20_bahamonde_sup_v3_pub/
Output:
  manuscript/figures/corner_t5_d20_bahamonde.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29507332" / "d20_bahamonde_amp_v3_pub"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29471255" / "d20_bahamonde_sup_v3_pub"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t5_d20_bahamonde.pdf"

PARAMS = ["beta1", "beta2", "beta3", "xi", "delta1"]
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "xi": r"$\xi$",
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
