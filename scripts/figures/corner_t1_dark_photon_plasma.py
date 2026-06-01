r"""Corner plot for T1 dark-photon-plasma — Stage A v3 publication chains.

Amplification (29205968) and suppression (29205982) posteriors overlaid in
the App J \\cref{NestedScore} sign convention (red = amplifying-region
likelihood, blue = suppressing-region likelihood). Four-dimensional prior
$\\{m_A^2, \\delta_m, \\xi, \\alpha_3\\}$.

Data sources:
  hpc_results/29205968/stage_a_amp_v3_pub/
  hpc_results/29205982/stage_a_sup_v3_pub/
Output:
  manuscript/figures/corner_t1_dark_photon_plasma.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Ensure sibling _corner_style.py is importable when invoked from anywhere.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29205968" / "stage_a_amp_v3_pub"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29205982" / "stage_a_sup_v3_pub"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t1_dark_photon_plasma.pdf"

PARAMS = ["mA2", "deltam", "xi", "alpha3"]
PARAM_LABELS = {
    "mA2": r"$m_A^2$",
    "deltam": r"$\delta_m$",
    "xi": r"$\xi$",
    "alpha3": r"$\alpha_3$",
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
