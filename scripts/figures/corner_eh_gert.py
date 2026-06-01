r"""Corner plot for Euler--Heisenberg + Gertsenshtein (2D, $\rho$, $\sigma$).

Amplification (job 29700083) and suppression posteriors overlaid in the
App J \cref{NestedScore} sign convention. Two-dimensional prior over the
two physical QED-EFT F$^4$ couplings:
  $\rho$ — coefficient of $(F\cdot F)^2$
  $\sigma$ — coefficient of $(F\cdot\tilde F)^2$
with the QED reference point $\sigma/\rho = 7/16$ (Adler 1971, Dunne 2012).

The chain returns log Z_+ = +0.184, log Z_- = -0.185 (log B = +0.37):
clean, symmetric null. Euler--Heisenberg F$^4$ corrections are inert for
Gertsenshtein conversion at the surveyed field strength $B_0 = 10^{-2}$.

The single F·F̃ Pontryagin topological term is correctly identified as a
total derivative at the symbolic level by the Wolfram E-L pipeline
(gertsenshtein_eh_top.json: lambdatop absent from the equations of motion,
27 H-terms identical to the base EH derivation); the symbolic-derivation
outcome is the publishable validation and no extended-prior chain is needed.

Data sources:
  hpc_results/29700083/eh_gert_amp_v1/
  hpc_results/29700083/eh_gert_sup_v1/
Output:
  manuscript/figures/corner_eh_gert.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29700083" / "eh_gert_amp_v1"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29700083" / "eh_gert_sup_v1"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_eh_gert.pdf"

PARAMS = ["rho", "sigma"]
PARAM_LABELS = {
    "rho": r"$\rho$",
    "sigma": r"$\sigma$",
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
