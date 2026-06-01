r"""Corner plot for the non-propagating-torsion χ-closure control (17D).

Amplification (job 29705560) and suppression posteriors overlaid in the
App J \cref{NestedScore} sign convention. Seventeen-dimensional sub-space
$\{\beta_{1\!-\!3}, \delta_1, \zeta_{1\!-\!3}, \chi_{1\!-\!10}\}$
obtained from the eighteen-parameter complete parity-even (χ-closure)
sector by pinning the torsion kinetic coefficient $\xi = 0$ — the
χ-closure analogue of the NP control of \cref{CornerNP}.

The chain returns log Z_+ = +9.93 ± 0.08, log Z_- = -3.28 ± 0.07
(log B = +13.21): statistically indistinguishable from the
ξ-free T7 v2 chain (+10.30 ± 0.08), confirming that the χ-closure
amplification is driven by the χ_{1..10} R̃×∇T coupling sector
rather than by the propagating torsion kinetic mode.

Data sources:
  hpc_results/29705560/np_ceven_amp_v1/
  hpc_results/29705560/np_ceven_sup_v1/
Output:
  manuscript/figures/corner_np_ceven.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_amp_v1"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29705560" / "np_ceven_sup_v1"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_np_ceven.pdf"

PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "delta1",
    "zeta1",
    "zeta2",
    "zeta3",
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi5",
    "chi6",
    "chi7",
    "chi8",
    "chi9",
    "chi10",
]
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "delta1": r"$\delta_1$",
    "zeta1": r"$\zeta_1$",
    "zeta2": r"$\zeta_2$",
    "zeta3": r"$\zeta_3$",
    "chi1": r"$\chi_1$",
    "chi2": r"$\chi_2$",
    "chi3": r"$\chi_3$",
    "chi4": r"$\chi_4$",
    "chi5": r"$\chi_5$",
    "chi6": r"$\chi_6$",
    "chi7": r"$\chi_7$",
    "chi8": r"$\chi_8$",
    "chi9": r"$\chi_9$",
    "chi10": r"$\chi_{10}$",
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
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
