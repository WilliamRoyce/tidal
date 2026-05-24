r"""Corner plot for T7 complete-even sector (18D) — v1 landscape partial chains.

Amplification and suppression posteriors (both from job 29588680, partial
chains under PolyChord TIMEOUT) overlaid in the App J \\cref{NestedScore}
sign convention. Eighteen-dimensional prior over the complete parity-even
sector with the full $\\chi_1$--$\\chi_{10}$ R̃×∂T basis; 74% of dead points
at the soft-floor noise band, geometric incompatibility under the plane-wave
benchmark.

Data sources:
  hpc_results/29588680/t7_amp_v1/  (partial chain)
  hpc_results/29588680/t7_sup_v1/  (partial chain)
Output:
  manuscript/figures/corner_t7.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29588680" / "t7_amp_v1"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29588680" / "t7_sup_v1"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t7.pdf"

PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "xi",
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
    "xi": r"$\xi$",
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
