r"""Corner plot for T9 complete-even + full $\xi$-kinetic closure (32D) — v2 rescue chains.

Amplification and suppression posteriors (job 29694142, t_end=1 rescue under
the overflow-safe horizon) overlaid in the App J \cref{NestedScore} sign
convention. Thirty-two-dimensional prior over the complete parity-even
sector with the full $\xi_1$--$\xi_{16}$ kinetic-invariant basis ($\xi_{11}$
vanishes in plane-wave reduction; 15 surviving) on top of the complete
$\chi_1$--$\chi_{10}$ R̃×∂T basis.

The v2 rescue returns log Z_+ = −0.18, log Z_- = −9.61 (log B = +9.43);
both log Z negative indicates partial floor contamination from the
extended $\xi_{11}$–$\xi_{16}$ sub-sector, so the Bayes factor is a LOWER
BOUND on the genuine amplification of the parity-even kinetic closure.

T9 closes the $\xi_1$--$\xi_{16}$ kinetic-sector gap (Barker's single $\xi$
$=$ LC$\{\xi_6, \xi_7, \xi_{12}, \xi_{16}\}$ is now the intersection of the
full surveyed basis rather than the survey scope).

Data sources:
  hpc_results/29694142/t9_amp_v2/
  hpc_results/29694142/t9_sup_v2/
Output:
  manuscript/figures/corner_t9.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29694142" / "t9_amp_v2"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29694142" / "t9_sup_v2"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t9.pdf"

# Parameter order matches the derivation pipeline's column convention
# (verified against examples/torsion_gertsenshtein/theory_complete_even_full_xi.toml).
# xi11 vanishes in the plane-wave reduction and is omitted from the
# chain columns; the remaining 32 = 3 (beta) + 15 (xi minus xi11) + 1 (delta1)
# + 3 (zeta) + 10 (chi) parameters appear in the order below.
PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "xi1",
    "xi2",
    "xi3",
    "xi4",
    "xi5",
    "xi6",
    "xi7",
    "xi8",
    "xi9",
    "xi10",
    "xi12",
    "xi13",
    "xi14",
    "xi15",
    "xi16",
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
    "xi1": r"$\xi_1$",
    "xi2": r"$\xi_2$",
    "xi3": r"$\xi_3$",
    "xi4": r"$\xi_4$",
    "xi5": r"$\xi_5$",
    "xi6": r"$\xi_6$",
    "xi7": r"$\xi_7$",
    "xi8": r"$\xi_8$",
    "xi9": r"$\xi_9$",
    "xi10": r"$\xi_{10}$",
    "xi12": r"$\xi_{12}$",
    "xi13": r"$\xi_{13}$",
    "xi14": r"$\xi_{14}$",
    "xi15": r"$\xi_{15}$",
    "xi16": r"$\xi_{16}$",
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
