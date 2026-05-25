r"""Corner plot for the non-propagating-torsion (NP) control sector (8D).

Amplification (job 29700462) and suppression posteriors overlaid in the
App J \cref{NestedScore} sign convention. Eight-dimensional sub-space
$\{\beta_{1\!-\!3}, \delta_1, \chi, \zeta_{1\!-\!3}\}$ obtained from the
nine-parameter Yang--Mills--PGT non-minimal sector by pinning the torsion
kinetic coefficient $\xi = 0$, leaving torsion algebraically determined by
the Cartan equation (Einstein--Cartan-like configuration with the
torsion--electromagnetism vertices retained).

The chain returns log Z_+ = +10.98 ± 0.15, log Z_- = +2.67 ± 0.10
(log B = +8.31): essentially the same amplification signal as the full
propagating theory. This is the load-bearing diagnostic against the
torsion-ghost interpretation of the amplification — with no torsion kinetic
operator, there is no torsion ghost mode, yet the amplification persists.

Data sources:
  hpc_results/29700462/np_amp_v1/
  hpc_results/29700462/np_sup_v1/
Output:
  manuscript/figures/corner_np.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29700462" / "np_amp_v1"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29700462" / "np_sup_v1"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_np.pdf"

PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "delta1",
    "chi",
    "zeta1",
    "zeta2",
    "zeta3",
]
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "delta1": r"$\delta_1$",
    "chi": r"$\chi$",
    "zeta1": r"$\zeta_1$",
    "zeta2": r"$\zeta_2$",
    "zeta3": r"$\zeta_3$",
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
