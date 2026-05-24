r"""Corner plot for T6 full parity-odd sector (20D) — v4 landscape partial chains.

Amplification and suppression posteriors (both from job 29567416, partial
chains under PolyChord TIMEOUT) overlaid in the App J \\cref{NestedScore}
sign convention. Twenty-dimensional prior over the full parity-odd sector;
both posteriors floor-bound (100% of dead points at $\\log L \\approx -100$),
geometrically forbidden by the plane-wave benchmark.

Data sources:
  hpc_results/29567416/t6_amp_v4/  (partial chain)
  hpc_results/29567416/t6_sup_v4/  (partial chain)
Output:
  manuscript/figures/corner_t6_full.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29567416" / "t6_amp_v4"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29567416" / "t6_sup_v4"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t6_full.pdf"

PARAMS = [
    "beta1",
    "beta2",
    "beta3",
    "xi",
    "delta1",
    "chi",
    "zeta1",
    "zeta2",
    "zeta3",
    "d14",
    "d15",
    "d17",
    "d19",
    "d20",
    "d21",
    "zt1",
    "zt2",
    "zt3",
    "zt5",
    "zt6",
]
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "xi": r"$\xi$",
    "delta1": r"$\delta_1$",
    "chi": r"$\chi$",
    "zeta1": r"$\zeta_1$",
    "zeta2": r"$\zeta_2$",
    "zeta3": r"$\zeta_3$",
    "d14": r"$d_{14}$",
    "d15": r"$d_{15}$",
    "d17": r"$d_{17}$",
    "d19": r"$d_{19}$",
    "d20": r"$d_{20}$",
    "d21": r"$d_{21}$",
    "zt1": r"$\tilde\zeta_1$",
    "zt2": r"$\tilde\zeta_2$",
    "zt3": r"$\tilde\zeta_3$",
    "zt5": r"$\tilde\zeta_5$",
    "zt6": r"$\tilde\zeta_6$",
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
