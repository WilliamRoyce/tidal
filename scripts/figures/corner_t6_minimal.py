r"""Corner plot for T6-minimal parity-odd sector (4D) — v3 INTR landscape chains.

Amplification and suppression posteriors (both from job 29515407) overlaid in
the App J \\cref{NestedScore} sign convention. Four-dimensional prior
$\\{\\beta_1, \\beta_2, \\beta_3, d_{21}\\}$ — the minimal parity-odd sub-sector.

Data sources:
  hpc_results/29515407/t6_minimal_amp_v3/
  hpc_results/29515407/t6_minimal_sup_v3/
Output:
  manuscript/figures/corner_t6_minimal.pdf
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _corner_style import COLUMN_WIDTH, overlay_corner

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_AMP = REPO_ROOT / "hpc_results" / "29515407" / "t6_minimal_amp_v3"
DEFAULT_SUP = REPO_ROOT / "hpc_results" / "29515407" / "t6_minimal_sup_v3"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "corner_t6_minimal.pdf"

PARAMS = ["beta1", "beta2", "beta3", "d21"]
PARAM_LABELS = {
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "d21": r"$d_{21}$",
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
