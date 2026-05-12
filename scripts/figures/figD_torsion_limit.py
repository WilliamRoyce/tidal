r"""Figure D -- Einstein--Maxwell limit recovery from the torsion-extended theory.

Single-panel log--log plot of $|\\Delta P_\\mathrm{max}|$ vs $\\xi$, the
torsion-coupling parameter. The corresponding Einstein--Maxwell reference
$P_\\mathrm{max}^\\mathrm{EM}$ is annotated as a horizontal line. The
expected behaviour for a non-trivial limit is $|\\Delta| \\to 0$ as
$\\xi \\to 0$; the surveyed $h_\\times \\leftrightarrow a_x$ channel
exhibits a kinematic identity ($|\\Delta| = 0$ for all $\\xi$), which is
plotted at the machine-precision floor.

Data source:  benchmark_results/canonical/torsion_limit.json
Output:       manuscript/figures/figD_torsion_limit.pdf
Appendix ref: manuscript/sections/appendices/validation.tex (App D, calibration 2)
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "torsion_limit.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_torsion_limit.pdf"

FLOOR = 1e-16


def _plot(data: dict, out_path: Path) -> None:
    rows = sorted(data["xi_sweep"], key=operator.itemgetter("xi"))
    if not rows:
        msg = "no xi_sweep rows in JSON"
        raise ValueError(msg)
    xi = np.array([r["xi"] for r in rows])
    diff = np.array([r["abs_diff_P_max"] for r in rows])
    diff_plot = np.maximum(diff, FLOOR)
    em_pmax = data["em_reference"]["P_max"]

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.loglog(
        xi,
        diff_plot,
        marker="o",
        ms=5,
        color="#1f77b4",
        lw=1.0,
        label=r"$|\Delta P_{\rm max}|$",
    )
    ax.axhline(FLOOR, ls=":", lw=0.8, color="#999", label="machine-precision floor")
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel(r"$|P_{\rm max}^{\rm torsion} - P_{\rm max}^{\rm EM}|$")
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.set_title(rf"EM reference $P_{{\rm max}} = {em_pmax:.4f}$", fontsize=10)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    with args.data.open() as fh:
        data = json.load(fh)
    _plot(data, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
