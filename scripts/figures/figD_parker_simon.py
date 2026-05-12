"""Figure D §4 — Parker-Simon higher-derivative perturbative calibration.

Same structural plot as forced_oscillator: error vs ε with O(ε²) guide,
plus the improvement ratio. The correction term is biharmonic (∂_x⁴)
which gives a distinct k-mode scaling from the mass shift.

Data:   benchmark_results/canonical/parker_simon_flrw.json
Output: manuscript/figures/figD_parker_simon.pdf
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "parker_simon_flrw.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_parker_simon.pdf"
FLOOR = 1e-16


def _plot(data: dict, out_path: Path) -> None:
    rows = sorted(data["results"], key=operator.itemgetter("eps"))
    eps = np.array([r["eps"] for r in rows])
    err_p0 = np.array([max(r["err_pass0_vs_full"], FLOOR) for r in rows])
    err_c = np.array([max(r["err_combined_vs_full"], FLOOR) for r in rows])
    slope = data["summary"].get("fitted_slope_err_vs_eps")

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.loglog(
        eps, err_p0, marker="o", ms=5, lw=1.0, color="#d62728", label="Pass 0 only"
    )
    ax.loglog(
        eps, err_c, marker="s", ms=5, lw=1.0, color="#1f77b4", label="Pass 0 + Pass 1"
    )
    ax.loglog(
        eps,
        err_c[-1] * (eps / eps[-1]) ** 2,
        ls=":",
        lw=0.8,
        color="#1f77b4",
        alpha=0.5,
        label=r"$\mathcal{O}(\varepsilon^2)$ guide",
    )
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("error vs full-$\\varepsilon$ modal")
    title = r"$\varepsilon \cdot \partial_x^4 \phi$ correction"
    if slope is not None:
        title += f" — fitted slope ≈ {slope:.2f}"
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)

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
