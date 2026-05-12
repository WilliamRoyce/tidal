"""Figure D §3 — Forced-oscillator Pass 0+1 calibration.

Two-panel: (a) err_pass0 and err_combined vs ε on log-log axes with
the O(ε) and O(ε²) reference guides; (b) improvement ratio vs ε.

Data:   benchmark_results/canonical/forced_oscillator.json
Output: manuscript/figures/figD_forced_oscillator.pdf
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "forced_oscillator.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_forced_oscillator.pdf"
FLOOR = 1e-16


def _plot(data: dict, out_path: Path) -> None:
    rows = sorted(data["results"], key=operator.itemgetter("eps"))
    eps = np.array([r["eps"] for r in rows])
    err_p0 = np.array([max(r["err_pass0_vs_full"], FLOOR) for r in rows])
    err_c = np.array([max(r["err_combined_vs_full"], FLOOR) for r in rows])
    ratio = np.array([r["improvement_ratio"] for r in rows])
    slope = data["summary"].get("fitted_slope_err_vs_eps")

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

    ax = axes[0]
    ax.loglog(
        eps, err_p0, marker="o", ms=5, lw=1.0, color="#d62728", label="Pass 0 only"
    )
    ax.loglog(
        eps, err_c, marker="s", ms=5, lw=1.0, color="#1f77b4", label="Pass 0 + Pass 1"
    )
    # Reference guides O(ε) and O(ε²) anchored at max eps
    ax.loglog(
        eps,
        err_p0[-1] * (eps / eps[-1]) ** 1,
        ls=":",
        lw=0.8,
        color="#d62728",
        alpha=0.5,
    )
    ax.loglog(
        eps,
        err_c[-1] * (eps / eps[-1]) ** 2,
        ls=":",
        lw=0.8,
        color="#1f77b4",
        alpha=0.5,
    )
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("error vs full-$\\varepsilon$ modal")
    title = "(a) error vs $\\varepsilon$"
    if slope is not None:
        title += f" — slope ≈ {slope:.2f}"
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)

    ax = axes[1]
    ax.loglog(eps, ratio, marker="o", ms=5, lw=1.0, color="#2ca02c")
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel("Pass 0 / (Pass 0 + Pass 1)")
    ax.set_title("(b) improvement ratio", fontsize=10)
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
