r"""Figure D -- Energy conservation audit across the example library.

Bar chart of $|dE/E|_\\mathrm{max}$ per example, with the de-Sitter outlier
(physical Hubble friction) annotated. The conservation audit calibrates
the structural integrity of the E-L velocity pipeline; failing examples
indicate pipeline regressions rather than per-example precision floors.

Data source:  benchmark_results/canonical/conservation_audit.json
Output:       manuscript/figures/figD_conservation_audit.pdf
Appendix ref: manuscript/sections/appendices/validation.tex (App D, calibration 5)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "conservation_audit.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_conservation_audit.pdf"

FLOOR = 1e-16


def _plot(data: dict, out_path: Path) -> None:
    examples = data["examples"]
    labels = [r["label"] for r in examples]
    drifts = [r.get("max_abs_dE_over_E") or 0.0 for r in examples]
    ok = [bool(r.get("ok")) for r in examples]

    fig, ax = plt.subplots(figsize=(max(5.0, 0.7 * len(labels)), 3.6))
    x = np.arange(len(labels))
    colors = ["#1f77b4" if o else "#d62728" for o in ok]
    plotted = [max(d, FLOOR) for d in drifts]
    ax.bar(x, plotted, color=colors, alpha=0.85)
    ax.set_yscale("log")
    ax.axhline(
        1e-3, ls="--", lw=0.8, color="#999", label=r"default threshold $10^{-3}$"
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel(r"$|dE/E|_{\rm max}$")
    ax.set_ylim(FLOOR, 1e-1)
    ax.grid(visible=True, axis="y", which="both", ls=":", alpha=0.4)
    ax.legend(frameon=False, fontsize=9, loc="upper left")

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
