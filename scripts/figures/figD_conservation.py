"""Figure D §8 — energy-conservation audit across the example library.

Single-panel log-scale bar chart of |dE/E|_max per example with the
default 1e-3 conservation threshold drawn for reference. Failed
examples (e.g. due to missing IC params) are shown in red.

Data:   benchmark_results/canonical/conservation_audit_full.json
Output: manuscript/figures/figD_conservation.pdf
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "conservation_audit_full.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_conservation.pdf"
FLOOR = 1e-17


def _plot(data: dict, out_path: Path) -> None:
    rows = data["results"]
    labels = [r["label"] for r in rows]
    drifts = [r.get("max_abs_dE_over_E") or 0.0 for r in rows]
    ok = [bool(r.get("ok")) for r in rows]

    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(7.0, 0.55 * n), 4.2))
    x = np.arange(n)
    colors = ["#1f77b4" if o else "#d62728" for o in ok]
    plotted = [max(d, FLOOR) for d in drifts]
    ax.bar(x, plotted, color=colors, alpha=0.85)
    ax.set_yscale("log")
    ax.axhline(
        1e-3, ls="--", lw=0.9, color="#444", label=r"default threshold $10^{-3}$"
    )
    ax.axhline(2.22e-16, ls=":", lw=0.7, color="#888", label="IEEE-754 floor")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(r"$|dE/E|_{\rm max}$")
    ax.set_ylim(FLOOR, 1.0)
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
