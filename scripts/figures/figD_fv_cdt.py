"""Figure D §5 — FV ↔ TorsionCDT bit-exact equivalence.

Two-panel (figure*): (a) P_max(point) bar comparison for FV and CDT
across the five canonical points; (b) log-scale relative difference
|ΔP_max|/P_max per point with the IEEE round-off floor annotated.

Data:   benchmark_results/canonical/fv_cdt_equivalence.json
Output: manuscript/figures/figD_fv_cdt.pdf
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "fv_cdt_equivalence.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_fv_cdt.pdf"

IEEE = 2.22e-16
FLOOR = 1e-18


def _plot(data: dict, out_path: Path) -> None:
    rows = [r for r in data["results"] if r.get("ok")]
    if not rows:
        msg = "no fv_cdt rows"
        raise ValueError(msg)
    labels = [r["label"] for r in rows]
    p_fv = np.array([r["P_max_FV"] for r in rows])
    p_cdt = np.array([r["P_max_CDT"] for r in rows])
    rel = np.array([max(r["rel_diff"], FLOOR) for r in rows])

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.8))

    ax = axes[0]
    x = np.arange(len(labels))
    width = 0.4
    ax.bar(x - width / 2, p_fv, width, label="FV (10 fields)", color="#1f77b4")
    ax.bar(x + width / 2, p_cdt, width, label="CDT (18 fields)", color="#d62728")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel(r"$P_{\rm max}$")
    ax.set_title("(a) $P_{\\rm max}$ at canonical points", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(visible=True, axis="y", ls=":", alpha=0.4)

    ax = axes[1]
    ax.bar(x, rel, color="#2ca02c", alpha=0.8)
    ax.axhline(IEEE, ls="--", lw=0.9, color="#888", label="IEEE-754 floor")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel(r"$|\Delta P_{\rm max}| / P_{\rm max}$")
    ax.set_title("(b) cross-formulation relative difference", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(visible=True, axis="y", which="both", ls=":", alpha=0.4)
    ax.set_ylim(IEEE / 2, max(1e-10, rel.max() * 5))

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
