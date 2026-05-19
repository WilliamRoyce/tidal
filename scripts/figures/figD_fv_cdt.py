r"""Figure D §5 — FV ↔ TorsionCDT cross-formulation equivalence.

Single-panel bar chart of $|\Delta P_{\max}|/P_{\max}$ across the
canonical parameter points, with the IEEE-754 round-off floor drawn
for reference. The cross-formulation equivalence test is bit-exact
at every point — the criterion is that all bars sit at or just above
the floor.

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
EPS_MACH = 2.220446049250313e-16


def _plot(data: dict, out_path: Path) -> None:
    rows = [r for r in data["results"] if r.get("ok")]
    labels = [r["label"] for r in rows]
    rels = np.array([max(r["rel_diff"], EPS_MACH / 10.0) for r in rows])
    max_rel = float(rels.max())

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    x = np.arange(len(labels))
    ax.bar(
        x,
        rels,
        color="#2ca02c",
        alpha=0.85,
        label=rf"$|\Delta P_{{\max}}|/P_{{\max}}$  (max $= {max_rel:.2e}$)",
    )
    ax.axhline(
        EPS_MACH,
        ls="--",
        lw=1.0,
        color="#444",
        label=rf"IEEE-754 floor $\varepsilon_{{\mathrm{{mach}}}} = {EPS_MACH:.2e}$",
    )
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel(r"$|\Delta P_{\max}|/P_{\max}$")
    ax.set_ylim(EPS_MACH / 3.0, max(max_rel * 5.0, 1e-10))
    ax.grid(visible=True, axis="y", which="both", ls=":", alpha=0.3)
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
