"""Figure D §6 — multi-method agreement heatmap.

Two-panel (figure*): (a) Backend-pair relative-difference heatmap per
theory; (b) summary bar chart of pairwise rel_diff distribution
plus the median and max across all pairs.

Data:   benchmark_results/canonical/multi_method.json
Output: manuscript/figures/figD_multi_method.pdf
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "multi_method.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_multi_method.pdf"
FLOOR = 1e-16


def _matrix(pairs: list[dict], theory: str, backends: list[str]) -> np.ndarray:
    n = len(backends)
    m = np.full((n, n), np.nan)
    for p in pairs:
        if p["theory"] != theory:
            continue
        i = backends.index(p["backend_a"])
        j = backends.index(p["backend_b"])
        m[i, j] = m[j, i] = p["rel_diff"]
    for i in range(n):
        m[i, i] = 0.0
    return m


def _plot(data: dict, out_path: Path) -> None:
    backends = data["metadata"]["parameters"]["backends"]
    theories = data["metadata"]["parameters"]["theories"]
    pairs = data["pairwise"]

    n_theories = len(theories)
    fig, axes = plt.subplots(
        1, n_theories + 1, figsize=(5.0 * n_theories + 4, 4.0), squeeze=False
    )
    for k, theory in enumerate(theories):
        ax = axes[0, k]
        m = _matrix(pairs, theory, backends)
        im = ax.imshow(
            np.where(np.isnan(m), FLOOR, np.maximum(m, FLOOR)),
            norm="log",
            cmap="viridis",
            vmin=FLOOR,
            vmax=max(1e-2, np.nanmax(m) if np.any(~np.isnan(m)) else 1e-2),
        )
        ax.set_xticks(range(len(backends)))
        ax.set_yticks(range(len(backends)))
        ax.set_xticklabels(backends, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(backends, fontsize=7)
        ax.set_title(theory, fontsize=10)
        fig.colorbar(im, ax=ax, label="rel. diff." if k == n_theories - 1 else "")

    # Right-most panel: distribution of pairwise rel_diff
    ax = axes[0, -1]
    rels = [max(p["rel_diff"], FLOOR) for p in pairs]
    ax.hist(np.log10(rels), bins=20, color="#1f77b4", alpha=0.8)
    ax.set_xlabel(r"$\log_{10}(\Delta P / P)$")
    ax.set_ylabel("pair count")
    ax.set_title("(c) pairwise distribution", fontsize=10)
    ax.grid(visible=True, ls=":", alpha=0.4)

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
