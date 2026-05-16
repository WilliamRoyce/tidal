r"""Figure D §6 — multi-method agreement (App-C styled).

Two heatmaps side-by-side (`figure*`) of the pairwise relative
difference $|P_a - P_b|/\max(|P_a|, |P_b|)$ across backend pairs on
two representative theories. Log-scale viridis colormap; the IEEE
floor is included in the lower colour bound for anchoring.

Annotated `min / max / median` rel_diff per theory in the in-figure
text.

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
EPS_MACH = 2.220446049250313e-16

# Descriptive titles (the panel key is the JSON's repository-style
# identifier; the figure must show physics-descriptive names instead).
THEORY_TITLE: dict[str, str] = {
    "gertsenshtein": "Einstein–Maxwell Gertsenshtein",
    "coupled_scalars": "Coupled-scalar two-channel mixing",
    "gertsenshtein_proca": "Einstein–Maxwell + Proca-mass photon",
    "torsion_gertsenshtein": "Propagating-PGT torsion–Gertsenshtein",
    "torsion_gertsenshtein_nonminimal": "Nonminimal torsion–EM coupling",
}


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

    n_th = len(theories)
    fig, axes = plt.subplots(1, n_th, figsize=(5.2 * n_th, 4.4), squeeze=False)
    vmax = max(
        (p["rel_diff"] for p in pairs if p["rel_diff"] > 0),
        default=1e-2,
    )

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="white")

    for k, theory in enumerate(theories):
        ax = axes[0, k]
        m = _matrix(pairs, theory, backends)
        # Replace zeros with EPS_MACH for log-norm rendering.
        rendered = np.where(np.isnan(m), EPS_MACH, np.maximum(m, EPS_MACH))
        # Mask upper triangle (including diagonal) to avoid double-showing
        # the symmetric pair (a, b) ≡ (b, a) and the trivial self-comparisons.
        mask_upper = np.triu(np.ones_like(rendered, dtype=bool))
        rendered = np.ma.masked_array(rendered, mask=mask_upper)
        im = ax.imshow(rendered, norm="log", cmap=cmap, vmin=EPS_MACH, vmax=vmax)
        ax.set_xticks(range(len(backends)))
        ax.set_yticks(range(len(backends)))
        ax.set_xticklabels(backends, rotation=45, ha="right", fontsize=8)
        ax.set_yticklabels(backends, fontsize=8)
        ax.set_title(THEORY_TITLE.get(theory, theory), fontsize=10)

        # Annotate min / max / median in figure text
        theory_pairs = [p["rel_diff"] for p in pairs if p["theory"] == theory]
        if theory_pairs:
            mn, mx = min(theory_pairs), max(theory_pairs)
            md = float(np.median(theory_pairs))
            ax.text(
                0.02,
                -0.18,
                rf"min $= {mn:.1e}$,  median $= {md:.1e}$,  max $= {mx:.1e}$",
                transform=ax.transAxes,
                fontsize=8,
                va="top",
            )
        fig.colorbar(
            im, ax=ax, shrink=0.85, label="rel. diff." if k == n_th - 1 else ""
        )

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
