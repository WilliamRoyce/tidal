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
    fig, axes = plt.subplots(1, n_th, figsize=(5.2 * n_th, 4.2), squeeze=False)

    # vmin/vmax from lower-triangle (visible) cells only.  vmax was already
    # clipped this way, but vmin used to be hard-pinned to EPS_MACH ~ 2e-16,
    # which made the colorbar span ~14 decades of empty range below the
    # actual data and left every visible cell sitting in the upper half of
    # the viridis ramp.
    visible_diffs: list[float] = []
    for theory in theories:
        for p in pairs:
            if p["theory"] != theory:
                continue
            i = backends.index(p["backend_a"])
            j = backends.index(p["backend_b"])
            if i > j and p["rel_diff"] > 0:
                visible_diffs.append(p["rel_diff"])
    vmax = max(visible_diffs, default=1e-2)
    vmin = min(visible_diffs, default=EPS_MACH)

    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad(color="white")

    # The lower-triangle mask leaves row 0 (modal) and the last column
    # (scipy_DOP853) with zero visible cells.  Clip to the populated
    # 5x5 sub-matrix: rows backends[1:], cols backends[:-1].
    row_labels = backends[1:]
    col_labels = backends[:-1]
    n_row = len(row_labels)
    n_col = len(col_labels)

    im = None
    for k, theory in enumerate(theories):
        ax = axes[0, k]
        m_full = _matrix(pairs, theory, backends)
        rendered_full = np.where(
            np.isnan(m_full), EPS_MACH, np.maximum(m_full, EPS_MACH)
        )
        mask_upper_full = np.triu(np.ones_like(rendered_full, dtype=bool))
        # Clip: drop row 0 (modal) and last col (DOP853).
        sub = rendered_full[1:, :-1]
        sub_mask = mask_upper_full[1:, :-1]
        rendered = np.ma.masked_array(sub, mask=sub_mask)
        im = ax.imshow(rendered, norm="log", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(n_col))
        ax.set_yticks(range(n_row))
        ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=8)
        if k == 0:
            ax.set_yticklabels(row_labels, fontsize=8)
        else:
            ax.set_yticklabels([])
        ax.set_title(THEORY_TITLE.get(theory, theory), fontsize=10)

        # Outline each populated lower-triangle cell with a Rectangle patch
        # so masked upper-right cells get no surrounding grid lines.
        from matplotlib.patches import Rectangle

        for i_sub in range(n_row):
            for j_sub in range(n_col):
                if i_sub >= j_sub:  # populated (lower-triangle in sub-coords)
                    ax.add_patch(
                        Rectangle(
                            (j_sub - 0.5, i_sub - 0.5),
                            1,
                            1,
                            fill=False,
                            edgecolor="black",
                            linewidth=0.5,
                        )
                    )

        # Drop the rectangular outer frame so the triangle reads as a triangle.
        for spine in ax.spines.values():
            spine.set_visible(False)

    # Single shared colorbar across both panels.
    if im is not None:
        fig.colorbar(im, ax=axes[0, :], shrink=0.85, label="rel. diff.")

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
