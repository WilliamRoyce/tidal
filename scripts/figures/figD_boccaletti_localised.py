"""Figure D §2 — Localised Boccaletti path-integrated calibration.

Two-panel (figure*): heatmap of P_sim(Bpeak, R) on the left; overlay
of P_sim vs analytic at fixed R on the right with a residual band.

Data:   benchmark_results/canonical/boccaletti_localised.json
Output: manuscript/figures/figD_boccaletti_localised.pdf
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_localised.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_localised.pdf"


def _plot(data: dict, out_path: Path) -> None:
    rows = [r for r in data.get("results", []) if r.get("ok")]
    if not rows:
        msg = "no localised rows"
        raise ValueError(msg)
    bps = sorted({r["Bpeak"] for r in rows})
    rs = sorted({r["R"] for r in rows})
    p_sim = np.full((len(rs), len(bps)), np.nan)
    p_ana = np.full_like(p_sim, np.nan)
    for r in rows:
        i = rs.index(r["R"])
        j = bps.index(r["Bpeak"])
        p_sim[i, j] = r["P_sim"]
        p_ana[i, j] = r["P_analytic"]
    residual = p_sim - p_ana

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.8))

    ax = axes[0]
    im = ax.imshow(
        residual,
        origin="lower",
        aspect="auto",
        extent=(bps[0], bps[-1], rs[0], rs[-1]),
        cmap="RdBu_r",
        vmin=-max(abs(np.nanmin(residual)), abs(np.nanmax(residual)), 1e-6),
        vmax=+max(abs(np.nanmin(residual)), abs(np.nanmax(residual)), 1e-6),
    )
    ax.set_xlabel(r"$B_{\rm peak}$")
    ax.set_ylabel(r"$R$")
    ax.set_title("(a) residual P_sim − P_ana", fontsize=10)
    fig.colorbar(im, ax=ax, label="residual")

    # right-panel: line overlay at the middle R
    r_idx = len(rs) // 2
    r0 = rs[r_idx]
    ax = axes[1]
    ax.plot(
        bps,
        p_ana[r_idx, :],
        ls="--",
        lw=1.0,
        color="#666",
        label=r"$\sin^2(\kappa B_{\rm peak} R \sqrt{\pi/2})$",
    )
    ax.plot(
        bps, p_sim[r_idx, :], marker="o", ms=5, lw=1.0, color="#1f77b4", label="TIDAL"
    )
    ax.set_xlabel(r"$B_{\rm peak}$")
    ax.set_ylabel("P")
    ax.set_title(f"(b) line cut at $R={r0}$", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
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
