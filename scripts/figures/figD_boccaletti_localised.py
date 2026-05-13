r"""Figure D §2 — path-integrated Boccaletti reproduction (App-C styled).

Two-panel `figure*` showing the validation of the path-integrated
Boccaletti kernel on a Gaussian-localised B-field profile:

  (a) Residual heatmap $P_\mathrm{sim} - P_\mathrm{analytic}$ across
      the $(B_\mathrm{peak}, R)$ parameter grid, diverging colormap
      centred at zero.
  (b) Universal-collapse plot: all 40 simulation points plotted against
      the dimensionless variable $s \equiv \kappa\,B_\mathrm{peak}\,R\,
      \sqrt{\pi/2}$, overlaid with the analytic prediction $\sin^2(s)$.

The collapse plot is the strongest visual demonstration that the
path-integrated kernel is correct: all simulation points should fall
on the universal curve regardless of their individual
$(B_\mathrm{peak}, R)$ values.

Data:   benchmark_results/canonical/boccaletti_localised.json
Output: manuscript/figures/figD_boccaletti_localised.pdf
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_localised.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_localised.pdf"
KAPPA = 1.0


def _plot(data: dict, out_path: Path) -> None:
    rows = [r for r in data["results"] if r.get("ok")]
    bps = sorted({r["Bpeak"] for r in rows})
    rs = sorted({r["R"] for r in rows})
    nbp, nr = len(bps), len(rs)
    p_sim = np.full((nr, nbp), np.nan)
    p_ana = np.full_like(p_sim, np.nan)
    for r in rows:
        i = rs.index(r["R"])
        j = bps.index(r["Bpeak"])
        p_sim[i, j] = r["P_sim"]
        p_ana[i, j] = r["P_analytic"]
    residual = p_sim - p_ana

    s_vals = np.array(
        [KAPPA * r["Bpeak"] * r["R"] * math.sqrt(math.pi / 2.0) for r in rows]
    )
    p_sim_flat = np.array([r["P_sim"] for r in rows])

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))

    # Panel (a): residual heatmap
    ax = axes[0]
    vmax = max(abs(np.nanmin(residual)), abs(np.nanmax(residual)), 1e-6)
    im = ax.imshow(
        residual,
        origin="lower",
        aspect="auto",
        extent=(bps[0], bps[-1], rs[0], rs[-1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=+vmax,
    )
    ax.set_xlabel(r"$B_{\mathrm{peak}}$")
    ax.set_ylabel(r"$R$")
    ax.set_title(
        rf"(a) residual $P_{{\mathrm{{sim}}}} - P_{{\mathrm{{ana}}}}$ "
        rf"(max $|\Delta| = {vmax:.1e}$)",
        fontsize=10,
    )
    fig.colorbar(im, ax=ax, shrink=0.85)

    # Panel (b): universal collapse
    ax = axes[1]
    s_dense = np.linspace(0.0, max(s_vals.max() * 1.05, 0.1), 400)
    ax.plot(
        s_dense,
        np.sin(s_dense) ** 2,
        ls="--",
        lw=1.0,
        color="#666",
        label=r"$\sin^2(s)$",
    )
    n_cells = nbp * nr
    ax.plot(
        s_vals,
        p_sim_flat,
        marker="o",
        ms=4,
        lw=0,
        color="#1f77b4",
        label=rf"TIDAL ($N = {n_cells}$ cells)",
    )
    ax.set_xlabel(r"$s = \kappa\,B_{\mathrm{peak}}\,R\,\sqrt{\pi/2}$")
    ax.set_ylabel(r"$P_{g\gamma}$")
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.grid(visible=True, ls=":", alpha=0.3)
    ax.set_title("(b) universal-collapse plot", fontsize=10)

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
