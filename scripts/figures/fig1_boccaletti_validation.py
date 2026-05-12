r"""Figure 1 (main body) -- Boccaletti validation.

Two-panel main-text figure for §4.1 ValidationResult: (a) the simulated
$P_\\mathrm{final}(B_0)$ curve overlaid with the analytic baseline
$\\sin^2(\\kappa B_0 t/2)$ across the perturbative sweep regime;
(b) the residual $P_\\mathrm{final}^\\mathrm{sim} - P_\\mathrm{final}^\\mathrm{ana}$
on the same x-axis. The detailed multi-resolution convergence panel
lives in App D (figD_boccaletti_calibration); this main-text figure
keeps to the single-resolution validation summary.

Data source:  benchmark_results/canonical/boccaletti_calibration.json
Output:       manuscript/figures/fig1_boccaletti_validation.pdf
Section ref:  manuscript/sections/results.tex (§4.1 BoccalettiValidation)
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_calibration.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "fig1_boccaletti_validation.pdf"


def _plot(data: dict, out_path: Path) -> None:
    meta = data["metadata"]["parameters"]
    kappa = meta["kappa"]
    t_end = meta["t_end"]

    rows = sorted(data["b0_sweep"], key=operator.itemgetter("B0"))
    b0 = np.array([r["B0"] for r in rows])
    p_sim = np.array([r["P_final_sim"] for r in rows])
    residual = np.array([r["residual_final"] for r in rows])

    b0_dense = np.linspace(b0.min(), b0.max(), 400)
    p_analytic = np.sin(0.5 * kappa * b0_dense * t_end) ** 2

    fig, (ax_top, ax_res) = plt.subplots(
        2,
        1,
        figsize=(5.0, 4.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )

    ax_top.plot(
        b0_dense,
        p_analytic,
        ls="--",
        lw=1.0,
        color="#666",
        label=r"$\sin^2(\kappa B_0 t/2)$",
    )
    ax_top.plot(b0, p_sim, marker="o", ms=4, lw=1.0, color="#1f77b4", label="TIDAL")
    ax_top.set_ylabel(r"$P_{g\gamma}$")
    ax_top.legend(frameon=False, fontsize=9, loc="upper right")
    ax_top.grid(visible=True, ls=":", alpha=0.4)

    ax_res.axhline(0.0, color="#999", lw=0.6, ls=":")
    ax_res.plot(b0, residual, marker="o", ms=4, lw=1.0, color="#d62728")
    ax_res.set_xlabel(r"$B_0$")
    ax_res.set_ylabel(r"sim $-$ ana")
    ax_res.grid(visible=True, ls=":", alpha=0.4)
    if residual.size and np.any(residual != 0):
        scale = max(np.max(np.abs(residual)) * 1.3, 1e-8)
        ax_res.set_ylim(-scale, scale)

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
