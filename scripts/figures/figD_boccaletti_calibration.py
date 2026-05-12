r"""Figure D — Boccaletti-kernel calibration for Einstein--Maxwell.

Three-panel figure: (a) the simulated $P_\\mathrm{final}(B_0)$ curve with
the analytic $\\sin^2(\\kappa B_0 t/2)$ overlay; (b) the residual
$P_\\mathrm{final}^\\mathrm{sim} - P_\\mathrm{final}^\\mathrm{analytic}$
vs $B_0$; (c) the multi-resolution convergence of the absolute error
at the fixed regime point ($B_0 = 0.05$, $t_\\mathrm{end} = 50$).

Data source:  benchmark_results/canonical/boccaletti_calibration.json
Output:       manuscript/figures/figD_boccaletti_calibration.pdf
Appendix ref: manuscript/sections/appendices/validation.tex (App D, calibration 1)
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
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_calibration.pdf"


def _plot(data: dict, out_path: Path) -> None:
    meta = data["metadata"]["parameters"]
    kappa = meta["kappa"]
    t_end = meta["t_end"]

    b0_rows = sorted(data["b0_sweep"], key=operator.itemgetter("B0"))
    b0 = np.array([r["B0"] for r in b0_rows])
    p_final_sim = np.array([r["P_final_sim"] for r in b0_rows])
    residual = np.array([r["residual_final"] for r in b0_rows])

    b0_dense = np.linspace(b0.min(), b0.max(), 400)
    p_final_analytic = np.sin(0.5 * kappa * b0_dense * t_end) ** 2

    conv_rows = sorted(data["convergence"], key=operator.itemgetter("N"))
    n_grid = np.array([r["N"] for r in conv_rows])
    abs_err = np.array([r["abs_error_final"] for r in conv_rows])

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))

    ax = axes[0]
    ax.plot(
        b0_dense,
        p_final_analytic,
        ls="--",
        lw=1.0,
        color="#999",
        label=r"$\sin^2(\kappa B_0 t/2)$",
    )
    ax.plot(b0, p_final_sim, marker="o", ms=4, lw=1.0, color="#1f77b4", label="TIDAL")
    ax.set_xlabel(r"$B_0$")
    ax.set_ylabel(r"$P_{\mathrm{final}}$")
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("(a) calibration curve", fontsize=10)

    ax = axes[1]
    ax.axhline(0.0, color="#999", lw=0.6, ls=":")
    ax.plot(b0, residual, marker="o", ms=4, lw=1.0, color="#d62728")
    ax.set_xlabel(r"$B_0$")
    ax.set_ylabel(
        r"$P_{\mathrm{final}}^{\mathrm{sim}} - P_{\mathrm{final}}^{\mathrm{ana}}$"
    )
    ax.set_title("(b) residual vs $B_0$", fontsize=10)
    if residual.size and np.any(residual != 0):
        scale = max(np.max(np.abs(residual)) * 1.3, 1e-8)
        ax.set_ylim(-scale, scale)

    ax = axes[2]
    if n_grid.size:
        ax.loglog(
            n_grid,
            np.maximum(abs_err, 1e-16),
            marker="s",
            ms=5,
            color="#2ca02c",
            lw=1.0,
        )
        # spectral guide (modal solver expected to plateau at machine precision)
        ax.axhline(1e-14, ls=":", lw=0.8, color="#999", label="round-off floor")
        ax.legend(frameon=False, fontsize=9, loc="lower left")
    ax.set_xlabel(r"grid resolution $N$")
    ax.set_ylabel(
        r"$|P_{\mathrm{final}}^{\mathrm{sim}} - P_{\mathrm{final}}^{\mathrm{ana}}|$"
    )
    ax.set_title("(c) multi-resolution convergence", fontsize=10)
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)

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
