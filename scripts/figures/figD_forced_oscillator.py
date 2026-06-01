r"""Figure D §3 — forced-oscillator order-eps^{0+1} calibration (App-C styled).

Single-panel error-vs-$\varepsilon$ plot calibrating the closed-form
Duhamel kernel against the full-$\varepsilon$ modal solution on a driven
Klein–Gordon oscillator with an $\varepsilon$ mass-shift correction.

Two curves: $|\phi^{(0)} - \phi^\mathrm{exact}|$ (base modal solution at
$\varepsilon = 0$, showing $\mathcal{O}(\varepsilon)$ error against the
full-$\varepsilon$ ground truth) and
$|\phi^{(0)} + \phi^{(1)} - \phi^\mathrm{exact}|$ (closed-form Duhamel
correction, showing $\mathcal{O}(\varepsilon^2)$ error). Fitted log–log
slopes and $R^2$ values are reported in the legend.

Data:   benchmark_results/canonical/forced_oscillator.json
Output: manuscript/figures/figD_forced_oscillator.pdf
"""

from __future__ import annotations

import argparse
import json
import math
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _palette import IBM_PALETTE

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "forced_oscillator.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_forced_oscillator.pdf"
EPS_MACH = 2.220446049250313e-16
FLOOR_TOL = 1e3 * EPS_MACH  # consider points within 1000*eps of floor as saturated


def _fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    lx, ly = np.log10(x), np.log10(y)
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    return float(slope), float(1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan)


def _plot(data: dict, out_path: Path) -> None:
    rows = sorted(data["results"], key=operator.itemgetter("eps"))
    eps = np.array([r["eps"] for r in rows])
    err_p0 = np.array([r["err_pass0_vs_full"] for r in rows])
    err_c = np.array([r["err_combined_vs_full"] for r in rows])

    # Filter floor-saturated points from the combined fit
    mask_c = (err_c > FLOOR_TOL) & (eps <= 0.2)
    s_p0, r2_p0 = _fit_loglog(eps, np.maximum(err_p0, EPS_MACH))
    s_c, r2_c = (
        _fit_loglog(eps[mask_c], err_c[mask_c])
        if mask_c.sum() >= 3
        else (math.nan, math.nan)
    )

    fig, ax = plt.subplots(figsize=(5.6, 4.0))

    ax.loglog(
        eps,
        np.maximum(err_p0, EPS_MACH),
        marker="o",
        ms=5,
        lw=0,
        color=IBM_PALETTE["magenta"],
        label=rf"$|\phi^{{(0)}} - \phi^\mathrm{{exact}}|$ ($\hat{{\alpha}} = {s_p0:.2f}$, $R^2 = {r2_p0:.3f}$)",
    )
    ax.loglog(
        eps,
        np.maximum(err_c, EPS_MACH),
        marker="s",
        ms=5,
        lw=0,
        color=IBM_PALETTE["blue"],
        label=rf"$|\phi^{{(0)}} + \phi^{{(1)}} - \phi^\mathrm{{exact}}|$ ($\hat{{\alpha}} = {s_c:.2f}$, $R^2 = {r2_c:.3f}$)",
    )
    # O(eps^2) reference guide anchored at largest eps (distinct gray
    # so the line doesn't merge with the phi^(0)+phi^(1) trace).
    ax.loglog(
        eps,
        err_c[-1] * (eps / eps[-1]) ** 2,
        ls="--",
        lw=0.8,
        color="#444",
        alpha=0.7,
        label=r"$\mathcal{O}(\varepsilon^2)$ guide",
    )
    # O(eps) reference for the base-solution error |phi^(0) - phi^exact|
    ax.loglog(
        eps,
        err_p0[-1] * (eps / eps[-1]),
        ls="--",
        lw=0.8,
        color="#888",
        alpha=0.7,
        label=r"$\mathcal{O}(\varepsilon)$ guide",
    )

    ax.set_xlabel(r"perturbation parameter $\varepsilon$")
    ax.set_ylabel(r"$|\phi^\mathrm{approx} - \phi^\mathrm{exact}|$")
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)
    ax.legend(frameon=False, fontsize=8, loc="best")

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
