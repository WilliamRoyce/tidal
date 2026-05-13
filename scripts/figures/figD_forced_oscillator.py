r"""Figure D §3 — forced-oscillator Pass 0+1 calibration (App-C styled).

Single-panel error-vs-$\varepsilon$ plot calibrating the closed-form
Duhamel kernel against the full-$\varepsilon$ modal solution on a driven
Klein–Gordon oscillator with an $\varepsilon$ mass-shift correction.

Two curves: Pass 0 only (the bare modal solution at $\varepsilon = 0$,
showing $\mathcal{O}(\varepsilon)$ error against the full-$\varepsilon$
ground truth) and Pass 0 + Pass 1 (the closed-form Duhamel correction,
showing $\mathcal{O}(\varepsilon^2)$ error). Fitted log–log slopes and
$R^2$ values are reported in the legend.

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
        lw=1.0,
        color="#d62728",
        label=rf"Pass 0 only ($\hat{{\alpha}} = {s_p0:.2f}$, $R^2 = {r2_p0:.3f}$)",
    )
    ax.loglog(
        eps,
        np.maximum(err_c, EPS_MACH),
        marker="s",
        ms=5,
        lw=1.0,
        color="#1f77b4",
        label=rf"Pass 0 + Pass 1 ($\hat{{\alpha}} = {s_c:.2f}$, $R^2 = {r2_c:.3f}$)",
    )
    # O(eps^2) reference guide anchored at largest eps
    ax.loglog(
        eps,
        err_c[-1] * (eps / eps[-1]) ** 2,
        ls=":",
        lw=0.7,
        color="#1f77b4",
        alpha=0.5,
        label=r"$\mathcal{O}(\varepsilon^2)$ guide",
    )
    # O(eps) reference for Pass 0
    ax.loglog(
        eps,
        err_p0[-1] * (eps / eps[-1]),
        ls=":",
        lw=0.7,
        color="#d62728",
        alpha=0.5,
        label=r"$\mathcal{O}(\varepsilon)$ guide",
    )
    ax.axhline(
        EPS_MACH, ls=":", lw=0.7, color="#666", label=r"$\varepsilon_{\mathrm{mach}}$"
    )

    ax.set_xlabel(r"perturbation parameter $\varepsilon$")
    ax.set_ylabel(r"error vs full-$\varepsilon$ modal solution")
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
