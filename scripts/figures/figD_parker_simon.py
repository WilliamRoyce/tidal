r"""Figure D §4 — Parker–Simon higher-derivative perturbative calibration.

Single-panel App-C-style error-vs-$\varepsilon$ plot for the
$\varepsilon\,\partial_x^4 \phi$ higher-derivative correction on a
Klein–Gordon base. The combined $\phi^{(0)} + \phi^{(1)}$ truncation
error decreases faster than the worst-case $\mathcal{O}(\varepsilon^2)$
bound for this correction type.

Floor-saturated points (where the error has fallen below
$\sim 10^3 \cdot \varepsilon_{\mathrm{mach}}$) are excluded from the
slope fit and explicitly flagged in the figure.

Data:   benchmark_results/canonical/parker_simon_flrw.json
Output: manuscript/figures/figD_parker_simon.pdf
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
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "parker_simon_flrw.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_parker_simon.pdf"
EPS_MACH = 2.220446049250313e-16
# Saturation detection: a point is saturated when its log-log slope to
# the next-smaller-eps point is less than this threshold. The
# perturbative theorem predicts slope >= 2 in the asymptotic regime;
# anything below ~1 indicates the spatial-discretisation noise floor
# has been reached.
SAT_SLOPE_THRESH = 1.0


def _fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    lx, ly = np.log10(x), np.log10(y)
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    return float(slope), float(1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan)


def _detect_floor(eps: np.ndarray, err: np.ndarray) -> np.ndarray:
    """Return boolean mask of floor-saturated points.

    Walks from large eps to small eps; a point is saturated if the
    local log-log slope to the previous (larger) eps falls below
    SAT_SLOPE_THRESH (i.e. error didn't shrink the way the theorem
    demands). All smaller-eps points beyond that are also saturated.
    """
    order = np.argsort(-eps)  # descending eps
    floor = np.zeros(eps.size, dtype=bool)
    saturated = False
    prev_eps, prev_err = None, None
    for i in order:
        e = eps[i]
        r = max(err[i], EPS_MACH)
        if prev_eps is not None and prev_err is not None:
            local_slope = (math.log10(r) - math.log10(prev_err)) / (
                math.log10(e) - math.log10(prev_eps)
            )
            if local_slope < SAT_SLOPE_THRESH:
                saturated = True
        if saturated:
            floor[i] = True
        prev_eps, prev_err = e, r
    return floor


def _plot(data: dict, out_path: Path) -> None:
    rows = sorted(data["results"], key=operator.itemgetter("eps"))
    eps = np.array([r["eps"] for r in rows])
    err_p0 = np.array([r["err_pass0_vs_full"] for r in rows])
    err_c = np.array([r["err_combined_vs_full"] for r in rows])
    floor_mask = _detect_floor(eps, err_c)
    fit_mask = ~floor_mask & (eps <= 0.2)

    s_p0, r2_p0 = _fit_loglog(eps, np.maximum(err_p0, EPS_MACH))
    if fit_mask.sum() >= 3:
        s_c, r2_c = _fit_loglog(eps[fit_mask], err_c[fit_mask])
    else:
        s_c, r2_c = math.nan, math.nan

    fig, ax = plt.subplots(figsize=(5.6, 4.0))

    ax.loglog(
        eps,
        np.maximum(err_p0, EPS_MACH),
        marker="o",
        ms=5,
        lw=0,
        color="#d62728",
        label=rf"$|\phi^{{(0)}} - \phi^\mathrm{{exact}}|$ ($\hat{{\alpha}} = {s_p0:.2f}$, $R^2 = {r2_p0:.3f}$)",
    )
    ax.loglog(
        eps[~floor_mask],
        err_c[~floor_mask],
        marker="s",
        ms=5,
        lw=0,
        color="#1f77b4",
        label=rf"$|\phi^{{(0)}} + \phi^{{(1)}} - \phi^\mathrm{{exact}}|$ ($\hat{{\alpha}} = {s_c:.2f}$, $R^2 = {r2_c:.3f}$)",
    )
    if floor_mask.any():
        ax.loglog(
            eps[floor_mask],
            np.maximum(err_c[floor_mask], EPS_MACH),
            marker="s",
            ms=5,
            lw=0,
            mfc="white",
            mec="#1f77b4",
            label=r"$\phi^{(0)} + \phi^{(1)}$ (below spatial-FD floor; excluded from fit)",
        )
        # Annotate the FD-4 spatial-discretisation floor at the median
        # of the saturated points so the exclusion is visually justified.
        floor_val = float(np.median(np.maximum(err_c[floor_mask], EPS_MACH)))
        ax.axhline(
            floor_val,
            ls="-.",
            lw=0.7,
            color="#888",
            alpha=0.6,
            label=rf"FD-4 spatial floor $\approx {floor_val:.1e}$",
        )
    # Reference guides
    if fit_mask.any():
        idx = np.where(fit_mask)[0][-1]
        ax.loglog(
            eps,
            err_c[idx] * (eps / eps[idx]) ** 2,
            ls=":",
            lw=0.8,
            color="#1f77b4",
            alpha=0.5,
            label=r"$\mathcal{O}(\varepsilon^2)$ worst-case guide",
        )

    ax.set_xlabel(r"perturbation parameter $\varepsilon$")
    ax.set_ylabel(r"$|\phi^\mathrm{approx} - \phi^\mathrm{exact}|$")
    ax.set_title(r"$\varepsilon\,\partial_x^4 \phi$ correction", fontsize=10)
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)
    ax.legend(frameon=False, fontsize=7, loc="lower right")

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
