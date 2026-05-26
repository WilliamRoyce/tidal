r"""Figure D §7 — spatial and temporal convergence (App-C styled).

Two-panel `figure*`:

  (a) Proca dispersion: $\omega(k)$ extracted by FFT from TIDAL at
      three values of $m^2 \in \{0.1, 1.0, 10.0\}$, overlaid with the
      analytic $\omega(k) = \sqrt{k^2 + m^2}$ at each $m^2$.
  (b) Time-integration order: $|P_\mathrm{final}^\mathrm{sim} -
      \sin^2(\kappa B_0 t/2)|$ versus $\Delta t$ for leapfrog
      Yoshida-2 and Yoshida-4. Fitted log–log slopes and $R^2$
      annotated in the legend. Implicit schemes are deliberately
      excluded (their flat-at-rtol behaviour is covered in App D §6).

The Rabi-frequency grid-convergence panel that previously sat in
position (b) has been dropped: the data is FFT-bin-width bounded
across the surveyed N range and adds no spatial-convergence
information at this resolution.

Data:   benchmark_results/canonical/convergence_rich.json
Output: manuscript/figures/figD_convergence.pdf
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
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "convergence_rich.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_convergence.pdf"
EPS_MACH = 2.220446049250313e-16
FLOOR_TOL = 1e3 * EPS_MACH


def _fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    lx, ly = np.log10(x), np.log10(y)
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    return float(slope), float(1.0 - ss_res / ss_tot if ss_tot > 0 else math.nan)


def _plot_dispersion(ax, data: dict) -> None:
    rows = [r for r in data.get("dispersion", []) if r.get("ok")]
    by_m: dict[float, list[tuple[float, float]]] = {}
    for r in rows:
        by_m.setdefault(r["mass2"], []).append((r["k_realized"], r["omega_sim"]))
    colors = [IBM_PALETTE["blue"], IBM_PALETTE["orange"], IBM_PALETTE["purple"]]
    for i, (m2, pts) in enumerate(sorted(by_m.items())):
        pts.sort()
        ks = np.array([p[0] for p in pts])
        oms = np.array([p[1] for p in pts])
        kd = np.linspace(ks.min() * 0.9, ks.max() * 1.05, 200)
        c = colors[i % len(colors)]
        ax.plot(kd, np.sqrt(kd**2 + m2), ls="--", lw=0.9, color=c, alpha=0.5)
        ax.plot(ks, oms, marker="o", ms=5, lw=0, color=c, label=rf"$m^2 = {m2}$")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$\omega$")
    ax.set_title("(a) Proca dispersion", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.grid(visible=True, ls=":", alpha=0.3)


def _plot_rabi(ax, data: dict) -> None:
    rows = sorted(
        [r for r in data.get("rabi_convergence", []) if r.get("ok")],
        key=operator.itemgetter("N"),
    )
    if not rows:
        ax.text(
            0.5, 0.5, "no rabi data", ha="center", va="center", transform=ax.transAxes
        )
        return
    kdx = np.array([r["k_dx"] for r in rows])
    err = np.array([max(abs(r["ratio"] - 1.0), EPS_MACH) for r in rows])
    ax.loglog(
        kdx,
        err,
        marker="o",
        ms=5,
        lw=1.0,
        color=IBM_PALETTE["blue"],
        label=r"$|\Omega_{\mathrm{eff}}/\Omega_{\mathrm{theory}} - 1|$",
    )
    # If the data is flat (FFT-extraction floor reached at every N),
    # annotate the plateau instead of fitting a meaningless slope.
    e_lo, e_hi = float(err.min()), float(err.max())
    is_flat = e_hi / max(e_lo, EPS_MACH) < 2.0
    if is_flat:
        plateau = float(np.median(err))
        ax.axhline(
            plateau,
            ls="-.",
            lw=0.9,
            color=IBM_PALETTE["blue"],
            alpha=0.5,
            label=rf"FFT-extraction floor $\approx {plateau:.1e}$",
        )
        # Anchor the y-range so the plateau reads as a clean
        # horizontal at the FFT floor rather than zoomed-in noise.
        ax.set_ylim(plateau / 30.0, plateau * 30.0)
    else:
        anchor_x = float(kdx.max())
        anchor_y = float(err[int(np.argmax(kdx))])
        ax.loglog(
            kdx,
            anchor_y * (kdx / anchor_x) ** 2,
            ls=":",
            lw=0.8,
            color="#888",
            label=r"$(k\Delta x)^2$ guide",
        )
    ax.set_xlabel(r"$k\,\Delta x$")
    ax.set_ylabel(r"$|\Omega_{\mathrm{eff}}/\Omega_{\mathrm{theory}} - 1|$")
    ax.set_title("(b) Rabi-frequency grid convergence", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)


def _plot_tio(ax, data: dict) -> None:
    # Prefer modal-reference error; spatial/IC/extraction floors cancel.
    err_key = "abs_error_vs_modal"
    rows = [
        r
        for r in data.get("time_integration_order", [])
        if r.get("ok") and r.get(err_key) is not None and math.isfinite(r[err_key])
    ]
    if not rows:
        err_key = "abs_error"
        rows = [
            r
            for r in data.get("time_integration_order", [])
            if r.get("ok") and r.get(err_key) is not None and math.isfinite(r[err_key])
        ]
    by_scheme: dict[str, list[dict]] = {}
    for r in rows:
        by_scheme.setdefault(r["scheme"], []).append(r)
    colors = {"leapfrog_Y2": IBM_PALETTE["blue"], "leapfrog_Y4": IBM_PALETTE["orange"]}
    markers = {"leapfrog_Y2": "o", "leapfrog_Y4": "s"}
    for scheme, srows in sorted(by_scheme.items()):
        srows.sort(key=operator.itemgetter("dt"))
        dts = np.array([r["dt"] for r in srows])
        errs = np.array([max(r[err_key], EPS_MACH) for r in srows])
        fit_mask = errs > FLOOR_TOL
        slope: float | None = None
        intercept: float | None = None
        if fit_mask.sum() >= 3:
            slope, r2 = _fit_loglog(dts[fit_mask], errs[fit_mask])
            # refit to capture intercept (np.polyfit returns slope only via _fit_loglog)
            lx = np.log10(dts[fit_mask])
            ly = np.log10(errs[fit_mask])
            _, intercept = np.polyfit(lx, ly, 1)
            label = rf"{scheme} ($\hat{{\alpha}} = {slope:.2f}$, $R^2 = {r2:.3f}$)"
        else:
            label = scheme
        ax.loglog(
            dts,
            errs,
            marker=markers.get(scheme, "."),
            ms=5,
            lw=0,
            color=colors.get(scheme, "#666"),
            label=label,
        )
        # Overlay the least-squares regression line as a dashed guide so the
        # reader can see at a glance what slope the legend reports. For Y4 this
        # also makes the size of the residual scatter visible against the fit.
        if slope is not None and intercept is not None:
            dt_line = np.array([dts.min(), dts.max()])
            err_line = 10.0 ** (slope * np.log10(dt_line) + intercept)
            ax.loglog(
                dt_line,
                err_line,
                ls="--",
                lw=1.0,
                color=colors.get(scheme, "#666"),
                alpha=0.6,
            )
    ax.set_xlabel(r"$\Delta t$")
    if err_key == "abs_error_vs_modal":
        ax.set_ylabel(
            r"$|P_{\mathrm{final}}^{\mathrm{sim}} - P_{\mathrm{final}}^{\mathrm{modal}}|$"
        )
    else:
        ax.set_ylabel(r"$|P_{\mathrm{final}}^{\mathrm{sim}} - \sin^2(\kappa B_0 t/2)|$")
    ax.set_title("(b) time-integration order (leapfrog only)", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)


def _plot(data: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    _plot_dispersion(axes[0], data)
    _plot_tio(axes[1], data)
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
