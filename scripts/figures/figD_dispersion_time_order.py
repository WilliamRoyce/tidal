r"""Figure D -- Dispersion + time-integration order calibration.

Two-panel figure: (a) Proca dispersion $\\omega(k)$ from FFT of TIDAL
simulations overlaid with the analytic baseline $\\omega = \\sqrt{k^2 + m^2}$;
(b) absolute $P_\\mathrm{final}$ error vs $\\Delta t$ for leapfrog
Yoshida-2, leapfrog Yoshida-4 and CVODE, with log--log slopes annotated
against the expected algebraic orders.

Data sources:
  benchmark_results/canonical/proca_dispersion.json
  benchmark_results/canonical/time_integration_order.json
Output:       manuscript/figures/figD_dispersion_time_order.pdf
Appendix ref: manuscript/sections/appendices/validation.tex (App D, calibration 4)
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DISP = REPO_ROOT / "benchmark_results" / "canonical" / "proca_dispersion.json"
DEFAULT_TIO = (
    REPO_ROOT / "benchmark_results" / "canonical" / "time_integration_order.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_dispersion_time_order.pdf"

FLOOR = 1e-16


def _plot_dispersion(ax: plt.Axes, data: dict) -> None:
    rows = [r for r in data.get("k_sweep", []) if r.get("ok")]
    if not rows:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        return
    mass2 = data["metadata"]["parameters"]["mass2"]
    ks = np.array([r["k_realised"] for r in rows])
    omegas = np.array([r["omega_sim"] for r in rows])
    k_dense = np.linspace(ks.min() * 0.9, ks.max() * 1.05, 200)
    ax.plot(
        k_dense,
        np.sqrt(k_dense**2 + mass2),
        ls="--",
        lw=1.0,
        color="#999",
        label=rf"$\sqrt{{k^2 + m^2}}$ ($m^2={mass2}$)",
    )
    ax.plot(ks, omegas, marker="o", ms=5, lw=0, color="#1f77b4", label="TIDAL")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$\omega$")
    ax.set_title("(a) Proca dispersion", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.grid(visible=True, ls=":", alpha=0.4)


def _plot_order(ax: plt.Axes, data: dict) -> None:
    runs = [r for r in data.get("runs", []) if r.get("ok")]
    if not runs:
        ax.text(0.5, 0.5, "no data", ha="center", va="center", transform=ax.transAxes)
        return
    groups = {
        ("leapfrog", 2): {"marker": "o", "color": "#1f77b4", "label": "leapfrog Y2"},
        ("leapfrog", 4): {"marker": "s", "color": "#ff7f0e", "label": "leapfrog Y4"},
        ("cvode", None): {"marker": "^", "color": "#2ca02c", "label": "cvode"},
    }
    slopes = data.get("summary", {}).get("fitted_slopes", {})
    for (scheme, lf), style in groups.items():
        rows = sorted(
            [
                r
                for r in runs
                if r["scheme"] == scheme and r.get("leapfrog_order") == lf
            ],
            key=operator.itemgetter("dt"),
        )
        if not rows:
            continue
        dts = np.array([r["dt"] for r in rows])
        errs = np.array([max(r["abs_error"], FLOOR) for r in rows])
        ax.loglog(
            dts,
            errs,
            marker=style["marker"],
            color=style["color"],
            lw=1.0,
            ms=5,
            label=style["label"],
        )

    # expected-rate guides anchored at the largest dt
    for order, color in [(2, "#1f77b4"), (4, "#ff7f0e")]:
        rows = sorted(
            [
                r
                for r in runs
                if r["scheme"] == "leapfrog" and r.get("leapfrog_order") == order
            ],
            key=operator.itemgetter("dt"),
        )
        if not rows:
            continue
        anchor = rows[-1]
        guide_dts = np.array([min(r["dt"] for r in rows), max(r["dt"] for r in rows)])
        guide = anchor["abs_error"] * (guide_dts / anchor["dt"]) ** order
        ax.loglog(guide_dts, guide, ls="--", color=color, alpha=0.3, lw=0.9)

    # annotate fitted slopes
    txts = []
    for key in ("leapfrog_2", "leapfrog_4"):
        s = slopes.get(key)
        if s is not None:
            txts.append(f"slope({key}) = {s:.2f}")
    if txts:
        ax.text(
            0.05,
            0.05,
            "\n".join(txts),
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
            ha="left",
            bbox={"facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
        )

    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"$|P_{\rm final}^{\rm sim} - P_{\rm final}^{\rm ana}|$")
    ax.set_title("(b) time-integration order", fontsize=10)
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)


def _plot(disp: dict, tio: dict, out_path: Path) -> None:
    fig, (ax_d, ax_t) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    _plot_dispersion(ax_d, disp)
    _plot_order(ax_t, tio)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dispersion", type=Path, default=DEFAULT_DISP)
    parser.add_argument("--time-order", type=Path, default=DEFAULT_TIO)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    with args.dispersion.open() as fh:
        disp = json.load(fh)
    with args.time_order.open() as fh:
        tio = json.load(fh)
    _plot(disp, tio, args.out)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
