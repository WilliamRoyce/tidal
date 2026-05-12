"""Figure D §7 — three-panel convergence (Proca m², Rabi N, multi-scheme dt).

(figure*) Three panels:
  (a) Proca dispersion omega(k) overlay vs sqrt(k²+m²) at three m².
  (b) Rabi-frequency grid convergence log-log with (k·Δx)² guide.
  (c) Time-integration error vs Δt for seven schemes with expected-order
      guides.

Data:   benchmark_results/canonical/convergence_rich.json
Output: manuscript/figures/figD_convergence.pdf
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "convergence_rich.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_convergence.pdf"
FLOOR = 1e-16


def _plot_dispersion(ax, data: dict) -> None:
    rows = [r for r in data.get("dispersion", []) if r.get("ok")]
    by_m = {}
    for r in rows:
        by_m.setdefault(r["mass2"], []).append((r["k_realised"], r["omega_sim"]))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for i, (m2, pts) in enumerate(sorted(by_m.items())):
        pts.sort()
        ks = np.array([p[0] for p in pts])
        oms = np.array([p[1] for p in pts])
        kd = np.linspace(ks.min() * 0.9, ks.max() * 1.05, 200)
        c = colors[i % len(colors)]
        ax.plot(kd, np.sqrt(kd**2 + m2), ls="--", lw=0.9, color=c, alpha=0.5)
        ax.plot(ks, oms, marker="o", ms=5, lw=0, color=c, label=rf"$m^2={m2}$")
    ax.set_xlabel(r"$k$")
    ax.set_ylabel(r"$\omega$")
    ax.set_title("(a) Proca dispersion", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(visible=True, ls=":", alpha=0.4)


def _plot_rabi(ax, data: dict) -> None:
    rows = [r for r in data.get("rabi_convergence", []) if r.get("ok")]
    if not rows:
        ax.text(
            0.5, 0.5, "no rabi data", ha="center", va="center", transform=ax.transAxes
        )
        return
    rows.sort(key=operator.itemgetter("N"))
    kdx = np.array([r["k_dx"] for r in rows])
    err = np.array([max(abs(r["ratio"] - 1.0), FLOOR) for r in rows])
    ax.loglog(kdx, err, marker="o", ms=5, lw=1.0, color="#1f77b4")
    if kdx.size >= 2:
        guide = err[0] * (kdx / kdx[0]) ** 2
        ax.loglog(
            kdx, guide, ls=":", lw=0.8, color="#888", label=r"$(k\Delta x)^2$ guide"
        )
        ax.legend(frameon=False, fontsize=9)
    ax.set_xlabel(r"$k\,\Delta x$")
    ax.set_ylabel(r"$|\Omega_{\rm eff}/\Omega_{\rm theory} - 1|$")
    ax.set_title("(b) Rabi-frequency grid convergence", fontsize=10)
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)


def _plot_tio(ax, data: dict) -> None:
    rows = [r for r in data.get("time_integration_order", []) if r.get("ok")]
    slopes = data["summary"].get("tio_fitted_slopes", {})
    by_scheme = {}
    for r in rows:
        by_scheme.setdefault(r["scheme"], []).append(r)
    colors = plt.get_cmap("tab10")
    for i, (scheme, srows) in enumerate(sorted(by_scheme.items())):
        srows.sort(key=operator.itemgetter("dt"))
        dts = np.array([r["dt"] for r in srows])
        errs = np.array([max(r["abs_error"], FLOOR) for r in srows])
        slope = slopes.get(scheme)
        label = scheme + (f" (slope={slope:.2f})" if slope is not None else "")
        ax.loglog(dts, errs, marker="o", ms=4, lw=1.0, color=colors(i), label=label)
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"$|P_{\rm final}^{\rm sim} - P_{\rm bare}|$")
    ax.set_title("(c) time-integration order", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best")
    ax.grid(visible=True, which="both", ls=":", alpha=0.4)


def _plot(data: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 3.8))
    _plot_dispersion(axes[0], data)
    _plot_rabi(axes[1], data)
    _plot_tio(axes[2], data)
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
