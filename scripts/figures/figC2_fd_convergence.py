"""Figure C.2 — FD-stencil and spectral convergence on the smooth gradient test.

Plots RMS error vs N on log-log axes for FD orders 2, 4, 6 and the FFT-based
pseudo-spectral path, on the test problem f(x) = exp(sin(x)) over a periodic
[0, 2pi] grid. Dashed fitted lines show the measured algebraic rates with
slope alpha-hat and R^2 from log-log linear regression. Spectral curve
descends exponentially to the floating-point floor then saturates.

Data source:  benchmark_results/canonical/fd_convergence.json
Output:       manuscript/figures/figC2_fd_convergence.pdf
Appendix ref: manuscript/sections/appendices/numerical.tex:236
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# Match the manuscript: revtex4-2/PRD, \usepackage{stix}, single-column 3.375 in.
mpl.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{stix}",
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "lines.linewidth": 1.0,
        "axes.linewidth": 0.5,
    }
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "fd_convergence.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figC2_fd_convergence.pdf"

# PRD column width = 3.375 in (85 mm). Height chosen for readability on log-log.
_FIG_WIDTH = 3.375
_FIG_HEIGHT = 2.65

SCHEME_STYLE = {
    "fd_o2": {"label": "FD order 2", "marker": "o", "color": "#1f77b4"},
    "fd_o4": {"label": "FD order 4", "marker": "s", "color": "#ff7f0e"},
    "fd_o6": {"label": "FD order 6", "marker": "^", "color": "#2ca02c"},
    "spectral": {"label": r"Spectral (f.p.\ floor)", "marker": "D", "color": "#d62728"},
}

_FLOOR = 1e-14  # exclude points at/near machine-precision floor


def _regress(
    ns: np.ndarray, errs: np.ndarray
) -> tuple[float | None, float | None, np.ndarray | None]:
    mask = errs > _FLOOR
    if mask.sum() < 2:
        return None, None, None
    slope, intercept, r, _, _ = linregress(np.log10(ns[mask]), np.log10(errs[mask]))
    fitted = 10.0 ** (intercept + slope * np.log10(ns))
    return float(slope), float(r**2), fitted


def _group_by_scheme(rows: list[dict]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    out: dict[str, list[tuple[int, float]]] = {}
    for row in rows:
        out.setdefault(row["scheme"], []).append((row["n"], row["l2_error"]))
    grouped: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for scheme, pairs in out.items():
        pairs.sort()
        ns = np.array([p[0] for p in pairs])
        errs = np.array([p[1] for p in pairs])
        grouped[scheme] = (ns, errs)
    return grouped


def _plot(data: dict, out_path: Path) -> None:
    grouped = _group_by_scheme(data["results"])
    n_vals = sorted({row["n"] for row in data["results"]})

    fig, ax = plt.subplots(figsize=(_FIG_WIDTH, _FIG_HEIGHT))

    # Build labels for FD curves from regression; spectral uses constant label.
    computed_labels: dict[str, str] = {}
    for order, color in [(2, "#1f77b4"), (4, "#ff7f0e"), (6, "#2ca02c")]:
        scheme = f"fd_o{order}"
        if scheme not in grouped:
            continue
        ns, errs = grouped[scheme]
        slope, r_sq, fitted = _regress(ns, errs)
        # Dashed regression line plotted first so markers sit on top.
        if fitted is not None:
            ax.plot(ns, fitted, ls="--", color=color, alpha=0.55, lw=0.8)
        if slope is not None:
            computed_labels[scheme] = (
                rf"FD-{order} ($\hat{{\alpha}}={slope:.2f},\,R^2={r_sq:.4f}$)"
            )
        else:
            computed_labels[scheme] = f"FD-{order}"

    # Data markers on top of regression lines.
    for scheme, (ns, errs) in grouped.items():
        style = SCHEME_STYLE[scheme]
        lbl = computed_labels.get(scheme, style["label"])
        mask = errs > 0
        ax.plot(
            ns[mask],
            errs[mask],
            marker=style["marker"],
            color=style["color"],
            label=lbl,
            lw=1.0,
            ms=4,
        )

    # Floating-point floor reference line.
    fp_eps = 2.0**-52
    ax.axhline(fp_eps, color="gray", ls=":", lw=0.6, alpha=0.6)
    ax.text(
        n_vals[-1] * 1.3,
        fp_eps * 2.5,
        r"$\varepsilon_{\rm mach}$",
        color="gray",
        fontsize=6.5,
        ha="right",
        va="bottom",
    )

    # Annotate the spectral drift tail so it is not mistaken for degradation.
    spec_ns, spec_errs = grouped.get("spectral", (np.array([]), np.array([])))
    drift_mask = spec_errs > fp_eps * 5  # points clearly above the floor
    if drift_mask.sum() >= 2:
        tail_n = spec_ns[drift_mask][-1]
        tail_e = spec_errs[drift_mask][-1]
        ax.annotate(
            r"FFT rounding",
            xy=(tail_n, tail_e),
            xytext=(tail_n * 0.28, tail_e * 8),
            color=SCHEME_STYLE["spectral"]["color"],
            fontsize=6,
            ha="center",
            va="bottom",
            arrowprops={
                "arrowstyle": "-",
                "color": SCHEME_STYLE["spectral"]["color"],
                "lw": 0.5,
                "alpha": 0.7,
            },
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Grid resolution $N$")
    ax.set_ylabel(r"RMS error", labelpad=4)
    pow2_ticks = [n for n in n_vals if (n & (n - 1)) == 0]
    ax.set_xticks(pow2_ticks)
    ax.set_xticklabels([str(n) for n in pow2_ticks])
    ax.set_xlim(n_vals[0] * 0.7, n_vals[-1] * 1.4)
    ax.grid(visible=True, which="both", ls=":", alpha=0.35, lw=0.5)
    ax.legend(
        loc="upper right",
        fontsize=7,
        frameon=False,
        handlelength=1.4,
        handletextpad=0.4,
        borderpad=0.0,
        labelspacing=0.25,
    )
    ax.tick_params(which="both", direction="in", top=False, right=False)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(pad=0.2)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, format="pdf", bbox_inches="tight", pad_inches=0.02)
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
