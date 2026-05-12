r"""Figure 1 (main body) -- Boccaletti validation.

Two-panel main-text figure for §4.1 ValidationResult: (a) the simulated
$P_{\rm final}(B_0)$ curve overlaid with the analytic baseline
$\sin^2(\kappa B_0 t/2)$ at the fixed validation $t_{\rm end}$;
(b) the residual $P_{\rm final}^{\rm sim} - P_{\rm final}^{\rm bare}$
on the same x-axis.

The detailed multi-overlay calibration (energy-cap and
Raffelt--Stodolsky comparisons) and the multi-N convergence panel live
in App D (figD_boccaletti_uniform); this main-text figure keeps to the
bare-formula validation summary.

Data source:  benchmark_results/canonical/boccaletti_uniform_rich.json
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
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_uniform_rich.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "fig1_boccaletti_validation.pdf"


def _plot(data: dict, out_path: Path) -> None:
    grid = data.get("grid") or data.get("results") or []
    if not grid:
        msg = "no grid rows"
        raise ValueError(msg)
    # Pick the first t_end as the validation slice.
    t_values = sorted({r["t_end"] for r in grid})
    t0 = t_values[0]
    rows = sorted([r for r in grid if r["t_end"] == t0], key=operator.itemgetter("B0"))
    b0 = np.array([r["B0"] for r in rows])
    p_sim = np.array([r["P_final_sim"] for r in rows])
    p_bare = np.array([r["P_bare"] for r in rows])

    kappa = data["metadata"]["parameters"]["kappa"]

    fig, (ax_top, ax_res) = plt.subplots(
        2,
        1,
        figsize=(5.0, 4.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )
    b0d = np.linspace(b0.min(), b0.max(), 400)
    ax_top.plot(
        b0d,
        np.sin(0.5 * kappa * b0d * t0) ** 2,
        ls="--",
        lw=1.0,
        color="#666",
        label=r"$\sin^2(\kappa B_0 t/2)$",
    )
    ax_top.plot(b0, p_sim, marker="o", ms=4, lw=1.0, color="#1f77b4", label="TIDAL")
    ax_top.set_ylabel(r"$P_{g\gamma}$")
    ax_top.legend(frameon=False, fontsize=9, loc="upper right")
    ax_top.grid(visible=True, ls=":", alpha=0.4)

    residual = p_sim - p_bare
    ax_res.axhline(0.0, color="#999", lw=0.6, ls=":")
    ax_res.plot(b0, residual, marker="o", ms=4, lw=1.0, color="#d62728")
    ax_res.set_xlabel(r"$B_0$")
    ax_res.set_ylabel("sim − ana")
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
