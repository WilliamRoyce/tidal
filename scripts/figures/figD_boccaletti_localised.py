r"""Figure D §2 — path-integrated Boccaletti reproduction (App-C styled).

Single-panel universal-collapse plot: every simulation cell of the
$(B_\mathrm{peak}, R)$ sweep is plotted against the dimensionless
variable $s \equiv \kappa\,B_\mathrm{peak}\,R\,\sqrt{\pi/2}$,
overlaid with the analytic prediction $\sin^2(s)$.

This is the strongest visual demonstration that the path-integrated
kernel is correct: every simulation point should fall on the
universal curve regardless of its individual $(B_\mathrm{peak}, R)$
values. Points at $s > \pi/2$ (past the first quarter-period) lie
below $\sin^2(s)$ as the bare-formula coherent approximation breaks
down.

Data:   benchmark_results/canonical/boccaletti_localised.json
Output: manuscript/figures/figD_boccaletti_localised.pdf
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from _palette import IBM_PALETTE

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_localised.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_localised.pdf"
KAPPA = 1.0


def _plot(data: dict, out_path: Path) -> None:
    rows = [r for r in data["results"] if r.get("ok")]
    s_vals = np.array(
        [KAPPA * r["Bpeak"] * r["R"] * math.sqrt(math.pi / 2.0) for r in rows]
    )
    p_sim_flat = np.array([r["P_sim"] for r in rows])
    len(rows)

    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    s_dense = np.linspace(0.0, max(s_vals.max() * 1.05, 0.1), 400)
    ax.plot(
        s_dense,
        np.sin(s_dense) ** 2,
        ls="--",
        lw=1.0,
        color="#666",
        label=r"$\sin^2(s)$",
    )
    ax.plot(
        s_vals,
        p_sim_flat,
        marker="o",
        ms=4,
        lw=0,
        color=IBM_PALETTE["blue"],
        label=r"$\mathit{TIDAL}$",
    )
    ax.set_xlabel(r"$s = \kappa\,B_{\mathrm{peak}}\,R\,\sqrt{\pi/2}$")
    ax.set_ylabel(r"$P_{g\gamma}$")
    ax.legend(frameon=False, fontsize=9, loc="best")
    ax.grid(visible=True, ls=":", alpha=0.3)

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
