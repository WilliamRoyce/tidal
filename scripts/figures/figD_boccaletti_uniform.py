"""Figure D §1 — Boccaletti uniform-B0 three-overlay calibration.

Three-panel figure (figure* width, two columns):
  (a) P_final(B0) vs three analytic predictions (bare, energy-cap,
      Raffelt-Stodolsky) at one t_end.
  (b) Residual band: P_final_sim − bare, P_max_sim − capped,
      P_max_sim − RS, on the same axis.
  (c) Multi-N convergence of |P_final_sim − bare| at the regime point.

Data:   benchmark_results/canonical/boccaletti_uniform_rich.json
Output: manuscript/figures/figD_boccaletti_uniform.pdf
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
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_uniform.pdf"
FLOOR = 1e-16


def _plot(data: dict, out_path: Path) -> None:
    grid = data["grid"]
    if not grid:
        msg = "no grid rows"
        raise ValueError(msg)
    t_values = sorted({r["t_end"] for r in grid})
    t0 = t_values[0]
    rows0 = sorted([r for r in grid if r["t_end"] == t0], key=operator.itemgetter("B0"))
    b0 = np.array([r["B0"] for r in rows0])
    p_final = np.array([r["P_final_sim"] for r in rows0])
    p_max = np.array([r["P_max_sim"] for r in rows0])
    p_bare = np.array([r["P_bare"] for r in rows0])
    p_capped = np.array([r["P_energy_cap"] for r in rows0])
    p_rs = np.array([r["P_raffelt_stodolsky"] for r in rows0])

    conv = sorted(data.get("convergence", []), key=operator.itemgetter("N"))
    n_grid = np.array([r["N"] for r in conv]) if conv else np.array([])
    abs_err = np.array([r["abs_error_final"] for r in conv]) if conv else np.array([])

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.6))

    # (a) calibration overlay
    ax = axes[0]
    b0d = np.linspace(b0.min(), b0.max(), 400)
    kappa = data["metadata"]["parameters"]["kappa"]
    ax.plot(
        b0d,
        np.sin(0.5 * kappa * b0d * t0) ** 2,
        ls=":",
        lw=1.0,
        color="#444",
        label="bare sin²",
    )
    ax.plot(b0, p_capped, ls="--", lw=0.9, color="#999", label="energy-cap")
    ax.plot(b0, p_rs, ls="-.", lw=0.9, color="#aa6", label="Raffelt–Stodolsky")
    ax.plot(
        b0,
        p_final,
        marker="o",
        ms=4,
        lw=0,
        color="#1f77b4",
        label=r"$P^{\rm sim}_{\rm final}$",
    )
    ax.plot(
        b0,
        p_max,
        marker="s",
        ms=4,
        lw=0,
        color="#d62728",
        label=r"$P^{\rm sim}_{\rm max}$",
    )
    ax.set_xlabel(r"$B_0$")
    ax.set_ylabel(r"$P$")
    ax.set_title(f"(a) calibration at $t={t0}$", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(visible=True, ls=":", alpha=0.4)

    # (b) residual band
    ax = axes[1]
    ax.axhline(0.0, color="#999", lw=0.6, ls=":")
    ax.plot(
        b0, p_final - p_bare, marker="o", ms=4, color="#1f77b4", label="final − bare"
    )
    ax.plot(b0, p_max - p_capped, marker="s", ms=4, color="#d62728", label="max − cap")
    ax.plot(b0, p_max - p_rs, marker="^", ms=4, color="#aa6", label="max − RS")
    ax.set_xlabel(r"$B_0$")
    ax.set_ylabel("sim − analytic")
    ax.set_title("(b) residual band", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(visible=True, ls=":", alpha=0.4)

    # (c) convergence
    ax = axes[2]
    if n_grid.size:
        ax.loglog(
            n_grid,
            np.maximum(abs_err, FLOOR),
            marker="o",
            ms=5,
            color="#2ca02c",
            lw=1.0,
        )
        ax.axhline(1e-14, ls=":", lw=0.8, color="#999", label="round-off floor")
        ax.legend(frameon=False, fontsize=8)
    ax.set_xlabel(r"$N$")
    ax.set_ylabel(r"$|P^{\rm sim}_{\rm final} - P_{\rm bare}|$")
    ax.set_title("(c) N-convergence", fontsize=10)
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
