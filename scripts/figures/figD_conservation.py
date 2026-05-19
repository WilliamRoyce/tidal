r"""Figure D §8 — energy-conservation audit (App-C styled).

Single-panel log-scale bar chart of $|dE/E|_{\max}$ per example,
**sorted by drift in descending order**, with the default $10^{-3}$
conservation threshold and the IEEE-754 machine-epsilon floor drawn
for reference. Failed examples are dropped (not silently rendered
as empty bars).

Data:   benchmark_results/canonical/conservation_audit_full.json
Output: manuscript/figures/figD_conservation.pdf
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
    REPO_ROOT / "benchmark_results" / "canonical" / "conservation_audit_full.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_conservation.pdf"
EPS_MACH = 2.220446049250313e-16

# Descriptive bar labels (the bar key is the JSON's repository-style
# identifier; the figure must show physics-descriptive names instead).
LABEL_MAP: dict[str, str] = {
    "gertsenshtein": "Einstein–Maxwell Gertsenshtein",
    "coupled_scalars": "Coupled-scalar two-channel mixing",
    "gertsenshtein_proca": "Einstein–Maxwell + Proca-mass photon",
    "coupled_scattering": "Coupled gradient scattering",
    "chern_simons_3d": "Chern–Simons EM (3+1D)",
    "conformal_kg_static": "Conformal Klein–Gordon (static)",
    "cylindrical_kg_1d": "Cylindrical Klein–Gordon (1+1D)",
    "navier_cauchy_2d": "Navier–Cauchy elasticity (2+1D)",
    "gw_plane_wave_1d": "Linearized gravitational wave (1+1D TT)",
    "torsion_gertsenshtein": "Propagating-PGT torsion–Gertsenshtein",
    "torsion_gertsenshtein_nonminimal": "Nonminimal torsion–EM coupling",
    "dark_photon_plasma": "Dark-photon plasma (TorsionCDT)",
    "torsion_dark_photon_fv": "Dark-photon plasma (fundamental-vector)",
}


def _plot(data: dict, out_path: Path) -> None:
    rows = [
        r
        for r in data["results"]
        if r.get("ok") and r.get("max_abs_dE_over_E") is not None
    ]
    # Sort by drift descending
    rows.sort(key=operator.itemgetter("max_abs_dE_over_E"), reverse=True)
    labels = [LABEL_MAP.get(r["label"], r["label"]) for r in rows]
    drifts = np.array([max(r["max_abs_dE_over_E"], EPS_MACH) for r in rows])

    n = len(labels)
    if n == 0:
        msg = "no successful conservation rows"
        raise ValueError(msg)

    fig, ax = plt.subplots(figsize=(max(7.5, 0.55 * n), 4.4))
    x = np.arange(n)
    ax.bar(x, drifts, color="#1f77b4", alpha=0.85)
    ax.set_yscale("log")
    ax.axhline(
        1e-3,
        ls="--",
        lw=1.0,
        color="#444",
        label=r"default conservation threshold $10^{-3}$",
    )
    ax.axhline(
        EPS_MACH,
        ls=":",
        lw=0.7,
        color="#888",
        label=rf"IEEE-754 floor $\varepsilon_{{\mathrm{{mach}}}} = {EPS_MACH:.2e}$",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=9)
    ax.set_ylabel(r"$|dE/E|_{\max}$")
    ax.set_ylim(EPS_MACH / 3, 1.0)
    ax.grid(visible=True, axis="y", which="both", ls=":", alpha=0.3)
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    n_failed = data["summary"].get("n_failed", 0)
    n_total = data["summary"].get("n_examples", n)
    title = (
        rf"energy-conservation audit ($n_{{\mathrm{{ok}}}} = {n}$ of {n_total}, "
        rf"$n_{{\mathrm{{failed}}}} = {n_failed}$ dropped)"
    )
    ax.set_title(title, fontsize=10)

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
