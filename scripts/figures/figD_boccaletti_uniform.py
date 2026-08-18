r"""Figure D §1 — Boccaletti uniform-B0 calibration (App-C styled).

Two-panel `figure*` (two-column spread) showing the photon–graviton
conversion probability $P_\mathrm{final}(B_0)$ on a uniform background
calibrated against the analytic Boccaletti kernel
$\sin^2(\kappa B_0 t/2)$ and the Raffelt–Stodolsky formula
with the graviton effective-mass detuning.

  (a) $P_\mathrm{final}^{\mathrm{sim}}(B_0)$ across the full sweep at
      fixed $t=t_\mathrm{end}$ with smooth analytic overlays at 400
      dense $B_0$ values. The bare sin² (dotted) and the Raffelt–
      Stodolsky two-mode formula (dashed) coincide at small
      $\kappa B_0 t/2$; the latter is shown for comparison.
  (b) Multi-resolution convergence
      $|P_\mathrm{final}^{\mathrm{sim}} - \sin^2(\kappa B_0 t/2)|$
      versus grid resolution $N$ on log–log axes, with the IEEE-754
      machine-epsilon floor drawn as a horizontal reference.

The figure makes a single coherent claim: the Boccaletti kernel is
reproduced uniformly across the perturbative regime and the residual
falls toward the round-off floor under grid refinement. No `P_max`
or energy-cap form is plotted — `P_max` saturates at the cap once a
Rabi peak fires and is not a clean perturbative-regime metric.

Data:   benchmark_results/canonical/boccaletti_uniform_rich.json
Output: manuscript/figures/figD_boccaletti_uniform.pdf
"""

from __future__ import annotations

import argparse
import json
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
from _palette import IBM_PALETTE

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_uniform_rich.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "figures" / "figD_boccaletti_uniform.pdf"
EPS_MACH = 2.220446049250313e-16


def _raffelt_stodolsky(
    b0: np.ndarray, kappa: float, t: float, omega: float
) -> np.ndarray:
    """Two-mode detuned formula from perturbative_regime.tex eq. 344-347.

    g = kappa B0 / 2; Delta = m_h^2 / (4 omega) = kappa^2 B0^2 / (8 omega);
    omega ~ k carrier. The 1/8 coefficient comes from the standard
    relativistic two-mode reduction Delta = (omega_h - omega_a)/2 with
    omega_h - omega_a = m_h^2/(2 omega) and m_h^2 = kappa^2 B0^2 / 2 from
    the B0^2 stress-energy backreaction on the graviton dispersion.
    """
    g = 0.5 * kappa * b0
    delta = -(kappa**2) * (b0**2) / (8.0 * omega)
    om = np.sqrt(g * g + delta * delta)
    return (g * g) / (g * g + delta * delta) * np.sin(om * t) ** 2


def _fit_loglog(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (slope, R²) on log-log axes."""
    lx = np.log10(x)
    ly = np.log10(y)
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(r2)


def _plot(data: dict, out_path: Path) -> None:
    params = data["metadata"]["parameters"]
    kappa = params["kappa"]
    kwave = params["kwave"]

    grid = data["grid"]
    t_values = sorted({r["t_end"] for r in grid})
    t0 = t_values[0]
    rows = sorted([r for r in grid if r["t_end"] == t0], key=operator.itemgetter("B0"))
    b0 = np.array([r["B0"] for r in rows])
    p_final = np.array([r["P_final_sim"] for r in rows])

    # Dense analytic curves
    b0_dense = np.linspace(b0.min(), b0.max(), 400)
    p_bare = np.sin(0.5 * kappa * b0_dense * t0) ** 2
    p_rs = _raffelt_stodolsky(b0_dense, kappa, t0, kwave)

    # Convergence rows, grouped by scheme. Older smoke runs may have
    # a single (modal) trace with no `scheme` field; treat that as
    # one scheme labeled 'modal'.
    conv_all = data.get("convergence", [])
    schemes: dict[str, list[dict]] = {}
    for row in conv_all:
        s = row.get("scheme", "modal")
        schemes.setdefault(s, []).append(row)
    for k in schemes:
        schemes[k].sort(key=operator.itemgetter("N"))

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.8))

    # Panel (a): TIDAL residual against RS + analytic sin²-RS gap as smooth curve.
    ax = axes[0]
    p_rs_at = _raffelt_stodolsky(b0, kappa, t0, kwave)
    data_residual = np.maximum(np.abs(p_final - p_rs_at), EPS_MACH)
    analytic_residual = np.maximum(np.abs(p_bare - p_rs), EPS_MACH)

    ax.semilogy(
        b0,
        data_residual,
        marker="o",
        ms=4,
        lw=0,
        color=IBM_PALETTE["blue"],
        label=r"$|P_\mathrm{TIDAL} - P_\mathrm{RS}|$",
    )
    ax.semilogy(
        b0_dense,
        analytic_residual,
        ls="--",
        lw=1.0,
        color="#888",
        alpha=0.7,
        label=r"$|\sin^2(\kappa B_0 t/2) - P_\mathrm{RS}|$",
    )
    ax.set_xlabel(r"$B_0$")
    ax.set_ylabel(r"$|P_{\mathrm{TIDAL}} - P_{\mathrm{RS}}|$")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)
    ax.set_title(rf"(a) residual at $t = {t0:g}$", fontsize=10)

    # Panel (b): multi-solver N-convergence at one regime point.
    # Errors are computed against the corrected two-mode RS prediction
    # (the best-available analytic baseline), NOT against the bare sin^2
    # formula. With the bare baseline, modal/CVODE would plateau at the
    # analytic 2-mode offset |sin^2 - P_RS| ~3.6e-5 and obscure the true
    # numerical precision of the backend; against RS, modal/CVODE drop
    # to ~4e-8 (the next-order (Delta/g)^4 truncation that the 2x2
    # reduction itself drops), giving a meaningful convergence diagnostic.
    ax = axes[1]
    scheme_styles = {
        "modal": (IBM_PALETTE["blue"], "o", "modal"),
        "cvode": (IBM_PALETTE["purple"], "s", "CVODE"),
        "leapfrog_Y2": (IBM_PALETTE["orange"], "^", "leapfrog $Y_2$"),
        "leapfrog_Y4": (IBM_PALETTE["magenta"], "D", "leapfrog $Y_4$"),
    }
    # Resolve the regime point (B0, t_end) once from any convergence row.
    any_conv_row = next(iter(schemes.values()))[0] if schemes else None
    if any_conv_row is None:
        b0_regime, t_regime = 0.0, 0.0
        p_rs_regime = 0.0
    else:
        b0_regime = float(any_conv_row["B0"])
        t_regime = float(any_conv_row["t_end"])
        p_rs_regime = float(
            _raffelt_stodolsky(np.array([b0_regime]), kappa, t_regime, kwave)[0]
        )
    for scheme_name in ("modal", "cvode", "leapfrog_Y2", "leapfrog_Y4"):
        rows = schemes.get(scheme_name)
        if not rows:
            continue
        color, marker, display = scheme_styles[scheme_name]
        ns = np.array([r["N"] for r in rows], dtype=float)
        # Recompute against corrected RS (the data field abs_error_final
        # stores |sim - bare| for backward compatibility with the JSON
        # schema; we want |sim - RS| instead).
        es = np.maximum(
            np.abs(np.array([r["P_final_sim"] for r in rows]) - p_rs_regime),
            EPS_MACH,
        )
        ax.loglog(
            ns,
            es,
            marker=marker,
            ms=5,
            color=color,
            lw=0,
            label=display,
        )
    ax.set_xlabel(r"grid resolution $N$")
    ax.set_ylabel(r"$|P_{\mathrm{final}}^{\mathrm{sim}} - P_{\mathrm{RS}}|$")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(visible=True, which="both", ls=":", alpha=0.3)
    ax.set_xticks([64, 128, 256, 512])
    ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b0_converge = params.get("b0_converge_point", "—")
    ax.set_title(rf"(b) $N$-convergence at $B_0 = {b0_converge:g}$", fontsize=10)

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
