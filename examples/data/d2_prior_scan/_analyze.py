"""Analyse the δ₁ perturbative-range scan and generate cutoffs + figure.

Reads two sweeps:
 * ``delta1_scan/results.csv``     — t_end=10 (campaign default)
 * ``delta1_scan_t20/results.csv`` — t_end=20 (post-hoc t-doubling check)

Emits:
 * ``delta1_perturbative_range.csv`` — three cutoff rows:
     - ``p_max_t10`` — largest |δ₁| with P_max(t=10) < 0.5
     - ``p_max_t20`` — largest |δ₁| with P_max(t=20) < 0.5
     - ``diverge``   — largest |δ₁| with simulation converging
 * ``delta1_scan.png`` — 2-panel P_max vs δ₁ (linear + log y).
"""

from __future__ import annotations

import csv
import math
import operator
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
T10_CSV = HERE / "delta1_scan" / "results.csv"
T20_CSV = HERE / "delta1_scan_t20" / "results.csv"
OUT_CSV = HERE / "delta1_perturbative_range.csv"
OUT_PNG = HERE / "delta1_scan.png"

P_MAX_LIMIT = 0.5  # Hwang-Noh perturbative bound


def _load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _to_float(s: str | None) -> float:
    if s is None or s == "":
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _largest_abs_below_threshold(
    rows: list[dict[str, str]],
    *,
    column: str,
    threshold: float,
) -> tuple[float, float]:
    """Return (delta_lower, delta_upper) bracketing the largest contiguous
    region around δ₁=0 where ``column`` stays below ``threshold``.

    Walks outward from δ₁=0; returns the first crossings on each side. If
    no crossing on a side, returns the table's min/max for that side.
    """
    parsed = [
        (float(r["delta1"]), _to_float(r.get(column)), r.get("run_status", ""))
        for r in rows
    ]
    parsed.sort(key=operator.itemgetter(0))
    deltas = [d for d, _, _ in parsed]
    vals = [v for _, v, _ in parsed]

    def _is_perturbative(v: float, status: str) -> bool:
        if status != "success":
            return False
        if not math.isfinite(v):
            return False
        return v < threshold

    zero_idx = min(range(len(deltas)), key=lambda i: abs(deltas[i]))

    # Walk right
    right = deltas[zero_idx]
    for i in range(zero_idx + 1, len(deltas)):
        if not _is_perturbative(vals[i], parsed[i][2]):
            break
        right = deltas[i]
    # Walk left
    left = deltas[zero_idx]
    for i in range(zero_idx - 1, -1, -1):
        if not _is_perturbative(vals[i], parsed[i][2]):
            break
        left = deltas[i]
    return left, right


def _largest_abs_converging(rows: list[dict[str, str]]) -> tuple[float, float]:
    """Same outward walk, but only checks ``run_status == 'success'``."""
    parsed = sorted(
        [(float(r["delta1"]), r.get("run_status", "")) for r in rows],
        key=operator.itemgetter(0),
    )
    deltas = [d for d, _ in parsed]
    statuses = [s for _, s in parsed]
    zero_idx = min(range(len(deltas)), key=lambda i: abs(deltas[i]))
    right = deltas[zero_idx]
    for i in range(zero_idx + 1, len(deltas)):
        if statuses[i] != "success":
            break
        right = deltas[i]
    left = deltas[zero_idx]
    for i in range(zero_idx - 1, -1, -1):
        if statuses[i] != "success":
            break
        left = deltas[i]
    return left, right


def main() -> None:
    rows10 = _load(T10_CSV)
    rows20 = _load(T20_CSV)

    lo10, hi10 = _largest_abs_below_threshold(
        rows10, column="P_max", threshold=P_MAX_LIMIT
    )
    lo20, hi20 = _largest_abs_below_threshold(
        rows20, column="P_max", threshold=P_MAX_LIMIT
    )
    lo_div, hi_div = _largest_abs_converging(rows10)
    # Take the tighter (innermost) of the two t_end converge windows.
    lo_div_t20, hi_div_t20 = _largest_abs_converging(rows20)
    lo_div = max(lo_div, lo_div_t20)  # least-negative lower bound
    hi_div = min(hi_div, hi_div_t20)  # smallest positive upper bound

    with OUT_CSV.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "delta1_lower", "delta1_upper"])
        w.writerow(["p_max_t10", f"{lo10:.6g}", f"{hi10:.6g}"])
        w.writerow(["p_max_t20", f"{lo20:.6g}", f"{hi20:.6g}"])
        w.writerow(["diverge", f"{lo_div:.6g}", f"{hi_div:.6g}"])

    print(
        f"p_max_t10 : delta1 in [{lo10:+.4f}, {hi10:+.4f}]  "
        f"(|delta*| = {min(abs(lo10), abs(hi10)):.4f})"
    )
    print(
        f"p_max_t20 : delta1 in [{lo20:+.4f}, {hi20:+.4f}]  "
        f"(|delta*| = {min(abs(lo20), abs(hi20)):.4f})"
    )
    print(
        f"diverge   : delta1 in [{lo_div:+.4f}, {hi_div:+.4f}]  "
        f"(|delta*| = {min(abs(lo_div), abs(hi_div)):.4f})"
    )

    # --- Figure ---
    deltas10 = np.asarray([float(r["delta1"]) for r in rows10])
    pmax10 = np.asarray([_to_float(r.get("P_max")) for r in rows10])
    sort10 = np.argsort(deltas10)
    deltas10 = deltas10[sort10]
    pmax10 = pmax10[sort10]
    statuses10 = np.asarray([r.get("run_status") for r in rows10])[sort10]

    deltas20 = np.asarray([float(r["delta1"]) for r in rows20])
    pmax20 = np.asarray([_to_float(r.get("P_max")) for r in rows20])
    sort20 = np.argsort(deltas20)
    deltas20 = deltas20[sort20]
    pmax20 = pmax20[sort20]

    _fig, (ax_lin, ax_log) = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax in (ax_lin, ax_log):
        ax.plot(deltas10, pmax10, "o-", label="t_end = 10", color="C0", markersize=4)
        ax.plot(deltas20, pmax20, "s-", label="t_end = 20", color="C3", markersize=4)
        ax.axhline(
            P_MAX_LIMIT,
            ls="--",
            color="0.4",
            lw=0.8,
            label=f"Hwang-Noh limit (P={P_MAX_LIMIT})",
        )
        ax.axvspan(
            lo10,
            hi10,
            alpha=0.10,
            color="C0",
            label=f"perturbative (t=10): δ₁ ∈ [{lo10:+.3f}, {hi10:+.3f}]",
        )
        ax.axvspan(
            lo20,
            hi20,
            alpha=0.10,
            color="C3",
            label=f"perturbative (t=20): δ₁ ∈ [{lo20:+.3f}, {hi20:+.3f}]",
        )
        # Mark diverged points
        diverged_mask = statuses10 != "success"
        if diverged_mask.any():
            ax.scatter(
                deltas10[diverged_mask],
                np.full(diverged_mask.sum(), 1.0),
                marker="x",
                color="red",
                s=40,
                label="diverged (sim guard fired)",
                zorder=5,
            )
        ax.set_xlabel(r"$\delta_1$")
        ax.set_ylabel(r"$P_\mathrm{max}$")
        ax.legend(loc="best", fontsize=7)
        ax.grid(alpha=0.3)

    ax_lin.set_title("Linear scale")
    ax_log.set_title("Log scale")
    ax_log.set_yscale("log")
    ax_log.set_ylim(1e-6, 5)
    plt.suptitle(
        r"$\delta_1$ perturbative-range scan at $(\alpha_1, \alpha_2, \alpha_3) = (0, -0.6, 0.5)$, "
        r"$\kappa B_0 = 0.01$, $N=256$, $L=100$",
        fontsize=10,
    )
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=120)
    print(f"\nWrote {OUT_CSV}")
    print(f"Wrote {OUT_PNG}")


if __name__ == "__main__":
    main()
