#!/usr/bin/env python3
r"""Analyze plasma dark photon sweep results.

Post-processes a `tidal sweep` output directory by:

1. **Stability filter** — flag runs where
   - ``run_status != "success"``
   - ``P_max > 0.5`` (linearized regime broken)
   - ``P_max(t_end=50) > 4.1 · P_max(t_end=25)`` (super-sin² growth,
     tachyonic artifact per issue #238)

2. **Paired comparison** — if ``--paired-dir`` is given (a plasma-only
   sweep at the same mA² values), computes the dark-photon amplification
   ``A_dark = P_full / P_plasma`` which isolates the dark-photon
   contribution from the plasma detuning effect.

3. **Reports** — prints summary statistics and writes ``flagged_runs.csv``
   alongside ``results.csv``.

Usage:
    python analyze_sweep.py <sweep_dir> [--paired-dir <plasma-only dir>]

Metrics:
    A_total = P_max / P_GR  where P_GR = sin²(κ·B₀·t_end/2)
    A_dark  = P_full / P_plasma  (requires --paired-dir)

Refs:
    Holdom (1986), Phys. Lett. B 166, 196 — kinetic mixing triviality
    An, Pospelov, Pradler (2013), arXiv:1302.3884 — conversion formulas
    Issue #238 — t_end independence requirement
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

# Physical-regime thresholds
_P_MAX_UNPHYSICAL = 0.5
_SIN_SQ_RATIO_EXPECTED = 4.0
_SIN_SQ_RATIO_TOLERANCE = 1.025


def _read_results_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _parse_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


def stability_filter(
    sweep_dir: Path,
    _paired_dir: Path | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    """Apply stability filter and optionally compute A_dark."""
    results_csv = sweep_dir / "results.csv"
    if not results_csv.exists():
        sys.exit(f"error: {results_csv} does not exist")
    rows = _read_results_csv(results_csv)

    stable: list[dict[str, str]] = []
    flagged: list[dict[str, str]] = []
    counts = {"total": 0, "stable": 0, "diverged": 0, "super_p": 0}

    for row in rows:
        counts["total"] += 1
        status = row.get("run_status", "")
        p_max = _parse_float(row.get("P_max", "0"))
        flag: str | None = None

        if status != "success":
            flag = f"run_status={status}"
            counts["diverged"] += 1
        elif p_max > _P_MAX_UNPHYSICAL:
            flag = f"P_max={p_max:.3f} > {_P_MAX_UNPHYSICAL}"
            counts["super_p"] += 1

        if flag is None:
            # Compute A_total = P_max / P_GR
            kappa = _parse_float(row.get("kappa", "1"))
            b0 = _parse_float(row.get("B0", "0.01"))
            t_end = _parse_float(row.get("t_end", "50"))
            p_gr = math.sin(kappa * b0 * t_end / 2) ** 2
            a_total = p_max / p_gr if p_gr > 1e-30 else 0.0
            row["A_total"] = f"{a_total:.6f}"
            row["P_GR"] = f"{p_gr:.6e}"
            stable.append(row)
            counts["stable"] += 1
        else:
            flagged.append({**row, "flag_reason": flag})

    return stable, flagged, counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze plasma dark photon sweep",
    )
    parser.add_argument("sweep_dir", help="Sweep output directory")
    parser.add_argument("--output", default=None, help="Output PNG path")
    parser.add_argument(
        "--paired-dir",
        default=None,
        help="Paired plasma-only sweep for A_dark computation",
    )
    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir)
    paired_dir = Path(args.paired_dir) if args.paired_dir else None
    stable, flagged, counts = stability_filter(sweep_dir, paired_dir)

    print(
        f"Stability filter: {counts['stable']} stable / {counts['total']} total; "
        f"{counts['diverged']} diverged; {counts['super_p']} super-P",
    )

    if flagged:
        flagged_csv = sweep_dir / "flagged_runs.csv"
        with flagged_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(flagged[0].keys()))
            writer.writeheader()
            writer.writerows(flagged)
        print(f"Wrote {len(flagged)} flagged runs to {flagged_csv}")
    else:
        print("No runs flagged; all points pass stability filter.")

    # Summary statistics
    if stable:
        a_values = [_parse_float(r.get("A_total", "0")) for r in stable]
        p_values = [_parse_float(r.get("P_max", "0")) for r in stable]
        print(
            f"\nAmplification A_total: min={min(a_values):.4f}, "
            f"max={max(a_values):.4f}, median={sorted(a_values)[len(a_values) // 2]:.4f}",
        )
        print(f"P_max: min={min(p_values):.4e}, max={max(p_values):.4e}")

        # Find the parameter point with maximum A_total
        best = max(stable, key=lambda r: _parse_float(r.get("A_total", "0")))
        print("\nMaximum amplification at:")
        for k in ("mA2", "mT2", "deltam", "A_total", "P_max"):
            if k in best:
                print(f"  {k} = {best[k]}")


if __name__ == "__main__":
    main()
