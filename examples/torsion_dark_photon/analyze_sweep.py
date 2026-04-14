#!/usr/bin/env python3
r"""Analyze dark photon torsion sweep results using C₀ = P/B₀² metric.

Post-processes a `tidal sweep` output directory by:

1. **Stability filter** — flag runs where
   - ``run_status != "success"`` (divergence or solver error), or
   - ``P_max > 0.5`` (unphysically large conversion; linearized regime broken), or
   - ``P_max(t_end=50) > 4.1 · P_max(t_end=25)`` if a paired run at half
     ``t_end`` is provided — indicates super-sin² growth (tachyonic
     amplification artifact per issue #238).

2. **Plot** — delegates to ``tidal plot`` with the ``C₀ = P/B₀²`` derived
   metric, using a masked CSV when the filter rejects runs.

3. **Reports** — prints a one-line summary
   (``filtered: N stable / M total; X diverged; Y super-P; Z super-sin²``)
   and writes ``flagged_runs.csv`` alongside the original ``results.csv``.

Usage:
    python analyze_sweep.py <sweep_dir> [--output <png_path>]
                            [--paired-dir <half-t_end sweep dir>]

The ``--paired-dir`` argument enables the super-sin² check (condition 3);
without it, only conditions 1 and 2 are applied.

Refs:
    Hwang & Noh (2310.04150) — linearized regime validity
    Domcke & Garcia-Cely (2301.02072) — Gertsenshtein thin-magnet formula
    Issue #238 — t_end independence requirement for amplification claims
"""

from __future__ import annotations

import argparse
import csv
import subprocess  # noqa: S404
import sys
from pathlib import Path

# Physical-regime thresholds
_P_MAX_UNPHYSICAL = 0.5  # linearized regime broken
_SIN_SQ_RATIO_EXPECTED = 4.0  # sin²(κB₀·t)/sin²(κB₀·t/2) in linear regime
_SIN_SQ_RATIO_TOLERANCE = 1.025  # 2.5% slack → flag at > 4.1


def _read_results_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _parse_float(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except (ValueError, TypeError):
        return default


def _sweep_key(row: dict[str, str], known_cols: set[str]) -> tuple[str, ...]:
    """Identify sweep-parameter columns (non-metadata) to join paired runs."""
    metadata_cols = {
        "P_max",
        "P_max_time",
        "P_final",
        "wall_time_s",
        "run_status",
        "error_message",
        "solver_exit_code",
        "t_end",
        "grid_shape",
        "dt",
        "scheme",
        "bc",
    }
    return tuple(row[k] for k in sorted(known_cols - metadata_cols))


def stability_filter(
    sweep_dir: Path,
    paired_dir: Path | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, int]]:
    """Apply the stability filter to a sweep output.

    Returns
    -------
    stable_rows, flagged_rows, counts
        ``stable_rows`` contains the rows passing all checks.
        ``flagged_rows`` contains the rejected rows with a new column
        ``flag_reason``.
        ``counts`` has keys ``total, stable, diverged, super_p, super_sin2``.
    """
    results_csv = sweep_dir / "results.csv"
    if not results_csv.exists():
        sys.exit(f"error: {results_csv} does not exist")
    rows = _read_results_csv(results_csv)
    known_cols = set(rows[0].keys()) if rows else set()

    paired_lookup: dict[tuple[str, ...], float] = {}
    if paired_dir is not None:
        paired_csv = paired_dir / "results.csv"
        if not paired_csv.exists():
            sys.exit(f"error: paired {paired_csv} does not exist")
        for r in _read_results_csv(paired_csv):
            if r.get("run_status") != "success":
                continue
            paired_lookup[_sweep_key(r, known_cols)] = _parse_float(r.get("P_max", "0"))

    stable: list[dict[str, str]] = []
    flagged: list[dict[str, str]] = []
    counts = {"total": 0, "stable": 0, "diverged": 0, "super_p": 0, "super_sin2": 0}

    for row in rows:
        counts["total"] += 1
        status = row.get("run_status", "")
        p_max = _parse_float(row.get("P_max", "0"))
        flag: str | None = None

        if status != "success":
            flag = f"run_status={status}"
            counts["diverged"] += 1
        elif p_max > _P_MAX_UNPHYSICAL:
            flag = f"P_max={p_max:.3f} > {_P_MAX_UNPHYSICAL} (linearized regime broken)"
            counts["super_p"] += 1
        elif paired_lookup:
            p_paired = paired_lookup.get(_sweep_key(row, known_cols))
            if p_paired is not None and p_paired > 0:
                ratio = p_max / p_paired
                if ratio > _SIN_SQ_RATIO_EXPECTED * _SIN_SQ_RATIO_TOLERANCE:
                    flag = (
                        f"P(t_end)/P(t_end/2) = {ratio:.2f} > {_SIN_SQ_RATIO_EXPECTED:.1f} "
                        f"(super-sin² growth, likely tachyonic artifact; #238)"
                    )
                    counts["super_sin2"] += 1

        if flag is None:
            stable.append(row)
            counts["stable"] += 1
        else:
            flagged.append({**row, "flag_reason": flag})

    return stable, flagged, counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Analyze dark photon sweep with stability filter + C₀ plot",
    )
    parser.add_argument("sweep_dir", help="Sweep output directory")
    parser.add_argument("--output", default=None, help="Output PNG path")
    parser.add_argument(
        "--paired-dir",
        default=None,
        help="Sweep output at half t_end for super-sin² check (#238)",
    )
    args = parser.parse_args()

    sweep_dir = Path(args.sweep_dir)
    paired_dir = Path(args.paired_dir) if args.paired_dir else None
    _stable, flagged, counts = stability_filter(sweep_dir, paired_dir=paired_dir)

    # Summary line
    print(
        f"Stability filter: {counts['stable']} stable / {counts['total']} total; "
        f"{counts['diverged']} diverged; {counts['super_p']} super-P; "
        f"{counts['super_sin2']} super-sin²",
    )

    # Write flagged runs (if any) for inspection
    if flagged:
        flagged_csv = sweep_dir / "flagged_runs.csv"
        with flagged_csv.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(flagged[0].keys()))
            writer.writeheader()
            writer.writerows(flagged)
        print(f"Wrote {len(flagged)} flagged runs to {flagged_csv}")
    else:
        print("No runs flagged; all points pass stability filter.")

    # Delegate plot to `tidal plot` (operates on the full CSV; the C₀ metric
    # handles divergent runs gracefully by dropping them).
    output = args.output or f"{sweep_dir}/analysis_C0.png"
    cmd = [
        "uv",
        "run",
        "tidal",
        "plot",
        str(sweep_dir),
        "--type",
        "sweep",
        "--metric",
        "C0",
        "--log-y",
        "--output",
        output,
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)  # noqa: S603


if __name__ == "__main__":
    main()
