"""``tidal analyze`` — post-hoc analysis of sweep results.

Provides sensitivity analysis (Sobol, Morris) on completed sweep data.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.measurement._sweep_results import SweepResults


def _auto_detect_metric(results: SweepResults) -> str | None:
    """Return the first recognized metric found in results rows."""
    for candidate in ("P_max", "max_energy_error", "L_mix"):
        if results.rows and candidate in results.rows[0]:
            return candidate
    return None


def _run_sensitivity(
    results: SweepResults,
    method: str,
    metric: str,
    n_bootstrap: int,
) -> int:
    """Run the chosen sensitivity analysis and print results."""
    from tidal.measurement._sensitivity import (
        compute_morris_screening,
        compute_sobol_indices,
        format_sensitivity_table,
    )

    try:
        if method == "sobol":
            result = compute_sobol_indices(results, metric, n_bootstrap=n_bootstrap)
        elif method == "morris":
            result = compute_morris_screening(results, metric)
        else:
            print(f"Error: unknown sensitivity method '{method}'", file=sys.stderr)
            return 1
    except (ImportError, ValueError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(format_sensitivity_table(result))
    return 0


def analyze_command(args: Namespace) -> int:
    """Run ``tidal analyze`` subcommand.

    Returns 0 on success, 1 on error.
    """
    from tidal.measurement._sweep_results import SweepResults

    data_path = Path(args.data_path)
    if not data_path.exists():
        print(f"Error: path not found: {data_path}", file=sys.stderr)
        return 1

    try:
        results = SweepResults.from_directory(data_path)
    except (FileNotFoundError, ValueError) as exc:
        print(f"Error loading sweep data: {exc}", file=sys.stderr)
        return 1

    method = getattr(args, "sensitivity", "sobol")
    metric = getattr(args, "metric", None)
    if metric is None:
        metric = _auto_detect_metric(results)
        if metric is None:
            print(
                "Error: --metric is required. "
                f"Available: {', '.join(results.metric_names)}",
                file=sys.stderr,
            )
            return 1

    n_bootstrap = getattr(args, "bootstrap", 100)
    return _run_sensitivity(results, method, metric, n_bootstrap)
