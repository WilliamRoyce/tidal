"""``tidal analyze`` — post-hoc analysis of sweep results.

Provides sensitivity analysis (Sobol, Morris) on completed sweep data.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.measurement._sweep_results import SweepResults


def _auto_detect_metric(results: SweepResults) -> str | None:
    """Return the first recognized metric found in results rows."""
    from tidal.measurement._sweep_results import DEFAULT_METRIC_CANDIDATES

    for candidate in DEFAULT_METRIC_CANDIDATES:
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
    from tidal.cli._console import error as _cerror
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
            from tidal.cli._console import error_with_hint

            error_with_hint(
                f"unknown sensitivity method '{method}'",
                hints=["Valid: sobol (global), morris (screening)"],
            )
            return 1
    except (ImportError, ValueError) as exc:
        _cerror(str(exc))
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
        from tidal.cli._console import error_with_hint

        error_with_hint(
            f"path not found: {data_path}",
            hints=["Use `tidal sweep --output` directory"],
        )
        return 1

    try:
        results = SweepResults.from_directory(data_path)
    except (FileNotFoundError, ValueError) as exc:
        from tidal.cli._console import error_with_hint

        error_with_hint(
            f"loading sweep data: {exc}",
            hints=["Ensure sweep.json or results.csv exists in directory"],
        )
        return 1

    method = getattr(args, "sensitivity", "sobol")
    metric = getattr(args, "metric", None)
    if metric is None:
        metric = _auto_detect_metric(results)
        if metric is None:
            from tidal.cli._console import error_with_hint

            error_with_hint(
                f"--metric is required. Available: {', '.join(results.metric_names)}",
                hints=[
                    "Example: `tidal analyze DIR --metric P_max --sensitivity sobol`"
                ],
            )
            return 1

    n_bootstrap = getattr(args, "bootstrap", 100)
    return _run_sensitivity(results, method, metric, n_bootstrap)
