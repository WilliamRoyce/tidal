"""``tidal analyze`` — post-hoc analysis of sweep results.

Provides sensitivity analysis (Sobol, Morris) and critical field extraction
on completed sweep data.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.measurement._critical_field import CriticalFieldResult
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


def _run_critical_field(results: SweepResults, args: Namespace) -> int:
    """Extract critical field from sweep data and write reduced results."""
    from tidal.cli._console import error as _cerror
    from tidal.cli._console import error_with_hint
    from tidal.measurement._critical_field import (
        compute_critical_field,
        compute_reference_threshold,
        critical_field_to_sweep_results,
    )

    field_param: str = args.critical_field
    metric: str = getattr(args, "metric", None) or "P_final"
    threshold: float = getattr(args, "threshold", 0.99)
    interpolate: bool = not getattr(args, "no_interpolate", False)
    output: str | None = getattr(args, "output", None)
    ref_formula: str | None = getattr(args, "reference_formula", None)
    ref_b: float | None = getattr(args, "reference_B", None)

    if output is None:
        error_with_hint(
            "--output is required for critical field analysis",
            hints=[
                "Example: `tidal analyze DIR --critical-field B0 --output results/`"
            ],
        )
        return 1

    # Compute threshold from analytical formula if requested
    if ref_formula is not None:
        if ref_b is None:
            error_with_hint(
                "--reference-B is required when using --reference-formula",
                hints=["Provide a reference B value: --reference-B 0.1"],
            )
            return 1
        try:
            threshold = compute_reference_threshold(
                ref_formula, ref_b, results.fixed_params, results.sim_settings
            )
        except ValueError as exc:
            _cerror(f"reference formula: {exc}")
            return 1
        print(f"Reference formula '{ref_formula}': P_EM(B={ref_b}) = {threshold:.6f}")

    # Run critical field extraction
    try:
        cf_result = compute_critical_field(
            results,
            field_param,
            metric=metric,
            threshold=threshold,
            interpolate=interpolate,
        )
    except ValueError as exc:
        _cerror(str(exc))
        return 1

    # Convert to SweepResults for serialization
    reduced = critical_field_to_sweep_results(cf_result, results)

    # Write output
    out_path = Path(output)
    out_path.mkdir(parents=True, exist_ok=True)
    reduced.to_csv(out_path / "results.csv")
    reduced.to_json(out_path / "results.json")
    reduced.save_sweep_json(out_path / "sweep.json")

    # Print summary
    _print_critical_field_summary(cf_result)
    print(f"\nResults written to {out_path}/")
    return 0


def _print_critical_field_summary(result: CriticalFieldResult) -> None:
    """Print a summary table of critical field results."""
    n_total = len(result.rows)
    n_good = sum(1 for r in result.rows if r.get("crossing_quality") == "good")
    n_coarse = sum(1 for r in result.rows if r.get("crossing_quality") == "coarse")
    n_edge = sum(1 for r in result.rows if r.get("crossing_quality") == "edge")
    n_none = sum(1 for r in result.rows if r.get("crossing_quality") == "none")

    print("\n--- Critical Field Analysis ---")
    print(f"Field parameter: {result.field_param}")
    print(f"Metric: {result.metric} >= {result.threshold}")
    print(f"Outer parameters: {list(result.outer_params.keys())}")
    print(f"Total points: {n_total}")
    print(f"  good: {n_good}  coarse: {n_coarse}  edge: {n_edge}  none: {n_none}")

    # Show B_min range for valid points
    valid_bmins = [
        r["B_min"]
        for r in result.rows
        if r.get("B_min") is not None
        and not (
            isinstance(r["B_min"], float) and r["B_min"] != r["B_min"]
        )  # NaN check
    ]
    if valid_bmins:
        print(f"B_min range: [{min(valid_bmins):.4g}, {max(valid_bmins):.4g}]")
        inv_bmins = [1.0 / b for b in valid_bmins if b > 0]
        if inv_bmins:
            print(f"1/B_min range: [{min(inv_bmins):.4g}, {max(inv_bmins):.4g}]")


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

    # Dispatch: critical field analysis or sensitivity analysis
    critical_field = getattr(args, "critical_field", None)
    if critical_field:
        return _run_critical_field(results, args)

    # Sensitivity analysis (original path)
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
