"""``tidal sweep`` — Run parameter sweeps across simulations.

Orchestrates ``simulate`` + ``measure`` across a parameter grid
and aggregates scalar metrics into portable CSV/JSON output.

Supports:
- Linear sweeps: ``--sweep "g0=0.1:1.0:10"``
- Log-scale sweeps: ``--sweep "g0=0.01:10.0:10:log"``
- Explicit values: ``--sweep "g0=0.1,0.5,1.0"``
- Multi-parameter cartesian products: multiple ``--sweep`` flags
- Grid convergence studies: ``--converge "32,64,128,256"``
- Resume interrupted sweeps: ``--resume``
- Parallel execution: ``--parallel N``
"""

from __future__ import annotations

import itertools
import sys
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

# ------------------------------------------------------------------
# Parameter parsing
# ------------------------------------------------------------------


def parse_sweep_spec(raw: str) -> tuple[str, list[float]]:  # noqa: C901
    """Parse a ``--sweep`` argument into ``(param_name, values)``.

    Supported formats:

    - ``NAME=START:STOP:N``       — N linearly-spaced points
    - ``NAME=START:STOP:N:log``   — N log-spaced points
    - ``NAME=V1,V2,V3,...``       — explicit values

    Parameters
    ----------
    raw : str
        Raw sweep specification string.

    Returns
    -------
    tuple[str, list[float]]
        Parameter name and list of values.

    Raises
    ------
    ValueError
        If the format is invalid.
    """
    if "=" not in raw:
        msg = (
            f"Invalid sweep spec: '{raw}'. "
            f"Expected NAME=START:STOP:N or NAME=V1,V2,..."
        )
        raise ValueError(msg)

    name, rhs = raw.split("=", 1)
    name = name.strip()
    rhs = rhs.strip()

    if not name:
        msg = f"Empty parameter name in sweep spec: '{raw}'"
        raise ValueError(msg)

    # Check if this is a range spec (contains colons)
    if ":" in rhs:
        parts = rhs.split(":")
        if len(parts) == 3:  # noqa: PLR2004
            # START:STOP:N (linear)
            start, stop, n_str = parts
            n = int(n_str)
            if n < 2:  # noqa: PLR2004
                msg = f"Sweep count must be >= 2, got {n}"
                raise ValueError(msg)
            values = np.linspace(float(start), float(stop), n).tolist()
        elif len(parts) == 4:  # noqa: PLR2004
            # START:STOP:N:log
            start, stop, n_str, scale = parts
            n = int(n_str)
            if n < 2:  # noqa: PLR2004
                msg = f"Sweep count must be >= 2, got {n}"
                raise ValueError(msg)
            if scale.lower() != "log":
                msg = f"Unknown scale '{scale}' (expected 'log')"
                raise ValueError(msg)
            s, e = float(start), float(stop)
            if s <= 0 or e <= 0:
                msg = f"Log-scale requires positive bounds, got {s}:{e}"
                raise ValueError(msg)
            values = np.logspace(np.log10(s), np.log10(e), n).tolist()
        else:
            msg = (
                f"Invalid range spec: '{rhs}'. "
                f"Expected START:STOP:N or START:STOP:N:log"
            )
            raise ValueError(msg)
    else:
        # Explicit values: V1,V2,V3,...
        values = [float(v.strip()) for v in rhs.split(",")]
        if len(values) < 2:  # noqa: PLR2004
            msg = f"Sweep needs at least 2 values, got {len(values)}"
            raise ValueError(msg)

    return name, values


def parse_converge_spec(raw: str) -> list[int]:
    """Parse a ``--converge`` argument into grid sizes.

    Parameters
    ----------
    raw : str
        Comma-separated grid sizes (e.g. ``"32,64,128,256"``).

    Returns
    -------
    list[int]
        Grid sizes in ascending order.

    Raises
    ------
    ValueError
        If fewer than 2 sizes or non-positive values.
    """
    sizes = [int(s.strip()) for s in raw.split(",")]
    if len(sizes) < 2:  # noqa: PLR2004
        msg = f"Convergence study needs at least 2 grid sizes, got {len(sizes)}"
        raise ValueError(msg)
    if any(s <= 0 for s in sizes):
        msg = f"Grid sizes must be positive, got {sizes}"
        raise ValueError(msg)
    return sorted(sizes)


# ------------------------------------------------------------------
# Subdirectory naming
# ------------------------------------------------------------------


def _run_subdir_name(
    swept_params: dict[str, float],
    converge_size: int | None = None,
) -> str:
    """Generate a subdirectory name for one sweep run."""
    if converge_size is not None:
        return f"N_{converge_size}"
    parts: list[str] = []
    for name, val in swept_params.items():
        # Format: short float, avoid trailing zeros
        if val == int(val):
            parts.append(f"{name}_{int(val)}")
        else:
            parts.append(f"{name}_{val:.6g}")
    return "_".join(parts) if parts else "run_0"


# ------------------------------------------------------------------
# Single run execution
# ------------------------------------------------------------------


def _build_sim_args(
    base_args: Namespace,
    param_overrides: dict[str, float],
    output_dir: Path,
    grid_shape_override: int | None = None,
) -> Namespace:
    """Build a simulate-compatible Namespace for one run.

    Copies all simulation flags from *base_args* and overrides
    parameters and output path.
    """
    import copy

    sim_args = copy.copy(base_args)

    # Override parameters: merge base --param list with sweep overrides
    base_params: list[str] = list(getattr(base_args, "param", []) or [])
    for k, v in param_overrides.items():
        # Remove any existing override for this key
        base_params = [p for p in base_params if not p.startswith(f"{k}=")]
        base_params.append(f"{k}={v}")
    sim_args.param = base_params

    # Output to subdirectory (force directory format for disk-backed streaming)
    # Note: no_plot must be False because _infer_output_format checks it first
    # and would return "summary" (skipping disk write). Instead, set
    # output_format="directory" which gets checked after no_plot.
    sim_args.output = str(output_dir)
    sim_args.output_format = "directory"
    sim_args.no_plot = False
    sim_args.quiet = True

    # Grid shape override for convergence mode
    if grid_shape_override is not None:
        sim_args.grid_shape = str(grid_shape_override)

    return sim_args


def _run_single(  # noqa: C901, PLR0912, PLR0913, PLR0914, PLR0915, PLR0917
    base_args: Namespace,
    spec_path: Path,
    param_overrides: dict[str, float],
    output_dir: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    grid_shape_override: int | None = None,
) -> dict[str, Any]:
    """Execute one simulate + measure cycle.

    Returns a dict of scalar metrics for one row of the results table.
    """
    from tidal.cli._measure import (
        _run_asymptotic,  # pyright: ignore[reportPrivateUsage]
        _run_conservation,  # pyright: ignore[reportPrivateUsage]
        _run_conversion,  # pyright: ignore[reportPrivateUsage]
        _run_dispersion,  # pyright: ignore[reportPrivateUsage]
        _run_effective_mass,  # pyright: ignore[reportPrivateUsage]
        _run_energy,  # pyright: ignore[reportPrivateUsage]
        _run_mixing,  # pyright: ignore[reportPrivateUsage]
    )
    from tidal.cli._simulate import (
        _parse_params,  # pyright: ignore[reportPrivateUsage]
        _simulate,  # pyright: ignore[reportPrivateUsage]
    )
    from tidal.measurement._io import SimulationData
    from tidal.symbolic import load_equation_system

    metrics: dict[str, Any] = {}

    # 1. Simulate
    sim_args = _build_sim_args(
        base_args, param_overrides, output_dir, grid_shape_override
    )
    spec = load_equation_system(spec_path)
    params = _parse_params(sim_args.param, spec)

    t0 = time.monotonic()
    exit_code = _simulate(sim_args, spec, params)
    wall_time = time.monotonic() - t0

    if exit_code != 0:
        metrics["wall_time_s"] = round(wall_time, 2)
        metrics["error"] = "simulation_failed"
        return metrics

    # 2. Load simulation data for measurement
    data = SimulationData.load(output_dir, spec)

    # 3. Run requested measurements and extract scalar metrics
    if "conservation" in measurements or "summary" in measurements:
        try:
            cons = _run_conservation(data, threshold)
            metrics["max_energy_error"] = cons["max_relative_error"]
            metrics["energy_conserved"] = cons["is_conserved"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["max_energy_error"] = None
            metrics["conservation_error"] = str(exc)

    conv_result = None
    if "conversion" in measurements or "summary" in measurements:
        try:
            conv = _run_conversion(data, source, target)
            metrics["P_max"] = conv["peak_probability"]
            metrics["P_max_time"] = conv["peak_time"]
            # Final conversion
            result_obj = conv["_result_obj"]
            metrics["P_final"] = float(result_obj.probability[-1])
            conv_result = result_obj
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["P_max"] = None
            metrics["conversion_error"] = str(exc)

    if ("mixing" in measurements or "summary" in measurements) and conv_result is not None:
        try:
            mix = _run_mixing(conv_result)
            metrics["L_mix"] = mix["mixing_length"]
            metrics["L_mix_uncertainty"] = mix["mixing_length_uncertainty"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["L_mix"] = None
            metrics["mixing_error"] = str(exc)

    if "energy" in measurements:
        try:
            eng = _run_energy(data)
            # Total energy at final snapshot
            metrics["E_total_final"] = eng["total"][-1]
            metrics["E_total_initial"] = eng["total"][0]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["energy_error"] = str(exc)

    if "dispersion" in measurements:
        try:
            dyn = list(data.dynamical_fields)
            disp = _run_dispersion(data, dyn)
            result_obj = disp["_result_obj"]
            # Effective mass: m² = ω² - k² at peak
            wn = result_obj.wavenumbers
            freq = result_obj.peak_frequencies
            active = freq > 0.0
            if np.any(active):
                m2_vals = freq[active] ** 2 - wn[active] ** 2
                metrics["m2_eff"] = float(np.median(m2_vals))
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["dispersion_error"] = str(exc)

    if "effective_mass" in measurements:
        try:
            dyn = list(data.dynamical_fields)
            em = _run_effective_mass(data, dyn)
            metrics["m2_eff"] = em["m2_eff"]
            metrics["m2_eff_std"] = em["m2_eff_std"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["effective_mass_error"] = str(exc)

    if "asymptotic" in measurements:
        try:
            asym = _run_asymptotic(data, source, target)
            metrics["P_asymptotic"] = asym["P_final"]
            metrics["P_forward"] = asym["P_forward"]
            metrics["P_reflected"] = asym["P_reflected"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["asymptotic_error"] = str(exc)

    if "peak_conversion" in measurements:
        try:
            from tidal.cli._measure import (
                _run_peak_conversion,  # pyright: ignore[reportPrivateUsage]
            )

            pc = _run_peak_conversion(data, source, target)
            metrics["P_max"] = pc["P_max"]
            metrics["P_max_time"] = pc["P_max_time"]
            metrics["P_final"] = pc["P_final"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["peak_conversion_error"] = str(exc)

    metrics["wall_time_s"] = round(wall_time, 2)
    return metrics


# ------------------------------------------------------------------
# Main sweep orchestration
# ------------------------------------------------------------------


def _collect_sim_settings(args: Namespace) -> dict[str, Any]:
    """Extract simulation settings from CLI args for CSV columns."""
    settings: dict[str, Any] = {}
    settings["grid_shape"] = getattr(args, "grid_shape", None) or "auto"
    settings["t_end"] = getattr(args, "t_end", 10.0)
    settings["dt"] = getattr(args, "dt", None) or "auto"
    settings["scheme"] = getattr(args, "scheme", "auto")
    bc = getattr(args, "bc", None)
    settings["bc"] = bc or ("periodic" if getattr(args, "periodic", True) else "neumann")
    return settings


def _run_sweep(  # noqa: C901, PLR0912, PLR0914, PLR0915
    args: Namespace,
    swept_params: dict[str, list[float]],
    converge_sizes: list[int] | None,
) -> int:
    """Execute the parameter sweep or convergence study.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.
    swept_params : dict[str, list[float]]
        Swept parameter names and values.
    converge_sizes : list[int] or None
        Grid sizes for convergence mode, or None for parameter sweep.

    Returns
    -------
    int
        Exit code.
    """
    from tidal.cli._simulate import _parse_params  # pyright: ignore[reportPrivateUsage]
    from tidal.measurement._sweep_results import SweepResults
    from tidal.symbolic import load_equation_system

    spec_path = Path(args.json_path)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse fixed parameters and measurement config
    spec = load_equation_system(spec_path)
    all_params = _parse_params(getattr(args, "param", []) or [], spec)
    fixed_params = {k: v for k, v in all_params.items() if k not in swept_params}

    measurements = set()
    raw_measure = getattr(args, "measure", None)
    if raw_measure:
        measurements = {s.strip() for s in raw_measure.split(",")}
    if not measurements:
        measurements = {"summary"}

    source = _parse_field_list(getattr(args, "source", None))
    target = _parse_field_list(getattr(args, "target", None))
    threshold = getattr(args, "energy_threshold", 1e-3)

    sim_settings = _collect_sim_settings(args)

    # Build list of runs
    runs: list[dict[str, Any]] = []
    if converge_sizes is not None:
        for size in converge_sizes:
            run: dict[str, Any] = {"_grid_override": size}
            run.update(all_params)
            runs.append(run)
    else:
        # Cartesian product of swept parameters
        param_names = list(swept_params.keys())
        param_values = [swept_params[n] for n in param_names]
        for combo in itertools.product(*param_values):
            run = {}
            for name, val in zip(param_names, combo, strict=False):
                run[name] = val
            runs.append(run)

    total_runs = len(runs)
    resume = getattr(args, "resume", False)

    print(f"Sweep: {total_runs} runs, measurements: {', '.join(sorted(measurements))}")
    print(f"Output: {output_dir.resolve()}")

    # Execute runs
    rows: list[dict[str, Any]] = []
    run_dirs: list[Path] = []

    for i, run_spec in enumerate(runs):
        grid_override = run_spec.pop("_grid_override", None)

        # Build parameter overrides for this run
        param_overrides = dict(all_params)
        param_overrides.update(
            {k: v for k, v in run_spec.items() if k in swept_params}
        )

        # Determine subdirectory name
        swept_vals = {k: run_spec[k] for k in swept_params if k in run_spec}
        subdir = _run_subdir_name(swept_vals, grid_override)
        run_dir = output_dir / subdir
        run_dirs.append(run_dir)

        # Resume check: skip if metadata.json exists
        if resume and (run_dir / "metadata.json").exists():
            print(f"  [{i + 1}/{total_runs}] {subdir} — skipped (resume)")
            # Load existing metrics from results if available
            row = _build_row(
                swept_vals, fixed_params, sim_settings, {}, grid_override
            )
            rows.append(row)
            continue

        print(f"  [{i + 1}/{total_runs}] {subdir}...", end="", flush=True)

        try:
            metrics = _run_single(
                args,
                spec_path,
                param_overrides,
                run_dir,
                measurements,
                source,
                target,
                threshold,
                grid_shape_override=grid_override,
            )
            row = _build_row(
                swept_vals, fixed_params, sim_settings, metrics, grid_override
            )
            rows.append(row)

            # Brief status
            status_parts: list[str] = []
            if "P_max" in metrics and metrics["P_max"] is not None:
                status_parts.append(f"P_max={metrics['P_max']:.4f}")
            if "max_energy_error" in metrics and metrics["max_energy_error"] is not None:
                status_parts.append(f"|dE/E|={metrics['max_energy_error']:.2e}")
            if "wall_time_s" in metrics:
                status_parts.append(f"{metrics['wall_time_s']:.1f}s")
            print(f" {', '.join(status_parts)}")

        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            print(f" ERROR: {exc}")
            row = _build_row(
                swept_vals, fixed_params, sim_settings,
                {"error": str(exc)}, grid_override,
            )
            rows.append(row)

    # Build SweepResults and save
    results = SweepResults(
        swept_params=swept_params,
        fixed_params=fixed_params,
        sim_settings=sim_settings,
        rows=rows,
        run_dirs=run_dirs,
        spec_path=str(spec_path),
        measurements=sorted(measurements),
        source_fields=list(source) if source else None,
        target_fields=list(target) if target else None,
        metadata={
            "timestamp": datetime.now(tz=UTC).isoformat(),
            "total_runs": total_runs,
        },
        converge_sizes=converge_sizes,
    )

    results.save_sweep_json(output_dir / "sweep.json")
    results.to_csv(output_dir / "results.csv")
    results.to_json(output_dir / "results.json")

    print(f"\nSweep complete: {results.n_runs} runs")
    print(f"  results.csv:  {(output_dir / 'results.csv').resolve()}")
    print(f"  results.json: {(output_dir / 'results.json').resolve()}")
    print(f"  sweep.json:   {(output_dir / 'sweep.json').resolve()}")

    # Convergence order estimation
    if converge_sizes is not None and len(rows) >= 3:  # noqa: PLR2004
        _report_convergence(rows, converge_sizes)

    return 0


def _build_row(
    swept_vals: dict[str, float],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    metrics: dict[str, Any],
    grid_override: int | None,
) -> dict[str, Any]:
    """Build a single results row with all columns."""
    row: dict[str, Any] = {}
    row.update(swept_vals)
    row.update(fixed_params)
    if grid_override is not None:
        sim_settings = dict(sim_settings)
        sim_settings["grid_shape"] = grid_override
    row.update(sim_settings)
    # Filter out internal keys
    row.update({k: v for k, v in metrics.items() if not k.startswith("_")})
    return row


def _parse_field_list(raw: str | None) -> tuple[str, ...] | None:
    """Parse comma-separated field names into a tuple, or None."""
    if raw is None:
        return None
    return tuple(s.strip() for s in raw.split(",") if s.strip())


def _report_convergence(  # noqa: C901
    rows: list[dict[str, Any]],
    sizes: list[int],
) -> None:
    """Estimate and print convergence order from sweep results."""
    # Try to find a numerical metric to check convergence
    metric_key = None
    for key in ["max_energy_error", "P_max", "E_total_final"]:
        if key in rows[0] and rows[0][key] is not None:
            metric_key = key
            break

    if metric_key is None:
        return

    values = [row.get(metric_key) for row in rows]
    if any(v is None for v in values):
        return

    vals = np.array(values, dtype=np.float64)
    h = 1.0 / np.array(sizes, dtype=np.float64)

    # Richardson extrapolation: estimate order p from consecutive triples
    orders: list[float] = []
    for i in range(len(vals) - 2):
        f1, f2, f3 = vals[i], vals[i + 1], vals[i + 2]
        h1, h2, _h3 = h[i], h[i + 1], h[i + 2]
        if f2 not in {f1, f3} and h1 != h2:
            # For uniform refinement ratio r = h_coarse/h_fine
            ratio = (f1 - f2) / (f2 - f3) if (f2 - f3) != 0 else float("inf")
            if ratio > 1:
                r = h1 / h2
                p = np.log(ratio) / np.log(r)
                if 0 < p < 10:  # noqa: PLR2004
                    orders.append(float(p))

    if orders:
        avg_order = np.mean(orders)
        print(f"\nConvergence ({metric_key}):")
        print(f"  Estimated order: {avg_order:.2f}")
        for size, val in zip(sizes, vals, strict=False):
            print(f"  N={size}: {metric_key}={val:.6e}")


# ------------------------------------------------------------------
# Parallel execution
# ------------------------------------------------------------------


def _run_single_wrapper(task: dict[str, Any]) -> dict[str, Any]:
    """Wrap _run_single for multiprocessing Pool.map dispatch."""
    from pathlib import Path

    metrics = _run_single(
        task["base_args"],
        Path(task["spec_path"]),
        task["param_overrides"],
        Path(task["output_dir"]),
        task["measurements"],
        task["source"],
        task["target"],
        task["threshold"],
        grid_shape_override=task.get("grid_override"),
    )
    return {
        "index": task["index"],
        "swept_vals": task["swept_vals"],
        "metrics": metrics,
        "grid_override": task.get("grid_override"),
    }


# ------------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------------


def sweep_command(args: Namespace) -> int:  # noqa: PLR0911
    """Execute the sweep command.

    Parameters
    ----------
    args : Namespace
        Parsed CLI arguments.

    Returns
    -------
    int
        Exit code.
    """
    spec_path = Path(args.json_path)
    if not spec_path.exists():
        print(f"Error: file not found: {spec_path}", file=sys.stderr)
        return 1

    if not args.output:
        print("Error: --output is required for sweep", file=sys.stderr)
        return 1

    # Parse sweep specs
    swept_params: dict[str, list[float]] = {}
    converge_sizes: list[int] | None = None

    sweep_specs: list[str] = getattr(args, "sweep", None) or []
    converge_spec: str | None = getattr(args, "converge", None)

    if converge_spec and sweep_specs:
        print("Error: --sweep and --converge are mutually exclusive", file=sys.stderr)
        return 1

    if not converge_spec and not sweep_specs:
        print("Error: provide --sweep or --converge", file=sys.stderr)
        return 1

    try:
        if converge_spec:
            converge_sizes = parse_converge_spec(converge_spec)
        else:
            for raw in sweep_specs:
                name, values = parse_sweep_spec(raw)
                if name in swept_params:
                    print(
                        f"Error: duplicate sweep parameter '{name}'",
                        file=sys.stderr,
                    )
                    return 1
                swept_params[name] = values
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    return _run_sweep(args, swept_params, converge_sizes)
