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
- Ensemble replicates: ``--n-replicates N --base-seed S``
- IC perturbation: ``--ic-perturbation SCALE``
- Parameter noise UQ: ``--param-noise "PARAM=SCALE"``

References
----------
Smith, R.C. (2013) *Uncertainty Quantification: Theory, Implementation,
and Applications*, SIAM. Ch. 3 (sample statistics, SEM, CI).

Palmer, T.N. et al. (1993) "Ensemble prediction", ECMWF Tech Memo 188.
(IC perturbation as ensemble generation strategy.)

ISO 5725-2:1994 — Accuracy of measurement methods and results.
(CV threshold for repeatability QC.)
"""

from __future__ import annotations

import itertools
import time
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace

    from tidal.symbolic.json_loader import EquationSystem

# Measurement types supported by sweep (subset of tidal measure)
_SWEEP_MEASUREMENTS = frozenset(
    {
        "summary",
        "energy",
        "conversion",
        "mixing",
        "dispersion",
        "conservation",
        "effective_mass",
        "asymptotic",
        "peak_conversion",
        "velocity",
        "resonance",
        "spectrum",
        "spectral_conversion",
    }
)

# Safety limits for sweep grid size
_WARN_SWEEP_SIZE = 1_000
_MAX_SWEEP_SIZE = 10_000

# ------------------------------------------------------------------
# Ensemble / replicate helpers
# ------------------------------------------------------------------

# CV threshold for repeatability QC warning.
# Ref: ISO 5725-2:1994 — Accuracy of measurement methods.
_CV_THRESHOLD = 0.5


def _warn_high_cv(
    agg: Any,  # noqa: ANN401  # SweepResults (circular import)
    swept_params: dict[str, Any],
) -> None:
    """Emit warnings for metrics with high coefficient of variation.

    CV = std/|mean|. A CV > 0.5 indicates large variability relative
    to the mean, suggesting more replicates may be needed.

    References
    ----------
    ISO 5725-2:1994 — Accuracy of measurement methods (repeatability).
    """
    for agg_row in agg.rows:
        for key in agg_row:
            if not key.endswith("_mean") or agg_row[key] is None:
                continue
            base = key[: -len("_mean")]
            mean_val = float(agg_row[key])
            std_val = agg_row.get(f"{base}_std")
            if std_val is None or abs(mean_val) == 0:
                continue
            cv = float(std_val) / abs(mean_val)
            if cv > _CV_THRESHOLD:
                param_info = ", ".join(
                    f"{p}={agg_row.get(p, '?')}" for p in swept_params
                )
                from tidal.cli._console import warn as _cwarn

                _cwarn(
                    f"high variability for {base} "
                    f"(CV={cv:.2f}) at {param_info} "
                    f"— consider more replicates."
                )


def _parse_param_noise(raw_list: list[str] | None) -> dict[str, float]:
    """Parse ``--param-noise PARAM=SCALE`` entries into ``{name: scale}``.

    Raises
    ------
    ValueError
        If an entry doesn't contain '='.
    """
    if not raw_list:
        return {}
    result: dict[str, float] = {}
    for entry in raw_list:
        if "=" not in entry:
            msg = f"--param-noise must be PARAM=SCALE, got {entry!r}"
            raise ValueError(msg)
        name, val = entry.split("=", 1)
        result[name.strip()] = float(val.strip())
    return result


def _apply_param_noise_to_plan(
    rp: dict[str, Any],
    param_noise: dict[str, float],
    seed: int,
) -> None:
    """Perturb swept_vals and param_overrides in-place for one replicate plan.

    Each parameter gets an independent RNG derived from the replicate seed
    via ``numpy.random.SeedSequence.spawn()``, so adding/removing a noised
    parameter does not change the draws for other parameters.

    References
    ----------
    numpy.random.SeedSequence — used for reproducible, independent streams.
    """
    # Spawn one independent child seed per parameter (sorted for determinism)
    sorted_params = sorted(param_noise.keys())
    children = np.random.SeedSequence(seed).spawn(len(sorted_params))
    param_rngs = {
        pname: np.random.default_rng(child)
        for pname, child in zip(sorted_params, children, strict=True)
    }

    noised_swept = dict(rp["swept_vals"])
    noised_overrides = dict(rp["param_overrides"])
    nominal_vals: dict[str, float] = {}
    for pname, scale in param_noise.items():
        draw = scale * param_rngs[pname].standard_normal()
        if pname in noised_swept:
            nominal_vals[pname] = noised_swept[pname]
            noised_swept[pname] += draw
        if pname in noised_overrides:
            if pname not in nominal_vals:
                nominal_vals[pname] = noised_overrides[pname]
            noised_overrides[pname] = noised_swept.get(
                pname, noised_overrides[pname] + draw
            )
    rp["swept_vals"] = noised_swept
    rp["param_overrides"] = noised_overrides
    rp["nominal_vals"] = nominal_vals


def _expand_replicates(
    run_plans: list[dict[str, Any]],
    n_replicates: int,
    base_seed: int,
    param_noise: dict[str, float],
) -> list[dict[str, Any]]:
    """Expand run plans with replicate dimension.

    Each original plan becomes *n_replicates* plans, each with a unique
    ``replicate`` index and ``seed``.  If *param_noise* is set, parameter
    values are perturbed per-replicate.

    Returns a new list of run plans (originals are not mutated).
    """
    expanded: list[dict[str, Any]] = []
    for rp in run_plans:
        for rep in range(n_replicates):
            seed = base_seed + rep
            rep_plan = dict(rp)
            rep_plan["replicate"] = rep
            rep_plan["seed"] = seed
            rep_plan["subdir"] = f"{rp['subdir']}/rep_{rep}"
            rep_plan["run_dir"] = rp["run_dir"] / f"rep_{rep}"
            rep_plan["index"] = len(expanded)

            if param_noise:
                _apply_param_noise_to_plan(rep_plan, param_noise, seed)

            expanded.append(rep_plan)
    return expanded


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
            f"Invalid sweep spec: '{raw}'. Expected NAME=START:STOP:N or NAME=V1,V2,..."
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
            values: list[float] = cast(
                "list[float]", np.linspace(float(start), float(stop), n).tolist()
            )
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
            values = cast(
                "list[float]", np.logspace(np.log10(s), np.log10(e), n).tolist()
            )
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
# Multi-D sampling strategies
# ------------------------------------------------------------------


def _generate_samples(
    param_bounds: dict[str, tuple[float, float]],
    n_samples: int,
    strategy: str,
    *,
    seed: int | None = None,
) -> list[dict[str, float]]:
    """Generate parameter samples using space-filling designs.

    Parameters
    ----------
    param_bounds : dict[str, tuple[float, float]]
        Parameter name -> (lower, upper) bounds.
    n_samples : int
        Number of samples to generate.
    strategy : str
        ``"latin_hypercube"`` or ``"sobol"``.
    seed : int or None
        Random seed for reproducibility.

    Returns
    -------
    list[dict[str, float]]
        Each dict maps parameter names to sampled values.

    Raises
    ------
    ValueError
        If *strategy* is unknown or *n_samples* < 1.
    """
    from scipy.stats.qmc import LatinHypercube, Sobol  # type: ignore[import-untyped]

    names = list(param_bounds.keys())
    d = len(names)
    lower = np.array([param_bounds[n][0] for n in names])
    upper = np.array([param_bounds[n][1] for n in names])

    if strategy == "latin_hypercube":
        sampler = LatinHypercube(d=d, seed=seed)  # pyright: ignore[reportCallIssue]
        unit_samples = sampler.random(n=n_samples)
    elif strategy == "sobol":
        sampler = Sobol(d=d, seed=seed)  # pyright: ignore[reportCallIssue]
        unit_samples = sampler.random(n=n_samples)
    else:
        msg = (
            f"Unknown sampling strategy '{strategy}' (expected: latin_hypercube, sobol)"
        )
        raise ValueError(msg)

    # Scale from [0, 1]^d to parameter bounds
    scaled = lower + unit_samples * (upper - lower)

    return [dict(zip(names, row, strict=True)) for row in scaled]


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


def _build_sim_args(  # noqa: PLR0913
    base_args: Namespace,
    param_overrides: dict[str, float],
    output_dir: Path,
    grid_shape_override: int | None = None,
    *,
    replicate_seed: int | None = None,
    ic_perturbation: float | None = None,
) -> Namespace:
    """Build a simulate-compatible Namespace for one run.

    Copies all simulation flags from *base_args* and overrides
    parameters and output path.

    Parameters
    ----------
    replicate_seed : int, optional
        If set, overrides ``ic_noise_seed`` and ``ic_perturbation_seed``
        for ensemble variation across replicates.
    ic_perturbation : float, optional
        If set, enables IC perturbation with this scale.
    """
    import copy

    sim_args = copy.copy(base_args)

    # Clear simulate-specific resume flags — sweep has its own boolean
    # --resume (resume interrupted sweep), which must not leak into
    # _simulate() where --resume expects a directory path string.
    sim_args.resume = None
    sim_args.snapshot = None
    sim_args.t_additional = None

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

    # Ensemble seed injection — each replicate gets independent seeds for
    # IC noise and IC perturbation via SeedSequence.spawn() to avoid
    # correlated randomness when both are active simultaneously.
    if replicate_seed is not None:
        children = np.random.SeedSequence(replicate_seed).spawn(2)
        sim_args.ic_noise_seed = int(children[0].generate_state(1)[0])
        sim_args.ic_perturbation_seed = int(children[1].generate_state(1)[0])

    # IC perturbation scale
    if ic_perturbation is not None:
        sim_args.ic_perturbation = ic_perturbation

    return sim_args


def _simulate_run(  # noqa: PLR0913
    base_args: Namespace,
    spec_path: Path,
    param_overrides: dict[str, float],
    output_dir: Path,
    grid_shape_override: int | None = None,
    *,
    replicate_seed: int | None = None,
    ic_perturbation: float | None = None,
    spec: EquationSystem | None = None,
) -> tuple[int, float, EquationSystem]:
    """Execute a single simulation and return (exit_code, wall_time_s, spec).

    Returning the loaded ``EquationSystem`` allows callers to reuse it
    for measurement without a redundant file read + JSON parse.
    """
    from tidal.cli._simulate import (
        _parse_params,  # pyright: ignore[reportPrivateUsage]
        _simulate,  # pyright: ignore[reportPrivateUsage]
    )

    sim_args = _build_sim_args(
        base_args,
        param_overrides,
        output_dir,
        grid_shape_override,
        replicate_seed=replicate_seed,
        ic_perturbation=ic_perturbation,
    )
    if spec is None:
        from tidal.symbolic import load_equation_system

        spec = load_equation_system(spec_path)
    params = _parse_params(sim_args.param, spec)

    t0 = time.monotonic()
    exit_code = _simulate(sim_args, spec, params)
    wall_time = time.monotonic() - t0
    return exit_code, wall_time, spec


def _measure_run(  # noqa: C901, PLR0912, PLR0913, PLR0914, PLR0915, PLR0917
    run_dir: Path,
    spec_path: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    spec: EquationSystem | None = None,
) -> dict[str, Any]:
    """Run all requested measurements on an existing simulation output.

    This is the single source of truth for measurement extraction in sweeps.
    Called both after a fresh simulation and when resuming a completed run.

    Returns a dict of scalar metrics.
    """
    from tidal.cli._measure import (
        _run_asymptotic,  # pyright: ignore[reportPrivateUsage]
        _run_conservation,  # pyright: ignore[reportPrivateUsage]
        _run_conversion,  # pyright: ignore[reportPrivateUsage]
        _run_dispersion,  # pyright: ignore[reportPrivateUsage]
        _run_effective_mass,  # pyright: ignore[reportPrivateUsage]
        _run_energy,  # pyright: ignore[reportPrivateUsage]
        _run_mixing,  # pyright: ignore[reportPrivateUsage]
        _run_peak_conversion,  # pyright: ignore[reportPrivateUsage]
        _run_resonance,  # pyright: ignore[reportPrivateUsage]
        _run_spectral_conversion,  # pyright: ignore[reportPrivateUsage]
        _run_spectrum,  # pyright: ignore[reportPrivateUsage]
        _run_velocity,  # pyright: ignore[reportPrivateUsage]
    )
    from tidal.measurement._io import SimulationData

    if spec is None:
        from tidal.symbolic import load_equation_system

        spec = load_equation_system(spec_path)
    data = SimulationData.load(run_dir, spec)
    metrics: dict[str, Any] = {}

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
            result_obj = conv["_result_obj"]
            metrics["P_final"] = float(result_obj.probability[-1])
            conv_result = result_obj
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["P_max"] = None
            metrics["conversion_error"] = str(exc)

    if (
        "mixing" in measurements or "summary" in measurements
    ) and conv_result is not None:
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
            metrics["E_total_final"] = eng["total"][-1]
            metrics["E_total_initial"] = eng["total"][0]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["energy_error"] = str(exc)

    if "dispersion" in measurements:
        try:
            dyn = list(data.dynamical_fields)
            disp = _run_dispersion(data, dyn)
            result_obj = disp["_result_obj"]
            wn: np.ndarray[Any, np.dtype[np.floating[Any]]] = result_obj.wavenumbers
            freq: np.ndarray[Any, np.dtype[np.floating[Any]]] = (
                result_obj.peak_frequencies
            )
            active = freq > 0.0
            if np.any(active):
                m2_vals: np.ndarray[Any, np.dtype[np.floating[Any]]] = (
                    freq[active] ** 2 - wn[active] ** 2
                )
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
            metrics["P_transmitted"] = asym["P_transmitted"]
            metrics["P_reflected"] = asym["P_reflected"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["asymptotic_error"] = str(exc)

    if "peak_conversion" in measurements:
        try:
            if conv_result is not None:
                # Reuse conversion result already computed above
                peak_idx = int(np.argmax(conv_result.probability))
                metrics["P_max"] = float(conv_result.probability[peak_idx])
                metrics["P_max_time"] = float(conv_result.times[peak_idx])
                metrics["P_final"] = float(conv_result.probability[-1])
            else:
                pc = _run_peak_conversion(data, source, target)
                metrics["P_max"] = pc["P_max"]
                metrics["P_max_time"] = pc["P_max_time"]
                metrics["P_final"] = pc["P_final"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["peak_conversion_error"] = str(exc)

    if "velocity" in measurements:
        try:
            dyn = list(data.dynamical_fields)
            vel = _run_velocity(data, dyn)
            metrics["v_group_mean"] = vel["v_group_mean"]
            metrics["v_phase_mean"] = vel["v_phase_mean"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["velocity_error"] = str(exc)

    if "resonance" in measurements:
        try:
            res = _run_resonance(data, source, target)
            metrics["n_resonant_modes"] = res["n_resonant_modes"]
            metrics["conversion_bandwidth"] = res["conversion_bandwidth"]
            metrics["peak_conversion_k"] = res["peak_conversion_k"]
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["resonance_error"] = str(exc)

    if "spectrum" in measurements:
        try:
            spec_results = _run_spectrum(data)
            # Aggregate scalars from spectrum of each field
            for fname, field_spec in spec_results.items():
                if isinstance(field_spec, dict) and "final" in field_spec:
                    final = cast("dict[str, Any]", field_spec["final"])
                    power_data: list[float] = final["power"]
                    wn_data: list[float] = final["wavenumbers"]
                    power = np.array(power_data)
                    wn_arr = np.array(wn_data)
                    if power.max() > 0:
                        metrics[f"peak_k_{fname}"] = float(wn_arr[np.argmax(power)])
                        metrics[f"peak_power_{fname}"] = float(power.max())
                        threshold_val = 0.01 * power.max()
                        metrics[f"n_active_modes_{fname}"] = int(
                            np.sum(power > threshold_val)
                        )
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["spectrum_error"] = str(exc)

    if "spectral_conversion" in measurements:
        try:
            sc = _run_spectral_conversion(data, source, target)
            result_obj = sc["_result_obj"]
            final_p_k = result_obj.probability[-1]
            if len(final_p_k) > 0 and np.any(final_p_k > 0):
                metrics["P_k_max"] = float(final_p_k.max())
                peak_idx = int(np.argmax(final_p_k))
                metrics["k_max_conversion"] = float(result_obj.wavenumbers[peak_idx])
                half_max = final_p_k.max() / 2.0
                above_half = final_p_k > half_max
                if np.any(above_half):
                    k_range = result_obj.wavenumbers[above_half]
                    metrics["spectral_conversion_bandwidth"] = float(
                        k_range[-1] - k_range[0]
                    )
                else:
                    metrics["spectral_conversion_bandwidth"] = 0.0
            else:
                metrics["P_k_max"] = 0.0
                metrics["k_max_conversion"] = 0.0
                metrics["spectral_conversion_bandwidth"] = 0.0
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            metrics["spectral_conversion_error"] = str(exc)

    return metrics


def _run_single(  # noqa: PLR0913, PLR0917
    base_args: Namespace,
    spec_path: Path,
    param_overrides: dict[str, float],
    output_dir: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    grid_shape_override: int | None = None,
    *,
    replicate_seed: int | None = None,
    ic_perturbation: float | None = None,
) -> dict[str, Any]:
    """Execute one simulate + measure cycle.

    Returns a dict of scalar metrics for one row of the results table.
    Always includes ``run_status``, ``error_message``, and ``solver_exit_code``.
    """
    # 1. Simulate
    exit_code, wall_time, spec = _simulate_run(
        base_args,
        spec_path,
        param_overrides,
        output_dir,
        grid_shape_override,
        replicate_seed=replicate_seed,
        ic_perturbation=ic_perturbation,
    )

    if exit_code != 0:
        return {
            "run_status": "solver_error",
            "error_message": f"solver exit code {exit_code}",
            "solver_exit_code": exit_code,
            "wall_time_s": round(wall_time, 2),
            "error": "simulation_failed",
        }

    # 2. Measure (reuse spec from simulate — avoids redundant JSON parse)
    metrics = _measure_run(
        output_dir,
        spec_path,
        measurements,
        source,
        target,
        threshold,
        spec=spec,
    )
    metrics["wall_time_s"] = round(wall_time, 2)
    metrics["run_status"] = "success"
    metrics["error_message"] = None
    metrics["solver_exit_code"] = 0
    return metrics


def _measure_existing(  # noqa: PLR0913, PLR0917
    run_dir: Path,
    spec_path: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    spec: EquationSystem | None = None,
) -> dict[str, Any]:
    """Measure an existing run directory with error handling.

    Used by resume logic in sequential, parallel, and adaptive execution
    paths. Wraps ``_measure_run()`` with status tracking and error capture.

    Parameters
    ----------
    spec : EquationSystem or None
        Pre-loaded equation system.  When supplied the JSON is not re-parsed,
        avoiding redundant I/O on sweep resume (one parse instead of N).

    Returns metrics dict with ``run_status``, ``error_message``, and
    ``solver_exit_code`` always set.
    """
    try:
        metrics = _measure_run(
            run_dir,
            spec_path,
            measurements,
            source,
            target,
            threshold,
            spec=spec,
        )
    except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
        return {
            "error": f"resume_measure_failed: {exc}",
            "run_status": "measurement_error",
            "error_message": str(exc)[:200],
            "solver_exit_code": 0,
        }
    else:
        metrics.setdefault("run_status", "success")
        metrics.setdefault("error_message", None)
        metrics.setdefault("solver_exit_code", 0)
        return metrics


# ------------------------------------------------------------------
# Main sweep orchestration
# ------------------------------------------------------------------


def _save_incremental(  # noqa: PLR0913, PLR0917
    output_dir: Path,
    rows: list[dict[str, Any]],
    run_dirs: list[Path],
    swept_params: dict[str, list[float]],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    spec_path: Path,
    converge_sizes: list[int] | None,
) -> None:
    """Save partial results after each completed run (crash recovery)."""
    from tidal.measurement._sweep_results import SweepResults

    partial = SweepResults(
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
            "total_runs": len(rows),
            "partial": True,
        },
        converge_sizes=converge_sizes,
    )
    partial.to_csv(output_dir / "results.csv")
    partial.to_json(output_dir / "results.json")


def _collect_sim_settings(args: Namespace) -> dict[str, Any]:
    """Extract simulation settings from CLI args for CSV columns."""
    settings: dict[str, Any] = {}
    settings["grid_shape"] = getattr(args, "grid_shape", None) or "auto"
    settings["t_end"] = getattr(args, "t_end", 10.0)
    settings["dt"] = getattr(args, "dt", None) or "auto"
    settings["scheme"] = getattr(args, "scheme", "auto")
    bc = getattr(args, "bc", None)
    settings["bc"] = bc or (
        "periodic" if getattr(args, "periodic", True) else "neumann"
    )
    return settings


def _execute_sequential(  # noqa: PLR0913, PLR0914, PLR0917
    args: Namespace,
    spec_path: Path,
    run_plans: list[dict[str, Any]],
    swept_params: dict[str, list[float]],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    converge_sizes: list[int] | None,
    output_dir: Path,
    *,
    resume: bool,
) -> list[dict[str, Any]]:
    """Execute sweep runs sequentially with incremental saving."""
    total_runs = len(run_plans)
    # Write CSV every ~5% of runs (min 1), balancing crash recovery vs I/O
    save_interval = max(1, total_runs // 20)
    rows: list[dict[str, Any]] = []
    run_dirs: list[Path] = [rp["run_dir"] for rp in run_plans]
    sweep_start = time.monotonic()

    # Pre-load spec once for measurement reuse (avoids N re-parses on resume)
    from tidal.symbolic import load_equation_system

    cached_spec = load_equation_system(spec_path)

    for rp in run_plans:
        i = rp["index"]
        run_dir: Path = rp["run_dir"]
        subdir: str = rp["subdir"]
        swept_vals: dict[str, float] = rp["swept_vals"]
        grid_override: int | None = rp["grid_override"]
        param_overrides: dict[str, float] = rp["param_overrides"]
        rep = rp.get("replicate")
        rep_seed = rp.get("seed")
        rep_ic_pert: float | None = rp.get("ic_perturbation")
        nominal_vals: dict[str, float] | None = rp.get("nominal_vals")

        # Resume check: re-measure existing output instead of re-simulating
        if resume and (run_dir / "metadata.json").exists():
            print(
                f"  [{i + 1}/{total_runs}] {subdir} — measuring existing...",
                end="",
                flush=True,
            )
            metrics = _measure_existing(
                run_dir,
                spec_path,
                measurements,
                source,
                target,
                threshold,
                spec=cached_spec,
            )
            if metrics.get("run_status") == "measurement_error":
                print(f" measure error: {metrics.get('error_message', '')}")
            else:
                print(" ok")
            rows.append(
                _build_row(
                    swept_vals,
                    fixed_params,
                    sim_settings,
                    metrics,
                    grid_override,
                    replicate=rep,
                    seed=rep_seed,
                    nominal_vals=nominal_vals,
                )
            )
            _save_incremental(
                output_dir,
                rows,
                run_dirs,
                swept_params,
                fixed_params,
                sim_settings,
                measurements,
                source,
                target,
                spec_path,
                converge_sizes,
            )
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
                replicate_seed=rep_seed,
                ic_perturbation=rep_ic_pert,
            )
            rows.append(
                _build_row(
                    swept_vals,
                    fixed_params,
                    sim_settings,
                    metrics,
                    grid_override,
                    replicate=rep,
                    seed=rep_seed,
                    nominal_vals=nominal_vals,
                )
            )
            _print_status(metrics)
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            print(f" ERROR: {exc}")
            rows.append(
                _build_row(
                    swept_vals,
                    fixed_params,
                    sim_settings,
                    {
                        "error": str(exc),
                        "run_status": "diverged",
                        "error_message": str(exc)[:200],
                        "solver_exit_code": -1,
                    },
                    grid_override,
                    replicate=rep,
                    seed=rep_seed,
                    nominal_vals=nominal_vals,
                )
            )

        # ETA
        _print_eta(len(rows), total_runs, sweep_start)

        # Incremental save for crash recovery (every save_interval runs
        # + always after the last run).  Reduces O(N²) I/O for large sweeps.
        if len(rows) % save_interval == 0 or len(rows) == total_runs:
            _save_incremental(
                output_dir,
                rows,
                run_dirs,
                swept_params,
                fixed_params,
                sim_settings,
                measurements,
                source,
                target,
                spec_path,
                converge_sizes,
            )

    return rows


def _execute_parallel(  # noqa: PLR0913, PLR0914, PLR0917
    args: Namespace,
    spec_path: Path,
    run_plans: list[dict[str, Any]],
    swept_params: dict[str, list[float]],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    converge_sizes: list[int] | None,
    output_dir: Path,
    *,
    resume: bool,
    n_workers: int,
) -> list[dict[str, Any]]:
    """Execute sweep runs in parallel using multiprocessing.Pool."""
    from multiprocessing import Pool

    total_runs = len(run_plans)
    run_dirs = [rp["run_dir"] for rp in run_plans]

    # Build tasks for Pool, handling resume first (sequentially)
    tasks: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = [{} for _ in range(total_runs)]
    completed: set[int] = set()
    sweep_start = time.monotonic()

    # Pre-load spec once for measurement reuse (avoids N re-parses on resume)
    from tidal.symbolic import load_equation_system

    cached_spec = load_equation_system(spec_path)

    for rp in run_plans:
        i = rp["index"]
        run_dir: Path = rp["run_dir"]
        rep = rp.get("replicate")
        rep_seed = rp.get("seed")
        nominal_vals: dict[str, float] | None = rp.get("nominal_vals")

        if resume and (run_dir / "metadata.json").exists():
            # Resume runs measured sequentially (fast — no simulation)
            print(
                f"  [{i + 1}/{total_runs}] {rp['subdir']} — measuring existing...",
                end="",
                flush=True,
            )
            metrics = _measure_existing(
                run_dir,
                spec_path,
                measurements,
                source,
                target,
                threshold,
                spec=cached_spec,
            )
            if metrics.get("run_status") == "measurement_error":
                print(f" measure error: {metrics.get('error_message', '')}")
            else:
                print(" ok")
            rows[i] = _build_row(
                rp["swept_vals"],
                fixed_params,
                sim_settings,
                metrics,
                rp["grid_override"],
                replicate=rep,
                seed=rep_seed,
                nominal_vals=nominal_vals,
            )
            completed.add(i)
            continue

        tasks.append(
            {
                "index": i,
                "base_args": args,
                "spec_path": str(spec_path),
                "param_overrides": rp["param_overrides"],
                "output_dir": str(rp["run_dir"]),
                "measurements": measurements,
                "source": source,
                "target": target,
                "threshold": threshold,
                "grid_override": rp["grid_override"],
                "swept_vals": rp["swept_vals"],
                "subdir": rp["subdir"],
                "replicate": rep,
                "seed": rep_seed,
                "nominal_vals": nominal_vals,
                "ic_perturbation": rp.get("ic_perturbation"),
            }
        )

    if tasks:
        print(f"  Running {len(tasks)} simulations with {n_workers} workers...")
        with Pool(processes=n_workers) as pool:
            for result in pool.imap_unordered(_run_single_wrapper, tasks):
                idx = result["index"]
                metrics = result["metrics"]
                swept_vals = result["swept_vals"]
                grid_override = result["grid_override"]
                rows[idx] = _build_row(
                    swept_vals,
                    fixed_params,
                    sim_settings,
                    metrics,
                    grid_override,
                    replicate=result.get("replicate"),
                    seed=result.get("seed"),
                    nominal_vals=result.get("nominal_vals"),
                )
                completed.add(idx)
                print(f"  [{len(completed)}/{total_runs}] {result['subdir']}", end="")
                _print_status(metrics)
                _print_eta(len(completed), total_runs, sweep_start)
                # Incremental save after each parallel completion
                _save_incremental(
                    output_dir,
                    list(rows),
                    run_dirs,
                    swept_params,
                    fixed_params,
                    sim_settings,
                    measurements,
                    source,
                    target,
                    spec_path,
                    converge_sizes,
                )

    return list(rows)


# ------------------------------------------------------------------
# Adaptive interval refinement (1D)
# ------------------------------------------------------------------


def _interval_scores(
    values: list[float],
    metric_vals: list[float | None],
) -> list[float]:
    """Compute curvature-based interest score for each interval.

    For interval [a, b] with midpoint m, the score is:
        |f(a) - 2*f(m) + f(b)| / (|f(b) - f(a)| + eps)

    This approximates the second derivative normalized by the first,
    identifying intervals with sharp features worth resolving.
    """
    scores: list[float] = []
    eps = 1e-12
    for i in range(len(values) - 1):
        fa = metric_vals[i]
        fb = metric_vals[i + 1]
        if fa is None or fb is None:
            scores.append(0.0)
            continue
        # Curvature estimate from endpoints only (no midpoint yet)
        # Use abs difference as proxy for interest
        scores.append(abs(fb - fa) / (abs(max(fa, fb)) + eps))
    return scores


def _adaptive_run_point(  # noqa: PLR0913, PLR0917
    val: float,
    param_name: str,
    args: Namespace,
    spec_path: Path,
    all_params: dict[str, float],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    energy_threshold: float,
    output_dir: Path,
    run_dirs: list[Path],
    *,
    resume: bool,
    cached_spec: EquationSystem | None = None,
) -> dict[str, Any]:
    """Run a single adaptive point and return the results row."""
    swept_vals = {param_name: val}
    subdir = _run_subdir_name(swept_vals, None)
    run_dir = output_dir / subdir
    run_dirs.append(run_dir)

    param_overrides = dict(all_params)
    param_overrides[param_name] = val

    if resume and (run_dir / "metadata.json").exists():
        print(f"  [adaptive] {subdir} — measuring existing...", end="", flush=True)
        metrics = _measure_existing(
            run_dir,
            spec_path,
            measurements,
            source,
            target,
            energy_threshold,
            spec=cached_spec,
        )
        if metrics.get("run_status") == "measurement_error":
            print(f" measure error: {metrics.get('error_message', '')}")
        else:
            print(" ok")
    else:
        print(f"  [adaptive] {subdir}...", end="", flush=True)
        try:
            metrics = _run_single(
                args,
                spec_path,
                param_overrides,
                run_dir,
                measurements,
                source,
                target,
                energy_threshold,
            )
            _print_status(metrics)
        except (ValueError, TypeError, KeyError, OSError, RuntimeError) as exc:
            print(f" ERROR: {exc}")
            metrics = {
                "error": str(exc),
                "run_status": "diverged",
                "error_message": str(exc)[:200],
                "solver_exit_code": -1,
            }

    return _build_row(swept_vals, fixed_params, sim_settings, metrics, None)


def _detect_adaptive_metric(rows: list[dict[str, Any]]) -> str | None:
    """Auto-detect a suitable metric key from the first row."""
    from tidal.measurement._sweep_results import DEFAULT_METRIC_CANDIDATES

    for candidate in DEFAULT_METRIC_CANDIDATES:
        if candidate in rows[0]:
            return candidate
    return None


def _execute_adaptive(  # noqa: PLR0913, PLR0917
    args: Namespace,
    spec_path: Path,
    param_name: str,
    initial_values: list[float],
    max_count: int,
    metric_key: str | None,
    threshold: float,
    all_params: dict[str, float],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    energy_threshold: float,
    output_dir: Path,
    *,
    resume: bool,
) -> tuple[list[dict[str, Any]], list[Path], dict[str, list[float]]]:
    """Execute adaptive refinement sweep for a single parameter.

    Returns
    -------
    tuple[list[dict], list[Path], dict[str, list[float]]]
        (rows sorted by parameter value, run_dirs, updated swept_params)
    """
    sweep_start = time.monotonic()
    points = sorted(initial_values)
    rows: list[dict[str, Any]] = []
    run_dirs: list[Path] = []

    # Pre-load spec once for measurement reuse (avoids N re-parses on resume)
    from tidal.symbolic import load_equation_system

    cached_spec = load_equation_system(spec_path)

    # Phase 1: Run initial grid
    print(f"  Adaptive: initial grid ({len(points)} points)")
    for val in points:
        rows.append(  # noqa: PERF401 — side effects
            _adaptive_run_point(
                val,
                param_name=param_name,
                args=args,
                spec_path=spec_path,
                all_params=all_params,
                fixed_params=fixed_params,
                sim_settings=sim_settings,
                measurements=measurements,
                source=source,
                target=target,
                energy_threshold=energy_threshold,
                output_dir=output_dir,
                run_dirs=run_dirs,
                resume=resume,
                cached_spec=cached_spec,
            )
        )

    swept_snapshot = {param_name: list(points)}
    _save_incremental(
        output_dir,
        rows,
        run_dirs,
        swept_snapshot,
        fixed_params,
        sim_settings,
        measurements,
        source,
        target,
        spec_path,
        None,
    )

    # Phase 2: Iterative refinement
    if metric_key is None:
        metric_key = _detect_adaptive_metric(rows)
    if metric_key is None:
        print(
            "  Warning: no metric found for adaptive refinement, skipping refinement phase"
        )
        return rows, run_dirs, swept_snapshot

    print(f"  Adaptive: refining on '{metric_key}' (budget: {max_count} total)")

    iteration = 0
    while len(points) < max_count:
        iteration += 1
        scores = _interval_scores(points, [r.get(metric_key) for r in rows])

        if not scores or max(scores) < threshold:
            print(
                f"  Adaptive: converged after {iteration} iterations (max score < {threshold})"
            )
            break

        idx = int(np.argmax(scores))
        midpoint = (points[idx] + points[idx + 1]) / 2.0

        row = _adaptive_run_point(
            midpoint,
            param_name=param_name,
            args=args,
            spec_path=spec_path,
            all_params=all_params,
            fixed_params=fixed_params,
            sim_settings=sim_settings,
            measurements=measurements,
            source=source,
            target=target,
            energy_threshold=energy_threshold,
            output_dir=output_dir,
            run_dirs=run_dirs,
            resume=resume,
            cached_spec=cached_spec,
        )
        points.insert(idx + 1, midpoint)
        rows.insert(idx + 1, row)

        _print_eta(len(points), max_count, sweep_start)
        swept_snapshot = {param_name: list(points)}
        _save_incremental(
            output_dir,
            rows,
            run_dirs,
            swept_snapshot,
            fixed_params,
            sim_settings,
            measurements,
            source,
            target,
            spec_path,
            None,
        )

    print(f"  Adaptive: {len(points)} points total")
    return rows, run_dirs, {param_name: list(points)}


def _print_status(metrics: dict[str, Any]) -> None:
    """Print brief metric status for a completed run."""
    status_parts: list[str] = []
    if "P_max" in metrics and metrics["P_max"] is not None:
        status_parts.append(f"P_max={metrics['P_max']:.4f}")
    if "max_energy_error" in metrics and metrics["max_energy_error"] is not None:
        status_parts.append(f"|dE/E|={metrics['max_energy_error']:.2e}")
    if "wall_time_s" in metrics:
        status_parts.append(f"{metrics['wall_time_s']:.1f}s")
    print(f" {', '.join(status_parts)}")


def _get_provenance() -> dict[str, str]:
    """Collect provenance metadata for sweep.json."""
    import platform
    import subprocess  # noqa: S404

    meta: dict[str, str] = {
        "hostname": platform.node(),
        "platform": platform.platform(),
    }
    try:
        from importlib.metadata import version

        meta["tidal_version"] = version("tidal")
    except Exception:  # noqa: BLE001
        meta["tidal_version"] = "unknown"
    try:
        git_hash = (
            subprocess.check_output(
                ["git", "rev-parse", "--short", "HEAD"],  # noqa: S607
                stderr=subprocess.DEVNULL,
                timeout=5,
            )
            .decode()
            .strip()
        )
        meta["git_hash"] = git_hash
    except (subprocess.SubprocessError, FileNotFoundError):
        pass
    return meta


def _print_eta(completed: int, total: int, start_time: float) -> None:
    """Print ETA if there are remaining runs."""
    remaining = total - completed
    if remaining <= 0 or completed <= 0:
        return
    elapsed = time.monotonic() - start_time
    avg_per_run = elapsed / completed
    eta_s = remaining * avg_per_run
    if eta_s < 120:  # noqa: PLR2004
        print(f"    ETA: ~{eta_s:.0f}s remaining")
    else:
        print(f"    ETA: ~{eta_s / 60:.1f}min remaining")


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

    # Validate swept parameter names against spec
    if swept_params:
        from tidal.cli._inspect import (
            discover_parameters,  # pyright: ignore[reportPrivateUsage]
        )

        try:
            known_params: set[str] = set(discover_parameters(spec).keys())
        except TypeError:
            known_params = set()
        meta_params_raw = spec.metadata.get("parameters", {})
        if isinstance(meta_params_raw, dict):
            known_params |= set(cast("dict[str, Any]", meta_params_raw).keys())
        for name in sorted(set(swept_params.keys()) - known_params):
            from tidal.cli._console import warn as _cwarn

            _cwarn(
                f"swept parameter '{name}' not found in equation spec. Possible typo?"
            )

    measurements: set[str] = set()
    raw_measure = getattr(args, "measure", None)
    if raw_measure:
        measurements = {s.strip() for s in raw_measure.split(",")}
        unknown = measurements - _SWEEP_MEASUREMENTS
        if unknown:
            from tidal.cli._console import error_with_hint

            error_with_hint(
                f"unknown measurement(s): {', '.join(sorted(unknown))}. "
                f"Valid: {', '.join(sorted(_SWEEP_MEASUREMENTS))}",
                [
                    "Check spelling. Available: energy, conversion, mixing, spectrum, conservation"
                ],
            )
            return 1
    if not measurements:
        measurements = {"summary"}

    source = _parse_field_list(getattr(args, "source", None))
    target = _parse_field_list(getattr(args, "target", None))
    threshold = getattr(args, "energy_threshold", 1e-3)

    sim_settings = _collect_sim_settings(args)

    # Build list of runs
    sweep_strategy = getattr(args, "sweep_strategy", None) or "grid"
    n_samples = getattr(args, "n_samples", None)
    runs: list[dict[str, Any]] = []
    if converge_sizes is not None:
        for size in converge_sizes:
            run: dict[str, Any] = {"_grid_override": size}
            run.update(all_params)
            runs.append(run)
    elif (
        sweep_strategy in {"latin_hypercube", "sobol"}
        and n_samples
        and len(swept_params) >= 1
    ):
        # Space-filling design (LHS or Sobol)
        param_bounds = {
            name: (min(vals), max(vals)) for name, vals in swept_params.items()
        }
        samples = _generate_samples(param_bounds, n_samples, sweep_strategy)
        runs.extend(dict(sample) for sample in samples)
    else:
        # Cartesian product of swept parameters
        param_names = list(swept_params.keys())
        param_values = [swept_params[n] for n in param_names]
        for combo in itertools.product(*param_values):
            run = {}
            for name, val in zip(param_names, combo, strict=True):
                run[name] = val
            runs.append(run)

    total_runs = len(runs)

    # Safety limit for large sweeps
    force_large = getattr(args, "force_large_sweep", False)
    if total_runs > _MAX_SWEEP_SIZE and not force_large:
        from tidal.cli._console import error_with_hint

        error_with_hint(
            f"sweep has {total_runs} runs (limit: {_MAX_SWEEP_SIZE}). "
            f"Use --force-large-sweep to override.",
            [
                "Reduce parameter ranges or grid points",
                "Override with `--force-large-sweep`",
            ],
        )
        return 1
    if total_runs > _WARN_SWEEP_SIZE:
        from tidal.cli._console import warn as _cwarn

        _cwarn(f"{total_runs} runs scheduled. This may take a long time.")

    resume = getattr(args, "resume", False)

    parallel = getattr(args, "parallel", None)
    mode_label = f"parallel={parallel}" if parallel and parallel > 1 else "sequential"
    print(
        f"Sweep: {total_runs} runs ({mode_label}), measurements: {', '.join(sorted(measurements))}"
    )
    print(f"Output: {output_dir.resolve()}")

    # Pre-compute run metadata (subdirs, overrides, etc.)
    run_plans: list[dict[str, Any]] = []
    for i, run_spec in enumerate(runs):
        grid_override = run_spec.pop("_grid_override", None)
        param_overrides = dict(all_params)
        param_overrides.update({k: v for k, v in run_spec.items() if k in swept_params})
        swept_vals = {k: run_spec[k] for k in swept_params if k in run_spec}
        subdir = _run_subdir_name(swept_vals, grid_override)
        run_dir = output_dir / subdir
        run_plans.append(
            {
                "index": i,
                "param_overrides": param_overrides,
                "swept_vals": swept_vals,
                "subdir": subdir,
                "run_dir": run_dir,
                "grid_override": grid_override,
            }
        )

    # Ensemble: expand run_plans with replicates
    n_replicates = getattr(args, "n_replicates", None) or 1
    base_seed = getattr(args, "base_seed", None) or 42
    ic_perturbation = getattr(args, "ic_perturbation", None)
    param_noise = _parse_param_noise(getattr(args, "param_noise", None))

    # Validate param-noise names against known parameters
    if param_noise:
        known_noise_targets = set(swept_params.keys()) | set(all_params.keys())
        meta_params = spec.metadata.get("parameters", {})
        if isinstance(meta_params, dict):
            known_noise_targets |= set(meta_params.keys())  # type: ignore[reportUnknownArgumentType]
        for pn_name in sorted(param_noise.keys()):
            if pn_name not in known_noise_targets:
                from tidal.cli._console import warn as _cwarn

                _cwarn(
                    f"--param-noise parameter '{pn_name}' not found "
                    f"in spec or swept params. Noise will have no effect."
                )

    if n_replicates > 1:
        n_points = len(run_plans)
        run_plans = _expand_replicates(run_plans, n_replicates, base_seed, param_noise)
        # Thread ic_perturbation into each plan
        if ic_perturbation is not None:
            for rp in run_plans:
                rp["ic_perturbation"] = ic_perturbation
        total_runs = len(run_plans)
        run_dirs = [rp["run_dir"] for rp in run_plans]
        print(
            f"  Ensemble: {n_points} points x {n_replicates} replicates = {total_runs} runs"
        )

        # Warn if replicates with no variation source
        ic_type = getattr(args, "ic", "gaussian")
        if ic_type != "noise" and ic_perturbation is None and not param_noise:
            from tidal.cli._console import warn as _cwarn

            _cwarn(
                "--n-replicates > 1 with deterministic ICs and no "
                "--ic-perturbation or --param-noise — all replicates will be identical."
            )

        # Re-check safety limit after expansion
        force_large = getattr(args, "force_large_sweep", False)
        if total_runs > _MAX_SWEEP_SIZE and not force_large:
            from tidal.cli._console import error_with_hint

            error_with_hint(
                f"sweep has {total_runs} runs (limit: {_MAX_SWEEP_SIZE}). "
                f"Use --force-large-sweep to override.",
                [
                    "Reduce parameter ranges or grid points",
                    "Override with `--force-large-sweep`",
                ],
            )
            return 1

    # Dry-run: print plan and exit
    if getattr(args, "dry_run", False):
        print(f"\nDry run: {total_runs} runs planned")
        print(f"  Swept: {list(swept_params.keys())}")
        if n_replicates > 1:
            print(f"  Replicates: {n_replicates} (base_seed={base_seed})")
        print(f"  Measurements: {', '.join(sorted(measurements))}")
        preview = run_plans[:10]
        for rp in preview:
            print(f"    {rp['subdir']}: {rp['swept_vals']}")
        if total_runs > len(preview):
            print(f"    ... and {total_runs - len(preview)} more")
        return 0

    # Check for adaptive mode (from TOML --config)
    adaptive_config: dict[str, Any] = getattr(args, "_adaptive_config", None) or {}
    adaptive_param: str | None = None
    adaptive_metric: str | None = None
    adaptive_budget: int = 20
    adaptive_threshold: float = 0.01
    if adaptive_config and len(swept_params) == 1:
        adaptive_param = next(iter(swept_params))
        if adaptive_param in adaptive_config:
            ac: dict[str, Any] = adaptive_config[adaptive_param]
            adaptive_metric = getattr(args, "adaptive_metric", None) or ac.get("metric")
            adaptive_budget = getattr(args, "adaptive_budget", None) or ac.get(
                "max_count", 20
            )
            adaptive_threshold = getattr(args, "adaptive_threshold", None) or ac.get(
                "threshold", 0.01
            )
        else:
            adaptive_param = None

    # Also check CLI-only adaptive (--sweep "g0=0.01:1.0:adaptive")
    if (
        adaptive_param is None
        and getattr(args, "adaptive_metric", None)
        and len(swept_params) == 1
    ):
        adaptive_param = next(iter(swept_params))
        adaptive_metric = args.adaptive_metric
        adaptive_budget = getattr(args, "adaptive_budget", 20)
        adaptive_threshold = getattr(args, "adaptive_threshold", 0.01)

    # Execute runs
    rows: list[dict[str, Any]] = []
    run_dirs: list[Path] = [rp["run_dir"] for rp in run_plans]

    if adaptive_param is not None:
        if n_replicates > 1:
            from tidal.cli._console import warn as _cwarn

            _cwarn(
                "adaptive refinement does not yet support --n-replicates. "
                "Running single realization per point."
            )
        initial_values = swept_params[adaptive_param]
        rows, run_dirs, swept_params = _execute_adaptive(
            args,
            spec_path,
            adaptive_param,
            initial_values,
            adaptive_budget,
            adaptive_metric,
            adaptive_threshold,
            all_params,
            fixed_params,
            sim_settings,
            measurements,
            source,
            target,
            threshold,
            output_dir,
            resume=resume,
        )
    elif parallel and parallel > 1:
        rows = _execute_parallel(
            args,
            spec_path,
            run_plans,
            swept_params,
            fixed_params,
            sim_settings,
            measurements,
            source,
            target,
            threshold,
            converge_sizes,
            output_dir,
            resume=resume,
            n_workers=parallel,
        )
    else:
        rows = _execute_sequential(
            args,
            spec_path,
            run_plans,
            swept_params,
            fixed_params,
            sim_settings,
            measurements,
            source,
            target,
            threshold,
            converge_sizes,
            output_dir,
            resume=resume,
        )

    # Build SweepResults and save
    meta: dict[str, Any] = {
        "timestamp": datetime.now(tz=UTC).isoformat(),
        "total_runs": len(rows),
        "sampling_strategy": (
            "adaptive"
            if adaptive_param
            else sweep_strategy
            if sweep_strategy != "grid"
            else "grid"
        ),
        **_get_provenance(),
    }
    if n_replicates > 1:
        meta["n_replicates"] = n_replicates
        meta["base_seed"] = base_seed
        if ic_perturbation is not None:
            meta["ic_perturbation"] = ic_perturbation
        if param_noise:
            meta["param_noise"] = param_noise

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
        metadata=meta,
        converge_sizes=converge_sizes,
    )

    results.save_sweep_json(output_dir / "sweep.json")
    results.to_csv(output_dir / "results.csv")
    results.to_json(output_dir / "results.json")

    # Save aggregated results if replicates exist
    if results.has_replicates:
        agg = results.aggregate()
        agg.to_csv(output_dir / "results_agg.csv")
        agg.to_json(output_dir / "results_agg.json")

        # QC: warn on high coefficient of variation
        # Ref: ISO 5725-2:1994 — repeatability/reproducibility QC thresholds
        _warn_high_cv(agg, swept_params)

    print(f"\nSweep complete: {results.n_runs} runs")
    print(f"  results.csv:  {(output_dir / 'results.csv').resolve()}")
    print(f"  results.json: {(output_dir / 'results.json').resolve()}")
    if results.has_replicates:
        print(f"  results_agg.csv: {(output_dir / 'results_agg.csv').resolve()}")
    print(f"  sweep.json:   {(output_dir / 'sweep.json').resolve()}")

    # Convergence order estimation
    if converge_sizes is not None and len(rows) >= 3:  # noqa: PLR2004
        _report_convergence(rows, converge_sizes)

    return 0


def _build_row(  # noqa: PLR0913
    swept_vals: dict[str, float],
    fixed_params: dict[str, float],
    sim_settings: dict[str, Any],
    metrics: dict[str, Any],
    grid_override: int | None,
    *,
    replicate: int | None = None,
    seed: int | None = None,
    nominal_vals: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Build a single results row with all columns.

    Column order: swept params → nominal vals → replicate/seed →
    fixed params → sim settings → metrics → status.
    """
    row: dict[str, Any] = {}
    # 1. Swept parameter values
    row.update(swept_vals)
    # 2. Nominal values for param noise (right after swept vals)
    if nominal_vals:
        for pname, nval in nominal_vals.items():
            row[f"{pname}_nominal"] = nval
    # 3. Replicate metadata (before metrics for readability)
    if replicate is not None:
        row["replicate"] = replicate
    if seed is not None:
        row["seed"] = seed
    # 4. Fixed params and sim settings
    row.update(fixed_params)
    if grid_override is not None:
        sim_settings = dict(sim_settings)
        sim_settings["grid_shape"] = grid_override
    row.update(sim_settings)
    # 5. Metrics (filter out internal keys)
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
        for size, val in zip(sizes, vals, strict=True):
            print(f"  N={size}: {metric_key}={val:.6e}")


# ------------------------------------------------------------------
# Parallel execution
# ------------------------------------------------------------------


def _run_single_wrapper(task: dict[str, Any]) -> dict[str, Any]:
    """Wrap _run_single for multiprocessing Pool.map dispatch.

    Must be a module-level function (picklable) for Pool.map.
    """
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
        replicate_seed=task.get("seed"),
        ic_perturbation=task.get("ic_perturbation"),
    )
    return {
        "index": task["index"],
        "swept_vals": task["swept_vals"],
        "subdir": task["subdir"],
        "metrics": metrics,
        "grid_override": task.get("grid_override"),
        "replicate": task.get("replicate"),
        "seed": task.get("seed"),
        "nominal_vals": task.get("nominal_vals"),
    }


# ------------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------------


def sweep_command(args: Namespace) -> int:  # noqa: C901, PLR0911
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
    # --config: load TOML and merge into args
    config_path = getattr(args, "config", None)
    config_swept: dict[str, list[float]] = {}
    config_converge: list[int] | None = None
    if config_path:
        try:
            from tidal.cli._sweep_config import apply_config_to_args, load_sweep_config

            config = load_sweep_config(Path(config_path))
            config_swept, config_converge = apply_config_to_args(config, args)
        except (FileNotFoundError, ValueError) as exc:
            from tidal.cli._console import error as _cerror

            _cerror(str(exc))
            return 1

    from tidal.cli._console import error as _cerror
    from tidal.cli._console import error_with_hint

    spec_path = Path(args.json_path) if getattr(args, "json_path", None) else None
    if spec_path is None:
        error_with_hint(
            "json_path is required (via positional arg or TOML spec)",
            [
                "Example: `tidal sweep spec.json --sweep 'm2=0.1:1:10' --output results/`"
            ],
        )
        return 1
    if not spec_path.exists():
        error_with_hint(
            f"file not found: {spec_path}",
            ["Use `tidal list` to find specs, or `tidal derive` to generate"],
        )
        return 1

    if not getattr(args, "output", None):
        error_with_hint(
            "--output is required for sweep",
            ["Required: `--output results/`"],
        )
        return 1

    # Parse sweep specs from CLI (these override/extend TOML)
    swept_params: dict[str, list[float]] = dict(config_swept)
    converge_sizes: list[int] | None = config_converge

    sweep_specs: list[str] = getattr(args, "sweep", None) or []
    converge_spec: str | None = getattr(args, "converge", None)

    if converge_spec and (sweep_specs or swept_params):
        error_with_hint(
            "--sweep and --converge are mutually exclusive",
            [
                "Choose one: parameter sweep (`--sweep`) or convergence study (`--converge`)"
            ],
        )
        return 1

    try:
        if converge_spec:
            converge_sizes = parse_converge_spec(converge_spec)
            swept_params = {}  # converge overrides TOML sweeps
        elif sweep_specs:
            for raw in sweep_specs:
                name, values = parse_sweep_spec(raw)
                swept_params[name] = values  # CLI overrides TOML for same param
    except ValueError as exc:
        error_with_hint(
            str(exc),
            [
                "Range: `--sweep 'm2=0.1:10:5'`",
                "List: `--sweep 'm2=0.1,1,10'`",
                "Log: `--sweep 'm2=0.01:10:5:log'`",
            ],
        )
        return 1

    if not swept_params and converge_sizes is None:
        error_with_hint(
            "provide --sweep, --converge, or --config with [sweep.*] sections",
            ["Specify `--sweep 'param=start:stop:N'` or `--converge '32,64,128'`"],
        )
        return 1

    return _run_sweep(args, swept_params, converge_sizes)
