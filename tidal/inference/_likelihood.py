"""Simulation-based likelihood function for Bayesian inference.

Wraps the existing ``_run_single`` + ``_measure_run`` pipeline from
:mod:`tidal.cli._sweep` as a callable ``log_likelihood(theta) -> float``.

Each likelihood evaluation runs one PDE simulation and extracts a scalar
metric (e.g. ``P_max``).  Failed simulations return ``-inf``.
"""

from __future__ import annotations

import math
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from argparse import Namespace


# Hwang-Noh perturbativity limit on the *absolute* conversion
# probability ``P_max``.  A linearised simulation that produces
# ``P_max > 0.5`` has crossed the perturbative validity boundary
# (Hwang & Noh 2002, Phys. Rev. D 65, 124010); the linearised dynamics
# no longer describe the underlying physics and the value is
# unreliable.  This is independent of the *physics quantity* the user
# is fitting: the configured likelihood may be on the raw ``P_max`` or
# on the amplification ratio ``A = P_max / P_GR`` (when
# ``--baseline-formula`` is supplied — see ``compute_log_likelihood``
# below, which returns ``log(P_max / baseline)`` in that case).
# Either way, the gate condition is on the raw ``P_max`` because the
# Hwang-Noh limit is on the absolute conversion probability, not on
# its ratio to baseline.  When a sample's raw ``P_max`` exceeds the
# limit, ``logL = -inf`` is returned with
# ``run_status="non_perturbative"`` so the sample never enters the
# PolyChord live-points pool.  Gaussian / threshold likelihoods are
# exempt because they may legitimately probe parameter regions where
# ``P_max > 0.5`` (e.g.\ fitting a measured ``P_max ≈ 0.6``).
P_MAX_LIMIT: float = 0.5


def _is_non_perturbative(
    metric: str,
    likelihood_type: str,
    metric_value: float,
) -> bool:
    """Inline Hwang-Noh perturbativity gate.

    Returns True iff the metric is ``P_max``, the likelihood is
    ``maximize`` or ``minimize`` (Gaussian / threshold likelihoods
    may legitimately probe ``P_max > 0.5``), and the *raw* value
    exceeds :data:`P_MAX_LIMIT`.  The check is on the raw ``P_max``
    even when the likelihood is on the baseline-normalised
    amplification ratio: the Hwang-Noh limit is on the absolute
    conversion probability, not on its ratio to baseline.
    """
    return (
        metric == "P_max"
        and likelihood_type in {"maximize", "minimize"}
        and float(metric_value) > P_MAX_LIMIT
    )


@dataclass(frozen=True)
class LikelihoodConfig:
    """Configuration for the simulation likelihood.

    Parameters
    ----------
    metric : str
        Name of the metric to extract (e.g. ``"P_max"``, ``"L_mix"``).
    likelihood_type : str
        One of ``"maximize"``, ``"minimize"``, ``"extremize"``,
        ``"gaussian"``, ``"threshold"``.
    target : float
        Target value for gaussian likelihood.
    sigma : float
        Uncertainty for gaussian likelihood.
    min_value : float
        Minimum value for threshold likelihood.
    baseline_formula : str | None
        Baseline formula for ``extremize`` type.  Evaluated per-point
        using the current parameter values + math functions.  Example:
        ``"sin(kappa * B0 * t_end / 2)**2"``.
    """

    metric: str
    likelihood_type: str = "maximize"
    target: float = 0.0
    sigma: float = 1.0
    min_value: float = 0.0
    baseline_formula: str | None = None

    def __post_init__(self) -> None:
        valid = {"gaussian", "threshold", "maximize", "minimize", "extremize"}
        if self.likelihood_type not in valid:
            msg = f"Unknown likelihood type '{self.likelihood_type}'. Must be one of {sorted(valid)}."
            raise ValueError(msg)
        if self.likelihood_type == "gaussian" and self.sigma <= 0:
            msg = f"Gaussian sigma must be positive, got {self.sigma}"
            raise ValueError(msg)
        if self.likelihood_type == "extremize" and not self.baseline_formula:
            msg = (
                "extremize likelihood requires --baseline-formula "
                '(e.g. "sin(kappa * B0 * t_end / 2)**2")'
            )
            raise ValueError(msg)


def parse_likelihood(
    spec: str,
    *,
    baseline_formula: str | None = None,
) -> LikelihoodConfig:
    """Parse a CLI likelihood specification string.

    Format: ``METRIC:TYPE[:ARGS]``

    Examples::

        "P_max:maximize"
        "P_max:minimize"
        "P_max:extremize"          # requires baseline_formula
        "P_max:gaussian:0.5:0.1"
        "P_max:threshold:0.01"

    Parameters
    ----------
    spec : str
        Likelihood specification string.
    baseline_formula : str | None
        Baseline formula for ``extremize`` type (passed separately via
        ``--baseline-formula`` CLI flag).

    Raises
    ------
    ValueError
        If the specification string is malformed.
    """
    parts = spec.split(":")
    if len(parts) < 2:
        msg = f"Likelihood spec needs METRIC:TYPE, got '{spec}'"
        raise ValueError(msg)

    metric = parts[0]
    ltype = parts[1]

    if ltype == "gaussian":
        if len(parts) < 4:
            msg = (
                f"Gaussian likelihood needs METRIC:gaussian:TARGET:SIGMA, got '{spec}'"
            )
            raise ValueError(msg)
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="gaussian",
            target=float(parts[2]),
            sigma=float(parts[3]),
        )
    if ltype == "threshold":
        if len(parts) < 3:
            msg = f"Threshold likelihood needs METRIC:threshold:MIN_VALUE, got '{spec}'"
            raise ValueError(msg)
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="threshold",
            min_value=float(parts[2]),
        )
    if ltype == "maximize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="maximize",
            baseline_formula=baseline_formula,
        )
    if ltype == "minimize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="minimize",
            baseline_formula=baseline_formula,
        )
    if ltype == "extremize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="extremize",
            baseline_formula=baseline_formula,
        )

    msg = (
        f"Unknown likelihood type '{ltype}'. "
        "Use: maximize, minimize, extremize, gaussian, threshold."
    )
    raise ValueError(msg)


def compute_log_likelihood(
    metric_value: float,
    config: LikelihoodConfig,
    eval_params: dict[str, float] | None = None,
) -> float:
    """Compute log-likelihood from a metric value and config.

    Parameters
    ----------
    metric_value : float
        The measured value of the metric from simulation.
    config : LikelihoodConfig
        Likelihood configuration.
    eval_params : dict | None
        Current parameter values for formula evaluation (needed by
        ``extremize`` type to compute the per-point baseline).
    """
    if math.isnan(metric_value) or math.isinf(metric_value):
        return -math.inf

    # Physics sanity: conversion probability P is bounded by 1. A
    # finite-but-huge P_max (e.g. 1e100) indicates the linearized
    # simulation has crossed into the non-perturbative regime (tachyon,
    # ghost, overflow). This is a numerical artefact, not physics — if
    # the modal solver's stability guard didn't catch it, reject here
    # before the unphysical value drives the posterior. (#298 follow-up)
    if config.metric == "P_max" and metric_value > 2.0:
        return -math.inf

    if config.likelihood_type == "gaussian":
        return -0.5 * ((metric_value - config.target) / config.sigma) ** 2

    if config.likelihood_type == "threshold":
        return 0.0 if metric_value >= config.min_value else -math.inf

    # If a baseline formula is provided, use the log-amplification
    # log(A) = log(metric / baseline) as the physics quantity.  This is
    # the Gertsenshtein-normalized conversion (dimensionless, centered
    # on 0 for no-amplification, positive for enhancement, negative for
    # suppression).  Using log-ratio tames the dynamic range so evidence
    # isn't dominated by a few large outliers.
    if config.baseline_formula:
        baseline = _eval_baseline(config.baseline_formula, eval_params)
        if baseline is None or baseline <= 0 or metric_value <= 0:
            return -math.inf
        log_amplification = math.log(metric_value / baseline)
        if config.likelihood_type == "maximize":
            return log_amplification  # find max amplification
        if config.likelihood_type == "minimize":
            return -log_amplification  # find max suppression
        if config.likelihood_type == "extremize":
            return abs(log_amplification)  # find max deviation (both sides)
        # gaussian/threshold fall through to raw-metric branches below

    if config.likelihood_type == "minimize":
        return -metric_value

    if config.likelihood_type == "extremize":
        # extremize requires a baseline (enforced by __post_init__);
        # this branch is unreachable but kept for safety.
        return -math.inf

    # maximize: use metric value directly as log-likelihood
    return metric_value


def _eval_baseline(
    formula: str | None,
    params: dict[str, float] | None,
) -> float | None:
    """Evaluate a baseline formula with the given parameter values.

    Reuses the same ``FORMULA_NAMESPACE`` (sin, cos, sqrt, pi, etc.)
    as the sweep-results derived-columns computation.
    """
    if formula is None:
        return None

    from tidal.cli._simulate import (
        FORMULA_NAMESPACE,  # pyright: ignore[reportPrivateUsage]
    )

    ns: dict[str, object] = {**FORMULA_NAMESPACE}
    if params:
        ns.update(params)

    try:
        return float(
            eval(formula, {"__builtins__": {}}, ns),  # noqa: S307
        )
    except Exception as exc:  # noqa: BLE001
        import logging

        logging.getLogger("tidal.inference").warning(
            "Baseline formula '%s' failed: %s (params: %s)",
            formula,
            exc,
            sorted(ns.keys()) if params else "none",
        )
        return None


class SimulationLikelihood:
    """Callable likelihood that runs a simulation per evaluation.

    This wraps the existing sweep infrastructure to produce a single
    scalar log-likelihood from a parameter vector.

    Parameters
    ----------
    base_args : Namespace
        Base simulation arguments (grid, IC, time, etc.).
    spec_path : Path
        Path to the JSON equation specification.
    param_names : list[str]
        Parameter names matching the theta vector ordering.
    measurements : set[str]
        Measurement types to compute.
    source : tuple[str, ...] | None
        Source field(s) for conversion.
    target : tuple[str, ...] | None
        Target field(s) for conversion.
    threshold : float
        Energy conservation threshold.
    likelihood_config : LikelihoodConfig
        How to convert metrics to log-likelihood.
    temp_dir : Path | None
        Base temp directory. If None, uses system temp.
    keep_sims : bool
        If True, keep simulation output directories.
    output_dir : Path | None
        Inference output directory.  When set, tachyonic-rejected samples
        are appended to ``output_dir/_tachyonic_samples.csv`` so they can
        be shown on corner plots.  Multiple MPI ranks write safely via
        ``fcntl.flock``.  On resume, rows accumulate across rounds.
    """

    def __init__(
        self,
        base_args: Namespace,
        spec_path: Path,
        param_names: list[str],
        measurements: set[str],
        source: tuple[str, ...] | None,
        target: tuple[str, ...] | None,
        threshold: float,
        likelihood_config: LikelihoodConfig,
        temp_dir: Path | None = None,
        *,
        keep_sims: bool = False,
        output_dir: Path | None = None,
    ) -> None:
        self.base_args = base_args
        self.spec_path = spec_path
        self.param_names = param_names
        self.measurements = measurements
        self.source = source
        self.target = target
        self.threshold = threshold
        self.likelihood_config = likelihood_config
        self.temp_dir = temp_dir
        self.keep_sims = keep_sims
        self._tachyonic_path: Path | None = (
            output_dir / "_tachyonic_samples.csv" if output_dir is not None else None
        )
        self._call_count = 0
        # Metadata from the most recent __call__: read by the sequential MC
        # path to populate per-sample run_status / tachyonic_excess.  Parallel
        # workers don't share state — they return the metadata from
        # `_likelihood_worker` directly instead.
        self._last_metadata: dict[str, Any] = {}

    def __call__(self, theta: Any) -> float:
        """Evaluate log-likelihood at parameter vector theta.

        Parameters
        ----------
        theta : array-like
            Parameter values in the order of ``param_names``.

        Returns
        -------
        float
            Log-likelihood value. Returns ``-inf`` for failed simulations.
            Per-sample metadata (run_status, tachyonic_excess, etc.) is
            stashed on ``self._last_metadata`` for the caller to read.
        """
        call_idx = self._call_count
        self._call_count += 1
        logl, meta = _evaluate_likelihood(
            theta=theta,
            base_args=self.base_args,
            spec_path=self.spec_path,
            param_names=self.param_names,
            measurements=self.measurements,
            source=self.source,
            target=self.target,
            threshold=self.threshold,
            likelihood_config=self.likelihood_config,
            temp_dir=self.temp_dir,
            keep_sims=self.keep_sims,
            call_index=call_idx,
        )
        self._last_metadata = meta
        if (
            not math.isfinite(logl)
            and meta.get("run_status") == "tachyonic"
            and self._tachyonic_path is not None
        ):
            _append_tachyonic_row(self._tachyonic_path, self.param_names, theta, meta)
        return logl


def _append_tachyonic_row(
    path: Path,
    param_names: list[str],
    theta: Any,
    meta: dict[str, Any],
) -> None:
    """Append one tachyonic-rejected sample to the sidecar CSV.

    Safe for concurrent MPI ranks writing to a shared RDS path: uses
    ``fcntl.flock`` for exclusive locking while checking/writing the header
    and appending the data row.  On resume the file already exists with
    non-zero size, so the header is skipped and rows accumulate.
    """
    import csv
    import fcntl

    row: list[Any] = [
        *list(theta),
        meta.get("tachyonic_excess", float("nan")),
        meta.get("k_tachyonic", float("nan")),
        meta.get("n_tachyonic_modes", 0),
        meta.get("stability_profile", ""),
        int(bool(meta.get("borderline_stability"))),
    ]
    header = [
        *param_names,
        "tachyonic_excess",
        "k_tachyonic",
        "n_tachyonic_modes",
        "stability_profile",
        "borderline_stability",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", newline="", encoding="utf-8") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            if path.stat().st_size == 0:
                csv.writer(f).writerow(header)
            csv.writer(f).writerow(row)
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def _evaluate_likelihood(
    *,
    theta: Any,
    base_args: Namespace,
    spec_path: Path,
    param_names: list[str],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    likelihood_config: LikelihoodConfig,
    temp_dir: Path | None,
    keep_sims: bool,
    call_index: int,
) -> tuple[float, dict[str, Any]]:
    """Core likelihood evaluation (module-level for pickling).

    This function is the integration point with the sweep infrastructure.

    Returns
    -------
    tuple[float, dict]
        ``(log_likelihood, metadata)`` where ``metadata`` always contains
        a ``run_status`` key with one of:
        ``"success"``, ``"tachyonic"``, ``"simulation_failed"``,
        ``"metric_missing"``, ``"exception"``.  When the rejection is
        tachyonic, ``metadata["tachyonic_excess"]`` carries the maximum
        Re(λ) found (a positive growth rate).  For successful evaluations,
        ``metadata`` may also include the measured metric value.
    """
    from tidal.cli._sweep import (
        _measure_run,  # pyright: ignore[reportPrivateUsage]
        _simulate_run,  # pyright: ignore[reportPrivateUsage]
    )

    # Build parameter overrides
    param_overrides = {name: float(theta[i]) for i, name in enumerate(param_names)}

    # Backend selection: in-memory default (fast path), fall back to disk
    # when the user explicitly requests --likelihood-backend=disk.  The
    # disk path is kept for bisectability / rollback (see issue #269).
    backend = getattr(base_args, "likelihood_backend", "memory")

    # The disk backend still needs a temp dir; in-memory does not.
    run_dir: Path | None = None
    if backend == "disk":
        # Create a unique temp directory for this evaluation.
        # Include PID to avoid collisions between multiprocessing workers
        # (each fork gets its own counter starting at 0).
        import os

        base = temp_dir or Path(tempfile.gettempdir())
        run_dir = base / f"inference_run_{os.getpid()}_{call_index:06d}"
        run_dir.mkdir(parents=True, exist_ok=True)

    # Pre-flight tachyonic guard — mirrors _run_row_inner in _sweep.py.
    # Rejects parameter points where the source-containing block has growing
    # eigenvalues (Re(λ) > threshold) before running any simulation.
    # Uses conservative=True: when cond(V) is too high to trust the IC-coupling
    # filter (typical for CDT/PGT models), counts the growing mode as tachyonic
    # rather than risk a false-negative.  Cost: microseconds per evaluation.
    stability_meta: dict[str, Any] = {}
    if source and target and ({"conversion", "peak_conversion"} & measurements):
        try:
            from tidal.cli._simulate import (
                _parse_params,  # pyright: ignore[reportPrivateUsage]
            )
            from tidal.measurement._stability import check_conversion_stability
            from tidal.solver.grid import GridInfo
            from tidal.symbolic._spec_cache import load_spec_cached

            raw_spec = load_spec_cached(spec_path)
            base_p = _parse_params(
                list(getattr(base_args, "param", []) or []), raw_spec
            )
            params = {**base_p, **param_overrides}
            grid_n_raw = getattr(base_args, "grid_shape", None)
            grid_n = int(grid_n_raw) if grid_n_raw is not None else 64
            bounds_str = getattr(base_args, "bounds", None)
            if isinstance(bounds_str, str):
                bparts = bounds_str.split(":")
                bounds_tuple: tuple[tuple[float, float], ...] = (
                    (float(bparts[0]), float(bparts[1])),
                )
            elif bounds_str is not None:
                bounds_tuple = tuple(bounds_str)
            else:
                bounds_tuple = ((0.0, 50.0),)
            grid = GridInfo(shape=(grid_n,), bounds=bounds_tuple, periodic=(True,))
            ic_wavevector_str = getattr(base_args, "ic_wavevector", None)
            ic_k: float | None = None
            if ic_wavevector_str:
                try:
                    ic_k = float(str(ic_wavevector_str).split(",")[0])
                except (ValueError, IndexError):
                    ic_k = None
            stability = check_conversion_stability(
                raw_spec,
                grid,
                params,
                source=source[0],
                target=target[0],
                ic_wavevector=ic_k,
            )
            if not stability.stable:
                return -math.inf, {
                    "run_status": "tachyonic",
                    "tachyonic_excess": float(stability.max_excess),
                    "k_tachyonic": (
                        float(stability.k_tachyonic)
                        if stability.k_tachyonic is not None
                        else float("nan")
                    ),
                    "n_tachyonic_modes": int(stability.n_tachyonic_modes),
                    "stability_profile": stability.profile_name,
                    "borderline_stability": False,
                }
            # Stable but possibly borderline — surface for the post-hoc
            # audit to prioritise this sample for re-simulation.
            stability_meta = {
                "stability_profile": stability.profile_name,
                "borderline_stability": bool(stability.borderline),
                "tachyonic_excess": float(stability.max_excess),
            }
        except Exception:  # noqa: BLE001
            import logging

            logging.getLogger("tidal.inference").debug(
                "Pre-flight tachyonic check failed (non-critical), "
                "falling through to simulation",
                exc_info=True,
            )

    try:
        if backend == "memory":
            from tidal.cli._sweep import (
                _measure_from_sim_data,  # pyright: ignore[reportPrivateUsage]
                run_inference_step,
            )

            sim_data = run_inference_step(base_args, spec_path, param_overrides)
            spec = sim_data.spec
            metrics = _measure_from_sim_data(
                sim_data,
                measurements,
                source,
                target,
                threshold,
            )
        else:
            # Run simulation (disk path — preserved for rollback)
            assert run_dir is not None  # set when backend == "disk"
            exit_code, _wall_time, spec = _simulate_run(
                base_args,
                spec_path,
                param_overrides,
                run_dir,
            )

            if exit_code != 0:
                return -math.inf, {
                    "run_status": "simulation_failed",
                    "exit_code": int(exit_code),
                }

            # Extract metrics
            metrics = _measure_run(
                run_dir,
                spec_path,
                measurements,
                source,
                target,
                threshold,
                spec=spec,
            )

        # Get the target metric
        metric_value = metrics.get(likelihood_config.metric)
        if metric_value is None:
            return -math.inf, {
                "run_status": "metric_missing",
                "missing_metric": likelihood_config.metric,
            }

        # Inline Hwang-Noh perturbativity gate (#323 Stage C).  See
        # :func:`_is_non_perturbative` for the gate condition.
        if _is_non_perturbative(
            likelihood_config.metric,
            likelihood_config.likelihood_type,
            float(metric_value),
        ):
            return -math.inf, {
                "run_status": "non_perturbative",
                "P_max": float(metric_value),
                **stability_meta,
            }

        # Build eval_params for formula-based likelihoods.
        # Merge order (same as _simulate_run): spec.metadata["parameters"]
        # defaults -> CLI --param overrides -> swept theta.  Without the
        # first two, formulas like "sin(kappa * B0 * t_end / 2)**2" fail
        # with NameError, _eval_baseline returns None, and every logL
        # collapses to -inf (see #270).
        eval_params: dict[str, float] | None = None
        if likelihood_config.baseline_formula:
            from tidal.cli._simulate import (
                _parse_params,  # pyright: ignore[reportPrivateUsage]
            )

            eval_params = _parse_params(
                list(getattr(base_args, "param", []) or []),
                spec,
            )
            eval_params.update(param_overrides)
            for attr in ("t_end", "dt"):
                val = getattr(base_args, attr, None)
                if val is not None:
                    eval_params[attr] = float(val)

        logl = compute_log_likelihood(
            float(metric_value),
            likelihood_config,
            eval_params,
        )
        success_meta: dict[str, Any] = {
            "run_status": "success",
            likelihood_config.metric: float(metric_value),
            **stability_meta,
        }
        # If P_max guard in compute_log_likelihood demoted the metric to
        # -inf (e.g. P_max > 2.0 — non-perturbative regime), record that
        # explicitly so it shows up as a distinct rejection on the corner plot.
        if math.isinf(logl) and logl < 0:
            success_meta["run_status"] = "non_perturbative"
        return logl, success_meta  # noqa: TRY300

    except Exception:  # noqa: BLE001
        import logging

        logging.getLogger("tidal.inference").debug(
            "Likelihood evaluation failed at theta=%s",
            theta,
            exc_info=True,
        )
        return -math.inf, {"run_status": "exception"}

    finally:
        if run_dir is not None and not keep_sims and run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)


# Module-level wrapper for multiprocessing Pool.map
# Follows the same pattern as _run_single_wrapper in _sweep.py
_LIKELIHOOD_CONFIG: dict[str, Any] = {}


def _likelihood_worker_init(config: dict[str, Any]) -> None:  # pyright: ignore[reportUnusedFunction]
    """Initialize worker process with likelihood config."""
    from tidal.cli._sweep import (
        _set_single_thread_blas,  # pyright: ignore[reportPrivateUsage]
    )

    _set_single_thread_blas()
    _LIKELIHOOD_CONFIG.update(config)


def _likelihood_worker(theta: Any) -> tuple[float, dict[str, Any]]:  # pyright: ignore[reportUnusedFunction]
    """Evaluate likelihood in worker process (picklable module-level function).

    Returns ``(log_likelihood, metadata)``.  The metadata dict carries the
    per-sample ``run_status`` and rejection details (e.g.
    ``tachyonic_excess``) so the parent process can populate per-sample
    columns in the result table.
    """
    import threading

    # Thread-safe call counter
    counter_lock = threading.Lock()
    with counter_lock:
        idx = _LIKELIHOOD_CONFIG.get("_counter", 0)
        _LIKELIHOOD_CONFIG["_counter"] = idx + 1

    return _evaluate_likelihood(
        theta=theta,
        base_args=_LIKELIHOOD_CONFIG["base_args"],
        spec_path=Path(_LIKELIHOOD_CONFIG["spec_path"]),
        param_names=_LIKELIHOOD_CONFIG["param_names"],
        measurements=_LIKELIHOOD_CONFIG["measurements"],
        source=_LIKELIHOOD_CONFIG["source"],
        target=_LIKELIHOOD_CONFIG["target"],
        threshold=_LIKELIHOOD_CONFIG["threshold"],
        likelihood_config=_LIKELIHOOD_CONFIG["likelihood_config"],
        temp_dir=Path(_LIKELIHOOD_CONFIG["temp_dir"])
        if _LIKELIHOOD_CONFIG.get("temp_dir")
        else None,
        keep_sims=_LIKELIHOOD_CONFIG.get("keep_sims", False),
        call_index=idx,
    )
