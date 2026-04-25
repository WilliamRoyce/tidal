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
        self._call_count = 0

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
        """
        call_idx = self._call_count
        self._call_count += 1
        return _evaluate_likelihood(
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
) -> float:
    """Core likelihood evaluation (module-level for pickling).

    This function is the integration point with the sweep infrastructure.
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
                return -math.inf

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
            return -math.inf

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

        return compute_log_likelihood(
            float(metric_value),
            likelihood_config,
            eval_params,
        )

    except Exception:  # noqa: BLE001
        import logging

        logging.getLogger("tidal.inference").debug(
            "Likelihood evaluation failed at theta=%s",
            theta,
            exc_info=True,
        )
        return -math.inf

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


def _likelihood_worker(theta: Any) -> float:  # pyright: ignore[reportUnusedFunction]
    """Evaluate likelihood in worker process (picklable module-level function)."""
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
