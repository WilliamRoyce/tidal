"""Simple Monte Carlo sampling for parameter exploration.

Samples from prior distributions, evaluates constraints, runs simulations,
and stores log-likelihood for posterior analysis.  This is the first step
in the inference pipeline — proper Bayesian inference via nested sampling
comes later.
"""

from __future__ import annotations

import math
import time
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace
    from pathlib import Path

    from tidal.inference._constraints import ConstraintSet
    from tidal.inference._likelihood import LikelihoodConfig
    from tidal.inference._prior import Prior
    from tidal.inference._results import InferenceResult


def run_monte_carlo(  # noqa: PLR0915
    priors: list[Prior],
    likelihood_config: LikelihoodConfig,
    base_args: Namespace,
    spec_path: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    n_samples: int,
    constraints: ConstraintSet | None = None,
    n_workers: int | None = None,
    seed: int = 42,
    temp_dir: Path | None = None,
    *,
    keep_sims: bool = False,
    quiet: bool = False,
) -> InferenceResult:
    """Run simple Monte Carlo sampling.

    1. Draw ``n_samples`` from the joint prior.
    2. Reject samples that violate constraints (fast, no simulation).
    3. Evaluate log-likelihood for surviving samples.
    4. Return all samples with log-prior and log-likelihood.

    Parameters
    ----------
    priors : list[Prior]
        Marginal prior distributions, one per parameter.
    likelihood_config : LikelihoodConfig
        How to convert simulation metrics to log-likelihood.
    base_args : Namespace
        Base simulation arguments.
    spec_path : Path
        Path to JSON equation specification.
    measurements : set[str]
        Measurement types to compute.
    source, target : tuple[str, ...] | None
        Source/target fields for conversion measurement.
    threshold : float
        Energy conservation threshold.
    n_samples : int
        Number of Monte Carlo samples to draw.
    constraints : ConstraintSet | None
        Hard prior constraints (reject before simulation).
    n_workers : int | None
        Number of parallel workers. None = sequential.
    seed : int
        RNG seed for reproducibility.
    temp_dir : Path | None
        Base directory for simulation output.
    keep_sims : bool
        If True, keep individual simulation directories.
    quiet : bool
        If True, suppress progress output.
    """
    from tidal.inference._likelihood import (
        SimulationLikelihood,
    )
    from tidal.inference._results import InferenceResult

    param_names = [p.name for p in priors]
    rng = np.random.default_rng(seed)

    # 1. Draw samples from joint prior
    all_samples = np.column_stack([p.sample(rng, n_samples) for p in priors])

    # 2. Compute log-prior and check constraints
    log_prior = np.zeros(n_samples)
    constraint_mask = np.ones(n_samples, dtype=bool)

    for i in range(n_samples):
        lp = sum(p.log_prob(float(all_samples[i, j])) for j, p in enumerate(priors))
        log_prior[i] = lp

        if math.isinf(lp):
            constraint_mask[i] = False
            continue

        if constraints:
            params = dict(zip(param_names, all_samples[i], strict=False))
            if not constraints.check(params):
                constraint_mask[i] = False
                log_prior[i] = -math.inf

    n_valid = int(np.sum(constraint_mask))
    n_rejected = n_samples - n_valid

    if not quiet:
        _print_mc_status(n_samples, n_rejected, n_valid)

    # 3. Evaluate log-likelihood for valid samples
    log_likelihood = np.full(n_samples, -math.inf)
    all_metrics: dict[str, np.ndarray] = {}

    if n_valid == 0:
        if not quiet:
            print("  All samples rejected by constraints — no simulations to run.")
        return InferenceResult(
            samples=all_samples,
            log_likelihood=log_likelihood,
            log_prior=log_prior,
            param_names=param_names,
            method="mc",
            metrics=None,
            metadata={"seed": seed, "n_rejected": n_rejected},
        )

    valid_indices = np.where(constraint_mask)[0]
    valid_samples = all_samples[constraint_mask]

    t0 = time.monotonic()

    metadata_valid: list[dict[str, Any]]
    if n_workers and n_workers > 1:
        # Parallel evaluation
        log_likelihoods_valid, metadata_valid = _run_parallel(
            valid_samples=valid_samples,
            base_args=base_args,
            spec_path=spec_path,
            param_names=param_names,
            measurements=measurements,
            source=source,
            target=target,
            threshold=threshold,
            likelihood_config=likelihood_config,
            temp_dir=temp_dir,
            keep_sims=keep_sims,
            n_workers=n_workers,
            quiet=quiet,
        )
    else:
        # Sequential evaluation
        likelihood_fn = SimulationLikelihood(
            base_args=base_args,
            spec_path=spec_path,
            param_names=param_names,
            measurements=measurements,
            source=source,
            target=target,
            threshold=threshold,
            likelihood_config=likelihood_config,
            temp_dir=temp_dir,
            keep_sims=keep_sims,
        )
        log_likelihoods_valid = np.zeros(n_valid)
        metadata_valid = []
        for k, sample in enumerate(valid_samples):
            likelihood_fn._call_count = k  # pyright: ignore[reportPrivateUsage]
            log_likelihoods_valid[k] = likelihood_fn(sample)
            metadata_valid.append(dict(likelihood_fn._last_metadata))  # pyright: ignore[reportPrivateUsage]
            if not quiet and (k + 1) % max(1, n_valid // 10) == 0:
                print(f"  MC progress: {k + 1}/{n_valid} simulations")

    wall_time = time.monotonic() - t0

    # Place valid results back into full arrays
    for k, idx in enumerate(valid_indices):
        log_likelihood[idx] = log_likelihoods_valid[k]

    # Build per-sample metric arrays from collected metadata.  Constraint-
    # rejected samples (not in valid_indices) get a "constraint" status.
    # Numeric fields use object arrays so we can mix floats with NaN.
    run_status_arr: np.ndarray = np.full(n_samples, "constraint", dtype=object)
    tachyonic_excess_arr: np.ndarray = np.full(n_samples, np.nan, dtype=float)
    for k, idx in enumerate(valid_indices):
        meta = metadata_valid[k] if k < len(metadata_valid) else {}
        run_status_arr[idx] = meta.get("run_status", "success")
        if "tachyonic_excess" in meta:
            tachyonic_excess_arr[idx] = float(meta["tachyonic_excess"])
    all_metrics["run_status"] = run_status_arr
    all_metrics["tachyonic_excess"] = tachyonic_excess_arr

    if not quiet:
        n_success = int(np.sum(np.isfinite(log_likelihoods_valid)))
        n_tachy = int(np.sum(run_status_arr == "tachyonic"))
        msg = f"  Simulations: {n_success}/{n_valid} succeeded ({wall_time:.1f}s)"
        if n_tachy > 0:
            msg += f"; {n_tachy} pre-rejected (tachyonic)"
        print(msg)

    return InferenceResult(
        samples=all_samples,
        log_likelihood=log_likelihood,
        log_prior=log_prior,
        param_names=param_names,
        method="mc",
        metrics=all_metrics or None,
        metadata={
            "seed": seed,
            "n_rejected": n_rejected,
            "n_simulated": n_valid,
            "wall_time_s": round(wall_time, 2),
        },
    )


def _run_parallel(
    *,
    valid_samples: np.ndarray,
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
    n_workers: int,
    quiet: bool,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Run likelihood evaluations in parallel using multiprocessing.

    Returns ``(log_likelihoods, per_sample_metadata)`` — workers cannot
    share state, so the metadata dicts are returned as values from each
    likelihood call rather than read from instance attributes.
    """
    from multiprocessing import Pool

    from tidal.inference._likelihood import (
        _likelihood_worker,  # pyright: ignore[reportPrivateUsage]
        _likelihood_worker_init,  # pyright: ignore[reportPrivateUsage]
    )

    config: dict[str, Any] = {
        "base_args": base_args,
        "spec_path": str(spec_path),
        "param_names": param_names,
        "measurements": measurements,
        "source": source,
        "target": target,
        "threshold": threshold,
        "likelihood_config": likelihood_config,
        "temp_dir": str(temp_dir) if temp_dir else None,
        "keep_sims": keep_sims,
        "_counter": 0,
    }

    n_valid = len(valid_samples)
    results = np.zeros(n_valid)
    metadata: list[dict[str, Any]] = [{} for _ in range(n_valid)]

    with Pool(
        n_workers,
        initializer=_likelihood_worker_init,
        initargs=(config,),
    ) as pool:
        for k, (logl, meta) in enumerate(pool.imap(_likelihood_worker, valid_samples)):
            results[k] = logl
            metadata[k] = meta
            if not quiet and (k + 1) % max(1, n_valid // 10) == 0:
                print(f"  MC progress: {k + 1}/{n_valid} simulations")

    return results, metadata


def _print_mc_status(n_total: int, n_rejected: int, n_valid: int) -> None:
    """Print Monte Carlo sampling status."""
    print(f"  Samples drawn: {n_total}")
    if n_rejected > 0:
        pct = 100 * n_rejected / n_total
        print(f"  Rejected by constraints: {n_rejected} ({pct:.0f}%)")
    print(f"  Simulations to run: {n_valid}")
