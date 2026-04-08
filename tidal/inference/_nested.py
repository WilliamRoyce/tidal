"""Nested sampling backends for Bayesian evidence computation.

Supports two backends via the ``--sampler`` flag:

- **dynesty** (default): Pure Python, pip-installable, good for development.
- **polychord**: Fortran-compiled, MPI-parallel, better for HPC production.

Both use the same interface: ``log_likelihood(theta) -> float`` and
``prior_transform(u) -> theta``.  Visualization uses **anesthetic**
(Handley 2019) for both backends.

References
----------
Speagle, J.S. (2020) "dynesty: a dynamic nested sampling package",
    MNRAS 493(3), 3132-3158.
Handley, W. et al. (2015) "PolyChord: next-generation nested sampling",
    MNRAS 453(4), 4384-4398.
Skilling, J. (2004) "Nested Sampling", AIP Conference Proceedings 735, 395-405.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.inference._results import InferenceResult


def run_nested_sampling(
    log_likelihood: Callable[..., float],
    prior_transform: Callable[..., Any],
    ndim: int,
    param_names: list[str],
    sampler: str = "dynesty",
    nlive: int = 100,
    dlogz: float = 0.01,
    n_workers: int | None = None,
    seed: int | None = None,
    *,
    dynamic: bool = False,
    quiet: bool = False,
    **kwargs: Any,
) -> InferenceResult:
    """Run nested sampling with the specified backend.

    Parameters
    ----------
    log_likelihood : callable
        Log-likelihood function: ``(theta: ndarray) -> float``.
    prior_transform : callable
        Prior transform: ``(u: ndarray) -> ndarray`` mapping [0,1]^d to
        physical parameter space.
    ndim : int
        Number of parameters.
    param_names : list[str]
        Parameter names.
    sampler : str
        Backend: ``"dynesty"`` or ``"polychord"``.
    nlive : int
        Number of live points.
    dlogz : float
        Evidence tolerance (stopping criterion).
    n_workers : int | None
        Number of parallel workers. None = sequential.
    seed : int | None
        Random seed.
    dynamic : bool
        Use dynamic nested sampling (dynesty only).
    quiet : bool
        Suppress progress output.
    **kwargs
        Additional backend-specific options.
    """
    if sampler == "dynesty":
        return _run_dynesty(
            log_likelihood=log_likelihood,
            prior_transform=prior_transform,
            ndim=ndim,
            param_names=param_names,
            nlive=nlive,
            dlogz=dlogz,
            n_workers=n_workers,
            seed=seed,
            dynamic=dynamic,
            quiet=quiet,
            **kwargs,
        )
    if sampler == "polychord":
        return _run_polychord(
            log_likelihood=log_likelihood,
            prior_transform=prior_transform,
            ndim=ndim,
            param_names=param_names,
            nlive=nlive,
            n_workers=n_workers,
            quiet=quiet,
            **kwargs,
        )
    msg = f"Unknown sampler '{sampler}'. Use 'dynesty' or 'polychord'."
    raise ValueError(msg)


def _run_dynesty(
    *,
    log_likelihood: Callable[..., float],
    prior_transform: Callable[..., Any],
    ndim: int,
    param_names: list[str],
    nlive: int,
    dlogz: float,
    n_workers: int | None,
    seed: int | None,
    dynamic: bool,
    quiet: bool,
    **kwargs: Any,
) -> InferenceResult:
    """Run nested sampling with dynesty."""
    try:
        import dynesty
    except ImportError:
        msg = (
            "dynesty is required for nested sampling. "
            "Install with: pip install tidal[inference]"
        )
        raise ImportError(msg) from None

    from tidal.inference._results import InferenceResult

    pool = None
    queue_size = 1

    if n_workers and n_workers > 1:
        from multiprocessing import Pool

        from tidal.cli._sweep import (
            _set_single_thread_blas,  # pyright: ignore[reportPrivateUsage]
        )

        pool = Pool(n_workers, initializer=_set_single_thread_blas)
        queue_size = n_workers

    rstate = np.random.default_rng(seed) if seed is not None else None

    try:
        if dynamic:
            sampler = dynesty.DynamicNestedSampler(
                log_likelihood,
                prior_transform,
                ndim,
                pool=pool,
                queue_size=queue_size,
                rstate=rstate,
            )
            sampler.run_nested(
                dlogz_init=dlogz,
                print_progress=not quiet,
                **kwargs,
            )
        else:
            sampler = dynesty.NestedSampler(
                log_likelihood,
                prior_transform,
                ndim,
                nlive=nlive,
                bound="multi",
                pool=pool,
                queue_size=queue_size,
                rstate=rstate,
            )
            sampler.run_nested(
                dlogz=dlogz,
                print_progress=not quiet,
                **kwargs,
            )

        results = sampler.results

        # Extract importance weights
        try:
            weights = np.exp(results.logwt - results.logz[-1])
            weights /= weights.sum()
        except (AttributeError, IndexError):
            weights = None

        return InferenceResult(
            samples=results.samples,
            log_likelihood=results.logl,
            log_prior=np.zeros(len(results.logl)),  # uniform in unit cube
            param_names=param_names,
            method="nested",
            log_evidence=float(results.logz[-1]),
            log_evidence_err=float(results.logzerr[-1]),
            weights=weights,
            metadata={
                "sampler": "dynesty",
                "nlive": nlive,
                "dlogz": dlogz,
                "dynamic": dynamic,
                "n_iterations": results.niter,
                "n_calls": sum(results.ncall),
            },
        )
    finally:
        if pool is not None:
            pool.close()
            pool.join()


def _run_polychord(
    *,
    log_likelihood: Callable[..., float],
    prior_transform: Callable[..., Any],
    ndim: int,
    param_names: list[str],
    nlive: int,
    n_workers: int | None,
    quiet: bool,
    **kwargs: Any,
) -> InferenceResult:
    """Run nested sampling with PolyChord.

    Requires ``pypolychord`` to be installed (Fortran-compiled).
    Typically used on HPC systems (Newton server, CSD3).
    """
    try:
        from pypolychord import run_polychord
        from pypolychord.settings import PolyChordSettings
    except ImportError:
        msg = (
            "pypolychord is required for PolyChord backend. "
            "Install from source: pip install git+https://github.com/PolyChord/PolyChordLite "
            "(requires gfortran). For development, use --sampler dynesty instead."
        )
        raise ImportError(msg) from None

    from tidal.inference._results import InferenceResult

    output_dir = kwargs.get("output_dir", "polychord_output")

    settings = PolyChordSettings(ndim, 0)  # 0 derived parameters
    settings.nlive = nlive
    settings.base_dir = str(output_dir)
    settings.file_root = "tidal"
    settings.do_clustering = True
    settings.read_resume = False
    settings.feedback = 0 if quiet else 1

    # PolyChord expects (theta, ndim, nderived) -> (logL, [derived])
    def polychord_loglike(theta: list[float]) -> tuple[float, list[float]]:
        logl = log_likelihood(theta)
        return logl, []

    def polychord_prior(hypercube: list[float]) -> list[float]:
        return prior_transform(hypercube)

    output = run_polychord(polychord_loglike, ndim, 0, settings, polychord_prior)

    # Read results via anesthetic for consistent interface
    try:
        from anesthetic import NestedSamples

        ns = NestedSamples(root=f"{output_dir}/tidal")
        samples = ns.to_numpy()[:, :ndim]
        logl = ns.logL.to_numpy()
        weights = ns.get_weights().to_numpy()
        logz = float(ns.logZ())
        logz_err = float(ns.logZ(100).std())
    except ImportError:
        # Fallback: read PolyChord output files directly
        samples = np.loadtxt(f"{output_dir}/tidal_equal_weights.txt")[:, :ndim]
        logl = np.loadtxt(f"{output_dir}/tidal_equal_weights.txt")[:, -1]
        weights = None
        logz = float(output.logZ)
        logz_err = float(output.logZerr)

    return InferenceResult(
        samples=samples,
        log_likelihood=logl,
        log_prior=np.zeros(len(logl)),
        param_names=param_names,
        method="nested",
        log_evidence=logz,
        log_evidence_err=logz_err,
        weights=weights,
        metadata={
            "sampler": "polychord",
            "nlive": nlive,
        },
    )
