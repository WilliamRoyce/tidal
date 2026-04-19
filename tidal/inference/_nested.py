"""Nested sampling backend for Bayesian evidence computation.

Uses **PolyChord** (Will Handley et al., 2015) — Fortran-compiled,
slice-sampling, scales to high dimensions, native anesthetic integration.

Install via ``bash scripts/install_polychord.sh`` (requires gfortran).

Visualization uses **anesthetic** (Handley 2019).

References
----------
Handley, W. et al. (2015) "PolyChord: next-generation nested sampling",
    MNRAS 453(4), 4384-4398.
Handley, W. (2019) "anesthetic: nested sampling visualization",
    JOSS 4(37), 1414.
Skilling, J. (2004) "Nested Sampling", AIP Conference Proceedings 735, 395-405.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.inference._results import InferenceResult


# ---------------------------------------------------------------------------
# nlive recommendation
# ---------------------------------------------------------------------------


def recommend_nlive(ndim: int, precision: str = "standard") -> int:
    """Recommend number of live points for a given dimensionality.

    Scaling follows Handley et al. (2015), Section 3.2.

    Parameters
    ----------
    ndim : int
        Number of parameters.
    precision : str
        One of ``"fast"``, ``"standard"``, ``"production"``.

    Returns
    -------
    int
        Recommended ``nlive``.
    """
    if precision == "fast":
        return max(25 * ndim, 50)
    if precision == "production":
        return max(50 * ndim, 250)
    return max(25 * ndim, 100)  # standard


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------


def run_nested_sampling(
    log_likelihood: Callable[..., float],
    prior_transform: Callable[..., Any],
    ndim: int,
    param_names: list[str],
    sampler: str = "polychord",
    nlive: int = 100,
    n_workers: int | None = None,
    *,
    quiet: bool = False,
    **kwargs: Any,
) -> InferenceResult:
    """Run nested sampling with PolyChord.

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
        Backend — only ``"polychord"`` supported (validated for clarity).
    nlive : int
        Number of live points.
    n_workers : int | None
        Number of parallel workers (unused by PolyChord's serial binding;
        retained for future MPI path).
    quiet : bool
        Suppress progress output.
    **kwargs
        PolyChord-specific settings: ``num_repeats``, ``precision_criterion``,
        ``do_clustering``, ``output_dir``, ``file_root``, etc.
    """
    if sampler != "polychord":
        msg = f"Unknown sampler '{sampler}'. Only 'polychord' is supported."
        raise ValueError(msg)

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


# ---------------------------------------------------------------------------
# PolyChord backend
# ---------------------------------------------------------------------------


def _run_polychord(  # noqa: PLR0915
    *,
    log_likelihood: Callable[..., float],
    prior_transform: Callable[..., Any],
    ndim: int,
    param_names: list[str],
    nlive: int,
    n_workers: int | None,  # PolyChord uses MPI, not Pool
    quiet: bool,
    **kwargs: Any,
) -> InferenceResult:
    """Run nested sampling with PolyChord (Handley et al. 2015).

    Requires ``pypolychord`` to be installed (Fortran-compiled).
    Install via: ``bash scripts/install_polychord.sh``
    """
    # --- Pre-flight: detect all-(-inf) likelihoods before PolyChord hangs.
    # PolyChord rejects non-finite initial draws and will spin forever if
    # the likelihood is uniformly -inf (see #270 for the original silent
    # 60-minute hang this guards against).  Runs before the pypolychord
    # import so it still fires in environments without Fortran.
    n_probe = int(kwargs.pop("n_probe", 10))
    probe_rng = np.random.default_rng(0)
    probe_logls = [
        float(log_likelihood(list(prior_transform(list(probe_rng.random(ndim))))))
        for _ in range(n_probe)
    ]
    if all(not math.isfinite(x) for x in probe_logls):
        msg = (
            f"All {n_probe} prior-sample likelihoods returned non-finite logL. "
            "PolyChord cannot build live points. Check: "
            "(a) --baseline-formula references only parameters in the spec or "
            "--param/--prior; (b) simulations succeed across the prior range; "
            "(c) --constraint isn't over-restrictive."
        )
        raise RuntimeError(msg)

    try:
        from pypolychord import run_polychord
    except ImportError as e:
        msg = (
            "pypolychord is required for nested sampling. "
            "Install with: bash scripts/install_polychord.sh "
            f"(requires gfortran). Original error: {e}"
        )
        raise ImportError(msg) from e
    try:
        from pypolychord import PolyChordSettings
    except ImportError:
        from pypolychord.settings import PolyChordSettings

    from tidal.inference._results import InferenceResult

    # --- Output directory ---
    output_dir = kwargs.pop("output_dir", "polychord_output")
    file_root = kwargs.pop("file_root", "tidal")

    # --- PolyChord settings ---
    nderived = kwargs.pop("nderived", 0)
    settings = PolyChordSettings(ndim, nderived)
    settings.nlive = nlive
    settings.base_dir = str(output_dir)
    settings.file_root = file_root
    settings.do_clustering = kwargs.pop("do_clustering", True)
    settings.read_resume = kwargs.pop("read_resume", False)
    settings.feedback = 0 if quiet else 1

    # Slice sampling repetitions (default 5*ndim per PolyChord convention)
    num_repeats = kwargs.pop("num_repeats", None)
    if num_repeats is not None:
        settings.num_repeats = num_repeats

    # Evidence precision criterion
    precision_criterion = kwargs.pop("precision_criterion", None)
    if precision_criterion is not None:
        settings.precision_criterion = precision_criterion

    # Apply any remaining kwargs directly to settings
    for key, val in kwargs.items():
        if hasattr(settings, key):
            setattr(settings, key, val)

    # --- Wrap likelihood/prior for PolyChord interface ---
    def polychord_loglike(theta: list[float]) -> tuple[float, list[float]]:
        logl = log_likelihood(theta)
        return logl, []

    def polychord_prior(hypercube: list[float]) -> list[float]:
        return list(prior_transform(hypercube))

    # --- Run ---
    output = run_polychord(
        polychord_loglike,
        ndim,
        nderived,
        settings,
        polychord_prior,
    )

    # --- Read results via anesthetic (native PolyChord integration) ---
    chain_root = f"{output_dir}/{file_root}"
    try:
        from anesthetic import read_chains

        ns = read_chains(chain_root)
        samples = np.asarray(ns.to_numpy())[:, :ndim]
        logl = np.asarray(ns.logL)
        # `get_weights()` returns a numpy array in anesthetic>=2.0 and a
        # pandas Series in earlier versions; np.asarray handles both.
        weights = np.asarray(ns.get_weights())
        logz = float(ns.logZ())
        logz_err = float(ns.logZ(100).std())
    except ImportError:
        # Fallback: read PolyChord output files directly
        ns = None
        eq_file = f"{chain_root}_equal_weights.txt"
        data = np.loadtxt(eq_file)
        if data.ndim != 2 or data.shape[1] < ndim + 1:
            msg = (
                f"Unexpected PolyChord output shape {data.shape} in {eq_file}; "
                f"expected at least {ndim + 1} columns"
            )
            raise ValueError(msg) from None
        samples = data[:, :ndim]
        logl = data[:, -1]
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
            "chain_root": chain_root,
            # Cache the anesthetic NestedSamples object for importance analysis
            "_anesthetic_samples": ns,
        },
    )
