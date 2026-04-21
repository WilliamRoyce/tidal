"""Parameter importance analysis via Kullback-Leibler divergence.

Uses **anesthetic** (Handley 2019) to compute:

- **D_KL**: KL divergence from prior to posterior (total information gain).
- **d_G**: Bayesian model dimensionality (effective constrained parameters).
- **Marginal D_KL**: Per-parameter information gain (which parameters matter).

References
----------
Handley, W. (2019) "anesthetic: nested sampling visualization",
    JOSS 4(37), 1414.
Handley, W. et al. (2015) "PolyChord: next-generation nested sampling",
    MNRAS 453(4), 4384-4398.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from tidal.inference._results import InferenceResult


@dataclass(frozen=True)
class ParameterImportanceResult:
    """Results from parameter importance analysis.

    Parameters
    ----------
    param_names : list[str]
        Parameter names.
    d_kl : float
        Total KL divergence (nats) — information gained from prior to
        posterior.
    d_kl_err : float
        Bootstrap uncertainty on D_KL.
    d_g : float
        Bayesian model dimensionality — effective number of constrained
        parameters.
    d_g_err : float
        Bootstrap uncertainty on d_G.
    marginal_d_kl : dict[str, float]
        Per-parameter marginal KL divergence (nats).  High D_KL means
        the data strongly constrains that parameter.
    log_evidence : float
        Log Bayesian evidence (log Z).
    log_evidence_err : float
        Bootstrap uncertainty on log Z.
    """

    param_names: list[str]
    d_kl: float
    d_kl_err: float
    d_g: float
    d_g_err: float
    marginal_d_kl: dict[str, float]
    log_evidence: float
    log_evidence_err: float


def to_anesthetic_samples(result: InferenceResult) -> Any:
    """Convert an InferenceResult to anesthetic NestedSamples.

    For PolyChord results, returns the cached ``NestedSamples`` object
    that was read from the chain files (preserving full dead-birth
    information).  When the cache is unavailable (e.g. a result loaded
    from CSV via ``InferenceResult.from_directory``) a ``NestedSamples``
    object is reconstructed from the sample arrays.

    Parameters
    ----------
    result : InferenceResult
        Inference results (must have ``method="nested"``).

    Returns
    -------
    anesthetic.NestedSamples
        Anesthetic samples object with D_KL, d_G, logZ methods.

    Raises
    ------
    ImportError
        If anesthetic is not installed.
    ValueError
        If called on non-nested results.
    """
    try:
        from anesthetic import NestedSamples, read_chains
    except ImportError:
        msg = (
            "anesthetic is required for parameter importance analysis. "
            "Install with: pip install tidal[inference]"
        )
        raise ImportError(msg) from None

    if result.method != "nested":
        msg = (
            "Parameter importance requires nested sampling results "
            f"(got method='{result.method}'). Use --method nested."
        )
        raise ValueError(msg)

    # Check for cached anesthetic object from PolyChord
    cached = result.metadata.get("_anesthetic_samples")
    if cached is not None:
        return cached

    # Check for PolyChord chain root to read from disk
    chain_root = result.metadata.get("chain_root")
    if chain_root is not None:
        try:
            return read_chains(chain_root)
        except (FileNotFoundError, OSError):
            pass  # Fall through to manual construction

    # Reconstruct from sample arrays (e.g. loaded from CSV).
    # Sort by logL and assign births assuming constant nlive.
    import numpy as np

    n = len(result.log_likelihood)
    nlive = result.metadata.get("nlive") or 100
    if nlive >= n:
        msg = (
            f"Cannot reconstruct logl_birth: nlive ({nlive}) >= "
            f"n_samples ({n}). Need more dead points than live points."
        )
        raise ValueError(msg)
    logl_birth = np.full(n, -np.inf)
    sorted_idx = np.argsort(result.log_likelihood)
    sorted_logl = result.log_likelihood[sorted_idx]
    for i in range(nlive, n):
        logl_birth[sorted_idx[i]] = sorted_logl[i - nlive]

    return NestedSamples(
        data=result.samples,
        logL=result.log_likelihood,
        logL_birth=logl_birth,
        columns=result.param_names,
    )


def compute_parameter_importance(
    result: InferenceResult,
    n_bootstrap: int = 100,
) -> ParameterImportanceResult:
    """Compute parameter importance from nested sampling results.

    Uses anesthetic to compute the total and per-parameter KL divergence
    (information gain from prior to posterior), the Bayesian model
    dimensionality, and the evidence.

    Parameters
    ----------
    result : InferenceResult
        Nested sampling results.
    n_bootstrap : int
        Number of bootstrap samples for uncertainty estimation.

    Returns
    -------
    ParameterImportanceResult
        Parameter importance metrics.
    """
    from anesthetic import NestedSamples

    ns = to_anesthetic_samples(result)

    # --- Total statistics with bootstrap ---
    stats = ns.stats(nsamples=n_bootstrap)

    d_kl_samples = stats["D_KL"].to_numpy()
    d_g_samples = stats["d_G"].to_numpy()
    logz_samples = stats["logZ"].to_numpy()

    d_kl = float(d_kl_samples.mean())
    d_kl_err = float(d_kl_samples.std())
    d_g = float(d_g_samples.mean())
    d_g_err = float(d_g_samples.std())
    log_evidence = float(logz_samples.mean())
    log_evidence_err = float(logz_samples.std())

    # --- Per-parameter marginal D_KL ---
    # Anesthetic's `read_chains` returns integer-indexed parameter columns
    # (``[0, 1, 'logL', ...]``) in v2.0+, while the manual-reconstruction
    # fallback in to_anesthetic_samples uses named columns.  Index by
    # POSITION to work for both paths — the column order is guaranteed
    # by both the read_chains reader and the reconstruction (see #287).
    import logging

    import numpy as np

    marginal_d_kl: dict[str, float] = {}
    logl_arr = np.asarray(
        ns.logL.to_numpy() if hasattr(ns.logL, "to_numpy") else ns.logL,
    )
    logl_birth_arr = None
    if hasattr(ns, "logL_birth"):
        logl_birth_arr = np.asarray(
            ns.logL_birth.to_numpy()
            if hasattr(ns.logL_birth, "to_numpy")
            else ns.logL_birth,
        )
    for i, name in enumerate(result.param_names):
        try:
            # ns.iloc[:, i] raises IndexError if out of range; let it
            # propagate to the except below so the warning names the
            # offending column.
            col = np.asarray(ns.iloc[:, i])
            marginal = NestedSamples(
                data=col.reshape(-1, 1),
                logL=logl_arr,
                logL_birth=logl_birth_arr,
                columns=[name],
            )
            marginal_d_kl[name] = float(marginal.D_KL())
        except (ValueError, ZeroDivisionError, AttributeError, IndexError) as exc:
            logging.getLogger("tidal.inference").warning(
                "marginal D_KL failed for '%s' (col %d): %s", name, i, exc,
            )
            marginal_d_kl[name] = float("nan")

    return ParameterImportanceResult(
        param_names=result.param_names,
        d_kl=d_kl,
        d_kl_err=d_kl_err,
        d_g=d_g,
        d_g_err=d_g_err,
        marginal_d_kl=marginal_d_kl,
        log_evidence=log_evidence,
        log_evidence_err=log_evidence_err,
    )


def format_importance_table(result: ParameterImportanceResult) -> str:
    """Format parameter importance as a human-readable table.

    Parameters are ranked by marginal D_KL (most constrained first).
    """
    import math

    lines: list[str] = []
    lines.extend(("", "--- Parameter Importance (KL Divergence) ---", ""))

    # Rank by marginal D_KL (descending)
    ranked = sorted(
        result.marginal_d_kl.items(),
        key=lambda x: x[1] if math.isfinite(x[1]) else -1,
        reverse=True,
    )

    lines.extend(
        (
            f"  {'Parameter':<20} {'D_KL (nats)':>12}  {'Importance':>12}",
            f"  {'─' * 20} {'─' * 12}  {'─' * 12}",
        ),
    )

    for name, dkl in ranked:
        if not math.isfinite(dkl):
            importance = "N/A"
        elif dkl > 1.0:
            importance = "STRONG"
        elif dkl > 0.1:
            importance = "moderate"
        else:
            importance = "weak"
        lines.append(f"  {name:<20} {dkl:>12.4f}  {importance:>12}")

    lines.extend(
        (
            "",
            f"  Total D_KL:  {result.d_kl:.4f} ± {result.d_kl_err:.4f} nats",
            f"  Model dim (d_G): {result.d_g:.2f} ± {result.d_g_err:.2f}  (of {len(result.param_names)} parameters)",
            f"  log Z:  {result.log_evidence:.2f} ± {result.log_evidence_err:.2f}",
            "",
        ),
    )

    return "\n".join(lines)
