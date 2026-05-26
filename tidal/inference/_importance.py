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

    # Per-parameter marginal D_KL: histogram the posterior samples for
    # each parameter (in its prior's natural space, so the prior is
    # uniform), then compute the discrete KL against the uniform prior.
    # We deliberately don't use anesthetic's NestedSamples.D_KL() on a
    # single-column slice — that method uses only logL/weights and
    # returns the full-joint D_KL regardless of columns, so every
    # parameter would get the same value. See #308.
    marginal_d_kl: dict[str, float] = {}
    weights_arr = np.asarray(
        ns.get_weights() if hasattr(ns, "get_weights") else ns.weights,
    )
    weights_arr = weights_arr / weights_arr.sum()  # noqa: PLR6104

    # Map parameter name to its prior config, if available, so we can
    # transform log_uniform parameters into log space (where the prior
    # is uniform and the histogram-based KL is meaningful).  Accept
    # priors from either result.priors (direct call) or result.metadata
    # ["priors"] (loaded from disk).
    prior_map: dict[str, tuple[str, float, float]] = {}
    priors_iter = getattr(result, "priors", None)
    if priors_iter is None and hasattr(result, "metadata"):
        priors_iter = (result.metadata or {}).get("priors")
    if priors_iter:
        for p in priors_iter:
            if isinstance(p, dict):
                pname = p.get("name")
                pdist = p.get("distribution") or p.get("dist") or p.get("kind")
                plo = p.get("low") if "low" in p else p.get("lo")
                phi = p.get("high") if "high" in p else p.get("hi")
            else:
                pname = getattr(p, "name", None)
                pdist = (
                    getattr(p, "distribution", None)
                    or getattr(p, "dist", None)
                    or getattr(p, "kind", None)
                )
                plo = getattr(p, "low", None) or getattr(p, "lo", None)
                phi = getattr(p, "high", None) or getattr(p, "hi", None)
            if pname and pdist and plo is not None and phi is not None:
                prior_map[str(pname)] = (str(pdist), float(plo), float(phi))

    n_bins = 40
    for i, name in enumerate(result.param_names):
        try:
            col = np.asarray(ns.iloc[:, i])
            kind, lo, hi = prior_map.get(
                name, ("uniform", float(col.min()), float(col.max()))
            )
            # Transform to space where prior is uniform
            if kind == "log_uniform":
                x = np.log(col)
                a, b = float(np.log(lo)), float(np.log(hi))
            else:
                x = col
                a, b = lo, hi
            hist, edges = np.histogram(
                x,
                bins=n_bins,
                weights=weights_arr,
                range=(a, b),
            )
            p_post = hist / hist.sum() / (edges[1] - edges[0])
            q_prior = 1.0 / (b - a)
            mask = p_post > 0
            kl = float(
                np.sum(
                    p_post[mask]
                    * np.log(p_post[mask] / q_prior)
                    * (edges[1] - edges[0]),
                )
            )
            marginal_d_kl[name] = kl
        except (ValueError, ZeroDivisionError, AttributeError, IndexError) as exc:
            logging.getLogger("tidal.inference").warning(
                "marginal D_KL failed for '%s' (col %d): %s",
                name,
                i,
                exc,
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


def compute_cross_kl(
    amp_samples: Any,
    sup_samples: Any,
    params: list[str],
    n_bins: int = 80,
) -> dict[str, float]:
    r"""Per-parameter KL divergence between amp and sup posteriors.

    Computes :math:`D_{\\mathrm{KL}}(P_{+}(\\theta_i) \\| P_{-}(\\theta_i))` for
    each parameter by histogramming the 1D marginal of each posterior over a
    common binning grid spanning the union of both supports, then summing
    :math:`P_{+} \\log (P_{+}/P_{-})\\,d\\theta`.  Useful for identifying
    which parameters distinguish the amplification posterior from the
    suppression posterior — complementing the per-direction marginal D_KL
    against the prior computed by :func:`compute_parameter_importance`.

    Parameters
    ----------
    amp_samples
        Anesthetic ``Samples``/``NestedSamples`` for the amplification
        likelihood (positive ``+log A`` direction).
    sup_samples
        Anesthetic ``Samples``/``NestedSamples`` for the suppression
        likelihood (negative ``-log A`` direction).
    params
        Parameter names to compute the cross-KL for.  Each name must be a
        column in both samples objects.
    n_bins
        Number of histogram bins shared between the two posteriors.

    Returns
    -------
    dict[str, float]
        ``{param_name: cross_kl_nats}``.  NaN entries indicate a parameter
        whose histogram could not be computed (e.g. zero weight or
        degenerate support).
    """
    import numpy as np

    amp_weights = np.asarray(
        amp_samples.get_weights()
        if hasattr(amp_samples, "get_weights")
        else amp_samples.weights,
    )
    amp_weights = amp_weights / amp_weights.sum()  # noqa: PLR6104  (anesthetic weights are read-only)
    sup_weights = np.asarray(
        sup_samples.get_weights()
        if hasattr(sup_samples, "get_weights")
        else sup_samples.weights,
    )
    sup_weights = sup_weights / sup_weights.sum()  # noqa: PLR6104  (anesthetic weights are read-only)

    cross: dict[str, float] = {}
    for name in params:
        try:
            amp_col = np.asarray(amp_samples[name])
            sup_col = np.asarray(sup_samples[name])
            # Shared binning grid: union of both supports so both
            # marginals live on the same axis when we compare.
            lo = float(min(amp_col.min(), sup_col.min()))
            hi = float(max(amp_col.max(), sup_col.max()))
            if not (hi > lo):
                cross[name] = float("nan")
                continue
            edges = np.linspace(lo, hi, n_bins + 1)
            p_hist, _ = np.histogram(amp_col, bins=edges, weights=amp_weights)
            q_hist, _ = np.histogram(sup_col, bins=edges, weights=sup_weights)
            # Normalise to densities.
            dx = edges[1] - edges[0]
            p = p_hist / max(p_hist.sum(), 1e-300) / dx
            q = q_hist / max(q_hist.sum(), 1e-300) / dx
            # Smooth zero bins of q so D_KL stays finite: replace zeros
            # with the smallest non-zero density in q, scaled by 1/n_bins.
            q_floor = max(q[q > 0].min() if np.any(q > 0) else 1e-12, 1e-12) / n_bins
            q = np.where(q > 0, q, q_floor)
            mask = p > 0
            cross[name] = float(np.sum(p[mask] * np.log(p[mask] / q[mask]) * dx))
        except (ValueError, KeyError, AttributeError, IndexError) as exc:
            import logging

            logging.getLogger("tidal.inference").warning(
                "cross-KL failed for '%s': %s", name, exc
            )
            cross[name] = float("nan")

    return cross


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
