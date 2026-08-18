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

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.inference._results import InferenceResult

# ---------------------------------------------------------------------------
# Marginal D_KL estimator constants
# ---------------------------------------------------------------------------

#: Histogram bins for scalar-prior marginals.  Kept at the pre-#420 value so
#: uniform/log_uniform results are bit-identical to earlier releases.
_MARGINAL_BINS = 40

#: Bins for the empirical-reference (radial_angular) marginals.  Matches
#: :func:`compute_cross_kl`'s default resolution.
_RA_BINS = 80

#: Number of prior reference samples drawn per radial_angular record.  The
#: histogram-KL bias is ~(bins - 1) / (2 N) ≈ 4e-4 nats at this size.
_RA_REFERENCE_SAMPLES = 100_000

#: Fixed seed for the radial_angular reference draw so importance.json is
#: deterministic across re-runs of the same chain.
_RA_REFERENCE_SEED = 420


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
    consistency : dict[str, Any]
        Self-check diagnostics for the marginal estimates (#420):
        ``sum_marginals``, ``superadditivity_ok`` (for a product prior
        the chain rule gives ``D_KL(joint) >= sum of marginals`` exactly,
        so a violating sum means the estimator is broken),
        ``product_prior`` (whether the bound is exact for this run),
        ``saturated_params`` (marginals within 90% of the histogram
        ceiling ``log(n_bins)`` — resolution-limited values), and
        ``note``.  Empty dict when marginals could not be computed.
    """

    param_names: list[str]
    d_kl: float
    d_kl_err: float
    d_g: float
    d_g_err: float
    marginal_d_kl: dict[str, float]
    log_evidence: float
    log_evidence_err: float
    consistency: dict[str, Any] = field(default_factory=dict)


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


# ---------------------------------------------------------------------------
# Marginal D_KL helpers (#420)
# ---------------------------------------------------------------------------


def _identity(c: NDArray[np.float64]) -> NDArray[np.float64]:
    return c


def _uniformizing_transform(
    kind: str,
    lo: float,
    hi: float,
) -> tuple[Callable[[NDArray[np.float64]], NDArray[np.float64]], float, float] | None:
    """Return ``(fn, a, b)`` mapping samples into the space where the prior is uniform.

    ``fn(samples)`` is uniformly distributed on ``(a, b)`` when ``samples``
    are drawn from the prior — the precondition for the histogram-vs-uniform
    KL in :func:`_hist_kl_vs_uniform`.  Returns ``None`` for unknown kinds
    (or invalid hyperparameters) so the caller can fall back.

    Pre-#420 only ``log_uniform`` was transformed; ``arctan_uniform`` and
    ``normal`` fell through to identity-vs-uniform, inflating D_KL by up to
    ~2.5 nats for posteriors identical to the prior.
    """
    if kind == "uniform":
        if not hi > lo:
            return None
        return _identity, float(lo), float(hi)
    if kind == "log_uniform":
        if not (0 < lo < hi):
            return None
        return np.log, float(np.log(lo)), float(np.log(hi))
    if kind == "arctan_uniform":
        # Prior.sample/transform ignore the recorded low/high entirely:
        # theta is uniform on the FIXED eps-truncated range and x = tan(theta)
        # (see _prior.py and the follow-up issue referenced in #420).  The
        # uniformizing transform must therefore use the same fixed range —
        # using the recorded (lo, hi) here would reintroduce the bug.
        from tidal.inference._prior import _ARCTAN_EPS

        theta_hi = float(np.pi / 2 - _ARCTAN_EPS)
        return np.arctan, -theta_hi, theta_hi
    if kind == "normal":
        # low = mean, high = std (Prior's convention for normal).
        if not hi > 0:
            return None
        from scipy.special import ndtr

        mu, sigma = float(lo), float(hi)
        return (lambda c: np.asarray(ndtr((c - mu) / sigma))), 0.0, 1.0
    return None


def _hist_kl_vs_uniform(
    x: NDArray[np.float64],
    weights: NDArray[np.float64],
    a: float,
    b: float,
    n_bins: int = _MARGINAL_BINS,
) -> float:
    """Discrete KL of a weighted histogram of ``x`` against uniform on ``(a, b)``.

    This is the pre-#420 estimator body verbatim, so uniform/log_uniform
    marginals are unchanged by the #420 fix.
    """
    hist, edges = np.histogram(x, bins=n_bins, weights=weights, range=(a, b))
    p_post = hist / hist.sum() / (edges[1] - edges[0])
    q_prior = 1.0 / (b - a)
    mask = p_post > 0
    return float(
        np.sum(p_post[mask] * np.log(p_post[mask] / q_prior) * (edges[1] - edges[0])),
    )


def _ra_field(record: Any, key: str) -> Any:
    """Read a RadialAngularPrior field from a dict record or a prior object."""
    if isinstance(record, dict):
        return record[key]
    return getattr(record, key)


def _sample_radial_angular(
    record: Any,
    m: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Draw ``m`` coupling vectors from a radial_angular prior record, vectorized.

    Reproduces :meth:`tidal.inference._prior.RadialAngularPrior.transform`
    (log-uniform radius x cubed-sphere tile direction) for a whole batch of
    unit-cube points at once — the per-sample method call is too slow for
    the ~1e5-sample reference draw this module needs on every
    ``InferenceResult.save()``.  Equivalence with ``transform`` is pinned by
    a test (atol=1e-12).

    Returns an ``(m, N)`` array; column ``j`` is the prior marginal sample
    for ``record.names[j]``.
    """
    from tidal.inference._sphere import (
        face_to_axis_sign,
        tile_bounds,
    )

    names = list(_ra_field(record, "names"))
    n = len(names)
    r_lo = float(_ra_field(record, "r_lo"))
    r_hi = float(_ra_field(record, "r_hi"))
    sub_tile = tuple(int(t) for t in _ra_field(record, "sub_tile"))
    subdiv = int(_ra_field(record, "M"))
    q_arr = np.asarray(_ra_field(record, "Q"), dtype=np.float64)

    u = rng.random((m, n))
    if r_lo == r_hi:  # fixed-radius mode: u[:, 0] consumed but ignored
        r = np.full(m, r_lo)
    else:
        log_lo, log_hi = np.log(r_lo), np.log(r_hi)
        r = np.exp(log_lo + u[:, 0] * (log_hi - log_lo))

    u_lo, u_hi = tile_bounds(sub_tile, subdiv)
    u_face = u_lo + u[:, 1:] * (u_hi - u_lo)

    k, s = face_to_axis_sign(int(_ra_field(record, "face_idx")))
    v = np.empty((m, n))
    v[:, k - 1] = s
    other_axes = [i for i in range(n) if i != k - 1]
    v[:, other_axes] = u_face
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    # face_to_direction applies Q^T x per sample; row form is v @ Q.
    theta_hat = v @ q_arr
    return r[:, None] * theta_hat


def _marginal_dkl_empirical(
    col: NDArray[np.float64],
    weights: NDArray[np.float64],
    prior_col: NDArray[np.float64],
    r_lo: float,
    n_bins: int = _RA_BINS,
) -> float:
    """KL(posterior || prior) for one coupling against an empirical prior reference.

    radial_angular components have no closed-form 1-D marginal CDF, so the
    reference is a histogram of prior samples.  Both sides are binned in
    ``asinh(c / r_lo)`` space: KL is invariant under a shared monotone
    transform, and asinh spreads the log-uniform radius decades evenly
    (linear binning would pile all mass near zero over a 6-decade r range,
    collapsing every KL toward 0) while staying smooth through zero for
    sign-crossing tiles.  Zero prior bins are floored as in
    :func:`compute_cross_kl` so the estimate stays finite.
    """
    x_post = np.arcsinh(np.asarray(col, dtype=np.float64) / r_lo)
    x_prior = np.arcsinh(np.asarray(prior_col, dtype=np.float64) / r_lo)
    lo = float(min(x_post.min(), x_prior.min()))
    hi = float(max(x_post.max(), x_prior.max()))
    if not hi > lo:
        return float("nan")
    edges = np.linspace(lo, hi, n_bins + 1)
    dx = edges[1] - edges[0]
    p_hist, _ = np.histogram(x_post, bins=edges, weights=weights)
    q_hist, _ = np.histogram(x_prior, bins=edges)
    p = p_hist / max(p_hist.sum(), 1e-300) / dx
    q = q_hist / max(q_hist.sum(), 1e-300) / dx
    q_floor = max(q[q > 0].min() if np.any(q > 0) else 1e-12, 1e-12) / n_bins
    q = np.where(q > 0, q, q_floor)
    mask = p > 0
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask]) * dx))


def compute_parameter_importance(
    result: InferenceResult,
    n_bootstrap: int = 100,
) -> ParameterImportanceResult:
    """Compute parameter importance from nested sampling results.

    Uses anesthetic to compute the total and per-parameter KL divergence
    (information gain from prior to posterior), the Bayesian model
    dimensionality, and the evidence.

    Each marginal is estimated in the space where its recorded prior is
    uniform (#420); prior records come from ``result.priors`` or
    ``result.metadata["priors"]``:

    ============== =====================================================
    prior kind     uniformizing treatment
    ============== =====================================================
    uniform        identity on ``(low, high)``
    log_uniform    ``log(x)`` on ``(log low, log high)``
    arctan_uniform ``arctan(x)`` on the fixed eps-truncated theta range
                   (recorded bounds are ignored, matching Prior.sample)
    normal         Gaussian CDF ``Phi((x - mu) / sigma)`` on ``(0, 1)``
    radial_angular empirical prior reference per coupling: KL of the
                   posterior histogram against a fixed-seed prior-sample
                   histogram, binned in ``asinh(c / r_lo)`` space
    (none)         histogram vs uniform over the empirical sample range
                   — legacy fallback; a warning fires if priors metadata
                   exists but a parameter has no usable record
    ============== =====================================================

    Parameters
    ----------
    result : InferenceResult
        Nested sampling results.
    n_bootstrap : int
        Number of bootstrap samples for uncertainty estimation.

    Returns
    -------
    ParameterImportanceResult
        Parameter importance metrics, including the ``consistency``
        self-check block (superadditivity of marginals vs the joint,
        histogram-ceiling saturation — see #420).
    """
    ns = to_anesthetic_samples(result)

    # --- Total statistics with bootstrap ---
    # anesthetic's bootstrap draws from NumPy's legacy global RNG
    # (anesthetic.samples.logX uses a bare np.random.rand), so logZ and
    # D_KL error bars vary run to run on identical input.  Seed and
    # restore around it so repeated calls agree.  See issue #388.
    from tidal.inference._visualize import _deterministic_render

    with _deterministic_render():
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

    # Map each parameter name to its prior config, if available, so each
    # column can be transformed into the space where its prior is uniform
    # (#420: previously only log_uniform was transformed).  Accept priors
    # from either result.priors (direct call) or result.metadata["priors"]
    # (loaded from disk).  Joint radial_angular records cover several names
    # at once and get an empirical reference instead of a CDF transform.
    logger = logging.getLogger("tidal.inference")
    prior_map: dict[str, tuple[str, float, float]] = {}
    ra_map: dict[str, tuple[int, Any, int]] = {}  # name -> (record id, record, col j)
    priors_iter = getattr(result, "priors", None)
    if priors_iter is None and hasattr(result, "metadata"):
        priors_iter = (result.metadata or {}).get("priors")
    if priors_iter:
        for rec_idx, p in enumerate(priors_iter):
            is_joint = (
                p.get("kind") == "radial_angular"
                if isinstance(p, dict)
                else getattr(p, "names", None) is not None
            )
            if is_joint:
                for j, joint_name in enumerate(_ra_field(p, "names")):
                    ra_map[str(joint_name)] = (rec_idx, p, j)
                continue
            if isinstance(p, dict):
                pname = p.get("name")
                pdist = p.get("distribution") or p.get("dist") or p.get("kind")
                plo = p.get("low") if "low" in p else p.get("lo")
                phi = p.get("high") if "high" in p else p.get("hi")
            else:
                pname = getattr(p, "name", None)
                pdist = getattr(p, "distribution", None) or getattr(p, "dist", None)
                # Explicit None checks: `or`-chaining drops a legitimate
                # low == 0.0 (e.g. uniform:0:X, normal:0:1) into the
                # empirical fallback (#420).
                plo = getattr(p, "low", None)
                if plo is None:
                    plo = getattr(p, "lo", None)
                phi = getattr(p, "high", None)
                if phi is None:
                    phi = getattr(p, "hi", None)
            if pname and pdist and plo is not None and phi is not None:
                prior_map[str(pname)] = (str(pdist), float(plo), float(phi))

    # One reference draw per radial_angular record, shared by its names.
    ra_samples: dict[int, NDArray[np.float64]] = {}

    saturated: list[str] = []
    bias_numerator = 0.0  # accumulates (n_bins - 1) / 2 per finite marginal
    for i, name in enumerate(result.param_names):
        try:
            col = np.asarray(ns.iloc[:, i])
            if name in ra_map:
                rec_idx, record, j = ra_map[name]
                if rec_idx not in ra_samples:
                    rng = np.random.default_rng(_RA_REFERENCE_SEED)
                    ra_samples[rec_idx] = _sample_radial_angular(
                        record,
                        _RA_REFERENCE_SAMPLES,
                        rng,
                    )
                kl = _marginal_dkl_empirical(
                    col,
                    weights_arr,
                    ra_samples[rec_idx][:, j],
                    float(_ra_field(record, "r_lo")),
                )
                bins_used = _RA_BINS
            else:
                entry = prior_map.get(name)
                transform = (
                    _uniformizing_transform(*entry) if entry is not None else None
                )
                if transform is None:
                    # No usable prior record: histogram vs uniform over the
                    # empirical sample range (pre-#420 fallback).  Warn only
                    # when priors metadata exists but this parameter fell
                    # through — that is a schema mismatch, the failure mode
                    # that silently mis-scored radial_angular runs (#420).
                    if priors_iter:
                        logger.warning(
                            "marginal D_KL for '%s': no usable prior record "
                            "(kind=%s); falling back to a uniform reference "
                            "over the sample range — value is unreliable",
                            name,
                            entry[0] if entry is not None else None,
                        )
                    fn, a, b = _identity, float(col.min()), float(col.max())
                else:
                    fn, a, b = transform
                kl = _hist_kl_vs_uniform(fn(col), weights_arr, a, b)
                bins_used = _MARGINAL_BINS
            marginal_d_kl[name] = kl
            if np.isfinite(kl):
                bias_numerator += (bins_used - 1) / 2.0
                if kl > 0.9 * float(np.log(bins_used)):
                    saturated.append(name)
        except (ValueError, ZeroDivisionError, AttributeError, IndexError) as exc:
            logger.warning(
                "marginal D_KL failed for '%s' (col %d): %s",
                name,
                i,
                exc,
            )
            marginal_d_kl[name] = float("nan")

    # --- Consistency self-checks (#420) ---
    # For a product prior the chain rule gives
    #   D_KL(joint) = sum_i D_KL(marginal_i) + multi-information(posterior),
    # and multi-information >= 0, so sum(marginals) <= D_KL(joint) exactly.
    # The pre-#420 bug announced itself by violating this bound (e.g. sum
    # 7.55 vs joint 4.64); check it on every run so a broken estimator
    # self-reports instead of waiting to be noticed in a table.
    finite_marginals = [v for v in marginal_d_kl.values() if np.isfinite(v)]
    sum_marginals = float(sum(finite_marginals))
    product_prior = not ra_map  # radial_angular is not a product across names
    # Each histogram marginal is biased UP by ~(n_bins - 1) / (2 N_eff)
    # (chi-square mean of the multinomial KL estimator), and these biases
    # add across parameters while the joint (anesthetic) estimate does
    # not share them.  Allow 2x the summed expected bias plus bootstrap
    # noise before declaring the bound violated — the pre-#420 signature
    # exceeded the joint by whole nats, far beyond this allowance.
    n_eff = float(1.0 / np.sum(weights_arr**2))  # weights are normalized
    bias_allowance = 2.0 * bias_numerator / n_eff
    tolerance = 3.0 * d_kl_err + bias_allowance + 0.05
    superadditivity_ok = sum_marginals <= d_kl + tolerance
    consistency: dict[str, Any] = {
        "sum_marginals": sum_marginals,
        "joint_d_kl": d_kl,
        "superadditivity_ok": bool(superadditivity_ok),
        "product_prior": product_prior,
        "saturated_params": saturated,
        # Kish effective sample size of the posterior weights: the marginal
        # estimator's noise floor is ~(n_bins - 1)/(2 n_eff) PER PARAMETER,
        # so tiny n_eff (e.g. floor-contaminated chains) makes every
        # marginal read as spuriously informative.  Recorded so readers
        # can judge whether small marginals are signal or floor.
        "n_eff": n_eff,
        "note": (
            "sum(marginals) <= joint D_KL is exact for product priors; "
            "approximate when a radial_angular joint prior is present. "
            "Saturated marginals are within 90% of the histogram ceiling "
            "log(n_bins) and are resolution-limited."
        ),
    }
    if not superadditivity_ok:
        logger.warning(
            "marginal D_KL consistency: sum of marginals (%.3f) exceeds "
            "joint D_KL (%.3f) beyond tolerance (%.3f)%s — the marginal "
            "estimator is inconsistent with the joint; do not quote these "
            "marginals (see #420)",
            sum_marginals,
            d_kl,
            tolerance,
            "" if product_prior else " [bound approximate: joint prior present]",
        )
    if saturated:
        logger.warning(
            "marginal D_KL consistency: %s within 90%% of the estimator "
            "ceiling log(n_bins) — resolution-limited, treat as lower bounds",
            ", ".join(saturated),
        )

    return ParameterImportanceResult(
        param_names=result.param_names,
        d_kl=d_kl,
        d_kl_err=d_kl_err,
        d_g=d_g,
        d_g_err=d_g_err,
        marginal_d_kl=marginal_d_kl,
        log_evidence=log_evidence,
        log_evidence_err=log_evidence_err,
        consistency=consistency,
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
            # Normalize to densities.
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

    # Consistency self-checks (#420): make a broken estimator visible in the
    # table itself, not just in the log.
    consistency = result.consistency or {}
    if consistency and not consistency.get("superadditivity_ok", True):
        lines.extend(
            (
                f"  WARNING: sum of marginals ({consistency['sum_marginals']:.3f}) "
                f"exceeds joint D_KL ({consistency['joint_d_kl']:.3f})",
                "           — impossible for a product prior; marginals are "
                "unreliable (see #420)",
                "",
            ),
        )
    saturated = consistency.get("saturated_params") or []
    if saturated:
        lines.extend(
            (
                f"  WARNING: {', '.join(saturated)} at the histogram "
                "resolution ceiling — treat as lower bounds",
                "",
            ),
        )

    return "\n".join(lines)
