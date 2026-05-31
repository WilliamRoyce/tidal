"""Posterior visualization for inference results.

Uses **anesthetic** (Handley 2019) when available for publication-quality
corner plots that integrate natively with nested sampling output.
Falls back to matplotlib for basic visualization.

References
----------
Handley, W. (2019) "anesthetic: nested sampling visualization",
JOSS 4(37), 1414. https://doi.org/10.21105/joss.01414
"""

from __future__ import annotations

import contextlib
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    import numpy as np

    from tidal.inference._importance import ParameterImportanceResult
    from tidal.inference._results import InferenceResult


def plot_corner(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
    show_rejected_inchain: bool = True,
    show_rejected_prior: bool = False,
    full_prior_bounds: bool = False,
) -> None:
    """Generate a corner plot (2D marginals + 1D histograms).

    Uses anesthetic if available, otherwise falls back to matplotlib.

    Two independent rejected-sample overlays are available:

    - ``show_rejected_inchain``: samples in ``results.csv`` with
      ``run_status='tachyonic'`` (MC pipeline only — always empty for
      PolyChord, which discards -inf samples upfront).  Default on.
    - ``show_rejected_prior``: post-hoc uniform prior-survey samples
      from ``_rejected_prior.csv`` (5000 prior draws evaluated by the
      stability guard only).  Default off — dense and obscures the
      posterior.

    The ``full_prior_bounds`` flag (v3 architecture) overrides the
    anesthetic per-panel auto-scale with the full prior range, so that
    compactified ``arctan_uniform`` priors visibly span their full
    ±tan(89°) ≈ ±57.3 sample range.  Useful when the user needs to
    confirm posterior compactness rather than read marginal shape.
    Only supported by the anesthetic backend.
    """
    try:
        _plot_corner_anesthetic(
            result,
            output_path,
            show=show,
            show_rejected_inchain=show_rejected_inchain,
            show_rejected_prior=show_rejected_prior,
            full_prior_bounds=full_prior_bounds,
        )
    except ImportError:
        _plot_corner_matplotlib(
            result,
            output_path,
            show=show,
            show_rejected_inchain=show_rejected_inchain,
            show_rejected_prior=show_rejected_prior,
        )


def plot_trace(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Generate trace plots (parameter values vs sample index).

    Parameters
    ----------
    result : InferenceResult
        Inference results.
    output_path : Path | None
        If provided, save figure to this path.
    show : bool
        If True, display the figure interactively.
    """
    import matplotlib.pyplot as plt

    n_params = result.n_params
    fig, axes = plt.subplots(n_params + 1, 1, figsize=(10, 2.5 * (n_params + 1)))
    if n_params == 0:
        return

    axes_list = [axes] if n_params == 0 else list(axes)

    for i, name in enumerate(result.param_names):
        ax = axes_list[i]
        ax.plot(result.samples[:, i], ".", markersize=1, alpha=0.5)
        ax.set_ylabel(name)
        ax.set_xlabel("")

    # Log-likelihood trace
    ax = axes_list[-1]
    ax.plot(result.log_likelihood, ".", markersize=1, alpha=0.5, color="red")
    ax.set_ylabel("log L")
    ax.set_xlabel("Sample index")

    fig.suptitle(f"Trace plot ({result.method}, n={result.n_samples})")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def plot_importance(
    result: ParameterImportanceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
) -> None:
    """Horizontal bar chart of per-parameter KL divergence.

    Parameters ranked by marginal D_KL (most constrained first).
    Color-coded: strong (red), moderate (orange), weak (blue).

    Parameters
    ----------
    result : ParameterImportanceResult
        Output from ``compute_parameter_importance``.
    output_path : Path | None
        If provided, save figure to this path.
    show : bool
        If True, display the figure interactively.
    """
    import math

    import matplotlib.pyplot as plt
    import numpy as np

    # Rank by D_KL
    ranked = sorted(
        result.marginal_d_kl.items(),
        key=lambda x: x[1] if math.isfinite(x[1]) else -1,
    )
    names = [name for name, _ in ranked]
    values = [dkl if math.isfinite(dkl) else 0.0 for _, dkl in ranked]

    # Color by importance
    colors = []
    for v in values:
        if v > 1.0:
            colors.append("#d62728")  # red = strong
        elif v > 0.1:
            colors.append("#ff7f0e")  # orange = moderate
        else:
            colors.append("#1f77b4")  # blue = weak

    fig, ax = plt.subplots(figsize=(8, max(3, 0.6 * len(names))))
    y_pos = np.arange(len(names))
    ax.barh(y_pos, values, color=colors, edgecolor="none", height=0.6)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.set_xlabel("Marginal $D_{\\mathrm{KL}}$ (nats)")
    ax.set_title(
        f"Parameter Importance  ($d_G = {result.d_g:.1f}$ of {len(names)} params)",
    )

    # Threshold line
    ax.axvline(0.1, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.text(0.1, len(names) - 0.3, "weak", fontsize=8, color="gray", ha="left")

    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _compute_log10_amplification(
    result: InferenceResult,
    samples: object,
) -> np.ndarray | None:
    """Derive log10(A) per sample where A = P_max / P_GR_baseline.

    A is the physical amplification factor — same definition for amplify
    and suppress runs, so the corner-plot column has a consistent meaning
    across paired figures.  For the Gertsenshtein channel
    P_GR = sin²(κ B₀ t_end / 2)² depends only on the fixed simulation
    parameters (κ, B₀, t_end), so it can be computed once and shared
    across all rows.

    Returns None if the derivation can't be done (no P_max metric, no
    baseline formula in metadata, or arithmetic failure).  The caller
    falls back to plain logL in that case.
    """
    import numpy as np

    metrics = getattr(result, "metrics", None) or {}
    metric_name = (result.metadata or {}).get("likelihood_metric") or "P_max"
    sim = (result.metadata or {}).get("simulation_params") or {}
    formula = sim.get("baseline_formula")

    # Path 1 (preferred, requires P3 metadata): per-sample metric value
    # plus baseline formula → exact derivation.
    if metric_name in metrics and formula:
        fixed = dict(sim.get("fixed_params") or {})
        t_end = sim.get("t_end")
        if t_end:
            fixed["t_end"] = float(t_end)
        if fixed:
            from tidal.inference._likelihood import (
                _eval_baseline,  # pyright: ignore[reportPrivateUsage]
            )

            baseline = _eval_baseline(formula, fixed)
            if baseline is not None and baseline > 0:
                p_arr = np.asarray(metrics[metric_name], dtype=float)
                with np.errstate(divide="ignore", invalid="ignore"):
                    ratio = np.where(p_arr > 0, p_arr / baseline, np.nan)
                    return np.log10(ratio)

    # Path 2 (fallback, for chains predating the metric metadata):
    # likelihood_type tells us how logL relates to log(A):
    #   maximize → logL = log(A)        → log10(A) = logL / ln 10
    #   minimize → logL = -log(A)       → log10(A) = -logL / ln 10
    # Skip for gaussian / threshold (no A definition there).
    ltype = (result.metadata or {}).get("likelihood_type")
    if ltype not in {"maximize", "minimize"}:
        return None
    logl = np.asarray(result.log_likelihood, dtype=float)
    sign = 1.0 if ltype == "maximize" else -1.0
    with np.errstate(invalid="ignore"):
        log10a = sign * logl / np.log(10.0)
        return np.where(np.isfinite(log10a), log10a, np.nan)


def _rejected_samples_array(result: InferenceResult) -> np.ndarray | None:
    """Return an Nx(n_params) array of tachyonic-rejected samples, or None.

    Samples flagged with ``run_status='tachyonic'`` by the pre-flight
    stability guard.  Returns ``None`` if no metadata or no rejections.
    Following anesthetic conventions for visualising inaccessible regions
    of parameter space (Handley 2019 JOSS).  See also
    ``_rejected_prior.csv`` produced by Phase 4 for nested-sampling output.
    """
    import numpy as np

    metrics = getattr(result, "metrics", None)
    if not metrics or "run_status" not in metrics:
        return None
    statuses = metrics["run_status"]
    mask = np.asarray([str(s) == "tachyonic" for s in statuses])
    if not mask.any():
        return None
    return result.samples[mask]


def _clip_to_posterior_range(
    values: np.ndarray,
    weights: np.ndarray | None,
    q_low: float = 0.001,
    q_high: float = 0.999,
) -> np.ndarray:
    """Clip a derived-column array to the posterior-weighted quantile range.

    Values outside [q_low, q_high] are replaced with NaN so anesthetic's
    KDE is calibrated for the posterior scale rather than the full prior-
    exploration chain.  NaN rows are skipped by anesthetic's kde_plot_1d
    and kde_contour_plot_2d (data_compressed step).  Weights are intact.

    Default range is [0.001, 0.999] (i.e. drops 0.1% from each tail) so
    KDE bandwidth stays calibrated for the posterior peak while the clip
    boundaries fall at asymptotic density — anesthetic's 2D contours
    therefore do not produce visibly-cut horizontal lines at the panel
    boundaries.  An earlier [0.02, 0.98] default caused obvious
    horizontal cuts in derived-column rows of the corner plot because
    anesthetic drew contours that ended sharply at the 2%/98% clip
    boundaries (where density is non-negligible).

    Guard: if q_low == q_high (perfectly degenerate posterior) the original
    array is returned unchanged so anesthetic handles the edge case at its
    own scale.
    """
    import numpy as np

    finite_mask = np.isfinite(values)
    if not finite_mask.any():
        return values.copy()

    v = values[finite_mask]

    if weights is not None:
        w = np.asarray(weights, dtype=float)[finite_mask]
        total = w.sum()
        if total > 0:
            w /= total
            sort_idx = np.argsort(v)
            cdf = np.cumsum(w[sort_idx])
            v_s = v[sort_idx]
            qlo = float(v_s[np.searchsorted(cdf, q_low, side="left")])
            qhi = float(v_s[np.searchsorted(cdf, q_high, side="left")])
        else:
            qlo, qhi = float(np.quantile(v, q_low)), float(np.quantile(v, q_high))
    else:
        qlo, qhi = float(np.quantile(v, q_low)), float(np.quantile(v, q_high))

    if qlo >= qhi:
        return values.copy()

    out = values.astype(float).copy()
    out[~np.isfinite(values)] = np.nan
    out[(values < qlo) | (values > qhi)] = np.nan
    return out


def _set_derived_axis_limits(
    axes_df: object,
    plot_params: list[str],
    col_name: str,
    values: np.ndarray,
    weights: np.ndarray | None = None,
) -> None:
    """Post-hoc fix axis limits for a derived column to the posterior-weighted range.

    anesthetic's plot_2d() sets axis limits from weighted quantiles, but
    low-weight dead-point samples in nested-sampling chains can still shift
    the 5th percentile outside the posterior support, making the posterior
    appear as a thin stripe.  This resets limits to the posterior-weighted
    2.5%-97.5% (95% credible) interval +/- 5% padding so the visible
    contour boundary sits at the panel edge without large whitespace
    margins.
    """
    import numpy as np

    finite_mask = np.isfinite(values)
    if not finite_mask.any():
        return

    v = values[finite_mask]

    if weights is not None:
        w = np.asarray(weights, dtype=float)[finite_mask]
        total = w.sum()
        if total > 0:
            w /= total
            sort_idx = np.argsort(v)
            v_s = v[sort_idx]
            cdf = np.cumsum(w[sort_idx])
            q05 = float(v_s[np.searchsorted(cdf, 0.025, side="left")])
            q95 = float(v_s[np.searchsorted(cdf, 0.975, side="left")])
        else:
            q05, q95 = float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))
    else:
        q05, q95 = float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))

    spread = q95 - q05
    pad = max(0.05 * spread, 1e-3)
    lo, hi = q05 - pad, q95 + pad

    col_idx = plot_params.index(col_name)

    # 1D marginal diagonal: col_name on x-axis
    diag = _diagonal_axis(axes_df, col_idx, col_name)
    if diag is not None:
        diag.set_xlim(lo, hi)

    # Lower-triangle cells where col_name is the column (x-axis): row index > col_idx
    for i, row_name in enumerate(plot_params):
        if i <= col_idx:
            continue
        ax = _cell_axis(axes_df, i, row_name, col_idx, col_name)
        if ax is not None:
            ax.set_xlim(lo, hi)

    # Lower-triangle cells where col_name is the row (y-axis): col index < col_idx
    for j, col in enumerate(plot_params):
        if j >= col_idx:
            continue
        ax = _cell_axis(axes_df, col_idx, col_name, j, col)
        if ax is not None:
            ax.set_ylim(lo, hi)


def _set_credible_axis_limits(
    axes_df: object,
    plot_params: list[str],
    chain_data: dict[str, np.ndarray],
    weights: np.ndarray | None,
    *,
    skip: tuple[str, ...] = (),
    q_low: float = 0.025,
    q_high: float = 0.975,
    pad_frac: float = 0.05,
) -> None:
    """Tighten every panel's axes to the weighted credible interval of each parameter.

    anesthetic's default ``plot_2d`` auto-scale uses a wide quantile range
    (close to the chain extent), which makes highly-concentrated posteriors
    appear as thin lines occupying < 10% of the panel.  This helper sets
    each panel's x/y axis to the weighted [q_low, q_high] credible interval
    of the underlying parameter samples, plus a ``pad_frac`` fractional
    padding.

    Default ``[0.025, 0.975]`` (95% credible) + 5% padding fits the
    outer-contour boundary tightly without leaving large whitespace
    margins.  The visible 2D contours are drawn at 95% / 68%, so axes
    matched to the 95%-credible bounds enclose them without padding
    becoming visible whitespace beyond the contour outline.

    Skips parameters not present in ``chain_data`` (e.g. derived columns
    handled by :func:`_set_derived_axis_limits` separately).
    """
    import numpy as np

    for col_idx, col_name in enumerate(plot_params):
        if col_name not in chain_data:
            continue
        if col_name in skip:
            # log_uniform-prior params get log axes separately; using
            # additive padding here would push the lower bound negative
            # and matplotlib silently reverts the axis to linear.
            continue
        values = chain_data[col_name]
        finite_mask = np.isfinite(values)
        if not finite_mask.any():
            continue
        v = values[finite_mask]

        if weights is not None:
            w = np.asarray(weights, dtype=float)[finite_mask]
            total = w.sum()
            if total > 0:
                w /= total
                sort_idx = np.argsort(v)
                cdf = np.cumsum(w[sort_idx])
                v_s = v[sort_idx]
                lo_idx = min(
                    int(np.searchsorted(cdf, q_low, side="left")), len(v_s) - 1
                )
                hi_idx = min(
                    int(np.searchsorted(cdf, q_high, side="left")), len(v_s) - 1
                )
                qlo = float(v_s[lo_idx])
                qhi = float(v_s[hi_idx])
            else:
                qlo = float(np.quantile(v, q_low))
                qhi = float(np.quantile(v, q_high))
        else:
            qlo = float(np.quantile(v, q_low))
            qhi = float(np.quantile(v, q_high))

        spread = qhi - qlo
        if spread <= 0:
            continue
        pad = max(pad_frac * spread, 1e-9)
        lo, hi = qlo - pad, qhi + pad

        # 1D diagonal: x-axis
        diag = _diagonal_axis(axes_df, col_idx, col_name)
        if diag is not None:
            with contextlib.suppress(AttributeError, ValueError):
                diag.set_xlim(lo, hi)

        # Lower-triangle x-axis cells (row > col_idx)
        for i, row_name in enumerate(plot_params):
            if i <= col_idx:
                continue
            ax = _cell_axis(axes_df, i, row_name, col_idx, col_name)
            if ax is not None:
                with contextlib.suppress(AttributeError, ValueError):
                    ax.set_xlim(lo, hi)

        # Lower-triangle y-axis cells (col < col_idx)
        for j, col in enumerate(plot_params):
            if j >= col_idx:
                continue
            ax = _cell_axis(axes_df, col_idx, col_name, j, col)
            if ax is not None:
                with contextlib.suppress(AttributeError, ValueError):
                    ax.set_ylim(lo, hi)


def _set_log_axes_for_log_uniform_priors(
    axes_df: object,
    plot_params: list[str],
    priors: dict[str, tuple[str, float, float]],
    chain_data: dict[str, np.ndarray] | None = None,
    weights: np.ndarray | None = None,
    *,
    min_decades: float = 1.0,
) -> set[str]:
    """Set log axes for parameters with log_uniform priors *AND* multi-decade posterior support.

    A posterior over a ``log_uniform`` prior may or may not be naturally
    log-distributed in practice:

    * If the posterior spans multiple decades (e.g. mA² ∈ [0.001, 100]),
      a log axis is essential — linear squashes the bulk near zero.
    * If the posterior is concentrated within < 1 decade (e.g.
      ξ ∈ [0.5, 2.0]), log axes stretch the small-value side and
      compress the large-value side, making the structure *less* clear.

    This helper applies log scale only when the 95%-credible posterior
    range spans at least ``min_decades`` decades.  Otherwise the axis
    stays linear (anesthetic's default), letting the credible-interval
    axis tightening take over.

    For each qualifying parameter:

    * Sets the diagonal panel's x-axis to log scale.
    * Sets each lower-triangle cell containing this parameter on either
      the x- or y-axis to the corresponding log scale.

    Returns the set of parameter names that received log-scale axes,
    so the caller can skip them in subsequent linear-credible-interval
    tightening (which would otherwise compute a negative-lo bound and
    cause matplotlib to revert the scale to linear).
    """
    import numpy as np

    logged: set[str] = set()
    for col_idx, col_name in enumerate(plot_params):
        if col_name not in priors:
            continue
        dist, _lo, _hi = priors[col_name]
        if dist != "log_uniform":
            continue

        # Decade-span check: skip log scale if posterior is concentrated
        # within < min_decades decades of dynamic range.
        if chain_data is not None and col_name in chain_data:
            values = chain_data[col_name]
            finite_mask = np.isfinite(values) & (values > 0)
            if finite_mask.any():
                v = values[finite_mask]
                if weights is not None:
                    w = np.asarray(weights, dtype=float)[finite_mask]
                    total = w.sum()
                    if total > 0:
                        w /= total
                        sort_idx = np.argsort(v)
                        cdf = np.cumsum(w[sort_idx])
                        v_s = v[sort_idx]
                        qlo_idx = min(int(np.searchsorted(cdf, 0.025)), len(v_s) - 1)
                        qhi_idx = min(int(np.searchsorted(cdf, 0.975)), len(v_s) - 1)
                        qlo, qhi = float(v_s[qlo_idx]), float(v_s[qhi_idx])
                    else:
                        qlo, qhi = (
                            float(np.quantile(v, 0.025)),
                            float(np.quantile(v, 0.975)),
                        )
                else:
                    qlo, qhi = (
                        float(np.quantile(v, 0.025)),
                        float(np.quantile(v, 0.975)),
                    )
                if qlo > 0 and qhi > qlo:
                    decades = math.log10(qhi / qlo)
                    if decades < min_decades:
                        # Posterior is concentrated within < min_decades —
                        # log axis would stretch things unhelpfully.
                        continue

        # 1D diagonal: x-axis only
        diag = _diagonal_axis(axes_df, col_idx, col_name)
        if diag is not None:
            with contextlib.suppress(AttributeError, ValueError):
                diag.set_xscale("log")

        # Lower-triangle x-axis cells (row > col_idx)
        for i, row_name in enumerate(plot_params):
            if i <= col_idx:
                continue
            ax = _cell_axis(axes_df, i, row_name, col_idx, col_name)
            if ax is not None:
                with contextlib.suppress(AttributeError, ValueError):
                    ax.set_xscale("log")

        # Lower-triangle y-axis cells (col < col_idx)
        for j, col in enumerate(plot_params):
            if j >= col_idx:
                continue
            ax = _cell_axis(axes_df, col_idx, col_name, j, col)
            if ax is not None:
                with contextlib.suppress(AttributeError, ValueError):
                    ax.set_yscale("log")

        logged.add(col_name)

    return logged


def _set_full_prior_axis_limits(
    axes_df: object,
    plot_params: list[str],
    priors: dict[str, tuple[str, float, float]],
) -> None:
    """Set every panel's x/y limits to the *full prior range* per parameter.

    Counter-part to :func:`_set_derived_axis_limits` which uses the
    posterior-weighted 5%-95% interval.  This mode is invoked by
    ``tidal plot --type corner --full-prior-bounds`` to make the
    compactified-prior coverage visually honest: arctan_uniform:-89:89
    panels show their full ±tan(89°) ≈ ±57.3 range so the user can
    confirm the posterior is genuinely concentrated rather than
    artificially auto-scaled to a sub-region.

    For ``arctan_uniform`` priors the bounds are degrees of the angle
    variable and the physical sample range is ``tan(deg · π/180)``.
    For ``log_uniform`` and ``uniform`` priors the ``[low, high]``
    bounds map directly to the sample axis.
    """
    import math

    def _physical_bounds(dist: str, lo: float, hi: float) -> tuple[float, float]:
        if dist == "arctan_uniform":
            # Map angle-degrees bounds → physical sample bounds via tan.
            return math.tan(math.radians(lo)), math.tan(math.radians(hi))
        return lo, hi

    for col_idx, col_name in enumerate(plot_params):
        if col_name not in priors:
            continue
        dist, lo, hi = priors[col_name]
        x_lo, x_hi = _physical_bounds(dist, lo, hi)
        if not (x_lo < x_hi):
            continue

        # 1D marginal diagonal
        diag = _diagonal_axis(axes_df, col_idx, col_name)
        if diag is not None:
            diag.set_xlim(x_lo, x_hi)

        # Lower-triangle x-axis cells (row > col_idx)
        for i, row_name in enumerate(plot_params):
            if i <= col_idx:
                continue
            ax = _cell_axis(axes_df, i, row_name, col_idx, col_name)
            if ax is not None:
                ax.set_xlim(x_lo, x_hi)

        # Lower-triangle y-axis cells (col < col_idx)
        for j, col in enumerate(plot_params):
            if j >= col_idx:
                continue
            ax = _cell_axis(axes_df, col_idx, col_name, j, col)
            if ax is not None:
                ax.set_ylim(x_lo, x_hi)


def _render_anesthetic_corner_into(
    samples: object,
    param_names: list[str],
    *,
    target_fig: object | None = None,
    labels: dict[str, str] | None = None,
    ticks: str | None = "inner",
    show_diagonal: bool = True,
    fill_two_tone: bool = True,
) -> tuple[object, object]:
    """Render a two-tone anesthetic corner plot.

    Used both by :func:`_plot_corner_anesthetic` (standalone full figure;
    ``target_fig=None`` — anesthetic's ``plot_2d`` builds its own
    Figure, preserving the standalone path bit-identically) and by the
    cubed-sphere atlas (each face panel passes its matplotlib
    :class:`~matplotlib.figure.SubFigure`; the helper then explicitly
    builds axes via :func:`anesthetic.make_2d_axes` so they live inside
    the SubFigure rather than in a fresh top-level Figure).

    The function does ONLY the minimal common core: invoke anesthetic
    plotting (with the degenerate-posterior fallbacks) and apply the
    two-tone Planck-blue solid fill patch.  Standalone-specific
    decoration (prior overlay, MAP marker, CI band, log axes, axis
    limits, rejected overlay, legend, suptitle) is layered on by the
    caller.

    Parameters
    ----------
    samples : anesthetic.Samples-like
        A NestedSamples / MCMCSamples object exposing ``plot_2d``.
    param_names : list of str
        Column names in ``samples`` to plot, in the desired corner-plot
        order.
    target_fig : Figure or SubFigure, optional
        A matplotlib Figure / SubFigure to render the axes into.  None
        means the standalone path — ``samples.plot_2d(param_names)``
        builds its own figure (bit-identical to the pre-refactor
        behaviour).  When provided, axes are built explicitly via
        :func:`anesthetic.make_2d_axes` inside ``target_fig`` (passed
        through anesthetic's ``fig=`` fig_kw); the axes fill the entire
        SubFigure.
    labels : dict, optional
        Mapping ``param_name -> LaTeX label`` for axis labels.  Only
        consulted on the composable (target_fig) path; the standalone
        path uses anesthetic's defaults.
    ticks : {"inner", "outer", None}, optional
        Anesthetic ``ticks`` mode for the composable path.  ``None``
        strips tick marks entirely (the atlas's setting).  Ignored on
        the standalone path (anesthetic uses "inner" by default).
    show_diagonal : bool, optional
        Whether the diagonal 1D KDE marginals are drawn.  True by
        default.  Only honoured on the composable path.
    fill_two_tone : bool, optional
        Apply :func:`_force_solid_credible_fills` after rendering to
        replace anesthetic's per-panel-normalized gradient with the
        cosmology-standard solid Planck-blue two-tone palette.

    Returns
    -------
    fig : Figure
        The matplotlib Figure / SubFigure containing the rendered axes.
    axes_df : AxesDataFrame
        The anesthetic axes DataFrame, suitable for further overlay
        operations (priors, MAP markers, log axes, etc.).
    """
    import numpy as np

    def _plot_with_fallback(axes_arg: object) -> object:
        """Run ``samples.plot_2d`` with the degenerate-posterior fallbacks.

        Three failure modes from anesthetic's KDE machinery on near-flat /
        lower-dimensional posteriors:
          - qhull Delaunay triangulation: 'singular input data'
          - scipy gaussian_kde: 'singular data covariance matrix'
          - matplotlib tricontour: 'Triangulation is invalid' (samples too
            concentrated for matplotlib's Delaunay; observed on Phase E
            atlas T5.1/T5.2 corner panels where one chain hit a sharp peak)
        All fall back to KDE diagonals + scatter cross-panels.
        """
        try:
            return samples.plot_2d(axes_arg, label="68% / 95% CL")
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as e:
            msg = str(e).lower()
            if (
                "qhull" not in msg
                and "singular" not in msg
                and "lower-dimensional" not in msg
                and "triangulation" not in msg
            ):
                raise
            try:
                return samples.plot_2d(
                    axes_arg,
                    label="samples",
                    kinds={
                        "diagonal": "kde_1d",
                        "lower": "scatter_2d",
                        "upper": "scatter_2d",
                    },
                )
            except (RuntimeError, ValueError, np.linalg.LinAlgError):
                return samples.plot_2d(
                    axes_arg,
                    label="samples",
                    kinds={
                        "diagonal": "hist_1d",
                        "lower": "scatter_2d",
                        "upper": "scatter_2d",
                    },
                )

    fig: object
    axes_df: object
    if target_fig is not None:
        # Composable path: axes live inside the caller-provided figure
        # (typically a SubFigure of an outer atlas grid).  Passing
        # ``fig=target_fig`` to ``make_2d_axes`` short-circuits the
        # ``plt.figure(**fig_kw)`` fallback at anesthetic plot.py:748,
        # so the new axes are added to ``target_fig`` rather than a
        # fresh hidden top-level Figure.
        from anesthetic import make_2d_axes

        fig, axes_df = make_2d_axes(
            param_names,
            labels=labels,
            ticks=ticks,
            lower=True,
            diagonal=show_diagonal,
            upper=False,
            fig=target_fig,
        )
        _plot_with_fallback(axes_df)
        if labels is not None:
            # ``plot_2d`` overwrites the LaTeX labels set by
            # ``make_2d_axes(labels=...)`` with the raw column names
            # (anesthetic plot.py: AxesDataFrame.__init__ calls set_labels
            # but plot_2d re-applies index/columns as labels).  Re-apply
            # via the public ``set_labels`` method post-render.
            axes_df.set_labels(labels)  # type: ignore[attr-defined]
    else:
        # Standalone path: bit-identical to the pre-refactor behaviour —
        # let anesthetic.plot_2d build its own figure and axes from the
        # parameter list.  Extract fig + axes_df from anesthetic's
        # return value: either ``(fig, axes_df)`` (older versions) or
        # just the AxesDataFrame (>= 2.0).
        ret = _plot_with_fallback(param_names)
        if isinstance(ret, tuple) and len(ret) == 2:
            fig, axes_df = ret[0], ret[1]
        else:
            axes_df = ret
            cell = ret.iloc[0, 0] if hasattr(ret, "iloc") else next(iter(ret))
            fig = cell.get_figure()

    if fill_two_tone:
        _force_solid_credible_fills(fig)
    return fig, axes_df


def _plot_corner_anesthetic(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
    show_rejected_inchain: bool = True,
    show_rejected_prior: bool = False,
    full_prior_bounds: bool = False,
) -> None:
    """Corner plot using anesthetic (Handley 2019).

    Layout follows the cosmology corner-plot convention:

    - 2D credible regions at 68% / 95% on every off-diagonal cell
      (anesthetic default ``levels=[0.95, 0.68]`` — same credible mass
      across all panels)
    - 1D marginals on the diagonal with a red prior overlay (see
      :func:`_overlay_priors`), the MAP marker, and a 68% credible-
      interval shaded band
    - ``logL`` appended as a final corner row/column so the distribution
      of log-amplification is shown directly in parameter-space context;
      this reuses anesthetic's weight-aware KDE rather than building a
      separate (and weight-unaware) histogram
    - Tachyonic-rejected samples overlaid in red on the upper triangle
      (see :func:`_overlay_rejected_anesthetic`)
    - Headline numbers (log Z, max A, D_KL) printed in the suptitle
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np

    from tidal.inference._importance import to_anesthetic_samples

    try:
        samples = to_anesthetic_samples(result)
    except ValueError:
        # Fall back to MCMCSamples for non-nested results
        from anesthetic import MCMCSamples

        samples = MCMCSamples(
            data=result.samples,
            columns=result.param_names,
        )

    # Corner plot shows the SAMPLED parameters only — no derived columns.
    # Per supervisor direction (2026-05-29): the prior log10(A) overlay was
    # confusing and non-standard; remove it from all corner plots.
    plot_params: list[str] = list(result.param_names)

    # Core render: make_2d_axes + plot_2d + two-tone fill patch.  The
    # standalone path passes ``subplot_spec=None`` so anesthetic builds
    # its own Figure; the atlas variant passes a SubplotSpec for one
    # panel of the outer 2 x N grid.  Both paths share the same KDE +
    # fallback + fill logic.
    fig, axes_df = _render_anesthetic_corner_into(
        samples,
        plot_params,
        target_fig=None,
        labels=None,
        ticks="inner",
        show_diagonal=True,
        fill_two_tone=True,
    )
    # Set log axes for log_uniform-prior parameters FIRST — before
    # prior/MAP overlays — so the overlays are drawn in the correct
    # log space, not against linear axes that get re-scaled later.
    # Only params whose posterior spans ≥ 1 decade get log scale; tightly
    # concentrated posteriors (< 1 decade) stay linear so log-stretching
    # doesn't make them harder to read.
    prior_map = _extract_prior_map(result)
    chain_data_for_log = {
        name: result.samples[:, i] for i, name in enumerate(result.param_names)
    }
    logged_params = _set_log_axes_for_log_uniform_priors(
        axes_df,
        plot_params,
        prior_map,
        chain_data_for_log,
        result.weights,
    )
    _overlay_priors(axes_df, result)
    _overlay_map_and_ci(axes_df, result)
    # v3 architecture: the upper-right triangle is hidden because at large N
    # (v3 chains have higher ESS than v2) anesthetic's default scatter
    # rendering becomes a solid mess with no extra information beyond the
    # lower-triangle KDE contours.  See docs/V3_ARCHITECTURE.md and
    # https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/347.
    _hide_upper_triangle(axes_df, plot_params)
    if full_prior_bounds:
        _set_full_prior_axis_limits(
            axes_df,
            plot_params,
            prior_map,
        )
    else:
        # Default: tighten per-panel axes to the weighted 95%-credible
        # interval of each parameter, plus 5% padding.  anesthetic's
        # default auto-scale is closer to the full chain extent, which
        # makes highly-concentrated posteriors (e.g. Stage A sup
        # deltam 98%-cred [-1.4, +1.3] vs chain extent [-19, +20])
        # appear as thin lines occupying < 10% of the panel.
        # Skip log-scale params (those whose posterior spanned ≥ 1
        # decade and got log axes) since additive-padding around their
        # lo-bound goes negative and matplotlib reverts to linear.
        _set_credible_axis_limits(
            axes_df,
            plot_params,
            chain_data_for_log,
            result.weights,
            skip=tuple(logged_params),
        )
        # log10_A handling removed per supervisor direction (2026-05-29):
        # derived amplification column was confusing and non-standard.
    # Re-apply log scale at the very end — overlays (priors, MAP, CI)
    # internally save/restore xlim around their plot() calls which on
    # some matplotlib versions silently reverts the scale to linear.
    # This final pass guarantees the rendered axis is log for every
    # log_uniform-prior parameter, regardless of intermediate changes.
    _set_log_axes_for_log_uniform_priors(
        axes_df,
        plot_params,
        prior_map,
        chain_data_for_log,
        result.weights,
    )
    has_any_rejected = _overlay_rejected_anesthetic(
        axes_df,
        result,
        show_inchain=show_rejected_inchain,
        show_prior=show_rejected_prior,
    )

    # Legend in the conventional top-right slot (the empty upper-corner
    # cell of the triangle).  Combine anesthetic's auto-added contour
    # patch with manual entries for MAP, prior, and rejected overlay.
    _add_corner_legend(axes_df, plot_params, has_rejected=has_any_rejected)

    # Suptitle with headline numbers.
    title_lines = [
        f"Posterior ({result.method}, n={result.n_samples}, "
        f"ESS={result.effective_sample_size():.0f})"
    ]
    headline_parts: list[str] = []
    if result.log_evidence is not None:
        if result.log_evidence_err is not None:
            headline_parts.append(
                f"log Z = {result.log_evidence:.3f} ± {result.log_evidence_err:.3f}"
            )
        else:
            headline_parts.append(f"log Z = {result.log_evidence:.3f}")
    # Headline: report mode-aware max(exp(logL)) headline as fallback;
    # the log10_A-based summary was removed per supervisor direction.
    if False:
        pass
    else:
        finite_logl = result.log_likelihood[np.isfinite(result.log_likelihood)]
        if finite_logl.size > 0:
            ltype = (result.metadata or {}).get("likelihood_type")
            max_factor = float(np.exp(np.max(finite_logl)))
            labels = {
                "maximize": "max A",
                "minimize": "max suppression P_GR/P",
                "extremize": "max |log A|",
            }
            label = labels.get(ltype) if ltype else "max exp(log L)"
            headline_parts.append(f"{label} = {max_factor:.3f}")
    pi = (result.metadata or {}).get("parameter_importance") or {}
    if isinstance(pi, dict) and "d_kl" in pi:
        headline_parts.append(f"D_KL = {pi['d_kl']:.3f} nats")
    if headline_parts:
        title_lines.append("   |   ".join(headline_parts))
    fig.suptitle("\n".join(title_lines), y=1.02)
    # Suppress an unused-import warning for mpatches if no rejected samples
    # forced the legend entry path; keep the import top-level for clarity.
    _ = mpatches

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _force_solid_credible_fills(fig: object) -> None:
    """Render solid-fill credible regions per Planck/DES/getdist convention.

    Replaces anesthetic's per-panel-normalized gradient fills with solid
    colours per credible level.

    **Why**: anesthetic's default ``kde_contour_plot_2d`` renders contour
    fills via ``ax.contourf(P, levels=[c95, c68], cmap=cmap, vmin=0,
    vmax=P.max())``.  Because ``vmax`` is the *per-panel* peak density,
    the same iso-probability levels (68%, 95%) land at different cmap
    positions in different panels, producing visually inconsistent fill
    colours across the corner plot.  See
    ``.venv/.../anesthetic/plot.py:1311``.

    The cosmology community standard (Planck 2018, DES Y3, ACT DR6,
    SPT-3G, Euclid forecasts; produced via ``getdist`` with ``filled=
    True``) uses **solid two-tone fills** — one colour for the 95%
    credible region, a darker colour for the 68% region — identical
    across every panel.  This is what we want: a reader extracting
    credible-region boundaries is not misled by panel-to-panel colour
    drift.

    **How**: walk each axis's ``ContourSet`` collections (anesthetic
    creates exactly one per panel for the 2D KDE fill), and override
    the face colours to a fixed two-tone Planck-blue palette.  The
    associated contour *line* collections (also drawn by anesthetic)
    are left untouched so the boundaries remain crisp.

    This post-hoc patch keeps anesthetic's KDE computation and contour
    *placement* intact; we only override the rendered fill colours.
    """
    import matplotlib.contour as mcontour
    from matplotlib import colors as mcolors

    # Two-tone palette: light blue for the 95% credible band, darker for
    # the 68%.  Matches the standard Planck-blue scheme via getdist.
    # anesthetic uses `iso_probability_contours([0.95, 0.68])` which
    # returns levels in *increasing* density order [c95, c68], producing
    # `contourf` bands [c95, c68] (95% band: between c95 and c68 — the
    # "outer" annular ring) and [c68, max] (68% core).  So palette index
    # 0 is the 95% band (light), index 1 is the 68% core (dark).
    fill_palette = ["#aac8e9", "#3877b8"]
    if not hasattr(fig, "axes"):
        return
    for ax in fig.axes:
        contour_sets = [
            obj
            for obj in ax.findobj(match=mcontour.ContourSet)
            if getattr(obj, "filled", False)
        ]
        for cs in contour_sets:
            _apply_solid_fill(cs, fill_palette, mcolors)


def _apply_solid_fill(
    cs: object,
    palette: list[str],
    mcolors: object,
) -> None:
    """Override one ContourSet's fill colours with a solid palette.

    Handles both old matplotlib (< 3.8: ``.collections`` attribute holds
    one PathCollection per band) and new matplotlib (≥ 3.8: ContourSet
    *is* a single Collection, ``.set_facecolor`` takes one entry per band).
    """

    def colour_for(j: int) -> tuple[float, float, float, float]:
        return mcolors.to_rgba(  # type: ignore[attr-defined]
            palette[min(j, len(palette) - 1)],
            alpha=0.85,
        )

    try:
        colls = list(getattr(cs, "collections", []) or [])
        if colls:
            for j, col in enumerate(colls):
                if hasattr(col, "set_facecolor"):
                    col.set_facecolor(colour_for(j))
                    col.set_edgecolor("none")
            return
        n_bands = max(len(getattr(cs, "levels", [])) - 1, 1)
        cs.set_facecolor([colour_for(j) for j in range(n_bands)])  # type: ignore[attr-defined]
        cs.set_edgecolor("none")  # type: ignore[attr-defined]
    except (AttributeError, TypeError):
        return


def _overlay_map_and_ci(axes_df: object, result: InferenceResult) -> None:
    """Mark the MAP estimate and 68% credible interval on each 1D marginal.

    Anesthetic does not provide a built-in MAP/CI annotation
    (Handley 2019 leaves these as user-side decisions, see
    `anesthetic/plot.py:867` `kde_plot_1d`).  This helper draws:

    - solid black vertical line at the MAP value
    - 68% credible interval as a translucent black axvspan

    The existing red prior overlay (:func:`_overlay_priors`) stays
    independent.
    """
    try:
        map_estimate = result.best()
    except Exception:  # noqa: BLE001
        return
    try:
        ci_68 = result.credible_interval(0.68)
    except Exception:  # noqa: BLE001
        ci_68 = {}

    for i, name in enumerate(result.param_names):
        ax = _diagonal_axis(axes_df, i, name)
        if ax is None:
            continue
        # Anesthetic's diagonal cells share an axis with a "twin" — KDE is
        # drawn on the twin.  Annotations look correct on the cell axis
        # itself (the visible one).  Pick whichever has data.
        target = ax
        if hasattr(ax, "twin") and ax.twin is not None:
            target = ax.twin
        try:
            xmap = float(map_estimate.get(name, float("nan")))
            if not math.isnan(xmap):
                target.axvline(xmap, color="black", lw=1.0, alpha=0.85, zorder=3)
        except (TypeError, ValueError):
            pass
        try:
            if name in ci_68:
                lo, hi = ci_68[name]
                target.axvspan(
                    float(lo),
                    float(hi),
                    alpha=0.12,
                    color="black",
                    lw=0,
                    zorder=0,
                )
        except (TypeError, ValueError):
            pass


def _add_corner_legend(
    axes_df: object,
    plot_params: list[str],
    *,
    has_rejected: bool,
) -> None:
    """Render a corner-plot legend in the conventional top-right slot.

    Combines anesthetic's auto-added 68%/95% CL Rectangle patch
    (created when ``label=`` is passed to ``plot_2d``) with manual
    entries for MAP, prior, and (optionally) the rejected overlay.
    """
    import matplotlib.patches as mpatches
    from matplotlib.lines import Line2D

    if not plot_params:
        return

    # Top-right cell of the triangle is the conventional legend slot.
    target = None
    if hasattr(axes_df, "loc"):
        try:
            target = axes_df.loc[plot_params[0], plot_params[-1]]
        except (KeyError, AttributeError):
            target = None
    if target is None and hasattr(axes_df, "iloc"):
        try:
            target = axes_df.iloc[0, -1]
        except (IndexError, AttributeError):
            target = None
    if target is None:
        return

    # Pull anesthetic's auto-added contour patches (one per cell where
    # plot_2d ran) by collecting all labelled handles from the figure.
    fig = target.get_figure() if hasattr(target, "get_figure") else None
    handles: list[object] = []
    seen_labels: set[str] = set()
    if fig is not None:
        for ax in fig.axes:
            for h, lbl in zip(*ax.get_legend_handles_labels(), strict=False):
                if lbl and lbl not in seen_labels:
                    handles.append(h)
                    seen_labels.add(lbl)

    # Manual MAP and prior entries
    if "MAP" not in seen_labels:
        handles.append(Line2D([], [], color="black", lw=1.0, label="MAP + 68% CI"))
    if "prior" not in seen_labels:
        handles.append(Line2D([], [], color="red", lw=1.2, label="prior"))
    if has_rejected and "tachyonic" not in seen_labels:
        handles.append(
            mpatches.Patch(color="C3", alpha=0.7, label="Unstable Re(λ) > 0")
        )

    if not handles:
        return
    target.legend(
        handles=handles,
        loc="upper right",
        fontsize=8,
        frameon=True,
        framealpha=0.9,
    )


def _hide_upper_triangle(axes_df: object, plot_params: list[str]) -> None:
    """Hide the upper-right triangle of an anesthetic corner plot.

    v3 default: only the lower triangle (KDE 2D contours) and the diagonal
    (1D marginals) are shown.  The upper triangle is set invisible so the
    scatter rendering doesn't clutter the figure at large N.

    Iterates the off-diagonal cells with ``j > i`` (upper triangle) and
    calls ``ax.set_visible(False)`` on each.  ``_cell_axis`` already
    handles the AxesDataFrame indexing convention and returns ``None``
    for missing cells (defensive).  Tolerates anesthetic API variations.
    """
    for i, ni in enumerate(plot_params):
        for j, nj in enumerate(plot_params):
            if j <= i:
                continue
            ax = _cell_axis(axes_df, i, ni, j, nj)
            if ax is None:
                continue
            try:
                ax.set_visible(False)
            except (AttributeError, ValueError):
                # Defensive: unfamiliar axes object (e.g.\ a subfigure
                # wrapper); skip rather than fail the whole plot.
                continue


def _overlay_rejected_anesthetic(
    axes_df: object,
    result: InferenceResult,
    *,
    show_inchain: bool = True,
    show_prior: bool = False,
) -> bool:
    """Overlay tachyonic-rejected samples on the corner plot.

    Adds a red scatter on upper-triangle panels of an anesthetic
    AxesDataFrame.  The two overlay layers are independently toggleable —
    the prior survey contains thousands of points and obscures the
    posterior contours when rendered, so it defaults off.

    Returns True iff at least one layer rendered samples — used by the
    caller to decide whether to add a legend entry.
    """
    rejected = _rejected_samples_array(result) if show_inchain else None
    rejected_prior = _load_rejected_prior_overlay(result) if show_prior else None
    if rejected is None and rejected_prior is None:
        return False

    # Upper triangle (j > i): anesthetic's scatter convention.  Lower
    # triangle (j < i): KDE contours — overlaying scatter there would
    # clash with the contour fill, so we leave it alone.
    names = result.param_names
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            if j <= i:
                continue
            ax = _cell_axis(axes_df, i, ni, j, nj)
            if ax is None:
                continue
            if rejected is not None:
                ax.scatter(
                    rejected[:, j],
                    rejected[:, i],
                    s=6,
                    alpha=0.7,
                    color="C3",
                    label="tachyonic",
                    zorder=2,
                )
            if rejected_prior is not None:
                ax.scatter(
                    rejected_prior[:, j],
                    rejected_prior[:, i],
                    s=3,
                    alpha=0.3,
                    color="C3",
                    marker=".",
                    zorder=1,
                )
    return True


def _cell_axis(axes_df: object, i: int, ni: str, j: int, nj: str) -> object | None:
    """Retrieve an off-diagonal (i, j) cell axis.

    Works for an anesthetic AxesDataFrame or a 2D array of axes, falling
    back through name and integer indexing.
    """
    if hasattr(axes_df, "loc"):
        try:
            return axes_df.loc[ni, nj]
        except (KeyError, AttributeError):
            pass
    if hasattr(axes_df, "iloc"):
        try:
            return axes_df.iloc[i, j]
        except (IndexError, AttributeError):
            pass
    try:
        return axes_df[i, j]  # type: ignore[index]
    except (IndexError, TypeError):
        return None


def _load_rejected_prior_overlay(result: InferenceResult) -> np.ndarray | None:
    """Read the post-hoc prior-only stability sweep file if present.

    Phase 4 writes ``_rejected_prior.csv`` in the inference output dir;
    ``InferenceResult.metadata['rejected_prior_path']`` carries the path
    when populated.  Returns an ndarray of rejected prior samples in the
    same column order as ``result.param_names``, or ``None`` if absent.
    """
    import numpy as np

    meta = getattr(result, "metadata", None) or {}
    from pathlib import Path as _PathType

    # Try (1) the explicit path stored in metadata, then (2) a sibling
    # _rejected_prior.csv next to inference.json (typical when results have
    # been pulled from HPC and the original absolute path is unreachable).
    candidates: list[_PathType] = []
    abs_path = meta.get("rejected_prior_path")
    if abs_path:
        candidates.append(_PathType(abs_path))
    loaded_from = meta.get("loaded_from")
    if loaded_from:
        candidates.append(_PathType(loaded_from) / "_rejected_prior.csv")

    p: _PathType | None = next((c for c in candidates if c.exists()), None)
    if p is None:
        return None
    try:
        import csv

        with p.open("r", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            cols = result.param_names
            rows: list[list[float]] = []
            for r in reader:
                if r.get("run_status") != "tachyonic":
                    continue
                try:
                    rows.append([float(r[c]) for c in cols])
                except (KeyError, ValueError):
                    continue
        if not rows:
            return None
        return np.asarray(rows, dtype=float)
    except OSError:
        return None


def _overlay_priors(axes_df: object, result: InferenceResult) -> None:
    """Overlay prior density on each 1D marginal; log-scale for log_uniform.

    Accepts either an anesthetic AxesDataFrame (row/col indexed by
    parameter name) or an ordinary 2D numpy array of matplotlib Axes.
    Diagonal cells hold the 1D marginal posteriors.

    No-op if no priors are stored in ``result.metadata["priors"]``.
    """
    priors = _extract_prior_map(result)
    if not priors:
        return

    for i, name in enumerate(result.param_names):
        if name not in priors:
            continue
        ax = _diagonal_axis(axes_df, i, name)
        if ax is None:
            continue
        dist, lo, hi = priors[name]
        # Overlay the prior density in the SAME space as the posterior
        # that anesthetic drew (linear x-axis). For log_uniform this is
        # the 1/(x·ln(hi/lo)) curve; for uniform it is a flat horizontal
        # line. The visual "null" signature is then "blue posterior
        # curve tracks the red prior curve" — no more mistaking a
        # log_uniform prior shape for a physical signal (#309).
        import math

        import numpy as np

        # Preserve the original axis limits so the log_uniform prior's
        # large density at x→lo doesn't squash the posterior via
        # autoscale.
        xlim = ax.get_xlim() if hasattr(ax, "get_xlim") else None
        ylim = ax.get_ylim() if hasattr(ax, "get_ylim") else None

        if dist == "log_uniform" and lo > 0 and hi > lo:
            xs = np.geomspace(lo, hi, 200)
            density = 1.0 / (xs * math.log(hi / lo))
            ax.plot(
                xs,
                density,
                color="red",
                lw=1.2,
                alpha=0.8,
                label="prior",
            )
        elif dist == "uniform" and hi > lo:
            ax.axhline(
                1.0 / (hi - lo),
                color="red",
                lw=1.2,
                alpha=0.8,
                label="prior",
            )

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)


def _extract_prior_map(
    result: InferenceResult,
) -> dict[str, tuple[str, float, float]]:
    """Pull (distribution, low, high) for each parameter from the result.

    Priors land in metadata via ``tidal/cli/_sample.py`` (commit that
    added ``result.metadata["priors"]``). When replotting old chains
    that pre-date that commit, use ``tidal plot --priors "..."`` to
    inject them via :func:`tidal.inference._prior.parse_prior`.
    """
    meta = getattr(result, "metadata", None) or {}
    entries = meta.get("priors") or []
    out: dict[str, tuple[str, float, float]] = {}
    for p in entries:
        if isinstance(p, dict):
            name = p.get("name")
            dist = p.get("distribution") or p.get("dist") or p.get("kind")
            low = p.get("low", p.get("lo"))
            high = p.get("high", p.get("hi"))
        else:
            name = getattr(p, "name", None)
            dist = (
                getattr(p, "distribution", None)
                or getattr(p, "dist", None)
                or getattr(p, "kind", None)
            )
            low = getattr(p, "low", None) or getattr(p, "lo", None)
            high = getattr(p, "high", None) or getattr(p, "hi", None)
        if name and dist and low is not None and high is not None:
            out[str(name)] = (str(dist), float(low), float(high))
    return out


def _diagonal_axis(axes_df: object, i: int, name: str) -> object | None:
    """Retrieve the 1D marginal axis for parameter ``name``.

    Works for an anesthetic AxesDataFrame (indexed by parameter name)
    or a plain 2D ndarray-like of matplotlib Axes.
    """
    # AxesDataFrame: row and column labels are parameter names
    if hasattr(axes_df, "loc"):
        try:
            return axes_df.loc[name, name]
        except (KeyError, AttributeError):
            pass
    if hasattr(axes_df, "iloc"):
        try:
            return axes_df.iloc[i, i]
        except (IndexError, AttributeError):
            pass
    # Plain 2D array
    try:
        return axes_df[i, i]  # type: ignore[index]
    except (IndexError, TypeError):
        return None


def _plot_corner_matplotlib(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
    show_rejected_inchain: bool = True,
    show_rejected_prior: bool = False,
) -> None:
    """Basic corner plot using matplotlib (fallback)."""
    import matplotlib.pyplot as plt
    import numpy as np

    n = result.n_params
    if n == 0:
        return

    fig, axes = plt.subplots(n, n, figsize=(3 * n, 3 * n))
    if n == 1:
        axes = np.array([[axes]])

    # Use finite samples only
    mask = np.isfinite(result.log_posterior)
    samples = result.samples[mask]
    weights = result.weights[mask] if result.weights is not None else None

    rejected = _rejected_samples_array(result) if show_rejected_inchain else None
    rejected_prior = (
        _load_rejected_prior_overlay(result) if show_rejected_prior else None
    )

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                continue
            if i == j:
                ax.hist(
                    samples[:, i],
                    bins=30,
                    weights=weights,
                    density=True,
                    alpha=0.7,
                )
                ax.set_xlabel(result.param_names[i])
            else:
                # Rejected first (zorder=0) then accepted on top.
                if rejected_prior is not None:
                    ax.scatter(
                        rejected_prior[:, j],
                        rejected_prior[:, i],
                        s=2,
                        alpha=0.2,
                        color="C3",
                        marker=".",
                    )
                if rejected is not None:
                    ax.scatter(
                        rejected[:, j],
                        rejected[:, i],
                        s=4,
                        alpha=0.4,
                        color="C3",
                        label="tachyonic" if (i == 1 and j == 0) else None,
                    )
                ax.scatter(
                    samples[:, j],
                    samples[:, i],
                    c=result.log_likelihood[mask],
                    s=2,
                    alpha=0.5,
                    cmap="viridis",
                )
                if j == 0:
                    ax.set_ylabel(result.param_names[i])
                if i == n - 1:
                    ax.set_xlabel(result.param_names[j])

    fig.suptitle(
        f"Posterior ({result.method}, n={result.n_samples}, "
        f"ESS={result.effective_sample_size():.0f})",
    )
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
