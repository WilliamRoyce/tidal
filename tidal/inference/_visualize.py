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
    show_rejected: bool = True,
) -> None:
    """Generate a corner plot (2D marginals + 1D histograms).

    Uses anesthetic if available, otherwise falls back to matplotlib.

    Parameters
    ----------
    result : InferenceResult
        Inference results with samples and weights.
    output_path : Path | None
        If provided, save figure to this path.
    show : bool
        If True, display the figure interactively.
    show_rejected : bool
        If True (default), overlay samples rejected by the pre-flight
        stability guard (``run_status='tachyonic'``) as a red scatter on
        each 2D panel.  Convey real physics — the unstable region of
        parameter space is part of the result.  Following anesthetic
        conventions (Handley 2019 JOSS) for showing inaccessible regions.
    """
    try:
        _plot_corner_anesthetic(
            result, output_path, show=show, show_rejected=show_rejected
        )
    except ImportError:
        _plot_corner_matplotlib(
            result, output_path, show=show, show_rejected=show_rejected
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


def _plot_corner_anesthetic(
    result: InferenceResult,
    output_path: Path | None = None,
    *,
    show: bool = False,
    show_rejected: bool = True,
) -> None:
    """Corner plot using anesthetic (Handley 2019).

    If ``result.metadata["priors"]`` is present (see issue #308 and
    ``tidal/cli/_sample.py``), this overlay-renders the prior density
    on each 1D marginal axis:

    - ``uniform`` prior → flat red line at 1/(hi-lo)
    - ``log_uniform`` prior → flat red line at 1/(log hi - log lo),
      with the x-axis switched to log scale so a flat posterior =
      flat prior visually. Resolves #309.
    """
    import matplotlib.pyplot as plt

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

    # anesthetic >= 2.0: plot_2d returns an AxesDataFrame whose cells
    # expose .get_figure().  Older versions returned (fig, axes).  Handle
    # both without version-gating.
    ret = samples.plot_2d(result.param_names)
    if isinstance(ret, tuple) and len(ret) == 2:
        fig = ret[0]
        axes_df = ret[1]
    else:
        # AxesDataFrame: row/col indexed by parameter name; diagonal
        # cells are the 1D marginals.
        axes_df = ret
        cell = ret.iloc[0, 0] if hasattr(ret, "iloc") else next(iter(ret))
        fig = cell.get_figure()

    _overlay_priors(axes_df, result)
    if show_rejected:
        _overlay_rejected_anesthetic(axes_df, result)

    fig.suptitle(
        f"Posterior ({result.method}, n={result.n_samples}, "
        f"ESS={result.effective_sample_size():.0f})",
        y=1.02,
    )

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def _overlay_rejected_anesthetic(
    axes_df: object,
    result: InferenceResult,
) -> None:
    """Overlay tachyonic-rejected samples on the corner plot.

    Adds a red scatter on lower-triangle panels of an anesthetic
    AxesDataFrame.  Uses raw scatter rather than ``Samples.plot_2d(
    kind='scatter_2d')`` so we don't depend on anesthetic version-specific
    overload semantics; the visual convention (red, lower triangle, low
    alpha) matches Handley-group corner plots showing inaccessible
    regions.  Side files written by Phase 4 (``_rejected_prior.csv``) are
    also picked up here when available — see
    ``_load_rejected_prior_overlay``.
    """
    rejected = _rejected_samples_array(result)
    rejected_prior = _load_rejected_prior_overlay(result)
    if rejected is None and rejected_prior is None:
        return

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
    show_rejected: bool = True,
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

    rejected = _rejected_samples_array(result) if show_rejected else None
    rejected_prior = _load_rejected_prior_overlay(result) if show_rejected else None

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
