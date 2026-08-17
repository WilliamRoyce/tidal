"""Render functions for sweep-specific plot types in ``tidal plot``.

Provides:

- ``render_sweep_1d``   — metric vs single swept parameter (line plot)
- ``render_sweep_1d_grouped`` — 2-param sweep, one line per group value
- ``render_sweep_2d``   — metric vs two swept parameters (heatmap)
- ``render_sweep_2d_with_overlay`` — 3-panel comparison: TIDAL | analytical | |error|
- ``render_sweep_compare`` — overlay timeseries from multiple runs
- ``render_convergence`` — log-log error vs resolution
- ``render_replicate_convergence`` — SEM vs replicate count diagnostic

References
----------
Smith, R.C. (2013) *Uncertainty Quantification: Theory, Implementation,
and Applications*, SIAM. Ch. 3 (sample statistics, SEM, CI).
"""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.colors import Normalize
    from matplotlib.figure import Figure
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData
    from tidal.measurement._sweep_results import SweepResults


# ------------------------------------------------------------------
# Greek parameter label mapping
# ------------------------------------------------------------------

_GREEK_LABEL_MAP: dict[str, str] = {
    "alpha": r"$\alpha$",
    "alpha1": r"$\alpha_1$",
    "alpha2": r"$\alpha_2$",
    "alpha3": r"$\alpha_3$",
    "beta1": r"$\beta_1$",
    "beta2": r"$\beta_2$",
    "beta3": r"$\beta_3$",
    "gamma": r"$\gamma$",
    "delta1": r"$\delta_1$",
    "delta2": r"$\delta_2$",
    "deltam": r"$\delta_m$",
    "xi": r"$\xi$",
    "chi": r"$\chi$",
    "zeta1": r"$\zeta_1$",
    "zeta2": r"$\zeta_2$",
    "zeta3": r"$\zeta_3$",
    "kappa": r"$\kappa$",
    "lambda": r"$\lambda$",
    "mu": r"$\mu$",
    "sigma": r"$\sigma$",
    "omega": r"$\omega$",
    "phi": r"$\phi$",
    "psi": r"$\psi$",
    "epsilon": r"$\epsilon$",
    "eta": r"$\eta$",
    "theta": r"$\theta$",
    "B0": r"$B_0$",
    "mPhi2": r"$m_\phi^2$",
}


def _greek_label(name: str) -> str:
    """Return a LaTeX label for a parameter name, or the name unchanged."""
    return _GREEK_LABEL_MAP.get(name, name)


# ------------------------------------------------------------------
# Arctan axis transform helpers
# ------------------------------------------------------------------


def _apply_arctan_axis(ax: Axes, axis: str = "x") -> None:
    """Apply arctan transform to an axis for unbounded parameter display.

    Compresses (-∞, +∞) into a bounded range while showing actual
    parameter values as tick labels at human-readable positions.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes to transform.
    axis : str
        ``"x"`` or ``"y"`` — which axis to transform.
    """
    from matplotlib.ticker import FuncFormatter

    # Tick positions in real-parameter space
    tick_values = [-100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100]

    def _fmt(val: float, _pos: object) -> str:
        """Format arctan-space tick as the original parameter value."""
        # Inverse arctan: val in arctan space → original value
        original = np.tan(val)
        if abs(original) < 0.01:  # noqa: PLR2004
            return "0"
        if abs(original) >= 100:  # noqa: PLR2004
            return f"{original:.0f}"
        if abs(original) >= 1:
            return f"{original:.1f}"
        return f"{original:.2f}"

    formatter = FuncFormatter(_fmt)
    arctan_ticks = [float(np.arctan(v)) for v in tick_values]

    if axis == "x":
        ax.set_xticks(arctan_ticks)
        ax.xaxis.set_major_formatter(formatter)
    else:
        ax.set_yticks(arctan_ticks)
        ax.yaxis.set_major_formatter(formatter)


def _arctan_transform_data(
    values: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Transform parameter values to arctan space."""
    return np.arctan(values)


# ------------------------------------------------------------------
# Analytical overlay helpers
# ------------------------------------------------------------------


def _build_overlay_scalars(results: SweepResults) -> dict[str, float]:
    """Build scalar namespace from fixed params and numeric sim settings.

    Returns a dict suitable for use in overlay formula evaluation,
    containing all fixed sweep parameters (e.g. ``kappa``) and any
    numeric simulation settings (e.g. ``t_end``, ``grid_shape``).
    """
    ns: dict[str, float] = dict(results.fixed_params)
    for k, v in results.sim_settings.items():
        with contextlib.suppress(TypeError, ValueError):
            ns[k] = float(v)
    return ns


def _evaluate_sweep_overlay(
    formula: str,
    param_arrays: dict[str, NDArray[np.float64]],
    scalar_params: dict[str, float],
) -> NDArray[np.float64]:
    """Evaluate an overlay formula over swept parameter arrays.

    Parameters
    ----------
    formula : str
        Python expression using swept parameter names, fixed parameter
        names, sim settings (``t_end``, ``grid_shape``), and standard
        math functions (``sin``, ``cos``, ``sqrt``, ``pi``, ``exp``, ...).
        Example: ``'sin(kappa * B0 * t_end / 2)**2'``.
    param_arrays : dict
        Swept parameters as numpy arrays (broadcast-ready).  For 1D
        sweeps this is a 1-D array; for 2D sweeps use meshgrid outputs.
    scalar_params : dict
        Fixed parameters and numeric sim settings as Python floats.

    Returns
    -------
    NDArray
        Evaluated formula values with the same shape as the input arrays.

    Raises
    ------
    ValueError
        If the formula contains a syntax/name error or evaluation fails.
    """
    from tidal.cli._simulate import (
        FORMULA_NAMESPACE,  # pyright: ignore[reportPrivateUsage]
    )

    ns: dict[str, object] = {**FORMULA_NAMESPACE, **scalar_params, **param_arrays}
    try:
        from tidal.cli._simulate import safe_formula_eval

        result = safe_formula_eval(formula, ns)
    except Exception as exc:
        msg = f"Error evaluating overlay formula {formula!r}: {exc}"
        raise ValueError(msg) from exc
    out = np.asarray(result, dtype=np.float64)
    # Broadcast scalar result (e.g. formula is constant, doesn't reference swept param)
    # to match the shape of the first swept-parameter array.
    if out.ndim == 0 and param_arrays:
        ref = next(iter(param_arrays.values()))
        out = np.broadcast_to(out, np.asarray(ref).shape).copy()
    return out


# ------------------------------------------------------------------
# Colour norm and quality overlay helpers
# ------------------------------------------------------------------


def _build_norm(
    values: NDArray[np.float64],
    *,
    log_scale: bool = False,
    divergent_center: float | None = None,
) -> Normalize:
    """Build a matplotlib Normalize, LogNorm, or TwoSlopeNorm.

    Returns a ``Normalize``, ``LogNorm``, or ``TwoSlopeNorm`` instance.
    Returns a default linear ``Normalize(0, 1)`` if *values* is empty or
    all-zero with ``log_scale=True``.

    Parameters
    ----------
    divergent_center : float or None
        If provided, use ``TwoSlopeNorm`` centered at this value.
        Useful for colormaps where a physical baseline (e.g. P_EM)
        should map to the colormap midpoint.
    """
    from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm

    if len(values) == 0:
        return Normalize(0, 1)
    vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
    if log_scale:
        # Clamp vmin to smallest positive value
        pos = values[values > 0]
        if len(pos) == 0:
            return Normalize(0, 1)
        vmin = float(pos.min())
        vmax = max(vmax, vmin * 1.01)  # avoid vmin == vmax
        return LogNorm(vmin=vmin, vmax=vmax)
    if divergent_center is not None and vmin < divergent_center < vmax:
        return TwoSlopeNorm(vcenter=divergent_center, vmin=vmin, vmax=vmax)
    return Normalize(vmin=vmin, vmax=vmax)


def _overlay_quality_hatching(
    ax: Axes,
    p1_vals: NDArray[np.float64],
    p2_vals: NDArray[np.float64],
    quality_grid: NDArray[np.object_],
) -> None:
    """Overlay hatching on 2D grid cells with non-'good' crossing quality.

    Hatching patterns:
    - ``"coarse"`` — light diagonal hatching (``"/"``)
    - ``"edge"``   — cross hatching (``"x"``)
    - ``"none"``   — dense diagonal (``"///"``)
    """
    import matplotlib.patches as mpatches

    n2, n1 = quality_grid.shape
    if n1 < 2 or n2 < 2:  # noqa: PLR2004
        return

    # Cell widths (use spacing between sorted unique values)
    dp1 = np.diff(p1_vals)
    dp2 = np.diff(p2_vals)

    hatch_map = {"coarse": "/", "edge": "x", "none": "///"}

    for i2 in range(n2):
        for i1 in range(n1):
            q = str(quality_grid[i2, i1])
            if q not in hatch_map:
                continue
            # Cell bounds: nearest shading means centred on (p1[i1], p2[i2])
            w = float(dp1[min(i1, len(dp1) - 1)])
            h = float(dp2[min(i2, len(dp2) - 1)])
            x0 = float(p1_vals[i1]) - w / 2
            y0 = float(p2_vals[i2]) - h / 2
            rect = mpatches.Rectangle(
                (x0, y0),
                w,
                h,
                hatch=hatch_map[q],
                fill=False,
                edgecolor="gray",
                linewidth=0.3,
                zorder=2,
            )
            ax.add_patch(rect)


def render_sweep_1d(  # noqa: PLR0913, PLR0912, PLR0915, C901
    ax: Axes,
    results: SweepResults,
    metric: str,
    *,
    overlay: str | None = None,
    log_y: bool = False,
    thresholds: list[str] | None = None,
    scatter: bool = False,
    annotate_extremes: bool = False,
    arctan_axes: bool = False,
) -> None:
    """Plot a scalar metric vs a single swept parameter.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes to draw on.
    results : SweepResults
        Loaded sweep data.
    metric : str
        Column name in results.rows to plot on the y-axis.
    overlay : str or None
        Optional analytical formula to overlay as a dashed reference
        curve.  The formula may reference the swept parameter by name,
        fixed parameters (e.g. ``kappa``), and simulation settings
        (e.g. ``t_end``).  Standard math functions (``sin``, ``cos``,
        ``sqrt``, ``pi``, ``exp``) are available.
        Example: ``'sin(kappa * B0 * t_end / 2)**2'``.
    log_y : bool
        Use logarithmic y-axis scale.
    thresholds : list[str] or None
        Horizontal threshold lines in ``"VALUE[:LABEL[:COLOR]]"`` format.
        Example: ``["1.0:P=1:red", "0.1::orange"]``.
    scatter : bool
        Use scatter plot instead of connected line plot.
    annotate_extremes : bool
        Annotate max/min values with arrows and labels.

    Raises
    ------
    ValueError
        If not exactly 1 swept parameter, or if the overlay formula fails.
    """
    param_names = list(results.swept_params.keys())
    if len(param_names) != 1:
        msg = (
            f"render_sweep_1d expects exactly 1 swept parameter, "
            f"got {len(param_names)}: {param_names}"
        )
        raise ValueError(msg)

    param_name = param_names[0]

    if results.has_replicates:
        agg = results.aggregate()
        x = agg.column(param_name)
        if arctan_axes:
            x = _arctan_transform_data(np.asarray(x, dtype=np.float64))
        y_mean = agg.column(f"{metric}_mean")
        y_std = agg.column(f"{metric}_std")
        num_label = "TIDAL" if overlay else None
        if scatter:
            ax.scatter(
                x,
                y_mean,
                c="tab:blue",
                s=20,
                alpha=0.8,
                zorder=3,
                edgecolors="none",
                label=num_label,
            )
        else:
            ax.plot(
                x,
                y_mean,
                "o-",
                color="tab:blue",
                linewidth=1.5,
                markersize=5,
                label=num_label,
            )
        ax.fill_between(
            x,
            y_mean - y_std,
            y_mean + y_std,
            alpha=0.2,
            color="tab:blue",
            label=r"$\pm 1\sigma$",
        )
        y_for_annotate = y_mean
    else:
        x = results.column(param_name)
        if arctan_axes:
            x = _arctan_transform_data(np.asarray(x, dtype=np.float64))
        y = results.column(metric)
        num_label = "TIDAL" if overlay else None
        if scatter:
            ax.scatter(
                x,
                y,
                c="tab:blue",
                s=20,
                alpha=0.8,
                zorder=3,
                edgecolors="none",
                label=num_label,
            )
        else:
            ax.plot(
                x,
                y,
                "o-",
                color="tab:blue",
                linewidth=1.5,
                markersize=5,
                label=num_label,
            )
        y_for_annotate = y

    # Annotate max/min values
    if annotate_extremes and len(y_for_annotate) > 0:
        y_arr = np.asarray(y_for_annotate, dtype=np.float64)
        x_arr = np.asarray(x, dtype=np.float64)
        valid = np.isfinite(y_arr)
        if np.any(valid):
            i_max = int(np.nanargmax(y_arr))
            i_min = int(np.nanargmin(y_arr))
            ax.annotate(
                f"max={y_arr[i_max]:.3g}",
                (x_arr[i_max], y_arr[i_max]),
                textcoords="offset points",
                xytext=(8, -12),
                fontsize=8,
                color="#B71C1C",
                fontweight="bold",
                arrowprops={"arrowstyle": "->", "color": "#B71C1C", "lw": 0.8},
            )
            if i_min != i_max:
                ax.annotate(
                    f"min={y_arr[i_min]:.3g}",
                    (x_arr[i_min], y_arr[i_min]),
                    textcoords="offset points",
                    xytext=(8, 8),
                    fontsize=8,
                    color="#1565C0",
                    arrowprops={"arrowstyle": "->", "color": "#1565C0", "lw": 0.8},
                )

    # Analytical overlay curve
    if overlay is not None:
        x_arr = np.asarray(x, dtype=np.float64)
        scalar_ns = _build_overlay_scalars(results)
        y_ref = _evaluate_sweep_overlay(overlay, {param_name: x_arr}, scalar_ns)
        ax.plot(
            x_arr,
            y_ref,
            "--",
            color="tab:orange",
            linewidth=1.5,
            label="analytical",
            zorder=3,
        )
        ax.legend(fontsize="small")

    # Threshold lines
    if thresholds:
        for spec in thresholds:
            parts = spec.split(":")
            val = float(parts[0])
            label = parts[1] if len(parts) > 1 and parts[1] else None
            color = parts[2] if len(parts) > 2 and parts[2] else "red"  # noqa: PLR2004
            ax.axhline(val, color=color, linestyle=":", alpha=0.7, label=label)

    if overlay or thresholds:
        ax.legend(fontsize="small")

    ax.set_xlabel(_greek_label(param_name))
    ax.set_ylabel(metric)
    ax.set_title(f"{metric} vs {_greek_label(param_name)}")
    ax.grid(visible=True, alpha=0.3)

    if log_y:
        ax.set_yscale("log")

    if arctan_axes:
        _apply_arctan_axis(ax, "x")


def _plot_1d_metric(  # noqa: PLR0913
    ax: Axes,
    results: SweepResults,
    param_name: str,
    metric: str,
    color: str,
    *,
    log_y: bool = False,
) -> None:
    """Plot a single metric on an axis, with error bands if replicates exist."""
    if results.has_replicates:
        agg = results.aggregate()
        x = agg.column(param_name)
        y_mean = agg.column(f"{metric}_mean")
        y_std = agg.column(f"{metric}_std")
        ax.plot(x, y_mean, "o-", color=color, linewidth=1.5, markersize=5)
        ax.fill_between(x, y_mean - y_std, y_mean + y_std, alpha=0.2, color=color)
    else:
        x = results.column(param_name)
        y = results.column(metric)
        ax.plot(x, y, "o-", color=color, linewidth=1.5, markersize=5)
    if log_y:
        ax.set_yscale("log")


def render_sweep_1d_multi(
    fig: Figure,
    results: SweepResults,
    metrics: list[str],
    *,
    log_y: bool = False,
) -> None:
    """Plot multiple metrics vs a single swept parameter.

    Uses subplots stacked vertically, sharing the x-axis.
    When replicates exist, shows mean +/- std error bands.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure.
    results : SweepResults
        Loaded sweep data.
    metrics : list[str]
        Column names to plot.
    log_y : bool
        Use logarithmic y-axis scale on all subplots.

    Raises
    ------
    ValueError
        If not exactly 1 swept parameter.
    """
    param_names = list(results.swept_params.keys())
    if len(param_names) != 1:
        msg = f"Expected 1 swept parameter, got {len(param_names)}"
        raise ValueError(msg)

    param_name = param_names[0]
    n = len(metrics)

    axes = fig.subplots(n, 1, sharex=True, squeeze=False)
    colors = ("tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple")

    for i, metric in enumerate(metrics):
        ax = axes[i, 0]
        _plot_1d_metric(
            ax,
            results,
            param_name,
            metric,
            colors[i % len(colors)],
            log_y=log_y,
        )
        ax.set_ylabel(metric)
        ax.grid(visible=True, alpha=0.3)

    axes[-1, 0].set_xlabel(param_name)
    fig.suptitle(f"Sweep: {', '.join(metrics)} vs {param_name}")


def render_sweep_1d_grouped(  # noqa: PLR0913
    ax: Axes,
    results: SweepResults,
    metric: str,
    *,
    x_param: str,
    group_param: str,
    overlay: str | None = None,
    log_y: bool = False,
    log_x: bool = False,
) -> None:
    """Plot a scalar metric vs one swept parameter, one line per value of another.

    For a 2-parameter sweep, group rows by ``group_param`` and draw one
    connected line (markers + line) per unique group value, sharing a
    common x-axis ``x_param``.  This is the canonical figure for
    "collapse across B₀" tests of linearized-regime validity.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes to draw on.
    results : SweepResults
        Loaded sweep data with **exactly 2** swept parameters.
    metric : str
        Column name in ``results.rows`` for the y-values.
    x_param : str
        Swept parameter name for the x-axis.
    group_param : str
        Swept parameter name to group lines by.
    overlay : str or None
        Optional analytical formula in ``x_param`` (and fixed params /
        sim settings).  Drawn as a single dashed black reference curve
        — assumed to be independent of ``group_param``.
    log_y : bool
        Use logarithmic y-axis scale.
    log_x : bool
        Use logarithmic x-axis scale.

    Raises
    ------
    ValueError
        If the sweep does not have exactly 2 swept parameters, or if
        ``{x_param, group_param}`` is not the set of swept parameters.
    """
    import matplotlib.pyplot as plt

    param_names = list(results.swept_params.keys())
    if len(param_names) != 2:  # noqa: PLR2004
        msg = (
            f"render_sweep_1d_grouped expects exactly 2 swept parameters, "
            f"got {len(param_names)}: {param_names}"
        )
        raise ValueError(msg)
    if {x_param, group_param} != set(param_names):
        msg = (
            f"x_param={x_param!r} and group_param={group_param!r} must match "
            f"the swept parameters {param_names}"
        )
        raise ValueError(msg)

    # Extract all rows as parallel arrays
    x_all = np.asarray(results.column(x_param), dtype=np.float64)
    g_all = np.asarray(results.column(group_param), dtype=np.float64)
    y_all = np.asarray(results.column(metric), dtype=np.float64)

    # Unique group values, sorted numerically
    group_vals = np.unique(g_all[np.isfinite(g_all)])
    n_groups = len(group_vals)
    if n_groups == 0:
        msg = f"no finite values found for group_param={group_param!r}"
        raise ValueError(msg)

    cmap = plt.colormaps["viridis"]
    # Span [0.1, 0.9] to avoid the lightest/darkest extremes
    color_positions = np.linspace(0.1, 0.9, max(n_groups, 1))

    for idx, gv in enumerate(group_vals):
        mask = np.isclose(g_all, gv) & np.isfinite(y_all)
        if not np.any(mask):
            continue
        xs = x_all[mask]
        ys = y_all[mask]
        # Sort by x for a clean connected line
        order = np.argsort(xs)
        xs, ys = xs[order], ys[order]
        label = f"{_greek_label(group_param)} = {gv:.3g}"
        ax.plot(
            xs,
            ys,
            "o-",
            color=cmap(float(color_positions[idx])),
            linewidth=1.5,
            markersize=5,
            label=label,
            zorder=3,
        )

    # Analytical overlay (single curve, group-independent)
    if overlay is not None:
        x_sorted = np.sort(np.unique(x_all[np.isfinite(x_all)]))
        scalar_ns = _build_overlay_scalars(results)
        y_ref = _evaluate_sweep_overlay(overlay, {x_param: x_sorted}, scalar_ns)
        ax.plot(
            x_sorted,
            y_ref,
            "--",
            color="black",
            linewidth=1.5,
            label="analytical",
            zorder=4,
        )

    ax.set_xlabel(_greek_label(x_param))
    ax.set_ylabel(metric)
    ax.set_title(
        f"{metric} vs {_greek_label(x_param)} (grouped by {_greek_label(group_param)})",
    )
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize="small", loc="best")

    if log_y:
        ax.set_yscale("log")
    if log_x:
        ax.set_xscale("log")


def render_sweep_2d(  # noqa: PLR0913
    ax: Axes,
    results: SweepResults,
    metric: str,
    *,
    log_scale: bool = False,
    divergent_center: float | None = None,
    cmap_name: str = "viridis",
    clamp_color: tuple[float, float] | None = None,
) -> None:
    """Plot a scalar metric over two swept parameters.

    For grid-aligned data (Cartesian sweeps), renders a pcolormesh heatmap.
    For scattered data (LHS/Sobol), renders a colored scatter plot with
    optional interpolation background.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes to draw on.
    results : SweepResults
        Loaded sweep data.
    metric : str
        Column name for the color values.
    log_scale : bool
        Use logarithmic colorbar (useful for ``inv_B_min``).
    divergent_center : float or None
        Center colormap at this value using ``TwoSlopeNorm``.
    clamp_color : tuple[float, float] or None
        Clamp colorbar range to ``(vmin, vmax)``.

    Raises
    ------
    ValueError
        If not exactly 2 swept parameters.
    """
    param_names = list(results.swept_params.keys())
    if len(param_names) != 2:  # noqa: PLR2004
        msg = (
            f"render_sweep_2d expects exactly 2 swept parameters, "
            f"got {len(param_names)}: {param_names}"
        )
        raise ValueError(msg)

    p1_name, p2_name = param_names

    if _is_scattered_data(results):
        _render_2d_scattered(
            ax,
            results,
            p1_name,
            p2_name,
            metric,
            log_scale=log_scale,
            divergent_center=divergent_center,
            cmap_name=cmap_name,
        )
    else:
        _render_2d_grid(
            ax,
            results,
            p1_name,
            p2_name,
            metric,
            log_scale=log_scale,
            divergent_center=divergent_center,
            cmap_name=cmap_name,
            clamp_color=clamp_color,
        )


def _render_2d_scattered(  # noqa: PLR0913
    ax: Axes,
    results: SweepResults,
    p1_name: str,
    p2_name: str,
    metric: str,
    *,
    log_scale: bool = False,
    divergent_center: float | None = None,
    cmap_name: str = "viridis",
) -> None:
    """Render 2D sweep as scatter + interpolation background."""
    p1 = np.array(results.column(p1_name), dtype=np.float64)
    p2 = np.array(results.column(p2_name), dtype=np.float64)
    metric_vals = np.array(results.column(metric), dtype=np.float64)

    valid = np.isfinite(metric_vals)
    if log_scale:
        valid &= metric_vals > 0
    p1, p2, metric_vals = p1[valid], p2[valid], metric_vals[valid]

    if len(metric_vals) == 0:
        ax.text(
            0.5,
            0.5,
            f"No valid data for {metric}",
            transform=ax.transAxes,
            ha="center",
        )
        return

    norm = _build_norm(
        metric_vals,
        log_scale=log_scale,
        divergent_center=divergent_center,
    )

    # Interpolation background (if enough points)
    if len(metric_vals) >= 10:  # noqa: PLR2004
        try:
            from scipy.interpolate import griddata  # type: ignore[import-untyped]

            n_grid = 100
            p1g = np.linspace(float(p1.min()), float(p1.max()), n_grid)
            p2g = np.linspace(float(p2.min()), float(p2.max()), n_grid)
            p1m, p2m = np.meshgrid(p1g, p2g)
            z = np.asarray(
                griddata(
                    np.column_stack([p1, p2]),
                    metric_vals,
                    (p1m, p2m),
                    method="cubic",
                ),
                dtype=np.float64,
            )
            ax.pcolormesh(
                p1g,
                p2g,
                z,
                shading="auto",
                cmap=cmap_name,
                alpha=0.3,
                norm=norm,
            )
        except (ValueError, ImportError):
            pass  # Fall back to scatter-only

    sc = ax.scatter(
        p1,
        p2,
        c=metric_vals,
        cmap=cmap_name,
        s=60,
        edgecolors="k",
        linewidths=0.5,
        norm=norm,
        zorder=5,
    )
    ax.set_xlabel(p1_name)
    ax.set_ylabel(p2_name)
    ax.set_title(metric)
    ax.figure.colorbar(sc, ax=ax, label=metric)  # type: ignore[union-attr]


def _render_2d_grid(  # noqa: PLR0913, C901
    ax: Axes,
    results: SweepResults,
    p1_name: str,
    p2_name: str,
    metric: str,
    *,
    log_scale: bool = False,
    divergent_center: float | None = None,
    cmap_name: str = "viridis",
    clamp_color: tuple[float, float] | None = None,
) -> None:
    """Render 2D sweep as pcolormesh heatmap (for grid-aligned data).

    Diverged or invalid runs are overlaid with cross-hatching.
    """
    p1_vals = np.sort(results.swept_params[p1_name])
    p2_vals = np.sort(results.swept_params[p2_name])

    n1, n2 = len(p1_vals), len(p2_vals)
    grid = np.full((n2, n1), np.nan)
    diverged_grid = np.zeros((n2, n1), dtype=bool)
    quality_grid: NDArray[np.object_] | None = None

    # Check for crossing_quality column (critical field results)
    has_quality = results.rows and "crossing_quality" in results.rows[0]
    if has_quality:
        quality_grid = np.full((n2, n1), "", dtype=object)

    for row in results.rows:
        v1 = row.get(p1_name)
        v2 = row.get(p2_name)
        val = row.get(metric)
        if v1 is None or v2 is None:
            continue
        i1 = min(int(np.searchsorted(p1_vals, float(v1))), n1 - 1)
        i2 = min(int(np.searchsorted(p2_vals, float(v2))), n2 - 1)

        status = row.get("run_status", "success")
        if status != "success" or val is None:
            diverged_grid[i2, i1] = True
            continue

        try:
            grid[i2, i1] = float(val)
        except (TypeError, ValueError):
            diverged_grid[i2, i1] = True

        if quality_grid is not None:
            quality_grid[i2, i1] = row.get("crossing_quality", "")

    # Build norm (linear or log), with optional clamping
    valid_vals = grid[np.isfinite(grid)]
    if log_scale:
        valid_vals = valid_vals[valid_vals > 0]

    if clamp_color is not None and len(valid_vals) > 0:
        vmin_c, vmax_c = clamp_color
        valid_vals = np.clip(valid_vals, vmin_c, vmax_c)
        grid = np.where(np.isfinite(grid), np.clip(grid, vmin_c, vmax_c), grid)

    norm = (
        _build_norm(valid_vals, log_scale=log_scale, divergent_center=divergent_center)
        if len(valid_vals) > 0
        else None
    )

    im = ax.pcolormesh(
        p1_vals,
        p2_vals,
        grid,
        shading="nearest",
        cmap=cmap_name,
        norm=norm,
    )
    ax.set_xlabel(_greek_label(p1_name))
    ax.set_ylabel(_greek_label(p2_name))
    ax.set_title(metric)
    cb = ax.figure.colorbar(im, ax=ax, label=metric)  # type: ignore[union-attr]

    # Reference line on colorbar at divergent center
    if divergent_center is not None:
        cb.ax.axhline(divergent_center, color="k", linewidth=0.8)

    # Diverged-run hatching overlay
    if np.any(diverged_grid) and n1 >= 2 and n2 >= 2:  # noqa: PLR2004
        from matplotlib.patches import Patch

        p1m, p2m = np.meshgrid(p1_vals, p2_vals)
        ax.contourf(
            p1m,
            p2m,
            diverged_grid.astype(float),
            levels=[0.5, 1.5],
            colors="none",
            hatches=["///"],
        )
        ax.legend(
            handles=[
                Patch(
                    facecolor="white",
                    edgecolor="grey",
                    hatch="///",
                    label="Diverged / invalid",
                ),
            ],
            loc="lower left",
            fontsize=8,
            framealpha=0.9,
        )

    # Quality-flag hatching overlay for critical field results
    if quality_grid is not None:
        _overlay_quality_hatching(ax, p1_vals, p2_vals, quality_grid)


def render_sweep_2d_with_overlay(
    fig: Figure,
    results: SweepResults,
    metric: str,
    overlay: str,
    *,
    log_scale: bool = False,
) -> None:
    """3-panel comparison: TIDAL numerical | analytical formula | absolute error.

    Creates three side-by-side pcolormesh panels so that the TIDAL
    simulation result can be directly compared against an analytical
    reference formula.  The TIDAL and analytical panels share the same
    colour scale; the error panel uses a separate sequential colourmap.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure (caller should size it to ~(15, 5)).
    results : SweepResults
        Loaded 2D sweep data (exactly 2 swept parameters required).
    metric : str
        Column name in results.rows for the y-values (e.g. ``"P_final"``).
    overlay : str
        Analytical formula string evaluated on a meshgrid of the two swept
        parameters.  Fixed parameters and numeric sim settings are available
        as scalar variables.  Example:
        ``'sin(kappa * Bpeak * R * sqrt(pi/2))**2'``.
    log_scale : bool
        Use logarithmic colour scale for numerical and analytical panels.

    Raises
    ------
    ValueError
        If not exactly 2 swept parameters or overlay evaluation fails.
    """
    param_names = list(results.swept_params.keys())
    if len(param_names) != 2:  # noqa: PLR2004
        msg = (
            f"render_sweep_2d_with_overlay expects exactly 2 swept parameters, "
            f"got {len(param_names)}: {param_names}"
        )
        raise ValueError(msg)

    p1_name, p2_name = param_names
    p1_vals = np.sort(np.asarray(results.swept_params[p1_name], dtype=np.float64))
    p2_vals = np.sort(np.asarray(results.swept_params[p2_name], dtype=np.float64))
    n1, n2 = len(p1_vals), len(p2_vals)

    # Build numerical grid
    z_num = np.full((n2, n1), np.nan)
    for row in results.rows:
        v1 = row.get(p1_name)
        v2 = row.get(p2_name)
        val = row.get(metric)
        if v1 is None or v2 is None or val is None:
            continue
        i1 = min(int(np.searchsorted(p1_vals, float(v1))), n1 - 1)
        i2 = min(int(np.searchsorted(p2_vals, float(v2))), n2 - 1)
        z_num[i2, i1] = float(val)

    # Evaluate analytical formula on meshgrid
    p1m, p2m = np.meshgrid(p1_vals, p2_vals)
    scalar_ns = _build_overlay_scalars(results)
    z_anal = _evaluate_sweep_overlay(overlay, {p1_name: p1m, p2_name: p2m}, scalar_ns)

    # Absolute error
    z_err = np.abs(z_num - z_anal)

    # Shared colour scale for numerical + analytical panels
    valid = z_num[np.isfinite(z_num)]
    if log_scale:
        valid = valid[valid > 0]
    norm = _build_norm(valid, log_scale=log_scale) if len(valid) > 0 else None

    axes = fig.subplots(1, 3)

    # Panel 0: TIDAL numerical
    ax0 = axes[0]
    im0 = ax0.pcolormesh(
        p1_vals,
        p2_vals,
        z_num,
        shading="nearest",
        cmap="viridis",
        norm=norm,
    )
    ax0.set_xlabel(p1_name)
    ax0.set_ylabel(p2_name)
    ax0.set_title(f"{metric} (TIDAL)")
    fig.colorbar(im0, ax=ax0, label=metric)

    # Panel 1: analytical formula
    ax1 = axes[1]
    im1 = ax1.pcolormesh(
        p1_vals,
        p2_vals,
        z_anal,
        shading="nearest",
        cmap="viridis",
        norm=norm,
    )
    ax1.set_xlabel(p1_name)
    ax1.set_ylabel(p2_name)
    ax1.set_title(f"{metric} (analytical)")
    fig.colorbar(im1, ax=ax1, label=metric)

    # Panel 2: absolute error
    ax2 = axes[2]
    im2 = ax2.pcolormesh(p1_vals, p2_vals, z_err, shading="nearest", cmap="YlOrRd")
    ax2.set_xlabel(p1_name)
    ax2.set_ylabel(p2_name)
    ax2.set_title("|error|")
    fig.colorbar(im2, ax=ax2, label="|error|")

    fig.suptitle(f"{metric}: TIDAL vs analytical", y=1.02)
    fig.tight_layout()


def render_sweep_compare(
    ax: Axes,
    results: SweepResults,
    measurement: str,
    spec_path: str | None = None,
) -> None:
    """Overlay timeseries (e.g. P(t)) from each sweep run, color-coded.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Loaded sweep data.
    measurement : str
        Measurement type to overlay (``"conversion"`` or ``"conservation"``).
    spec_path : str or None
        Override spec path (default: from sweep metadata).

    Raises
    ------
    ValueError
        If not exactly 1 swept parameter or unsupported measurement.
    """
    import matplotlib.pyplot as plt

    param_names = list(results.swept_params.keys())
    if len(param_names) != 1:
        msg = (
            f"sweep-compare requires exactly 1 swept parameter, got {len(param_names)}"
        )
        raise ValueError(msg)

    param_name = param_names[0]
    param_vals = results.column(param_name)

    # Color scale
    norm = _safe_normalize(param_vals)
    cmap = plt.colormaps["viridis"]

    resolved_spec = spec_path or results.spec_path
    if not resolved_spec:
        msg = "No spec path available — provide --spec"
        raise ValueError(msg)

    from pathlib import Path

    from tidal.measurement._io import SimulationData
    from tidal.symbolic import load_equation_system

    spec = load_equation_system(Path(resolved_spec))

    for i, run_dir in enumerate(results.run_dirs):
        if not run_dir.is_dir():
            continue
        try:
            data = SimulationData.load(run_dir, spec)
        except (FileNotFoundError, ValueError, OSError):
            continue

        color: tuple[float, ...] = cmap(norm[i])  # type: ignore[assignment]
        label = f"{param_name}={param_vals[i]:.3g}"

        if measurement == "conversion":
            _overlay_conversion(ax, data, results, color, label)  # type: ignore[arg-type]
        elif measurement == "conservation":
            _overlay_conservation(ax, data, color, label)  # type: ignore[arg-type]
        else:
            msg = f"Unsupported measurement for compare: {measurement}"
            raise ValueError(msg)

    ax.legend(fontsize=7, loc="best")
    ax.grid(visible=True, alpha=0.3)


def _overlay_conversion(
    ax: Axes,
    data: SimulationData,
    results: SweepResults,
    color: tuple[float, ...],
    label: str,
) -> None:
    """Overlay conversion P(t) for one run."""
    from tidal.measurement import compute_group_conversion

    source = results.source_fields
    target = results.target_fields
    try:
        conv = compute_group_conversion(data, source or [], target)
        ax.plot(conv.times, conv.probability, color=color, label=label, linewidth=1)
        ax.set_xlabel("t")
        ax.set_ylabel("P(t)")
        ax.set_title("Conversion probability")
    except (ValueError, TypeError, KeyError):
        pass  # Skip runs where conversion cannot be computed


def _overlay_conservation(
    ax: Axes,
    data: SimulationData,
    color: tuple[float, ...],
    label: str,
) -> None:
    """Overlay energy conservation for one run."""
    from tidal.measurement import check_energy_conservation

    try:
        diag = check_energy_conservation(data)
        ax.plot(
            data.times,
            np.abs(diag.relative_error),
            color=color,
            label=label,
            linewidth=1,
        )
        ax.set_yscale("log")
        ax.set_xlabel("t")
        ax.set_ylabel("|dE/E|")
        ax.set_title("Energy conservation")
    except (ValueError, TypeError, KeyError):
        pass  # Skip runs where conservation cannot be computed


def render_convergence(
    ax: Axes,
    results: SweepResults,
    metric: str,
) -> None:
    """Log-log plot of a metric vs grid resolution.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Convergence sweep data.
    metric : str
        Column name for the y-axis (e.g. ``"max_energy_error"``).
    """
    if results.has_replicates:
        agg = results.aggregate()
        sizes = agg.column("grid_shape")
        values = agg.column(f"{metric}_mean")
        y_std = agg.column(f"{metric}_std")
    else:
        sizes = results.column("grid_shape")
        values = results.column(metric)
        y_std = None

    # Filter valid points
    mask = np.isfinite(values) & (values > 0)
    if not np.any(mask):
        ax.text(
            0.5,
            0.5,
            f"No valid data for {metric}",
            transform=ax.transAxes,
            ha="center",
        )
        return

    h = 1.0 / sizes[mask]
    y = values[mask]
    n = sizes[mask]

    if y_std is not None:
        yerr = y_std[mask]
        ax.errorbar(
            n,
            y,
            yerr=yerr,
            fmt="o-",
            color="tab:blue",
            linewidth=1.5,
            markersize=6,
            capsize=3,
            label="measured",
        )
        ax.set_yscale("log")
        ax.set_xscale("log")
    else:
        ax.loglog(
            n,
            y,
            "o-",
            color="tab:blue",
            linewidth=1.5,
            markersize=6,
            label="measured",
        )

    # Fit convergence order
    if len(h) >= 2:  # noqa: PLR2004
        coeffs = np.polyfit(np.log(h), np.log(y), 1)
        order = coeffs[0]
        ax.set_title(f"{metric} — convergence order: {order:.2f}")

        # Reference lines
        h_ref = np.array([h.min(), h.max()])
        n_ref = 1.0 / h_ref
        for p, style, lbl in [(2, "--", "O(h²)"), (4, ":", "O(h⁴)")]:
            y_ref = y[0] * (h_ref / h_ref[0]) ** p
            ax.loglog(n_ref, y_ref, style, color="gray", alpha=0.5, label=lbl)
    else:
        ax.set_title(metric)

    ax.set_xlabel("N (grid points)")
    ax.set_ylabel(metric)
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3, which="both")


# -- Helpers -------------------------------------------------------------------


def _safe_normalize(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Normalize values to [0, 1] range for colormap mapping.

    NaN values are preserved (not included in min/max computation).
    Returns 0.5 for all-NaN or constant arrays.
    """
    if len(values) == 0 or np.all(np.isnan(values)):
        return np.full_like(values, 0.5)
    vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
    if vmax == vmin:
        return np.full_like(values, 0.5)
    return (values - vmin) / (vmax - vmin)


def _is_scattered_data(results: SweepResults) -> bool:
    """Detect whether sweep data is scattered (LHS/Sobol) vs grid-aligned."""
    strategy = results.metadata.get("sampling_strategy", "")
    if strategy in {"latin_hypercube", "sobol"}:
        return True
    if strategy == "grid":
        return False
    # Heuristic fallback for old data without metadata
    n = results.n_runs
    if n < 2:  # noqa: PLR2004
        return False
    for name in results.swept_params:
        vals = np.array(results.column(name), dtype=np.float64)
        n_unique = len(np.unique(vals[np.isfinite(vals)]))
        if n_unique < n * 0.8:
            return False
    return True


# -- Advanced visualization (F8) ---------------------------------------------


def render_sweep_parallel(
    ax: Axes,
    results: SweepResults,
    metric: str,
    *,
    cmap_name: str = "viridis",
    divergent_center: float | None = None,
) -> None:
    """Parallel coordinates: each axis is a parameter + metric, color = metric.

    Lines are sorted by metric value (low-to-high) so high-metric runs
    draw on top.  Alpha and linewidth scale with distance from median to
    emphasize extreme (most interesting) values.

    Diverged runs (NaN metric) are drawn as dashed grey lines so they
    still show which parameter regions are unstable.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Sweep data with 2+ swept parameters.
    metric : str
        Metric for coloring the polylines.
    cmap_name : str
        Matplotlib colormap name (default ``"viridis"``).
    divergent_center : float or None
        Center colormap at this value using ``TwoSlopeNorm``.

    Raises
    ------
    ValueError
        If fewer than 2 swept parameters.
    """
    import matplotlib.pyplot as plt

    param_names = list(results.swept_params.keys())
    if len(param_names) < 2:  # noqa: PLR2004
        msg = f"sweep-parallel requires 2+ swept parameters, got {len(param_names)}"
        raise ValueError(msg)

    metric_vals = np.array(results.column(metric), dtype=np.float64)
    ok_mask = np.isfinite(metric_vals)

    # Normalize only successful runs for color mapping
    ok_metric = metric_vals[ok_mask]
    norm_vals = np.full_like(metric_vals, np.nan)
    if len(ok_metric) > 0:
        norm_vals[ok_mask] = _safe_normalize(ok_metric)

    cmap = plt.colormaps[cmap_name]

    # Include metric as the final axis
    axis_names = [*param_names, metric]

    # Normalize each axis to [0, 1]
    axis_data: list[NDArray[np.float64]] = []
    for name in axis_names:
        raw = np.array(results.column(name), dtype=np.float64)
        axis_data.append(_safe_normalize(raw))

    # Draw diverged runs first (background, dashed grey)
    for idx in np.where(~ok_mask)[0]:
        coords = [float(axis_data[a][idx]) for a in range(len(axis_names))]
        # Replace NaN coordinate (metric axis) with 0 so line is visible
        coords = [0.0 if np.isnan(c) else c for c in coords]
        ax.plot(
            range(len(axis_names)),
            coords,
            color="0.6",
            alpha=0.3,
            linewidth=0.6,
            linestyle="--",
        )

    # Draw successful runs sorted by metric (low first, high on top)
    ok_indices = np.where(ok_mask)[0]
    for idx in ok_indices[np.argsort(metric_vals[ok_mask])]:
        coords = [float(axis_data[a][idx]) for a in range(len(axis_names))]
        nv = float(norm_vals[idx])
        emphasis = abs(nv - 0.5) * 2  # 0 at median, 1 at extremes
        ax.plot(
            range(len(axis_names)),
            coords,
            color=cmap(nv),
            alpha=0.15 + 0.6 * emphasis,
            linewidth=0.8 + 0.8 * emphasis,
        )

    ax.set_xticks(range(len(axis_names)))
    ax.set_xticklabels(axis_names, rotation=30, ha="right")
    ax.set_ylabel("Normalized value")
    ax.set_title(f"Parallel Coordinates (color = {metric})")

    from matplotlib.cm import ScalarMappable

    if len(ok_metric) > 0:
        cb_norm = _build_norm(ok_metric, divergent_center=divergent_center)
    else:
        from matplotlib.colors import Normalize

        cb_norm = Normalize(0, 1)
    sm = ScalarMappable(cmap=cmap, norm=cb_norm)
    sm.set_array(ok_metric)
    ax.figure.colorbar(sm, ax=ax, label=metric)  # type: ignore[union-attr]


def render_sweep_tornado(
    ax: Axes,
    results: SweepResults,
    metric: str,
) -> None:
    """Tornado chart: horizontal bars showing parameter impact on metric.

    For scattered data (LHS/Sobol), shows Spearman rank correlation.
    For grid-aligned data, shows min-to-max metric range per parameter.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Sweep data.
    metric : str
        Metric to analyze.

    Raises
    ------
    ValueError
        If no swept parameters found.
    """
    param_names = list(results.swept_params.keys())
    if not param_names:
        msg = "sweep-tornado requires at least 1 swept parameter"
        raise ValueError(msg)

    metric_vals = np.array(results.column(metric), dtype=np.float64)

    if _is_scattered_data(results):
        _render_tornado_correlation(ax, results, param_names, metric_vals, metric)
    else:
        _render_tornado_range(ax, results, param_names, metric_vals, metric)


def _render_tornado_correlation(
    ax: Axes,
    results: SweepResults,
    param_names: list[str],
    metric_vals: NDArray[np.float64],
    metric: str,
) -> None:
    """Tornado chart using Spearman rank correlation (for scattered data).

    Bars are colored by statistical significance:
    - Deep red (p < 0.001): highly significant
    - Orange (p < 0.05): significant
    - Grey: not significant

    Significance stars are shown next to bar labels.
    """
    from scipy.stats import spearmanr  # type: ignore[import-untyped]

    # (name, |corr|, sign, p-value)
    impacts: list[tuple[str, float, float, float]] = []
    for name in param_names:
        param_vals = np.array(results.column(name), dtype=np.float64)
        valid = np.isfinite(param_vals) & np.isfinite(metric_vals)
        if valid.sum() < 3:  # noqa: PLR2004
            impacts.append((name, 0.0, 0.0, 1.0))
            continue
        result = spearmanr(param_vals[valid], metric_vals[valid])
        corr = float(result.statistic)  # type: ignore[union-attr]
        pval = float(result.pvalue)  # type: ignore[union-attr]
        impacts.append((name, abs(corr), float(np.sign(corr)), pval))

    impacts.sort(key=lambda t: t[1], reverse=True)  # noqa: FURB118

    names = [t[0] for t in impacts]
    abs_corrs = [t[1] for t in impacts]
    signs = [t[2] for t in impacts]
    pvals = [t[3] for t in impacts]
    y_pos = np.arange(len(names))

    # Color by p-value significance level
    colors = [
        "#D84315" if p < 0.001 else "#FFA726" if p < 0.05 else "#BDBDBD"  # noqa: PLR2004
        for p in pvals
    ]
    bar_container = ax.barh(
        y_pos,
        abs_corrs,
        color=colors,
        alpha=0.85,
        height=0.6,
        edgecolor="white",
        linewidth=0.5,
    )

    for rect, corr, sign, pval in zip(
        bar_container.patches,
        abs_corrs,
        signs,
        pvals,
        strict=True,
    ):
        # Significance stars
        sig = (
            "***"
            if pval < 0.001  # noqa: PLR2004
            else "**"
            if pval < 0.01  # noqa: PLR2004
            else "*"
            if pval < 0.05  # noqa: PLR2004
            else ""
        )
        label = f"{corr * sign:+.2f}{sig}"
        ax.text(
            float(rect.get_width()) + 0.02,
            float(rect.get_y()) + float(rect.get_height()) / 2,
            label,
            va="center",
            fontsize=8,
            fontweight="bold" if sig else "normal",
            color="#D84315" if sig else "#666666",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels([_greek_label(n) for n in names])
    ax.set_xlabel(f"|Spearman r| with {metric}")
    ax.set_xlim(0, 1.15)
    ax.set_title(f"Parameter Sensitivity ({metric})")
    ax.invert_yaxis()

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor="#D84315", label="p < 0.001"),
            Patch(facecolor="#FFA726", label="p < 0.05"),
            Patch(facecolor="#BDBDBD", label="not significant"),
        ],
        loc="lower right",
        fontsize=8,
    )


def _render_tornado_range(
    ax: Axes,
    results: SweepResults,
    param_names: list[str],
    metric_vals: NDArray[np.float64],
    metric: str,
) -> None:
    """Tornado chart using median-of-unique-values (for grid-aligned data)."""
    global_median = float(np.median(metric_vals))

    impacts: list[tuple[str, float, float]] = []
    for name in param_names:
        param_vals = np.array(results.column(name), dtype=np.float64)
        unique_vals = np.unique(param_vals)
        if len(unique_vals) < 2:  # noqa: PLR2004
            impacts.append((name, global_median, global_median))
            continue
        medians: list[float] = []
        for v in unique_vals:
            mask = np.isclose(param_vals, v)
            medians.append(float(np.median(metric_vals[mask])))
        impacts.append((name, min(medians), max(medians)))

    impacts.sort(key=lambda t: t[2] - t[1], reverse=True)

    names = [t[0] for t in impacts]
    lows = [t[1] for t in impacts]
    highs = [t[2] for t in impacts]
    y_pos = np.arange(len(names))

    ax.barh(
        y_pos,
        [hi - lo for lo, hi in zip(lows, highs, strict=True)],
        left=lows,
        color="tab:blue",
        alpha=0.7,
        height=0.6,
    )
    ax.axvline(global_median, color="gray", linestyle="--", linewidth=0.8)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.set_xlabel(metric)
    ax.set_title(f"Parameter Impact on {metric}")
    ax.invert_yaxis()


def render_sweep_scatter(  # noqa: C901, PLR0912, PLR0913, PLR0915
    fig: Figure,
    results: SweepResults,
    metric: str,
    *,
    cmap_name: str = "viridis",
    divergent_center: float | None = None,
    log_diagonal: bool = False,
    params: list[str] | None = None,
    clamp_color: tuple[float, float] | None = None,
) -> None:
    """Pairwise scatter matrix with marginal param-vs-metric on diagonal.

    Diagonal cells show each parameter vs the metric (marginal relationship).
    Off-diagonal cells show pairwise parameter scatter, colored by metric.

    Diverged runs (NaN metric) are plotted as grey 'x' markers so they
    still map the instability boundary without dominating the color scale.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure.
    results : SweepResults
        Sweep data with 2+ swept parameters.
    metric : str
        Metric for coloring.
    cmap_name : str
        Matplotlib colormap name (default ``"viridis"``).
    divergent_center : float or None
        Center colormap at this value using ``TwoSlopeNorm``.
    clamp_color : tuple of (vmin, vmax) or None
        Clamp metric values to this range for coloring.
    log_diagonal : bool
        Use logarithmic y-axis on diagonal panels.
    params : list[str] or None
        Subset of swept parameters to include.  If *None*, use all.

    Raises
    ------
    ValueError
        If fewer than 2 (selected) swept parameters.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable

    all_param_names = list(results.swept_params.keys())
    param_names = (
        [p for p in params if p in all_param_names] if params else all_param_names
    )
    n = len(param_names)
    if n < 2:  # noqa: PLR2004
        msg = f"sweep-scatter requires 2+ swept parameters, got {n}"
        raise ValueError(msg)

    metric_vals = np.array(results.column(metric), dtype=np.float64)
    ok_mask = np.isfinite(metric_vals)
    # Clamp metric values before norm computation so outliers don't
    # dominate the colorbar (e.g. a single divergent A=250 washes out
    # the 0.5-1.1 range where the physics lives).
    if clamp_color is not None:
        vmin_c, vmax_c = clamp_color
        metric_vals = np.where(
            np.isfinite(metric_vals),
            np.clip(metric_vals, vmin_c, vmax_c),
            metric_vals,
        )
    ok_metric = metric_vals[ok_mask]

    cmap = plt.colormaps[cmap_name]
    if len(ok_metric) > 0:
        norm = _build_norm(ok_metric, divergent_center=divergent_center)
    else:
        from matplotlib.colors import Normalize

        norm = Normalize(0, 1)

    # Pre-extract all parameter arrays
    param_data = {p: np.array(results.column(p), dtype=np.float64) for p in param_names}

    axes = fig.subplots(n, n, squeeze=False)

    for i in range(n):
        for j in range(n):
            ax = axes[i][j]
            if i == j:
                # Diagonal: marginal scatter (param vs metric)
                vals = param_data[param_names[i]]
                # Successful runs: colored by metric
                if np.any(ok_mask):
                    ax.scatter(
                        vals[ok_mask],
                        metric_vals[ok_mask],
                        c=metric_vals[ok_mask],
                        cmap=cmap,
                        norm=norm,
                        s=12,
                        alpha=0.7,
                        edgecolors="none",
                    )
                # Diverged runs: grey x markers
                if np.any(~ok_mask):
                    ax.scatter(
                        vals[~ok_mask],
                        np.zeros(np.sum(~ok_mask)),
                        marker="x",
                        c="0.6",
                        s=15,
                        alpha=0.5,
                        zorder=5,
                    )
                if log_diagonal:
                    ax.set_yscale("log")
                # GR baseline for amplification-like metrics
                if metric in {"A", "log10_A", "A_paired", "log10_A_paired", "P_max"}:
                    baseline = 0.0 if metric == "log10_A" else 1.0
                    ax.axhline(
                        baseline,
                        color="#2E7D32",
                        linestyle="--",
                        alpha=0.5,
                        linewidth=1,
                    )
                ax.set_ylabel(metric, fontsize=7)
            else:
                # Off-diagonal: pairwise parameter scatter
                x = param_data[param_names[j]]
                y = param_data[param_names[i]]
                # Successful runs: colored by metric
                if np.any(ok_mask):
                    ax.scatter(
                        x[ok_mask],
                        y[ok_mask],
                        c=metric_vals[ok_mask],
                        cmap=cmap,
                        norm=norm,
                        s=10,
                        alpha=0.7,
                        edgecolors="none",
                    )
                # Diverged runs: grey x markers
                if np.any(~ok_mask):
                    ax.scatter(
                        x[~ok_mask],
                        y[~ok_mask],
                        marker="x",
                        c="0.6",
                        s=15,
                        alpha=0.5,
                        zorder=5,
                    )

            if i == n - 1:
                ax.set_xlabel(_greek_label(param_names[j]), fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0 and i != j:
                ax.set_ylabel(_greek_label(param_names[i]), fontsize=8)
            elif j != 0:
                ax.set_yticklabels([])

    fig.suptitle(f"Scatter Matrix (color = {metric})", fontsize=12)

    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array(ok_metric)
    fig.colorbar(sm, ax=axes[:, -1].ravel().tolist(), label=metric, shrink=0.8)


# -- Replicate convergence diagnostic ----------------------------------------


def render_replicate_convergence(
    ax: Axes,
    results: SweepResults,
    metric: str,
) -> None:
    """Plot running SEM vs replicate count (convergence diagnostic).

    For each parameter point, computes the running mean and SEM as
    replicates are added (k = 1, 2, ..., N). The SEM should decrease
    as 1/sqrt(k). A theoretical reference line is shown for comparison.

    This is the standard "have I run enough replicates?" diagnostic
    used in Monte Carlo simulation.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Ensemble sweep data (must have replicates).
    metric : str
        Metric column to analyze.

    References
    ----------
    Smith, R.C. (2013) *Uncertainty Quantification*, SIAM. Ch. 3.2
    (standard error convergence).
    """
    import matplotlib.pyplot as plt

    groups = results.group_by_point()
    param_names = list(results.swept_params.keys())

    cmap = plt.colormaps["viridis"]
    n_groups = len(groups)

    for gi, (point_key, rep_rows) in enumerate(groups.items()):
        # Extract metric values in replicate order
        sorted_rows = sorted(rep_rows, key=lambda r: r.get("replicate", 0))
        vals = [
            r[metric]
            for r in sorted_rows
            if r.get(metric) is not None and r.get("run_status") == "success"
        ]
        if len(vals) < 2:  # noqa: PLR2004
            continue

        arr = np.array(vals, dtype=np.float64)
        n_reps = len(arr)
        ks = np.arange(2, n_reps + 1)  # SEM undefined for k=1

        # Running SEM: std(first k values) / sqrt(k), for k=2..N
        running_sem = np.array(
            [float(np.std(arr[:k], ddof=1) / np.sqrt(k)) for k in ks],
            dtype=np.float64,
        )

        # Label from parameter point
        label_parts = [
            f"{p}={v:.3g}" for p, v in zip(param_names, point_key, strict=True)
        ]
        label = ", ".join(label_parts)

        color = cmap(gi / max(n_groups - 1, 1))
        ax.plot(
            ks,
            running_sem,
            "o-",
            color=color,
            markersize=4,
            linewidth=1.2,
            label=label,
        )

    # Reference: 1/sqrt(k) scaled to typical SEM at k=2
    # (shows expected convergence rate)
    if groups:
        all_sems: list[float] = []
        for rep_rows in groups.values():
            sorted_rows = sorted(rep_rows, key=lambda r: r.get("replicate", 0))
            vals = [
                r[metric]
                for r in sorted_rows
                if r.get(metric) is not None and r.get("run_status") == "success"
            ]
            if len(vals) >= 2:  # noqa: PLR2004
                arr = np.array(vals[:2], dtype=np.float64)
                all_sems.append(float(np.std(arr, ddof=1) / np.sqrt(2)))
        if all_sems:
            sem_at_2 = np.median(all_sems)
            max_k = max(
                len(
                    [
                        r
                        for r in rr
                        if r.get(metric) is not None
                        and r.get("run_status") == "success"
                    ],
                )
                for rr in groups.values()
            )
            k_ref = np.arange(2, max_k + 1)
            sem_ref = sem_at_2 * np.sqrt(2) / np.sqrt(k_ref)
            ax.plot(
                k_ref,
                sem_ref,
                "--",
                color="gray",
                alpha=0.6,
                linewidth=1.5,
                label=r"$\propto 1/\sqrt{k}$",
            )

    ax.set_xlabel("Number of replicates (k)")
    ax.set_ylabel(f"SEM({metric})")
    ax.set_title(f"Replicate Convergence — {metric}")
    ax.legend(fontsize=7, loc="best")
    ax.grid(visible=True, alpha=0.3)


# -- Histogram / distribution -------------------------------------------------


def render_sweep_histogram(
    ax: Axes,
    results: SweepResults,
    metric: str,
) -> None:
    """Histogram of metric values with summary statistics.

    Shows distribution of the metric across all successful runs,
    with a GR baseline reference line and a statistics annotation box.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Sweep data.
    metric : str
        Metric column to histogram.
    """
    metric_vals = np.array(results.column(metric), dtype=np.float64)
    valid = metric_vals[np.isfinite(metric_vals)]

    if len(valid) == 0:
        ax.text(
            0.5,
            0.5,
            f"No valid data for {metric}",
            transform=ax.transAxes,
            ha="center",
        )
        return

    ax.hist(
        valid,
        bins=40,
        color="#1976D2",
        alpha=0.8,
        edgecolor="white",
        linewidth=0.5,
    )

    # GR baseline reference
    if metric == "log10_A":
        ax.axvline(0, color="#2E7D32", linestyle="--", linewidth=2, label="A = 1 (GR)")
    elif metric == "A":
        ax.axvline(1, color="#2E7D32", linestyle="--", linewidth=2, label="A = 1 (GR)")

    ax.set_xlabel(metric)
    ax.set_ylabel("Count")
    ax.set_title(f"Distribution of {metric} (n={len(valid)})")

    # Statistics box
    n_total = len(valid)
    if metric in {"A", "log10_A"}:
        a_vals = 10.0**valid if metric == "log10_A" else valid
        n_near_gr = int(np.sum((a_vals > 0.9) & (a_vals < 1.1)))  # noqa: PLR2004
        n_gt2 = int(np.sum(a_vals > 2))  # noqa: PLR2004
        n_gt10 = int(np.sum(a_vals > 10))  # noqa: PLR2004
        stats_text = (
            f"Near GR: {n_near_gr} ({100 * n_near_gr / n_total:.0f}%)\n"
            f"A > 2: {n_gt2} ({100 * n_gt2 / n_total:.0f}%)\n"
            f"A > 10: {n_gt10} ({100 * n_gt10 / n_total:.0f}%)"
        )
    else:
        stats_text = (
            f"median: {np.median(valid):.3g}\n"
            f"mean: {np.mean(valid):.3g}\n"
            f"std: {np.std(valid):.3g}"
        )
    ax.text(
        0.97,
        0.95,
        stats_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "#E3F2FD", "alpha": 0.8},
    )

    ax.legend(fontsize=10)


# -- Divergence rate -----------------------------------------------------------


def render_sweep_divergence(  # noqa: C901
    ax: Axes,
    results: SweepResults,
) -> None:
    """Horizontal bar chart of divergence rate per swept parameter.

    For each swept parameter, computes the fraction of runs that
    diverged when that parameter was non-zero (or varied).

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Sweep data.
    """
    param_names = list(results.swept_params.keys())
    if not param_names:
        ax.text(0.5, 0.5, "No swept parameters", transform=ax.transAxes, ha="center")
        return

    names: list[str] = []
    diverge_rates: list[float] = []

    for name in param_names:
        param_vals = np.array(results.column(name), dtype=np.float64)
        unique_vals = np.unique(param_vals[np.isfinite(param_vals)])
        if len(unique_vals) < 2:  # noqa: PLR2004
            continue

        n_total = 0
        n_diverged = 0
        for row in results.rows:
            val = row.get(name)
            if val is None:
                continue
            n_total += 1
            if row.get("run_status", "success") != "success":
                n_diverged += 1

        rate = 100.0 * n_diverged / n_total if n_total > 0 else 0.0
        names.append(_greek_label(name))
        diverge_rates.append(rate)

    if not names:
        ax.text(0.5, 0.5, "No varying parameters", transform=ax.transAxes, ha="center")
        return

    # Sort by divergence rate
    order = sorted(range(len(names)), key=lambda i: diverge_rates[i])
    names = [names[i] for i in order]
    diverge_rates = [diverge_rates[i] for i in order]

    y_pos = np.arange(len(names))
    colors = [
        "#D32F2F" if r > 50 else "#FFA726" if r > 10 else "#1976D2"  # noqa: PLR2004
        for r in diverge_rates
    ]

    bars = ax.barh(
        y_pos,
        diverge_rates,
        color=colors,
        alpha=0.85,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel("Divergence Rate (%)")
    ax.set_xlim(0, max(105, max(diverge_rates) * 1.1))
    ax.set_title("Divergence Rate by Parameter")

    # Add percentage labels
    for bar, rate in zip(bars, diverge_rates, strict=False):  # type: ignore[reportUnknownVariableType]
        if rate > 5:  # noqa: PLR2004
            ax.text(
                bar.get_width() - 3,  # type: ignore[reportUnknownMemberType]
                bar.get_y() + bar.get_height() / 2,  # type: ignore[reportUnknownMemberType]
                f"{rate:.0f}%",
                ha="right",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )
        elif rate > 0:
            ax.text(
                bar.get_width() + 1,  # type: ignore[reportUnknownMemberType]
                bar.get_y() + bar.get_height() / 2,  # type: ignore[reportUnknownMemberType]
                f"{rate:.1f}%",
                ha="left",
                va="center",
                fontsize=7,
                color="#333",
            )
