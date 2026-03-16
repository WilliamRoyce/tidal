"""Render functions for sweep-specific plot types in ``tidal plot``.

Provides:

- ``render_sweep_1d``   — metric vs single swept parameter (line plot)
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
    from matplotlib.figure import Figure
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData
    from tidal.measurement._sweep_results import SweepResults


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
        result = eval(formula, {"__builtins__": {}}, ns)  # noqa: S307
    except Exception as exc:
        msg = f"Error evaluating overlay formula {formula!r}: {exc}"
        raise ValueError(msg) from exc
    return np.asarray(result, dtype=np.float64)


def render_sweep_1d(
    ax: Axes,
    results: SweepResults,
    metric: str,
    *,
    overlay: str | None = None,
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
        y_mean = agg.column(f"{metric}_mean")
        y_std = agg.column(f"{metric}_std")
        num_label = "TIDAL" if overlay else None
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
    else:
        x = results.column(param_name)
        y = results.column(metric)
        num_label = "TIDAL" if overlay else None
        ax.plot(
            x, y, "o-", color="tab:blue", linewidth=1.5, markersize=5, label=num_label
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

    ax.set_xlabel(param_name)
    ax.set_ylabel(metric)
    ax.set_title(f"{metric} vs {param_name}")
    ax.grid(visible=True, alpha=0.3)


def _plot_1d_metric(
    ax: Axes,
    results: SweepResults,
    param_name: str,
    metric: str,
    color: str,
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


def render_sweep_1d_multi(
    fig: Figure,
    results: SweepResults,
    metrics: list[str],
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
        _plot_1d_metric(ax, results, param_name, metric, colors[i % len(colors)])
        ax.set_ylabel(metric)
        ax.grid(visible=True, alpha=0.3)

    axes[-1, 0].set_xlabel(param_name)
    fig.suptitle(f"Sweep: {', '.join(metrics)} vs {param_name}")


def render_sweep_2d(
    ax: Axes,
    results: SweepResults,
    metric: str,
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
        _render_2d_scattered(ax, results, p1_name, p2_name, metric)
    else:
        _render_2d_grid(ax, results, p1_name, p2_name, metric)


def _render_2d_scattered(
    ax: Axes,
    results: SweepResults,
    p1_name: str,
    p2_name: str,
    metric: str,
) -> None:
    """Render 2D sweep as scatter + interpolation background."""
    p1 = np.array(results.column(p1_name), dtype=np.float64)
    p2 = np.array(results.column(p2_name), dtype=np.float64)
    metric_vals = np.array(results.column(metric), dtype=np.float64)

    valid = np.isfinite(metric_vals)
    p1, p2, metric_vals = p1[valid], p2[valid], metric_vals[valid]

    if len(metric_vals) == 0:
        ax.text(
            0.5, 0.5, f"No valid data for {metric}", transform=ax.transAxes, ha="center"
        )
        return

    vmin, vmax = float(metric_vals.min()), float(metric_vals.max())

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
                cmap="viridis",
                alpha=0.3,
                vmin=vmin,
                vmax=vmax,
            )
        except (ValueError, ImportError):
            pass  # Fall back to scatter-only

    sc = ax.scatter(
        p1,
        p2,
        c=metric_vals,
        cmap="viridis",
        s=60,
        edgecolors="k",
        linewidths=0.5,
        vmin=vmin,
        vmax=vmax,
        zorder=5,
    )
    ax.set_xlabel(p1_name)
    ax.set_ylabel(p2_name)
    ax.set_title(metric)
    ax.figure.colorbar(sc, ax=ax, label=metric)  # type: ignore[union-attr]


def _render_2d_grid(
    ax: Axes,
    results: SweepResults,
    p1_name: str,
    p2_name: str,
    metric: str,
) -> None:
    """Render 2D sweep as pcolormesh heatmap (for grid-aligned data)."""
    p1_vals = np.sort(results.swept_params[p1_name])
    p2_vals = np.sort(results.swept_params[p2_name])

    n1, n2 = len(p1_vals), len(p2_vals)
    grid = np.full((n2, n1), np.nan)

    for row in results.rows:
        v1 = row.get(p1_name)
        v2 = row.get(p2_name)
        val = row.get(metric)
        if v1 is None or v2 is None or val is None:
            continue
        i1 = min(int(np.searchsorted(p1_vals, float(v1))), n1 - 1)
        i2 = min(int(np.searchsorted(p2_vals, float(v2))), n2 - 1)
        grid[i2, i1] = float(val)

    im = ax.pcolormesh(
        p1_vals,
        p2_vals,
        grid,
        shading="nearest",
        cmap="viridis",
    )
    ax.set_xlabel(p1_name)
    ax.set_ylabel(p2_name)
    ax.set_title(metric)
    ax.figure.colorbar(im, ax=ax, label=metric)  # type: ignore[union-attr]


def render_sweep_2d_with_overlay(  # noqa: PLR0914
    fig: Figure,
    results: SweepResults,
    metric: str,
    overlay: str,
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
    vmin = float(valid.min()) if len(valid) else 0.0
    vmax = float(valid.max()) if len(valid) else 1.0

    axes = fig.subplots(1, 3)

    # Panel 0: TIDAL numerical
    ax0 = axes[0]
    im0 = ax0.pcolormesh(
        p1_vals, p2_vals, z_num, shading="nearest", cmap="viridis", vmin=vmin, vmax=vmax
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
        vmin=vmin,
        vmax=vmax,
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
            0.5, 0.5, f"No valid data for {metric}", transform=ax.transAxes, ha="center"
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
            n, y, "o-", color="tab:blue", linewidth=1.5, markersize=6, label="measured"
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
    """Normalize values to [0, 1] range for colormap mapping."""
    vmin, vmax = values.min(), values.max()
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
) -> None:
    """Parallel coordinates: each axis is a parameter + metric, color = metric.

    Lines are sorted by metric value (low-to-high) so high-metric runs
    draw on top.  Alpha and linewidth scale with distance from median to
    emphasize extreme (most interesting) values.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    results : SweepResults
        Sweep data with 2+ swept parameters.
    metric : str
        Metric for coloring the polylines.

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
    norm = _safe_normalize(metric_vals)
    cmap = plt.colormaps["viridis"]

    # Include metric as the final axis
    axis_names = [*param_names, metric]

    # Normalize each axis to [0, 1]
    axis_data: list[NDArray[np.float64]] = []
    for name in axis_names:
        raw = np.array(results.column(name), dtype=np.float64)
        axis_data.append(_safe_normalize(raw))

    # Sort by metric: draw low values first, high values on top
    for idx in np.argsort(metric_vals):
        coords = [axis_data[a][idx] for a in range(len(axis_names))]
        emphasis = abs(float(norm[idx]) - 0.5) * 2  # 0 at median, 1 at extremes
        ax.plot(
            range(len(axis_names)),
            coords,
            color=cmap(norm[idx]),
            alpha=0.15 + 0.6 * emphasis,
            linewidth=0.8 + 0.8 * emphasis,
        )

    ax.set_xticks(range(len(axis_names)))
    ax.set_xticklabels(axis_names, rotation=30, ha="right")
    ax.set_ylabel("Normalized value")
    ax.set_title(f"Parallel Coordinates (color = {metric})")

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    sm = ScalarMappable(
        cmap=cmap,
        norm=Normalize(vmin=float(metric_vals.min()), vmax=float(metric_vals.max())),
    )
    sm.set_array(metric_vals)
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
    """Tornado chart using Spearman rank correlation (for scattered data)."""
    from scipy.stats import spearmanr  # type: ignore[import-untyped]

    impacts: list[tuple[str, float, float]] = []  # (name, |corr|, sign)
    for name in param_names:
        param_vals = np.array(results.column(name), dtype=np.float64)
        valid = np.isfinite(param_vals) & np.isfinite(metric_vals)
        if valid.sum() < 3:  # noqa: PLR2004
            impacts.append((name, 0.0, 0.0))
            continue
        result = spearmanr(param_vals[valid], metric_vals[valid])
        corr = float(result.statistic)  # type: ignore[union-attr]
        impacts.append((name, abs(corr), float(np.sign(corr))))

    impacts.sort(key=lambda t: t[1], reverse=True)  # noqa: FURB118

    names = [t[0] for t in impacts]
    abs_corrs = [t[1] for t in impacts]
    signs = [t[2] for t in impacts]
    y_pos = np.arange(len(names))

    colors = ["tab:blue" if s >= 0 else "tab:red" for s in signs]
    bar_container = ax.barh(y_pos, abs_corrs, color=colors, alpha=0.7, height=0.6)

    for rect, corr, sign in zip(bar_container.patches, abs_corrs, signs, strict=True):
        label = f"{corr * sign:+.2f}"
        ax.text(
            float(rect.get_width()) + 0.02,
            float(rect.get_y()) + float(rect.get_height()) / 2,
            label,
            va="center",
            fontsize=8,
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names)
    ax.set_xlabel(f"|Spearman r| with {metric}")
    ax.set_xlim(0, 1.15)
    ax.set_title(f"Parameter Sensitivity ({metric})")
    ax.invert_yaxis()

    from matplotlib.patches import Patch

    ax.legend(
        handles=[
            Patch(facecolor="tab:blue", alpha=0.7, label="positive"),
            Patch(facecolor="tab:red", alpha=0.7, label="negative"),
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


def render_sweep_scatter(
    fig: Figure,
    results: SweepResults,
    metric: str,
) -> None:
    """Pairwise scatter matrix with marginal param-vs-metric on diagonal.

    Diagonal cells show each parameter vs the metric (marginal relationship).
    Off-diagonal cells show pairwise parameter scatter, colored by metric.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure.
    results : SweepResults
        Sweep data with 2+ swept parameters.
    metric : str
        Metric for coloring.

    Raises
    ------
    ValueError
        If fewer than 2 swept parameters.
    """
    import matplotlib.pyplot as plt

    param_names = list(results.swept_params.keys())
    n = len(param_names)
    if n < 2:  # noqa: PLR2004
        msg = f"sweep-scatter requires 2+ swept parameters, got {n}"
        raise ValueError(msg)

    metric_vals = np.array(results.column(metric), dtype=np.float64)
    cmap = plt.colormaps["viridis"]
    vmin, vmax = float(np.nanmin(metric_vals)), float(np.nanmax(metric_vals))

    axes = fig.subplots(n, n, squeeze=False)

    for i in range(n):
        for j in range(n):
            ax = axes[i][j]
            if i == j:
                # Diagonal: marginal scatter (param vs metric)
                vals = np.array(results.column(param_names[i]), dtype=np.float64)
                ax.scatter(
                    vals,
                    metric_vals,
                    c=metric_vals,
                    cmap=cmap,
                    s=12,
                    alpha=0.7,
                    vmin=vmin,
                    vmax=vmax,
                )
                ax.set_ylabel(metric, fontsize=7)
            else:
                # Off-diagonal: pairwise parameter scatter
                x = np.array(results.column(param_names[j]), dtype=np.float64)
                y = np.array(results.column(param_names[i]), dtype=np.float64)
                ax.scatter(
                    x,
                    y,
                    c=metric_vals,
                    cmap=cmap,
                    s=10,
                    alpha=0.7,
                    vmin=vmin,
                    vmax=vmax,
                )

            if i == n - 1:
                ax.set_xlabel(param_names[j], fontsize=8)
            else:
                ax.set_xticklabels([])
            if j == 0 and i != j:
                ax.set_ylabel(param_names[i], fontsize=8)
            elif j != 0:
                ax.set_yticklabels([])

    fig.suptitle(f"Scatter Matrix (color = {metric})", fontsize=12)

    # Colorbar beside rightmost column only (avoids overlap)
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    sm = ScalarMappable(
        cmap=cmap,
        norm=Normalize(vmin=vmin, vmax=vmax),
    )
    sm.set_array(metric_vals)
    fig.colorbar(sm, ax=axes[:, -1].ravel().tolist(), label=metric, shrink=0.8)


# -- Replicate convergence diagnostic ----------------------------------------


def render_replicate_convergence(  # noqa: PLR0914
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
            ks, running_sem, "o-", color=color, markersize=4, linewidth=1.2, label=label
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
                    ]
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
