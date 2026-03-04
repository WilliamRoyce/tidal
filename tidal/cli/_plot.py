"""Visualization helpers for ``tidal simulate``."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from tidal.cli._panels import single_field

if TYPE_CHECKING:
    from pathlib import Path

    from matplotlib.axes import Axes

    from tidal.measurement._io import SimulationData
    from tidal.solver.grid import GridInfo

from tidal.cli._simulate import (
    DPI,
    SPATIAL_DIM_2D,
    VMAX_FLOOR,
)

# --- Shared helpers ---


def _amplitude_peaks(field_history: np.ndarray) -> np.ndarray:
    """Vectorized peak amplitude per snapshot.

    Parameters
    ----------
    field_history : ndarray, shape (n_snapshots, *grid_shape)

    Returns
    -------
    ndarray, shape (n_snapshots,)
    """
    reduce_axes = tuple(range(1, field_history.ndim))
    return np.max(np.abs(field_history), axis=reduce_axes)


# --- Dimension-specific plot functions ---


def _plot_1d(path: Path, sd: SimulationData, gi: GridInfo) -> None:
    """1D: spacetime heatmap + amplitude decay (vectorized)."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    name = single_field(sd, None)
    field_hist = sd.fields[name]  # (n_snapshots, n_x)
    x = gi.axes_coords(0)
    times = sd.times

    # Panel 1: spacetime heatmap — already a 2D array
    ax = axes[0]
    ax.imshow(
        field_hist,
        aspect="auto",
        origin="lower",
        extent=[float(x[0]), float(x[-1]), float(times[0]), float(times[-1])],
        cmap="RdBu_r",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(f"{name} spacetime evolution")

    # Panel 2: vectorized amplitude decay
    ax = axes[1]
    peaks = _amplitude_peaks(field_hist)
    ax.plot(times, peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    param_str = ", ".join(f"{k}={v}" for k, v in sd.parameters.items())
    fig.suptitle(f"{name} ({param_str})" if sd.parameters else name)
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path.resolve()}")


def _plot_2d(path: Path, sd: SimulationData, gi: GridInfo) -> None:
    """2D: initial + final snapshots + amplitude decay (vectorized)."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    name = single_field(sd, None)
    field_hist = sd.fields[name]  # (n_snapshots, nx, ny)
    times = sd.times
    bounds = gi.bounds

    # Panel 1: initial x-y
    ax = axes[0]
    init_data = field_hist[0]
    vmax = max(float(np.max(np.abs(init_data))), VMAX_FLOOR)
    ax.imshow(
        init_data.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} (t=0)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 2: final x-y
    ax = axes[1]
    final_data = field_hist[-1]
    vmax_f = max(float(np.max(np.abs(final_data))), VMAX_FLOOR)
    ax.imshow(
        final_data.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax_f,
        vmax=vmax_f,
    )
    ax.set_title(f"{name} (t={float(times[-1]):.1f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 3: vectorized amplitude decay
    ax = axes[2]
    peaks = _amplitude_peaks(field_hist)
    ax.plot(times, peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    param_str = ", ".join(f"{k}={v}" for k, v in sd.parameters.items())
    fig.suptitle(f"{name} ({param_str})" if sd.parameters else name)
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path.resolve()}")


def _plot_component_independence(
    ax: Axes,
    sd: SimulationData,
) -> None:
    """Plot amplitude of non-primary components over time (vectorized)."""
    colors = ["red", "green", "purple", "orange", "brown"]
    for i in range(1, sd.spec.n_components):
        comp_name = sd.spec.component_names[i]
        if comp_name not in sd.fields:
            continue
        comp_peaks = _amplitude_peaks(sd.fields[comp_name])
        ax.plot(
            sd.times,
            comp_peaks,
            color=colors[(i - 1) % len(colors)],
            linewidth=2,
            label=comp_name,
        )
    ax.set_xlabel("Time")
    ax.set_ylabel("max |field|")
    ax.set_title("Other components")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))


def _plot_3d(path: Path, sd: SimulationData, gi: GridInfo) -> None:
    """3D: z-profile + x-y slice + amplitude decay + component check."""
    import matplotlib.pyplot as plt

    name = single_field(sd, None)
    field_hist = sd.fields[name]  # (n_snapshots, nx, ny, nz)
    times = sd.times
    bounds = gi.bounds

    n_panels = 4 if sd.spec.n_components > 1 else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

    # Panel 1: z-profile at center
    ax = axes[0]
    ic = gi.shape[0] // 2
    z_1d = gi.axes_coords(2)
    ax.plot(z_1d, field_hist[0, ic, ic, :], "b-", linewidth=2, label="t=0")
    ax.plot(
        z_1d,
        field_hist[-1, ic, ic, :],
        "r-",
        linewidth=2,
        label=f"t={float(times[-1]):.1f}",
    )
    ax.set_xlabel("z")
    ax.set_ylabel(name)
    ax.set_title(f"{name} z-profile (x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # Panel 2: z-t spacetime heatmap at x=y=center
    ax = axes[1]
    zt_slice = field_hist[:, ic, ic, :]  # (n_snapshots, nz)
    vmax = max(float(np.max(np.abs(zt_slice))), VMAX_FLOOR)
    ax.imshow(
        zt_slice.T,
        aspect="auto",
        origin="lower",
        extent=(float(times[0]), float(times[-1]), bounds[2][0], bounds[2][1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} z-t (x=y=center)")
    ax.set_xlabel("Time")
    ax.set_ylabel("z")

    # Panel 3: vectorized amplitude decay
    ax = axes[2]
    peaks = _amplitude_peaks(field_hist)
    ax.plot(times, peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)

    # Panel 4: other components (multi-field only)
    if sd.spec.n_components > 1:
        _plot_component_independence(axes[3], sd)

    param_str = ", ".join(f"{k}={v}" for k, v in sd.parameters.items())
    fig.suptitle(f"{name} ({param_str})" if sd.parameters else name)
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path.resolve()}")


def save_plot(path: Path, sd: SimulationData, gi: GridInfo) -> None:
    """Generate and save visualization from SimulationData."""
    import matplotlib as mpl

    mpl.use("Agg")

    spatial_dim = sd.spec.spatial_dimension

    if spatial_dim == 1:
        _plot_1d(path, sd, gi)
    elif spatial_dim == SPATIAL_DIM_2D:
        _plot_2d(path, sd, gi)
    else:
        _plot_3d(path, sd, gi)
