"""Visualization helpers for ``tg simulate``."""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from pathlib import Path

    from matplotlib.axes import Axes
    from pde import FieldCollection

from torsion_gertsenshtein.cli._simulate import (
    DPI,
    SPATIAL_DIM_2D,
    VMAX_FLOOR,
    PlotContext,
    field_slots,
)

# --- Shared plot helpers ---


def plot_amplitude_decay(
    ax: Axes,
    ctx: PlotContext,
    slot: int,
    name: str,
    times: list[float],
) -> None:
    """Plot peak amplitude decay over time (shared by 1D/2D/3D)."""
    peaks = [float(np.max(np.abs(ctx.initial_state[slot].data)))]
    for t_idx in range(len(ctx.storage)):
        snap = cast("FieldCollection", ctx.storage[t_idx])
        peaks.append(float(np.max(np.abs(snap[slot].data))))
    ax.plot([0.0, *times], peaks, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(f"max |{name}|")
    ax.set_title("Peak amplitude")
    ax.grid(visible=True, alpha=0.3)


def plot_1d(path: Path, ctx: PlotContext) -> None:
    """1D visualization: spacetime heatmap + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    times = list(ctx.storage.times)
    x_coords = cast("np.ndarray", ctx.grid.cell_coords[..., 0])
    slots = field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    # Panel 1: spacetime heatmap for first field
    ax = axes[0]
    heatmap_data: list[np.ndarray] = []
    for t_idx in range(len(ctx.storage)):
        snap = cast("FieldCollection", ctx.storage[t_idx])
        heatmap_data.append(snap[slot].data)

    heatmap = np.array(heatmap_data)
    ax.imshow(
        heatmap,
        aspect="auto",
        origin="lower",
        extent=[float(x_coords[0]), float(x_coords[-1]), 0, times[-1]],
        cmap="RdBu_r",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(f"{name} spacetime evolution")

    # Panel 2: amplitude decay
    plot_amplitude_decay(axes[1], ctx, slot, name, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(
        f"{ctx.spec.component_names[0]} ({param_str})"
        if ctx.params
        else ctx.spec.component_names[0]
    )
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


def plot_2d(path: Path, ctx: PlotContext) -> None:
    """2D visualization: initial + final snapshots + amplitude decay."""
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    times = list(ctx.storage.times)
    final = cast("FieldCollection", ctx.storage[-1])
    slots = field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    bounds = ctx.grid.axes_bounds

    # Panel 1: initial x-y
    ax = axes[0]
    init_data = ctx.initial_state[slot].data
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
    final_data = final[slot].data
    vmax_f = max(float(np.max(np.abs(final_data))), VMAX_FLOOR)
    ax.imshow(
        final_data.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax_f,
        vmax=vmax_f,
    )
    ax.set_title(f"{name} (t={times[-1]:.1f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    # Panel 3: amplitude decay
    plot_amplitude_decay(axes[2], ctx, slot, name, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(f"{name} ({param_str})" if ctx.params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


# --- 3D plot sub-helpers ---


def _plot_z_profile(
    ax: Axes,
    ctx: PlotContext,
    slot: int,
    name: str,
    times: list[float],
) -> None:
    """Plot z-profile at center of x-y plane."""
    n = ctx.grid.shape[0]
    ic = n // 2
    z_1d = cast("np.ndarray", ctx.grid.cell_coords[ic, ic, :, 2])
    final = cast("FieldCollection", ctx.storage[-1])

    ax.plot(z_1d, ctx.initial_state[slot].data[ic, ic, :], "b-", linewidth=2, label="t=0")
    ax.plot(
        z_1d, final[slot].data[ic, ic, :], "r-", linewidth=2, label=f"t={times[-1]:.1f}"
    )
    ax.set_xlabel("z")
    ax.set_ylabel(name)
    ax.set_title(f"{name} z-profile (x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)


def _plot_xy_slice(
    ax: Axes,
    ctx: PlotContext,
    slot: int,
    name: str,
) -> None:
    """Plot x-y slice at z=center."""
    n = ctx.grid.shape[0]
    ic = n // 2
    bounds = ctx.grid.axes_bounds

    init_slice = ctx.initial_state[slot].data[:, :, ic]
    vmax = max(float(np.max(np.abs(init_slice))), VMAX_FLOOR)
    ax.imshow(
        init_slice.T,
        origin="lower",
        extent=(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]),
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title(f"{name} x-y (t=0, z=center)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")


def _plot_component_independence(
    ax: Axes,
    ctx: PlotContext,
    slots: dict[str, int],
    times: list[float],
) -> None:
    """Plot amplitude of other components over time."""
    colors = ["red", "green", "purple", "orange", "brown"]
    for i in range(1, ctx.spec.n_components):
        comp_name = ctx.spec.component_names[i]
        comp_slot = slots[comp_name]
        comp_peaks = [0.0]
        for t_idx in range(len(ctx.storage)):
            snap = cast("FieldCollection", ctx.storage[t_idx])
            comp_peaks.append(float(np.max(np.abs(snap[comp_slot].data))))
        ax.plot(
            [0.0, *times],
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


def plot_3d(path: Path, ctx: PlotContext) -> None:
    """3D visualization: z-profile + x-y slice + amplitude decay + component check."""
    import matplotlib.pyplot as plt

    times = list(ctx.storage.times)
    slots = field_slots(ctx.spec)
    name = ctx.spec.component_names[0]
    slot = slots[name]

    n_panels = 4 if ctx.spec.n_components > 1 else 3
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

    # Panel 1: z-profile
    _plot_z_profile(axes[0], ctx, slot, name, times)

    # Panel 2: x-y slice initial
    _plot_xy_slice(axes[1], ctx, slot, name)

    # Panel 3: amplitude decay
    plot_amplitude_decay(axes[2], ctx, slot, name, times)

    # Panel 4: component independence (multi-field only)
    if ctx.spec.n_components > 1:
        _plot_component_independence(axes[3], ctx, slots, times)

    param_str = ", ".join(f"{k}={v}" for k, v in ctx.params.items())
    fig.suptitle(f"{name} ({param_str})" if ctx.params else name)
    plt.tight_layout()
    plt.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved plot to: {path}")


def save_plot(path: Path, ctx: PlotContext) -> None:
    """Generate and save visualization — dispatches to 1D/2D/3D."""
    import matplotlib as mpl

    mpl.use("Agg")

    spatial_dim = ctx.spec.spatial_dimension

    if spatial_dim == 1:
        plot_1d(path, ctx)
    elif spatial_dim == SPATIAL_DIM_2D:
        plot_2d(path, ctx)
    else:
        plot_3d(path, ctx)
