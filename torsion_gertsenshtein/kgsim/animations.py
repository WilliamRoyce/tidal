"""Unified plotting utilities for Klein-Gordon simulations.

This module provides high-level functions for creating visualizations and animations
from simulation snapshots, reducing code duplication across examples.
"""

from __future__ import annotations

import logging
import pathlib
import shutil
from typing import TYPE_CHECKING, Any, Literal, cast

import matplotlib as mpl

mpl.use("Agg")
import operator

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.colors import Normalize, TwoSlopeNorm
from tqdm import tqdm

if TYPE_CHECKING:
    from pde import CartesianGrid

logger = logging.getLogger(__name__)


def choose_writer_and_out(
    snap_count: int, t_end: float, out_path: str | pathlib.Path, fps: int | None = None
) -> tuple[FFMpegWriter | PillowWriter, str, str]:
    """
    Select an appropriate matplotlib animation writer and output filename.

    The caller provides out_path (without extension) and this function will
    append the appropriate extension (".mp4" for ffmpeg, ".gif" for pillow).

    Parameters
    ----------
    snap_count : int
        Number of frames/snapshots available for the animation.
    t_end : float
        Total simulation time (same units as snapshot spacing). Used to compute a
        target frames-per-second (fps) for playback.
    out_path : str | pathlib.Path
        Base path for output file (extension will be added).
    fps : int | None, optional
        Frames per second. If None, computed from snap_count and t_end.

    Returns
    -------
    tuple[FFMpegWriter | PillowWriter, str, str]
        A tuple containing:
        - writer: an instance of a matplotlib.animation writer (FFMpegWriter or PillowWriter).
        - out: str, path to the output file with extension appended.
        - use_writer: str, identifier of the writer used ("ffmpeg" or "pillow").

    Notes
    -----
    - Prefers ffmpeg if available; falls back to PillowWriter when ffmpeg is not available.
    - Chooses output filename and a string label for which writer was chosen.
    """
    if fps is None:
        fps = max(1, int(snap_count / max(1.0, t_end / 5.0)))

    p = pathlib.Path(out_path)
    base = str(
        p.with_suffix("")
    )  # strip any existing suffix, we'll append the chosen one

    # Check if ffmpeg is available
    if shutil.which("ffmpeg") is not None:
        writer = FFMpegWriter(fps=fps, metadata={"artist": "kgsim"}, bitrate=2000)
        out = base + ".mp4"
        use_writer = "ffmpeg"
    else:
        writer = PillowWriter(fps=fps)
        out = base + ".gif"
        use_writer = "pillow"
    return writer, out, use_writer


def create_spacetime_plot(  # noqa: PLR0913
    snapshots: list[tuple[float, np.ndarray]],
    grid: CartesianGrid,
    output_path: str | pathlib.Path,
    *,
    title: str = r"Klein-Gordon evolution: $\phi(x,t)$",
    xlabel: str = r"$x$",
    ylabel: str = r"$t$",
    cbar_label: str = r"$\phi$",
    cmap: str = "bwr",
    use_twoslope_norm: bool = True,
    dpi: int = 200,
    figsize: tuple[float, float] = (4, 3),
) -> None:
    """
    Create a spacetime heatmap from 1D simulation snapshots.

    Parameters
    ----------
    snapshots : list[tuple[float, np.ndarray]]
        List of (time, field_data) tuples where field_data is 1D.
    grid : CartesianGrid
        Grid providing spatial bounds.
    output_path : str | pathlib.Path
        Output file path.
    title : str, optional
        Plot title.
    xlabel : str, optional
        X-axis label (spatial coordinate).
    ylabel : str, optional
        Y-axis label (time coordinate).
    cbar_label : str, optional
        Colorbar label.
    cmap : str, optional
        Colormap name.
    use_twoslope_norm : bool, optional
        If True, use TwoSlopeNorm centered at zero.
    dpi : int, optional
        Output DPI.
    figsize : tuple[float, float], optional
        Figure size (width, height) in inches.

    Raises
    ------
    ValueError
        If no snapshots provided.
    """
    if not snapshots:
        msg = "No snapshots provided"
        raise ValueError(msg)

    # Sort by time
    snapshots = sorted(snapshots, key=operator.itemgetter(0))

    # Build evolution array: shape (nt, nx)
    times = [t for t, _ in snapshots]
    phi_rows = [
        arr.reshape(grid.shape) if hasattr(grid, "shape") else arr.ravel()
        for _, arr in snapshots
    ]
    data = np.vstack([row.ravel() for row in phi_rows])

    # Setup colormap normalization
    vmin = float(data.min())
    vmax = float(data.max())

    if use_twoslope_norm:
        # Ensure symmetric range for TwoSlopeNorm
        abs_max = max(abs(vmin), abs(vmax))
        norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Create plot
    pathlib.Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        data,
        aspect="auto",
        origin="lower",
        extent=(
            grid.axes_bounds[0][0],
            grid.axes_bounds[0][1],
            min(times),
            max(times),
        ),
        cmap=cmap,
        norm=norm,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    fig.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def create_spacetime_plot_adjacent(  # noqa: PLR0913, PLR0914
    snapshots: list[tuple[float, np.ndarray, np.ndarray]],
    grid: CartesianGrid,
    output_path: str | pathlib.Path,
    *,
    titles: tuple[str, str] = (
        r"$\phi_{0}(x,t)$",
        r"$\phi_{1}(x,t)$",
    ),
    xlabel: str = r"$x$",
    ylabel: str = r"$t$",
    cbar_label: str = r"$\phi$",
    cmap: str = "bwr",
    use_twoslope_norm: bool = True,
    dpi: int = 200,
    figsize: tuple[float, float] = (12, 6),
) -> None:
    """
    Create side-by-side spacetime heatmaps for coupled 1D fields.

    Parameters
    ----------
    snapshots : list[tuple[float, np.ndarray, np.ndarray]]
        List of (time, field0_data, field1_data) tuples.
    grid : CartesianGrid
        Grid providing spatial bounds.
    output_path : str | pathlib.Path
        Output file path.
    titles : tuple[str, str], optional
        Titles for left and right panels.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    cbar_label : str, optional
        Colorbar label.
    cmap : str, optional
        Colormap name.
    use_twoslope_norm : bool, optional
        If True, use TwoSlopeNorm centered at zero.
    dpi : int, optional
        Output DPI.
    figsize : tuple[float, float], optional
        Figure size.

    Raises
    ------
    ValueError
        If no snapshots provided.
    RuntimeError
        If no image was created.
    """
    if not snapshots:
        msg = "No snapshots provided"
        raise ValueError(msg)

    # Sort by time
    snapshots = sorted(snapshots, key=operator.itemgetter(0))

    # Build evolution arrays for both fields
    times = [t for t, *_ in snapshots]

    def _reshape(arr: np.ndarray) -> np.ndarray:
        return arr.reshape(grid.shape) if hasattr(grid, "shape") else arr.ravel()

    phi0_rows = [_reshape(arr0) for _, arr0, _ in snapshots]
    phi1_rows = [_reshape(arr1) for _, _, arr1 in snapshots]

    data = [
        np.vstack([row.ravel() for row in phi0_rows]),
        np.vstack([row.ravel() for row in phi1_rows]),
    ]

    # Setup colormap normalization across both fields
    vmin = min(data[0].min(), data[1].min())
    vmax = max(data[0].max(), data[1].max())

    if use_twoslope_norm:
        abs_max = max(abs(vmin), abs(vmax))
        norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Create plot with GridSpec for better colorbar placement
    pathlib.Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.03], wspace=0.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

    extent = (
        grid.axes_bounds[0][0],
        grid.axes_bounds[0][1],
        min(times),
        max(times),
    )

    im = None
    for i, (ax, title) in enumerate([(ax0, titles[0]), (ax1, titles[1])]):
        im = ax.imshow(
            data[i],
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=cmap,
            norm=norm,
        )
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        if i == 0:
            ax.set_ylabel(ylabel)
        else:
            ax.tick_params(labelleft=False, left=False)

    if im is None:
        msg = "No image was created"
        raise RuntimeError(msg)

    # Add colorbar
    cbar = fig.colorbar(im, cax=fig.add_subplot(gs[0, 2]), orientation="vertical")
    cbar.set_label(cbar_label)

    fig.savefig(str(output_path), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def create_1d_line_animation(  # noqa: PLR0913
    snapshots: list[tuple[float, np.ndarray]],
    grid: CartesianGrid,
    output_path: str | pathlib.Path,
    *,
    title_template: str = r"Klein-Gordon evolution: $\phi(x)$ at t = {t:.2f}",
    xlabel: str = r"$x$",
    ylabel: str = r"$\phi$",
    fps: int | None = None,
    dpi: int = 150,
    figsize: tuple[float, float] = (8, 5),
    ylim: tuple[float, float] | None = None,
    line_color: str = "C0",
    grid_alpha: float = 0.3,
) -> None:
    """
    Create an animation of 1D field evolution as a line plot.

    Parameters
    ----------
    snapshots : list[tuple[float, np.ndarray]]
        List of (time, field_data) tuples where field_data is 1D.
    grid : CartesianGrid
        Grid providing spatial coordinates.
    output_path : str | pathlib.Path
        Output file path.
    title_template : str, optional
        Title template with {t} placeholder for time.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    fps : int | None, optional
        Frames per second. If None, auto-computed.
    dpi : int, optional
        Output DPI.
    figsize : tuple[float, float], optional
        Figure size.
    ylim : tuple[float, float] | None, optional
        Y-axis limits. If None, auto-computed from data.
    line_color : str, optional
        Line color.
    grid_alpha : float, optional
        Grid alpha transparency.

    Raises
    ------
    ValueError
        If no snapshots provided.
    """
    if not snapshots:
        msg = "No snapshots provided"
        raise ValueError(msg)

    # Sort by time
    snapshots = sorted(snapshots, key=operator.itemgetter(0))

    # Get spatial coordinates
    x = cast("np.ndarray", grid.axes_coords[0])

    # Compute y-axis limits if not provided
    if ylim is None:
        all_min = min(frame.min() for _, frame in snapshots)
        all_max = max(frame.max() for _, frame in snapshots)
        y_padding = 0.1 * (all_max - all_min) if all_max > all_min else 0.1
        ylim = (all_min - y_padding, all_max + y_padding)

    # Setup figure
    fig, ax = plt.subplots(figsize=figsize)
    (line,) = ax.plot(x, snapshots[0][1].ravel(), color=line_color)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(ylim)
    ax.grid(visible=True, alpha=grid_alpha)
    title = ax.set_title(title_template.format(t=snapshots[0][0]))

    # Select writer and create animation
    pathlib.Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    # Extract t_end from snapshots
    t_end = snapshots[-1][0]
    writer, out, use_writer = choose_writer_and_out(
        len(snapshots), t_end, str(output_path), fps=fps
    )

    with writer.saving(fig, str(out), dpi=dpi):
        for t, frame in tqdm(snapshots, desc="Writing frames", total=len(snapshots)):
            line.set_ydata(frame.ravel())
            title.set_text(title_template.format(t=t))
            fig.canvas.draw()
            writer.grab_frame()

    plt.close(fig)
    logger.info("Saved %s using %s writer", out, use_writer)


def create_2d_heatmap_animation(  # noqa: PLR0913
    snapshots: list[tuple[float, np.ndarray]],
    grid: CartesianGrid,
    output_path: str | pathlib.Path,
    *,
    title_template: str = r"$\phi(x,y)$ at t = {t:.3f}",
    xlabel: str = r"$x$",
    ylabel: str = r"$y$",
    cbar_label: str = r"$\phi$",
    cmap: str = "bwr",
    fps: int | None = None,
    dpi: int = 150,
    figsize: tuple[float, float] = (6, 5),
    interpolation: str = "nearest",
    aspect: float | Literal["equal", "auto"] | None = "equal",
) -> None:
    """
    Create an animation of 2D field evolution as a heatmap.

    Parameters
    ----------
    snapshots : list[tuple[float, np.ndarray]]
        List of (time, field_data) tuples where field_data is 2D.
    grid : CartesianGrid
        Grid providing spatial bounds.
    output_path : str | pathlib.Path
        Output file path.
    title_template : str, optional
        Title template with {t} placeholder for time.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    cbar_label : str, optional
        Colorbar label.
    cmap : str, optional
        Colormap name.
    fps : int | None, optional
        Frames per second.
    dpi : int, optional
        Output DPI.
    figsize : tuple[float, float], optional
        Figure size.
    interpolation : str, optional
        Image interpolation method.
    aspect : str, optional
        Aspect ratio.

    Raises
    ------
    ValueError
        If no snapshots provided.
    """
    if not snapshots:
        msg = "No snapshots provided"
        raise ValueError(msg)

    # Sort by time
    snapshots = sorted(snapshots, key=operator.itemgetter(0))

    # Compute global color scale, centered at zero
    all_min = min(frame.min() for _, frame in snapshots)
    all_max = max(frame.max() for _, frame in snapshots)
    # Force symmetric range around zero for better visualization
    abs_max = max(abs(all_min), abs(all_max))
    all_min = -abs_max
    all_max = abs_max

    # Setup figure
    fig, ax = plt.subplots(figsize=figsize)
    first_frame = snapshots[0][1].reshape(grid.shape)
    # imshow expects (rows, columns) = (y, x) with origin='lower', but our data is (x, y)
    # So we need to transpose for correct visualization
    im = ax.imshow(
        first_frame.T,  # Transpose to convert from (x, y) to (y, x) for imshow
        origin="lower",
        extent=(
            grid.axes_bounds[0][0],
            grid.axes_bounds[0][1],
            grid.axes_bounds[1][0],
            grid.axes_bounds[1][1],
        ),
        cmap=cmap,
        vmin=all_min,
        vmax=all_max,
        interpolation=interpolation,
        aspect=aspect,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title_template.format(t=snapshots[0][0]))
    fig.colorbar(im, ax=ax, shrink=0.8, label=cbar_label)

    # Select writer and create animation
    pathlib.Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    t_end = snapshots[-1][0]
    writer, out, use_writer = choose_writer_and_out(
        len(snapshots), t_end, str(output_path), fps=fps
    )

    with writer.saving(fig, str(out), dpi=dpi):
        for t, frame in snapshots:
            frame2d = frame.reshape(grid.shape)
            im.set_data(frame2d.T)  # Transpose to match imshow's (y, x) expectation
            ax.set_title(title_template.format(t=t))
            fig.canvas.draw()
            writer.grab_frame()

    plt.close(fig)
    logger.info("Saved %s using %s writer", out, use_writer)


def create_2d_coupled_animation(  # noqa: PLR0913, PLR0914
    snapshots: list[tuple[float, np.ndarray, np.ndarray]],
    grid: CartesianGrid,
    output_path: str | pathlib.Path,
    *,
    titles: tuple[str, str] = ("Field 0: φ₀(x, y, t)", "Field 1: φ₁(x, y, t)"),
    xlabel: str = "x",
    ylabel: str = "y",
    cbar_labels: tuple[str, str] = ("φ₀", "φ₁"),
    cmap: str = "RdBu_r",
    use_twoslope_norm: bool = True,
    fps: int = 30,
    dpi: int = 150,
    figsize: tuple[float, float] = (12, 5),
    interpolation: str = "bilinear",
    aspect: str = "equal",
    time_text_y: float = 0.95,
) -> None:
    """
    Create a side-by-side animation of two coupled 2D fields.

    Parameters
    ----------
    snapshots : list[tuple[float, np.ndarray, np.ndarray]]
        List of (time, field0_data, field1_data) tuples.
    grid : CartesianGrid
        Grid providing spatial bounds.
    output_path : str | pathlib.Path
        Output file path.
    titles : tuple[str, str], optional
        Titles for left and right panels.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    cbar_labels : tuple[str, str], optional
        Colorbar labels for each field.
    cmap : str, optional
        Colormap name.
    use_twoslope_norm : bool, optional
        If True, use symmetric TwoSlopeNorm.
    fps : int, optional
        Frames per second.
    dpi : int, optional
        Output DPI.
    figsize : tuple[float, float], optional
        Figure size.
    interpolation : str, optional
        Image interpolation method.
    aspect : str, optional
        Aspect ratio.
    time_text_y : float, optional
        Y position for time display text (in figure coordinates).

    Raises
    ------
    ValueError
        If no snapshots provided.
    """
    if not snapshots:
        msg = "No snapshots provided"
        raise ValueError(msg)

    # Sort by time
    snapshots = sorted(snapshots, key=operator.itemgetter(0))

    # Determine global color range
    all_phi0 = np.array([phi0 for _, phi0, _ in snapshots])
    all_phi1 = np.array([phi1 for _, _, phi1 in snapshots])
    vmin = min(all_phi0.min(), all_phi1.min())
    vmax = max(all_phi0.max(), all_phi1.max())

    # Setup normalization
    if use_twoslope_norm:
        abs_max = max(abs(vmin), abs(vmax))
        norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    # Create figure
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=figsize)

    extent = [
        grid.axes_bounds[0][0],
        grid.axes_bounds[0][1],
        grid.axes_bounds[1][0],
        grid.axes_bounds[1][1],
    ]

    first_time, first_phi0, first_phi1 = snapshots[0]

    # Field 0 (left panel)
    im0 = ax0.imshow(
        first_phi0.T,
        origin="lower",
        extent=extent,
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        aspect=aspect,
    )
    ax0.set_xlabel(xlabel)
    ax0.set_ylabel(ylabel)
    ax0.set_title(titles[0])
    plt.colorbar(im0, ax=ax0, label=cbar_labels[0])

    # Field 1 (right panel)
    im1 = ax1.imshow(
        first_phi1.T,
        origin="lower",
        extent=extent,
        cmap=cmap,
        norm=norm,
        interpolation=interpolation,
        aspect=aspect,
    )
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    ax1.set_title(titles[1])
    plt.colorbar(im1, ax=ax1, label=cbar_labels[1])

    fig.tight_layout()

    # Add time display
    time_text = fig.text(
        0.5, time_text_y, f"t = {first_time:.2f}", ha="center", va="top", fontsize=12
    )

    # Create animation using FuncAnimation
    def update_frame(frame_idx: int) -> tuple[Any, Any, Any]:
        """Update function for animation."""
        t, phi0, phi1 = snapshots[frame_idx]
        im0.set_data(phi0.T)
        im1.set_data(phi1.T)
        time_text.set_text(f"t = {t:.2f}")
        return im0, im1, time_text

    anim = FuncAnimation(
        fig, update_frame, frames=len(snapshots), interval=1000 / fps, blit=True
    )

    # Save animation with writer selection
    pathlib.Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    if shutil.which("ffmpeg") is not None:
        writer = FFMpegWriter(fps=fps, bitrate=2000, codec="libx264")
        anim.save(str(output_path), writer=writer, dpi=dpi)
        logger.info("Saved %s using ffmpeg writer", output_path)
    else:
        # Fallback to GIF
        gif_path = pathlib.Path(output_path).with_suffix(".gif")
        writer = PillowWriter(fps=fps)
        anim.save(str(gif_path), writer=writer, dpi=dpi)
        logger.info("Saved %s using pillow writer", gif_path)

    plt.close(fig)
