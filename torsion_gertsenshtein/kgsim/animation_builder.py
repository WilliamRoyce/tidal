"""Animation builder for simplified animation creation.

This module provides AnimationBuilder, a fluent API for creating animations
and plots from Klein-Gordon simulations with reduced boilerplate.
"""

from __future__ import annotations

import pathlib
import shutil
from dataclasses import dataclass
from typing import TYPE_CHECKING

import matplotlib as mpl
from pde import FieldCollection

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import lines, text
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.colors import Normalize, TwoSlopeNorm

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.image import AxesImage
    from pde import CartesianGrid, MemoryStorage


@dataclass
class AnimationConfig:
    """Configuration for animation creation."""

    output_path: str | pathlib.Path
    title: str = "Klein-Gordon Evolution"
    cmap: str = "RdBu_r"
    fps: int | None = None
    dpi: int = 200
    figsize: tuple[float, float] = (8, 6)
    xlabel: str = "x"
    ylabel: str = "t"
    cbar_label: str = "φ"
    use_twoslope_norm: bool = True


class AnimationBuilder:
    """Fluent API for creating Klein-Gordon animations and plots.

    This class consolidates common animation patterns and eliminates boilerplate
    by providing a unified interface for different visualization types.

    Examples
    --------
    >>> from torsion_gertsenshtein.kgsim import AnimationBuilder, AnimationConfig
    >>> # After running simulation with run_with_snapshots()
    >>> builder = AnimationBuilder(storage, grid)
    >>> config = AnimationConfig(output_path="output/animation.mp4", title="KG Evolution")
    >>> builder.create_2d_heatmap_animation(config)

    Or use method chaining:

    >>> AnimationBuilder(storage, grid).create_spacetime_1d(
    ...     AnimationConfig(output_path="output/spacetime.png")
    ... )
    """

    def __init__(self, storage: MemoryStorage, grid: CartesianGrid) -> None:
        """Initialize animation builder.

        Parameters
        ----------
        storage : MemoryStorage
            Storage containing simulation snapshots
        grid : CartesianGrid
            Grid used in simulation
        """
        self.storage = storage
        self.grid = grid
        self.times = storage.times

    @staticmethod
    def _choose_writer(
        snap_count: int, t_end: float, fps: int | None = None
    ) -> tuple[FFMpegWriter | PillowWriter, str]:
        """Choose animation writer (ffmpeg or pillow)."""
        if fps is None:
            fps = max(1, int(snap_count / max(1.0, t_end / 5.0)))

        if shutil.which("ffmpeg") is not None:
            writer = FFMpegWriter(fps=fps, metadata={"artist": "kgsim"}, bitrate=2000)
            ext = ".mp4"
        else:
            writer = PillowWriter(fps=fps)
            ext = ".gif"

        return writer, ext

    @staticmethod
    def _setup_colormap_norm(
        data: np.ndarray, *, use_twoslope: bool = True
    ) -> Normalize | TwoSlopeNorm:
        """Set up colormap normalization."""
        vmin, vmax = data.min(), data.max()

        if use_twoslope and vmin < 0 < vmax:
            return TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
        return Normalize(vmin=vmin, vmax=vmax)

    def create_spacetime_1d(self, config: AnimationConfig) -> None:
        """Create spacetime heatmap for 1D simulations.

        Parameters
        ----------
        config : AnimationConfig
            Animation configuration

        Raises
        ------
        ValueError
            If grid is not 1D
        """
        if self.grid.dim != 1:
            msg = f"spacetime_1d requires 1D grid, got {self.grid.dim}D"
            raise ValueError(msg)

        # Extract phi data from all snapshots
        phi_data = np.array([self.storage[i][0].data for i in range(len(self.storage))])  # type: ignore[index,misc]

        # Get spatial coordinates
        x = self.grid.axes_coords[0]

        # Create figure
        fig, ax = plt.subplots(figsize=config.figsize)

        # Setup extent and plot
        extent = (
            float(x[0]),
            float(x[-1]),
            float(self.times[0]),
            float(self.times[-1]),
        )
        norm = self._setup_colormap_norm(
            phi_data, use_twoslope=config.use_twoslope_norm
        )

        im = ax.imshow(
            phi_data,
            aspect="auto",
            origin="lower",
            extent=extent,
            cmap=config.cmap,
            norm=norm,
            interpolation="bilinear",
        )

        ax.set_xlabel(config.xlabel, fontsize=12)
        ax.set_ylabel(config.ylabel, fontsize=12)
        ax.set_title(config.title, fontsize=13, fontweight="bold")
        fig.colorbar(im, ax=ax, label=config.cbar_label)

        # Save
        output_path = pathlib.Path(config.output_path)
        output_path.parent.mkdir(exist_ok=True, parents=True)
        fig.savefig(output_path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)

    def create_2d_heatmap_animation(self, config: AnimationConfig) -> None:
        """Create animated 2D heatmap.

        Parameters
        ----------
        config : AnimationConfig
            Animation configuration

        Raises
        ------
        ValueError
            If grid is not 2D
        """
        if self.grid.dim != 2:  # noqa: PLR2004
            msg = f"2d_heatmap_animation requires 2D grid, got {self.grid.dim}D"
            raise ValueError(msg)

        # Get spatial extent
        bounds = self.grid.axes_bounds
        extent = (
            float(bounds[0][0]),
            float(bounds[0][1]),
            float(bounds[1][0]),
            float(bounds[1][1]),
        )

        # Setup figure
        fig, ax = plt.subplots(figsize=config.figsize)

        # Get first frame and setup colormap
        first_frame = self.storage[0][0].data  # type: ignore[index,misc]
        all_data = np.array([self.storage[i][0].data for i in range(len(self.storage))])  # type: ignore[index,misc]
        norm = self._setup_colormap_norm(
            all_data, use_twoslope=config.use_twoslope_norm
        )

        # Create initial image
        im = ax.imshow(
            first_frame.T,  # type: ignore[arg-type]  # Transpose for correct orientation
            origin="lower",
            extent=extent,
            cmap=config.cmap,
            norm=norm,
            animated=True,
        )

        ax.set_xlabel(config.xlabel, fontsize=12)
        ax.set_ylabel(config.ylabel if config.ylabel != "t" else "y", fontsize=12)
        ax.set_aspect("equal")

        title_text = ax.text(
            0.5,
            1.02,
            "",
            transform=ax.transAxes,
            ha="center",
            fontsize=13,
            fontweight="bold",
        )

        fig.colorbar(im, ax=ax, label=config.cbar_label)

        def update(frame_idx: int) -> tuple[AxesImage, text.Text]:
            """Update function for animation."""
            frame2d = self.storage[frame_idx][0].data  # type: ignore[index,misc]
            im.set_data(frame2d.T)  # type: ignore[arg-type]  # Transpose for correct orientation
            title_text.set_text(f"{config.title} (t={self.times[frame_idx]:.2f})")
            return im, title_text

        # Create animation
        anim = FuncAnimation(
            fig, update, frames=len(self.storage), interval=50, blit=True
        )

        # Choose writer and save
        writer, ext = self._choose_writer(len(self.storage), self.times[-1], config.fps)

        # Ensure output path has correct extension
        out_path = pathlib.Path(config.output_path)
        if out_path.suffix not in {".mp4", ".gif"}:
            out_path = out_path.with_suffix(ext)

        out_path.parent.mkdir(exist_ok=True, parents=True)
        anim.save(str(out_path), writer=writer, dpi=config.dpi)
        plt.close(fig)

    def create_1d_line_animation(self, config: AnimationConfig) -> None:
        """Create animated 1D line plot.

        Parameters
        ----------
        config : AnimationConfig
            Animation configuration

        Raises
        ------
        ValueError
            If grid is not 1D
        """
        if self.grid.dim != 1:
            msg = f"1d_line_animation requires 1D grid, got {self.grid.dim}D"
            raise ValueError(msg)

        x = self.grid.axes_coords[0]

        # Get data range for y-limits
        all_data = np.array([self.storage[i][0].data for i in range(len(self.storage))])  # type: ignore[index,misc]
        ymin, ymax = all_data.min(), all_data.max()
        y_margin = 0.1 * (ymax - ymin)

        # Setup figure
        fig, ax = plt.subplots(figsize=config.figsize)
        (line,) = ax.plot([], [], "b-", linewidth=2)
        ax.set_xlim(x[0], x[-1])
        ax.set_ylim(ymin - y_margin, ymax + y_margin)
        ax.set_xlabel(config.xlabel, fontsize=12)
        ax.set_ylabel(config.cbar_label, fontsize=12)
        ax.grid(visible=True, alpha=0.3)

        title_text = ax.set_title("", fontsize=13, fontweight="bold")

        def update(frame_idx: int) -> tuple[lines.Line2D, text.Text]:
            """Update function for animation."""
            y_data = self.storage[frame_idx][0].data  # type: ignore[index,misc]
            line.set_data(x, y_data)  # type: ignore[arg-type]
            title_text.set_text(f"{config.title} (t={self.times[frame_idx]:.2f})")
            return line, title_text

        # Create animation
        anim = FuncAnimation(
            fig, update, frames=len(self.storage), interval=50, blit=True
        )

        # Save
        writer, ext = self._choose_writer(len(self.storage), self.times[-1], config.fps)
        out_path = pathlib.Path(config.output_path)
        if out_path.suffix not in {".mp4", ".gif"}:
            out_path = out_path.with_suffix(ext)

        out_path.parent.mkdir(exist_ok=True, parents=True)
        anim.save(str(out_path), writer=writer, dpi=config.dpi)
        plt.close(fig)

    def create_dual_field_animation(
        self, config: AnimationConfig, field_labels: tuple[str, str] = ("φ₀", "φ₁")
    ) -> None:
        """Create side-by-side animation for two-field systems.

        Parameters
        ----------
        config : AnimationConfig
            Animation configuration
        field_labels : tuple[str, str]
            Labels for the two fields

        Raises
        ------
        ValueError
            If grid is not 2D or storage does not contain at least 2 fields
        """
        if self.grid.dim != 2:  # noqa: PLR2004
            msg = f"dual_field_animation requires 2D grid, got {self.grid.dim}D"
            raise ValueError(msg)

        # Check we have at least 2 fields

        first_snapshot = self.storage[0]
        if isinstance(first_snapshot, FieldCollection) and len(first_snapshot) < 2:  # noqa: PLR2004
            msg = "dual_field_animation requires at least 2 fields in storage"
            raise ValueError(msg)

        # Setup figure and extract data for both fields
        fig, axes = self._setup_dual_field_figure(config, field_labels)
        im0, im1 = self._create_dual_field_images(axes, config)

        def update(frame_idx: int) -> tuple[AxesImage, AxesImage]:
            """Update function for animation."""
            im0.set_data(self.storage[frame_idx][0].data.T)  # type: ignore[index,misc]
            im1.set_data(self.storage[frame_idx][1].data.T)  # type: ignore[index,misc]
            fig.suptitle(f"{config.title} (t={self.times[frame_idx]:.2f})")
            return im0, im1

        # Create and save animation
        self._save_dual_animation(fig, update, config)

    @staticmethod
    def _setup_dual_field_figure(
        config: AnimationConfig, field_labels: tuple[str, str]
    ) -> tuple[Figure, tuple[Axes, Axes]]:
        """Set up figure with two subplots for dual field animation."""
        fig, (ax0, ax1) = plt.subplots(
            1, 2, figsize=(config.figsize[0] * 1.8, config.figsize[1])
        )

        for ax, label in zip([ax0, ax1], field_labels, strict=False):
            ax.set_xlabel(config.xlabel, fontsize=11)
            ax.set_ylabel("y", fontsize=11)
            ax.set_title(label, fontsize=12)
            ax.set_aspect("equal")

        fig.suptitle("", fontsize=13, fontweight="bold")
        return fig, (ax0, ax1)

    def _create_dual_field_images(
        self, axes: tuple[Axes, Axes], config: AnimationConfig
    ) -> tuple[AxesImage, AxesImage]:
        """Create initial images for dual field animation."""
        ax0, ax1 = axes
        bounds = self.grid.axes_bounds
        extent = (
            float(bounds[0][0]),
            float(bounds[0][1]),
            float(bounds[1][0]),
            float(bounds[1][1]),
        )

        # Get data for normalization
        all_data0 = np.array(
            [self.storage[i][0].data for i in range(len(self.storage))]  # type: ignore[index,misc]
        )
        all_data1 = np.array(
            [self.storage[i][1].data for i in range(len(self.storage))]  # type: ignore[index,misc]
        )

        norm0 = self._setup_colormap_norm(
            all_data0, use_twoslope=config.use_twoslope_norm
        )
        norm1 = self._setup_colormap_norm(
            all_data1, use_twoslope=config.use_twoslope_norm
        )

        # Create initial images
        im0 = ax0.imshow(
            self.storage[0][0].data.T,  # type: ignore[index,misc]
            origin="lower",
            extent=extent,
            cmap=config.cmap,
            norm=norm0,
            animated=True,
        )
        im1 = ax1.imshow(
            self.storage[0][1].data.T,  # type: ignore[index,misc]
            origin="lower",
            extent=extent,
            cmap=config.cmap,
            norm=norm1,
            animated=True,
        )

        # Add colorbars
        ax0.figure.colorbar(im0, ax=ax0, label="", shrink=0.8)  # type: ignore[union-attr]
        ax1.figure.colorbar(im1, ax=ax1, label="", shrink=0.8)  # type: ignore[union-attr]

        return im0, im1

    def _save_dual_animation(
        self,
        fig: Figure,
        update_func: object,
        config: AnimationConfig,
    ) -> None:
        """Save dual field animation to file."""
        anim = FuncAnimation(
            fig,
            update_func,  # type: ignore[arg-type]
            frames=len(self.storage),
            interval=50,
            blit=True,  # type: ignore[arg-type]
        )

        writer, ext = self._choose_writer(len(self.storage), self.times[-1], config.fps)
        out_path = pathlib.Path(config.output_path)
        if out_path.suffix not in {".mp4", ".gif"}:
            out_path = out_path.with_suffix(ext)

        out_path.parent.mkdir(exist_ok=True, parents=True)
        anim.save(str(out_path), writer=writer, dpi=config.dpi)
        plt.close(fig)  # type: ignore[arg-type]
