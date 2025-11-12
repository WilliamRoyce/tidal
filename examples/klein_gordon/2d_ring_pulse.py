from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING, Any

import matplotlib as mpl

mpl.use("Agg")
import operator
import shutil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from matplotlib.animation import PillowWriter

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
    ring_pulse_2d,
    run,
)

if TYPE_CHECKING:
    from pde import FieldCollection


def main() -> None:
    # Setup grid / PDE / initial state
    grid_config = GridConfig(
        dim=2,
        shape=(256, 256),
        bounds=((-50.0, 50.0), (-50.0, 50.0)),
        periodic=True,
    )
    grid = make_grid(grid_config)

    params = KGParameters(mass=1.0)
    pde = KleinGordonPDE(params)

    state = ring_pulse_2d(grid, amplitude=1.0, initial_radius=15.0, sigma=2.0)

    # collect snapshots (time, 2D array)
    snapshots: list[tuple[float, np.ndarray]] = []
    snapshots.append((0.0, np.asarray(state[0].data).copy()))

    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        phi_field = state_coll[0]
        snapshots.append((float(t), np.asarray(phi_field.data).copy()))
        return {}

    simulation_config = SimulationConfig(
        t_end=50.0,
        dt=None,  # adaptive
        solver="explicit",
        method="RK45",
        backend="numpy",
        progress=True,
    )

    # run simulation (fills snapshots via observer)
    run(
        pde=pde,
        state=state,
        config=simulation_config,
        extra_observer=record_phi,
        snapshot_interval=0.5,
    )

    if not snapshots:
        msg = "No snapshots recorded"
        raise RuntimeError(msg)

    # Sort by time
    snapshots.sort(key=operator.itemgetter(0))

    # prepare output dir
    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out = "outputs/phi_evolution_2d.mp4"

    # Prepare figure and writer
    first_frame = snapshots[0][1].reshape(grid.shape)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        first_frame,
        origin="lower",
        extent=(
            grid.axes_bounds[0][0],
            grid.axes_bounds[0][1],
            grid.axes_bounds[1][0],
            grid.axes_bounds[1][1],
        ),
        cmap="bwr",
        vmin=np.min(first_frame),
        vmax=np.max(first_frame),
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im, ax=ax, shrink=0.8, label=r"$\phi$")

    # choose writer (FFMpegWriter preferred)
    fps = max(
        1, int(len(snapshots) / max(1.0, simulation_config.t_end / 5.0))
    )  # simple heuristic

    # prefer ffmpeg if available; animation.writers.is_available checks writer backend
    if animation.writers.is_available("ffmpeg") and shutil.which("ffmpeg") is not None:
        writer = animation.FFMpegWriter
        writer = writer(fps=fps, metadata={"artist": "kgsim"}, bitrate=2000)
        use_writer = "ffmpeg"
    else:
        # fallback to PillowWriter (GIF) if ffmpeg not available
        writer = PillowWriter(fps=fps)
        out = "outputs/phi_evolution_2d.gif"
        use_writer = "pillow"

    # normalize color scale across all frames for visual consistency
    all_min = min(frame.min() for _, frame in snapshots)
    all_max = max(frame.max() for _, frame in snapshots)
    im.set_clim(all_min, all_max)

    with writer.saving(fig, str(out), dpi=150):
        for t, frame in snapshots:
            frame2d = frame.reshape(grid.shape)
            im.set_data(frame2d)
            ax.set_title(rf"$\phi(x,y)$ at t = {t:.3f}")
            # draw + grab
            fig.canvas.draw()
            writer.grab_frame()

    plt.close(fig)
    print(f"Saved {out} using {use_writer} writer (fps={fps})")


if __name__ == "__main__":
    main()
