from __future__ import annotations

from typing import TYPE_CHECKING, Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

import operator

import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
    ring_pulse_2d,
    run,
)

# Create 2D heatmap animation using unified plotting module
from torsion_gertsenshtein.kgsim.animations import create_2d_heatmap_animation

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import CartesianGrid, FieldCollection


def build_grid_and_state() -> tuple[CartesianGrid, FieldCollection]:
    """
    Build and return a 2D Cartesian grid together with an initial field state containing a ring-shaped pulse.

    Returns
    -------
        tuple[CartesianGrid, FieldCollection]:
            A tuple where the first element is a 2D CartesianGrid configured as a periodic domain
            with shape (256, 256) and bounds ((-50.0, 50.0), (-50.0, 50.0)), and the second element
            is a FieldCollection initialized by ring_pulse_2d on that grid. The pulse in the returned
            state has amplitude=1.0, initial_radius=15.0, and sigma=2.0.
    """
    grid_config = GridConfig(
        dim=2,
        shape=(256, 256),
        bounds=((-50.0, 50.0), (-50.0, 50.0)),
        periodic=True,
    )
    grid = make_grid(grid_config)
    state = ring_pulse_2d(grid, amplitude=1.0, initial_radius=15.0, sigma=2.0)
    return grid, state


def make_recorder() -> tuple[
    Callable[[FieldCollection, float], dict[str, Any]], list[tuple[float, np.ndarray]]
]:
    """
    Create a simple recorder callback and its backing storage for snapshots.

    The returned recorder appends a time-stamped copy of the first field's data
    (from state_coll[0].data) to the snapshots list each time it is called.

    Returns
    -------
        tuple[Callable[[FieldCollection, float], dict[str, Any]], list[tuple[float, numpy.ndarray]]]:
            A pair (recorder, snapshots) where:
            - recorder(state_coll, t) -> dict: a callable that records a snapshot.
              It appends (float(t), np.asarray(state_coll[0].data).copy()) to `snapshots`
              and returns an empty dict (useful for callback APIs expecting a dict).
            - snapshots: a list of tuples (time: float, data: numpy.ndarray) containing
              copies of the recorded field data in the order they were recorded.
    """
    snapshots: list[tuple[float, np.ndarray]] = []

    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        snapshots.append((float(t), np.asarray(state_coll[0].data).copy()))
        return {}

    return record_phi, snapshots


def main() -> None:
    """
    Run a Klein-Gordon PDE simulation, record snapshots, and save an animation.

    This top-level routine orchestrates a time-dependent simulation of the Klein-Gordon
    equation on a 2D computational grid and produces a movie of the primary field φ(x,y)
    over time. The function performs the following high-level steps:
    - Build the computational grid and initial state by calling the project-specific
        builder helper.
    - Create a recorder/observer that collects field snapshots during the simulation.
    - Configure and execute the PDE solver (adaptive time stepping is used by default).
    - Validate that snapshots were recorded; if none were captured, raise RuntimeError.
    - Sort collected snapshots by simulation time.
    - Prepare an output directory ("outputs") and a matplotlib figure for rendering frames.
    - Choose an appropriate matplotlib animation writer (e.g., FFMpegWriter when available)
        and normalize the color scale across all frames for visual consistency.
    - Iterate through snapshots, update the figure with each 2D frame, and have the writer
        record the frames into a single output file.
    - Close the figure and print a message reporting the final saved filename and writer.

    Notes
    -----
    - Side effects: creates/uses an "outputs" directory and writes an animation file there.
        It also prints a confirmation message identifying the saved file and writer used.
    - The function expects several project-specific helpers and types to be available
        in the module scope (e.g., build_grid_and_state, make_recorder, SimulationConfig,
        KleinGordonPDE, KGParameters, run, prepare_figure, choose_writer_and_out).
    - A suitable animation writer backend (such as ffmpeg) may be required on the host
        system for some writers to function.

    Raises
    ------
    RuntimeError
        If no snapshots were recorded during the simulation.
    """
    # Setup grid / PDE / initial state
    grid, state = build_grid_and_state()
    record_phi, snapshots = make_recorder()

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
        pde=KleinGordonPDE(KGParameters(mass=1.0)),
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

    create_2d_heatmap_animation(
        snapshots,
        grid,
        "outputs/phi_evolution_2d.mp4",
        title_template=r"$\phi(x,y)$ at t = {t:.3f}",
        xlabel=r"$x$",
        ylabel=r"$y$",
        cbar_label=r"$\phi$",
        cmap="bwr",
        interpolation="nearest",
    )


if __name__ == "__main__":
    main()
