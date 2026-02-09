from __future__ import annotations

import logging
from typing import Any

from tidal.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

from tidal.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)


def _build_simulation_components() -> dict[str, Any]:
    """
    Build and return all components needed to run the simulation as a single dict.

    Returning a single container lets `main` keep few local variables.
    """
    # --- grid ---
    grid_config = GridConfig(
        dim=1, shape=(1024,), bounds=((0.0, 200.0),), periodic=True
    )
    grid = make_grid(grid_config)

    pde = KleinGordonPDE(params=KGParameters(mass=0.5))

    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # Adaptive time step
        solver="scipy",  # or "explicit"
        method="RK45",
        backend="numba",  # prefer 'numpy' here for portability unless numba RHS is provided
        progress=True,
    )

    # --- any observers / snapshots ---
    observers: list[Any] = []  # keep as list to be appended inside builder if needed

    return {
        "grid_config": grid_config,
        "grid": grid,
        "state": state,
        "pde": pde,
        "simulation_config": simulation_config,
        "observers": observers,
    }


def main() -> None:
    """
    Run a 1D Klein-Gordon simulation of a Gaussian pulse.

    This function performs the following steps:
    - Construct a 1D periodic computational grid (default: 1024 points over [0, 200]).
    - Instantiate Klein-Gordon PDE parameters (mass=0.5) and the PDE object.
    - Initialize the field as a Gaussian pulse with specified amplitude and width.
    - Collect snapshots of φ at the initial time and at each tracker interrupt via
        an observer callback. Snapshots are stored as a list of (time, numpy.ndarray)
        tuples where each array has shape (nx,).
    - Configure and run the time integrator (adaptive dt with the "RK45" SciPy
        solver by default, using the "numpy" backend). The run call fills the
        snapshots list as a side effect.
    - Assemble the recorded snapshots into a 2D array of shape (nt, nx) where
        rows correspond to times and columns to spatial grid points.
    - Create and save an image of the evolution using matplotlib.imshow with the
        horizontal axis as x and the vertical axis as t. The output image is written
        to "outputs/phi_evolution.png" (the outputs directory is created if needed).
    - Print the saved file path to stdout.
    - Assemble the recorded snapshots and create an animation showing φ(x) vs x
        for each time step.
    - Save the animation as either MP4 (if ffmpeg available) or GIF.

    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    simulation_components = _build_simulation_components()

    # Run simulation with automatic snapshot storage
    _result, storage = run_with_snapshots(
        pde=simulation_components["pde"],
        state=simulation_components["state"],
        config=simulation_components["simulation_config"],
        snapshot_interval=0.5,
    )

    # Create 1D line animation using AnimationBuilder
    builder = AnimationBuilder(storage, simulation_components["grid"])
    config = AnimationConfig(
        output_path="outputs/KG_evolution.mp4",
        title=r"Klein-Gordon evolution: $\phi(x)$",
        xlabel=r"$x$",
        cbar_label=r"$\phi$",
        fps=30,
    )
    builder.create_1d_line_animation(config)
    print("Saved outputs/KG_evolution.[mp4|gif]")


if __name__ == "__main__":
    main()
