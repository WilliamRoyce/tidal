from __future__ import annotations

import logging
from typing import Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

from torsion_gertsenshtein.kgsim import (
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
    - Run simulation with automatic snapshot collection using MemoryStorage.
    - Create spacetime heatmap using AnimationBuilder.
    - Save the output to "outputs/KG_evolution.pdf".
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    simulation_components = _build_simulation_components()

    # Run simulation with automatic snapshot storage
    _result, storage = run_with_snapshots(
        pde=simulation_components["pde"],
        state=simulation_components["state"],
        config=simulation_components["simulation_config"],
        snapshot_interval=1.0,
    )

    # Create spacetime plot using AnimationBuilder
    builder = AnimationBuilder(storage, simulation_components["grid"])
    config = AnimationConfig(
        output_path="outputs/KG_evolution.pdf",
        title=r"Klein-Gordon evolution: $\phi(x,t)$",
        xlabel=r"$x$",
        ylabel=r"$t$",
        cbar_label=r"$\phi$",
    )
    builder.create_spacetime_1d(config)
    print("Saved outputs/KG_evolution.pdf")


if __name__ == "__main__":
    main()
