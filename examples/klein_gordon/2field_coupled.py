from __future__ import annotations

import logging
from typing import Any

from tidal.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

from tidal.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GridConfig,
    MultiFieldParams,
    SimulationConfig,
    make_coupled_kg_pde,
    make_grid,
    multi_gaussian,
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

    # Two fields with masses and symmetric coupling
    masses = [0.25, 1.0]
    g = 0.2  # off-diagonal bilinear coupling
    coupling = [[0.0, g], [g, 0.0]]
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # IC: excite only field 0
    state = multi_gaussian(
        grid, amplitudes=[1.0, 0.0], widths=[5.0, 5.0], velocities=[0.0, 0.0]
    )

    # Sim config
    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,
        solver="scipy",
        method="RK45",
        backend="numba",
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
    """Run the coupled Klein-Gordon simulation.

    This function collects snapshots of both fields and saves a figure showing the
    space-time evolution.

    This function builds the simulation components, records initial and tracked
    snapshots during the run, assembles the evolution arrays for both fields,
    and writes an image file with two subplots (one per field). It does not
    return a value.
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

    # Create dual field spacetime plot using AnimationBuilder
    builder = AnimationBuilder(storage, simulation_components["grid"])
    config = AnimationConfig(
        output_path="outputs/coupled_fields.pdf",
        title="Coupled Klein-Gordon Fields",
        xlabel=r"$x$",
        ylabel=r"$t$",
        cbar_label=r"$\phi$",
    )
    builder.create_spacetime_1d_dual(config, field_labels=(r"$\phi_0$", r"$\phi_1$"))
    print("Saved outputs/coupled_fields.pdf")


if __name__ == "__main__":
    main()
