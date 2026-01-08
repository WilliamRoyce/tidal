from __future__ import annotations

import operator
from typing import TYPE_CHECKING, Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    MultiFieldParams,
    SimulationConfig,
    make_coupled_kg_pde,
    make_grid,
    multi_gaussian,
    run,
)

# Create side-by-side spacetime heatmaps using unified plotting module
from torsion_gertsenshtein.kgsim.animations import create_spacetime_plot_adjacent

if TYPE_CHECKING:
    from pde import FieldCollection


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

    Raises
    ------
    RuntimeError
        If no snapshots are recorded during the simulation or if a non-finite
        field value is encountered during observation.
    """
    simulation_components = _build_simulation_components()

    # Collector for snapshots (time, phi_array)
    snapshots: list[tuple[float, np.ndarray, np.ndarray]] = []

    # Record initial condition
    snapshots.append(
        (
            0.0,
            np.asarray(simulation_components["state"][0].data).copy(),
            np.asarray(simulation_components["state"][2].data).copy(),
        )
    )

    # Observer that records φ at each tracker interrupt
    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        # fields order: phi0, pi0, phi1, pi1, ...
        arr0 = np.asarray(state_coll[0].data)
        arr1 = np.asarray(state_coll[2].data)
        if not (np.isfinite(arr0).all() and np.isfinite(arr1).all()):
            msg = f"Non-finite field at t={t}"
            raise RuntimeError(msg)
        # record both field profiles with the same timestamp
        snapshots.append((t, arr0.copy(), arr1.copy()))
        return {}

    run(
        pde=simulation_components["pde"],
        state=simulation_components["state"],
        config=simulation_components["simulation_config"],
        extra_observer=record_phi,
    )

    # Ensure there is at least one snapshot
    if not snapshots:
        msg = "No snapshots were recorded during the simulation."
        raise RuntimeError(msg)

    # Sort by time (observer callbacks might not be strictly increasing)
    snapshots.sort(key=operator.itemgetter(0))

    create_spacetime_plot_adjacent(
        snapshots,
        simulation_components["grid"],
        "outputs/KG_coupled_evolution.pdf",
        titles=(
            r"Coupled Klein-Gordon evolution: $\phi_{0}(x,t)$",
            r"Coupled Klein-Gordon evolution: $\phi_{1}(x,t)$",
        ),
        xlabel=r"$x$",
        ylabel=r"$t$",
        cbar_label=r"$\phi$",
    )
    print("Saved outputs/KG_coupled_evolution.pdf")


if __name__ == "__main__":
    main()
