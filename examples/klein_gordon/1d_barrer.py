from __future__ import annotations

from typing import TYPE_CHECKING, Any

# Create spacetime heatmap using unified plotting module
from torsion_gertsenshtein.kgsim.animations import create_spacetime_plot
from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

import operator

import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    InhomogeneousKGPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
    step_region_1d,
)

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

    # --- coefficients ---
    m_out = 0.25
    m_in = 5.0
    x0, x1 = 150.0, 150.5
    # Build m^2(x)
    m2_field = step_region_1d(
        grid,
        x0=x0,
        x1=x1,
        inside_value=m_in**2,
        outside_value=m_out**2,
    )

    # --- initial state ---
    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    # --- PDE / solver config ---
    pde = InhomogeneousKGPDE(m2_field=m2_field)
    # --- simulation config ---
    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # adaptive
        solver="scipy",
        method="RK45",
        backend="numpy",  # will auto-fallback to numpy for InhomogeneousKGPDE
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
    """Run a 1D Klein-Gordon simulation with a mass step, collect snapshots, and save an x-t image.

    This function delegates grid and field construction, PDE/state setup, time evolution
    and snapshot collection, and plotting to smaller helpers so the public function is
    documented and keeps a small number of local variables.

    Raises
    ------
    RuntimeError
        If no snapshots are recorded during the simulation or if a non-finite field is encountered during evolution.
    """
    simulation_components = _build_simulation_components()

    # Collector for snapshots (time, phi_array)
    snapshots: list[tuple[float, np.ndarray]] = []

    # Record initial condition
    snapshots.append((0.0, np.asarray(simulation_components["state"][0].data).copy()))

    # Observer that records φ at each tracker interrupt
    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        arr = np.asarray(state_coll[0].data)
        if not np.isfinite(arr).all():
            msg = f"Non-finite field at t={t}"
            raise RuntimeError(msg)
        snapshots.append((t, arr.copy()))
        return {}

    # run accepts an extra_observer callback (side-effect: fills snapshots)
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

    create_spacetime_plot(
        snapshots,
        simulation_components["grid"],
        "outputs/KG_barrier.pdf",
        title=r"Klein-Gordon $\phi(x,t)$ with Potential Barrier",
        xlabel=r"$x$",
        ylabel=r"$t$",
        cbar_label=r"$\phi$",
        figsize=(8, 6),
    )
    print("Saved outputs/KG_barrier.pdf")


if __name__ == "__main__":
    main()
