from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING, Any

import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm

mpl.use("Agg")
import operator

import matplotlib.pyplot as plt
import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
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

    params = KGParameters(mass=0.5)
    pde = KleinGordonPDE(params)

    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # adaptive
        solver="scipy",  # or "explicit"
        method="RK45",
        backend="numpy",  # prefer 'numpy' here for portability unless numba RHS is provided
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

    Raises
    ------
    RuntimeError
        If no snapshots were recorded during the simulation.
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

    # Build evolution array: shape (nt, nx) with time as vertical axis
    times = [t for t, _ in snapshots]
    phi_rows = [
        arr.reshape(simulation_components["grid"].shape)
        if hasattr(simulation_components["grid"], "shape")
        else arr.ravel()
        for _, arr in snapshots
    ]
    data = np.vstack([row.ravel() for row in phi_rows])  # (nt, nx)

    # Center colormap on zero even if data is asymmetric
    vmin = float(data.min())
    vmax = float(data.max())
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out = "outputs/KG_evolution.png"

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        data,
        aspect="auto",
        origin="lower",
        extent=(
            simulation_components["grid"].axes_bounds[0][0],
            simulation_components["grid"].axes_bounds[0][1],
            min(times),
            max(times),
        ),
        cmap="bwr",
        norm=norm,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(r"Klein-Gordon evolution: $\phi(x,t)$")
    fig.colorbar(im, ax=ax, label=r"$\phi$")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
