from __future__ import annotations

import operator
import pathlib
from typing import TYPE_CHECKING, Any

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    MultiFieldParams,
    SimulationConfig,
    make_coupled_kg_pde,
    make_grid,
    multi_gaussian,
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


def main() -> None:  # noqa: PLR0914
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

    # Build evolution arrays for both fields: shape (nt, nx) with time as vertical axis
    # snapshots are tuples (t, phi0_array, phi1_array)
    times = [t for t, *_ in snapshots]

    def _reshape(arr: np.ndarray) -> np.ndarray:
        return (
            arr.reshape(simulation_components["grid"].shape)
            if hasattr(simulation_components["grid"], "shape")
            else arr.ravel()
        )

    phi0_rows = [_reshape(arr0) for _, arr0, _ in snapshots]
    phi1_rows = [_reshape(arr1) for _, _, arr1 in snapshots]

    data = [
        np.vstack([row.ravel() for row in phi0_rows]),
        np.vstack([row.ravel() for row in phi1_rows]),
    ]  # (nt, nx) for field 1

    # Center colormap on zero even if data is asymmetric
    vmin = np.min([data[0].min(), data[1].min()])
    vmax = np.max([data[0].max(), data[1].max()])
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)

    # Two panels + dedicated colorbar column via GridSpec to avoid overlap with tight_layout
    fig = plt.figure(figsize=(12, 6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.03], wspace=0.05)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1], sharey=ax0)

    im = None
    for i, a in enumerate([ax0, ax1]):
        im = a.imshow(
            data[i],
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
        a.set_title(r"Coupled Klein-Gordon evolution: " + rf"$\phi_{{{i}}}(x,t)$")
        a.set_xlabel("x")
        if i == 0:
            a.set_ylabel("t")
        else:
            a.tick_params(labelleft=False, left=False)

    if im is None:
        msg = "No image was created to attach a colorbar to."
        raise RuntimeError(msg)

    # place colorbar into the small third column created above
    cbar = fig.colorbar(im, cax=fig.add_subplot(gs[0, 2]), orientation="vertical")
    cbar.set_label(r"$\phi$")

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out = "outputs/KG_coupled_evolution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
