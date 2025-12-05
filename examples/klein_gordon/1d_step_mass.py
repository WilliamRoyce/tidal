from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING, Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

import operator

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

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
    m_out = 0.5
    m_in = 0.55
    x0, x1 = 125.0, 200.0
    # Build m^2(x)
    m2_field = step_region_1d(
        grid,
        x0=x0,
        x1=x1,
        inside_value=m_in**2,
        outside_value=m_out**2,
    )

    # --- initial state ---
    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=50.0)

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
    out = "outputs/KG_mass_step.pdf"

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
    ax.set_title(r"Klein-Gordon $\phi(x,t)$ with mass step")
    fig.colorbar(im, ax=ax, label=r"$\phi$")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
