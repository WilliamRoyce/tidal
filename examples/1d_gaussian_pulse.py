from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING, Any

import matplotlib as mpl

mpl.use("Agg")
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


def main() -> None:
    grid_config = GridConfig(
        dim=1, shape=(1024,), bounds=((0.0, 200.0),), periodic=True
    )
    grid = make_grid(grid_config)

    params = KGParameters(mass=0.5)
    pde = KleinGordonPDE(params)

    state = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    # Collector for snapshots (time, phi_array)
    snapshots: list[tuple[float, np.ndarray]] = []

    # Record initial condition
    phi0 = state[0]
    snapshots.append((0.0, np.asarray(phi0.data).copy()))

    # Observer that records φ at each tracker interrupt
    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        phi_field = state_coll[0]
        snapshots.append((float(t), np.asarray(phi_field.data).copy()))
        return {}

    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # adaptive
        solver="scipy",  # or "explicit"
        method="RK45",
        backend="numba",
        progress=True,
    )

    # run accepts an extra_observer callback (side-effect: fills snapshots)
    run(pde=pde, state=state, config=simulation_config, extra_observer=record_phi)

    # Build evolution array: shape (nt, nx) with time as vertical axis
    times = [t for t, _ in snapshots]
    phi_rows = [arr.reshape(-1) for _, arr in snapshots]
    data = np.vstack(phi_rows)  # (nt, nx)

    x0, x1 = grid.axes_bounds[0]
    t0, t1 = min(times), max(times)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        data,
        aspect="auto",
        origin="lower",
        extent=(x0, x1, t0, t1),
    )
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_title(r"Klein-Gordon evolution: $\phi(x,t)$")
    fig.colorbar(im, ax=ax, label=r"$\phi$")

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out = "outputs/phi_evolution.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
