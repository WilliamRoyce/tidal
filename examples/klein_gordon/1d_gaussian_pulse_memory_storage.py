"""Klein-Gordon simulation demonstrating MemoryStorage-based snapshot collection.

This example shows the improved snapshot collection pattern using py-pde's
built-in MemoryStorage tracker instead of manual list management. This is
the recommended approach for new code.

Compare with 1d_gaussian_pulse.py to see the difference.
"""

from __future__ import annotations

import logging
import pathlib

import matplotlib.pyplot as plt
import numpy as np

from tidal.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)

logger = logging.getLogger(__name__)


def main() -> None:
    """Run Klein-Gordon simulation using MemoryStorage for snapshot collection."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    # Grid setup
    grid_config = GridConfig(dim=1, shape=(512,), bounds=((0.0, 200.0),), periodic=True)
    grid = make_grid(grid_config)

    # Initial condition
    state = gaussian_pulse(
        grid, amplitude=1.0, width=5.0, center=[100.0], initial_velocity=0.0
    )

    # PDE and simulation config
    pde = KleinGordonPDE(params=KGParameters(mass=0.5))
    simulation_config = SimulationConfig(
        t_end=200.0,
        dt=None,  # Adaptive time step
        solver="scipy",
        method="RK45",
        backend="numpy",
        progress=True,
    )

    # Run simulation with automatic snapshot storage
    # This replaces the manual snapshot list pattern
    result, storage = run_with_snapshots(
        pde=pde,
        state=state,
        config=simulation_config,
        snapshot_interval=1.0,  # Snapshot every 1.0 time units
    )

    logger.info("Simulation complete")
    logger.info("  Collected %d snapshots", len(storage))
    logger.info("  Time range: [%.2f, %.2f]", storage.times[0], storage.times[-1])
    logger.info("  Final state: φ_max = %.4f", result[0].data.max())

    # Extract data from MemoryStorage
    # Access times and field data
    times = storage.times
    positions = grid.axes_coords[0]

    # Build spacetime array: phi_data[time_idx, space_idx]
    phi_data = np.array([storage[i][0].data for i in range(len(storage))])  # type: ignore[index]

    logger.info("  Spacetime array shape: %s", phi_data.shape)

    # Create spacetime heatmap
    fig, ax = plt.subplots(figsize=(10, 6))
    extent = (
        float(positions[0]),
        float(positions[-1]),
        float(times[0]),
        float(times[-1]),
    )
    im = ax.imshow(
        phi_data,
        aspect="auto",
        origin="lower",
        extent=extent,
        cmap="RdBu_r",
        interpolation="bilinear",
    )
    ax.set_xlabel(r"Position $x$", fontsize=12)
    ax.set_ylabel(r"Time $t$", fontsize=12)
    ax.set_title(
        r"Klein-Gordon $\phi(x,t)$ using MemoryStorage", fontsize=13, fontweight="bold"
    )
    fig.colorbar(im, ax=ax, label=r"$\phi$")

    # Save output
    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out_path = "outputs/kg_memory_storage_demo.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved %s", out_path)

    # Demonstrate accessing individual snapshots
    logger.info("--- MemoryStorage Access Examples ---")
    logger.info("First snapshot time: %.2f", storage.times[0])
    logger.info("First snapshot φ max: %.4f", storage[0][0].data.max())  # type: ignore[index]
    logger.info("Last snapshot time: %.2f", storage.times[-1])
    logger.info("Last snapshot φ max: %.4f", storage[-1][0].data.max())  # type: ignore[index]

    # Can also iterate over storage
    logger.info("Snapshot times (first 5): %s", list(storage.times[:5]))


if __name__ == "__main__":
    main()
