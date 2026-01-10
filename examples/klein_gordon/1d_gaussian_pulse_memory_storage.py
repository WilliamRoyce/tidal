"""Klein-Gordon simulation demonstrating MemoryStorage-based snapshot collection.

This example shows the improved snapshot collection pattern using py-pde's
built-in MemoryStorage tracker instead of manual list management. This is
the recommended approach for new code.

Compare with 1d_gaussian_pulse.py to see the difference.
"""

from __future__ import annotations

import pathlib

import matplotlib.pyplot as plt
import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)


def main() -> None:
    """Run Klein-Gordon simulation using MemoryStorage for snapshot collection."""
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

    print("\n✓ Simulation complete")
    print(f"  Collected {len(storage)} snapshots")
    print(f"  Time range: [{storage.times[0]:.2f}, {storage.times[-1]:.2f}]")
    print(f"  Final state: φ_max = {result[0].data.max():.4f}")

    # Extract data from MemoryStorage
    # Access times and field data
    times = storage.times
    positions = grid.axes_coords[0]

    # Build spacetime array: phi_data[time_idx, space_idx]
    phi_data = np.array([storage[i][0].data for i in range(len(storage))])  # type: ignore[index]

    print(f"  Spacetime array shape: {phi_data.shape}")

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
    print(f"\n✓ Saved {out_path}")

    # Demonstrate accessing individual snapshots
    print("\n--- MemoryStorage Access Examples ---")
    print(f"First snapshot time: {storage.times[0]:.2f}")
    print(f"First snapshot φ max: {storage[0][0].data.max():.4f}")  # type: ignore[index]
    print(f"Last snapshot time: {storage.times[-1]:.2f}")
    print(f"Last snapshot φ max: {storage[-1][0].data.max():.4f}")  # type: ignore[index]

    # Can also iterate over storage
    print(f"\nSnapshot times (first 5): {list(storage.times[:5])}")


if __name__ == "__main__":
    main()
