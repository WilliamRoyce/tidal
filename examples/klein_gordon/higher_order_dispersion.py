"""Higher-order Klein-Gordon: 1D wave with fourth-order dispersion.

This example demonstrates the effect of fourth-order spatial derivatives on
wave propagation. Standard Klein-Gordon (alpha_4=0) is compared with fourth-order
dispersion (alpha_4>0) to show how higher-order terms modify the dispersion relation.

The fourth-order term acts as dispersive correction, causing wave packets to spread
differently than in the standard case. This is relevant for:
- Quantum corrections to classical field theories
- Beam/plate equations in mechanics
- Modified dispersion relations in condensed matter
"""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)
from torsion_gertsenshtein.kgsim.advanced_equations import HigherOrderKGPDE


def build_grid_and_state() -> dict[str, Any]:
    """Build 1D periodic grid with centered Gaussian pulse.

    Returns
    -------
    dict[str, Any]
        Dictionary containing grid and initial state.
    """
    grid_config = GridConfig(
        dim=1,
        shape=(512,),
        bounds=((0.0, 200.0),),
        periodic=True,
    )
    grid = make_grid(grid_config)
    state = gaussian_pulse(
        grid,
        amplitude=1.0,
        width=5.0,
        center=[100.0],
    )
    return {"grid": grid, "state": state}


def run_simulation(
    grid: Any,
    state: Any,
    alpha_4: float,
    label: str,
) -> Any:
    """Run simulation with specified fourth-order coefficient.

    Parameters
    ----------
    grid : Any
        Simulation grid.
    state : Any
        Initial state.
    alpha_4 : float
        Fourth-order dispersion coefficient.
    label : str
        Label for progress output.

    Returns
    -------
    Any
        Storage containing simulation snapshots.
    """
    pde = HigherOrderKGPDE(
        mass=0.5,
        alpha_2=1.0,  # Standard Laplacian coefficient
        alpha_4=alpha_4,  # Fourth-order dispersion
        alpha_6=0.0,  # No sixth-order terms
    )

    sim_config = SimulationConfig(
        t_end=200.0,
        dt=None,
        backend="numba",
        solver="scipy",
        method="RK45",
        progress=True,
    )

    print(f"\n{label}:")
    print(f"  alpha_2 = {pde.alpha_2:.3f}, alpha_4 = {pde.alpha_4:.6f}")

    _result, storage = run_with_snapshots(
        pde=pde,
        state=state,
        config=sim_config,
        snapshot_interval=2.0,
    )

    print(f"  Recorded {len(storage.times)} snapshots")
    return storage


def main() -> None:
    """Run standard and higher-order simulations and compare."""
    print("Higher-Order Klein-Gordon Dispersion Comparison")
    print("=" * 60)

    setup = build_grid_and_state()
    grid = setup["grid"]
    state_initial = setup["state"]

    print(f"Grid: {grid.shape} points, bounds {grid.axes_bounds}")
    print("Initial: Gaussian pulse, width=5.0, center=100.0")

    # Run standard Klein-Gordon (no fourth-order term)
    storage_standard = run_simulation(
        grid,
        state_initial.copy(),
        alpha_4=0.0,
        label="Standard Klein-Gordon (alpha_4=0)",
    )

    # Run with fourth-order dispersion
    storage_dispersive = run_simulation(
        grid,
        state_initial.copy(),
        alpha_4=10,  # Fourth-order coefficient
        label="Higher-Order Klein-Gordon (alpha_4=10)",
    )

    # Create side-by-side spacetime comparison
    print("\nCreating comparison plot...")

    _fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Standard KG
    times_standard = storage_standard.times
    x_coords = np.asarray(grid.axes_coords[0])
    data_standard = np.array(
        [storage_standard[i][0].data for i in range(len(storage_standard))]
    )

    # Get global vmax for both plots
    data_dispersive = np.array(
        [storage_dispersive[i][0].data for i in range(len(storage_dispersive))]
    )
    vmax = max(
        np.abs(data_standard.min()),
        np.abs(data_standard.max()),
        np.abs(data_dispersive.min()),
        np.abs(data_dispersive.max()),
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    im1 = ax1.imshow(
        data_standard,
        aspect="auto",
        extent=[x_coords[0], x_coords[-1], times_standard[0], times_standard[-1]],
        origin="lower",
        cmap="bwr",
        norm=norm,
    )
    ax1.set_xlabel(r"$x$")
    ax1.set_ylabel(r"$t$")
    ax1.set_title(r"Standard KG ($\alpha_4=0$)")
    plt.colorbar(im1, ax=ax1, label=r"$\phi$")

    # Higher-order KG
    times_dispersive = storage_dispersive.times

    im2 = ax2.imshow(
        data_dispersive,
        aspect="auto",
        extent=[x_coords[0], x_coords[-1], times_dispersive[0], times_dispersive[-1]],
        origin="lower",
        cmap="bwr",
        norm=norm,
    )
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"$t$")
    ax2.set_title(r"Higher-Order KG ($\alpha_4=10$)")
    plt.colorbar(im2, ax=ax2, label=r"$\phi$")

    plt.tight_layout()
    output_path = "outputs/higher_order_comparison.pdf"
    plt.savefig(output_path, dpi=200)
    plt.close()

    print(f"Saved comparison plot: {output_path}")
    print("\nExpected differences:")
    print("  - Standard: Clean propagation with minimal spreading")
    print("  - Higher-order: Additional dispersion from fourth-order term")
    print("  - Fourth-order term modifies high-k modes more strongly")


if __name__ == "__main__":
    main()
