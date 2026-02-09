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

import logging
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from tidal.plot_pgf import enable_pgf

enable_pgf("xelatex")

from tidal.kgsim import (
    GridConfig,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)
from tidal.kgsim.advanced_equations import HigherOrderKGPDE

if TYPE_CHECKING:
    from pde import CartesianGrid, FieldCollection, MemoryStorage

logger = logging.getLogger(__name__)


def build_grid_and_state() -> tuple[CartesianGrid, FieldCollection]:
    """Build 1D periodic grid with centered Gaussian pulse.

    Returns
    -------
    tuple[CartesianGrid, FieldCollection]
        Grid and initial state.
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
    return grid, state


def run_simulation(
    state: FieldCollection,
    alpha_4: float,
    label: str,
) -> MemoryStorage:
    """Run simulation with specified fourth-order coefficient.

    Parameters
    ----------
    state : FieldCollection
        Initial state.
    alpha_4 : float
        Fourth-order dispersion coefficient.
    label : str
        Label for progress output.

    Returns
    -------
    MemoryStorage
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

    logger.info("%s:", label)
    logger.info("  alpha_2 = %.3f, alpha_4 = %.6f", pde.alpha_2, pde.alpha_4)

    _result, storage = run_with_snapshots(
        pde=pde,
        state=state,
        config=sim_config,
        snapshot_interval=2.0,
    )

    logger.info("  Recorded %d snapshots", len(storage.times))
    return storage


def create_comparison_plot(
    grid: CartesianGrid,
    storage_standard: MemoryStorage,
    storage_dispersive: MemoryStorage,
    output_path: str = "outputs/higher_order_comparison.pdf",
) -> None:
    """Create side-by-side spacetime comparison plot.

    Parameters
    ----------
    grid : CartesianGrid
        Simulation grid.
    storage_standard : MemoryStorage
        Storage from standard KG simulation.
    storage_dispersive : MemoryStorage
        Storage from higher-order KG simulation.
    output_path : str
        Path to save the plot.
    """
    logger.info("Creating comparison plot...")

    _fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Extract data
    x_coords = np.asarray(grid.axes_coords[0])
    data_standard = np.array(
        [storage_standard[i][0].data for i in range(len(storage_standard))]  # type: ignore[index]
    )
    data_dispersive = np.array(
        [storage_dispersive[i][0].data for i in range(len(storage_dispersive))]  # type: ignore[index]
    )

    # Setup shared colormap
    vmax = max(
        np.abs(data_standard.min()),
        np.abs(data_standard.max()),
        np.abs(data_dispersive.min()),
        np.abs(data_dispersive.max()),
    )
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Plot standard KG
    times_standard = storage_standard.times
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

    # Plot higher-order KG
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
    plt.savefig(output_path, dpi=200)
    plt.close()

    logger.info("Saved comparison plot: %s", output_path)


def main() -> None:
    """Run standard and higher-order simulations and compare."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    logger.info("Higher-Order Klein-Gordon Dispersion Comparison")
    logger.info("=" * 60)

    grid, state_initial = build_grid_and_state()

    logger.info("Grid: %s points, bounds %s", grid.shape, grid.axes_bounds)
    logger.info("Initial: Gaussian pulse, width=5.0, center=100.0")

    # Run standard Klein-Gordon (no fourth-order term)
    storage_standard = run_simulation(
        state_initial.copy(),
        alpha_4=0.0,
        label="Standard Klein-Gordon (alpha_4=0)",
    )

    # Run with fourth-order dispersion
    storage_dispersive = run_simulation(
        state_initial.copy(),
        alpha_4=10,  # Fourth-order coefficient
        label="Higher-Order Klein-Gordon (alpha_4=10)",
    )

    # Create comparison plot
    create_comparison_plot(grid, storage_standard, storage_dispersive)

    logger.info("\nExpected differences:")
    logger.info("  - Standard: Clean propagation with minimal spreading")
    logger.info("  - Higher-order: Additional dispersion from fourth-order term")
    logger.info("  - Fourth-order term modifies high-k modes more strongly")


if __name__ == "__main__":
    main()
