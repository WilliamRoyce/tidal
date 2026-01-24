"""Anisotropic Klein-Gordon: 2D Gaussian pulse with directional wave speeds.

This example demonstrates anisotropic wave propagation where the wave speed
differs in x and y directions. The initial Gaussian pulse evolves into an
elliptical wavefront, illustrating the direction-dependent propagation.

Comparing isotropic (c_x = c_y = 1.0) vs anisotropic (c_x = 2.0, c_y = 0.5)
cases shows how anisotropy affects wavefront shape and arrival times.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")

from torsion_gertsenshtein.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GridConfig,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run_with_snapshots,
)
from torsion_gertsenshtein.kgsim.advanced_equations import AnisotropicKGPDE

if TYPE_CHECKING:
    from pde import CartesianGrid, FieldCollection


def build_grid_and_state() -> tuple[CartesianGrid, FieldCollection]:
    """Build 2D grid with centered Gaussian pulse.

    Returns
    -------
    tuple[CartesianGrid, FieldCollection]
        Grid and initial state [phi, pi] with Gaussian pulse at origin.
    """
    grid_config = GridConfig(
        dim=2,
        shape=(256, 256),
        bounds=((0.0, 200.0), (0.0, 200.0)),
        periodic=True,
    )
    grid = make_grid(grid_config)
    state = gaussian_pulse(
        grid,
        amplitude=1.0,
        width=5.0,
        center=[100.0, 100.0],
    )
    return grid, state


def main() -> None:
    """Run anisotropic Klein-Gordon simulation and create animation."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    grid, state = build_grid_and_state()

    # Anisotropic PDE: wave propagates faster in x (c_x=2.0) than y (c_y=0.5)
    pde = AnisotropicKGPDE(mass=0.5, speeds=[2.0, 0.5])

    simulation_config = SimulationConfig(
        t_end=100.0,
        dt=None,
        backend="numba",
        solver="scipy",
        method="RK45",
        progress=True,
    )

    # Run simulation with automatic snapshot storage
    _result, storage = run_with_snapshots(
        pde=pde,
        state=state,
        config=simulation_config,
        snapshot_interval=1.0,
    )

    # Create animation using AnimationBuilder
    builder = AnimationBuilder(storage, grid)
    config = AnimationConfig(
        output_path="outputs/anisotropic_2d_pulse.mp4",
        title=r"Anisotropic KG: $\phi(x,y)$ at $t =$",
        xlabel=r"$x$",
        ylabel=r"$y$",
        cbar_label=r"$\phi$",
        fps=30,
        dpi=150,
    )
    builder.create_2d_heatmap_animation(config)
    print("Saved outputs/anisotropic_2d_pulse.[mp4|gif]")


if __name__ == "__main__":
    main()
