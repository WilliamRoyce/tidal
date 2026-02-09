"""
2D Coupled Klein-Gordon Simulation Example.

This script demonstrates a coupled Klein-Gordon simulation in 2D with two fields
(phi0, phi1) having different masses and off-diagonal coupling. The fields start
with spatially separated Gaussian pulses and evolve, showing energy transfer
between fields due to coupling.

The simulation creates an animation showing the evolution of both fields and
their interaction over time.
"""

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)

from typing import TYPE_CHECKING

from tidal.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GridConfig,
    MultiFieldParams,
    SimulationConfig,
    make_coupled_kg_pde,
    make_grid,
    multi_gaussian_2d,
    run_with_snapshots,
)

if TYPE_CHECKING:
    from pde import CartesianGrid, FieldCollection


def build_grid_and_state() -> tuple[CartesianGrid, FieldCollection]:
    """
    Build a 2D grid and initial state with two spatially separated Gaussian pulses.

    Returns
    -------
    tuple[CartesianGrid, FieldCollection]
        Grid and initial state for coupled 2D Klein-Gordon simulation.
        Field 0 starts with amplitude 1.0 at (-15, 0), field 1 with amplitude 0.0 at (15, 0).
    """
    grid_config = GridConfig(
        dim=2,
        shape=(128, 128),
        bounds=((0.0, 200.0), (0.0, 200.0)),
        periodic=True,
    )
    grid = make_grid(grid_config)

    # Initialize two fields: field 0 excited on the left, field 1 dormant on the right
    state = multi_gaussian_2d(
        grid,
        amplitudes=[1.0, 1.0],
        widths=[4.0, 4.0],
        centers=[(100.0, 100.0), (125.0, 125.0)],
        velocities=[0.0, 0.0],
    )
    return grid, state


def main() -> None:
    """Run 2D coupled Klein-Gordon simulation and create dual field animation."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    grid, state = build_grid_and_state()

    # Create coupled PDE
    masses = [0.3, 0.7]  # lighter field (0) and heavier field (1)
    g = 0.2  # off-diagonal coupling strength
    coupling = [[0.0, g], [g, 0.0]]  # symmetric off-diagonal coupling
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    simulation_config = SimulationConfig(
        t_end=150.0,
        dt=None,
        solver="explicit",
        backend="numpy",
        progress=True,
    )

    logger.info("Starting 2D coupled KG simulation...")
    logger.info("  Masses: %s, Coupling strength: %s", masses, g)

    # Run simulation with automatic snapshot storage
    _result, storage = run_with_snapshots(
        pde=pde,
        state=state,
        config=simulation_config,
        snapshot_interval=1.0,
    )

    # Create dual field animation using AnimationBuilder
    builder = AnimationBuilder(storage, grid)
    config = AnimationConfig(
        output_path="outputs/coupled_2d_fields.mp4",
        title="Coupled 2D Klein-Gordon Fields",
        xlabel=r"$x$",
        ylabel=r"$y$",
        fps=30,
        dpi=150,
    )
    builder.create_dual_field_animation(config, field_labels=(r"$\phi_0$", r"$\phi_1$"))
    print("Saved outputs/coupled_2d_fields.[mp4|gif]")


if __name__ == "__main__":
    main()
