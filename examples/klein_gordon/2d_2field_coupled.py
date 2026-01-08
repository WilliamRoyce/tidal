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
import pathlib
from typing import TYPE_CHECKING, Any

logger = logging.getLogger(__name__)

# Disable PGF/LaTeX rendering for this example (causes issues with Unicode in titles)
import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    MultiFieldParams,
    SimulationConfig,
    make_coupled_kg_pde,
    make_grid,
    multi_gaussian_2d,
    run,
)

# Create animation using unified plotting module
from torsion_gertsenshtein.kgsim.animations import create_2d_coupled_animation

if TYPE_CHECKING:
    from collections.abc import Callable

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


def make_recorder() -> tuple[
    Callable[[FieldCollection, float], dict[str, Any]],
    list[tuple[float, np.ndarray, np.ndarray]],
]:
    """
    Create a recorder callback and storage for snapshots of both fields.

    Returns
    -------
    tuple[Callable, list]
        recorder: callback function to record field snapshots
        snapshots: list of (time, phi0_data, phi1_data) tuples
    """
    snapshots: list[tuple[float, np.ndarray, np.ndarray]] = []

    def record_fields(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        # Field order: phi0, pi0, phi1, pi1
        phi0_data = np.asarray(state_coll[0].data).copy()
        phi1_data = np.asarray(state_coll[2].data).copy()

        if not (np.isfinite(phi0_data).all() and np.isfinite(phi1_data).all()):
            msg = f"Non-finite field values at t={t}"
            raise RuntimeError(msg)

        snapshots.append((float(t), phi0_data, phi1_data))
        return {}

    return record_fields, snapshots


def main() -> None:
    """
    Run the 2D coupled Klein-Gordon simulation and create visualization.

    This function:
    1. Sets up a 2D grid with two coupled fields having different masses
    2. Initializes spatially separated Gaussian pulses
    3. Runs the simulation with off-diagonal coupling
    4. Creates an MP4 animation showing energy transfer between fields
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    # Build simulation components
    grid, state = build_grid_and_state()

    # Two fields with different masses and symmetric coupling
    masses = [0.3, 0.7]  # lighter field (0) and heavier field (1)
    g = 0.2  # off-diagonal bilinear coupling strength
    coupling = [[0.0, g], [g, 0.0]]
    pde = make_coupled_kg_pde(MultiFieldParams(masses=masses, coupling=coupling))

    # Simulation configuration
    # Use explicit solver with adaptive time-stepping for 2D
    simulation_config = SimulationConfig(
        t_end=150.0,
        dt=None,  # adaptive
        solver="explicit",
        backend="numpy",  # numpy backend for multi-field 2D
        progress=True,
    )

    # Create recorder
    recorder, snapshots = make_recorder()

    # Record initial condition
    snapshots.append(
        (
            0.0,
            np.asarray(state[0].data).copy(),
            np.asarray(state[2].data).copy(),
        )
    )

    logger.info("Starting 2D coupled KG simulation...")
    logger.info("  Grid: %s cells, domain: %s", grid.shape, grid.axes_bounds)
    logger.info("  Masses: %s", masses)
    logger.info("  Coupling strength: %s", g)
    logger.info("  Simulation time: 0 to %s", simulation_config.t_end)

    # Run simulation
    run(
        pde=pde,
        state=state,
        config=simulation_config,
        extra_observer=recorder,
    )

    logger.info("Simulation complete. Recorded %s snapshots.", len(snapshots))

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out_base = pathlib.Path("outputs") / "phi_evolution_2d_2field_coupled.mp4"

    logger.info("Creating animation: %s", out_base)
    create_2d_coupled_animation(
        snapshots,
        grid,
        out_base,
        titles=(r"Field 0: $\phi_0(x, y, t)$", r"Field 1: $\phi_1(x, y, t)$"),
        xlabel=r"$x$",
        ylabel=r"$y$",
        cbar_labels=(r"$\phi_0$", r"$\phi_1$"),
        cmap="bwr",
        fps=15,
        dpi=150,
    )

    logger.info("Animation saved to: %s", out_base)
    logger.info("Results:")
    logger.info("  Field 0 initially excited at (100, 100) with mass = %s", masses[0])
    logger.info("  Field 1 initially dormant at (150, 150) with mass = %s", masses[1])
    logger.info(
        "  Off-diagonal coupling g = %s causes energy transfer between fields", g
    )


if __name__ == "__main__":
    main()
