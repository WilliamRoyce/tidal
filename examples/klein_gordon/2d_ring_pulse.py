"""2D Klein-Gordon Ring Pulse Simulation.

This example demonstrates creating a 2D ring-shaped initial condition using either:

1. **Function-based API**:
   ```python
   state = ring_pulse_2d(grid, amplitude=1.0, initial_radius=15.0, sigma=2.0)
   ```

2. **Class-based API**:
   ```python
   ic = RingPulse2D(amplitude=1.0, initial_radius=15.0, sigma=2.0)
   state = ic.build(grid)
   ```

The class-based approach offers reusability and easier customization
for scenarios involving multiple ring pulses or parameter sweeps.
"""

from __future__ import annotations

import logging
from typing import Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

from torsion_gertsenshtein.kgsim import (
    AnimationBuilder,
    AnimationConfig,
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    make_grid,
    ring_pulse_2d,
    run_with_snapshots,
)


def _build_simulation_components() -> dict[str, Any]:
    """Build simulation components for 2D ring pulse.

    Note: This example uses the function-based API (ring_pulse_2d).
    For the class-based alternative, use:
        ic = RingPulse2D(amplitude=1.0, initial_radius=15.0, sigma=2.0)
        state = ic.build(grid)
    """
    grid_config = GridConfig(
        dim=2,
        shape=(256, 256),
        bounds=((-50.0, 50.0), (-50.0, 50.0)),
        periodic=True,
    )
    grid = make_grid(grid_config)

    # Function-based API (convenient wrapper)
    state = ring_pulse_2d(grid, amplitude=1.0, initial_radius=15.0, sigma=2.0)

    simulation_config = SimulationConfig(
        t_end=50.0,
        dt=None,
        solver="explicit",
        method="RK45",
        backend="numpy",
        progress=True,
    )

    pde = KleinGordonPDE(params=KGParameters(mass=1.0))

    return {
        "grid": grid,
        "state": state,
        "pde": pde,
        "simulation_config": simulation_config,
    }


def main() -> None:
    """Run 2D Klein-Gordon simulation for ring pulse and create animation."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    simulation_components = _build_simulation_components()

    # Run simulation with automatic snapshot storage
    _result, storage = run_with_snapshots(
        pde=simulation_components["pde"],
        state=simulation_components["state"],
        config=simulation_components["simulation_config"],
        snapshot_interval=0.5,
    )

    # Create animation using AnimationBuilder
    builder = AnimationBuilder(storage, simulation_components["grid"])
    config = AnimationConfig(
        output_path="outputs/phi_evolution_2d.mp4",
        title=r"2D Klein-Gordon: $\phi(x,y)$",
        xlabel=r"$x$",
        ylabel=r"$y$",
        cbar_label=r"$\phi$",
        cmap="bwr",
    )
    builder.create_2d_heatmap_animation(config)
    print("Saved outputs/phi_evolution_2d.[mp4|gif]")


if __name__ == "__main__":
    main()
