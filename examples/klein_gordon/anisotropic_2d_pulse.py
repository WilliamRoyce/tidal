"""Anisotropic Klein-Gordon: 2D Gaussian pulse with directional wave speeds.

This example demonstrates anisotropic wave propagation where the wave speed
differs in x and y directions. The initial Gaussian pulse evolves into an
elliptical wavefront, illustrating the direction-dependent propagation.

Comparing isotropic (c_x = c_y = 1.0) vs anisotropic (c_x = 2.0, c_y = 0.5)
cases shows how anisotropy affects wavefront shape and arrival times.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")

import numpy as np

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
)
from torsion_gertsenshtein.kgsim.advanced_equations import AnisotropicKGPDE
from torsion_gertsenshtein.kgsim.animations import create_2d_heatmap_animation

if TYPE_CHECKING:
    from collections.abc import Callable

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


def make_recorder() -> tuple[
    Callable[[FieldCollection, float], dict[str, Any]],
    list[tuple[float, np.ndarray]],
]:
    """Create recorder callback for snapshots.

    Returns
    -------
    tuple[callable, list]
        Recorder function and snapshots list containing (time, phi_data) tuples.
    """
    snapshots: list[tuple[float, np.ndarray]] = []

    def record_phi(state_coll: FieldCollection, t: float) -> dict[str, Any]:
        snapshots.append((float(t), np.asarray(state_coll[0].data).copy()))
        return {}

    return record_phi, snapshots


def main() -> None:
    """Run anisotropic Klein-Gordon simulation and create animation."""
    grid, state = build_grid_and_state()

    # Anisotropic PDE: wave propagates faster in x (c_x=2.0) than y (c_y=0.5)
    # This creates elliptical wavefronts stretched along the x-direction
    pde = AnisotropicKGPDE(
        mass=0.5,
        speeds=[2.0, 0.5],  # [c_x, c_y] - anisotropy ratio of 4:1
    )

    sim_config = SimulationConfig(
        t_end=100.0,
        dt=None,
        backend="numba",  # Numba-accelerated for improved performance
        solver="scipy",
        method="RK45",
        progress=True,
    )

    recorder, snapshots = make_recorder()

    print("Running anisotropic Klein-Gordon simulation...")
    print(f"  Grid: {grid.shape} points, bounds {grid.axes_bounds}")
    print(f"  Anisotropic speeds: c_x = {pde.speeds[0]:.1f}, c_y = {pde.speeds[1]:.1f}")
    print(f"  Mass: m = {np.sqrt(pde.m2):.2f}")
    print(f"  Time range: 0 to {sim_config.t_end}")

    run(
        pde=pde,
        state=state,
        config=sim_config,
        extra_observer=recorder,
    )

    print(f"Simulation complete. Recorded {len(snapshots)} snapshots.")

    # Create animation showing elliptical wavefront expansion
    output_path = "outputs/anisotropic_2d_pulse.mp4"
    print(f"\nCreating animation: {output_path}")

    create_2d_heatmap_animation(
        snapshots=snapshots,
        grid=grid,
        output_path=output_path,
        title_template="Anisotropic KG (c_x=2.0, c_y=0.5): t={t:.1f}",
        cmap="bwr",
        fps=15,
    )

    print("Animation saved successfully.")
    print("\nExpected behavior:")
    print("  - Wavefront expands as an ellipse (not a circle)")
    print("  - Faster propagation along x-axis (horizontal)")
    print("  - Slower propagation along y-axis (vertical)")
    print("  - Anisotropy ratio c_x/c_y = 4.0 determines ellipse aspect ratio")


if __name__ == "__main__":
    main()
