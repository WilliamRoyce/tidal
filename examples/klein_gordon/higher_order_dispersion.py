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
from torsion_gertsenshtein.kgsim.advanced_equations import HigherOrderKGPDE
from torsion_gertsenshtein.kgsim.animations import create_spacetime_plot_adjacent

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import CartesianGrid, FieldCollection


def build_grid_and_state() -> tuple[CartesianGrid, FieldCollection]:
    """Build 1D periodic grid with centered Gaussian pulse.

    Returns
    -------
    tuple[CartesianGrid, FieldCollection]
        Grid and initial state [phi, pi] with Gaussian pulse.
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


def run_simulation(
    grid: CartesianGrid,
    state: FieldCollection,
    alpha_4: float,
    label: str,
) -> list[tuple[float, np.ndarray]]:
    """Run simulation with specified fourth-order coefficient.

    Parameters
    ----------
    grid : CartesianGrid
        Simulation grid.
    state : FieldCollection
        Initial state.
    alpha_4 : float
        Fourth-order dispersion coefficient.
    label : str
        Label for progress output.

    Returns
    -------
    list[tuple[float, np.ndarray]]
        List of (time, phi_data) snapshots.
    """
    pde = HigherOrderKGPDE(
        mass=0.5,
        alpha_2=1.0,  # Standard Laplacian coefficient
        alpha_4=alpha_4,  # Fourth-order dispersion
        alpha_6=0.0,  # No sixth-order terms
    )

    # Note: Fourth-order terms require smaller time steps!
    # Rule of thumb: dt ~ O(dx^4) for fourth-order
    grid.discretization[0] if hasattr(grid, "discretization") else 0.2

    sim_config = SimulationConfig(
        t_end=200.0,
        dt=None,
        backend="numpy",  # Required for custom evolution_rate
        solver="scipy",
        method="RK45",
        progress=True,
    )

    recorder, snapshots = make_recorder()

    print(f"\n{label}:")
    print(f"  alpha_2 = {pde.alpha_2:.3f}, alpha_4 = {pde.alpha_4:.6f}")

    run(
        pde=pde,
        state=state,
        config=sim_config,
        extra_observer=recorder,
    )

    print(f"  Recorded {len(snapshots)} snapshots")
    return snapshots


def main() -> None:
    """Run standard and higher-order simulations and compare."""
    print("Higher-Order Klein-Gordon Dispersion Comparison")
    print("=" * 60)

    grid, state_initial = build_grid_and_state()
    print(f"Grid: {grid.shape} points, bounds {grid.axes_bounds}")
    print("Initial: Gaussian pulse, width=5.0, center=50.0")

    # Run standard Klein-Gordon (no fourth-order term)
    snapshots_standard = run_simulation(
        grid,
        state_initial.copy(),
        alpha_4=0.0,
        label="Standard Klein-Gordon (alpha_4=0)",
    )

    # Run with fourth-order dispersion
    snapshots_dispersive = run_simulation(
        grid,
        state_initial.copy(),
        alpha_4=10,  # Small fourth-order coefficient
        label="Higher-Order Klein-Gordon (alpha_4=10)",
    )

    # Create side-by-side spacetime comparison plot
    print("\nCreating comparison plot...")
    output_path = "outputs/higher_order_comparison.pdf"

    # Combine snapshots into the format expected by create_spacetime_plot_adjacent
    # It expects list of (time, field0, field1) tuples
    time_tolerance = 1e-6
    combined_snapshots: list[tuple[float, np.ndarray, np.ndarray]] = []
    for (t0, data0), (t1, data1) in zip(
        snapshots_standard, snapshots_dispersive, strict=False
    ):
        # Verify times match (they should, same t_end and similar dt)
        if abs(t0 - t1) > time_tolerance:
            print(f"Warning: Time mismatch at t0={t0:.3f}, t1={t1:.3f}")
        combined_snapshots.append((t0, data0, data1))

    create_spacetime_plot_adjacent(
        snapshots=combined_snapshots,
        grid=grid,
        output_path=output_path,
        titles=(
            r"Standard KG ($\alpha_4=0$)",
            r"Higher-Order KG ($\alpha_4=10$)",
        ),
        xlabel=r"$x$",
        ylabel=r"$t$",
        cbar_label=r"$\phi$",
        cmap="bwr",
        use_twoslope_norm=True,
        dpi=200,
    )

    print(f"Saved comparison plot: {output_path}")
    print("\nExpected differences:")
    print("  - Standard: Clean propagation with minimal spreading")
    print("  - Higher-order: Additional dispersion from fourth-order term")
    print("  - Fourth-order term modifies high-k modes more strongly")
    print("  - Wave packet may develop oscillatory tails")


if __name__ == "__main__":
    main()
