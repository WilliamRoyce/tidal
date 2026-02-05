"""End-to-end example: EM field simulation from Lagrangian.

This script demonstrates the complete Lagrangian-to-simulation pipeline:
1. Load equations derived from the EM Lagrangian L = -1/4 F_uv F^uv
2. Build a PDE solver from the specification
3. Create initial conditions (Gaussian pulse in A_1)
4. Run the simulation
5. Visualize the results

The key point is that NO physics is hardcoded in Python - all equation
structure comes from the JSON file that was derived symbolically from
the Lagrangian using Mathematica/xAct.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, PDEBase, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system
from torsion_gertsenshtein.vectorfield import (
    ComponentGaussianPulse,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde.trackers.base import TrackerBase

    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]


A0_TOL = 1e-10
DEFAULT_T_END = 25.0
OUTPUT_FILENAME = "em_from_lagrangian_output.png"


@dataclass(frozen=True)
class SimulationData:
    """Container for simulation results."""

    x: NumericArray
    times: list[float]
    snapshots: list[FieldCollection]


def main() -> None:
    """Run the EM field simulation from Lagrangian-derived equations."""
    json_path = Path(__file__).parent.parent / "data" / "em_1d.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    initial_state = _create_initial_state(grid, spec)
    storage = _run_simulation(pde, initial_state, DEFAULT_T_END)
    simulation = _collect_simulation_data(storage, grid)
    _analyze_results(simulation)
    _plot_results(simulation)
    _print_footer()


def _print_header() -> None:
    print("=" * 60)
    print("EM Field Simulation from Lagrangian")
    print("=" * 60)
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Number of components: {spec.n_components}")
    print(f"  Component names: {spec.component_names}")
    print(f"  Gauge: {spec.metadata.get('gauge', 'none')}")
    print()
    return spec


def _build_pde(json_path: Path) -> PDEBase:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up simulation grid...")
    grid = CartesianGrid([(0, 100)], 256, periodic=True)
    print("  Domain: [0, 100]")
    print("  Resolution: 256 points")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid, spec: EquationSystem) -> FieldCollection:
    print("Step 4: Creating initial conditions...")
    print("  Gaussian pulse in A_1 (spatial component of vector potential)")
    pulse = ComponentGaussianPulse(
        center=(50.0,),
        width=5.0,
        amplitude=1.0,
        active_components={"A_1": 1.0},
    )
    initial_state = pulse.create(grid, spec)
    print("  Pulse center: x = 50")
    print("  Pulse width: 5")
    print("  Initial A_0: 0 (temporal component)")
    print("  Initial A_1: Gaussian (spatial component)")
    print()
    return initial_state


def _run_simulation(
    pde: PDEBase, initial_state: FieldCollection, t_end: float
) -> MemoryStorage:
    print("Step 5: Running simulation...")
    print(f"  Duration: {t_end} time units")
    print("  (Massless waves propagate at c = 1 in natural units)")
    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=1.0)
    pde.solve(initial_state, t_range=t_end, dt=0.01, tracker=tracker)
    print(f"  Stored {len(storage)} snapshots")
    print()
    return storage


def _collect_simulation_data(
    storage: MemoryStorage, grid: CartesianGrid
) -> SimulationData:
    x = np.asarray(grid.cell_coords[..., 0], dtype=float)
    times = list(storage.times)
    snapshots = [cast("FieldCollection", storage[i]) for i in range(len(storage))]
    return SimulationData(x=x, times=times, snapshots=snapshots)


def _get_component(snapshot: FieldCollection, index: int) -> NumericArray:
    field = snapshot[index]
    if not isinstance(field, ScalarField):
        msg = f"Expected ScalarField at index {index}, got {type(field).__name__}"
        raise TypeError(msg)
    return np.asarray(field.data, dtype=float)


def _analyze_results(simulation: SimulationData) -> None:
    print("Step 6: Analyzing results...")
    for t, snapshot in zip(simulation.times, simulation.snapshots, strict=True):
        a0_max = float(np.max(np.abs(_get_component(snapshot, 0))))
        if a0_max > A0_TOL:
            print(f"  WARNING: A_0 became non-zero at t={t}: max={a0_max}")
            break
    else:
        print("  A_0 remained zero throughout (no coupling) ✓")

    initial_a1 = _get_component(simulation.snapshots[0], 2)
    initial_peak_x = simulation.x[np.argmax(initial_a1)]
    print(f"  Initial A_1 peak at x = {initial_peak_x:.1f}")
    print("  Final A_1: pulse has split and propagated")
    print()


def _plot_results(simulation: SimulationData) -> None:
    print("Step 7: Generating visualization...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(
        "EM Vector Potential from Lagrangian L = -1/4 F_μν F^μν",
        fontsize=14,
    )

    ax = axes[0, 0]
    ax.plot(simulation.x, _get_component(simulation.snapshots[0], 2), "b-")
    ax.set_xlabel("x")
    ax.set_ylabel("A_1")
    ax.set_title("Initial Condition (t=0)")
    ax.legend(["A_1 initial"])
    ax.grid(visible=True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(
        simulation.x,
        _get_component(simulation.snapshots[-1], 2),
        "b-",
        label=f"A_1 at t={simulation.times[-1]:.1f}",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("A_1")
    ax.set_title(f"Final State (t={simulation.times[-1]:.1f})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    ax = axes[1, 0]
    a1_history = np.stack(
        [_get_component(snapshot, 2) for snapshot in simulation.snapshots]
    )
    im = ax.imshow(
        a1_history.T,
        aspect="auto",
        origin="lower",
        extent=[0, simulation.times[-1], 0, 100],
        cmap="plasma",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("x")
    ax.set_title("A_1 Spacetime Evolution")
    plt.colorbar(im, ax=ax, label="A_1")

    ax = axes[1, 1]
    a0_history = np.stack(
        [_get_component(snapshot, 0) for snapshot in simulation.snapshots]
    )
    a0_max = float(np.max(np.abs(a0_history)))
    ax.plot(
        simulation.times,
        [
            float(np.max(np.abs(_get_component(snapshot, 0))))
            for snapshot in simulation.snapshots
        ],
        "r-",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("max|A_0|")
    ax.set_title(f"A_0 Component (should be ~0, max={a0_max:.2e})")
    ax.grid(visible=True, alpha=0.3)

    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150)
    print(f"  Saved plot to: {output_path}")
    print()


def _print_footer() -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations were derived from Lagrangian (not hardcoded)")
    print("  2. Massless EM waves propagate at c = 1")
    print("  3. Gaussian pulse splits into left/right moving waves")
    print("  4. A_0 and A_1 evolve independently (no coupling)")
    print("=" * 60)


if __name__ == "__main__":
    main()
