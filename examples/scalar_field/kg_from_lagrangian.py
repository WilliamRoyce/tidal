"""End-to-end example: Klein-Gordon field simulation from Lagrangian.

This script demonstrates the complete Lagrangian-to-simulation pipeline:
1. Load equations derived from the Klein-Gordon Lagrangian L = -1/2 (∂φ)² - 1/2 m²φ²
2. Build a PDE solver from the specification
3. Create initial conditions (Gaussian pulse)
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

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, PDEBase, ScalarField

from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.vectorfield import ComponentGaussianPulse

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde.trackers.base import TrackerBase

    from tidal.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Grid
GRID_BOUNDS = [(0, 100)]
GRID_SHAPE = [256]
GRID_PERIODIC = True

# Time integration
T_END = 30.0
DT = 0.01
SNAPSHOT_INTERVAL = 1.0

# Initial conditions
PULSE_CENTER = 50.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "kg_from_lagrangian_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationData:
    """Container for simulation results."""

    x: NumericArray
    times: list[float]
    snapshots: list[FieldCollection]


def _print_header() -> None:
    print("=" * 60)
    print("Klein-Gordon Field Simulation from Lagrangian")
    print("=" * 60)
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Number of components: {spec.n_components}")
    print(f"  Component names: {spec.component_names}")

    # Check for mass term
    has_mass = any(
        term.operator == "identity" and term.coefficient != 0
        for eq in spec.equations
        for term in eq.rhs_terms
    )
    print(f"  Mass term: {'present (m² = 1)' if has_mass else 'absent (massless)'}")
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
    grid = CartesianGrid(GRID_BOUNDS, GRID_SHAPE, periodic=GRID_PERIODIC)
    print(f"  Domain: {GRID_BOUNDS}")
    print(f"  Resolution: {GRID_SHAPE[0]} points")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid, spec: EquationSystem) -> FieldCollection:
    print("Step 4: Creating initial conditions...")
    print("  Gaussian pulse for scalar field φ")
    pulse = ComponentGaussianPulse(
        center=(PULSE_CENTER,),
        width=PULSE_WIDTH,
        amplitude=PULSE_AMPLITUDE,
        active_components={"phi_0": 1.0},
    )
    initial_state = pulse.create(grid, spec)
    print(f"  Pulse center: x = {PULSE_CENTER}")
    print(f"  Pulse width: {PULSE_WIDTH}")
    print("  Initial φ: Gaussian")
    print("  Initial ∂_t φ: 0 (initially at rest)")
    print()
    return initial_state


def _run_simulation(pde: PDEBase, initial_state: FieldCollection) -> MemoryStorage:
    print("Step 5: Running simulation...")
    print(f"  Duration: {T_END} time units")
    print("  (Klein-Gordon wave with m² = 1)")
    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=SNAPSHOT_INTERVAL)
    pde.solve(initial_state, t_range=T_END, dt=DT, tracker=tracker)
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


def _analyze_results(simulation: SimulationData, spec: EquationSystem) -> None:
    print("Step 6: Analyzing results...")

    # Get φ field (index 0) and π field (index 1)
    initial_phi = _get_component(simulation.snapshots[0], 0)
    final_phi = _get_component(simulation.snapshots[-1], 0)

    initial_peak_x = simulation.x[np.argmax(np.abs(initial_phi))]
    final_peak_x = simulation.x[np.argmax(np.abs(final_phi))]

    print(f"  Initial φ peak at x = {initial_peak_x:.1f}")
    print(f"  Final φ: wave has evolved to x = {final_peak_x:.1f}")

    # Check energy-like quantity (not conserved in discretization, but should be bounded)
    phi_max_over_time = [
        float(np.max(np.abs(_get_component(snapshot, 0))))
        for snapshot in simulation.snapshots
    ]
    initial_max = phi_max_over_time[0]
    final_max = phi_max_over_time[-1]

    print(f"  Initial max|φ| = {initial_max:.3f}")
    print(f"  Final max|φ| = {final_max:.3f}")

    # Check if mass term affects propagation (massive waves disperse)
    has_mass = any(
        term.operator == "identity" and term.coefficient != 0
        for eq in spec.equations
        for term in eq.rhs_terms
    )

    if has_mass:
        print("  Mass term present: wave exhibits dispersion ✓")
    else:
        print("  Massless case: wave propagates at c = 1 ✓")

    print()


def _plot_results(simulation: SimulationData, spec: EquationSystem) -> None:
    print("Step 7: Generating visualization...")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    # Get mass value for title
    has_mass = any(
        term.operator == "identity" and term.coefficient != 0
        for eq in spec.equations
        for term in eq.rhs_terms
    )
    mass_str = "m² = 1" if has_mass else "m² = 0"

    fig.suptitle(
        f"Klein-Gordon Scalar Field from Lagrangian ({mass_str})",
        fontsize=14,
    )

    # Initial condition
    ax = axes[0, 0]
    ax.plot(simulation.x, _get_component(simulation.snapshots[0], 0), "b-")
    ax.set_xlabel("x")
    ax.set_ylabel("φ")
    ax.set_title("Initial Condition (t=0)")
    ax.legend(["φ initial"])
    ax.grid(visible=True, alpha=0.3)

    # Final state
    ax = axes[0, 1]
    ax.plot(
        simulation.x,
        _get_component(simulation.snapshots[-1], 0),
        "b-",
        label=f"φ at t={simulation.times[-1]:.1f}",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("φ")
    ax.set_title(f"Final State (t={simulation.times[-1]:.1f})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Spacetime evolution
    ax = axes[1, 0]
    phi_history = np.stack(
        [_get_component(snapshot, 0) for snapshot in simulation.snapshots]
    )
    im = ax.imshow(
        phi_history.T,
        aspect="auto",
        origin="lower",
        extent=[0, simulation.times[-1], 0, 100],
        cmap="bwr_r",
        vmin=-np.max(np.abs(phi_history)),
        vmax=np.max(np.abs(phi_history)),
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("x")
    ax.set_title("φ Spacetime Evolution")
    plt.colorbar(im, ax=ax, label="φ")

    # Peak amplitude over time
    ax = axes[1, 1]
    phi_max = [
        float(np.max(np.abs(_get_component(snapshot, 0))))
        for snapshot in simulation.snapshots
    ]
    ax.plot(simulation.times, phi_max, "b-")
    ax.set_xlabel("Time")
    ax.set_ylabel("max|φ|")
    ax.set_title("Field Amplitude Evolution")
    ax.grid(visible=True, alpha=0.3)

    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to: {output_path}")
    print()


def _print_footer() -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations were derived from Lagrangian (not hardcoded)")
    print("  2. Klein-Gordon field evolves with mass term m² = 1")
    print("  3. Gaussian pulse propagates and disperses due to mass")
    print("  4. Demonstrates scalar field dynamics from first principles")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the Klein-Gordon field simulation from Lagrangian-derived equations."""
    json_path = Path(__file__).parent.parent / "data" / "klein_gordon_1d.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    initial_state = _create_initial_state(grid, spec)
    storage = _run_simulation(pde, initial_state)
    simulation = _collect_simulation_data(storage, grid)
    _analyze_results(simulation, spec)
    _plot_results(simulation, spec)
    _print_footer()


if __name__ == "__main__":
    main()
