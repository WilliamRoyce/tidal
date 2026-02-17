"""Proca field (massive vector) 1+1D simulation from Lagrangian-derived equations.

This script demonstrates the Proca theory (massive electrodynamics) in 1+1D:
  L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a

The mass term breaks gauge invariance and gives the photon mass m.
This leads to dispersive wave propagation with:
  - Dispersion relation: omega^2 = k^2 + m^2
  - Group velocity: v_g = k/omega < 1 (subluminal)
  - Phase velocity: v_p = omega/k > 1 (superluminal)

Equations of motion (with automatic Lorenz condition):
  d2_t(A_0) = d2_x(A_0) - m^2 A_0  (massive Klein-Gordon)
  d2_t(A_1) = d2_x(A_1) - m^2 A_1  (massive Klein-Gordon)

Unlike Maxwell where pulses maintain shape, Proca pulses spread due to dispersion.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from tidal.symbolic import build_pde_from_json, load_equation_system
from tidal.utils import normalize_solve_result

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
MASS_SQUARED = 1.0

# Grid
X_MIN = 0.0
X_MAX = 100.0
GRID_POINTS = 512

# Time integration
T_END = 30.0
DT = 0.005  # Smaller timestep for better energy conservation
SNAPSHOT_INTERVAL = 1.0

# Initial conditions
PULSE_CENTER_X = 50.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Analysis
DISPERSION_THRESHOLD = 1.2

# Output
OUTPUT_FILENAME = "proca_output.png"

# State indices (2 components, each with field + momentum)
IDX_A0 = 0
IDX_PI0 = 1
IDX_A1 = 2
IDX_PI1 = 3

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for Proca simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    x_coords: NumericArray


def _print_header() -> None:
    print("=" * 60)
    print("Proca Field (Massive Vector) 1+1D Simulation")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Dimension: {spec.dimension} (1+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    # Show equation structure
    print()
    print("  Equation structure:")
    for eq in spec.equations:
        terms_str = " + ".join(
            f"{term.coefficient:.1f}*{term.operator}({term.field})"
            for term in eq.rhs_terms
        )
        print(f"    d2_t({eq.field_name}) = {terms_str}")
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    # The mass parameter (procaMassSquared) is stored symbolically in JSON
    # and can be set at simulation time
    pde = build_pde_from_json(json_path, parameters={"procaMassSquared": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED} (runtime-specified)")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 1D simulation grid...")
    grid = CartesianGrid(
        bounds=[(X_MIN, X_MAX)],
        shape=[GRID_POINTS],
        periodic=True,
    )
    print(f"  Domain: [{X_MIN}, {X_MAX}]")
    print(f"  Resolution: {grid.shape[0]} points")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    # Initial state: 4 fields = 2 components * (field + momentum)
    # State layout: procaA_0, pi_0, procaA_1, pi_1
    a0 = ScalarField(grid, data=0.0, label="procaA_0")
    pi0 = ScalarField(grid, data=0.0, label="pi_0")
    a1 = ScalarField(grid, data=0.0, label="procaA_1")
    pi1 = ScalarField(grid, data=0.0, label="pi_1")

    # Initialize a Gaussian pulse in A_1 (spatial component)
    # This demonstrates dispersive propagation for massive vector
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    gaussian = PULSE_AMPLITUDE * np.exp(
        -((x - PULSE_CENTER_X) ** 2) / (2 * PULSE_WIDTH**2)
    )
    a1.data[:] = gaussian

    state = FieldCollection([a0, pi0, a1, pi1])
    print(f"  Gaussian pulse in A_1 at center x={PULSE_CENTER_X}")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print("  A_0 and momenta initialized to zero")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")
    # Note: Using smaller dt for better energy conservation.
    # The explicit Euler integrator does not exactly conserve the Hamiltonian;
    # energy drift decreases with smaller dt (dt=0.01 gives ~36% drift,
    # dt=0.005 gives ~16%).

    storage = MemoryStorage()
    result = pde.solve(  # type: ignore[union-attr]
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",  # RK4 for better energy conservation
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    x = cast("np.ndarray", grid.cell_coords[..., 0])

    print(f"  Duration: {T_END} time units")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        x_coords=x,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])
    x = result.x_coords

    # Extract field amplitudes
    initial_a0 = np.max(np.abs(initial[IDX_A0].data))
    initial_a1 = np.max(np.abs(initial[IDX_A1].data))

    final_a0 = np.max(np.abs(final[IDX_A0].data))
    final_a1 = np.max(np.abs(final[IDX_A1].data))

    print(f"  Initial: max|A_0| = {initial_a0:.4f}, max|A_1| = {initial_a1:.4f}")
    print(f"  Final:   max|A_0| = {final_a0:.4f}, max|A_1| = {final_a1:.4f}")

    # Check for dispersion (pulse spreading)
    initial_width = np.sqrt(
        np.sum((x - PULSE_CENTER_X) ** 2 * initial[IDX_A1].data ** 2)
        / np.sum(initial[IDX_A1].data ** 2)
    )
    # Find center of mass of final A_1
    final_com = np.sum(x * final[IDX_A1].data ** 2) / (
        np.sum(final[IDX_A1].data ** 2) + 1e-10
    )
    final_width = np.sqrt(
        np.sum((x - final_com) ** 2 * final[IDX_A1].data ** 2)
        / (np.sum(final[IDX_A1].data ** 2) + 1e-10)
    )

    print(f"  Initial pulse width: {initial_width:.2f}")
    print(f"  Final pulse width: {final_width:.2f}")

    if final_width > initial_width * DISPERSION_THRESHOLD:
        print("  Dispersion observed: pulse has spread (as expected for massive field)")
    else:
        print("  Note: Pulse spread may be subtle at this mass/time scale")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    x = result.x_coords

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Get time snapshots
    times = [0, len(storage) // 3, 2 * len(storage) // 3, len(storage) - 1]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times)))

    # A_0 evolution
    ax = axes[0, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            x,
            snapshot[IDX_A0].data,
            color=colors[i],
            label=f"t={t_idx:.0f}",
            alpha=0.8,
        )
    ax.set_xlabel("x")
    ax.set_ylabel("A_0")
    ax.set_title("A_0 (temporal component) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # A_1 evolution
    ax = axes[0, 1]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            x,
            snapshot[IDX_A1].data,
            color=colors[i],
            label=f"t={t_idx:.0f}",
            alpha=0.8,
        )
    ax.set_xlabel("x")
    ax.set_ylabel("A_1")
    ax.set_title("A_1 (spatial component) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Momentum pi_1 evolution
    ax = axes[1, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            x,
            snapshot[IDX_PI1].data,
            color=colors[i],
            label=f"t={t_idx:.0f}",
            alpha=0.8,
        )
    ax.set_xlabel("x")
    ax.set_ylabel("\u03c0_1")
    ax.set_title("\u03c0_1 (momentum of A_1) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Energy (Hamiltonian) over time
    # Note: The Hamiltonian should be conserved for the continuous Klein-Gordon
    # equation. Using RK4 scheme (scheme='runge-kutta') provides much better
    # energy conservation than the default explicit Euler, which is
    # non-symplectic and causes linear drift.
    ax = axes[1, 1]
    energies = []
    time_values = []

    dx = float(x[1] - x[0])

    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])

        # Extract fields and momenta
        a0_data = snapshot[IDX_A0].data
        pi0_data = snapshot[IDX_PI0].data
        a1_data = snapshot[IDX_A1].data
        pi1_data = snapshot[IDX_PI1].data

        # Compute spatial gradients
        grad_a0 = np.gradient(a0_data, dx)
        grad_a1 = np.gradient(a1_data, dx)

        # Klein-Gordon Hamiltonian: H = 1/2 int [pi^2 + (d_x phi)^2 + m^2 phi^2] dx
        energy_a0 = np.sum(
            0.5 * pi0_data**2 + 0.5 * grad_a0**2 + 0.5 * MASS_SQUARED * a0_data**2
        )
        energy_a1 = np.sum(
            0.5 * pi1_data**2 + 0.5 * grad_a1**2 + 0.5 * MASS_SQUARED * a1_data**2
        )

        energies.append(energy_a0 + energy_a1)
        time_values.append(t_idx)

    ax.plot(time_values, energies, "b-", linewidth=2)
    ax.set_xlabel("Snapshot index")
    ax.set_ylabel("Total Energy")
    ax.set_title(
        r"Hamiltonian $H = \frac{1}{2}\int[\pi^2 + (\nabla\phi)^2 + m^2\phi^2]dx$"
    )
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        f"Proca Field 1+1D: Massive Vector (m^2={MASS_SQUARED})\n"
        r"Dispersion relation: $\omega^2 = k^2 + m^2$",
        fontsize=12,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved plot to: {output_path}")
    plt.close()
    print()


def _print_footer() -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations derived from Proca Lagrangian (Maxwell + mass term)")
    print("  2. Each component satisfies massive Klein-Gordon equation")
    print("  3. Dispersion relation: omega^2 = k^2 + m^2 (subluminal propagation)")
    print("  4. Unlike massless Maxwell, pulses spread due to dispersion")
    print("  5. Lorenz condition d_mu A^mu = 0 is automatic for m != 0")
    print(
        f"  6. Mass parameter m^2={MASS_SQUARED} specified at runtime (not pre-evaluated)"
    )
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the Proca field simulation."""
    json_path = Path(__file__).parent.parent / "data" / "proca_1d.json"
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid)
    result = _run_simulation(pde, grid, state)
    _analyze_results(result)
    _plot_results(result)
    _print_footer()


if __name__ == "__main__":
    main()
