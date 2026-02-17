"""End-to-end example: Scalar field scattering from a finite potential well.

This script demonstrates the background-field feature of the TIDAL pipeline:

    L = -1/2 (partial phi)^2 - V(x)/2 phi^2

where V(x) is a *background* scalar field (not varied in Euler-Lagrange),
spatially localized via Heaviside step functions:

    V(x) = V0  for 30 <= x <= 70   (inside the well)
    V(x) = 0   otherwise

The equation of motion is:

    d^2 phi / dt^2 = d^2 phi / dx^2 - V(x) phi

Key observations:
- Inside the well, V(x) = V0 acts as a mass term: omega^2 = k^2 + V0
- Outside the well, V(x) = 0: massless wave equation omega = |k|
- A pulse approaching the well boundary partially reflects and partially transmits
- Transmitted part disperses inside the well due to the mass gap
- Demonstrates coordinate-dependent (Piecewise) coefficient evaluation
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

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde.pdes.base import PDEBase
    from pde.trackers.base import TrackerBase

    from tidal.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
V0 = 4.0
WELL_LEFT = 30.0
WELL_RIGHT = 60.0

# Grid
GRID_BOUNDS = [(0, 100)]
GRID_SHAPE = [512]
GRID_PERIODIC = True

# Time integration
T_END = 100.0
DT = 0.005
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
PULSE_CENTER = 15.0
PULSE_WIDTH = 3.0
PULSE_K0 = 2.0  # Initial wavevector (rightward)

# Output
OUTPUT_FILENAME = "scalar_potential_well_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResults:
    """Container for simulation results."""

    x: NumericArray
    times: list[float]
    snapshots: list[FieldCollection]


def _print_header() -> None:
    print("=" * 60)
    print("Scalar Potential Well — Scattering Simulation")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/2 (d phi)^2 - V(x)/2 phi^2")
    print(f"  V(x) = V0  for {WELL_LEFT:.0f} <= x <= {WELL_RIGHT:.0f}")
    print("  V(x) = 0   otherwise")
    print(f"EOM:  d2_t(phi) = d2_x(phi) - V(x) phi  [V0={V0}]")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Components: {spec.component_names}")
    print(f"  Coordinates: {spec.coordinates}")

    for eq in spec.equations:
        print(f"  Equation for {eq.field_name}:")
        for term in eq.rhs_terms:
            sym = (
                f" (symbolic: {term.coefficient_symbolic})"
                if term.coefficient_symbolic
                else ""
            )
            dep = (
                f" [coord-dependent: {term.coordinate_dependent}]"
                if term.coordinate_dependent
                else ""
            )
            print(
                f"    {term.coefficient:+.1f} * {term.operator}({term.field}){sym}{dep}"
            )

    print()
    return spec


def _build_pde(json_path: Path) -> PDEBase:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters={"V0": V0})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Parameters: V0 = {V0}")
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


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")
    # Place Gaussian pulse OUTSIDE the well, give it rightward momentum
    print(f"  Gaussian pulse in phi at x={PULSE_CENTER:.0f}, width={PULSE_WIDTH:.0f}")
    print(f"  Initial momentum: k0={PULSE_K0:.0f} (rightward, toward well)")

    x = np.asarray(grid.cell_coords[..., 0], dtype=float)
    envelope = np.exp(-((x - PULSE_CENTER) ** 2) / (2 * PULSE_WIDTH**2))
    phi_data = envelope * np.cos(PULSE_K0 * (x - PULSE_CENTER))
    # pi_phi = d(phi)/dt at t=0 — for a rightward-propagating wavepacket:
    # pi = omega * envelope * sin(k0 * (x - x0)) approximated by k0 * ...
    pi_data = PULSE_K0 * envelope * np.sin(PULSE_K0 * (x - PULSE_CENTER))

    phi = ScalarField(grid, data=phi_data, label="phi_0")
    pi_phi = ScalarField(grid, data=pi_data, label="pi_phi_0")
    initial_state = FieldCollection([phi, pi_phi])
    print(f"  Created {len(initial_state)} fields: [phi_0, pi_phi_0]")
    print()
    return initial_state


def _run_simulation(pde: PDEBase, initial_state: FieldCollection) -> MemoryStorage:
    print("Step 5: Running simulation...")
    print(f"  Duration: {T_END} time units")
    print(f"  Scheme: runge-kutta, dt={DT}")
    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=SNAPSHOT_INTERVAL)
    pde.solve(
        initial_state, t_range=T_END, dt=DT, scheme="runge-kutta", tracker=tracker
    )
    print(f"  Stored {len(storage)} snapshots")
    print()
    return storage


def _collect_simulation_data(
    storage: MemoryStorage,
    grid: CartesianGrid,
) -> SimulationResults:
    x = np.asarray(grid.cell_coords[..., 0], dtype=float)
    times = list(storage.times)
    snapshots = [cast("FieldCollection", storage[i]) for i in range(len(storage))]
    return SimulationResults(x=x, times=times, snapshots=snapshots)


def _get_field(snapshot: FieldCollection, index: int) -> NumericArray:
    field = snapshot[index]
    if not isinstance(field, ScalarField):
        msg = f"Expected ScalarField at index {index}, got {type(field).__name__}"
        raise TypeError(msg)
    return np.asarray(field.data, dtype=float)


def _analyze_results(sim: SimulationResults, _spec: EquationSystem) -> None:
    print("Step 6: Analyzing results...")

    initial_phi = _get_field(sim.snapshots[0], 0)
    final_phi = _get_field(sim.snapshots[-1], 0)

    print(f"  Initial max|phi| = {np.max(np.abs(initial_phi)):.4f}")
    print(f"  Final   max|phi| = {np.max(np.abs(final_phi)):.4f}")

    omega_min = np.sqrt(V0)
    print(f"  V0 = {V0}, omega_min = sqrt(V0) = {omega_min:.2f}")
    print(f"  Well region: x in [{WELL_LEFT:.0f}, {WELL_RIGHT:.0f}]")
    print("  Pulse scatters off well boundary: partial reflection + transmission")

    # Check energy inside vs outside the well at final time
    x = sim.x
    inside = (x >= WELL_LEFT) & (x <= WELL_RIGHT)
    energy_inside = float(np.sum(final_phi[inside] ** 2))
    energy_outside = float(np.sum(final_phi[~inside] ** 2))
    total = energy_inside + energy_outside
    if total > 0:
        print(
            f"  Final |phi|^2 fraction: "
            f"{energy_inside / total:.1%} inside, {energy_outside / total:.1%} outside"
        )

    print()


def _plot_results(sim: SimulationResults) -> None:
    print("Step 7: Generating visualization...")
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))

    fig.suptitle(
        "Scalar Field Scattering from Finite Potential Well\n"
        r"$\mathcal{L} = -\frac{1}{2}(\partial\phi)^2"
        r" - \frac{V(x)}{2}\phi^2$"
        r"$,\quad V(x) = V_0\,\theta(x-30)\,\theta(60-x)$",
        fontsize=13,
    )

    # (0,0): Spacetime evolution
    ax = axes[0, 0]
    phi_history = np.stack([_get_field(s, 0) for s in sim.snapshots])
    vmax = np.max(np.abs(phi_history))
    im = ax.imshow(
        phi_history.T,
        aspect="auto",
        origin="lower",
        extent=[0, sim.times[-1], 0, 100],
        cmap="bwr_r",
        vmin=-vmax,
        vmax=vmax,
    )
    # Mark well boundaries
    ax.axhline(WELL_LEFT, color="k", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(WELL_RIGHT, color="k", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("Time")
    ax.set_ylabel("x")
    ax.set_title(r"$\phi(t, x)$ Spacetime Evolution")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    # (0,1): Field profiles at selected times
    ax = axes[0, 1]
    n_snaps = len(sim.snapshots)
    indices = [0, n_snaps // 4, n_snaps // 2, 3 * n_snaps // 4, n_snaps - 1]
    colors = ["b", "g", "orange", "r", "purple"]
    for idx, color in zip(indices, colors, strict=True):
        ax.plot(
            sim.x,
            _get_field(sim.snapshots[idx], 0),
            color=color,
            linewidth=1.2,
            label=rf"$t={sim.times[idx]:.0f}$",
        )
    # Shade the well region
    ax.axvspan(WELL_LEFT, WELL_RIGHT, alpha=0.08, color="gray")
    ax.axvline(WELL_LEFT, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(WELL_RIGHT, color="k", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title("Field Profiles")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(visible=True, alpha=0.3)

    # (1,0): Peak amplitude over time
    ax = axes[1, 0]
    phi_max = [float(np.max(np.abs(_get_field(s, 0)))) for s in sim.snapshots]
    ax.plot(sim.times, phi_max, "b-", linewidth=1.5)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"max$|\phi|$")
    ax.set_title("Peak Amplitude Evolution")
    ax.grid(visible=True, alpha=0.3)

    # (1,1): Physics info panel
    ax = axes[1, 1]
    omega_min = np.sqrt(V0)
    info_lines = [
        ("System Parameters:", 0.92, 12, "bold", "black"),
        (
            rf"$V_0 = {V0:.1f}$,  well: $x \in [{WELL_LEFT:.0f}, {WELL_RIGHT:.0f}]$",
            0.80,
            11,
            "normal",
            "black",
        ),
        (
            rf"$\omega_{{\min}} = \sqrt{{V_0}} = {omega_min:.1f}$ (inside well)",
            0.68,
            11,
            "normal",
            "black",
        ),
        ("Physics:", 0.52, 12, "bold", "black"),
        (
            r"Outside: $\omega = |k|$ (massless)",
            0.40,
            10,
            "normal",
            "black",
        ),
        (
            r"Inside: $\omega^2 = k^2 + V_0$ (massive)",
            0.28,
            10,
            "normal",
            "black",
        ),
        (
            "Pulse reflects/transmits at well boundary",
            0.14,
            10,
            "normal",
            "gray",
        ),
    ]
    for text, y, size, weight, color in info_lines:
        ax.text(
            0.5,
            y,
            text,
            transform=ax.transAxes,
            fontsize=size,
            ha="center",
            weight=weight,
            color=color,
        )
    ax.axis("off")

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
    print("  1. V(x) is a spatially-localized background (Piecewise in JSON)")
    print("  2. Python evaluates Piecewise[...] via _mathematica_to_python -> np.where")
    print("  3. Pulse partially reflects at x=30, partially transmits into well")
    print("  4. Inside well: massive dispersion. Outside: massless propagation")
    print("  5. Demonstrates [[background_fields]] + coordinate_dependent coefficients")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the scalar potential well scattering simulation."""
    json_path = Path(__file__).parent.parent / "data" / "scalar_potential_well.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    initial_state = _create_initial_state(grid)
    storage = _run_simulation(pde, initial_state)
    sim = _collect_simulation_data(storage, grid)
    _analyze_results(sim, spec)
    _plot_results(sim)
    _print_footer()


if __name__ == "__main__":
    main()
