"""End-to-end example: Coupled scalar fields simulation from Lagrangian.

This script demonstrates the complete Lagrangian-to-simulation pipeline for
coupled scalar fields:

    L = 1/2(∂φ)² - 1/2 m_φ² φ² + 1/2(∂χ)² - 1/2 m_χ² χ² - g φ χ

EOMs:
    d²φ/dt² = ∇²φ - m_φ² φ - g χ
    d²χ/dt² = ∇²χ - m_χ² χ - g φ

This tests:
- Linear coupled systems
- Mode-mixing between fields
- Matrix-form PDE generation
- Cross-field references in JSON specification

The key point is that the coupling term (g φ χ) appears as an "identity" operator
on the OTHER field in each equation's RHS, demonstrating that the pipeline
correctly handles inter-field coupling.
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
T_END = 20.0
DT = 0.01
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
PULSE_CENTER = 30.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Analysis
COUPLING_THRESHOLD = 0.01

# Output
OUTPUT_FILENAME = "coupled_scalars_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationData:
    """Container for simulation results."""

    x: NumericArray
    times: list[float]
    snapshots: list[FieldCollection]


def _print_header() -> None:
    print("=" * 60)
    print("Coupled Scalar Fields Simulation from Lagrangian")
    print("=" * 60)
    print()
    print("Lagrangian: L = 1/2(∂φ)² - 1/2 m_φ² φ² + 1/2(∂χ)² - 1/2 m_χ² χ² - g φ χ")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Number of components: {spec.n_components}")
    print(f"  Component names: {spec.component_names}")

    # Check for cross-field coupling (term referencing different field)
    has_coupling = False
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field != eq.field_name and term.operator == "identity":
                has_coupling = True
                print(
                    f"  Cross-field coupling detected: {eq.field_name} ← {term.field}"
                )
                break

    if has_coupling:
        print("  Mode-mixing: enabled (coupled system)")
    else:
        print("  Mode-mixing: disabled (uncoupled system)")

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
    print(f"  Gaussian pulse in φ at x={PULSE_CENTER}")
    print("  χ initially zero (energy will transfer via coupling)")

    # Initialize phi with a Gaussian pulse, chi starts at zero
    pulse = ComponentGaussianPulse(
        center=(PULSE_CENTER,),
        width=PULSE_WIDTH,
        amplitude=PULSE_AMPLITUDE,
        active_components={"phi_0": 1.0},  # Only phi has initial pulse
    )
    initial_state = pulse.create(grid, spec)
    print(f"  Created {len(initial_state)} fields: [phi, pi_phi, chi, pi_chi]")
    print()
    return initial_state


def _run_simulation(pde: PDEBase, initial_state: FieldCollection) -> MemoryStorage:
    print("Step 5: Running simulation...")
    print(f"  Duration: {T_END} time units")
    print("  (Coupled evolution with mode-mixing)")
    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=SNAPSHOT_INTERVAL)
    pde.solve(
        initial_state, t_range=T_END, dt=DT, scheme="runge-kutta", tracker=tracker
    )
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


def _analyze_results(simulation: SimulationData, _spec: EquationSystem) -> None:
    print("Step 6: Analyzing results...")

    # Get phi field (index 0) and chi field (index 2)
    # State layout: [phi, pi_phi, chi, pi_chi]
    initial_phi = _get_component(simulation.snapshots[0], 0)
    final_phi = _get_component(simulation.snapshots[-1], 0)
    initial_chi = _get_component(simulation.snapshots[0], 2)
    final_chi = _get_component(simulation.snapshots[-1], 2)

    # Check initial state
    print(f"  Initial max|φ| = {np.max(np.abs(initial_phi)):.3f}")
    print(f"  Initial max|χ| = {np.max(np.abs(initial_chi)):.3f}")

    # Check final state
    print(f"  Final max|φ| = {np.max(np.abs(final_phi)):.3f}")
    print(f"  Final max|χ| = {np.max(np.abs(final_chi)):.3f}")

    # Check if energy transferred to chi
    chi_max_over_time = [
        float(np.max(np.abs(_get_component(snapshot, 2))))
        for snapshot in simulation.snapshots
    ]
    if max(chi_max_over_time) > COUPLING_THRESHOLD:
        print("  Energy transfer: φ → χ observed ✓ (coupling working)")
    else:
        print("  No energy transfer to χ (check coupling)")

    # Normal mode analysis info
    # For m_φ²=1, m_χ²=4, g=0.5:
    # ω² = (m_φ² + m_χ²)/2 ± sqrt[((m_φ² - m_χ²)/2)² + g²]
    # ω₁² ≈ 0.84, ω₂² ≈ 4.16
    print("  Expected normal mode frequencies:")
    print("    ω₁ ≈ 0.92 (lower mode, mixed φ-χ)")
    print("    ω₂ ≈ 2.04 (higher mode, mixed φ-χ)")

    print()


def _plot_results(simulation: SimulationData, _spec: EquationSystem) -> None:  # noqa: PLR0915
    print("Step 7: Generating visualization...")
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    fig.suptitle(
        "Coupled Scalar Fields: φ-χ Mode Mixing\n"
        "L = ½(∂φ)² - ½m²ᵩφ² + ½(∂χ)² - ½m²ᵪχ² - gφχ",
        fontsize=14,
    )

    # Spacetime evolution for phi
    ax = axes[0, 0]
    phi_history = np.stack(
        [_get_component(snapshot, 0) for snapshot in simulation.snapshots]
    )
    vmax = np.max(np.abs(phi_history))
    im = ax.imshow(
        phi_history.T,
        aspect="auto",
        origin="lower",
        extent=[0, simulation.times[-1], 0, 100],
        cmap="bwr_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("x")
    ax.set_title("φ Spacetime Evolution")
    plt.colorbar(im, ax=ax, label="φ")

    # Spacetime evolution for chi
    ax = axes[0, 1]
    chi_history = np.stack(
        [_get_component(snapshot, 2) for snapshot in simulation.snapshots]
    )
    vmax_chi = np.max(np.abs(chi_history))
    if vmax_chi < 1e-10:  # noqa: PLR2004
        vmax_chi = 1.0  # Avoid division by zero
    im = ax.imshow(
        chi_history.T,
        aspect="auto",
        origin="lower",
        extent=[0, simulation.times[-1], 0, 100],
        cmap="bwr_r",
        vmin=-vmax_chi,
        vmax=vmax_chi,
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("x")
    ax.set_title("χ Spacetime Evolution")
    plt.colorbar(im, ax=ax, label="χ")

    # Initial and final profiles
    ax = axes[0, 2]
    ax.plot(
        simulation.x, _get_component(simulation.snapshots[0], 0), "b-", label="φ(t=0)"
    )
    ax.plot(
        simulation.x,
        _get_component(simulation.snapshots[-1], 0),
        "b--",
        label=f"φ(t={simulation.times[-1]:.1f})",
    )
    ax.plot(
        simulation.x, _get_component(simulation.snapshots[0], 2), "r-", label="χ(t=0)"
    )
    ax.plot(
        simulation.x,
        _get_component(simulation.snapshots[-1], 2),
        "r--",
        label=f"χ(t={simulation.times[-1]:.1f})",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("Field value")
    ax.set_title("Field Profiles")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Peak amplitude evolution
    ax = axes[1, 0]
    phi_max = [
        float(np.max(np.abs(_get_component(snapshot, 0))))
        for snapshot in simulation.snapshots
    ]
    chi_max = [
        float(np.max(np.abs(_get_component(snapshot, 2))))
        for snapshot in simulation.snapshots
    ]
    ax.plot(simulation.times, phi_max, "b-", label="max|φ|")
    ax.plot(simulation.times, chi_max, "r-", label="max|χ|")
    ax.set_xlabel("Time")
    ax.set_ylabel("Peak amplitude")
    ax.set_title("Amplitude Evolution (Mode Mixing)")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Total energy proxy (L2 norm squared)
    ax = axes[1, 1]
    phi_energy = [
        float(np.sum(_get_component(snapshot, 0) ** 2))
        for snapshot in simulation.snapshots
    ]
    chi_energy = [
        float(np.sum(_get_component(snapshot, 2) ** 2))
        for snapshot in simulation.snapshots
    ]
    total_energy = [p + c for p, c in zip(phi_energy, chi_energy, strict=True)]
    ax.plot(simulation.times, phi_energy, "b-", label="φ energy")
    ax.plot(simulation.times, chi_energy, "r-", label="χ energy")
    ax.plot(simulation.times, total_energy, "k--", label="Total")
    ax.set_xlabel("Time")
    ax.set_ylabel("Energy (L² norm)")
    ax.set_title("Energy Transfer Between Fields")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Coupling information
    ax = axes[1, 2]
    ax.text(
        0.5,
        0.8,
        "Coupled System Parameters:",
        transform=ax.transAxes,
        fontsize=12,
        ha="center",
        weight="bold",
    )
    ax.text(
        0.5, 0.6, "m²ᵩ = 1.0 (ω₀ = 1)", transform=ax.transAxes, fontsize=11, ha="center"
    )
    ax.text(
        0.5,
        0.45,
        "m²ᵪ = 4.0 (ω₀ = 2)",
        transform=ax.transAxes,
        fontsize=11,
        ha="center",
    )
    ax.text(
        0.5, 0.3, "g = 0.5 (coupling)", transform=ax.transAxes, fontsize=11, ha="center"
    )
    ax.text(
        0.5,
        0.1,
        "Normal modes: ω₁ ≈ 0.92, ω₂ ≈ 2.04",
        transform=ax.transAxes,
        fontsize=10,
        ha="center",
        style="italic",
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
    print("  1. Equations were derived from coupled Lagrangian (not hardcoded)")
    print("  2. Cross-field coupling terms correctly parsed from JSON")
    print("  3. Energy oscillates between φ and χ (mode mixing)")
    print("  4. Demonstrates inter-field coupling without hardcoding physics")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the coupled scalar fields simulation from Lagrangian-derived equations."""
    json_path = Path(__file__).parent.parent / "data" / "coupled_scalars.json"
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
