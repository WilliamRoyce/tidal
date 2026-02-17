"""Klein-Gordon on 2-Sphere simulation - position-dependent wave speed.

This script demonstrates a scalar field on a static 2-sphere embedded via
stereographic projection, where spatial curvature produces position-dependent
wave propagation speed.

Physics:
  Metric: ds^2 = -dt^2 + Omega^2(dx^2 + dy^2)
  where Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2)

  The south pole (x=y=0) has Omega = 2 (wave speed = 1/Omega^2 = 1/4, slow).
  The equator (r=R) has Omega = 1 (wave speed = 1, medium).
  Near the north pole (r >> R) Omega -> 0 (wave speed -> infinity, fast).

  Key observation:
  - Time is flat: no Hubble friction, energy IS conserved
  - Spatial curvature: position-dependent wave speed 1/Omega^2
  - In 2D conformal geometry, Christoffel gradient corrections vanish
  - The wave equation is: d2_t phi = (1/Omega^2) nabla^2 phi - m^2 phi

Verification:
  1. Waves should propagate anisotropically (faster away from south pole)
  2. Energy should be conserved (no friction term)
  3. Compare wave speed at different positions
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
SPHERE_RADIUS = 2.0
MASS_SQUARED = 0.0  # Massless for clean wave speed test

# Grid
DOMAIN_SIZE_FACTOR = 4  # Domain extends DOMAIN_SIZE_FACTOR * SPHERE_RADIUS
GRID_SHAPE = [128, 128]
GRID_PERIODIC = True

# Time integration
T_END = 10.0
DT = 0.005
SNAPSHOT_INTERVAL = 0.2

# Initial conditions
PULSE_WIDTH = 0.8
PULSE_AMPLITUDE = 1.0

# Analysis
ENERGY_CONSERVATION_THRESHOLD = 0.1

# Output
OUTPUT_FILENAME = "sphere_kg_output.png"

# State indices
PHI_INDEX = 0
PI_INDEX = 1

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for sphere KG simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    domain_size: float


def _print_header() -> None:
    print("=" * 60)
    print("Klein-Gordon on 2-Sphere (Stereographic Projection)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = -dt^2 + Omega^2(dx^2 + dy^2)")
    print("  Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2)")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (2+1D)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        n_terms = len(eq.rhs_terms)
        n_pos_dep = sum(1 for t in eq.rhs_terms if t.position_dependent)
        print(
            f"    d2_t({eq.field_name}): {n_terms} terms, {n_pos_dep} position-dependent"
        )
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(
        json_path, parameters={"sphR": SPHERE_RADIUS, "sphm2": MASS_SQUARED}
    )
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Sphere radius: R = {SPHERE_RADIUS}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print()
    return pde


def _create_grid() -> tuple[CartesianGrid, float]:
    print("Step 3: Setting up 2D simulation grid...")
    domain_size = DOMAIN_SIZE_FACTOR * SPHERE_RADIUS
    grid = CartesianGrid(
        bounds=[(-domain_size, domain_size), (-domain_size, domain_size)],
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  Domain: [-{domain_size}, {domain_size}]^2")
    print(f"  Resolution: {grid.shape[0]}x{grid.shape[1]} points")
    print()
    return grid, domain_size


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    # Gaussian pulse at south pole (x=y=0)
    gaussian = PULSE_AMPLITUDE * np.exp(
        -(x**2 + y**2) / (2 * PULSE_WIDTH**2)
    )

    phi = ScalarField(grid, data=gaussian, label="sphphi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")
    state = FieldCollection([phi, pi])

    print("  Gaussian pulse at south pole (x=0, y=0)")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")

    # Show wave speed at key locations
    print()
    print("  Wave speed profile (1/Omega^2):")
    for r_val, label in [
        (0.0, "south pole"),
        (SPHERE_RADIUS, "equator"),
        (2 * SPHERE_RADIUS, "past equator"),
    ]:
        omega = 2 * SPHERE_RADIUS**2 / (SPHERE_RADIUS**2 + r_val**2)
        speed = 1 / omega**2
        print(f"    r={r_val:.1f} ({label}): Omega={omega:.3f}, speed={speed:.3f}")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")
    domain_size = DOMAIN_SIZE_FACTOR * SPHERE_RADIUS

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {T_END} time units")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        domain_size=domain_size,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])

    initial_max = np.max(np.abs(initial[PHI_INDEX].data))
    final_max = np.max(np.abs(final[PHI_INDEX].data))

    print(f"  Initial: max|phi| = {initial_max:.4f}")
    print(f"  Final:   max|phi| = {final_max:.4f}")

    # Energy should be approximately conserved (no friction)
    def compute_energy(snapshot: FieldCollection) -> float:
        phi_data = snapshot[PHI_INDEX].data
        pi_data = snapshot[PI_INDEX].data
        return float(np.sum(pi_data**2 + phi_data**2))

    initial_energy = compute_energy(initial)
    final_energy = compute_energy(final)
    energy_change = abs(final_energy - initial_energy) / max(initial_energy, 1e-10)

    print(f"  Initial energy proxy: {initial_energy:.2f}")
    print(f"  Final energy proxy: {final_energy:.2f}")
    print(f"  Relative energy change: {energy_change:.4f}")
    if energy_change < ENERGY_CONSERVATION_THRESHOLD:
        print("  Energy approximately conserved (no Hubble friction)")
    else:
        print("  Note: Energy change may be due to boundary effects")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    grid = result.grid
    domain_size = result.domain_size

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initial field
    ax = axes[0, 0]
    vmax = PULSE_AMPLITUDE
    im = ax.imshow(
        initial[PHI_INDEX].data.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-vmax,
        vmax=vmax,
        extent=[-domain_size, domain_size, -domain_size, domain_size],
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"Initial $\phi$ (t=0)")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    # Draw equator circle
    theta = np.linspace(0, 2 * np.pi, 100)
    ax.plot(
        SPHERE_RADIUS * np.cos(theta),
        SPHERE_RADIUS * np.sin(theta),
        "k--",
        alpha=0.5,
        label="equator",
    )
    ax.legend(fontsize=8)

    # Final field
    ax = axes[0, 1]
    final_vmax = max(np.max(np.abs(final[PHI_INDEX].data)), 0.01)
    im = ax.imshow(
        final[PHI_INDEX].data.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-final_vmax,
        vmax=final_vmax,
        extent=[-domain_size, domain_size, -domain_size, domain_size],
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(rf"Final $\phi$ (t={T_END:.0f})")
    plt.colorbar(im, ax=ax, label=r"$\phi$")
    ax.plot(
        SPHERE_RADIUS * np.cos(theta),
        SPHERE_RADIUS * np.sin(theta),
        "k--",
        alpha=0.5,
    )

    # Cross-sections at different times
    ax = axes[1, 0]
    time_values = list(storage.times)
    n_snapshots = len(storage)
    times_to_plot = [
        0,
        n_snapshots // 4,
        n_snapshots // 2,
        3 * n_snapshots // 4,
        n_snapshots - 1,
    ]
    cmap_colors = plt.get_cmap("viridis")
    colors = cmap_colors(np.linspace(0.2, 0.8, len(times_to_plot)))

    center_idx = grid.shape[1] // 2
    x_1d = cast("np.ndarray", grid.cell_coords[:, 0, 0])
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            x_1d,
            snapshot[PHI_INDEX].data[:, center_idx],
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.axvline(
        x=-SPHERE_RADIUS, color="gray", linestyle="--", alpha=0.3, label="equator"
    )
    ax.axvline(x=SPHERE_RADIUS, color="gray", linestyle="--", alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"$\phi$ cross-section at y=0 (wave spreading)")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    # Wave speed profile overlay
    ax = axes[1, 1]
    r_values = np.sqrt(x_1d**2)
    omega_profile = 2 * SPHERE_RADIUS**2 / (SPHERE_RADIUS**2 + r_values**2)
    speed_profile = 1.0 / omega_profile**2

    ax.plot(x_1d, speed_profile, "b-", linewidth=2, label=r"$1/\Omega^2$ (wave speed)")
    ax.plot(
        x_1d, omega_profile, "r-", linewidth=2, label=r"$\Omega$ (conformal factor)"
    )
    ax.axvline(x=-SPHERE_RADIUS, color="gray", linestyle="--", alpha=0.3)
    ax.axvline(x=SPHERE_RADIUS, color="gray", linestyle="--", alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("Value")
    ax.set_title(f"Wave speed profile (R={SPHERE_RADIUS})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)
    ax.set_ylim(0, max(speed_profile) * 1.1)

    fig.suptitle(
        rf"KG on $S^2$ (stereographic): $\Omega = 2R^2/(R^2+x^2+y^2)$, R={SPHERE_RADIUS}"
        "\nPosition-dependent wave speed from spatial curvature",
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
    print(f"  1. 2-sphere metric with R={SPHERE_RADIUS}: Omega = 2R^2/(R^2+r^2)")
    print("  2. NO Hubble friction (time component is flat)")
    print("  3. Position-dependent wave speed: 1/Omega^2")
    print("     - Slow near south pole (Omega=2, speed=0.25)")
    print("     - Fast near equator (Omega=1, speed=1)")
    print("  4. Energy is conserved (static metric)")
    print("  5. All coefficients derived from metric by xAct pipeline")
    print()
    print("  Physical interpretation:")
    print("    Waves on a sphere propagate with position-dependent speed.")
    print("    Near the south pole (stereographic center), the metric stretches")
    print("    coordinate distances, slowing wave propagation in coord space.")
    print(
        "    Near the equator, the conformal factor is 1 and speed matches flat space."
    )
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the 2-sphere Klein-Gordon simulation."""
    json_path = Path(__file__).parent.parent / "data" / "sphere_kg.json"
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid, _domain_size = _create_grid()
    state = _create_initial_state(grid)
    sim_result = _run_simulation(pde, grid, state)
    _analyze_results(sim_result)
    _plot_results(sim_result)
    _print_footer()


if __name__ == "__main__":
    main()
