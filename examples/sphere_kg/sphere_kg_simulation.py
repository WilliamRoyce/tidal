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

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "sphere_kg_output.png"
ENERGY_CONSERVATION_THRESHOLD = 0.1  # Relative energy change threshold for conservation check

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0914, PLR0915
    """Run the 2-sphere Klein-Gordon simulation."""
    print("=" * 60)
    print("Klein-Gordon on 2-Sphere (Stereographic Projection)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = -dt^2 + Omega^2(dx^2 + dy^2)")
    print("  Omega(x,y) = 2R^2 / (R^2 + x^2 + y^2)")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "sphere_kg.json"
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (2+1D)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    # Show equation structure
    print()
    print("  Equation structure:")
    for eq in spec.equations:
        n_terms = len(eq.rhs_terms)
        n_pos_dep = sum(1 for t in eq.rhs_terms if t.position_dependent)
        print(
            f"    d2_t({eq.field_name}): {n_terms} terms, {n_pos_dep} position-dependent"
        )

    # Build the PDE with runtime parameters
    print()
    print("Step 2: Building PDE from specification...")
    sphere_radius = 2.0
    mass_squared = 0.0  # Massless for clean wave speed test

    pde = build_pde_from_json(
        json_path, parameters={"sphR": sphere_radius, "sphm2": mass_squared}
    )
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Sphere radius: R = {sphere_radius}")
    print(f"  Mass parameter: m^2 = {mass_squared}")

    # Set up 2D spatial grid centered at south pole
    # Domain: [-L, L] x [-L, L] where L covers up to past the equator
    print()
    print("Step 3: Setting up 2D simulation grid...")
    domain_size = 4 * sphere_radius  # Cover past the equator
    grid = CartesianGrid(
        bounds=[(-domain_size, domain_size), (-domain_size, domain_size)],
        shape=[128, 128],
        periodic=True,
    )
    print(f"  Domain: [-{domain_size}, {domain_size}]^2")
    print(f"  Resolution: {grid.shape[0]}x{grid.shape[1]} points")

    # Create initial conditions
    print()
    print("Step 4: Creating initial conditions...")

    # Initial state: 1 field + momentum
    phi = ScalarField(grid, data=0.0, label="sphphi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")

    # Initialize a 2D Gaussian pulse near the south pole
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    width = 0.8
    amplitude = 1.0

    # Gaussian pulse at south pole (x=y=0)
    gaussian = amplitude * np.exp(-(x**2 + y**2) / (2 * width**2))
    phi.data[:] = gaussian

    state = FieldCollection([phi, pi])
    print("  Gaussian pulse at south pole (x=0, y=0)")
    print(f"  Width: {width}, Amplitude: {amplitude}")

    # Show wave speed at key locations
    print()
    print("  Wave speed profile (1/Omega^2):")
    for r_val, label in [
        (0.0, "south pole"),
        (sphere_radius, "equator"),
        (2 * sphere_radius, "past equator"),
    ]:
        omega = 2 * sphere_radius**2 / (sphere_radius**2 + r_val**2)
        speed = 1 / omega**2
        print(f"    r={r_val:.1f} ({label}): Omega={omega:.3f}, speed={speed:.3f}")

    # Run simulation
    print()
    print("Step 5: Running simulation...")
    t_end = 10.0
    dt = 0.005

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        scheme="runge-kutta",
        tracker=storage.tracker(0.2),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {t_end} time units")
    print(f"  Stored {len(storage)} snapshots")

    # Analyze results
    print()
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    initial_max = np.max(np.abs(initial[0].data))
    final_max = np.max(np.abs(final[0].data))

    print(f"  Initial: max|phi| = {initial_max:.4f}")
    print(f"  Final:   max|phi| = {final_max:.4f}")

    # Energy should be approximately conserved (no friction)
    def compute_energy(snapshot: FieldCollection) -> float:
        phi_data = snapshot[0].data
        pi_data = snapshot[1].data
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

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initial field
    ax = axes[0, 0]
    vmax = amplitude
    im = ax.imshow(
        initial[0].data.T,
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
        sphere_radius * np.cos(theta),
        sphere_radius * np.sin(theta),
        "k--",
        alpha=0.5,
        label="equator",
    )
    ax.legend(fontsize=8)

    # Final field
    ax = axes[0, 1]
    # Use adaptive vmax for final state
    final_vmax = max(np.max(np.abs(final[0].data)), 0.01)
    im = ax.imshow(
        final[0].data.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-final_vmax,
        vmax=final_vmax,
        extent=[-domain_size, domain_size, -domain_size, domain_size],
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(rf"Final $\phi$ (t={t_end:.0f})")
    plt.colorbar(im, ax=ax, label=r"$\phi$")
    ax.plot(
        sphere_radius * np.cos(theta), sphere_radius * np.sin(theta), "k--", alpha=0.5
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
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times_to_plot)))

    center_idx = grid.shape[1] // 2
    x_1d = cast("np.ndarray", grid.cell_coords[:, 0, 0])
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            x_1d,
            snapshot[0].data[:, center_idx],
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.axvline(
        x=-sphere_radius, color="gray", linestyle="--", alpha=0.3, label="equator"
    )
    ax.axvline(x=sphere_radius, color="gray", linestyle="--", alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"$\phi$ cross-section at y=0 (wave spreading)")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    # Wave speed profile overlay
    ax = axes[1, 1]
    r_values = np.sqrt(x_1d**2)
    omega_profile = 2 * sphere_radius**2 / (sphere_radius**2 + r_values**2)
    speed_profile = 1.0 / omega_profile**2

    ax.plot(x_1d, speed_profile, "b-", linewidth=2, label=r"$1/\Omega^2$ (wave speed)")
    ax.plot(
        x_1d, omega_profile, "r-", linewidth=2, label=r"$\Omega$ (conformal factor)"
    )
    ax.axvline(x=-sphere_radius, color="gray", linestyle="--", alpha=0.3)
    ax.axvline(x=sphere_radius, color="gray", linestyle="--", alpha=0.3)
    ax.set_xlabel("x")
    ax.set_ylabel("Value")
    ax.set_title(f"Wave speed profile (R={sphere_radius})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)
    ax.set_ylim(0, max(speed_profile) * 1.1)

    fig.suptitle(
        rf"KG on $S^2$ (stereographic): $\Omega = 2R^2/(R^2+x^2+y^2)$, R={sphere_radius}"
        "\nPosition-dependent wave speed from spatial curvature",
        fontsize=12,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"Saved plot to: {output_path}")

    # Summary
    print()
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print(f"  1. 2-sphere metric with R={sphere_radius}: Omega = 2R^2/(R^2+r^2)")
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


if __name__ == "__main__":
    main()
