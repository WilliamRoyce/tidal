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

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "proca_output.png"

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0915, PLR0914
    """Run the Proca field simulation."""
    print("=" * 60)
    print("Proca Field (Massive Vector) 1+1D Simulation")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/4 F_ab F^ab - 1/2 m^2 A_a A^a")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "proca_1d.json"
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

    # Build the PDE with runtime-specified mass parameter
    # The mass parameter (procaMassSquared) is stored symbolically in JSON
    # and can be set at simulation time
    print()
    print("Step 2: Building PDE from specification...")
    mass_squared = 1.0  # Can be changed at runtime without re-deriving equations
    pde = build_pde_from_json(json_path, parameters={"procaMassSquared": mass_squared})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Mass parameter: m² = {mass_squared} (runtime-specified)")

    # Set up 1D spatial grid
    print()
    print("Step 3: Setting up 1D simulation grid...")
    grid = CartesianGrid(
        bounds=[(0, 100)],  # x domain
        shape=[512],  # 512 grid points
        periodic=True,
    )
    print("  Domain: [0, 100]")
    print(f"  Resolution: {grid.shape[0]} points")
    print("  Periodic boundary conditions")

    # Create initial conditions
    print()
    print("Step 4: Creating initial conditions...")

    # Initial state: 4 fields = 2 components * (field + momentum)
    # State layout: procaA_0, pi_0, procaA_1, pi_1

    # Start with zero fields
    a0 = ScalarField(grid, data=0.0, label="procaA_0")
    pi0 = ScalarField(grid, data=0.0, label="pi_0")
    a1 = ScalarField(grid, data=0.0, label="procaA_1")
    pi1 = ScalarField(grid, data=0.0, label="pi_1")

    # Initialize a Gaussian pulse in A_1 (spatial component)
    # This demonstrates dispersive propagation for massive vector
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    center_x = 50.0
    width = 5.0
    amplitude = 1.0

    gaussian = amplitude * np.exp(-((x - center_x) ** 2) / (2 * width**2))
    a1.data[:] = gaussian

    state = FieldCollection([a0, pi0, a1, pi1])
    print(f"  Gaussian pulse in A_1 at center x={center_x}")
    print(f"  Width: {width}, Amplitude: {amplitude}")
    print("  A_0 and momenta initialized to zero")

    # Run simulation
    # Note: Using smaller dt for better energy conservation.
    # The explicit Euler integrator does not exactly conserve the Hamiltonian;
    # energy drift decreases with smaller dt (dt=0.01 gives ~36% drift, dt=0.005 gives ~16%).
    print()
    print("Step 5: Running simulation...")
    t_end = 30.0
    dt = 0.005  # Smaller timestep for better energy conservation

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        scheme="runge-kutta",  # RK4 for better energy conservation
        tracker=storage.tracker(1.0),  # Store every 1.0 time units
    )
    result = normalize_solve_result(result)  # Ensure consistent output format

    print(f"  Duration: {t_end} time units")
    print(f"  Stored {len(storage)} snapshots")

    # Analyze results
    print()
    print("Step 6: Analyzing results...")

    # Get initial and final states
    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    # Extract field amplitudes
    initial_a0 = np.max(np.abs(initial[0].data))
    initial_a1 = np.max(np.abs(initial[2].data))

    final_a0 = np.max(np.abs(final[0].data))
    final_a1 = np.max(np.abs(final[2].data))

    print(f"  Initial: max|A_0| = {initial_a0:.4f}, max|A_1| = {initial_a1:.4f}")
    print(f"  Final:   max|A_0| = {final_a0:.4f}, max|A_1| = {final_a1:.4f}")

    # Check for dispersion (pulse spreading)
    initial_width = np.sqrt(
        np.sum((x - center_x) ** 2 * initial[2].data ** 2)
        / np.sum(initial[2].data ** 2)
    )
    # Find center of mass of final A_1
    final_com = np.sum(x * final[2].data ** 2) / (np.sum(final[2].data ** 2) + 1e-10)
    final_width = np.sqrt(
        np.sum((x - final_com) ** 2 * final[2].data ** 2)
        / (np.sum(final[2].data ** 2) + 1e-10)
    )

    print(f"  Initial pulse width: {initial_width:.2f}")
    print(f"  Final pulse width: {final_width:.2f}")

    if final_width > initial_width * 1.2:
        print("  Dispersion observed: pulse has spread (as expected for massive field)")
    else:
        print("  Note: Pulse spread may be subtle at this mass/time scale")

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Top row: A_0 and A_1 evolution
    # Get time snapshots
    times = [0, len(storage) // 3, 2 * len(storage) // 3, len(storage) - 1]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times)))

    # A_0 evolution
    ax = axes[0, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[0].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("A_0")
    ax.set_title("A_0 (temporal component) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # A_1 evolution
    ax = axes[0, 1]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[2].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("A_1")
    ax.set_title("A_1 (spatial component) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Bottom left: Momentum pi_1 evolution
    ax = axes[1, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[3].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel("π_1")
    ax.set_title("π_1 (momentum of A_1) evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Bottom right: Energy (Hamiltonian) over time
    # Note: The Hamiltonian should be conserved for the continuous Klein-Gordon equation.
    # Using RK4 scheme (scheme='runge-kutta') provides much better energy conservation
    # than the default explicit Euler, which is non-symplectic and causes linear drift.
    ax = axes[1, 1]
    energies = []
    time_values = []

    # Grid spacing for gradient calculation
    dx = float(x[1] - x[0])

    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])

        # Extract fields and momenta
        a0, pi0 = snapshot[0].data, snapshot[1].data
        a1, pi1 = snapshot[2].data, snapshot[3].data

        # Compute spatial gradients
        grad_a0 = np.gradient(a0, dx)
        grad_a1 = np.gradient(a1, dx)

        # Klein-Gordon Hamiltonian: H = ½∫[π² + (∂_x φ)² + m²φ²] dx
        energy_a0 = np.sum(0.5 * pi0**2 + 0.5 * grad_a0**2 + 0.5 * mass_squared * a0**2)
        energy_a1 = np.sum(0.5 * pi1**2 + 0.5 * grad_a1**2 + 0.5 * mass_squared * a1**2)

        energies.append(energy_a0 + energy_a1)
        time_values.append(t_idx)

    ax.plot(time_values, energies, "b-", linewidth=2)
    ax.set_xlabel("Snapshot index")
    ax.set_ylabel("Total Energy")
    ax.set_title(r"Hamiltonian $H = \frac{1}{2}\int[\pi^2 + (\nabla\phi)^2 + m^2\phi^2]dx$")
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        f"Proca Field 1+1D: Massive Vector (m²={mass_squared})\nDispersion relation: ω² = k² + m²",
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
    print("  1. Equations derived from Proca Lagrangian (Maxwell + mass term)")
    print("  2. Each component satisfies massive Klein-Gordon equation")
    print("  3. Dispersion relation: ω² = k² + m² (subluminal propagation)")
    print("  4. Unlike massless Maxwell, pulses spread due to dispersion")
    print("  5. Lorenz condition ∂_μ A^μ = 0 is automatic for m ≠ 0")
    print(f"  6. Mass parameter m²={mass_squared} specified at runtime (not pre-evaluated)")
    print("=" * 60)


if __name__ == "__main__":
    main()
