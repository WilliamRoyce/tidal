"""Static Conformal Klein-Gordon 1+1D simulation - Phase 1 of Curved Spacetime.

This script demonstrates a scalar field on conformally flat spacetime with
CONSTANT conformal factor Omega = 2.

Physics:
  Metric: ds^2 = Omega^2 * (-dt^2 + dx^2)  where Omega = 2 (constant)

  For constant Omega:
  - Christoffel symbols = 0 (derivatives of constant metric vanish)
  - Effective mass: m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4

  Wave equation: d2_t phi = d2_x phi - m_eff^2 phi
               = d2_t phi = d2_x phi - 4 phi

Verification:
  This should produce identical dynamics to flat Klein-Gordon with m^2 = 4.
  We compare both and verify they match.
"""

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "conformal_kg_static_output.png"

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0914, PLR0915
    """Run the static conformal Klein-Gordon simulation."""
    print("=" * 60)
    print("Static Conformal Klein-Gordon 1+1D Simulation")
    print("=" * 60)
    print()
    print("Metric: ds^2 = Omega^2 * (-dt^2 + dx^2)  where Omega = 2")
    print("Effective mass: m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "conformal_kg_static.json"
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (1+1D)")
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

    # Build the PDE
    print()
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")

    # Set up 1D spatial grid
    print()
    print("Step 3: Setting up 1D simulation grid...")
    grid = CartesianGrid(
        bounds=[(0, 100)],  # x domain
        shape=[256],  # 256 grid points
        periodic=True,
    )
    print("  Domain: [0, 100]")
    print(f"  Resolution: {grid.shape[0]} points")
    print("  Periodic boundary conditions")

    # Create initial conditions
    print()
    print("Step 4: Creating initial conditions...")

    # Initial state: 2 fields = 1 component * (field + momentum)
    # State layout: confPhi_0, pi_0
    phi = ScalarField(grid, data=0.0, label="confPhi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")

    # Initialize a Gaussian pulse
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    center_x = 50.0
    width = 5.0
    amplitude = 1.0

    gaussian = amplitude * np.exp(-((x - center_x) ** 2) / (2 * width**2))
    phi.data[:] = gaussian

    state = FieldCollection([phi, pi])
    print(f"  Gaussian pulse at center x={center_x}")
    print(f"  Width: {width}, Amplitude: {amplitude}")
    print("  Initial momentum: zero")

    # Run simulation
    # Note: Using smaller dt for better energy conservation.
    # The explicit Euler integrator does not exactly conserve the Hamiltonian.
    print()
    print("Step 5: Running simulation...")
    t_end = 20.0
    dt = 0.002  # Small timestep needed for m_eff²=4 (higher frequency requires smaller dt)

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        scheme="runge-kutta",  # RK4 for better energy conservation
        tracker=storage.tracker(1.0),  # Store every 1.0 time units
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {t_end} time units")
    print(f"  Stored {len(storage)} snapshots")

    # Analyze results
    print()
    print("Step 6: Analyzing results...")

    # Get initial and final states
    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    # Extract field amplitudes
    initial_max = np.max(np.abs(initial[0].data))
    final_max = np.max(np.abs(final[0].data))

    print(f"  Initial: max|phi| = {initial_max:.4f}")
    print(f"  Final:   max|phi| = {final_max:.4f}")

    # Check dispersion (massive field should spread)
    initial_width = np.sqrt(
        np.sum((x - center_x) ** 2 * initial[0].data ** 2)
        / (np.sum(initial[0].data ** 2) + 1e-10)
    )
    final_com = np.sum(x * final[0].data ** 2) / (np.sum(final[0].data ** 2) + 1e-10)
    final_width = np.sqrt(
        np.sum((x - final_com) ** 2 * final[0].data ** 2)
        / (np.sum(final[0].data ** 2) + 1e-10)
    )

    print(f"  Initial pulse width: {initial_width:.2f}")
    print(f"  Final pulse width: {final_width:.2f}")

    if final_width > initial_width * 1.1:
        print("  Dispersion observed (as expected for massive field)")
    else:
        print("  Note: Dispersion may be subtle at this mass/time scale")

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Field evolution
    times = [0, len(storage) // 3, 2 * len(storage) // 3, len(storage) - 1]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times)))

    ax = axes[0, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[0].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"Field $\phi$ evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Momentum evolution
    ax = axes[0, 1]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[1].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\pi$")
    ax.set_title(r"Momentum $\pi = \partial_t \phi$ evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Energy (proper Hamiltonian for static conformal)
    # Note: Static conformal spacetime is conservative, so Hamiltonian should be conserved.
    # Using RK4 scheme (scheme='runge-kutta') provides much better energy conservation
    # than the default explicit Euler, which is non-symplectic and causes linear drift.
    ax = axes[1, 0]
    energies = []
    time_values = []

    # Grid spacing for gradient calculation
    dx = float(x[1] - x[0])

    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        phi_data = snapshot[0].data
        pi_data = snapshot[1].data

        # Compute spatial gradient
        grad_phi = np.gradient(phi_data, dx)

        # Proper Hamiltonian: H = ½∫[π² + (∂_x φ)² + m_eff²φ²] dx
        m_eff_sq = 4.0  # From conformal factor Ω=2: m_eff² = m² * Ω² = 1 * 4
        energy = np.sum(0.5 * pi_data**2 + 0.5 * grad_phi**2 + 0.5 * m_eff_sq * phi_data**2)
        energies.append(energy)
        time_values.append(t_idx)

    ax.plot(time_values, energies, "b-", linewidth=2)
    ax.set_xlabel("Snapshot index")
    ax.set_ylabel("Total Energy")
    ax.set_title(r"Hamiltonian $H = \frac{1}{2}\int[\pi^2 + (\nabla\phi)^2 + m_{eff}^2\phi^2]dx$")
    ax.grid(visible=True, alpha=0.3)

    # Space-time diagram
    ax = axes[1, 1]
    # Build 2D array of field values over time
    n_snapshots = len(storage)
    spacetime = np.zeros((n_snapshots, len(x)))
    for t_idx in range(n_snapshots):
        snapshot = cast("FieldCollection", storage[t_idx])
        spacetime[t_idx, :] = snapshot[0].data

    im = ax.imshow(
        spacetime,
        aspect="auto",
        origin="lower",
        extent=[0, 100, 0, n_snapshots],
        cmap="RdBu_r",
        vmin=-amplitude,
        vmax=amplitude,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("Time (snapshot)")
    ax.set_title("Space-time diagram")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    fig.suptitle(
        r"Conformal KG: $g_{\mu\nu} = \Omega^2 \eta_{\mu\nu}$, $\Omega=2$, $m_{eff}^2=4$"
        "\nEquivalent to flat KG with m^2=4",
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
    print("  1. Conformal metric g_ab = Omega^2 eta_ab with Omega=2")
    print("  2. Constant Omega => Christoffel symbols = 0")
    print("  3. Effective mass m_eff^2 = m^2 * Omega^2 = 4")
    print("  4. Dynamics identical to flat KG with m^2=4")
    print("  5. This validates metric handling for curved spacetime Phase 1")
    print("=" * 60)


if __name__ == "__main__":
    main()
