"""De Sitter Klein-Gordon simulation - Hubble friction in expanding spacetime.

This script demonstrates a scalar field on de Sitter (exponentially expanding) spacetime.

Physics:
  Metric: ds^2 = Omega(t)^2 * (-dt^2 + dx^2)  where Omega(t) = exp(H*t)

  For de Sitter expansion:
  - Christoffel symbols are non-zero (computed from metric, not hardcoded)
  - Hubble friction term: -n*H * d_t phi where n = number of spatial dimensions
  - For 1+1D (n=1): -H * d_t phi
  - Effective mass: m_eff^2 = m^2 * Omega^2 = m^2 * exp(2Ht) (time-dependent)

  Wave equation (1+1D): d2_t phi = d2_x phi - H * d_t phi - m^2 * phi

  Key observation: Energy is NOT conserved - waves lose energy to expansion.
  Amplitude decays approximately as exp(-H*t/2) for 1+1D due to Hubble friction.
  Energy decays as exp(-H*t) for 1+1D.

Verification:
  Compare with flat space KG - de Sitter should show:
  1. Amplitude decay (Hubble damping)
  2. Slower group velocity (expansion stretches wavelengths)
"""

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "de_sitter_kg_output.png"
AMPLITUDE_DECAY_THRESHOLD = 0.9  # Threshold to detect Hubble friction effect

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0914, PLR0915
    """Run the de Sitter Klein-Gordon simulation."""
    print("=" * 60)
    print("De Sitter Klein-Gordon Simulation (Hubble Friction)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = exp(2Ht) * (-dt^2 + dx^2)")
    print("Wave equation (1+1D): d2_t phi = d2_x phi - H * d_t phi - m^2 * phi")
    print("                                            ^^^ n=1 spatial dimension")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "de_sitter_kg.json"
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

    # Build the PDE with runtime parameters
    # These override the symbolic coefficients from JSON
    print()
    print("Step 2: Building PDE from specification...")
    hubble_param = 0.1  # Hubble constant (expansion rate)
    mass_squared = 1.0  # Mass parameter

    pde = build_pde_from_json(
        json_path, parameters={"H": hubble_param, "m2": mass_squared}
    )
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  Hubble parameter: H = {hubble_param}")
    print(f"  Mass parameter: m^2 = {mass_squared}")

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
    phi = ScalarField(grid, data=0.0, label="dSphi_0")
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
    print()
    print("Step 5: Running simulation...")
    t_end = 30.0
    dt = 0.01

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
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

    # Calculate amplitude decay ratio
    # For 1+1D (n=1), amplitude decays as exp(-H*t/2) = exp(-0.5*H*t)
    decay_ratio = final_max / initial_max if initial_max > 0 else 0
    expected_decay = np.exp(-0.5 * hubble_param * t_end)  # 1+1D Hubble damping

    print(f"  Amplitude ratio (final/initial): {decay_ratio:.4f}")
    print(f"  Expected from Hubble damping exp(-H*t/2) [1+1D]: {expected_decay:.4f}")

    if decay_ratio < AMPLITUDE_DECAY_THRESHOLD:
        print("  Hubble friction observed: amplitude has decayed")
    else:
        print("  Note: For small H or short time, damping may be subtle")

    # Energy-like quantity (will decrease due to Hubble friction)
    def compute_energy(snapshot: FieldCollection) -> float:
        phi_data = snapshot[0].data
        pi_data = snapshot[1].data
        # Approximate energy: kinetic + mass term (gradient term omitted for simplicity)
        return float(np.sum(pi_data**2 + mass_squared * phi_data**2))

    initial_energy = compute_energy(initial)
    final_energy = compute_energy(final)

    print(f"  Initial 'energy': {initial_energy:.2f}")
    print(f"  Final 'energy': {final_energy:.2f}")
    print(f"  Energy loss: {100 * (1 - final_energy / initial_energy):.1f}%")

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
    ax.set_title(r"Field $\phi$ evolution (note amplitude decay)")
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

    # Energy decay
    ax = axes[1, 0]
    energies = []
    time_values = []
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        energies.append(compute_energy(snapshot))
        time_values.append(t_idx)

    ax.plot(time_values, energies, "b-", linewidth=2, label="Simulation")

    # Expected exponential decay from Hubble friction
    # For 1+1D (n=1), energy decays as exp(-H*t)
    # Note: time_values and energies types come from MemoryStorage iteration
    t_array = np.array(time_values)  # pyright: ignore[reportUnknownArgumentType, reportUnknownVariableType]
    expected_energies = energies[0] * np.exp(-hubble_param * t_array)  # pyright: ignore[reportUnknownVariableType]
    ax.plot(
        time_values, expected_energies, "r--", linewidth=1.5, label=r"$\sim e^{-Ht}$ [1+1D]"
    )

    ax.set_xlabel("Snapshot index")
    ax.set_ylabel("Total 'energy'")
    ax.set_title(f"Energy decay from Hubble friction (H={hubble_param})")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Space-time diagram
    ax = axes[1, 1]
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
        rf"De Sitter KG: $g_{{\mu\nu}} = e^{{2Ht}} \eta_{{\mu\nu}}$, H={hubble_param}, m$^2$={mass_squared}"
        "\nHubble friction causes energy loss to expansion",
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
    print(f"  1. De Sitter metric: g_ab = exp(2Ht) * eta_ab with H={hubble_param}")
    print("  2. Hubble friction term: -n*H * d_t phi where n = spatial dimensions")
    print("     For this 1+1D simulation (n=1): -H * d_t phi")
    print("  3. Christoffel symbols computed from metric definition (not hardcoded)")
    print("  4. Energy is NOT conserved - lost to cosmic expansion")
    print(f"  5. Amplitude decayed by factor ~{decay_ratio:.3f} over t={t_end}")
    print(f"  6. Expected decay from Hubble damping exp(-H*t/2): ~{expected_decay:.3f}")
    print()
    print("  Physical interpretation:")
    print("    Waves in expanding spacetime lose energy to expansion.")
    print("    The Hubble friction term acts like a velocity-dependent damping.")
    print("    This is why the CMB photons have redshifted since the Big Bang.")
    print("=" * 60)


if __name__ == "__main__":
    main()
