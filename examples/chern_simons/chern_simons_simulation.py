"""Chern-Simons 2+1D simulation from Lagrangian-derived equations.

This script demonstrates the Maxwell-Chern-Simons theory in 2+1D:
  L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c

The Chern-Simons term gives the photon a topological mass, leading to
helical wave propagation patterns.

Equations of motion (in Lorenz gauge):
  d2_t(A_0) = laplacian(A_0) + kappa*(d_x A_2 - d_y A_1)
  d2_t(A_1) = laplacian(A_1) + kappa*(d_y A_0)
  d2_t(A_2) = laplacian(A_2) - kappa*(d_x A_0)
"""

from pathlib import Path
from typing import cast

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0915, PLR0914
    """Run the Chern-Simons simulation."""
    print("=" * 60)
    print("Chern-Simons 2+1D Simulation from Lagrangian")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "chern_simons_3d.json"
    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Dimension: {spec.dimension} (2+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    # Check for cross-field coupling (CS terms)
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field != eq.field_name:
                print(
                    f"  CS coupling: {eq.field_name} <- {term.field} ({term.operator})"
                )

    # Build the PDE
    print()
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")

    # Set up 2D spatial grid (x, y)
    print()
    print("Step 3: Setting up 2D simulation grid...")
    grid = CartesianGrid(
        bounds=[(0, 50), (0, 50)],  # x and y domain
        shape=[64, 64],  # 64x64 grid
        periodic=True,
    )
    print(f"  Domain: [{grid.axes_bounds[0]}, {grid.axes_bounds[1]}]")
    print(f"  Resolution: {grid.shape}")
    print("  Periodic boundary conditions")

    # Create initial conditions
    print()
    print("Step 4: Creating initial conditions...")

    # Initial state: 6 fields = 3 components * (field + momentum)
    # State layout: A_0, pi_0, A_1, pi_1, A_2, pi_2

    # Start with zero fields
    a0 = ScalarField(grid, data=0.0, label="A_0")
    pi0 = ScalarField(grid, data=0.0, label="pi_0")
    a1 = ScalarField(grid, data=0.0, label="A_1")
    pi1 = ScalarField(grid, data=0.0, label="pi_1")
    a2 = ScalarField(grid, data=0.0, label="A_2")
    pi2 = ScalarField(grid, data=0.0, label="pi_2")

    # Initialize a Gaussian pulse in A_1 (x-component of vector potential)
    # This will excite all components via the CS coupling
    x, y = (
        cast("np.ndarray", grid.cell_coords[..., 0]),
        cast("np.ndarray", grid.cell_coords[..., 1]),
    )
    center_x, center_y = 25.0, 25.0
    width = 5.0
    amplitude = 1.0

    gaussian = amplitude * np.exp(
        -((x - center_x) ** 2 + (y - center_y) ** 2) / (2 * width**2)
    )
    a1.data[:] = gaussian

    state = FieldCollection([a0, pi0, a1, pi1, a2, pi2])
    print(f"  Gaussian pulse in A_1 at center ({center_x}, {center_y})")
    print(f"  Width: {width}, Amplitude: {amplitude}")

    # Run simulation
    print()
    print("Step 5: Running simulation...")
    t_end = 10.0
    dt = 0.01

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        tracker=storage.tracker(0.5),  # Store every 0.5 time units
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
    initial_a2 = np.max(np.abs(initial[4].data))

    final_a0 = np.max(np.abs(final[0].data))
    final_a1 = np.max(np.abs(final[2].data))
    final_a2 = np.max(np.abs(final[4].data))

    print(
        f"  Initial: max|A_0| = {initial_a0:.3f}, max|A_1| = {initial_a1:.3f}, max|A_2| = {initial_a2:.3f}"
    )
    print(
        f"  Final:   max|A_0| = {final_a0:.3f}, max|A_1| = {final_a1:.3f}, max|A_2| = {final_a2:.3f}"
    )

    # Check for energy transfer (CS coupling effect)
    if final_a0 > 0.01 or final_a2 > 0.01:  # noqa: PLR2004
        print("  CS coupling effect observed: energy transferred to A_0 and A_2")
    else:
        print("  Note: CS coupling may be weak at this kappa value")

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Initial state
    vmax = max(initial_a1, 0.1)
    im0 = axes[0, 0].imshow(
        initial[0].data.T, origin="lower", cmap="RdBu", vmin=-vmax, vmax=vmax
    )
    axes[0, 0].set_title("Initial A_0")
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(
        initial[2].data.T, origin="lower", cmap="RdBu", vmin=-vmax, vmax=vmax
    )
    axes[0, 1].set_title("Initial A_1")
    plt.colorbar(im1, ax=axes[0, 1])

    im2 = axes[0, 2].imshow(
        initial[4].data.T, origin="lower", cmap="RdBu", vmin=-vmax, vmax=vmax
    )
    axes[0, 2].set_title("Initial A_2")
    plt.colorbar(im2, ax=axes[0, 2])

    # Final state
    vmax_f = max(final_a0, final_a1, final_a2, 0.1)
    im3 = axes[1, 0].imshow(
        final[0].data.T, origin="lower", cmap="RdBu", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 0].set_title(f"Final A_0 (t={t_end})")
    plt.colorbar(im3, ax=axes[1, 0])

    im4 = axes[1, 1].imshow(
        final[2].data.T, origin="lower", cmap="RdBu", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 1].set_title(f"Final A_1 (t={t_end})")
    plt.colorbar(im4, ax=axes[1, 1])

    im5 = axes[1, 2].imshow(
        final[4].data.T, origin="lower", cmap="RdBu", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 2].set_title(f"Final A_2 (t={t_end})")
    plt.colorbar(im5, ax=axes[1, 2])

    fig.suptitle("Chern-Simons 2+1D: Maxwell-CS with kappa=0.5")
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True)
    output_path = output_dir / "chern_simons_output.png"
    plt.savefig(output_path, dpi=150)
    print(f"  Saved plot to: {output_path}")

    # Summary
    print()
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations derived from CS Lagrangian (Maxwell + topological term)")
    print("  2. Cross-field coupling via gradient operators")
    print("  3. Energy can transfer between A components")
    print("  4. Chern-Simons term adds effective mass to photon")
    print("=" * 60)


if __name__ == "__main__":
    main()
