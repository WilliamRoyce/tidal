"""Linearized massive gravity in 2+1D simulation.

This script demonstrates massive gravity in 2+1D:
  G^(1)_ab[h] + m^2 h_ab = 0

The mass term m^2 h_ab is obtained by linearizing G_ab + m^2 g_ab = 0
around flat Minkowski. In 2+1D, pure GR has no local degrees of freedom;
the mass term creates a propagating massive mode with:
  - Dispersion relation: omega^2 = k^2 + m^2
  - Group velocity: v_g = k/omega < 1 (subluminal)

Component structure (symmetric rank-2 in 3D):
  h_0 (h_tt): constraint (no d2_t)
  h_1, h_2 (h_tx, h_ty): first-order in time
  h_3, h_4, h_5 (h_xx, h_xy, h_yy): second-order evolution
"""

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "massive_gravity_output.png"

import matplotlib as mpl

from tidal.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage, ScalarField

from tidal.symbolic import build_pde_from_json, load_equation_system


def main() -> None:  # noqa: PLR0915, PLR0914
    """Run the massive gravity simulation."""
    print("=" * 60)
    print("Linearized Massive Gravity 2+1D Simulation")
    print("=" * 60)
    print()
    print("Equation: G^(1)_ab[h] + m^2 h_ab = 0")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "massive_gravity_3d.json"
    mass_squared = 1.0
    params = {"m2": mass_squared}
    spec = load_equation_system(json_path, parameters=params)

    print(f"  Dimension: {spec.dimension} (2+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Mass parameter: m^2 = {mass_squared}")

    # Show equation structure
    print()
    print("  Equation structure:")
    for eq in spec.equations:
        lhs = eq.lhs
        terms_str = " + ".join(
            f"{term.coefficient:.2g}*{term.operator}({term.field})"
            for term in eq.rhs_terms
        )
        print(f"    {lhs.expression} = {terms_str}")

    # Show mass matrix (should be non-zero)
    print()
    print("  Mass matrix (should have m^2 on diagonal):")
    for row in spec.coupling.mass_matrix:
        print(f"    [{', '.join(f'{v:.2g}' for v in row)}]")

    # Build the PDE
    print()
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters=params)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")

    # Set up 2D spatial grid
    print()
    print("Step 3: Setting up 2D simulation grid...")
    lx, ly = 50.0, 50.0
    nx, ny = 64, 64
    grid = CartesianGrid(
        bounds=[(0, lx), (0, ly)],
        shape=[nx, ny],
        periodic=True,
    )
    print(f"  Domain: [0, {lx}] x [0, {ly}]")
    print(f"  Resolution: {nx} x {ny}")
    print("  Periodic boundary conditions")

    # Create initial conditions
    print()
    print("Step 4: Creating initial conditions...")

    # Build state from spec layout
    fields: list[ScalarField] = []
    x_coords = cast("np.ndarray", grid.cell_coords[..., 0])
    y_coords = cast("np.ndarray", grid.cell_coords[..., 1])

    # Gaussian pulse in h_3 (h_xx spatial component)
    center_x, center_y = lx / 2, ly / 2
    width = 3.0
    amplitude = 0.1
    gaussian = amplitude * np.exp(
        -((x_coords - center_x) ** 2 + (y_coords - center_y) ** 2) / (2 * width**2)
    )

    # Build state according to spec layout
    for slot in spec.state_layout:
        sf = ScalarField(grid, data=0.0, label=slot)
        if slot == "h_3":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)
    print(f"  Gaussian pulse in h_3 (h_xx) at ({center_x}, {center_y})")
    print(f"  Width: {width}, Amplitude: {amplitude}")
    print(f"  State layout: {spec.state_layout}")

    # Run simulation
    print()
    print("Step 5: Running simulation...")
    t_end = 5.0
    dt = 0.005
    snapshot_interval = 0.25

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        scheme="runge-kutta",
        tracker=storage.tracker(snapshot_interval),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {t_end} time units")
    print(f"  Stored {len(storage)} snapshots")

    # Analyze results
    print()
    print("Step 6: Analyzing results...")

    # Find h_3 index in state layout
    h3_idx = spec.state_layout.index("h_3")

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    initial_max = float(np.max(np.abs(initial[h3_idx].data)))
    final_max = float(np.max(np.abs(final[h3_idx].data)))

    print(f"  Initial max|h_3|: {initial_max:.6f}")
    print(f"  Final max|h_3|:   {final_max:.6f}")

    if final_max < initial_max * 0.95:
        print("  Dispersion observed: amplitude decreased (massive field spreading)")
    else:
        print("  Note: Pulse spread may be subtle at this mass/time scale")

    # Track amplitude over time
    max_amplitudes = []
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        max_amplitudes.append(float(np.max(np.abs(snapshot[h3_idx].data))))

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Top-left: h_3 initial profile (2D heatmap)
    ax = axes[0, 0]
    im = ax.imshow(
        initial[h3_idx].data.T,
        origin="lower",
        extent=[0, lx, 0, ly],
        cmap="RdBu_r",
        vmin=-amplitude,
        vmax=amplitude,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("h_3 (h_xx) — Initial")
    plt.colorbar(im, ax=ax)

    # Top-right: h_3 final profile (2D heatmap)
    ax = axes[0, 1]
    im = ax.imshow(
        final[h3_idx].data.T,
        origin="lower",
        extent=[0, lx, 0, ly],
        cmap="RdBu_r",
        vmin=-amplitude,
        vmax=amplitude,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"h_3 (h_xx) — t = {t_end}")
    plt.colorbar(im, ax=ax)

    # Bottom-left: max|h_3| over time
    ax = axes[1, 0]
    time_values = np.linspace(0, t_end, len(max_amplitudes))
    ax.plot(time_values, max_amplitudes, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("max|h_3|")
    ax.set_title("Peak amplitude evolution")
    ax.grid(visible=True, alpha=0.3)

    # Bottom-right: cross-section comparison (y = ly/2)
    ax = axes[1, 1]
    mid_y = ny // 2
    x_1d = np.linspace(0, lx, nx)
    ax.plot(x_1d, initial[h3_idx].data[:, mid_y], "b-", label="Initial", linewidth=2)
    ax.plot(x_1d, final[h3_idx].data[:, mid_y], "r--", label=f"t = {t_end}", linewidth=2)
    ax.set_xlabel("x")
    ax.set_ylabel("h_3")
    ax.set_title(f"Cross-section at y = {ly / 2}")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        f"Massive Gravity 2+1D: m$^2$ = {mass_squared}\n"
        r"$G^{(1)}_{ab}[h] + m^2 h_{ab} = 0$   |   $\omega^2 = k^2 + m^2$",
        fontsize=12,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved plot to: {output_path}")

    # Summary
    print()
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations derived by linearizing G_ab + m^2 g_ab = 0 via xPert")
    print("  2. Mass term m^2 h_ab creates propagating DOF in 2+1D")
    print("  3. Dispersion relation: omega^2 = k^2 + m^2 (subluminal)")
    print("  4. Gaussian pulse in h_xx disperses due to mass")
    print(f"  5. Mass parameter m^2 = {mass_squared} specified at runtime")
    print("=" * 60)


if __name__ == "__main__":
    main()
