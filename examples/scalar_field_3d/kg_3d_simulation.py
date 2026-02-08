"""3+1D Klein-Gordon simulation from Lagrangian-derived equations.

This script demonstrates Klein-Gordon wave propagation in full 3+1D spacetime
using the equation derived from the Lagrangian:

    L = -1/2 eta^{ab} d_a phi d_b phi - 1/2 m^2 phi^2

In flat Minkowski spacetime with signature (-,+,+,+), the equation of motion is:

    d^2_t(phi) = laplacian_x(phi) + laplacian_y(phi) + laplacian_z(phi) - m^2 * phi

The simulation starts with a 3D Gaussian pulse at rest, which radiates outward
spherically. For a massless field (m=0), amplitude decays as 1/r (the 3D Green's
function). For a massive field, dispersion causes faster decay.

Performance note:
    3D grids scale as N^3. A 32^3 grid has ~33K cells; 64^3 has ~262K cells
    (8x more). This example uses 32^3 for a practical demo (~60s runtime).
    Increase resolution for production runs at the cost of longer computation.
"""

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "kg_3d_output.png"

import matplotlib as mpl

from torsion_gertsenshtein.utils import normalize_solve_result

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage

from torsion_gertsenshtein.symbolic import (
    build_pde_from_json,
    create_initial_state,
    load_equation_system,
)


def main() -> None:  # noqa: PLR0914, PLR0915
    """Run the 3+1D Klein-Gordon simulation.

    Raises
    ------
    FileNotFoundError
        If the JSON specification file is not found.
    """
    print("=" * 60)
    print("Klein-Gordon 3+1D (Massive Scalar Field)")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/2 (nabla phi)^2 - 1/2 m^2 phi^2")
    print("Spacetime: 3+1D flat Minkowski, signature (-,+,+,+)")
    print()

    # ----------------------------------------------------------------
    # Step 1: Load the equation specification
    # ----------------------------------------------------------------
    print("Step 1: Loading equation specification...")
    json_path = (
        Path(__file__).parent.parent / "data" / "klein_gordon_3d.json"
    )

    if not json_path.exists():
        msg = (
            f"JSON specification not found: {json_path}\n"
            "Run the Wolfram derivation first: "
            "wolframscript -file examples/scalar_field_3d/klein_gordon_3d.wls"
        )
        raise FileNotFoundError(msg)

    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Spacetime dimension: {spec.dimension} (3+1D)")
    print(f"  Spatial dimension: {spec.spatial_dimension}")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Coordinates: {spec.effective_coordinates}")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(f"    d2_t({eq.field_name}) = {' + '.join(operators)}")

    print()
    print("  State layout:")
    for i, (name, slot_type) in enumerate(spec.state_layout):
        print(f"    slot {i}: {name} ({slot_type})")

    # ----------------------------------------------------------------
    # Step 2: Build the PDE
    # ----------------------------------------------------------------
    print()
    print("Step 2: Building PDE from specification...")

    m2_value = 1.0
    pde = build_pde_from_json(json_path, parameters={"m2": m2_value})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Mass parameter: m^2 = {m2_value}")
    print(f"  State size: {spec.state_size} fields (field + momentum)")

    # ----------------------------------------------------------------
    # Step 3: Create the 3D spatial grid
    # ----------------------------------------------------------------
    print()
    print("Step 3: Setting up 3D simulation grid...")

    n = 32  # Grid points per axis (32^3 = 32,768 cells)
    domain_size = 20.0

    grid = CartesianGrid(
        bounds=[(0, domain_size), (0, domain_size), (0, domain_size)],
        shape=[n, n, n],
        periodic=True,
    )
    print(f"  Domain: [0, {domain_size}]^3")
    print(f"  Resolution: {n} x {n} x {n} = {n**3:,} cells")
    print(f"  Grid spacing: dx = dy = dz = {domain_size / n:.3f}")
    print("  Periodic boundary conditions")

    # ----------------------------------------------------------------
    # Step 4: Set up initial conditions
    # ----------------------------------------------------------------
    print()
    print("Step 4: Creating initial conditions...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    z = cast("np.ndarray", grid.cell_coords[..., 2])

    center = domain_size / 2.0
    width = 2.0
    amplitude = 1.0

    # 3D Gaussian pulse centered in the domain, at rest (momentum = 0)
    phi_data = amplitude * np.exp(
        -((x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2)
        / (2 * width**2)
    )

    state = create_initial_state(
        grid, spec, field_data={"phi_0": phi_data}
    )

    initial_state = state.copy()
    initial_peak = float(np.max(np.abs(phi_data)))

    print(f"  Gaussian pulse: center=({center}, {center}, {center}), width={width}")
    print(f"  Amplitude: {amplitude}, initial peak = {initial_peak:.4f}")
    print("  Momentum: zero (pulse starts at rest)")

    # ----------------------------------------------------------------
    # Step 5: Run simulation
    # ----------------------------------------------------------------
    print()
    print("Step 5: Running simulation...")

    t_end = 8.0
    dt = 0.05
    snapshot_interval = 0.5

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
    print(f"  Time step: dt = {dt}")
    print("  Scheme: Runge-Kutta (RK4)")
    print(f"  Stored {len(storage)} snapshots")

    # ----------------------------------------------------------------
    # Step 6: Analyze results
    # ----------------------------------------------------------------
    print()
    print("Step 6: Analyzing results...")

    final = cast("FieldCollection", storage[-1])

    # Field is slot 0, momentum is slot 1
    field_slot = 0

    final_peak = float(np.max(np.abs(final[field_slot].data)))
    print(f"  Initial peak |phi|: {initial_peak:.4f}")
    print(f"  Final peak |phi|:   {final_peak:.4f}")
    print(f"  Decay ratio: {final_peak / initial_peak:.4f}")
    print("  (Expected: amplitude decreases as pulse spreads in 3D)")

    # ----------------------------------------------------------------
    # Step 7: Collect evolution data for plotting
    # ----------------------------------------------------------------
    print()
    print("Step 7: Generating visualization...")

    # Extract 1D z-profile at x=center, y=center
    ix_center = n // 2
    iy_center = n // 2
    iz_center = n // 2

    z_1d = cast("np.ndarray", grid.cell_coords[ix_center, iy_center, :, 2])

    # Collect max|phi| over time
    times = [0.0, *list(storage.times)]
    max_phi_over_time = [initial_peak]
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        max_phi_over_time.append(float(np.max(np.abs(snapshot[field_slot].data))))

    # ----------------------------------------------------------------
    # Step 8: Create plots
    # ----------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # --- Panel 1: phi z-profile (initial vs final) ---
    ax = axes[0, 0]
    initial_profile = initial_state[field_slot].data[ix_center, iy_center, :]
    final_profile = final[field_slot].data[ix_center, iy_center, :]

    ax.plot(z_1d, initial_profile, "b-", linewidth=2, label="t = 0.0 (initial)")
    ax.plot(
        z_1d, final_profile, "r-", linewidth=2, label=f"t = {times[-1]:.1f} (final)"
    )

    # Add intermediate snapshots
    n_snapshots = len(storage)
    storage_times = list(storage.times)
    for frac, color, alpha in [(0.25, "green", 0.5), (0.5, "orange", 0.5)]:
        idx = int(frac * (n_snapshots - 1))
        snap = cast("FieldCollection", storage[idx])
        profile = snap[field_slot].data[ix_center, iy_center, :]
        ax.plot(
            z_1d, profile, color=color, linewidth=1.5, alpha=alpha,
            label=f"t = {storage_times[idx]:.1f}",
        )

    ax.set_xlabel("z")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"$\phi$ along z-axis (at x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 2: 2D x-y slice (initial) ---
    ax = axes[0, 1]
    initial_xy_slice = initial_state[field_slot].data[:, :, iz_center]
    im = ax.imshow(
        initial_xy_slice.T,
        origin="lower",
        extent=[0, domain_size, 0, domain_size],
        cmap="RdBu_r",
        vmin=-amplitude,
        vmax=amplitude,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"$\phi$ in x-y plane (t=0, z=center)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # --- Panel 3: 2D x-y slice (final) ---
    ax = axes[1, 0]
    final_xy_slice = final[field_slot].data[:, :, iz_center]
    vmax_final = max(float(np.max(np.abs(final_xy_slice))), 0.01)
    im = ax.imshow(
        final_xy_slice.T,
        origin="lower",
        extent=[0, domain_size, 0, domain_size],
        cmap="RdBu_r",
        vmin=-vmax_final,
        vmax=vmax_final,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(rf"$\phi$ in x-y plane (t={times[-1]:.1f}, z=center)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # --- Panel 4: max|phi| over time ---
    ax = axes[1, 1]
    ax.plot(times, max_phi_over_time, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"max $|\phi|$")
    ax.set_title("Peak amplitude decay")
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        f"Klein-Gordon 3+1D: m$^2$={m2_value}, "
        f"grid={n}$^3$, dt={dt}",
        fontsize=14,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"  Saved plot to: {output_path}")

    print()
    print("*** 3+1D Klein-Gordon simulation complete! ***")


if __name__ == "__main__":
    main()
