"""Gravitational wave propagation simulation from linearized Einstein-Hilbert action.

This script demonstrates gravitational wave propagation in 3+1D using the
de Donder gauge equations derived from the linearized Einstein-Hilbert Lagrangian:
  L = R * sqrt(-g) + gauge fixing term

In the de Donder (harmonic) gauge, each component of the metric perturbation
h_{ab} satisfies a flat-space wave equation:
  d2_t(h_i) = laplacian_x(h_i) + laplacian_y(h_i) + laplacian_z(h_i)

The 10 independent components of the symmetric h_{ab} are indexed as:
  h_0 = h_{tt}, h_1 = h_{tx}, h_2 = h_{ty}, h_3 = h_{tz}
  h_4 = h_{xx}, h_5 = h_{xy}, h_6 = h_{xz}
  h_7 = h_{yy}, h_8 = h_{yz}
  h_9 = h_{zz}

For transverse-traceless (TT) gauge gravitational waves propagating along z:
  - h_+ polarization: h_{xx} = -h_{yy} (h_4 = -h_7)
  - h_x polarization: h_{xy} = h_{yx}  (h_5)
  - All other components vanish

The simulation uses TT-gauge-compatible initial data to demonstrate
wave packet propagation along the z-axis.
"""

from pathlib import Path
from typing import cast

OUTPUT_FILENAME = "gw_output.png"

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


def main() -> None:  # noqa: PLR0915, PLR0914
    """Run the gravitational wave propagation simulation."""
    print("=" * 60)
    print("Gravitational Wave Propagation (Linearized EH, de Donder Gauge)")
    print("=" * 60)
    print()
    print("Lagrangian: R * sqrt(-g) + gauge fixing term")
    print("Gauge: de Donder (harmonic)")
    print()

    # ----------------------------------------------------------------
    # Step 1: Load the equation specification
    # ----------------------------------------------------------------
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "linearized_gravity_dedonder.json"

    if not json_path.exists():
        msg = (
            f"JSON specification not found: {json_path}\n"
            "Run the Wolfram derivation first: "
            "wolframscript -file examples/gravitational_waves/linearized_gravity.wls"
        )
        raise FileNotFoundError(msg)

    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Gauge: {spec.metadata.get('gauge', 'N/A')}")
    print(f"  Spacetime dimension: {spec.dimension} (3+1D)")
    print(f"  Spatial dimension: {spec.spatial_dimension}")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Coordinates: {spec.effective_coordinates}")

    # Print equation structure
    print()
    print("  Equation structure:")
    for eq in spec.equations:
        time_order = eq.time_derivative_order
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(f"    d{time_order}_t({eq.field_name}) = {' + '.join(operators)}")

    # Print state layout
    print()
    print("  State layout:")
    for i, (name, slot_type) in enumerate(spec.state_layout):
        print(f"    slot {i}: {name} ({slot_type})")

    # ----------------------------------------------------------------
    # Step 2: Build the PDE
    # ----------------------------------------------------------------
    print()
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print(f"  State size: {spec.state_size} fields (each component has field + momentum)")

    # ----------------------------------------------------------------
    # Step 3: Create the 3D spatial grid
    # ----------------------------------------------------------------
    print()
    print("Step 3: Setting up 3D simulation grid...")

    # Coarse grid: narrow in x,y (transverse), extended in z (propagation direction)
    # 4x4x64 gives reasonable resolution for z-propagation demo
    nx, ny, nz = 4, 4, 64
    lx, ly, lz = 4.0, 4.0, 40.0

    grid = CartesianGrid(
        bounds=[(0, lx), (0, ly), (0, lz)],
        shape=[nx, ny, nz],
        periodic=True,
    )
    print(f"  Domain: [0, {lx}] x [0, {ly}] x [0, {lz}]")
    print(f"  Resolution: {nx} x {ny} x {nz} points")
    print(f"  Grid spacing: dx={lx / nx:.2f}, dy={ly / ny:.2f}, dz={lz / nz:.3f}")
    print("  Periodic boundary conditions")

    # ----------------------------------------------------------------
    # Step 4: Set up TT-gauge initial data
    # ----------------------------------------------------------------
    print()
    print("Step 4: Creating TT-gauge initial conditions...")

    # Coordinate arrays for the 3D grid
    z = cast("np.ndarray", grid.cell_coords[..., 2])

    # Gaussian wave packet propagating along z
    z_center = lz / 2.0
    z_width = 3.0
    amplitude = 1.0
    k_z = 2.0 * np.pi / 10.0  # wavenumber for oscillation within envelope

    # h_+ polarization: Gaussian envelope * sinusoidal carrier along z
    h_plus = amplitude * np.exp(-((z - z_center) ** 2) / (2 * z_width**2)) * np.cos(
        k_z * (z - z_center)
    )

    # h_x polarization: smaller amplitude, phase-shifted
    h_cross_amplitude = 0.3 * amplitude
    h_cross = h_cross_amplitude * np.exp(
        -((z - z_center) ** 2) / (2 * z_width**2)
    ) * np.sin(k_z * (z - z_center))

    # TT-gauge conditions:
    #   h_4 (h_xx) = h_+
    #   h_7 (h_yy) = -h_+ (tracelessness: h_xx + h_yy = 0)
    #   h_5 (h_xy) = h_x
    #   All other components = 0 (transverse: no h_{tz}, h_{xz}, etc.)
    field_data = {
        "h_4": h_plus,
        "h_7": -h_plus,
        "h_5": h_cross,
    }

    state = create_initial_state(grid, spec, field_data=field_data)

    # Keep a copy of the true initial state (before any time stepping)
    initial_state = state.copy()

    # Validate TT-gauge conditions on initial data
    _field_slots = {
        name: i
        for i, (name, stype) in enumerate(spec.state_layout)
        if stype == "field"
    }
    _init = cast("FieldCollection", initial_state)
    trace_violation = float(
        np.max(np.abs(_init[_field_slots["h_4"]].data + _init[_field_slots["h_7"]].data))
    )
    if trace_violation > 1e-12:
        msg = (
            f"Initial data violates TT tracelessness: |h_xx + h_yy| = {trace_violation:.2e}. "
            "For TT gauge, require h_4 + h_7 = 0."
        )
        raise ValueError(msg)

    print(f"  h_+ polarization: Gaussian(center={z_center}, width={z_width})")
    print(f"    amplitude = {amplitude}, wavenumber k_z = {k_z:.4f}")
    print(f"  h_x polarization: amplitude = {h_cross_amplitude}")
    print("  TT-gauge: h_4 (h_xx) = h_+, h_7 (h_yy) = -h_+, h_5 (h_xy) = h_x")
    print("  All other components initialized to zero")

    # ----------------------------------------------------------------
    # Step 5: Run simulation
    # ----------------------------------------------------------------
    print()
    print("Step 5: Running simulation...")

    t_end = 15.0
    dt = 0.01
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
    print(f"  Time step: {dt}")
    print(f"  Scheme: Runge-Kutta (RK4)")
    print(f"  Stored {len(storage)} snapshots")

    # ----------------------------------------------------------------
    # Step 6: Analyze results
    # ----------------------------------------------------------------
    print()
    print("Step 6: Analyzing results...")

    # Use the true initial state (not storage[0], which may be at t>0)
    initial = cast("FieldCollection", initial_state)
    final = cast("FieldCollection", storage[-1])

    # Find state slot indices for field data (not momentum)
    field_slot_map = {
        name: i
        for i, (name, stype) in enumerate(spec.state_layout)
        if stype == "field"
    }
    h4_slot = field_slot_map["h_4"]
    h7_slot = field_slot_map["h_7"]
    h5_slot = field_slot_map["h_5"]

    initial_h4_max = np.max(np.abs(initial[h4_slot].data))
    final_h4_max = np.max(np.abs(final[h4_slot].data))
    initial_h5_max = np.max(np.abs(initial[h5_slot].data))
    final_h5_max = np.max(np.abs(final[h5_slot].data))

    print(f"  Initial: max|h_+| = {initial_h4_max:.4f}, max|h_x| = {initial_h5_max:.4f}")
    print(f"  Final:   max|h_+| = {final_h4_max:.4f}, max|h_x| = {final_h5_max:.4f}")

    # Check tracelessness is maintained: h_4 + h_7 should stay near zero
    trace_initial = np.max(np.abs(initial[h4_slot].data + initial[h7_slot].data))
    trace_final = np.max(np.abs(final[h4_slot].data + final[h7_slot].data))
    print(f"  Tracelessness |h_xx + h_yy|: initial={trace_initial:.2e}, final={trace_final:.2e}")

    # Check that gauge-violating components stay small
    h0_final = np.max(np.abs(final[field_slot_map["h_0"]].data))
    h1_final = np.max(np.abs(final[field_slot_map["h_1"]].data))
    h2_final = np.max(np.abs(final[field_slot_map["h_2"]].data))
    h3_final = np.max(np.abs(final[field_slot_map["h_3"]].data))
    print(
        f"  Non-TT components (should be ~0): "
        f"|h_0|={h0_final:.2e}, |h_1|={h1_final:.2e}, "
        f"|h_2|={h2_final:.2e}, |h_3|={h3_final:.2e}"
    )

    # ----------------------------------------------------------------
    # Step 7: Collect evolution data for plotting
    # ----------------------------------------------------------------
    print()
    print("Step 7: Generating visualization...")

    # Extract z-profile at x=0, y=0 (transverse center)
    # For the 3D grid with shape (nx, ny, nz), take slice [0, 0, :]
    ix_center = 0
    iy_center = 0

    z_1d = cast("np.ndarray", grid.cell_coords[ix_center, iy_center, :, 2])

    # Collect max |h_+| over time for panel 3
    # Start with the true initial state at t=0
    times = [0.0, *list(storage.times)]
    max_hplus_over_time = [float(np.max(np.abs(initial[h4_slot].data)))]
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        max_hplus_over_time.append(float(np.max(np.abs(snapshot[h4_slot].data))))

    # ----------------------------------------------------------------
    # Step 8: Create plots
    # ----------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # --- Panel 1: h_+ z-profile (initial vs final) ---
    ax = axes[0]
    initial_profile = initial[h4_slot].data[ix_center, iy_center, :]
    final_profile = final[h4_slot].data[ix_center, iy_center, :]

    ax.plot(z_1d, initial_profile, "b-", linewidth=2, label="t = 0.0 (initial)")
    ax.plot(z_1d, final_profile, "r-", linewidth=2, label=f"t = {times[-1]:.1f} (final)")

    # Also plot a few intermediate snapshots from storage
    n_snapshots = len(storage)
    storage_times = list(storage.times)
    for frac, color, alpha in [(0.25, "green", 0.5), (0.5, "orange", 0.5)]:
        idx = int(frac * (n_snapshots - 1))
        snap = cast("FieldCollection", storage[idx])
        profile = snap[h4_slot].data[ix_center, iy_center, :]
        ax.plot(
            z_1d, profile, color=color, linewidth=1.5, alpha=alpha,
            label=f"t = {storage_times[idx]:.1f}",
        )

    ax.set_xlabel("z")
    ax.set_ylabel(r"$h_+$ ($h_{xx}$)")
    ax.set_title(r"$h_+$ polarization z-profile")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 2: h_x z-profile ---
    ax = axes[1]
    initial_cross = initial[h5_slot].data[ix_center, iy_center, :]
    final_cross = final[h5_slot].data[ix_center, iy_center, :]

    ax.plot(z_1d, initial_cross, "b-", linewidth=2, label="t = 0.0 (initial)")
    ax.plot(z_1d, final_cross, "r-", linewidth=2, label=f"t = {times[-1]:.1f} (final)")

    for frac, color, alpha in [(0.25, "green", 0.5), (0.5, "orange", 0.5)]:
        idx = int(frac * (n_snapshots - 1))
        snap = cast("FieldCollection", storage[idx])
        profile = snap[h5_slot].data[ix_center, iy_center, :]
        ax.plot(
            z_1d, profile, color=color, linewidth=1.5, alpha=alpha,
            label=f"t = {storage_times[idx]:.1f}",
        )

    ax.set_xlabel("z")
    ax.set_ylabel(r"$h_\times$ ($h_{xy}$)")
    ax.set_title(r"$h_\times$ polarization z-profile")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 3: Max |h_+| over time ---
    ax = axes[2]
    ax.plot(times, max_hplus_over_time, "b-", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"max $|h_+|$")
    ax.set_title(r"Peak $|h_+|$ amplitude over time")
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        "Gravitational Wave Propagation (Linearized EH, de Donder gauge)\n"
        r"TT-gauge initial data: $h_+ = h_{xx} = -h_{yy}$, "
        r"$h_\times = h_{xy}$, propagation along z",
        fontsize=11,
    )
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved plot to: {output_path}")

    # ----------------------------------------------------------------
    # Summary
    # ----------------------------------------------------------------
    print()
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations derived from linearized Einstein-Hilbert Lagrangian")
    print("     via xAct/Mathematica with de Donder gauge fixing")
    print("  2. Each metric perturbation component satisfies a massless wave equation")
    print("  3. TT-gauge initial data: h_+ and h_x polarizations propagate along z")
    print("  4. Tracelessness h_xx + h_yy = 0 is preserved by evolution")
    print("  5. Components initialized to zero remain zero (decoupled wave equations)")
    print(f"  6. Wave speed = 1 (speed of light), propagation distance ~ {t_end:.0f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
