"""Linearized massive gravity in 2+1D simulation.

This script demonstrates massive gravity in 2+1D:
  G^(1)_ab[h] - m^2 h_ab = 0

The mass term is obtained by linearizing G_ab - m^2 g_ab = 0 around flat
Minkowski. The MINUS sign is essential: with "+", the identity terms have
positive coefficient in the evolution equations (exponential growth); with
"-", we get d2_t(h) = -m^2 h + spatial (stable Klein-Gordon-like).

In 2+1D, pure GR has no local degrees of freedom; the mass term creates
a propagating massive mode with dispersion relation: omega^2 = k^2 + m^2.

Component structure (symmetric rank-2 in 3D, gauge-unfixed):
  h_0 (h_tt):         constraint (time_order=0, Helmholtz-type)
  h_1 (h_tx):         constraint (time_order=0, partial Helmholtz)
  h_2 (h_ty):         constraint (time_order=0, partial Helmholtz)
  h_3 (h_xx):         constraint (time_order=0, algebraic)
  h_4 (h_xy):         EVOLUTION (time_order=2, the only propagating DOF)
  h_5 (h_yy):         constraint (time_order=0, algebraic)

CONSTRAINT SOLVING (unified solver):
  All 5 constraint equations are solved at each timestep using the unified
  constraint solver with FFT block solve for coupled constraints:

  - Cluster {h_0, h_3, h_5}: mutually coupled via Laplacian cross-terms.
    Solved simultaneously in Fourier space as a 3x3 dense system at each
    wavenumber. h_0 is Helmholtz-type (-m2*I - L), h_3 and h_5 are
    algebraic (-m2*I), coupled through Laplacians.

  - Cluster {h_1, h_2}: mutually coupled. h_1 is partial Helmholtz
    (-m2*I - 0.5*L_y), h_2 is (-m2*I - 0.5*L_x), with cross_derivative
    coupling. Their sources are zero at t=0 (depend on pi_3, pi_4, pi_5).

  The constraint coupling feeds back into h_4 via cross_derivative_xy(h_0),
  adding spatial structure to the evolution. In this gauge, h_4 oscillates
  at approximately omega = sqrt(2*m^2) with small corrections from
  constraints. True wave propagation requires gauge-fixing (de Donder).
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
    print("Equation: G^(1)_ab[h] - m^2 h_ab = 0")
    print()

    # Load the equation specification
    print("Step 1: Loading equation specification...")
    json_path = Path(__file__).parent.parent / "data" / "massive_gravity_3d.json"
    mass_squared = 1.0
    params = {"m2": mass_squared}
    spec = load_equation_system(json_path)

    print(f"  Dimension: {spec.dimension} (2+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Mass parameter: m^2 = {mass_squared}")

    # Show equation structure
    print()
    print("  Equation structure:")
    for eq in spec.equations:
        time_order = eq.time_derivative_order
        eq_type = "evolution" if time_order >= 2 else "constraint"  # noqa: PLR2004
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(f"    {eq.field_name}: d{time_order}_t = {' + '.join(operators)}  [{eq_type}]")

    # Physics parameters
    omega = np.sqrt(2 * mass_squared)  # oscillation frequency
    period = 2 * np.pi / omega

    print()
    print(f"  Oscillation frequency: omega = sqrt(2*m^2) = {omega:.4f}")
    print(f"  Period: T = 2*pi/omega = {period:.4f}")

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

    # Gaussian pulse in h_4 (h_xy — the ONLY evolution equation)
    center_x, center_y = lx / 2, ly / 2
    width = 5.0
    amplitude = 0.5
    gaussian = amplitude * np.exp(
        -((x_coords - center_x) ** 2 + (y_coords - center_y) ** 2) / (2 * width**2)
    )

    # Build state according to spec layout
    for name, slot_type in spec.state_layout:
        label = f"{name}_{slot_type}"
        sf = ScalarField(grid, data=0.0, label=label)
        if name == "h_4" and slot_type == "field":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)
    print(f"  Gaussian pulse in h_4 (h_xy) at ({center_x}, {center_y})")
    print(f"  Width: {width}, Amplitude: {amplitude}")
    print(f"  State layout: {spec.state_layout}")

    # Run simulation — cover ~3 oscillation periods
    print()
    print("Step 5: Running simulation...")
    t_end = 3.5 * period  # ~3.5 oscillation periods
    dt = 0.01
    snapshot_interval = t_end / 200  # ~200 snapshots for smooth curves

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=t_end,
        dt=dt,
        scheme="runge-kutta",
        tracker=storage.tracker(snapshot_interval),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {t_end:.2f} time units ({t_end / period:.1f} periods)")
    print(f"  Stored {len(storage)} snapshots")

    # Analyze results
    print()
    print("Step 6: Analyzing results...")

    # Find h_4 field index in state layout
    field_slot_map = {
        name: i for i, (name, stype) in enumerate(spec.state_layout) if stype == "field"
    }
    h4_idx = field_slot_map["h_4"]

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])
    mid_x_idx = nx // 2
    mid_y_idx = ny // 2

    initial_center = float(initial[h4_idx].data[mid_x_idx, mid_y_idx])
    final_center = float(final[h4_idx].data[mid_x_idx, mid_y_idx])

    print(f"  Initial h_4 at center: {initial_center:.6f}")
    print(f"  Final h_4 at center:   {final_center:.6f}")

    # Track center value and max amplitude over time
    center_values: list[float] = []
    max_amplitudes: list[float] = []
    time_values: list[float] = []
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        center_values.append(float(snapshot[h4_idx].data[mid_x_idx, mid_y_idx]))
        max_amplitudes.append(float(np.max(np.abs(snapshot[h4_idx].data))))
        time_values.append(float(t_idx) * snapshot_interval)

    times = np.array(time_values)
    centers = np.array(center_values)

    # Verify oscillation matches expected frequency
    analytic = amplitude * np.cos(omega * times)
    print(f"  Expected: A*cos(omega*t), omega = {omega:.4f}")
    residual = float(np.max(np.abs(centers - analytic)))
    print(f"  Max residual vs analytic: {residual:.2e}")

    # Find snapshot indices for key oscillation phases
    def find_nearest_snapshot(target_t: float) -> int:
        return int(np.argmin(np.abs(np.array(time_values) - target_t)))

    # Generate visualization
    print()
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # ---- Row 1: Three h_4 heatmap snapshots at key oscillation phases ----
    snapshot_times = [0.0, period / 4, period / 2]
    snapshot_labels = ["t = 0 (peak +)", "t = T/4 (node)", "t = T/2 (peak -)"]

    for col, (snap_t, snap_label) in enumerate(zip(snapshot_times, snapshot_labels, strict=True)):
        ax = axes[0, col]
        snap_idx = find_nearest_snapshot(snap_t)
        snap_data = cast("FieldCollection", storage[snap_idx])
        im = ax.imshow(
            snap_data[h4_idx].data.T,
            origin="lower",
            extent=[0, lx, 0, ly],
            cmap="RdBu_r",
            vmin=-amplitude,
            vmax=amplitude,
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        actual_t = time_values[snap_idx]
        ax.set_title(f"h_4 (h_xy) — {snap_label}\n(t = {actual_t:.2f})")
        plt.colorbar(im, ax=ax, shrink=0.8)

    # ---- Row 2, Left: Center amplitude vs time with analytic overlay ----
    ax = axes[1, 0]
    ax.plot(times, centers, "b-", linewidth=2, label="Simulation")
    ax.plot(
        times,
        analytic,
        "r--",
        linewidth=1.5,
        alpha=0.7,
        label=rf"$A\cos(\omega t)$, $\omega$={omega:.3f}",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("h_4 at center")
    ax.set_title("Center-point oscillation")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)
    # Mark oscillation period
    for n in range(1, 4):
        ax.axvline(n * period, color="gray", linestyle=":", alpha=0.4)

    # ---- Row 2, Center: Cross-sections at multiple times ----
    ax = axes[1, 1]
    x_1d = np.linspace(0, lx, nx)
    cross_times = [0.0, period / 4, period / 2, 3 * period / 4, period]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0, 1, len(cross_times)))
    for ct, color in zip(cross_times, colors, strict=True):
        cross_idx = find_nearest_snapshot(ct)
        snap_data = cast("FieldCollection", storage[cross_idx])
        actual_t = time_values[cross_idx]
        ax.plot(
            x_1d,
            snap_data[h4_idx].data[:, mid_y_idx],
            color=color,
            linewidth=1.5,
            label=f"t = {actual_t:.1f}",
        )
    ax.set_xlabel("x")
    ax.set_ylabel("h_4")
    ax.set_title(f"Cross-section at y = {ly / 2}")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # ---- Row 2, Right: Peak amplitude envelope ----
    ax = axes[1, 2]
    ax.plot(times, max_amplitudes, "b-", linewidth=2, label="max|h_4|")
    ax.axhline(amplitude, color="gray", linestyle="--", alpha=0.5, label=f"A = {amplitude}")
    ax.set_xlabel("Time")
    ax.set_ylabel("max|h_4|")
    ax.set_title("Peak amplitude (stability check)")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)
    ax.set_ylim(0, amplitude * 1.3)

    fig.suptitle(
        f"Massive Gravity 2+1D: m$^2$ = {mass_squared}   |   "
        rf"$\omega \approx \sqrt{{2m^2}}$ = {omega:.3f}   |   "
        f"T = {period:.2f}\n"
        r"$\ddot{h}_{xy} = -2m^2 h_{xy} + \partial_{xy} h_{tt}$  "
        "(gauge-unfixed, constraints solved via FFT block)",
        fontsize=11,
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
    print("  1. h_4 (h_xy) is the ONLY evolution equation (time_order=2)")
    print("  2. h_0-h_3, h_5 are constraints (solved via unified solver)")
    print(f"  3. Oscillation at omega ~ sqrt(2*m^2) = {omega:.4f}")
    print(f"  4. Period T = {period:.4f}, simulation covers {t_end / period:.1f} periods")
    print(f"  5. Max residual vs cos(omega*t): {residual:.2e}")
    print()
    print("  NOTE: Constraints feed back via cross_derivative_xy(h_0),")
    print("  but in gauge-unfixed form, this does not produce spatial")
    print("  propagation. True wave propagation requires gauge-fixing.")
    print("=" * 60)


if __name__ == "__main__":
    main()
