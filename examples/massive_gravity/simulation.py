"""Linearized massive gravity in 2+1D simulation (Fierz-Pauli mass term).

This script demonstrates massive gravity in 2+1D with the Fierz-Pauli
mass term — the unique ghost-free linear mass term for spin-2:

  G^(1)_ab[h] - m^2 (h_ab - eta_ab h) = 0

where h = eta^cd h_cd is the trace.  The Fierz-Pauli combination
(h_ab h^ab - h^2) propagates 5 DOF; wrong relative coefficient activates
a Boulware-Deser ghost (6th DOF, negative energy).

The trace term creates off-diagonal mass coupling between diagonal metric
components (h_0, h_3, h_5) via h = -h_0 + h_3 + h_5.  Off-diagonal
components (h_1, h_2, h_4) keep simple -m^2 self-coupling.

In 2+1D, pure GR has no local degrees of freedom; the mass term creates
a propagating massive mode with dispersion relation: omega^2 = k^2 + m^2.

Component structure (symmetric rank-2 in 3D, gauge-unfixed):
  h_0 (h_tt):         constraint (time_order=0, trace-coupled to h_3, h_5)
  h_1 (h_tx):         constraint (time_order=0, partial Helmholtz)
  h_2 (h_ty):         constraint (time_order=0, partial Helmholtz)
  h_3 (h_xx):         constraint (time_order=0, trace-coupled to h_0, h_5)
  h_4 (h_xy):         EVOLUTION (time_order=2, the only propagating DOF)
  h_5 (h_yy):         constraint (time_order=0, trace-coupled to h_0, h_3)

CONSTRAINT SOLVING (unified solver):
  All 5 constraint equations are solved at each timestep using the unified
  constraint solver with FFT block solve for coupled constraints:

  - Cluster {h_0, h_3, h_5}: mutually coupled via Fierz-Pauli trace terms
    AND Laplacian cross-terms. Solved simultaneously in Fourier space as a
    3x3 dense system at each wavenumber.

  - Cluster {h_1, h_2}: mutually coupled. h_1 is partial Helmholtz
    (-m2*I - 0.5*L_y), h_2 is (-m2*I - 0.5*L_x), with cross_derivative
    coupling. Their sources are zero at t=0 (depend on pi_3, pi_4, pi_5).

  The constraint coupling feeds back into h_4 via cross_derivative_xy(h_0),
  adding spatial structure to the evolution. In this gauge, h_4 oscillates
  at approximately omega = sqrt(2*m^2) with small corrections from
  constraints. True wave propagation requires gauge-fixing (de Donder).
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
MASS_SQUARED = 1.0

# Grid
DOMAIN_X = 50.0
DOMAIN_Y = 50.0
NX = 64
NY = 64

# Time integration
OMEGA = np.sqrt(2 * MASS_SQUARED)
PERIOD = 2 * np.pi / OMEGA
T_END = 3.5 * PERIOD
DT = 0.01
N_SNAPSHOTS = 200
SNAPSHOT_INTERVAL = T_END / N_SNAPSHOTS

# Initial conditions
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 0.5

# Output
OUTPUT_FILENAME = "massive_gravity_output.png"

# State indices
EVOLVING_FIELD = "h_4"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for massive gravity simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    h4_idx: int
    time_values: list[float]
    center_values: list[float]
    max_amplitudes: list[float]
    analytic: NumericArray
    residual: float


def _print_header() -> None:
    print("=" * 60)
    print("Linearized Massive Gravity 2+1D Simulation")
    print("=" * 60)
    print()
    print("Equation: G^(1)_ab[h] - m^2 (h_ab - eta_ab h) = 0  (Fierz-Pauli)")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Dimension: {spec.dimension} (2+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        time_order = eq.time_derivative_order
        eq_type = "evolution" if time_order >= 2 else "constraint"  # noqa: PLR2004
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(
            f"    {eq.field_name}: d{time_order}_t = "
            f"{' + '.join(operators)}  [{eq_type}]"
        )

    print()
    print(f"  Oscillation frequency: omega = sqrt(2*m^2) = {OMEGA:.4f}")
    print(f"  Period: T = 2*pi/omega = {PERIOD:.4f}")
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters={"m2": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 2D simulation grid...")
    grid = CartesianGrid(
        bounds=[(0, DOMAIN_X), (0, DOMAIN_Y)],
        shape=[NX, NY],
        periodic=True,
    )
    print(f"  Domain: [0, {DOMAIN_X}] x [0, {DOMAIN_Y}]")
    print(f"  Resolution: {NX} x {NY}")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(
    grid: CartesianGrid, json_path: Path
) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    spec = load_equation_system(json_path)
    x_coords = cast("np.ndarray", grid.cell_coords[..., 0])
    y_coords = cast("np.ndarray", grid.cell_coords[..., 1])

    center_x, center_y = DOMAIN_X / 2, DOMAIN_Y / 2
    gaussian = PULSE_AMPLITUDE * np.exp(
        -(
            (x_coords - center_x) ** 2
            + (y_coords - center_y) ** 2
        )
        / (2 * PULSE_WIDTH**2)
    )

    fields: list[ScalarField] = []
    for name, slot_type in spec.state_layout:
        label = f"{name}_{slot_type}"
        sf = ScalarField(grid, data=0.0, label=label)
        if name == EVOLVING_FIELD and slot_type == "field":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)
    print(f"  Gaussian pulse in {EVOLVING_FIELD} (h_xy) at ({center_x}, {center_y})")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print(f"  State layout: {spec.state_layout}")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection, json_path: Path
) -> SimulationResult:
    print("Step 5: Running simulation...")

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {T_END:.2f} time units ({T_END / PERIOD:.1f} periods)")
    print(f"  Stored {len(storage)} snapshots")
    print()

    # Find h_4 field index in state layout
    spec = load_equation_system(json_path)
    field_slot_map = {
        name: i
        for i, (name, stype) in enumerate(spec.state_layout)
        if stype == "field"
    }
    h4_idx = field_slot_map[EVOLVING_FIELD]

    mid_x_idx = NX // 2
    mid_y_idx = NY // 2

    # Track center value and max amplitude over time
    center_values: list[float] = []
    max_amplitudes: list[float] = []
    time_values: list[float] = []
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        center_values.append(
            float(snapshot[h4_idx].data[mid_x_idx, mid_y_idx])
        )
        max_amplitudes.append(
            float(np.max(np.abs(snapshot[h4_idx].data)))
        )
        time_values.append(float(t_idx) * SNAPSHOT_INTERVAL)

    times = np.array(time_values)
    centers = np.array(center_values)

    # Verify oscillation matches expected frequency
    analytic = PULSE_AMPLITUDE * np.cos(OMEGA * times)
    residual = float(np.max(np.abs(centers - analytic)))

    return SimulationResult(
        grid=grid,
        storage=storage,
        h4_idx=h4_idx,
        time_values=time_values,
        center_values=center_values,
        max_amplitudes=max_amplitudes,
        analytic=analytic,
        residual=residual,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    mid_x_idx = NX // 2
    mid_y_idx = NY // 2

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])

    initial_center = float(
        initial[result.h4_idx].data[mid_x_idx, mid_y_idx]
    )
    final_center = float(
        final[result.h4_idx].data[mid_x_idx, mid_y_idx]
    )

    print(f"  Initial h_4 at center: {initial_center:.6f}")
    print(f"  Final h_4 at center:   {final_center:.6f}")
    print(f"  Expected: A*cos(omega*t), omega = {OMEGA:.4f}")
    print(f"  Max residual vs analytic: {result.residual:.2e}")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    h4_idx = result.h4_idx
    time_values = result.time_values
    times = np.array(time_values)
    centers = np.array(result.center_values)
    max_amplitudes = result.max_amplitudes
    analytic = result.analytic

    def find_nearest_snapshot(target_t: float) -> int:
        return int(np.argmin(np.abs(np.array(time_values) - target_t)))

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # ---- Row 1: Three h_4 heatmap snapshots at key oscillation phases ----
    snapshot_times = [0.0, PERIOD / 4, PERIOD / 2]
    snapshot_labels = ["t = 0 (peak +)", "t = T/4 (node)", "t = T/2 (peak -)"]

    for col, (snap_t, snap_label) in enumerate(
        zip(snapshot_times, snapshot_labels, strict=True)
    ):
        ax = axes[0, col]
        snap_idx = find_nearest_snapshot(snap_t)
        snap_data = cast("FieldCollection", storage[snap_idx])
        im = ax.imshow(
            snap_data[h4_idx].data.T,
            origin="lower",
            extent=[0, DOMAIN_X, 0, DOMAIN_Y],
            cmap="RdBu_r",
            vmin=-PULSE_AMPLITUDE,
            vmax=PULSE_AMPLITUDE,
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
        label=rf"$A\cos(\omega t)$, $\omega$={OMEGA:.3f}",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("h_4 at center")
    ax.set_title("Center-point oscillation")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)
    # Mark oscillation period
    for n in range(1, 4):
        ax.axvline(n * PERIOD, color="gray", linestyle=":", alpha=0.4)

    # ---- Row 2, Center: Cross-sections at multiple times ----
    ax = axes[1, 1]
    x_1d = np.linspace(0, DOMAIN_X, NX)
    mid_y_idx = NY // 2
    cross_times = [0.0, PERIOD / 4, PERIOD / 2, 3 * PERIOD / 4, PERIOD]
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
    ax.set_title(f"Cross-section at y = {DOMAIN_Y / 2}")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # ---- Row 2, Right: Peak amplitude envelope ----
    ax = axes[1, 2]
    ax.plot(times, max_amplitudes, "b-", linewidth=2, label="max|h_4|")
    ax.axhline(
        PULSE_AMPLITUDE,
        color="gray",
        linestyle="--",
        alpha=0.5,
        label=f"A = {PULSE_AMPLITUDE}",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel("max|h_4|")
    ax.set_title("Peak amplitude (stability check)")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)
    ax.set_ylim(0, PULSE_AMPLITUDE * 1.3)

    fig.suptitle(
        f"Massive Gravity 2+1D: m$^2$ = {MASS_SQUARED}   |   "
        rf"$\omega \approx \sqrt{{2m^2}}$ = {OMEGA:.3f}   |   "
        f"T = {PERIOD:.2f}\n"
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
    plt.close()
    print()


def _print_footer(residual: float) -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. h_4 (h_xy) is the ONLY evolution equation (time_order=2)")
    print("  2. h_0-h_3, h_5 are constraints (solved via unified solver)")
    print(f"  3. Oscillation at omega ~ sqrt(2*m^2) = {OMEGA:.4f}")
    print(
        f"  4. Period T = {PERIOD:.4f}, simulation covers "
        f"{T_END / PERIOD:.1f} periods"
    )
    print(f"  5. Max residual vs cos(omega*t): {residual:.2e}")
    print()
    print("  NOTE: Constraints feed back via cross_derivative_xy(h_0),")
    print("  but in gauge-unfixed form, this does not produce spatial")
    print("  propagation. True wave propagation requires gauge-fixing.")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the massive gravity simulation."""
    json_path = Path(__file__).parent.parent / "data" / "massive_gravity_3d.json"
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid, json_path)
    result = _run_simulation(pde, grid, state, json_path)
    _analyze_results(result)
    _plot_results(result)
    _print_footer(result.residual)


if __name__ == "__main__":
    main()
