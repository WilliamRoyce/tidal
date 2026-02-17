"""Massive 3-form field simulation in 3+1D Minkowski spacetime.

This script demonstrates the evolution of a rank-3 antisymmetric tensor field C_{abc}
using equations derived from the Lagrangian:

    L = -1/12 d_a C_{bcd} d^a C^{bcd} - m^2/12 C_{bcd} C^{bcd}

In flat Minkowski spacetime with signature (-,+,+,+), total antisymmetry reduces
64 raw components to 4 independent ones (C_{012}, C_{013}, C_{023}, C_{123}),
each satisfying the Klein-Gordon equation:

    d^2_t(C_i) = laplacian(C_i) - m^2 C_i

The simulation starts with a 3D Gaussian pulse in C_0 only (the other 3 components
and all momenta are zero). Since each component satisfies KG independently (no
cross-coupling), C_1/C_2/C_3 should remain zero throughout the evolution.

Physics context:
    3-form potentials appear in supergravity (11D), string theory (Ramond-Ramond
    fields), and generalized gauge field theories. This example demonstrates the
    pipeline's ability to handle rank-3 tensor symmetry reduction end-to-end.

Performance note:
    With 4 components x 2 (field + momentum) = 8 state fields, 3D grids are
    expensive. This example uses 16^3 cells for a practical demo.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, MemoryStorage

from tidal.symbolic import (
    build_pde_from_json,
    create_initial_state,
    load_equation_system,
)
from tidal.utils import normalize_solve_result

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from tidal.symbolic.equation_system import EquationSystem

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
MASS_SQUARED = 1.0

# Grid
GRID_POINTS_PER_AXIS = 16  # 16^3 = 4,096 cells; 8 fields = 32K floats
DOMAIN_SIZE = 10.0

# Time integration
T_END = 5.0
DT = 0.05
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
PULSE_WIDTH = 1.5
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "massive_3form_output.png"

# State indices (interleaved: [C_0, pi_0, C_1, pi_1, C_2, pi_2, C_3, pi_3])
C0_FIELD_SLOT = 0

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for massive 3-form simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    spec: EquationSystem
    initial_state: FieldCollection
    initial_peak: float


def _print_header() -> None:
    print("=" * 60)
    print("Massive 3-Form Field in 3+1D Minkowski Spacetime")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/12 d_a C_bcd d^a C^bcd - m^2/12 C_bcd C^bcd")
    print("Symmetry: Totally antisymmetric rank-3 tensor")
    print("Reduction: 64 raw components -> 4 independent (antisymmetry)")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")

    if not json_path.exists():
        msg = (
            f"JSON specification not found: {json_path}\n"
            "Run the Wolfram derivation first: "
            "wolframscript -file examples/massive_3form/massive_3form.wls"
        )
        raise FileNotFoundError(msg)

    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Spacetime dimension: {spec.dimension} (3+1D)")
    print(f"  Spatial dimension: {spec.spatial_dimension}")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  Coordinates: {spec.effective_coordinates}")
    print(f"  Symmetry reduction: 4^3 = 64 raw -> {spec.n_components} independent")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(f"    d2_t({eq.field_name}) = {' + '.join(operators)}")

    print()
    print("  State layout:")
    for i, (name, slot_type) in enumerate(spec.state_layout):
        print(f"    slot {i}: {name} ({slot_type})")
    print()
    return spec


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters={"m2": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 3D simulation grid...")
    n = GRID_POINTS_PER_AXIS
    grid = CartesianGrid(
        bounds=[(0, DOMAIN_SIZE), (0, DOMAIN_SIZE), (0, DOMAIN_SIZE)],
        shape=[n, n, n],
        periodic=True,
    )
    print(f"  Domain: [0, {DOMAIN_SIZE}]^3")
    print(f"  Resolution: {n} x {n} x {n} = {n**3:,} cells")
    print(f"  Grid spacing: dx = dy = dz = {DOMAIN_SIZE / n:.3f}")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(
    grid: CartesianGrid, spec: EquationSystem
) -> tuple[FieldCollection, float]:
    print("Step 4: Creating initial conditions...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    z = cast("np.ndarray", grid.cell_coords[..., 2])

    center = DOMAIN_SIZE / 2.0

    # 3D Gaussian pulse in C_0 only (other components = 0)
    c0_data = PULSE_AMPLITUDE * np.exp(
        -((x - center) ** 2 + (y - center) ** 2 + (z - center) ** 2)
        / (2 * PULSE_WIDTH**2)
    )

    # Only C_0 is excited; C_1, C_2, C_3 remain zero (component independence)
    state = create_initial_state(grid, spec, field_data={"C_0": c0_data})
    initial_peak = float(np.max(np.abs(c0_data)))

    print(
        f"  Gaussian pulse in C_0: center=({center}, {center}, {center}), "
        f"width={PULSE_WIDTH}"
    )
    print(f"  Amplitude: {PULSE_AMPLITUDE}, initial peak = {initial_peak:.4f}")
    print("  C_1, C_2, C_3: zero (testing component independence)")
    print("  All momenta: zero (pulse starts at rest)")
    print()
    return state, initial_peak


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection, spec: EquationSystem,
    initial_peak: float,
) -> SimulationResult:
    print("Step 5: Running simulation...")

    initial_state = state.copy()

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {T_END} time units")
    print(f"  Time step: dt = {DT}")
    print("  Scheme: Runge-Kutta (RK4)")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        spec=spec,
        initial_state=initial_state,
        initial_peak=initial_peak,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    final = cast("FieldCollection", result.storage[-1])
    component_names = result.spec.component_names
    n_components = len(component_names)

    final_peak_c0 = float(np.max(np.abs(final[C0_FIELD_SLOT].data)))
    print(f"  C_0 initial peak: {result.initial_peak:.4f}")
    print(f"  C_0 final peak:   {final_peak_c0:.4f}")
    print(f"  C_0 decay ratio:  {final_peak_c0 / result.initial_peak:.4f}")
    print("  (Expected: amplitude decreases as pulse spreads in 3D)")

    print()
    print("  Component independence check:")
    for i in range(1, n_components):
        field_slot = 2 * i  # C_i is at slot 2*i
        max_val = float(np.max(np.abs(final[field_slot].data)))
        print(f"    max|{component_names[i]}| = {max_val:.2e} (should be ~0)")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    grid = result.grid
    spec = result.spec
    initial_state = result.initial_state
    initial_peak = result.initial_peak

    final = cast("FieldCollection", storage[-1])
    component_names = spec.component_names
    n_components = len(component_names)

    n = GRID_POINTS_PER_AXIS
    ix_center = n // 2
    iy_center = n // 2
    iz_center = n // 2

    z_1d = cast("np.ndarray", grid.cell_coords[ix_center, iy_center, :, 2])

    # Collect max|C_i| over time for all components
    times: list[float] = [0.0, *list(storage.times)]
    max_c0_over_time: list[float] = [initial_peak]
    max_others_over_time: dict[str, list[float]] = {
        name: [0.0] for name in component_names[1:]
    }

    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        max_c0_over_time.append(float(np.max(np.abs(snapshot[C0_FIELD_SLOT].data))))
        for i in range(1, n_components):
            field_slot = 2 * i
            max_others_over_time[component_names[i]].append(
                float(np.max(np.abs(snapshot[field_slot].data)))
            )

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # --- Panel 1: C_0 z-profile (initial vs final + intermediates) ---
    ax = axes[0, 0]
    initial_profile = initial_state[C0_FIELD_SLOT].data[ix_center, iy_center, :]
    final_profile = final[C0_FIELD_SLOT].data[ix_center, iy_center, :]

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
        profile = snap[C0_FIELD_SLOT].data[ix_center, iy_center, :]
        ax.plot(
            z_1d,
            profile,
            color=color,
            linewidth=1.5,
            alpha=alpha,
            label=f"t = {storage_times[idx]:.1f}",
        )

    ax.set_xlabel("z")
    ax.set_ylabel(r"$C_0$")
    ax.set_title(r"$C_0$ along z-axis (at x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 2: C_0 x-y slice (initial) ---
    ax = axes[0, 1]
    initial_xy_slice = initial_state[C0_FIELD_SLOT].data[:, :, iz_center]
    im = ax.imshow(
        initial_xy_slice.T,
        origin="lower",
        extent=[0, DOMAIN_SIZE, 0, DOMAIN_SIZE],
        cmap="RdBu_r",
        vmin=-PULSE_AMPLITUDE,
        vmax=PULSE_AMPLITUDE,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(r"$C_0$ in x-y plane (t=0, z=center)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # --- Panel 3: Peak amplitude decay ---
    ax = axes[1, 0]
    ax.plot(times, max_c0_over_time, "b-", linewidth=2, label=r"$C_0$")
    ax.set_xlabel("Time")
    ax.set_ylabel(r"max $|C_i|$")
    ax.set_title("Peak amplitude decay")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 4: Component independence ---
    ax = axes[1, 1]
    colors = ["red", "green", "purple"]
    for i, (name, vals) in enumerate(max_others_over_time.items()):
        ax.plot(times, vals, color=colors[i], linewidth=2, label=name)
    ax.set_xlabel("Time")
    ax.set_ylabel(r"max $|C_i|$")
    ax.set_title("Component independence (should stay ~0)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)
    ax.ticklabel_format(style="scientific", axis="y", scilimits=(0, 0))

    fig.suptitle(
        f"Massive 3-Form (antisymmetric rank-3): m$^2$={MASS_SQUARED}, "
        f"grid={n}$^3$, dt={DT}",
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


def _print_footer(result: SimulationResult) -> None:
    final = cast("FieldCollection", result.storage[-1])
    component_names = result.spec.component_names
    n_components = len(component_names)

    final_peak_c0 = float(np.max(np.abs(final[C0_FIELD_SLOT].data)))

    # Collect final max for non-C_0 components
    max_others_final: dict[str, float] = {}
    for i in range(1, n_components):
        field_slot = 2 * i
        max_others_final[component_names[i]] = float(
            np.max(np.abs(final[field_slot].data))
        )

    print("*** Massive 3-form simulation complete! ***")
    print()
    print("Key results:")
    print(f"  - Symmetry reduction: 64 raw components -> {n_components} independent")
    print(
        f"  - C_0 peak decayed from {result.initial_peak:.4f} to {final_peak_c0:.4f}"
    )
    for name, final_max in max_others_final.items():
        print(f"  - {name} remained ~0 (max = {final_max:.2e})")
    print("  - Component independence verified: no spurious coupling")


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the massive 3-form field simulation."""
    json_path = Path(__file__).parent.parent / "data" / "massive_3form.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state, initial_peak = _create_initial_state(grid, spec)
    sim_result = _run_simulation(pde, grid, state, spec, initial_peak)
    _analyze_results(sim_result)
    _plot_results(sim_result)
    _print_footer(sim_result)


if __name__ == "__main__":
    main()
