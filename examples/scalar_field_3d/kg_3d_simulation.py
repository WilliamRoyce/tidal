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

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
MASS_SQUARED = 1.0

# Grid
DOMAIN_SIZE = 20.0
N_GRID = 32  # Grid points per axis (32^3 = 32,768 cells)

# Time integration
T_END = 8.0
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
PULSE_WIDTH = 2.0
PULSE_AMPLITUDE = 1.0

# Analysis
PROFILE_FRACTIONS = [0.25, 0.5]
PROFILE_COLORS = ["green", "orange"]
PROFILE_ALPHA = 0.5

# Output
OUTPUT_FILENAME = "kg_3d_output.png"

# State indices
FIELD_SLOT = 0

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for 3+1D Klein-Gordon simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    initial_state: FieldCollection
    initial_peak: float
    final_peak: float


def _print_header() -> None:
    print("=" * 60)
    print("Klein-Gordon 3+1D (Massive Scalar Field)")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/2 (nabla phi)^2 - 1/2 m^2 phi^2")
    print("Spacetime: 3+1D flat Minkowski, signature (-,+,+,+)")
    print()


def _load_spec(json_path: Path) -> None:
    """Load and display the equation specification.

    Raises
    ------
    FileNotFoundError
        If the JSON specification file is not found.
    """
    print("Step 1: Loading equation specification...")

    if not json_path.exists():
        msg = (
            f"JSON specification not found: {json_path}\n"
            "Run the Wolfram derivation first: "
            "tidal derive examples/scalar_field_3d/theory.toml"
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
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    spec = load_equation_system(json_path)
    pde = build_pde_from_json(json_path, parameters={"m2": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print(f"  State size: {spec.state_size} fields (field + momentum)")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 3D simulation grid...")
    grid = CartesianGrid(
        bounds=[
            (0, DOMAIN_SIZE),
            (0, DOMAIN_SIZE),
            (0, DOMAIN_SIZE),
        ],
        shape=[N_GRID, N_GRID, N_GRID],
        periodic=True,
    )
    print(f"  Domain: [0, {DOMAIN_SIZE}]^3")
    print(f"  Resolution: {N_GRID} x {N_GRID} x {N_GRID} = {N_GRID**3:,} cells")
    print(f"  Grid spacing: dx = dy = dz = {DOMAIN_SIZE / N_GRID:.3f}")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(
    grid: CartesianGrid, json_path: Path
) -> tuple[FieldCollection, float]:
    print("Step 4: Creating initial conditions...")

    spec = load_equation_system(json_path)

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    z = cast("np.ndarray", grid.cell_coords[..., 2])

    center = DOMAIN_SIZE / 2.0

    # 3D Gaussian pulse centered in the domain, at rest (momentum = 0)
    phi_data = PULSE_AMPLITUDE * np.exp(
        -(
            (x - center) ** 2
            + (y - center) ** 2
            + (z - center) ** 2
        )
        / (2 * PULSE_WIDTH**2)
    )

    state = create_initial_state(
        grid, spec, field_data={"phi_0": phi_data}
    )

    initial_peak = float(np.max(np.abs(phi_data)))

    print(
        f"  Gaussian pulse: center=({center}, {center}, {center}), "
        f"width={PULSE_WIDTH}"
    )
    print(f"  Amplitude: {PULSE_AMPLITUDE}, initial peak = {initial_peak:.4f}")
    print("  Momentum: zero (pulse starts at rest)")
    print()
    return state, initial_peak


def _run_simulation(
    pde: PDEFromSpec,
    grid: CartesianGrid,
    state: FieldCollection,
    initial_peak: float,
) -> SimulationResult:
    print("Step 5: Running simulation...")

    initial_state = state.copy()

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        solver="scipy",
        method="DOP853",
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    print(f"  Duration: {T_END} time units")
    print("  Solver: scipy/DOP853 (adaptive)")
    print(f"  Stored {len(storage)} snapshots")
    print()

    final = cast("FieldCollection", storage[-1])
    final_peak = float(np.max(np.abs(final[FIELD_SLOT].data)))

    return SimulationResult(
        grid=grid,
        storage=storage,
        initial_state=initial_state,
        initial_peak=initial_peak,
        final_peak=final_peak,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    print(f"  Initial peak |phi|: {result.initial_peak:.4f}")
    print(f"  Final peak |phi|:   {result.final_peak:.4f}")
    print(f"  Decay ratio: {result.final_peak / result.initial_peak:.4f}")
    print("  (Expected: amplitude decreases as pulse spreads in 3D)")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    grid = result.grid
    initial_state = result.initial_state

    final = cast("FieldCollection", storage[-1])

    ix_center = N_GRID // 2
    iy_center = N_GRID // 2
    iz_center = N_GRID // 2

    z_1d = cast("np.ndarray", grid.cell_coords[ix_center, iy_center, :, 2])

    # Collect max|phi| over time
    times = [0.0, *list(storage.times)]
    max_phi_over_time = [result.initial_peak]
    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        max_phi_over_time.append(
            float(np.max(np.abs(snapshot[FIELD_SLOT].data)))
        )

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # --- Panel 1: phi z-profile (initial vs final) ---
    ax = axes[0, 0]
    initial_profile = initial_state[FIELD_SLOT].data[ix_center, iy_center, :]
    final_profile = final[FIELD_SLOT].data[ix_center, iy_center, :]

    ax.plot(z_1d, initial_profile, "b-", linewidth=2, label="t = 0.0 (initial)")
    ax.plot(
        z_1d,
        final_profile,
        "r-",
        linewidth=2,
        label=f"t = {times[-1]:.1f} (final)",
    )

    # Add intermediate snapshots
    n_snapshots = len(storage)
    storage_times = list(storage.times)
    for frac, color in zip(PROFILE_FRACTIONS, PROFILE_COLORS, strict=True):
        idx = int(frac * (n_snapshots - 1))
        snap = cast("FieldCollection", storage[idx])
        profile = snap[FIELD_SLOT].data[ix_center, iy_center, :]
        ax.plot(
            z_1d,
            profile,
            color=color,
            linewidth=1.5,
            alpha=PROFILE_ALPHA,
            label=f"t = {storage_times[idx]:.1f}",
        )

    ax.set_xlabel("z")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"$\phi$ along z-axis (at x=y=center)")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    # --- Panel 2: 2D x-y slice (initial) ---
    ax = axes[0, 1]
    initial_xy_slice = initial_state[FIELD_SLOT].data[:, :, iz_center]
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
    ax.set_title(r"$\phi$ in x-y plane (t=0, z=center)")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # --- Panel 3: 2D x-y slice (final) ---
    ax = axes[1, 0]
    final_xy_slice = final[FIELD_SLOT].data[:, :, iz_center]
    vmax_final = max(float(np.max(np.abs(final_xy_slice))), 0.01)
    im = ax.imshow(
        final_xy_slice.T,
        origin="lower",
        extent=[0, DOMAIN_SIZE, 0, DOMAIN_SIZE],
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
        f"Klein-Gordon 3+1D: m$^2$={MASS_SQUARED}, "
        f"grid={N_GRID}$^3$, solver=scipy/DOP853 (adaptive)",
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


def _print_footer() -> None:
    print("*** 3+1D Klein-Gordon simulation complete! ***")


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the 3+1D Klein-Gordon simulation."""
    json_path = Path(__file__).parent.parent / "data" / "klein_gordon_3d.json"
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state, initial_peak = _create_initial_state(grid, json_path)
    result = _run_simulation(pde, grid, state, initial_peak)
    _analyze_results(result)
    _plot_results(result)
    _print_footer()


if __name__ == "__main__":
    main()
