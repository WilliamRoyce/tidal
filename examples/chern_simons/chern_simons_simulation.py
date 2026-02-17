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

# Grid
X_MIN = 0.0
X_MAX = 50.0
Y_MIN = 0.0
Y_MAX = 50.0
GRID_SHAPE = [64, 64]

# Time integration
T_END = 10.0
DT = 0.01
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
CENTER_X = 25.0
CENTER_Y = 25.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Analysis
CS_COUPLING_THRESHOLD = 0.01

# Output
OUTPUT_FILENAME = "chern_simons_output.png"

# State indices (A_0 is first-order, A_1/A_2 are second-order)
IDX_A0 = 0
IDX_A1 = 1
IDX_PI1 = 2
IDX_A2 = 3
IDX_PI2 = 4

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for Chern-Simons simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    x_coords: NumericArray
    y_coords: NumericArray


def _print_header() -> None:
    print("=" * 60)
    print("Chern-Simons 2+1D Simulation from Lagrangian")
    print("=" * 60)
    print()
    print("Lagrangian: L = -1/4 F_ab F^ab + (kappa/2) epsilon^abc A_a D_b A_c")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
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
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 2D simulation grid...")
    grid = CartesianGrid(
        bounds=[(X_MIN, X_MAX), (Y_MIN, Y_MAX)],
        shape=GRID_SHAPE,
        periodic=True,
    )
    print(f"  Domain: [{grid.axes_bounds[0]}, {grid.axes_bounds[1]}]")
    print(f"  Resolution: {grid.shape}")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    # A_0 is first-order (no momentum slot), A_1 and A_2 are second-order
    # State layout: A_0, A_1, pi_1, A_2, pi_2 (5 fields total)

    a0 = ScalarField(grid, data=0.0, label="A_0")
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

    gaussian = PULSE_AMPLITUDE * np.exp(
        -((x - CENTER_X) ** 2 + (y - CENTER_Y) ** 2) / (2 * PULSE_WIDTH**2)
    )
    a1.data[:] = gaussian

    state = FieldCollection([a0, a1, pi1, a2, pi2])
    print(f"  Gaussian pulse in A_1 at center ({CENTER_X}, {CENTER_Y})")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")

    storage = MemoryStorage()
    result = pde.solve(  # type: ignore[union-attr]
        state,
        t_range=T_END,
        dt=DT,
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    print(f"  Duration: {T_END} time units")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        x_coords=x,
        y_coords=y,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])

    # Extract field amplitudes (state layout: A_0, A_1, pi_1, A_2, pi_2)
    initial_a0 = np.max(np.abs(initial[IDX_A0].data))
    initial_a1 = np.max(np.abs(initial[IDX_A1].data))
    initial_a2 = np.max(np.abs(initial[IDX_A2].data))

    final_a0 = np.max(np.abs(final[IDX_A0].data))
    final_a1 = np.max(np.abs(final[IDX_A1].data))
    final_a2 = np.max(np.abs(final[IDX_A2].data))

    print(
        f"  Initial: max|A_0| = {initial_a0:.3f}, "
        f"max|A_1| = {initial_a1:.3f}, max|A_2| = {initial_a2:.3f}"
    )
    print(
        f"  Final:   max|A_0| = {final_a0:.3f}, "
        f"max|A_1| = {final_a1:.3f}, max|A_2| = {final_a2:.3f}"
    )

    # Check for energy transfer (CS coupling effect)
    if final_a0 > CS_COUPLING_THRESHOLD or final_a2 > CS_COUPLING_THRESHOLD:
        print("  CS coupling effect observed: energy transferred to A_0 and A_2")
    else:
        print("  Note: CS coupling may be weak at this kappa value")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914
    print("Step 7: Generating visualization...")

    storage = result.storage

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    initial_a1 = float(np.max(np.abs(initial[IDX_A1].data)))
    final_a0 = float(np.max(np.abs(final[IDX_A0].data)))
    final_a1 = float(np.max(np.abs(final[IDX_A1].data)))
    final_a2 = float(np.max(np.abs(final[IDX_A2].data)))

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Initial state (state layout: A_0, A_1, pi_1, A_2, pi_2)
    vmax = max(initial_a1, 0.1)
    im0 = axes[0, 0].imshow(
        initial[IDX_A0].data.T, origin="lower", cmap="bwr", vmin=-vmax, vmax=vmax
    )
    axes[0, 0].set_title("Initial A_0")
    plt.colorbar(im0, ax=axes[0, 0])

    im1 = axes[0, 1].imshow(
        initial[IDX_A1].data.T, origin="lower", cmap="bwr", vmin=-vmax, vmax=vmax
    )
    axes[0, 1].set_title("Initial A_1")
    plt.colorbar(im1, ax=axes[0, 1])

    im2 = axes[0, 2].imshow(
        initial[IDX_A2].data.T, origin="lower", cmap="bwr", vmin=-vmax, vmax=vmax
    )
    axes[0, 2].set_title("Initial A_2")
    plt.colorbar(im2, ax=axes[0, 2])

    # Final state
    vmax_f = max(final_a0, final_a1, final_a2, 0.1)
    im3 = axes[1, 0].imshow(
        final[IDX_A0].data.T, origin="lower", cmap="bwr", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 0].set_title(f"Final A_0 (t={T_END})")
    plt.colorbar(im3, ax=axes[1, 0])

    im4 = axes[1, 1].imshow(
        final[IDX_A1].data.T, origin="lower", cmap="bwr", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 1].set_title(f"Final A_1 (t={T_END})")
    plt.colorbar(im4, ax=axes[1, 1])

    im5 = axes[1, 2].imshow(
        final[IDX_A2].data.T, origin="lower", cmap="bwr", vmin=-vmax_f, vmax=vmax_f
    )
    axes[1, 2].set_title(f"Final A_2 (t={T_END})")
    plt.colorbar(im5, ax=axes[1, 2])

    fig.suptitle("Chern-Simons 2+1D: Maxwell-CS with kappa=0.5")
    plt.tight_layout()

    output_dir = Path(__file__).parent.parent.parent / "outputs"
    output_dir.mkdir(exist_ok=True, parents=True)
    output_path = output_dir / OUTPUT_FILENAME
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"  Saved plot to: {output_path}")
    plt.close()
    print()


def _print_footer() -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Equations derived from CS Lagrangian (Maxwell + topological term)")
    print("  2. Cross-field coupling via gradient operators")
    print("  3. Energy can transfer between A components")
    print("  4. Chern-Simons term adds effective mass to photon")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the Chern-Simons simulation."""
    json_path = Path(__file__).parent.parent / "data" / "chern_simons_3d.json"
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid)
    result = _run_simulation(pde, grid, state)
    _analyze_results(result)
    _plot_results(result)
    _print_footer()


if __name__ == "__main__":
    main()
