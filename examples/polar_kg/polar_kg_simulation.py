"""Klein-Gordon in polar coordinates - Christoffel corrections from metric.

Demonstrates the full Lagrangian->JSON->PDE pipeline with a non-Cartesian
coordinate metric.  The Wolfram script derives the wave equation on flat 2D
space in polar coordinates (r, theta), automatically producing 1/r gradient
and 1/r^2 angular Laplacian corrections from the Christoffel symbols.

Physics:
    Metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2   (x=r, y=theta)
    EOM:    d2_t phi = d2_r phi + (1/r) d_r phi + (1/r^2) d2_theta phi - m^2 phi

    The (1/r) gradient term causes cylindrical wave spreading.
    The angular Laplacian is suppressed at small r.

Boundary conditions:
    r direction:     Neumann (non-periodic) at r_min and r_max
    theta direction: periodic [0, 2*pi]
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

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system
from torsion_gertsenshtein.utils import normalize_solve_result

if TYPE_CHECKING:
    from numpy.typing import NDArray

    NumericArray = NDArray[np.float64]

OUTPUT_FILENAME = "polar_kg_output.png"
R_MIN = 0.5
R_MAX = 10.0
THETA_MIN = 0.0
THETA_MAX = 2 * np.pi
GRID_NR = 128
GRID_NTHETA = 128
MASS_SQUARED = 0.5
PULSE_RADIUS = 3.0
PULSE_WIDTH = 0.5
PULSE_AMPLITUDE = 1.0
T_END = 8.0
DT = 0.005


@dataclass(frozen=True)
class SimulationResult:
    """Container for polar KG simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    r_coords: NumericArray
    theta_coords: NumericArray


def main() -> None:
    """Run the polar coordinate Klein-Gordon simulation."""
    _print_header()
    json_path = Path(__file__).parent.parent / "data" / "polar_kg.json"
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid)
    result = _run_simulation(pde, grid, state)
    _analyze_results(result)
    _plot_results(result)
    _print_footer()


def _print_header() -> None:
    print("=" * 60)
    print("Klein-Gordon in Polar Coordinates (r, theta)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2")
    print("  using coordinates {t, x, y} = {t, r, theta}")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (2+1D)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        n_terms = len(eq.rhs_terms)
        n_pos_dep = sum(1 for t in eq.rhs_terms if t.position_dependent)
        operators = [t.operator for t in eq.rhs_terms]
        print(
            f"    d2_t({eq.field_name}): {n_terms} terms, "
            f"{n_pos_dep} position-dependent"
        )
        print(f"      Operators: {operators}")
    print()


def _build_pde(json_path: Path) -> object:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters={"polm2": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up polar coordinate grid...")
    grid = CartesianGrid(
        bounds=[(R_MIN, R_MAX), (THETA_MIN, THETA_MAX)],
        shape=[GRID_NR, GRID_NTHETA],
        periodic=[False, True],  # r: bounded, theta: periodic
    )
    print(f"  r range: [{R_MIN}, {R_MAX}]")
    print("  theta range: [0, 2*pi]")
    print(f"  Resolution: {GRID_NR} x {GRID_NTHETA}")
    print("  Boundary: r=Neumann, theta=periodic")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    r = cast("np.ndarray", grid.cell_coords[..., 0])
    cast("np.ndarray", grid.cell_coords[..., 1])

    # Gaussian ring at r = PULSE_RADIUS, uniform in theta
    gaussian_ring = PULSE_AMPLITUDE * np.exp(
        -((r - PULSE_RADIUS) ** 2) / (2 * PULSE_WIDTH**2)
    )

    phi = ScalarField(grid, data=gaussian_ring, label="polphi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")
    state = FieldCollection([phi, pi])

    print(f"  Gaussian ring at r={PULSE_RADIUS}, width={PULSE_WIDTH}")
    print(f"  Amplitude: {PULSE_AMPLITUDE}")
    print(f"  Initial max|phi| = {np.max(np.abs(gaussian_ring)):.4f}")
    print()
    return state


def _run_simulation(
    pde: object, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")

    storage = MemoryStorage()
    result = pde.solve(  # type: ignore[union-attr]
        state,
        t_range=T_END,
        dt=DT,
        solver="scipy",
        method="RK45",
        tracker=storage.tracker(0.2),
    )
    result = normalize_solve_result(result)

    r = cast("np.ndarray", grid.cell_coords[..., 0])
    theta = cast("np.ndarray", grid.cell_coords[..., 1])

    print(f"  Duration: {T_END} time units, solver=scipy/RK45")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        r_coords=r,
        theta_coords=theta,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])

    initial_max = float(np.max(np.abs(initial[0].data)))
    final_max = float(np.max(np.abs(final[0].data)))

    print(f"  Initial max|phi| = {initial_max:.4f}")
    print(f"  Final   max|phi| = {final_max:.4f}")

    # Energy proxy (should be roughly conserved for flat space, static metric)
    def compute_energy(snapshot: FieldCollection) -> float:
        phi_data = snapshot[0].data
        pi_data = snapshot[1].data
        return float(np.sum(pi_data**2 + phi_data**2))

    initial_energy = compute_energy(initial)
    final_energy = compute_energy(final)
    energy_change = abs(final_energy - initial_energy) / max(initial_energy, 1e-10)

    print(f"  Initial energy proxy: {initial_energy:.2f}")
    print(f"  Final energy proxy:   {final_energy:.2f}")
    print(f"  Relative change:      {energy_change:.4f}")
    if energy_change < 0.15:
        print("  Energy approximately conserved (flat space, static metric)")
    else:
        print("  Note: Energy change may be due to boundary effects")
    print()


def _plot_results(result: SimulationResult) -> None:
    print("Step 7: Generating visualization...")

    storage = result.storage
    grid = result.grid

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initial field in native (r, theta) coordinates using imshow
    ax = axes[0, 0]
    vmax = PULSE_AMPLITUDE
    im = ax.imshow(
        initial[0].data.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-vmax,
        vmax=vmax,
        extent=[R_MIN, R_MAX, THETA_MIN, THETA_MAX],
        aspect="auto",
        interpolation="bilinear",
    )
    ax.set_xlabel("r")
    ax.set_ylabel("theta")
    ax.set_title("Initial phi (t=0)")
    fig.colorbar(im, ax=ax, label="phi")

    # Final field in native (r, theta) coordinates
    ax = axes[0, 1]
    final_vmax = max(float(np.max(np.abs(final[0].data))), 0.01)
    im = ax.imshow(
        final[0].data.T,
        origin="lower",
        cmap="bwr_r",
        vmin=-final_vmax,
        vmax=final_vmax,
        extent=[R_MIN, R_MAX, THETA_MIN, THETA_MAX],
        aspect="auto",
        interpolation="bilinear",
    )
    ax.set_xlabel("r")
    ax.set_ylabel("theta")
    ax.set_title(f"Final phi (t={T_END:.0f})")
    fig.colorbar(im, ax=ax, label="phi")

    # Radial cross-section at theta=0 over time
    ax = axes[1, 0]
    time_values = list(storage.times)
    n_snapshots = len(storage)
    times_to_plot = [
        0,
        n_snapshots // 4,
        n_snapshots // 2,
        3 * n_snapshots // 4,
        n_snapshots - 1,
    ]
    cmap_colors = plt.get_cmap("viridis")
    colors = cmap_colors(np.linspace(0.2, 0.8, len(times_to_plot)))

    r_1d = cast("np.ndarray", grid.cell_coords[:, 0, 0])
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            r_1d,
            snapshot[0].data[:, 0],  # theta=0 slice
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.set_xlabel("r")
    ax.set_ylabel("phi")
    ax.set_title("Radial cross-section (theta=0)")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    # Angular cross-section at fixed r over time
    ax = axes[1, 1]
    r_idx = GRID_NR // 3  # Pick a representative radius
    r_value = r_1d[r_idx]
    theta_1d = cast("np.ndarray", grid.cell_coords[0, :, 1])
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            theta_1d / np.pi,
            snapshot[0].data[r_idx, :],
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.set_xlabel("theta / pi")
    ax.set_ylabel("phi")
    ax.set_title(f"Angular cross-section (r={r_value:.1f})")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        "KG in Polar Coordinates: "
        r"$ds^2 = -dt^2 + dr^2 + r^2 d\theta^2$"
        f"\nMass m^2={MASS_SQUARED}, ring pulse at r={PULSE_RADIUS}",
        fontsize=12,
    )
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
    print("Polar coordinate simulation complete!")
    print()
    print("Key observations:")
    print("  1. Metric ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2")
    print("  2. Christoffel correction: (1/r) gradient term")
    print("  3. Angular Laplacian scaled by 1/r^2")
    print("  4. All coefficients derived from metric by xAct pipeline")
    print("  5. Mixed BCs: Neumann (r) + periodic (theta)")
    print("=" * 60)


if __name__ == "__main__":
    main()
