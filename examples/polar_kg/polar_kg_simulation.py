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
from pde import CallbackTracker, CartesianGrid, FieldCollection, ScalarField

from tidal.measurement import (
    EnergyDiagnostics,
    SimulationData,
    check_energy_conservation,
    create_snapshot_callback,
)
from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
MASS_SQUARED = 0.5

# Grid
R_MIN = 0.5
R_MAX = 10.0
THETA_MIN = 0.0
THETA_MAX = 2 * np.pi
GRID_SHAPE = [128, 128]
GRID_PERIODIC: list[bool] = [False, True]

# Time integration
T_END = 8.0
DT = 0.005
SNAPSHOT_INTERVAL = 0.2

# Initial conditions
PULSE_RADIUS = 3.0
PULSE_WIDTH = 0.5
PULSE_AMPLITUDE = 1.0

# Analysis
ENERGY_THRESHOLD = 0.15

# Output
OUTPUT_FILENAME = "polar_kg_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for polar KG simulation results."""

    grid: CartesianGrid
    data: SimulationData
    r_coords: NumericArray
    theta_coords: NumericArray


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


def _build_pde(json_path: Path) -> PDEFromSpec:
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
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  r range: [{R_MIN}, {R_MAX}]")
    print("  theta range: [0, 2*pi]")
    print(f"  Resolution: {GRID_SHAPE[0]} x {GRID_SHAPE[1]}")
    print("  Boundary: r=Neumann, theta=periodic")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    r = cast("np.ndarray", grid.cell_coords[..., 0])

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
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection, json_path: Path
) -> SimulationResult:
    print("Step 5: Running simulation...")

    spec = load_equation_system(json_path)
    params = {"polm2": MASS_SQUARED}
    output_data_dir = Path(__file__).parent.parent / "data" / "polar_kg_output"

    writer, callback = create_snapshot_callback(
        output_dir=output_data_dir,
        spec=spec,
        grid=grid,
        t_end=T_END,
        snapshot_interval=SNAPSHOT_INTERVAL,
        parameters=params,
        spec_path=json_path,
    )
    tracker = CallbackTracker(callback, interrupts=SNAPSHOT_INTERVAL)
    pde.solve(  # type: ignore[union-attr]
        state,
        t_range=T_END,
        dt=DT,
        solver="scipy",
        method="RK45",
        tracker=tracker,
    )
    writer.close()

    data = SimulationData.from_directory(output_data_dir, spec)
    r = cast("np.ndarray", grid.cell_coords[..., 0])
    theta = cast("np.ndarray", grid.cell_coords[..., 1])

    print(f"  Duration: {T_END} time units, solver=scipy/RK45")
    print(f"  {data.n_snapshots} snapshots saved to {output_data_dir}")
    print()

    return SimulationResult(
        grid=grid,
        data=data,
        r_coords=r,
        theta_coords=theta,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    data = result.data
    field_name = data.field_names[0]

    initial_max = float(np.max(np.abs(data.fields[field_name][0])))
    final_max = float(np.max(np.abs(data.fields[field_name][-1])))

    print(f"  Initial max|phi| = {initial_max:.4f}")
    print(f"  Final   max|phi| = {final_max:.4f}")
    print()


def _plot_results(result: SimulationResult, diag: EnergyDiagnostics) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    data = result.data
    grid = result.grid
    field_name = data.field_names[0]
    n_snapshots = data.n_snapshots

    initial_phi = data.fields[field_name][0]
    final_phi = data.fields[field_name][-1]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Initial field in native (r, theta) coordinates using imshow
    ax = axes[0, 0]
    vmax = PULSE_AMPLITUDE
    im = ax.imshow(
        initial_phi.T,
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
    final_vmax = max(float(np.max(np.abs(final_phi))), 0.01)
    im = ax.imshow(
        final_phi.T,
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
    time_values = list(data.times)
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
        ax.plot(
            r_1d,
            data.fields[field_name][t_idx][:, 0],  # theta=0 slice
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
    r_idx = GRID_SHAPE[0] // 3  # Pick a representative radius
    r_value = r_1d[r_idx]
    theta_1d = cast("np.ndarray", grid.cell_coords[0, :, 1])
    for i, t_idx in enumerate(times_to_plot):
        ax.plot(
            theta_1d / np.pi,
            data.fields[field_name][t_idx][r_idx, :],
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    _ = diag  # energy diagnostics available for future use
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


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the polar coordinate Klein-Gordon simulation.

    Raises
    ------
    FileNotFoundError
        If the JSON spec has not been generated yet.
    """
    json_path = Path(__file__).parent.parent / "data" / "polar_kg.json"
    if not json_path.exists():
        msg = (
            f"{json_path} not found. "
            "Run 'tidal derive theory.toml' first (requires wolframscript)."
        )
        raise FileNotFoundError(msg)
    _print_header()
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid)
    result = _run_simulation(pde, grid, state, json_path)
    _analyze_results(result)

    print("Checking energy conservation...")
    diag = check_energy_conservation(result.data, threshold=ENERGY_THRESHOLD)
    print(f"  max |dE/E0| = {diag.max_relative_error:.2e}")
    print(f"  Conserved: {diag.is_conserved}")
    print()

    _plot_results(result, diag)
    _print_footer()


if __name__ == "__main__":
    main()
