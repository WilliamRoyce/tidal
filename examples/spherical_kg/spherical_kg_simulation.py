"""Klein-Gordon in spherical coordinates - 3+1D with trig Christoffel terms.

Demonstrates the full Lagrangian->JSON->PDE pipeline with a 3+1D spherical
coordinate metric.  The Wolfram script derives the wave equation on flat 3D
space in spherical coordinates (r, theta, phi), producing position-dependent
Christoffel corrections including trigonometric functions (Cot, Csc).

Physics:
    Metric: ds^2 = -dt^2 + dx^2 + x^2 dy^2 + x^2 sin^2(y) dz^2
            where x=r, y=theta, z=phi

    EOM: d2_t phi = d2_r phi + (2/r) d_r phi
                  + (1/r^2) d2_theta phi + (cot(theta)/r^2) d_theta phi
                  + (1/(r^2 sin^2(theta))) d2_phi phi - m^2 phi

    Key features:
    - 6 RHS terms with position-dependent coefficients
    - Trigonometric functions: Cot[y[]], Csc[y[]]^2
    - 1/r amplitude decay for radial propagation
    - Coordinate singularity avoided at theta=0 and theta=pi

Boundary conditions:
    r direction:     Neumann (non-periodic)
    theta direction: Neumann (non-periodic, avoiding poles)
    phi direction:   periodic [0, 2*pi]
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

    from torsion_gertsenshtein.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

OUTPUT_FILENAME = "spherical_kg_output.png"
R_MIN = 0.5
R_MAX = 8.0
THETA_MIN = 0.05  # Avoid theta=0 singularity (sin(0)=0)
THETA_MAX = np.pi - 0.05  # Avoid theta=pi singularity
PHI_MIN = 0.0
PHI_MAX = 2 * np.pi
GRID_NR = 64
GRID_NTHETA = 64
GRID_NPHI = 64
MASS_SQUARED = 0.0  # Massless for clean 1/r decay test
PULSE_RADIUS = 3.0
PULSE_WIDTH = 0.6
PULSE_AMPLITUDE = 1.0
T_END = 5.0
DT = 0.01


@dataclass(frozen=True)
class SimulationResult:
    """Container for spherical KG simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    r_coords: NumericArray
    theta_coords: NumericArray
    phi_coords: NumericArray


def main() -> None:
    """Run the spherical coordinate Klein-Gordon simulation."""
    _print_header()
    json_path = Path(__file__).parent.parent / "data" / "spherical_kg.json"
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
    print("Klein-Gordon in Spherical Coordinates (r, theta, phi)")
    print("=" * 60)
    print()
    print("Metric: ds^2 = -dt^2 + dr^2 + r^2 d(theta)^2 + r^2 sin^2(theta) d(phi)^2")
    print("  using coordinates {t, x, y, z} = {t, r, theta, phi}")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (3+1D)")
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
    pde = build_pde_from_json(json_path, parameters={"spm2": MASS_SQUARED})
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Mass parameter: m^2 = {MASS_SQUARED}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up spherical coordinate grid...")
    grid = CartesianGrid(
        bounds=[(R_MIN, R_MAX), (THETA_MIN, THETA_MAX), (PHI_MIN, PHI_MAX)],
        shape=[GRID_NR, GRID_NTHETA, GRID_NPHI],
        periodic=[False, False, True],  # r: bounded, theta: bounded, phi: periodic
    )
    print(f"  r range: [{R_MIN}, {R_MAX}]")
    print(f"  theta range: [{THETA_MIN:.1f}, {THETA_MAX:.1f}] (avoiding poles)")
    print(f"  phi range: [{PHI_MIN}, {PHI_MAX}]")
    print(f"  Resolution: {GRID_NR} x {GRID_NTHETA} x {GRID_NPHI}")
    print(f"  Total cells: {GRID_NR * GRID_NTHETA * GRID_NPHI}")
    print("  Boundary: r=Neumann, theta=Neumann, phi=periodic")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    r = cast("np.ndarray", grid.cell_coords[..., 0])

    # Gaussian shell at r = PULSE_RADIUS, uniform in theta and phi
    gaussian_shell = PULSE_AMPLITUDE * np.exp(
        -((r - PULSE_RADIUS) ** 2) / (2 * PULSE_WIDTH**2)
    )

    phi = ScalarField(grid, data=gaussian_shell, label="spphi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")
    state = FieldCollection([phi, pi])

    print(f"  Gaussian shell at r={PULSE_RADIUS}, width={PULSE_WIDTH}")
    print(f"  Amplitude: {PULSE_AMPLITUDE}")
    print(f"  Initial max|phi| = {np.max(np.abs(gaussian_shell)):.4f}")
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
        dt=DT,  # Let the solver choose adaptive time steps
        solver="scipy",
        method="RK45",
        tracker=storage.tracker(0.5),
    )
    result = normalize_solve_result(result)

    r = cast("np.ndarray", grid.cell_coords[..., 0])
    theta = cast("np.ndarray", grid.cell_coords[..., 1])
    phi_coord = cast("np.ndarray", grid.cell_coords[..., 2])

    print(f"  Duration: {T_END} time units, solver=scipy/RK45")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        r_coords=r,
        theta_coords=theta,
        phi_coords=phi_coord,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])

    initial_max = float(np.max(np.abs(initial[0].data)))
    final_max = float(np.max(np.abs(final[0].data)))

    print(f"  Initial max|phi| = {initial_max:.4f}")
    print(f"  Final   max|phi| = {final_max:.4f}")

    # Check for 1/r decay: amplitude should decrease as wave spreads
    if final_max < initial_max:
        print("  Amplitude decay observed (expected for 3D spherical spreading)")
    else:
        print("  Note: Amplitude did not decay (may be due to boundary reflections)")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    grid = result.grid

    initial = cast("FieldCollection", storage[0])
    final = cast("FieldCollection", storage[-1])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Take a phi=const slice (equatorial plane in phi)
    phi_mid = GRID_NPHI // 2

    # Initial field: r-theta slice
    ax = axes[0, 0]
    vmax = PULSE_AMPLITUDE
    r_1d = cast("np.ndarray", grid.cell_coords[:, 0, 0, 0])
    cast("np.ndarray", grid.cell_coords[0, :, 0, 1])
    im = ax.imshow(
        initial[0].data[:, :, phi_mid].T,
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
    ax.set_title("Initial phi (r-theta slice, t=0)")
    fig.colorbar(im, ax=ax, label="phi")

    # Final field: r-theta slice
    ax = axes[0, 1]
    final_vmax = max(float(np.max(np.abs(final[0].data[:, :, phi_mid]))), 0.01)
    im = ax.imshow(
        final[0].data[:, :, phi_mid].T,
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
    ax.set_title(f"Final phi (r-theta slice, t={T_END:.0f})")
    fig.colorbar(im, ax=ax, label="phi")

    # Radial profile at theta=pi/2 (equator), phi=0
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

    theta_mid = GRID_NTHETA // 2  # Near theta=pi/2
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(
            r_1d,
            snapshot[0].data[:, theta_mid, 0],
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.set_xlabel("r")
    ax.set_ylabel("phi")
    ax.set_title("Radial profile (theta=pi/2, phi=0)")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    # r * phi(r) to check 1/r decay → should flatten
    ax = axes[1, 1]
    for i, t_idx in enumerate(times_to_plot):
        snapshot = cast("FieldCollection", storage[t_idx])
        profile = snapshot[0].data[:, theta_mid, 0]
        ax.plot(
            r_1d,
            r_1d * profile,
            color=colors[i],
            label=f"t={time_values[t_idx]:.1f}",
            alpha=0.8,
        )
    ax.set_xlabel("r")
    ax.set_ylabel("r * phi(r)")
    ax.set_title("r * phi(r): should flatten for 1/r decay")
    ax.legend(fontsize=7)
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        "KG in Spherical Coordinates: "
        r"$ds^2 = -dt^2 + dr^2 + r^2 d\theta^2 + r^2\sin^2\theta\,d\phi^2$"
        f"\nMassless (m^2={MASS_SQUARED}), shell pulse at r={PULSE_RADIUS}",
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
    print("Spherical coordinate simulation complete!")
    print()
    print("Key observations:")
    print("  1. Metric ds^2 = -dt^2 + dr^2 + r^2 dtheta^2 + r^2 sin^2(theta) dphi^2")
    print("  2. Christoffel corrections: (2/r) gradient_x, cot(theta)/r^2 gradient_y")
    print("  3. Trigonometric coefficients: Cot[y[]], Csc[y[]]^2")
    print("  4. 3+1D: 3 spatial dimensions (r, theta, phi)")
    print("  5. All coefficients derived from metric by xAct pipeline")
    print("=" * 60)


if __name__ == "__main__":
    main()
