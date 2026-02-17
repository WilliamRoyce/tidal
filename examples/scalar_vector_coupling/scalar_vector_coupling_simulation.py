"""Scalar-vector coupling 2+1D simulation from Lagrangian-derived equations.

This script demonstrates the complete pipeline for a mixed-rank coupled system:

    L = -1/2 (d_a phi)(d^a phi) - phim2/2 phi^2
        -1/4 F_ab F^ab - Am2/2 A_a A^a
        + kCS/2 epsilon^abc A_a d_b A_c
        + gSV phi (d_a A^a)

Equations of motion:
    d2_t(phi_0) = laplacian(phi) - phim2*phi - gSV*d_t(A_0) + gSV*div(A)
    A_0         = constraint (Gauss law with coupling)
    d2_t(A_1)   = d2_y(A_1) - d_xy(A_2) + CS terms - Am2*A_1 - gSV*d_x(phi)
    d2_t(A_2)   = d2_x(A_2) - d_xy(A_1) + CS terms - Am2*A_2 - gSV*d_y(phi)

This stress-tests the pipeline with:
- Mixed-rank cross-field coupling (scalar phi + vector A)
- Cross-field first_derivative_t and gradient operators
- Epsilon tensor + divergence coupling in the same Lagrangian
- 4 symbolic constants simultaneously (phim2, Am2, kCS, gSV)
- 4x4 mass/coupling matrices
- Constraint equation (A_0) with cross-field terms

The key physics: a Gaussian pulse in phi transfers energy to the vector
field A via the divergence coupling gSV*phi*div(A), while the Chern-Simons
term creates helical mixing between A components.
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
    from pde import PDEBase
    from pde.trackers.base import TrackerBase

    from tidal.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
PARAMETERS = {"phim2": 1.0, "Am2": 0.5, "kCS": 0.3, "gSV": 0.2}

# Grid
GRID_BOUNDS = [(0, 50), (0, 50)]
GRID_SHAPE = [64, 64]
GRID_PERIODIC = True

# Time integration
T_END = 10.0
DT = 0.01
SNAPSHOT_INTERVAL = 0.5

# Initial conditions
PULSE_CENTER_X = 25.0
PULSE_CENTER_Y = 25.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Output
OUTPUT_FILENAME = "scalar_vector_coupling_output.png"

# State vector indices (from JSON spec time-derivative orders):
#   phi_0 (2nd order) -> field=0, momentum=1
#   A_0   (constraint) -> field=2
#   A_1   (2nd order) -> field=3, momentum=4
#   A_2   (2nd order) -> field=5, momentum=6
IDX_PHI = 0
IDX_PI_PHI = 1
IDX_A0 = 2
IDX_A1 = 3
IDX_PI_A1 = 4
IDX_A2 = 5
IDX_PI_A2 = 6

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationData:
    """Container for simulation results on a 2D grid."""

    x: NumericArray
    y: NumericArray
    times: list[float]
    snapshots: list[FieldCollection]


def _print_header() -> None:
    print("=" * 60)
    print("Scalar-Vector Coupling 2+1D Simulation")
    print("=" * 60)
    print()
    print("Lagrangian:")
    print("  L = -1/2 (d phi)^2 - phim2/2 phi^2")
    print("      -1/4 F_ab F^ab - Am2/2 A^2")
    print("      + kCS/2 eps^abc A_a d_b A_c")
    print("      + gSV phi div(A)")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Lagrangian: {spec.metadata.get('lagrangian_expr', 'N/A')}")
    print(f"  Dimension: {spec.dimension} (2+1D spacetime)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    # Report cross-field coupling terms
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field != eq.field_name:
                print(
                    f"  Cross-field: {eq.field_name} <- "
                    f"{term.operator}({term.field})"
                )

    print()
    return spec


def _build_pde(json_path: Path) -> PDEBase:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path, parameters=PARAMETERS)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Parameters: {PARAMETERS}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 2D simulation grid...")
    grid = CartesianGrid(
        bounds=GRID_BOUNDS,
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  Domain: {GRID_BOUNDS}")
    print(f"  Resolution: {grid.shape}")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(
    grid: CartesianGrid, spec: EquationSystem
) -> FieldCollection:
    """Create initial conditions: Gaussian in phi, all A components zero.

    This setup demonstrates energy transfer from the scalar field to
    the vector field via the gSV divergence coupling.
    """
    print("Step 4: Creating initial conditions...")

    # State layout: phi_0, pi_phi, A_0, A_1, pi_A1, A_2, pi_A2
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    # Gaussian pulse in the scalar field phi
    gaussian = PULSE_AMPLITUDE * np.exp(
        -((x - PULSE_CENTER_X) ** 2 + (y - PULSE_CENTER_Y) ** 2)
        / (2 * PULSE_WIDTH**2)
    )

    # Build the 7-field state vector
    phi = ScalarField(grid, data=gaussian, label="phi_0")
    pi_phi = ScalarField(grid, data=0.0, label="pi_phi")
    a0 = ScalarField(grid, data=0.0, label="A_0")
    a1 = ScalarField(grid, data=0.0, label="A_1")
    pi_a1 = ScalarField(grid, data=0.0, label="pi_A1")
    a2 = ScalarField(grid, data=0.0, label="A_2")
    pi_a2 = ScalarField(grid, data=0.0, label="pi_A2")

    state = FieldCollection([phi, pi_phi, a0, a1, pi_a1, a2, pi_a2])

    assert len(state) == spec.state_size, (
        f"State size mismatch: {len(state)} != {spec.state_size}"
    )

    print(f"  Gaussian pulse in phi at ({PULSE_CENTER_X}, {PULSE_CENTER_Y})")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print("  All A components initialized to zero")
    print(f"  State vector: {len(state)} fields")
    print()
    return state


def _run_simulation(pde: PDEBase, initial_state: FieldCollection) -> MemoryStorage:
    print("Step 5: Running simulation...")
    print(f"  Duration: {T_END} time units")

    storage = MemoryStorage()
    tracker: TrackerBase = storage.tracker(interrupts=SNAPSHOT_INTERVAL)

    result = pde.solve(
        initial_state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=tracker,
    )
    result = normalize_solve_result(result)

    print(f"  Stored {len(storage)} snapshots")
    print()
    return storage


def _collect_simulation_data(
    storage: MemoryStorage, grid: CartesianGrid
) -> SimulationData:
    x = np.asarray(grid.cell_coords[..., 0], dtype=float)
    y = np.asarray(grid.cell_coords[..., 1], dtype=float)
    times = list(storage.times)
    snapshots = [cast("FieldCollection", storage[i]) for i in range(len(storage))]
    return SimulationData(x=x, y=y, times=times, snapshots=snapshots)


def _get_component(snapshot: FieldCollection, index: int) -> NumericArray:
    field = snapshot[index]
    if not isinstance(field, ScalarField):
        msg = f"Expected ScalarField at index {index}, got {type(field).__name__}"
        raise TypeError(msg)
    return np.asarray(field.data, dtype=float)


def _analyze_results(simulation: SimulationData) -> None:
    print("Step 6: Analyzing results...")

    # Initial amplitudes
    initial = simulation.snapshots[0]
    i_phi = float(np.max(np.abs(_get_component(initial, IDX_PHI))))
    i_a0 = float(np.max(np.abs(_get_component(initial, IDX_A0))))
    i_a1 = float(np.max(np.abs(_get_component(initial, IDX_A1))))
    i_a2 = float(np.max(np.abs(_get_component(initial, IDX_A2))))

    print(f"  Initial: max|phi| = {i_phi:.3f}, max|A_0| = {i_a0:.3f}, "
          f"max|A_1| = {i_a1:.3f}, max|A_2| = {i_a2:.3f}")

    # Final amplitudes
    final = simulation.snapshots[-1]
    f_phi = float(np.max(np.abs(_get_component(final, IDX_PHI))))
    f_a0 = float(np.max(np.abs(_get_component(final, IDX_A0))))
    f_a1 = float(np.max(np.abs(_get_component(final, IDX_A1))))
    f_a2 = float(np.max(np.abs(_get_component(final, IDX_A2))))

    print(f"  Final:   max|phi| = {f_phi:.3f}, max|A_0| = {f_a0:.3f}, "
          f"max|A_1| = {f_a1:.3f}, max|A_2| = {f_a2:.3f}")

    # Check for energy transfer: phi -> A
    if f_a1 > 0.01 or f_a2 > 0.01:  # noqa: PLR2004
        print("  Coupling effect: energy transferred phi -> A (gSV coupling active)")
    else:
        print("  Note: coupling may be weak at these parameter values")

    # Check constraint A_0 stays small
    a0_max_over_time = [
        float(np.max(np.abs(_get_component(s, IDX_A0))))
        for s in simulation.snapshots
    ]
    print(f"  Constraint A_0: max over all time = {max(a0_max_over_time):.4f}")

    print()


def _plot_results(simulation: SimulationData) -> None:  # noqa: PLR0915, PLR0914
    """Generate a 3x2 visualization of the simulation results.

    Row 0: Initial and final phi_0 (2D heatmaps)
    Row 1: Final A_1 and A_2 (2D heatmaps showing coupling effect)
    Row 2: Amplitude evolution and energy transfer over time
    """
    print("Step 7: Generating visualization...")

    fig, axes = plt.subplots(3, 2, figsize=(12, 14))
    fig.suptitle(
        "Scalar-Vector Coupling 2+1D\n"
        r"$\mathcal{L} = -\frac{1}{2}(\partial\phi)^2"
        r" - \frac{1}{4}F^2 - \frac{Am2}{2}A^2"
        r" + \frac{kCS}{2}\epsilon A \partial A"
        r" + gSV\,\phi\,\nabla\!\cdot\!A$",
        fontsize=13,
    )

    initial = simulation.snapshots[0]
    final = simulation.snapshots[-1]
    x_max = GRID_BOUNDS[0][1]
    y_max = GRID_BOUNDS[1][1]

    # --- Row 0: phi_0 initial and final ---

    phi_initial = _get_component(initial, IDX_PHI)
    phi_final = _get_component(final, IDX_PHI)
    vmax_phi = max(np.max(np.abs(phi_initial)), np.max(np.abs(phi_final)), 0.1)

    im00 = axes[0, 0].imshow(
        phi_initial.T,
        origin="lower",
        cmap="bwr",
        vmin=-vmax_phi,
        vmax=vmax_phi,
        extent=[0, x_max, 0, y_max],
    )
    axes[0, 0].set_title("Initial phi_0")
    axes[0, 0].set_xlabel("x")
    axes[0, 0].set_ylabel("y")
    plt.colorbar(im00, ax=axes[0, 0])

    im01 = axes[0, 1].imshow(
        phi_final.T,
        origin="lower",
        cmap="bwr",
        vmin=-vmax_phi,
        vmax=vmax_phi,
        extent=[0, x_max, 0, y_max],
    )
    axes[0, 1].set_title(f"Final phi_0 (t={simulation.times[-1]:.1f})")
    axes[0, 1].set_xlabel("x")
    axes[0, 1].set_ylabel("y")
    plt.colorbar(im01, ax=axes[0, 1])

    # --- Row 1: Final A_1 and A_2 (shows coupling-induced excitation) ---

    a1_final = _get_component(final, IDX_A1)
    a2_final = _get_component(final, IDX_A2)
    vmax_a = max(np.max(np.abs(a1_final)), np.max(np.abs(a2_final)), 0.1)

    im10 = axes[1, 0].imshow(
        a1_final.T,
        origin="lower",
        cmap="bwr",
        vmin=-vmax_a,
        vmax=vmax_a,
        extent=[0, x_max, 0, y_max],
    )
    axes[1, 0].set_title(f"Final A_1 (t={simulation.times[-1]:.1f})")
    axes[1, 0].set_xlabel("x")
    axes[1, 0].set_ylabel("y")
    plt.colorbar(im10, ax=axes[1, 0])

    im11 = axes[1, 1].imshow(
        a2_final.T,
        origin="lower",
        cmap="bwr",
        vmin=-vmax_a,
        vmax=vmax_a,
        extent=[0, x_max, 0, y_max],
    )
    axes[1, 1].set_title(f"Final A_2 (t={simulation.times[-1]:.1f})")
    axes[1, 1].set_xlabel("x")
    axes[1, 1].set_ylabel("y")
    plt.colorbar(im11, ax=axes[1, 1])

    # --- Row 2: Time evolution ---

    # Peak amplitude evolution
    ax = axes[2, 0]
    phi_max = [
        float(np.max(np.abs(_get_component(s, IDX_PHI))))
        for s in simulation.snapshots
    ]
    a1_max = [
        float(np.max(np.abs(_get_component(s, IDX_A1))))
        for s in simulation.snapshots
    ]
    a2_max = [
        float(np.max(np.abs(_get_component(s, IDX_A2))))
        for s in simulation.snapshots
    ]
    a0_max = [
        float(np.max(np.abs(_get_component(s, IDX_A0))))
        for s in simulation.snapshots
    ]

    ax.plot(simulation.times, phi_max, "b-", linewidth=2, label=r"max|$\phi$|")
    ax.plot(simulation.times, a1_max, "r-", linewidth=1.5, label="max|A_1|")
    ax.plot(simulation.times, a2_max, "g-", linewidth=1.5, label="max|A_2|")
    ax.plot(simulation.times, a0_max, "k--", linewidth=1, label="max|A_0| (constraint)")
    ax.set_xlabel("Time")
    ax.set_ylabel("Peak amplitude")
    ax.set_title("Amplitude Evolution")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)

    # Energy transfer (L2 norms)
    ax = axes[2, 1]
    phi_energy = [
        float(np.sum(_get_component(s, IDX_PHI) ** 2))
        for s in simulation.snapshots
    ]
    a1_energy = [
        float(np.sum(_get_component(s, IDX_A1) ** 2))
        for s in simulation.snapshots
    ]
    a2_energy = [
        float(np.sum(_get_component(s, IDX_A2) ** 2))
        for s in simulation.snapshots
    ]
    vector_energy = [a1 + a2 for a1, a2 in zip(a1_energy, a2_energy, strict=True)]

    ax.plot(simulation.times, phi_energy, "b-", linewidth=2, label=r"$\phi$ energy")
    ax.plot(simulation.times, vector_energy, "r-", linewidth=2, label="A energy")
    ax.plot(
        simulation.times,
        [p + v for p, v in zip(phi_energy, vector_energy, strict=True)],
        "k--",
        linewidth=1,
        label="Total",
    )
    ax.set_xlabel("Time")
    ax.set_ylabel(r"Energy ($L^2$ norm)")
    ax.set_title("Energy Transfer: Scalar to Vector")
    ax.legend(fontsize=9)
    ax.grid(visible=True, alpha=0.3)

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
    print("  1. Equations derived from Lagrangian with 4 interaction terms")
    print("  2. Cross-field coupling: phi <-> A via gradient and first_derivative_t")
    print("  3. Energy transfers from scalar phi to vector A via gSV coupling")
    print("  4. CS term creates helical mixing between A components")
    print("  5. A_0 constraint solved at each step (not time-evolved)")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the scalar-vector coupling simulation."""
    json_path = Path(__file__).parent.parent / "data" / "scalar_vector_coupling.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    initial_state = _create_initial_state(grid, spec)
    storage = _run_simulation(pde, initial_state)
    simulation = _collect_simulation_data(storage, grid)
    _analyze_results(simulation)
    _plot_results(simulation)
    _print_footer()


if __name__ == "__main__":
    main()
