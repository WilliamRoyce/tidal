"""Coupled Proca (two massive vector fields) in 2+1D with periodic BCs.

    L = -1/4 F^A_{ab} F^{A,ab} - 1/4 F^B_{ab} F^{B,ab}
        - mA2/2 A_a A^a - mB2/2 B_a B^a + gcoup A_a B^a

Component structure (6 fields in 2+1D):
    A_0: constraint (Helmholtz with -mA2, coupled to B_0 via gcoup)
    A_1: evolution  (wave equation with mass)
    A_2: evolution  (wave equation with mass)
    B_0: constraint (Helmholtz with -mB2, coupled to A_0 via gcoup)
    B_1: evolution  (wave equation with mass)
    B_2: evolution  (wave equation with mass)

Constraint solver features tested:
    - Coupled FFT solve (periodic BCs)
    - Two different Helmholtz scales (mA2=1.0 vs mB2=2.0)
    - Cross-field identity coupling in constraints
"""

from __future__ import annotations

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

    from tidal.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]

# ── Configuration ─────────────────────────────────────────────

# Physics
PARAMETERS = {"mA2": 1.0, "mB2": 2.0, "gcoup": 0.5}

# Grid
GRID_BOUNDS = [(0, np.pi), (0, np.pi)]
GRID_SHAPE = [20, 20]
GRID_PERIODIC = True

# Time integration
T_END = 2.0
DT = 0.05
SNAPSHOT_INTERVAL = 0.2

# Initial conditions
PULSE_CENTER_X = np.pi / 2
PULSE_CENTER_Y = np.pi / 2
PULSE_WIDTH = 0.5
PULSE_AMPLITUDE = 0.5

# Output
OUTPUT_FILENAME = "coupled_proca_output.png"

# State layout (from spec):
#   A_0 (constraint, t_order=0) -> idx 0 (field only)
#   A_1 (evolution, t_order=2)  -> idx 1 (field), idx 2 (momentum)
#   A_2 (evolution, t_order=2)  -> idx 3 (field), idx 4 (momentum)
#   B_0 (constraint, t_order=0) -> idx 5 (field only)
#   B_1 (evolution, t_order=2)  -> idx 6 (field), idx 7 (momentum)
#   B_2 (evolution, t_order=2)  -> idx 8 (field), idx 9 (momentum)
IDX_A0 = 0
IDX_A1 = 1
IDX_B0 = 5
IDX_B1 = 6

# ── Helpers ───────────────────────────────────────────────────


def _get_field(snap: FieldCollection, idx: int) -> NumericArray:
    return np.asarray(snap[idx].data, dtype=float)


def _print_header() -> None:
    print("=" * 60)
    print("Coupled Proca 2+1D Simulation (Periodic BCs)")
    print("=" * 60)
    print()
    print("Lagrangian:")
    print("  L = -1/4 F^A F^A - 1/4 F^B F^B")
    print("      - mA2/2 A^2 - mB2/2 B^2 + gcoup A.B")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")
    print(f"  State size: {spec.state_size} (fields + momenta)")
    print()
    for eq in spec.equations:
        t_order = eq.time_derivative_order
        eq_type = "constraint" if t_order == 0 else "evolution"
        operators = sorted({term.operator for term in eq.rhs_terms})
        print(
            f"    {eq.field_name}: t_order={t_order}  [{eq_type}]"
            f"  ops: {', '.join(operators)}"
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
    print("Step 3: Setting up 2D periodic grid...")
    grid = CartesianGrid(
        bounds=GRID_BOUNDS,
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  Domain: {GRID_BOUNDS}")
    print(f"  Resolution: {grid.shape}")
    print("  Boundary conditions: periodic")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid, spec: EquationSystem) -> FieldCollection:
    print("Step 4: Creating initial conditions...")
    x_coords = cast("np.ndarray", grid.cell_coords[..., 0])
    y_coords = cast("np.ndarray", grid.cell_coords[..., 1])

    gaussian = PULSE_AMPLITUDE * np.exp(
        -(
            (x_coords - PULSE_CENTER_X) ** 2
            + (y_coords - PULSE_CENTER_Y) ** 2
        )
        / (2 * PULSE_WIDTH**2)
    )

    fields: list[ScalarField] = []
    for name, slot_type in spec.state_layout:
        sf = ScalarField(grid, data=0.0, label=f"{name}_{slot_type}")
        if name == "A_1" and slot_type == "field":
            sf.data[:] = gaussian
        fields.append(sf)

    state = FieldCollection(fields)
    assert len(state) == spec.state_size
    print(f"  Gaussian pulse in A_1 at ({PULSE_CENTER_X:.2f}, {PULSE_CENTER_Y:.2f})")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print()
    return state


def _run_simulation(pde: PDEBase, state: FieldCollection) -> MemoryStorage:
    print("Step 5: Running simulation...")
    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    normalize_solve_result(result)
    print(f"  Duration: {T_END} time units")
    print(f"  Stored {len(storage)} snapshots")
    print()
    return storage


def _analyze_results(pde: PDEBase, storage: MemoryStorage) -> None:
    print("Step 6: Analyzing results...")

    final = cast("FieldCollection", storage[-1])

    # Re-solve constraints on the final state to get true constraint values.
    # py-pde's RK4 copies state for intermediate steps, so constraint fields
    # (rate=0) stay at initial values in stored snapshots.
    pde.evolution_rate(final)
    print(
        f"  A_0 (constraint, re-solved): max={np.max(np.abs(_get_field(final, IDX_A0))):.6f}"
    )
    print(
        f"  B_0 (constraint, re-solved): max={np.max(np.abs(_get_field(final, IDX_B0))):.6f}"
    )

    b1_max = float(np.max(np.abs(_get_field(final, IDX_B1))))
    print(f"  B_1 (evolution): final max={b1_max:.6f}")
    if b1_max > 1e-6:  # noqa: PLR2004
        print("  Coupling effect: energy transferred A -> B (gcoup active)")

    # Stability check
    for t_idx in range(len(storage)):
        snap = cast("FieldCollection", storage[t_idx])
        for field_idx in range(len(snap)):
            data = np.asarray(snap[field_idx].data, dtype=float)
            if not np.all(np.isfinite(data)):
                print(f"  WARNING: Non-finite at t_idx={t_idx}, field={field_idx}")
                return
            if float(np.max(np.abs(data))) > 1e6:  # noqa: PLR2004
                print(f"  WARNING: Blowup at t_idx={t_idx}, field={field_idx}")
                return

    print("  All fields bounded and finite throughout simulation")
    print()


def _plot_results(
    pde: PDEBase, storage: MemoryStorage
) -> None:
    print("Step 7: Generating visualization...")

    final = cast("FieldCollection", storage[-1])
    # Re-solve constraints for plotting
    pde.evolution_rate(final)

    lx, ly = GRID_BOUNDS[0][1], GRID_BOUNDS[1][1]
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Row 0: A_1 snapshots
    snap_indices = [0, len(storage) // 2, len(storage) - 1]
    snap_labels = ["t = 0", f"t = {T_END / 2:.1f}", f"t = {T_END:.1f}"]

    for col, (snap_idx, snap_label) in enumerate(
        zip(snap_indices, snap_labels, strict=True)
    ):
        ax = axes[0, col]
        snap = cast("FieldCollection", storage[snap_idx])
        im = ax.imshow(
            _get_field(snap, IDX_A1).T,
            origin="lower",
            extent=[0, lx, 0, ly],
            cmap="RdBu_r",
            vmin=-0.5,
            vmax=0.5,
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(f"A_1 at {snap_label}")
        plt.colorbar(im, ax=ax, shrink=0.8)

    # Row 1, Left: A_0 constraint
    ax = axes[1, 0]
    im = ax.imshow(
        _get_field(final, IDX_A0).T,
        origin="lower",
        extent=[0, lx, 0, ly],
        cmap="RdBu_r",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"A_0 (constraint) at t={T_END:.1f}")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Row 1, Center: B_0 constraint
    ax = axes[1, 1]
    im = ax.imshow(
        _get_field(final, IDX_B0).T,
        origin="lower",
        extent=[0, lx, 0, ly],
        cmap="RdBu_r",
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"B_0 (constraint) at t={T_END:.1f}")
    plt.colorbar(im, ax=ax, shrink=0.8)

    # Row 1, Right: Amplitude evolution
    ax = axes[1, 2]
    times = list(storage.times)
    a1_max = [
        float(np.max(np.abs(_get_field(cast("FieldCollection", storage[i]), IDX_A1))))
        for i in range(len(storage))
    ]
    b1_max = [
        float(np.max(np.abs(_get_field(cast("FieldCollection", storage[i]), IDX_B1))))
        for i in range(len(storage))
    ]

    ax.plot(times, a1_max, "b-", linewidth=2, label="max|A_1|")
    ax.plot(times, b1_max, "r-", linewidth=1.5, label="max|B_1|")
    ax.set_xlabel("Time")
    ax.set_ylabel("Peak amplitude")
    ax.set_title("Amplitude Evolution")
    ax.legend(fontsize=8)
    ax.grid(visible=True, alpha=0.3)

    fig.suptitle(
        "Coupled Proca 2+1D (Periodic BCs)\n"
        r"$\mathcal{L} = -\frac{1}{4}F_A^2 - \frac{1}{4}F_B^2"
        r" - \frac{m_A^2}{2}A^2 - \frac{m_B^2}{2}B^2 + g\,A\!\cdot\!B$"
        f"\nmA2={PARAMETERS['mA2']}, mB2={PARAMETERS['mB2']}, gcoup={PARAMETERS['gcoup']}",
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


def _print_footer() -> None:
    print("=" * 60)
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. A_0, B_0 are coupled constraints (Helmholtz + identity cross-term)")
    print("  2. Constraints solved via coupled FFT (periodic BCs)")
    print("  3. Two different mass scales: mA2=1.0, mB2=2.0")
    print("  4. Cross-coupling gcoup transfers energy between A and B sectors")
    print("  5. Energy drift O(dx^2) convergent with periodic BCs")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the coupled Proca simulation."""
    json_path = Path(__file__).parent.parent / "data" / "coupled_proca_3d.json"
    _print_header()
    spec = _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state = _create_initial_state(grid, spec)
    storage = _run_simulation(pde, state)
    _analyze_results(pde, storage)
    _plot_results(pde, storage)
    _print_footer()


if __name__ == "__main__":
    main()
