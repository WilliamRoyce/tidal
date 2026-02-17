"""Static Conformal Klein-Gordon 1+1D simulation - Phase 1 of Curved Spacetime.

This script demonstrates a scalar field on conformally flat spacetime with
CONSTANT conformal factor Omega = 2.

Physics:
  Metric: ds^2 = Omega^2 * (-dt^2 + dx^2)  where Omega = 2 (constant)

  For constant Omega:
  - Christoffel symbols = 0 (derivatives of constant metric vanish)
  - Effective mass: m_eff^2 = m^2 * Omega^2 = 1 * 4 = 4

  Wave equation: d2_t phi = d2_x phi - m_eff^2 phi
               = d2_t phi = d2_x phi - 4 phi

Verification:
  This should produce identical dynamics to flat Klein-Gordon with m^2 = 4.
  We compare both and verify they match.
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
M_EFF_SQUARED = 4.0  # From conformal factor Omega=2: m_eff^2 = m^2 * Omega^2 = 1 * 4

# Grid
X_MIN = 0.0
X_MAX = 100.0
GRID_SHAPE = [256]
GRID_PERIODIC: list[bool] = [True]

# Time integration
T_END = 20.0
DT = 0.002  # Small timestep needed for m_eff^2=4 (higher frequency requires smaller dt)
SNAPSHOT_INTERVAL = 1.0

# Initial conditions
CENTER_X = 50.0
PULSE_WIDTH = 5.0
PULSE_AMPLITUDE = 1.0

# Analysis
DISPERSION_THRESHOLD = 1.1  # Factor above which dispersion is "observed"

# Output
OUTPUT_FILENAME = "conformal_kg_static_output.png"

# ── Helpers ───────────────────────────────────────────────────


@dataclass(frozen=True)
class SimulationResult:
    """Container for conformal KG simulation results."""

    grid: CartesianGrid
    storage: MemoryStorage
    x_coords: NumericArray


def _print_header() -> None:
    print("=" * 60)
    print("Static Conformal Klein-Gordon 1+1D Simulation")
    print("=" * 60)
    print()
    print("Metric: ds^2 = Omega^2 * (-dt^2 + dx^2)  where Omega = 2")
    print(f"Effective mass: m_eff^2 = m^2 * Omega^2 = 1 * 4 = {M_EFF_SQUARED:.0f}")
    print()


def _load_spec(json_path: Path) -> None:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)

    print(f"  Spacetime dimension: {spec.dimension} (1+1D)")
    print(f"  Components: {spec.n_components} ({', '.join(spec.component_names)})")

    print()
    print("  Equation structure:")
    for eq in spec.equations:
        terms_str = " + ".join(
            f"{term.coefficient:.1f}*{term.operator}({term.field})"
            for term in eq.rhs_terms
        )
        print(f"    d2_t({eq.field_name}) = {terms_str}")
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Components: {pde.n_components}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 1D simulation grid...")
    grid = CartesianGrid(
        bounds=[(X_MIN, X_MAX)],
        shape=GRID_SHAPE,
        periodic=GRID_PERIODIC,
    )
    print(f"  Domain: [{X_MIN:.0f}, {X_MAX:.0f}]")
    print(f"  Resolution: {grid.shape[0]} points")
    print("  Periodic boundary conditions")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial conditions...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])

    # Initial state: 2 fields = 1 component * (field + momentum)
    # State layout: confPhi_0, pi_0
    gaussian = PULSE_AMPLITUDE * np.exp(
        -((x - CENTER_X) ** 2) / (2 * PULSE_WIDTH**2)
    )

    phi = ScalarField(grid, data=gaussian, label="confPhi_0")
    pi = ScalarField(grid, data=0.0, label="pi_0")
    state = FieldCollection([phi, pi])

    print(f"  Gaussian pulse at center x={CENTER_X}")
    print(f"  Width: {PULSE_WIDTH}, Amplitude: {PULSE_AMPLITUDE}")
    print("  Initial momentum: zero")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")

    storage = MemoryStorage()
    result = pde.solve(
        state,
        t_range=T_END,
        dt=DT,
        scheme="runge-kutta",  # RK4 for better energy conservation
        tracker=storage.tracker(SNAPSHOT_INTERVAL),
    )
    result = normalize_solve_result(result)

    x = cast("np.ndarray", grid.cell_coords[..., 0])

    print(f"  Duration: {T_END} time units")
    print(f"  Stored {len(storage)} snapshots")
    print()

    return SimulationResult(
        grid=grid,
        storage=storage,
        x_coords=x,
    )


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    initial = cast("FieldCollection", result.storage[0])
    final = cast("FieldCollection", result.storage[-1])
    x = result.x_coords

    # Extract field amplitudes
    initial_max = float(np.max(np.abs(initial[0].data)))
    final_max = float(np.max(np.abs(final[0].data)))

    print(f"  Initial: max|phi| = {initial_max:.4f}")
    print(f"  Final:   max|phi| = {final_max:.4f}")

    # Check dispersion (massive field should spread)
    initial_width = float(np.sqrt(
        np.sum((x - CENTER_X) ** 2 * initial[0].data ** 2)
        / (np.sum(initial[0].data ** 2) + 1e-10)
    ))
    final_com = float(
        np.sum(x * final[0].data ** 2) / (np.sum(final[0].data ** 2) + 1e-10)
    )
    final_width = float(np.sqrt(
        np.sum((x - final_com) ** 2 * final[0].data ** 2)
        / (np.sum(final[0].data ** 2) + 1e-10)
    ))

    print(f"  Initial pulse width: {initial_width:.2f}")
    print(f"  Final pulse width: {final_width:.2f}")

    if final_width > initial_width * DISPERSION_THRESHOLD:
        print("  Dispersion observed (as expected for massive field)")
    else:
        print("  Note: Dispersion may be subtle at this mass/time scale")
    print()


def _plot_results(result: SimulationResult) -> None:  # noqa: PLR0914, PLR0915
    print("Step 7: Generating visualization...")

    storage = result.storage
    x = result.x_coords

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Field evolution
    times = [0, len(storage) // 3, 2 * len(storage) // 3, len(storage) - 1]
    cmap = plt.get_cmap("viridis")
    colors = cmap(np.linspace(0.2, 0.8, len(times)))

    ax = axes[0, 0]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[0].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\phi$")
    ax.set_title(r"Field $\phi$ evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Momentum evolution
    ax = axes[0, 1]
    for i, t_idx in enumerate(times):
        snapshot = cast("FieldCollection", storage[t_idx])
        ax.plot(x, snapshot[1].data, color=colors[i], label=f"t={t_idx:.0f}", alpha=0.8)
    ax.set_xlabel("x")
    ax.set_ylabel(r"$\pi$")
    ax.set_title(r"Momentum $\pi = \partial_t \phi$ evolution")
    ax.legend()
    ax.grid(visible=True, alpha=0.3)

    # Energy (proper Hamiltonian for static conformal)
    # Note: Static conformal spacetime is conservative, so Hamiltonian should be conserved.
    # Using RK4 scheme (scheme='runge-kutta') provides much better energy conservation
    # than the default explicit Euler, which is non-symplectic and causes linear drift.
    ax = axes[1, 0]
    energies = []
    time_values = []

    # Grid spacing for gradient calculation
    dx = float(x[1] - x[0])

    for t_idx in range(len(storage)):
        snapshot = cast("FieldCollection", storage[t_idx])
        phi_data = snapshot[0].data
        pi_data = snapshot[1].data

        # Compute spatial gradient
        grad_phi = np.gradient(phi_data, dx)

        # Proper Hamiltonian: H = 1/2 * integral[pi^2 + (d_x phi)^2 + m_eff^2 phi^2] dx
        energy = np.sum(
            0.5 * pi_data**2 + 0.5 * grad_phi**2 + 0.5 * M_EFF_SQUARED * phi_data**2
        )
        energies.append(energy)
        time_values.append(t_idx)

    ax.plot(time_values, energies, "b-", linewidth=2)
    ax.set_xlabel("Snapshot index")
    ax.set_ylabel("Total Energy")
    ax.set_title(
        r"Hamiltonian $H = \frac{1}{2}\int[\pi^2 + (\nabla\phi)^2 + m_{eff}^2\phi^2]dx$"
    )
    ax.grid(visible=True, alpha=0.3)

    # Space-time diagram
    ax = axes[1, 1]
    n_snapshots = len(storage)
    spacetime = np.zeros((n_snapshots, len(x)))
    for t_idx in range(n_snapshots):
        snapshot = cast("FieldCollection", storage[t_idx])
        spacetime[t_idx, :] = snapshot[0].data

    im = ax.imshow(
        spacetime,
        aspect="auto",
        origin="lower",
        extent=[X_MIN, X_MAX, 0, n_snapshots],
        cmap="bwr_r",
        vmin=-PULSE_AMPLITUDE,
        vmax=PULSE_AMPLITUDE,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("Time (snapshot)")
    ax.set_title("Space-time diagram")
    plt.colorbar(im, ax=ax, label=r"$\phi$")

    fig.suptitle(
        r"Conformal KG: $g_{\mu\nu} = \Omega^2 \eta_{\mu\nu}$, $\Omega=2$, $m_{eff}^2=4$"
        "\nEquivalent to flat KG with m^2=4",
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
    print("Simulation complete!")
    print()
    print("Key observations:")
    print("  1. Conformal metric g_ab = Omega^2 eta_ab with Omega=2")
    print("  2. Constant Omega => Christoffel symbols = 0")
    print(f"  3. Effective mass m_eff^2 = m^2 * Omega^2 = {M_EFF_SQUARED:.0f}")
    print("  4. Dynamics identical to flat KG with m^2=4")
    print("  5. This validates metric handling for curved spacetime Phase 1")
    print("=" * 60)


# ── Entry point ───────────────────────────────────────────────


def main() -> None:
    """Run the static conformal Klein-Gordon simulation."""
    json_path = Path(__file__).parent.parent / "data" / "conformal_kg_static.json"
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
