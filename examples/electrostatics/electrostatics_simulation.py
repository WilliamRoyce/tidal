"""Electrostatics example: Poisson equation via constraint solver.

Demonstrates the elliptic constraint solver (Issue #91) by solving
the Poisson equation for the electrostatic potential:

    nabla^2 phi = -rho

where rho is a Gaussian charge distribution centered in the domain.

The constraint solver (py-pde's solve_poisson_equation) is invoked
at each timestep during evolution_rate(). Since rho has no evolution
(frozen field), the potential phi is computed once and remains static.

Physics:
    L = 1/2 (nabla phi)^2 + rho * phi
    Euler-Lagrange:  nabla^2 phi + rho = 0
    Rearranged:      nabla^2 phi = -rho

Boundary conditions: Dirichlet phi=0 on all boundaries (grounded box).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde import PDEBase

    from torsion_gertsenshtein.symbolic.json_loader import EquationSystem

    NumericArray = NDArray[np.float64]


OUTPUT_FILENAME = "electrostatics_solution.png"
DOMAIN_HALF_SIZE = 5.0
GRID_POINTS_PER_AXIS = 64
CHARGE_WIDTH = 0.5  # Gaussian sigma
CHARGE_MAGNITUDE = 1.0


@dataclass(frozen=True)
class SimulationResult:
    """Container for electrostatics simulation results."""

    grid: CartesianGrid
    state: FieldCollection
    rho_data: NumericArray
    phi_data: NumericArray
    x: NumericArray
    y: NumericArray


def main() -> None:
    """Run the electrostatics Poisson equation example."""
    _print_header()
    json_path = Path(__file__).parent.parent / "data" / "electrostatics_2d.json"
    _load_spec(json_path)
    pde = _build_pde(json_path)
    grid = _create_grid()
    state, rho_data = _create_initial_state(grid)
    result = _solve_poisson(pde, grid, state, rho_data)
    _verify_solution(result)
    _plot_results(result)
    _print_footer()


def _print_header() -> None:
    print("=" * 60)
    print("Electrostatics: Poisson Equation via Constraint Solver")
    print("=" * 60)
    print()
    print("Physics: nabla^2 phi = -rho (Gauss's law)")
    print("  Lagrangian: L = 1/2 (nabla phi)^2 + rho * phi")
    print("  Boundary conditions: Dirichlet phi=0 (grounded box)")
    print()


def _load_spec(json_path: Path) -> EquationSystem:
    print("Step 1: Loading equation specification...")
    spec = load_equation_system(json_path)
    print(f"  Fields: {spec.component_names}")
    print(f"  Dimension: {spec.dimension}D ({spec.spatial_dimension}D spatial)")

    # Show constraint solver status
    for eq in spec.equations:
        order = eq.time_derivative_order
        solver_status = "enabled" if eq.constraint_solver.enabled else "disabled"
        print(
            f"    {eq.field_name}: time_order={order}, constraint_solver={solver_status}"
        )

    print()
    return spec


def _build_pde(json_path: Path) -> PDEBase:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(str(json_path))
    print(f"  PDE class: {type(pde).__name__}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up simulation grid...")
    grid = CartesianGrid(
        [(-DOMAIN_HALF_SIZE, DOMAIN_HALF_SIZE), (-DOMAIN_HALF_SIZE, DOMAIN_HALF_SIZE)],
        [GRID_POINTS_PER_AXIS, GRID_POINTS_PER_AXIS],
        periodic=False,
    )
    print(f"  Domain: [{-DOMAIN_HALF_SIZE}, {DOMAIN_HALF_SIZE}]^2")
    print(f"  Resolution: {GRID_POINTS_PER_AXIS}x{GRID_POINTS_PER_AXIS} points")
    print("  Boundary conditions: Dirichlet phi=0 (non-periodic)")
    print()
    return grid


def _create_initial_state(
    grid: CartesianGrid,
) -> tuple[FieldCollection, NumericArray]:
    print("Step 4: Creating charge distribution...")
    print(f"  Gaussian charge: sigma={CHARGE_WIDTH}, q={CHARGE_MAGNITUDE}")

    # Gaussian charge distribution: rho(x,y) = -q * exp(-(x^2+y^2)/(2*sigma^2))
    # Negative sign because equation is nabla^2 phi + rho = 0, so rho = -charge
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])
    rho_data = -CHARGE_MAGNITUDE * np.exp(-(x**2 + y**2) / (2 * CHARGE_WIDTH**2))

    print(f"  rho max: {np.max(np.abs(rho_data)):.4f}")

    # Fields: phi (potential, starts at zero, solved via constraint), rho (charge
    # density, frozen), pi_rho (momentum, zero evolution)
    state = FieldCollection(
        [
            ScalarField(grid, data=0.0, label="phi"),
            ScalarField(grid, data=rho_data, label="rho"),
            ScalarField(grid, data=0.0, label="pi_rho"),
        ]
    )
    print()
    return state, rho_data


def _solve_poisson(
    pde: PDEBase,
    grid: CartesianGrid,
    state: FieldCollection,
    rho_data: NumericArray,
) -> SimulationResult:
    print("Step 5: Solving Poisson equation via constraint solver...")

    # Calling evolution_rate() triggers the constraint solver in Pass 1
    # The Poisson equation is solved for phi given the source rho
    rates = pde.evolution_rate(state, t=0.0)

    phi_data = np.asarray(state[0].data, dtype=float)
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    print(f"  phi max: {np.max(phi_data):.6f}")
    print(f"  phi min: {np.min(phi_data):.6f}")
    print(
        f"  phi(0,0): {phi_data[GRID_POINTS_PER_AXIS // 2, GRID_POINTS_PER_AXIS // 2]:.6f}"
    )

    # Verify constraint evolution rate is zero
    if not np.allclose(rates[0].data, 0.0):
        msg = "Constraint field evolution rate should be zero"
        raise ValueError(msg)

    print()
    return SimulationResult(
        grid=grid,
        state=state,
        rho_data=rho_data,
        phi_data=phi_data,
        x=x,
        y=y,
    )


def _verify_solution(result: SimulationResult) -> None:
    print("Step 6: Verifying solution properties...")

    # Check boundary conditions (phi = 0 at boundaries for Dirichlet)
    phi = result.phi_data
    boundary_max = max(
        np.max(np.abs(phi[0, :])),  # Left edge
        np.max(np.abs(phi[-1, :])),  # Right edge
        np.max(np.abs(phi[:, 0])),  # Bottom edge
        np.max(np.abs(phi[:, -1])),  # Top edge
    )
    print(f"  Boundary max |phi|: {boundary_max:.6e} (should be ~0 for Dirichlet)")

    # Check that potential is negative (positive charge with grounded boundaries)
    print(
        f"  Potential sign: {'negative' if np.max(phi) < 0 else 'positive'} (expected: negative)"
    )

    # The potential should be minimum (most negative) at the charge center
    center_idx = GRID_POINTS_PER_AXIS // 2
    phi_center = phi[center_idx, center_idx]
    phi_min = np.min(phi)
    print(f"  Potential at center: {phi_center:.4f}")
    print(f"  Minimum potential: {phi_min:.4f}")

    print()


def _plot_results(result: SimulationResult) -> None:
    print("Step 7: Generating visualization...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    fig.suptitle(
        "Electrostatics: Poisson Equation nabla^2(phi) = -rho\n"
        "Constraint Solver Integration (Issue #91)",
        fontsize=14,
    )

    # Charge distribution
    ax = axes[0]
    im0 = ax.imshow(
        result.rho_data.T,
        origin="lower",
        extent=[
            -DOMAIN_HALF_SIZE,
            DOMAIN_HALF_SIZE,
            -DOMAIN_HALF_SIZE,
            DOMAIN_HALF_SIZE,
        ],
        cmap="RdBu_r",
    )
    ax.set_title("Charge density rho(x,y)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im0, ax=ax, label="rho")

    # Potential
    ax = axes[1]
    im1 = ax.imshow(
        result.phi_data.T,
        origin="lower",
        extent=[
            -DOMAIN_HALF_SIZE,
            DOMAIN_HALF_SIZE,
            -DOMAIN_HALF_SIZE,
            DOMAIN_HALF_SIZE,
        ],
        cmap="viridis",
    )
    ax.set_title("Potential φ(x,y)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    fig.colorbar(im1, ax=ax, label="φ")

    # Radial profile
    ax = axes[2]
    r = np.sqrt(
        result.x[GRID_POINTS_PER_AXIS // 2, :] ** 2
        + result.y[GRID_POINTS_PER_AXIS // 2, :] ** 2
    )
    phi_radial = result.phi_data[GRID_POINTS_PER_AXIS // 2, :]
    sort_idx = np.argsort(r)
    ax.plot(r[sort_idx], phi_radial[sort_idx], "b-", linewidth=2, label="φ(r, θ=0)")
    ax.set_xlabel("r")
    ax.set_ylabel("φ(r)")
    ax.set_title("Radial Profile (y=0 slice)")
    ax.legend()
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
    print("Electrostatics simulation complete!")
    print()
    print("Key observations:")
    print("  1. Constraint equation (time_order=0) solved via Poisson solver")
    print("  2. Boundary conditions: Dirichlet φ=0 (grounded box)")
    print("  3. Potential is negative (positive charge, grounded boundaries)")
    print("  4. Solution computed in a single evolution_rate() call")
    print("=" * 60)


if __name__ == "__main__":
    main()
