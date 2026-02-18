"""Linear elasticity simulation using Lagrangian-derived Navier-Cauchy equations.

This script demonstrates the full Lagrangian-to-PDE pipeline for vector PDEs:
1. The Navier-Cauchy equations were derived symbolically from the elasticity Lagrangian
2. The equations are loaded from JSON and converted to py-pde format
3. A simulation is run showing elastic wave propagation

The Lagrangian is:
    L = (1/2)*rho*|u_dot|^2 - (1/2)*lambda*(div u)^2 - mu*u_ij*u_ij

where u_ij = (1/2)(d_i u_j + d_j u_i) is the strain tensor.

The derived equations (Navier-Cauchy) with rho=lambda=mu=1 are:
    d2_t(ux) = 3 d2_x(ux) + 1 d2_y(ux) + 2 d_x d_y(uy)
    d2_t(uy) = 1 d2_x(uy) + 3 d2_y(uy) + 2 d_x d_y(ux)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, cast

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField

from tidal.symbolic import build_pde_from_json, load_equation_system

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.symbolic.pde_builder import PDEFromSpec

    NumericArray = NDArray[np.float64]

# -- Configuration --------------------------------------------------------

# Grid
GRID_SHAPE = (64, 64)
DOMAIN_SIZE = (10.0, 10.0)

# Time integration
T_FINAL = 3.0
DT = 0.005

# Initial conditions
GAUSSIAN_SIGMA = 1.0

# Output
OUTPUT_FILENAME = "elasticity_result.png"

# -- Helpers ---------------------------------------------------------------


@dataclass(frozen=True)
class SimulationResult:
    """Container for elasticity simulation results."""

    grid: CartesianGrid
    final_state: FieldCollection


def _print_header() -> None:
    print("=" * 60)
    print("Navier-Cauchy Elastic Wave Propagation")
    print("=" * 60)
    print()
    print("Lagrangian: L = (1/2) rho |u_dot|^2")
    print("              - (1/2) lambda (div u)^2 - mu u_ij u_ij")
    print("  with rho = lambda = mu = 1")
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
        operators = [t.operator for t in eq.rhs_terms]
        print(f"    d2_t({eq.field_name}): {n_terms} terms")
        print(f"      Operators: {operators}")
    print()


def _build_pde(json_path: Path) -> PDEFromSpec:
    print("Step 2: Building PDE from specification...")
    pde = build_pde_from_json(json_path)
    print(f"  PDE class: {type(pde).__name__}")
    print(f"  Fields: {pde.spec.component_names}")
    print(f"  Number of components: {pde.n_components}")
    print()
    return pde


def _create_grid() -> CartesianGrid:
    print("Step 3: Setting up 2D grid...")
    grid = CartesianGrid(
        [(0, DOMAIN_SIZE[0]), (0, DOMAIN_SIZE[1])],
        GRID_SHAPE,
        periodic=True,
    )
    print(f"  Grid: {GRID_SHAPE[0]}x{GRID_SHAPE[1]}, domain {DOMAIN_SIZE}")
    print("  Boundary: periodic")
    print()
    return grid


def _create_initial_state(grid: CartesianGrid) -> FieldCollection:
    print("Step 4: Creating initial displacement (Gaussian bump in ux)...")

    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    # Center of the grid
    x0 = grid.axes_bounds[0][1] / 2
    y0 = grid.axes_bounds[1][1] / 2

    # Gaussian bump in ux
    ux_init = np.exp(
        -((x - x0) ** 2 + (y - y0) ** 2) / (2 * GAUSSIAN_SIGMA**2)
    )

    # State layout: ux, d_t(ux), uy, d_t(uy)
    state = FieldCollection(
        [
            ScalarField(grid, data=ux_init),  # ux_0 (horizontal displacement)
            ScalarField(grid, data=0.0),  # pi_ux (velocity)
            ScalarField(grid, data=0.0),  # uy_0 (vertical displacement)
            ScalarField(grid, data=0.0),  # pi_uy (velocity)
        ]
    )

    print(f"  Gaussian sigma: {GAUSSIAN_SIGMA}")
    print(f"  Initial ux max: {np.max(np.abs(state[0].data)):.4f}")
    print(f"  Initial uy max: {np.max(np.abs(state[2].data)):.4f}")
    print()
    return state


def _run_simulation(
    pde: PDEFromSpec, grid: CartesianGrid, state: FieldCollection
) -> SimulationResult:
    print("Step 5: Running simulation...")

    result = pde.solve(state, t_range=T_FINAL, dt=DT)
    assert isinstance(result, FieldCollection)

    print(f"  Duration: {T_FINAL} time units, dt={DT}")
    print()

    return SimulationResult(grid=grid, final_state=result)


def _analyze_results(result: SimulationResult) -> None:
    print("Step 6: Analyzing results...")

    final = result.final_state
    print(f"  ux max:    {np.max(np.abs(final[0].data)):.4f}")
    print(f"  pi_ux max: {np.max(np.abs(final[1].data)):.4f}")
    print(f"  uy max:    {np.max(np.abs(final[2].data)):.4f}")
    print(f"  pi_uy max: {np.max(np.abs(final[3].data)):.4f}")
    print()


def _plot_results(result: SimulationResult) -> None:
    print("Step 7: Generating visualization...")

    final = result.final_state

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # ux displacement
    im0 = axes[0, 0].imshow(
        final[0].data.T, origin="lower", cmap="bwr", aspect="equal"
    )
    axes[0, 0].set_title("ux (horizontal displacement)")
    fig.colorbar(im0, ax=axes[0, 0])

    # ux velocity
    im1 = axes[0, 1].imshow(
        final[1].data.T, origin="lower", cmap="bwr", aspect="equal"
    )
    axes[0, 1].set_title("d_t(ux) (horizontal velocity)")
    fig.colorbar(im1, ax=axes[0, 1])

    # uy displacement
    im2 = axes[1, 0].imshow(
        final[2].data.T, origin="lower", cmap="bwr", aspect="equal"
    )
    axes[1, 0].set_title("uy (vertical displacement)")
    fig.colorbar(im2, ax=axes[1, 0])

    # uy velocity
    im3 = axes[1, 1].imshow(
        final[3].data.T, origin="lower", cmap="bwr", aspect="equal"
    )
    axes[1, 1].set_title("d_t(uy) (vertical velocity)")
    fig.colorbar(im3, ax=axes[1, 1])

    fig.suptitle(
        "Navier-Cauchy Elastic Wave Propagation\n(derived from Lagrangian)"
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
    print("Elasticity simulation complete!")
    print()
    print("Key observations:")
    print("  1. Navier-Cauchy equations derived from elasticity Lagrangian")
    print("  2. Cross-derivative coupling (d_x d_y) between ux and uy")
    print("  3. Anisotropic Laplacians: 3 along, 1 transverse (lambda+2mu vs mu)")
    print("  4. Elastic wave splitting into P-wave and S-wave components")
    print("=" * 60)


# -- Entry point -----------------------------------------------------------


def main() -> None:
    """Run the elasticity simulation."""
    json_path = Path(__file__).parent.parent / "data" / "navier_cauchy_2d.json"

    if not json_path.exists():
        print(f"JSON file not found: {json_path}")
        print("Run: tidal derive theory.toml (in examples/elasticity/)")
        return

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
