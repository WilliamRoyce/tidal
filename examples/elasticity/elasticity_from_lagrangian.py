"""Linear elasticity simulation using Lagrangian-derived Navier-Cauchy equations.

This script demonstrates the full Lagrangian-to-PDE pipeline for vector PDEs:
1. The Navier-Cauchy equations were derived symbolically from the elasticity Lagrangian
2. The equations are loaded from JSON and converted to py-pde format
3. A simulation is run showing elastic wave propagation

The Lagrangian is:
    L = (1/2)*rho*|u_dot|^2 - (1/2)*lambda*(div u)^2 - mu*u_ij*u_ij

where u_ij = (1/2)(d_i u_j + d_j u_i) is the strain tensor.

The derived equations (Navier-Cauchy) with rho=lambda=mu=1 are:
    d²_t(ux) = 3∂²_x(ux) + 1∂²_y(ux) + 2∂_x∂_y(uy)
    d²_t(uy) = 1∂²_x(uy) + 3∂²_y(uy) + 2∂_x∂_y(ux)
"""

from __future__ import annotations

from pathlib import Path
from typing import cast

import numpy as np
from pde import CartesianGrid, FieldCollection, ScalarField

OUTPUT_FILENAME = "elasticity_result.png"

from torsion_gertsenshtein.symbolic import build_pde_from_json


def create_initial_displacement(grid: CartesianGrid) -> FieldCollection:
    """Create initial state with a localized displacement in ux.

    Parameters
    ----------
    grid : CartesianGrid
        The 2D simulation grid.

    Returns
    -------
    FieldCollection
        Initial state: [ux, pi_ux, uy, pi_uy]
        where pi_ux = d_t(ux) and pi_uy = d_t(uy)
    """
    x = cast("np.ndarray", grid.cell_coords[..., 0])
    y = cast("np.ndarray", grid.cell_coords[..., 1])

    # Center of the grid
    _lx, _ly = grid.shape
    x0 = grid.axes_bounds[0][1] / 2
    y0 = grid.axes_bounds[1][1] / 2

    # Gaussian bump in ux
    sigma = 1.0
    ux_init = np.exp(-((x - x0) ** 2 + (y - y0) ** 2) / (2 * sigma**2))

    # State layout: ux, d_t(ux), uy, d_t(uy)
    return FieldCollection(
        [
            ScalarField(grid, data=ux_init),  # ux_0 (horizontal displacement)
            ScalarField(grid, data=0.0),  # pi_ux (velocity)
            ScalarField(grid, data=0.0),  # uy_0 (vertical displacement)
            ScalarField(grid, data=0.0),  # pi_uy (velocity)
        ]
    )


def run_simulation(
    json_path: str | Path,
    grid_size: tuple[int, int] = (64, 64),
    domain_size: tuple[float, float] = (10.0, 10.0),
    t_final: float = 5.0,
    dt: float = 0.01,
) -> FieldCollection:
    """Run the elasticity simulation.

    Parameters
    ----------
    json_path : str or Path
        Path to the JSON equation specification.
    grid_size : tuple[int, int]
        Number of grid points in each direction.
    domain_size : tuple[float, float]
        Physical size of the domain in each direction.
    t_final : float
        Final simulation time.
    dt : float
        Time step for the solver.

    Returns
    -------
    FieldCollection
        Final state after simulation.
    """
    # Load PDE from JSON
    print(f"Loading PDE from {json_path}...")
    pde = build_pde_from_json(json_path)
    print(f"  Fields: {pde.spec.component_names}")
    print(f"  Number of components: {pde.n_components}")

    # Create 2D grid with periodic boundaries
    grid = CartesianGrid(
        [(0, domain_size[0]), (0, domain_size[1])],
        grid_size,
        periodic=True,
    )
    print(f"  Grid: {grid_size[0]}x{grid_size[1]}, domain {domain_size}")

    # Create initial state
    print("Creating initial displacement (Gaussian bump in ux)...")
    state = create_initial_displacement(grid)

    # Print initial state statistics
    print(f"  Initial ux max: {np.max(np.abs(state[0].data)):.4f}")
    print(f"  Initial uy max: {np.max(np.abs(state[2].data)):.4f}")

    # Run simulation
    print(f"Running simulation to t={t_final}...")
    result = pde.solve(state, t_range=t_final, dt=dt)

    # solve() returns FieldCollection for our usage
    assert isinstance(result, FieldCollection)
    return result


def main() -> None:
    """Run the elasticity simulation."""
    # Path to the JSON file
    json_path = Path(__file__).parent.parent / "data" / "navier_cauchy_2d.json"

    if not json_path.exists():
        print(f"JSON file not found: {json_path}")
        print("Run the Wolfram script first: wolframscript -file navier_cauchy.wls")
        return

    # Run simulation
    final_state = run_simulation(
        json_path=json_path,
        grid_size=(64, 64),
        domain_size=(10.0, 10.0),
        t_final=3.0,
        dt=0.005,
    )

    # Print final statistics
    print("\nFinal state:")
    print(f"  ux max: {np.max(np.abs(final_state[0].data)):.4f}")
    print(f"  pi_ux max: {np.max(np.abs(final_state[1].data)):.4f}")
    print(f"  uy max: {np.max(np.abs(final_state[2].data)):.4f}")
    print(f"  pi_uy max: {np.max(np.abs(final_state[3].data)):.4f}")

    # Optional: Save visualization
    try:
        import matplotlib.pyplot as plt  # noqa: PLC0415

        _fig, axes = plt.subplots(2, 2, figsize=(12, 10))

        # ux displacement
        im0 = axes[0, 0].imshow(
            final_state[0].data.T, origin="lower", cmap="RdBu", aspect="equal"
        )
        axes[0, 0].set_title("ux (horizontal displacement)")
        plt.colorbar(im0, ax=axes[0, 0])

        # ux velocity
        im1 = axes[0, 1].imshow(
            final_state[1].data.T, origin="lower", cmap="RdBu", aspect="equal"
        )
        axes[0, 1].set_title("∂ux/∂t (horizontal velocity)")
        plt.colorbar(im1, ax=axes[0, 1])

        # uy displacement
        im2 = axes[1, 0].imshow(
            final_state[2].data.T, origin="lower", cmap="RdBu", aspect="equal"
        )
        axes[1, 0].set_title("uy (vertical displacement)")
        plt.colorbar(im2, ax=axes[1, 0])

        # uy velocity
        im3 = axes[1, 1].imshow(
            final_state[3].data.T, origin="lower", cmap="RdBu", aspect="equal"
        )
        axes[1, 1].set_title("∂uy/∂t (vertical velocity)")
        plt.colorbar(im3, ax=axes[1, 1])

        plt.suptitle(
            "Navier-Cauchy Elastic Wave Propagation\n(derived from Lagrangian)"
        )
        plt.tight_layout()

        output_dir = Path(__file__).parent.parent.parent / "outputs"
        output_dir.mkdir(exist_ok=True)
        output_path = output_dir / OUTPUT_FILENAME
        plt.savefig(output_path, dpi=150)
        print(f"\nVisualization saved to: {output_path}")

    except ImportError:
        print("\nMatplotlib not available for visualization")


if __name__ == "__main__":
    main()
