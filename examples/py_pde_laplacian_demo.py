# examples/py_pde_demo.py

from __future__ import annotations

import pathlib

from torsion_gertsenshtein.plot_pgf import enable_pgf

enable_pgf("xelatex")  # or "pdflatex"/"lualatex"

import matplotlib.pyplot as plt
import numpy as np
from pde import (
    CartesianGrid,
    solve_laplace_equation,  # pyright: ignore[reportUnknownVariableType]
)


def main() -> None:
    """Run a py-pde demonstration solving the Laplace equation on a 2D grid.

    Sets up a 64x64 CartesianGrid on [0, 2 pi] x [0, 2 pi], applies simple
    sinusoidal boundary conditions, solves the Laplace equation using py-pde,
    generates and saves a PNG figure and prints the saved path.
    """
    # 2D grid 64x64 on [0,2pi]x[0,2pi]
    grid = CartesianGrid([(0, 2 * np.pi), (0, 2 * np.pi)], 64)
    bcs = {"x": {"value": "sin(y)"}, "y": {"value": "sin(x)"}}

    # Diffusion equation
    result = solve_laplace_equation(grid, bcs)

    # Save a figure
    fig, ax = plt.subplots(figsize=(4, 3))
    im = ax.imshow(result.data.T, origin="lower", extent=(0, 2 * np.pi, 0, 2 * np.pi))
    ax.set_title("py-pde: Solution to Laplace equation")
    fig.colorbar(im, ax=ax, shrink=0.8)

    pathlib.Path("outputs").mkdir(exist_ok=True, parents=True)
    out = "outputs/py_pde_demo.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
