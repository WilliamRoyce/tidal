from __future__ import annotations

from typing import TYPE_CHECKING

from pde import CartesianGrid

if TYPE_CHECKING:
    from .config import GridConfig


def make_grid(cfg: GridConfig) -> CartesianGrid:
    """
    Create a CartesianGrid from a GridConfig.

    Parameters
    ----------
    cfg : GridConfig
        Configuration object describing the desired grid. Must provide valid
        'bounds' (per-dimension min/max pairs), 'shape' (number of cells or points
        per dimension), and 'periodic' (boolean or per-axis periodicity flags).

    Returns
    -------
    CartesianGrid
        A py-pde CartesianGrid instance constructed according to the configuration.
    """
    cfg.validate()
    return CartesianGrid(
        bounds=list(cfg.bounds), shape=list(cfg.shape), periodic=cfg.periodic
    )
