"""GridInfo — lightweight Cartesian grid descriptor for TIDAL.

Replaces py-pde's CartesianGrid with a minimal frozen dataclass that provides
only what TIDAL needs: grid spacing for FD stencils, cell coordinates for
position-dependent coefficient evaluation, and periodicity flags for BC dispatch.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class GridInfo:
    """Immutable Cartesian grid descriptor.

    Parameters
    ----------
    bounds : tuple of (lo, hi) pairs
        Domain bounds per spatial axis, e.g. ``((0, 10), (0, 10))``.
    shape : tuple of int
        Number of grid cells per axis, e.g. ``(64, 64)``.
    periodic : tuple of bool
        Whether each axis is periodic, e.g. ``(True, True)``.

    Examples
    --------
    >>> g = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(True,))
    >>> g.ndim
    1
    >>> g.dx
    (0.15625,)
    >>> g.num_points
    64
    """

    bounds: tuple[tuple[float, float], ...]
    shape: tuple[int, ...]
    periodic: tuple[bool, ...]

    def __post_init__(self) -> None:
        if len(self.bounds) != len(self.shape):
            msg = (
                f"bounds has {len(self.bounds)} axes but shape has "
                f"{len(self.shape)} — must match"
            )
            raise ValueError(msg)
        if len(self.bounds) != len(self.periodic):
            msg = (
                f"bounds has {len(self.bounds)} axes but periodic has "
                f"{len(self.periodic)} — must match"
            )
            raise ValueError(msg)
        for i, (lo, hi) in enumerate(self.bounds):
            if hi <= lo:
                msg = f"bounds[{i}] = ({lo}, {hi}): upper must exceed lower"
                raise ValueError(msg)
        for i, n in enumerate(self.shape):
            if n < 2:  # noqa: PLR2004
                msg = f"shape[{i}] = {n}: need at least 2 cells per axis"
                raise ValueError(msg)

    @property
    def ndim(self) -> int:
        """Number of spatial dimensions."""
        return len(self.shape)

    @property
    def num_points(self) -> int:
        """Total number of grid cells (product of shape)."""
        return math.prod(self.shape)

    @property
    def dx(self) -> tuple[float, ...]:
        """Grid spacing per axis (cell-centred: dx = (hi - lo) / N)."""
        return tuple(
            (hi - lo) / n for (lo, hi), n in zip(self.bounds, self.shape, strict=False)
        )

    @property
    def discretization(self) -> tuple[float, ...]:
        """Alias for ``dx`` — matches py-pde CartesianGrid API."""
        return self.dx

    @property
    def axes_bounds(self) -> tuple[tuple[float, float], ...]:
        """Alias for ``bounds`` — matches py-pde CartesianGrid API."""
        return self.bounds

    @property
    def num_axes(self) -> int:
        """Alias for ``ndim`` — matches py-pde CartesianGrid API."""
        return self.ndim

    def axes_coords(self, axis: int) -> np.ndarray:
        """1-D array of cell centres along *axis*.

        Returns shape ``(shape[axis],)`` — a single 1-D vector, not a
        broadcasted grid.  Useful for building per-axis coordinate arrays
        without the overhead of full meshgrid.

        Raises
        ------
        IndexError
            If *axis* is out of range.
        """
        if axis < 0 or axis >= self.ndim:
            msg = f"axis {axis} out of range for {self.ndim}-D grid"
            raise IndexError(msg)
        lo, hi = self.bounds[axis]
        n = self.shape[axis]
        d = (hi - lo) / n
        return np.linspace(lo + d / 2, hi - d / 2, n)

    @property
    def cell_coords(self) -> np.ndarray:
        """Cell-centre coordinates, shape ``(*shape, ndim)``.

        For a 1D grid with bounds (0, 10) and shape (4,), the cell centres
        are at [1.25, 3.75, 6.25, 8.75] (centred in each cell).

        Unlike py-pde's ``CartesianGrid.cell_coords`` (which returns a list
        of arrays with inconsistent shapes), this always returns a single
        ndarray with the coordinate dimension as the last axis.
        """
        grids = np.meshgrid(
            *(self.axes_coords(i) for i in range(self.ndim)),
            indexing="ij",
        )
        return np.stack(grids, axis=-1)

    def coord_arrays(self) -> tuple[np.ndarray, ...]:
        """Broadcasted coordinate arrays, one per axis.

        Returns ``ndim`` arrays, each of shape ``self.shape``, containing the
        coordinate value along that axis at every grid point.  Equivalent to
        ``np.meshgrid(*1d_axes, indexing="ij")`` but cached-friendly.

        This is the most efficient way to evaluate position-dependent
        coefficients that depend on a single coordinate (e.g. ``f(x)``).
        """
        return tuple(
            np.meshgrid(
                *(self.axes_coords(i) for i in range(self.ndim)),
                indexing="ij",
            )
        )
