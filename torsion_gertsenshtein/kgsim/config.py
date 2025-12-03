from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from collections.abc import Sequence


@dataclass(frozen=True)
class KGParameters:
    """Physical parameters for Klein-Gordon."""

    mass: float = 1.0  # m in natural units


@dataclass(frozen=True)
class MultiFieldParams:
    """
    Parameters for N coupled KG fields.

    masses: length-N list of bare masses m_i >= 0
    coupling: NxN symmetric matrix for bilinear mixing in the potential
              V = 1/2 * sum_{ij} phi_i (M2_ij) phi_j, where
              M2 = diag(m_i^2) + coupling
    """

    masses: Sequence[float]  # m in natural units
    coupling: Sequence[Sequence[float]]  # NxN

    def __post_init__(self) -> None:
        """
        Validate that masses are non-negative and coupling is NxN symmetric.

        Raises
        ------
        ValueError
            If any mass is negative.
            If coupling is not square or does not match length of masses.
            If coupling is not symmetric.
        """
        n = len(self.masses)
        # Check all masses are non-negative
        for i, m in enumerate(self.masses):
            if m < 0:
                msg = f"mass at index {i} is negative: {m}"
                raise ValueError(msg)
        # Check coupling is n x n
        if len(self.coupling) != n:
            msg = f"coupling matrix has {len(self.coupling)} rows, expected {n}"
            raise ValueError(msg)
        for i, row in enumerate(self.coupling):
            if len(row) != n:
                msg = f"coupling matrix row {i} has length {len(row)}, expected {n}"
                raise ValueError(msg)
        # Check symmetry
        for i in range(n):
            for j in range(i + 1, n):
                if self.coupling[i][j] != self.coupling[j][i]:
                    msg = f"coupling matrix is not symmetric at ({i},{j}) and ({j},{i}): {self.coupling[i][j]} != {self.coupling[j][i]}"
                    raise ValueError(msg)


@dataclass(frozen=True)
class GridConfig:
    """Grid configuration for py-pde CartesianGrid."""

    dim: int = 1
    shape: Sequence[int] = (512,)
    bounds: Sequence[tuple[float, float]] = ((0.0, 100.0),)
    periodic: bool | Sequence[bool] = True

    def __post_init__(self) -> None:
        """
        Validate that the object's shape and bounds match the configured dimensionality.

        Raises
        ------
        ValueError
            If the length of `shape` or `bounds` does not equal `dim`. The exception message
            will indicate which invariant failed, for example:
            - "shape length {len(self.shape)} must equal dim {self.dim}"
            - "bounds length {len(self.bounds)} must equal dim {self.dim}"
        """
        if len(self.shape) != self.dim:
            msg = f"shape length {len(self.shape)} must equal dim {self.dim}"
            raise ValueError(msg)
        if len(self.bounds) != self.dim:
            msg = f"bounds length {len(self.bounds)} must equal dim {self.dim}"
            raise ValueError(msg)


@dataclass(frozen=True)
class SimulationConfig:
    """Time stepping and runtime options."""

    t_end: float = 100.0
    dt: float | None = None  # If None, let solver choose (adaptive)
    backend: Literal["auto", "numpy", "numba"] = (
        "auto"  # py-pde backend: "numpy" or "numba"
    )
    solver: Literal["scipy", "explicit"] = "scipy"  # "scipy" (adaptive) or "explicit"
    method: str = "RK45"  # for scipy solver
    data_dir: str | None = None
    save_every: float | None = None  # save snapshots every X time units
    progress: bool = True

    def __post_init__(self) -> None:
        """
        Validate configuration after initialization.

        Raises
        ------
        ValueError
            If dt is not None and <= 0.0.
            If t_end <= 0.0.
            If solver is not 'scipy' or 'explicit'.
            If backend is not 'auto', 'numpy', or 'numba'.
        """
        if self.dt is not None and self.dt <= 0.0:
            msg = "dt must be positive or None for adaptive stepping"
            raise ValueError(msg)
        if self.t_end <= 0.0:
            msg = "t_end must be positive"
            raise ValueError(msg)
        valid_solvers = {"scipy", "explicit"}
        if self.solver not in valid_solvers:
            msg = f"invalid solver: {self.solver!r}"
            raise ValueError(msg)
        valid_backends = {"auto", "numpy", "numba"}
        if self.backend not in valid_backends:
            msg = f"invalid backend: {self.backend!r}"
            raise ValueError(msg)
