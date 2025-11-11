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
class GridConfig:
    """Grid configuration for py-pde CartesianGrid."""

    dim: int = 1
    shape: Sequence[int] = (512,)
    bounds: Sequence[tuple[float, float]] = ((0.0, 100.0),)
    periodic: bool | Sequence[bool] = True

    def validate(self) -> None:
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
    backend: str = "auto"  # py-pde backend: "numpy" or "numba"
    solver: Literal["scipy", "explicit"] = "scipy"  # "scipy" (adaptive) or "explicit"
    method: str = "RK45"  # for scipy solver
    data_dir: str | None = None
    save_every: float | None = None  # save snapshots every X time units
    progress: bool = True
