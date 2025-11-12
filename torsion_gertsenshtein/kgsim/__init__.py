"""Klein-Gordon equation simulation package."""

from torsion_gertsenshtein.kgsim.config import (
    GridConfig,
    KGParameters,
    SimulationConfig,
)
from torsion_gertsenshtein.kgsim.equations import KleinGordonPDE
from torsion_gertsenshtein.kgsim.grids import make_grid
from torsion_gertsenshtein.kgsim.initial_conditions import (
    gaussian_pulse,
    plane_wave,
    ring_pulse_2d,
)
from torsion_gertsenshtein.kgsim.observers import total_energy_observer
from torsion_gertsenshtein.kgsim.runners import run

__all__ = [
    "GridConfig",
    "KGParameters",
    "KleinGordonPDE",
    "SimulationConfig",
    "gaussian_pulse",
    "make_grid",
    "plane_wave",
    "ring_pulse_2d",
    "run",
    "total_energy_observer",
]
