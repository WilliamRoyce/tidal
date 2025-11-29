"""Klein-Gordon equation simulation package."""

from torsion_gertsenshtein.kgsim.config import (
    GridConfig,
    KGParameters,
    SimulationConfig,
)
from torsion_gertsenshtein.kgsim.equations import InhomogeneousKGPDE, KleinGordonPDE
from torsion_gertsenshtein.kgsim.grids import make_grid
from torsion_gertsenshtein.kgsim.initial_conditions import (
    gaussian_pulse,
    plane_wave,
    ring_pulse_2d,
)
from torsion_gertsenshtein.kgsim.observers import total_energy_observer
from torsion_gertsenshtein.kgsim.profiles import (
    constant_field,
    step_region_1d,
)
from torsion_gertsenshtein.kgsim.runners import run

__all__ = [
    "GridConfig",
    "InhomogeneousKGPDE",
    "KGParameters",
    "KleinGordonPDE",
    "SimulationConfig",
    "constant_field",
    "gaussian_pulse",
    "make_grid",
    "plane_wave",
    "ring_pulse_2d",
    "run",
    "step_region_1d",
    "total_energy_observer",
]
