"""Klein-Gordon equation simulation package."""

from torsion_gertsenshtein.kgsim.advanced_equations import (
    AnisotropicHigherOrderKGPDE,
    AnisotropicKGPDE,
    DirectionalKGPDE,
    HigherOrderKGPDE,
)
from torsion_gertsenshtein.kgsim.config import (
    AnisotropicKGParameters,
    GridConfig,
    HigherOrderKGParameters,
    KGParameters,
    MultiFieldParams,
    SimulationConfig,
)
from torsion_gertsenshtein.kgsim.equations import (
    InhomogeneousKGPDE,
    KleinGordonPDE,
    make_coupled_kg_pde,
)
from torsion_gertsenshtein.kgsim.grids import make_grid
from torsion_gertsenshtein.kgsim.initial_conditions import (
    gaussian_pulse,
    multi_gaussian,
    multi_gaussian_2d,
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
    "AnisotropicHigherOrderKGPDE",
    "AnisotropicKGPDE",
    "AnisotropicKGParameters",
    "DirectionalKGPDE",
    "GridConfig",
    "HigherOrderKGPDE",
    "HigherOrderKGParameters",
    "InhomogeneousKGPDE",
    "KGParameters",
    "KleinGordonPDE",
    "MultiFieldParams",
    "SimulationConfig",
    "constant_field",
    "gaussian_pulse",
    "make_coupled_kg_pde",
    "make_grid",
    "multi_gaussian",
    "multi_gaussian_2d",
    "plane_wave",
    "ring_pulse_2d",
    "run",
    "step_region_1d",
    "total_energy_observer",
]
