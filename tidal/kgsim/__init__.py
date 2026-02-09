"""Klein-Gordon equation simulation package."""

from tidal.kgsim.advanced_equations import (
    AnisotropicHigherOrderKGPDE,
    AnisotropicKGPDE,
    DirectionalKGPDE,
    HigherOrderKGPDE,
)
from tidal.kgsim.animation_builder import (
    AnimationBuilder,
    AnimationConfig,
)
from tidal.kgsim.config import (
    AnisotropicKGParameters,
    GridConfig,
    HigherOrderKGParameters,
    KGParameters,
    MultiFieldParams,
    SimulationConfig,
)
from tidal.kgsim.equations import (
    InhomogeneousKGPDE,
    KleinGordonPDE,
    make_coupled_kg_pde,
)
from tidal.kgsim.grids import make_grid
from tidal.kgsim.initial_conditions import (
    GaussianPulse,
    InitialCondition,
    RingPulse2D,
    gaussian_pulse,
    multi_gaussian,
    multi_gaussian_2d,
    plane_wave,
    ring_pulse_2d,
)
from tidal.kgsim.observers import total_energy_observer
from tidal.kgsim.profiles import (
    constant_field,
    step_region_1d,
)
from tidal.kgsim.runners import run, run_with_snapshots

__all__ = [
    "AnimationBuilder",
    "AnimationConfig",
    "AnisotropicHigherOrderKGPDE",
    "AnisotropicKGPDE",
    "AnisotropicKGParameters",
    "DirectionalKGPDE",
    "GaussianPulse",
    "GridConfig",
    "HigherOrderKGPDE",
    "HigherOrderKGParameters",
    "InhomogeneousKGPDE",
    "InitialCondition",
    "KGParameters",
    "KleinGordonPDE",
    "MultiFieldParams",
    "RingPulse2D",
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
    "run_with_snapshots",
    "step_region_1d",
    "total_energy_observer",
]
