"""Vector and tensor field simulation infrastructure.

This package provides utilities for simulating multi-component field systems
derived from Lagrangians. It builds on the symbolic package to provide
convenient initial condition creation and configuration.
"""

from torsion_gertsenshtein.vectorfield.config import (
    ComponentFieldParams,
)
from torsion_gertsenshtein.vectorfield.initial_conditions import (
    ComponentGaussianPulse,
    ComponentPlaneWave,
    create_gaussian_pulse_state,
)

__all__ = [
    "ComponentFieldParams",
    "ComponentGaussianPulse",
    "ComponentPlaneWave",
    "create_gaussian_pulse_state",
]
