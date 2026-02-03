"""Symbolic computation layer for Lagrangian-to-PDE pipeline.

This package provides the Python-side interface for loading and processing
field equations derived symbolically from Lagrangians via Mathematica/xAct.
"""

from torsion_gertsenshtein.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)
from torsion_gertsenshtein.symbolic.pde_builder import (
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
)

__all__ = [
    "ComponentEquation",
    "EquationSystem",
    "OperatorTerm",
    "PDEFromSpec",
    "build_pde_from_json",
    "create_initial_state",
    "load_equation_system",
]
