"""Symbolic computation layer for Lagrangian-to-PDE pipeline.

This package provides the Python-side interface for loading and processing
field equations derived symbolically from Lagrangians via Mathematica/xAct.
"""

from tidal.symbolic.json_loader import (
    AXIS_LETTERS,
    BoundaryCondition,
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)
from tidal.symbolic.pde_builder import (
    PDEFromSpec,
    build_pde_from_json,
    create_initial_state,
    register_operator,
)

__all__ = [
    "AXIS_LETTERS",
    "BoundaryCondition",
    "ComponentEquation",
    "ConstraintSolverConfig",
    "EquationSystem",
    "OperatorTerm",
    "PDEFromSpec",
    "build_pde_from_json",
    "create_initial_state",
    "load_equation_system",
    "register_operator",
]
