"""Symbolic computation layer for Lagrangian-to-PDE pipeline.

This package provides the Python-side interface for loading and processing
field equations derived symbolically from Lagrangians via Mathematica/xAct.
"""

from __future__ import annotations

from tidal.symbolic.json_loader import (
    AXIS_LETTERS,
    BoundaryCondition,
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
    load_equation_system,
)

__all__ = [
    "AXIS_LETTERS",
    "BoundaryCondition",
    "ComponentEquation",
    "ConstraintSolverConfig",
    "EquationSystem",
    "OperatorTerm",
    "load_equation_system",
]
