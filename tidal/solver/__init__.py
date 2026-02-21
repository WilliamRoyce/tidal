"""TIDAL solver package — grid, spatial operators, time integrators."""

from tidal.solver.grid import GridInfo
from tidal.solver.operators import OPERATOR_REGISTRY, apply_operator

__all__ = ["OPERATOR_REGISTRY", "GridInfo", "apply_operator"]
