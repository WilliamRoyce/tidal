"""TIDAL solver package — grid, spatial operators, time integrators."""

from tidal.solver.grid import GridInfo
from tidal.solver.operators import OPERATOR_REGISTRY, apply_operator
from tidal.solver.state import StateLayout

__all__ = ["OPERATOR_REGISTRY", "GridInfo", "StateLayout", "apply_operator"]
