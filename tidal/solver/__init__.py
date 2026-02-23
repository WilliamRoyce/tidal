"""TIDAL solver package — grid, spatial operators, time integrators."""

from tidal.solver._exceptions import SimulationDivergedError
from tidal.solver._types import SolverResult
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.fields import FieldSet
from tidal.solver.grid import GridInfo
from tidal.solver.operators import OPERATOR_REGISTRY, apply_operator, operator_min_dim
from tidal.solver.rhs import RHSEvaluator
from tidal.solver.state import StateLayout
from tidal.solver.validation import (
    check_cfl_stability,
    check_mass_sign,
    check_pointwise_mass_stability,
    validate_field_references,
    validate_operator_dimensions,
)

__all__ = [
    "OPERATOR_REGISTRY",
    "CoefficientEvaluator",
    "FieldSet",
    "GridInfo",
    "RHSEvaluator",
    "SimulationDivergedError",
    "SolverResult",
    "StateLayout",
    "apply_operator",
    "check_cfl_stability",
    "check_mass_sign",
    "check_pointwise_mass_stability",
    "operator_min_dim",
    "validate_field_references",
    "validate_operator_dimensions",
]
