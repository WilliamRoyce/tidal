"""TIDAL solver package — grid, spatial operators, time integrators."""

from tidal.solver._exceptions import (
    KineticEvaluationError,
    SimulationDivergedError,
)
from tidal.solver._types import SolverResult
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.fields import FieldSet
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_evolution_matrices,  # type: ignore[reportPrivateUsage]
    find_independent_blocks,
    solve_modal,
    solve_modal_pass1,
)
from tidal.solver.operators import OPERATOR_REGISTRY, apply_operator, operator_min_dim
from tidal.solver.perturbative_driver import (
    PerturbativeResult,
    PerturbativeSolver,
)
from tidal.solver.rhs import RHSEvaluator
from tidal.solver.state import StateLayout
from tidal.solver.validation import (
    StabilityResult,
    check_cfl_stability,
    check_mass_sign,
    validate_field_references,
    validate_operator_dimensions,
)

__all__ = [
    "OPERATOR_REGISTRY",
    "CoefficientEvaluator",
    "FieldSet",
    "GridInfo",
    "KineticEvaluationError",
    "PerturbativeResult",
    "PerturbativeSolver",
    "RHSEvaluator",
    "SimulationDivergedError",
    "SolverResult",
    "StabilityResult",
    "StateLayout",
    "_build_evolution_matrices",
    "apply_operator",
    "check_cfl_stability",
    "check_mass_sign",
    "find_independent_blocks",
    "operator_min_dim",
    "solve_modal",
    "solve_modal_pass1",
    "validate_field_references",
    "validate_operator_dimensions",
]
