"""Shared solver setup utilities.

Common initialization patterns extracted from cvode.py, ida.py,
scipy_solver.py, and leapfrog.py to eliminate code duplication.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.rhs import RHSEvaluator
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem


def warn_frozen_constraints(
    layout: StateLayout,
    solver_name: str,
) -> list[str]:
    """Detect constraint fields (time_order=0) and warn they are frozen.

    Returns the list of constraint field names.
    """
    constraint_fields = [
        s.field_name for s in layout.slots if s.kind == "field" and s.time_order == 0
    ]
    if constraint_fields:
        warnings.warn(
            f"{solver_name}: constraint fields {constraint_fields} frozen at initial "
            f"values (not evolved). Use --scheme ida for constraint systems.",
            stacklevel=3,
        )
    return constraint_fields


def build_rhs_evaluator(
    spec: EquationSystem,
    grid: GridInfo,
    parameters: dict[str, float] | None,
    bc: BCSpec | None,
) -> RHSEvaluator:
    """Build CoefficientEvaluator + RHSEvaluator from spec.

    Handles the lazy import to avoid circular dependencies.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
    from tidal.solver.rhs import RHSEvaluator as _RHSEvaluator  # noqa: PLC0415

    coeff_eval = CoefficientEvaluator(spec, grid, parameters or {})
    return _RHSEvaluator(spec, grid, coeff_eval, bc=bc)


def configure_linear_solver(  # noqa: PLR0913
    options: dict[str, Any],
    layout: StateLayout,
    spec: EquationSystem,
    grid: GridInfo,
    bc: BCSpec | None,
    *,
    parameters: dict[str, float] | None = None,
) -> None:
    """Choose and configure the linear solver based on system size.

    Mutates *options* in-place, adding ``linsolver`` (and ``sparsity``
    for the sparse tier, or ``jacfn``/``jactimes`` for constant-coefficient
    systems with an analytical Jacobian).

    Parameters
    ----------
    parameters : dict[str, float] or None
        Runtime parameter overrides.  When provided, enables the analytical
        Jacobian optimisation for constant-coefficient systems.
    """
    from tidal.solver._types import DENSE_THRESHOLD, SPARSE_THRESHOLD  # noqa: PLC0415

    # Analytical Jacobian for constant-coefficient systems
    if parameters is not None:
        from tidal.solver.analytical_jacobian import (  # noqa: PLC0415
            try_analytical_jacobian,
        )

        if try_analytical_jacobian(options, spec, layout, grid, bc, parameters):
            return

    # Existing tier system for non-constant or when parameters not provided
    n_state = layout.total_size
    if n_state <= DENSE_THRESHOLD:
        options["linsolver"] = "dense"
    elif n_state <= SPARSE_THRESHOLD:
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        pattern = build_jacobian_sparsity(spec, layout, grid, bc)
        options["linsolver"] = "sparse"
        options["sparsity"] = pattern
    else:
        options["linsolver"] = "gmres"
