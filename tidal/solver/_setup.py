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

    # Warn about non-periodic coefficients with periodic BCs — these break
    # the integration-by-parts identity and cause O(1) energy non-conservation.
    if bc is not None:
        from tidal.solver.operators import is_periodic_bc  # noqa: PLC0415

        periodic = tuple(is_periodic_bc(b) for b in bc)
        coeff_eval.check_periodic_coefficient_continuity(periodic)

    return _RHSEvaluator(spec, grid, coeff_eval, bc=bc)


def configure_linear_solver(  # noqa: PLR0913
    options: dict[str, Any],
    layout: StateLayout,
    spec: EquationSystem,
    grid: GridInfo,
    bc: BCSpec | None,
    *,
    parameters: dict[str, float] | None = None,
    solver: str = "ida",
) -> None:
    """Choose and configure the linear solver based on system size.

    Mutates *options* in-place, adding ``linsolver`` (and ``sparsity``
    for the sparse tier, or ``jacfn``/``jactimes`` for constant-coefficient
    systems with an analytical Jacobian).

    Parameters
    ----------
    parameters : dict[str, float] or None
        Runtime parameter overrides.  Passes ``{}`` when ``None`` to allow
        constant-coefficient detection without explicit overrides.
    solver : str
        ``"ida"`` or ``"cvode"``.  Controls the analytical Jacobian callback
        signature (IDA uses 6-arg ``jacfn`` with ``cj``; CVODE uses 4-arg).
    """
    from tidal.solver._types import (  # noqa: PLC0415
        DENSE_THRESHOLD,
        SPARSE_THRESHOLD,
        SUPERLU_NNZ_LIMIT,
    )
    from tidal.solver.analytical_jacobian import (  # noqa: PLC0415
        try_analytical_jacobian,
    )

    # Analytical Jacobian for constant-coefficient systems
    if try_analytical_jacobian(
        options, spec, layout, grid, bc, parameters or {}, solver=solver,
    ):
        return

    # Existing tier system for non-constant systems
    n_state = layout.total_size
    if n_state <= DENSE_THRESHOLD:
        options["linsolver"] = "dense"
    elif n_state <= SPARSE_THRESHOLD:
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        pattern = build_jacobian_sparsity(spec, layout, grid, bc)
        if pattern.nnz > SUPERLU_NNZ_LIMIT:
            # SuperLU_MT fill-in can reach 20-50x the sparsity pattern nnz for
            # 2D FD systems, exhausting memory. Fall back to FD GMRES instead.
            # FD GMRES (no jactimes callback) is safe with IDACalcIC; see plan
            # comment on analytical GMRES regression (commit 7f3df3e) which was
            # specific to the jactimes callback, not this FD path.
            warnings.warn(
                f"Sparse tier: sparsity pattern nnz={pattern.nnz} exceeds "
                f"SUPERLU_NNZ_LIMIT={SUPERLU_NNZ_LIMIT}. SuperLU_MT fill-in "
                f"may exhaust memory for this {n_state}-state system. "
                f"Falling back to FD GMRES (linsolver='gmres').",
                UserWarning,
                stacklevel=3,
            )
            options["linsolver"] = "gmres"
        else:
            options["linsolver"] = "sparse"
            options["sparsity"] = pattern
    else:
        options["linsolver"] = "gmres"
