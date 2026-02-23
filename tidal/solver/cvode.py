"""SUNDIALS/CVODE integration for TIDAL — adaptive ODE solver for wave systems.

Pure ODE solver using the CVODE module from SUNDIALS via scikit-sundae.
Supports both BDF (stiff, order 1-5) and Adams (non-stiff, order 1-12)
methods with tolerance-controlled adaptive time stepping.

For wave (second-order) equations, the system is reduced to first-order
ODE form:
    dq/dt = K^{-1}(pi - S)    (velocity from momenta)
    dpi/dt = F(q)               (force from spatial operators)

Constraint fields (time_order=0) are frozen at initial values — correct
for gauge-fixed systems (e.g. Coulomb gauge A_0 = 0).

Reference: Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers", ACM TOMS 31(3), 2005.
scikit-sundae: NREL, BSD-3 license.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver.fields import FieldSet
from tidal.solver.leapfrog import compute_force, compute_velocity
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm

# System size thresholds for linear solver selection (same as IDA).
_DENSE_THRESHOLD = 2_000
_SPARSE_THRESHOLD = 200_000

# Time-derivative order threshold for dynamical (wave) equations
_SECOND_ORDER = 2


def _build_rhsfn(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    kinetic: np.ndarray | None,
    spatial_momenta: Mapping[str, Sequence[OperatorTerm]] | None,
    rhs_eval: RHSEvaluator,
) -> Callable[[float, np.ndarray, np.ndarray], None]:
    """Build the CVODE RHS closure: ``rhsfn(t, y, yp)``."""
    n = grid.num_points
    eq_map = {eq.field_name: i for i, eq in enumerate(spec.equations)}

    def rhsfn(t: float, y: np.ndarray, yp: np.ndarray) -> None:
        force = compute_force(spec, layout, grid, bc, y, t, rhs_eval)
        fieldset = FieldSet.from_flat(layout, grid.shape, y)
        velocity = compute_velocity(
            layout, grid, kinetic, spatial_momenta, y, fieldset, bc, t, rhs_eval,
        )

        for slot_idx, slot in enumerate(layout.slots):
            s = slice(slot_idx * n, (slot_idx + 1) * n)
            if slot.kind == "momentum":
                yp[s] = force[s]
            elif slot.kind == "field" and slot.time_order >= _SECOND_ORDER:
                yp[s] = velocity[s]
            elif slot.kind == "field" and slot.time_order == 1:
                eq_idx = eq_map.get(slot.field_name)
                if eq_idx is not None:
                    result = rhs_eval.evaluate(eq_idx, fieldset, t)
                    yp[s] = result.ravel()
                else:
                    yp[s] = 0.0
            else:
                yp[s] = 0.0

    return rhsfn


def solve_cvode(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    method: str = "BDF",
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_step: float = 0.0,
    max_num_steps: int = 50000,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
) -> dict[str, Any]:
    """Solve a TIDAL equation system using SUNDIALS/CVODE.

    Parameters
    ----------
    spec : EquationSystem
        Parsed equation specification.
    grid : GridInfo
        Spatial grid.
    y0 : np.ndarray
        Initial state vector (flat).
    t_span : tuple[float, float]
        (t_start, t_end).
    bc : str or tuple, optional
        Boundary conditions.
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients.
    method : str
        Integration method: ``'BDF'`` (stiff, default) or ``'Adams'`` (non-stiff).
    rtol, atol : float
        Relative and absolute tolerances.
    max_step : float
        Maximum step size. ``0`` (default) = unbounded (SUNDIALS default).
    max_num_steps : int
        Maximum solver steps per output interval.
    num_snapshots : int
        Number of output time points.
    snapshot_callback : callable, optional
        Called as ``callback(t, y)`` at each output time.

    Returns
    -------
    dict
        Result dictionary with keys: ``t``, ``y``, ``success``, ``message``.

    Warns
    -----
    UserWarning
        If constraint fields (time_order=0) are present — they remain frozen
        at initial values.
    """
    from sksundae.cvode import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        CVODE,
    )

    layout = StateLayout.from_spec(spec, grid.num_points)
    canonical = spec.canonical

    # Warn about constraint fields (frozen at IC)
    constraint_fields = [
        s.field_name for s in layout.slots
        if s.kind == "field" and s.time_order == 0
    ]
    if constraint_fields:
        warnings.warn(
            f"CVODE: constraint fields {constraint_fields} frozen at initial "
            f"values (not evolved). Use --scheme ida for constraint systems.",
            stacklevel=2,
        )

    # Pre-extract kinetic matrix
    kinetic = None
    if canonical and canonical.kinetic_matrix:
        kinetic = np.array(canonical.kinetic_matrix.to_dense())
    spatial_momenta = canonical.spatial_momenta if canonical else None

    # Build RHSEvaluator for coefficient resolution
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
    from tidal.solver.rhs import RHSEvaluator  # noqa: PLC0415

    coeff_eval = CoefficientEvaluator(spec, grid, parameters or {})
    rhs_eval = RHSEvaluator(spec, grid, coeff_eval, bc=bc)

    # Build RHS closure
    rhsfn = _build_rhsfn(spec, layout, grid, bc, kinetic, spatial_momenta, rhs_eval)

    # Configure CVODE solver
    options: dict[str, Any] = {
        "method": method,
        "rtol": rtol,
        "atol": atol,
        "max_num_steps": max_num_steps,
    }
    if max_step > 0:
        options["max_step"] = max_step

    # Choose linear solver based on system size (same thresholds as IDA)
    n_state = layout.total_size
    if n_state <= _DENSE_THRESHOLD:
        options["linsolver"] = "dense"
    elif n_state <= _SPARSE_THRESHOLD:
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        pattern = build_jacobian_sparsity(spec, layout, grid, bc)
        options["linsolver"] = "sparse"
        options["sparsity"] = pattern
    else:
        options["linsolver"] = "gmres"

    solver = CVODE(rhsfn, **options)

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    result: Any = solver.solve(t_eval, y0)

    # Call snapshot callback at each output time
    if snapshot_callback is not None and result.success:
        for i in range(len(result.t)):
            snapshot_callback(result.t[i], result.y[i])

    return {
        "t": result.t,
        "y": result.y,
        "success": result.success,
        "message": result.message,
    }
