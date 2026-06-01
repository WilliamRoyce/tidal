"""SUNDIALS/CVODE integration for TIDAL — adaptive ODE solver for wave systems.

Pure ODE solver using the CVODE module from SUNDIALS via scikit-sundae.
Supports both BDF (stiff, order 1-5) and Adams (non-stiff, order 1-12)
methods with tolerance-controlled adaptive time stepping.

For wave (second-order) equations, the system is reduced to first-order
ODE form using E-L velocity formulation:
    dq/dt = v                   (trivial kinematic)
    dv/dt = E-L RHS             (from spatial operators)

Constraint fields (time_order=0) are frozen at initial values — correct
for gauge-fixed systems (e.g. Coulomb gauge A_0 = 0).

Reference: Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers", ACM TOMS 31(3), 2005.
scikit-sundae: NREL, BSD-3 license.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._setup import (
    build_rhs_evaluator,
    configure_linear_solver,
    warn_frozen_constraints,
)
from tidal.solver._sksundae import SundialsResult, call_cvode, call_cvode_stepwise
from tidal.solver.fields import FieldSet
from tidal.solver.leapfrog import compute_force, compute_velocity
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.progress import SimulationProgress
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem


# ---------------------------------------------------------------------------
# RHS builder
# ---------------------------------------------------------------------------


def _build_rhsfn(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    rhs_eval: RHSEvaluator,
    m_inv: dict[str, float | NDArray[np.float64]] | None,
) -> Callable[[float, np.ndarray, np.ndarray], None]:
    """Build the CVODE RHS closure: ``rhsfn(t, y, yp)``.

    When ``m_inv`` is a dict, apply the per-field inverse kinetic coefficient
    to the velocity slot update so ``d²ₜ q = M⁻¹·K(q)`` for theories with
    non-trivial ``kinetic_coefficient_symbolic`` (#301 / #302). ``None``
    triggers the fast path (M = I).
    """
    eq_map = spec.equation_map
    fs = FieldSet.zeros(layout, grid.shape)
    force_buf = np.zeros(layout.total_size)
    vel_buf = np.zeros(layout.total_size)

    def rhsfn(t: float, y: np.ndarray, yp: np.ndarray) -> None:
        force = compute_force(
            spec,
            layout,
            grid,
            bc,
            y,
            t,
            rhs_eval,
            out=force_buf,
            fieldset=fs,
        )
        velocity = compute_velocity(layout, y, out=vel_buf)

        if m_inv is None:
            for _si, s, _fn in layout.velocity_slot_groups:
                yp[s] = force[s]
        else:
            for _si, s, fn in layout.velocity_slot_groups:
                scale = m_inv.get(fn, 1.0)
                yp[s] = scale * force[s]
        for _si, s, _vs in layout.dynamical_field_slot_groups:
            yp[s] = velocity[s]
        for _si, s, field_name in layout.first_order_slot_groups:
            eq_idx = eq_map.get(field_name)
            if eq_idx is not None:
                result = rhs_eval.evaluate(eq_idx, fs, t)
                yp[s] = result.ravel()
            else:
                yp[s] = 0.0
        for _si, s, _fn in layout.constraint_slot_groups:
            yp[s] = 0.0

    return rhsfn


# ---------------------------------------------------------------------------
# Solver entry point
# ---------------------------------------------------------------------------


def solve_cvode(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    method: str = "BDF",
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    max_step: float = 0.0,
    max_num_steps: int = 50000,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    progress: SimulationProgress | None = None,
) -> SolverResult:
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
    layout = StateLayout.from_spec(spec, grid.num_points)

    warn_frozen_constraints(layout, "CVODE")
    rhs_eval = build_rhs_evaluator(spec, grid, parameters, bc, rtol=rtol)

    # #301 / #302: evaluate kinetic_coefficient_symbolic so d²ₜ q = M⁻¹·K(q)
    # for theories with non-trivial mass matrix. None triggers the M=I fast path.
    from tidal.solver._kinetic import build_inverse_kinetic_diag  # noqa: PLC0415

    m_inv = build_inverse_kinetic_diag(spec, parameters or {}, grid=grid)

    # Build RHS closure
    rhsfn = _build_rhsfn(spec, layout, grid, bc, rhs_eval, m_inv)

    # Configure CVODE solver
    options: dict[str, Any] = {
        "method": method,
        "rtol": rtol,
        "atol": atol,
        "max_num_steps": max_num_steps,
    }
    if max_step > 0:
        options["max_step"] = max_step

    configure_linear_solver(
        options,
        layout,
        spec,
        grid,
        bc,
        parameters=parameters,
        solver="cvode",
    )

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    if progress is not None:
        # Step-by-step mode: progress updates between solver steps (zero overhead)
        result: SundialsResult = call_cvode_stepwise(
            rhsfn,
            t_eval,
            y0,
            progress,
            snapshot_callback=snapshot_callback,
            **options,
        )
    else:
        result = call_cvode(rhsfn, t_eval, y0, **options)

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
