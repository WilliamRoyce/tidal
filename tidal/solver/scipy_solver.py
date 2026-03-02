"""scipy.integrate.solve_ivp wrapper for TIDAL — explicit adaptive ODE solver.

Provides DOP853 (8th-order Dormand-Prince) as default method, which excels
for smooth non-stiff wave problems with fewer function evaluations than
implicit BDF for the same accuracy.

Also supports implicit methods (Radau, BDF) with Jacobian sparsity from
the existing TIDAL sparsity builder.

Reference: Dormand & Prince, "A family of embedded Runge-Kutta formulae",
J. Comp. Appl. Math. 6, 1980.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._scipy_types import IVPResult, call_solve_ivp
from tidal.solver._setup import build_rhs_evaluator, warn_frozen_constraints
from tidal.solver.fields import FieldSet
from tidal.solver.leapfrog import compute_force, compute_velocity
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.progress import SimulationProgress
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem

# Implicit methods that benefit from Jacobian sparsity
_IMPLICIT_METHODS = {"Radau", "BDF"}


def _build_rhs_fn(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    rhs_eval: RHSEvaluator,
    progress: SimulationProgress | None = None,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Build the scipy RHS closure: ``rhs_fn(t, y) -> dydt``."""
    eq_map = spec.equation_map
    fs = FieldSet.zeros(layout, grid.shape)

    def rhs_fn(t: float, y: np.ndarray) -> np.ndarray:
        if progress is not None:
            progress.update(t)

        dydt = np.zeros_like(y)

        force = compute_force(spec, layout, grid, bc, y, t, rhs_eval, fieldset=fs)
        velocity = compute_velocity(layout, y)

        for _si, s, _fn in layout.velocity_slot_groups:
            dydt[s] = force[s]
        for _si, s, _vs in layout.dynamical_field_slot_groups:
            dydt[s] = velocity[s]
        for _si, s, field_name in layout.first_order_slot_groups:
            eq_idx = eq_map.get(field_name)
            if eq_idx is not None:
                result = rhs_eval.evaluate(eq_idx, fs, t)
                dydt[s] = result.ravel()

        return dydt

    return rhs_fn


def solve_scipy(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    method: str = "DOP853",
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    max_step: float = np.inf,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    progress: SimulationProgress | None = None,
) -> SolverResult:
    """Solve a TIDAL equation system using scipy.integrate.solve_ivp.

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
        Integration method: ``'DOP853'`` (default), ``'RK45'``, ``'Radau'``,
        ``'BDF'``, ``'RK23'``, ``'LSODA'``.
    rtol, atol : float
        Relative and absolute tolerances.
    max_step : float
        Maximum step size. ``np.inf`` (default) = unbounded.
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

    warn_frozen_constraints(layout, "scipy")
    rhs_eval = build_rhs_evaluator(spec, grid, parameters, bc)

    # Build RHS closure (with optional progress tracking)
    rhs_fn = _build_rhs_fn(spec, layout, grid, bc, rhs_eval, progress=progress)

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    # For implicit methods, provide Jacobian sparsity
    jac_sparsity = None
    if method in _IMPLICIT_METHODS:
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        jac_sparsity = build_jacobian_sparsity(spec, layout, grid, bc)

    result: IVPResult = call_solve_ivp(
        rhs_fn,
        t_span,
        y0,
        method=method,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        jac_sparsity=jac_sparsity,
    )

    if progress is not None:
        progress.finish()

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
