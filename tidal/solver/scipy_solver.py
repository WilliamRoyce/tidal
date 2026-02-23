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

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.integrate import (  # pyright: ignore[reportMissingTypeStubs]
    solve_ivp,  # pyright: ignore[reportUnknownVariableType]
)

from tidal.solver.fields import FieldSet
from tidal.solver.leapfrog import compute_force, compute_velocity
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm

# Time-derivative order threshold for dynamical (wave) equations
_SECOND_ORDER = 2

# Implicit methods that benefit from Jacobian sparsity
_IMPLICIT_METHODS = {"Radau", "BDF"}


def _build_rhs_fn(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    kinetic: np.ndarray | None,
    spatial_momenta: Mapping[str, Sequence[OperatorTerm]] | None,
    rhs_eval: RHSEvaluator,
) -> Callable[[float, np.ndarray], np.ndarray]:
    """Build the scipy RHS closure: ``rhs_fn(t, y) -> dydt``."""
    n = grid.num_points
    eq_map = {eq.field_name: i for i, eq in enumerate(spec.equations)}

    def rhs_fn(t: float, y: np.ndarray) -> np.ndarray:
        dydt = np.zeros_like(y)

        force = compute_force(spec, layout, grid, bc, y, t, rhs_eval)
        fieldset = FieldSet.from_flat(layout, grid.shape, y)
        velocity = compute_velocity(
            layout,
            grid,
            kinetic,
            spatial_momenta,
            y,
            fieldset,
            bc,
            t,
            rhs_eval,
        )

        for slot_idx, slot in enumerate(layout.slots):
            s = slice(slot_idx * n, (slot_idx + 1) * n)
            if slot.kind == "momentum":
                dydt[s] = force[s]
            elif slot.kind == "field" and slot.time_order >= _SECOND_ORDER:
                dydt[s] = velocity[s]
            elif slot.kind == "field" and slot.time_order == 1:
                eq_idx = eq_map.get(slot.field_name)
                if eq_idx is not None:
                    result = rhs_eval.evaluate(eq_idx, fieldset, t)
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
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_step: float = np.inf,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
) -> dict[str, Any]:
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
    canonical = spec.canonical

    # Warn about constraint fields (frozen at IC)
    constraint_fields = [
        s.field_name for s in layout.slots if s.kind == "field" and s.time_order == 0
    ]
    if constraint_fields:
        warnings.warn(
            f"scipy: constraint fields {constraint_fields} frozen at initial "
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
    rhs_fn = _build_rhs_fn(spec, layout, grid, bc, kinetic, spatial_momenta, rhs_eval)

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    # Prepare solve_ivp kwargs
    ivp_kwargs: dict[str, Any] = {
        "rtol": rtol,
        "atol": atol,
        "max_step": max_step,
        "t_eval": t_eval,
        "dense_output": False,
    }

    # For implicit methods, provide Jacobian sparsity
    if method in _IMPLICIT_METHODS:
        from tidal.solver.sparsity import build_jacobian_sparsity  # noqa: PLC0415

        pattern = build_jacobian_sparsity(spec, layout, grid, bc)
        ivp_kwargs["jac_sparsity"] = pattern

    result: Any = solve_ivp(rhs_fn, t_span, y0, method=method, **ivp_kwargs)  # pyright: ignore[reportUnknownVariableType]

    # solve_ivp returns y as (n_states, n_times) — transpose to (n_times, n_states)
    y_out: np.ndarray = result.y.T if result.y is not None else np.empty((0, len(y0)))  # pyright: ignore[reportUnknownVariableType]
    t_out: np.ndarray = result.t if result.t is not None else np.empty(0)  # pyright: ignore[reportUnknownVariableType]

    # Call snapshot callback at each output time
    if snapshot_callback is not None and result.success:
        for i in range(len(t_out)):  # pyright: ignore[reportUnknownArgumentType]
            snapshot_callback(t_out[i], y_out[i])  # pyright: ignore[reportUnknownArgumentType]

    msg: str = result.message or (  # pyright: ignore[reportUnknownVariableType]
        "scipy integration completed" if result.success else "scipy integration failed"
    )
    return {
        "t": t_out,
        "y": y_out,
        "success": result.success,
        "message": msg,
    }
