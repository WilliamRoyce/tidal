"""Stormer-Verlet (leapfrog) integrator for TIDAL Hamiltonian systems.

Symplectic integrator for separable Hamiltonians:
    H = T(pi) + V(q)
    T = 1/2 pi^T K^{-1} pi  (kinetic, includes spatial momenta corrections)
    V = V(q)                  (potential, from spatial operators)

The Stormer-Verlet scheme preserves a shadow Hamiltonian to machine
precision, giving zero secular energy drift for conservative systems.

**Not appropriate for**: dissipative systems, absorbing BCs, constraint
damping, or energy outflow -- use IDA instead.

For non-diagonal K: solve K v = (pi - S) per step via ``np.linalg.solve``.
K is typically 1x1 to 10x10 -- microseconds per solve.

Reference: Hairer, Lubich, Wanner, "Geometric Numerical Integration",
Springer, 2006. Chapter VI: Symplectic Integration of Hamiltonian Systems.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver.fields import FieldSet
from tidal.solver.operators import apply_operator
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver.grid import GridInfo
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem

# Time-derivative order threshold for dynamical (wave) equations
_SECOND_ORDER = 2


def _compute_force(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: str | tuple[str, ...] | None,
    y: np.ndarray,
    t: float = 0.0,
    rhs_eval: RHSEvaluator | None = None,
) -> np.ndarray:
    """Compute dpi/dt = F(q) for all momentum slots.

    Evaluates the RHS of each second-order equation (spatial operators
    applied to field values). Returns a flat array over all slots.
    """
    n = grid.num_points
    shape = grid.shape

    fieldset = FieldSet.from_flat(layout, shape, y)
    eq_map = {eq.field_name: i for i, eq in enumerate(spec.equations)}

    force = np.zeros(layout.total_size)
    for slot_idx, slot in enumerate(layout.slots):
        if slot.kind != "momentum":
            continue
        s = slice(slot_idx * n, (slot_idx + 1) * n)
        eq_idx = eq_map.get(slot.field_name)
        if eq_idx is None:
            continue

        if rhs_eval is not None:
            result = rhs_eval.evaluate(eq_idx, fieldset, t)
        else:
            # Legacy path: constant coefficients only
            eq = spec.equations[eq_idx]
            result = np.zeros(shape)
            fields = fieldset.as_dict()
            for term in eq.rhs_terms:
                target_data = fields.get(term.field, np.zeros(shape))
                operated = apply_operator(term.operator, target_data, grid, bc)
                result += term.coefficient * operated

        force[s] = result.ravel()

    return force


def _compute_velocity(  # noqa: PLR0913, PLR0917
    layout: StateLayout,
    grid: GridInfo,
    kinetic: np.ndarray | None,
    spatial_momenta: dict | None,
    pi_flat: np.ndarray,
    fieldset: FieldSet,
    bc: str | tuple[str, ...] | None,
    t: float = 0.0,
    rhs_eval: RHSEvaluator | None = None,
) -> np.ndarray:
    """Compute dq/dt = K^{-1}(pi - S) for all field slots.

    For identity K (most cases), this is just pi.
    For non-diagonal K, solves K v = (pi - S) via np.linalg.solve.
    """
    n = grid.num_points
    velocity = np.zeros(layout.total_size)

    dyn_fields = layout.dynamical_fields
    n_dyn = len(dyn_fields)
    if n_dyn == 0:
        return velocity

    # Build per-field RHS: rhs[i] = pi_i - S_i
    rhs_vecs = np.zeros((n_dyn, n))
    for i, fname in enumerate(dyn_fields):
        mom_slot = layout.momentum_slot_map[fname]
        rhs_vecs[i] = pi_flat[mom_slot * n : (mom_slot + 1) * n]

        if rhs_eval is not None:
            s_i = rhs_eval.evaluate_spatial_momentum(fname, fieldset, t)
            rhs_vecs[i] -= s_i
        elif spatial_momenta and fname in spatial_momenta:
            # Legacy path
            fields = fieldset.as_dict()
            s_i = np.zeros(grid.shape)
            for term in spatial_momenta[fname]:
                target_data = fields.get(term.field, np.zeros(grid.shape))
                operated = apply_operator(term.operator, target_data, grid, bc)
                s_i += term.coefficient * operated
            rhs_vecs[i] -= s_i.ravel()

    # Solve K v = rhs
    if kinetic is None or np.allclose(kinetic, np.eye(n_dyn)):
        vel_vecs = rhs_vecs
    else:
        vel_vecs = np.linalg.solve(kinetic, rhs_vecs)

    for i, fname in enumerate(dyn_fields):
        field_slot = layout.field_slot_map[fname]
        s = slice(field_slot * n, (field_slot + 1) * n)
        velocity[s] = vel_vecs[i]

    return velocity


def _half_kick(
    y: np.ndarray, force: np.ndarray, dt: float,
    layout: StateLayout, n: int,
) -> None:
    """Apply half-kick: pi += (dt/2) F(q), in-place."""
    for slot_idx, slot in enumerate(layout.slots):
        if slot.kind == "momentum":
            s = slice(slot_idx * n, (slot_idx + 1) * n)
            y[s] += 0.5 * dt * force[s]


def _drift(
    y: np.ndarray, velocity: np.ndarray, dt: float,
    layout: StateLayout, n: int,
) -> None:
    """Apply drift: q += dt v, in-place."""
    for slot_idx, slot in enumerate(layout.slots):
        if slot.kind == "field" and slot.time_order >= _SECOND_ORDER:
            s = slice(slot_idx * n, (slot_idx + 1) * n)
            y[s] += dt * velocity[s]


def solve_leapfrog(  # noqa: C901, PLR0913, PLR0914, PLR0915
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    dt: float,
    *,
    bc: str | tuple[str, ...] | None = None,
    parameters: dict[str, float] | None = None,
    snapshot_interval: float | None = None,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
) -> dict[str, Any]:
    """Solve a TIDAL Hamiltonian system using Stormer-Verlet (leapfrog).

    Works for second-order (wave) equations with optional constraint fields.
    Constraint fields (time_order=0) are frozen at their initial values —
    correct for gauge-fixed systems (e.g. Coulomb gauge A_0 = 0).

    Parameters
    ----------
    spec : EquationSystem
        Parsed equation specification.
    grid : GridInfo
        Spatial grid.
    y0 : np.ndarray
        Initial state vector (flat: [q_0, pi_0, q_1, pi_1, ...]).
    t_span : tuple[float, float]
        (t_start, t_end).
    dt : float
        Fixed time step.
    bc : str or tuple, optional
        Boundary conditions.
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients.
    snapshot_interval : float, optional
        Time between snapshots. If None, only initial and final states saved.
    snapshot_callback : callable, optional
        Called as ``callback(t, y)`` at each snapshot time.

    Returns
    -------
    dict
        Result dictionary with keys: ``t``, ``y``, ``success``, ``message``.

    Raises
    ------
    ValueError
        If the system contains first-order (diffusion/transport) equations.

    Warns
    -----
    UserWarning
        If constraint fields (time_order=0) are present — they remain frozen
        at initial values.
    """
    layout = StateLayout.from_spec(spec, grid.num_points)
    n = grid.num_points
    canonical = spec.canonical

    # Validate: leapfrog supports second-order (wave) + frozen constraints.
    # First-order (diffusion/transport) equations require IDA.
    constraint_fields: list[str] = []
    for slot in layout.slots:
        if slot.kind == "momentum":
            continue
        if slot.time_order == 0:
            constraint_fields.append(slot.field_name)
        elif slot.time_order == 1:
            msg = (
                f"Leapfrog does not support first-order equations. "
                f"Field '{slot.field_name}' has time_order={slot.time_order}. "
                f"Use --scheme ida for first-order (diffusion/transport) systems."
            )
            raise ValueError(msg)

    if constraint_fields:
        warnings.warn(
            f"Leapfrog: constraint fields {constraint_fields} frozen at initial "
            f"values (not evolved). Correct for gauge-fixed systems (e.g. A_0=0).",
            stacklevel=2,
        )

    # Pre-extract kinetic matrix
    kinetic = None
    if canonical and canonical.kinetic_matrix:
        kinetic = np.array(canonical.kinetic_matrix.to_dense())
    spatial_momenta = canonical.spatial_momenta if canonical else None

    # Build RHSEvaluator if parameters provided
    rhs_eval: RHSEvaluator | None = None
    if parameters is not None:
        from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
        from tidal.solver.rhs import RHSEvaluator as _RHSEvaluator  # noqa: PLC0415

        coeff_eval = CoefficientEvaluator(spec, grid, parameters)
        rhs_eval = _RHSEvaluator(spec, grid, coeff_eval, bc=bc)

    # Initialize state
    y = y0.copy()
    t = t_span[0]
    t_end = t_span[1]

    # Snapshot collection
    if snapshot_interval is None:
        snapshot_interval = t_end - t
    next_snapshot = t
    times: list[float] = []
    snapshots: list[np.ndarray] = []

    def _save(t_now: float) -> None:
        times.append(t_now)
        snapshots.append(y.copy())
        if snapshot_callback is not None:
            snapshot_callback(t_now, y)

    _save(t)

    # Leapfrog loop
    n_steps = int((t_end - t) / dt)
    for _step in range(n_steps):
        # Half-kick
        force = _compute_force(spec, layout, grid, bc, y, t, rhs_eval)
        _half_kick(y, force, dt, layout, n)

        # Drift
        fieldset = FieldSet.from_flat(layout, grid.shape, y)
        velocity = _compute_velocity(
            layout, grid, kinetic, spatial_momenta, y, fieldset, bc, t, rhs_eval,
        )
        _drift(y, velocity, dt, layout, n)

        # Half-kick
        force = _compute_force(spec, layout, grid, bc, y, t, rhs_eval)
        _half_kick(y, force, dt, layout, n)

        t += dt

        # Snapshot check
        if t >= next_snapshot + snapshot_interval - dt * 0.01:
            _save(t)
            next_snapshot += snapshot_interval

    # Ensure final state is saved
    if not times or abs(times[-1] - t) > dt * 0.01:
        _save(t)

    return {
        "t": np.array(times),
        "y": np.array(snapshots),
        "success": True,
        "message": "Leapfrog integration completed",
    }
