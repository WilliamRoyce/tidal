"""Stormer-Verlet (leapfrog) integrator for TIDAL systems.

Symplectic integrator for second-order E-L equations using velocity form:
    dq/dt = v          (trivial kinematic)
    dv/dt = E-L RHS    (from equations[] array)

The Stormer-Verlet scheme preserves a shadow Hamiltonian to machine
precision, giving zero secular energy drift for conservative systems
with velocity-independent forces.

**Not appropriate for**: dissipative systems, absorbing BCs, constraint
damping, or energy outflow -- use IDA instead.

Reference: Hairer, Lubich, Wanner, "Geometric Numerical Integration",
Springer, 2006. Chapter VI: Symplectic Integration of Hamiltonian Systems.
"""

from __future__ import annotations

import math
import warnings
from typing import TYPE_CHECKING

import numpy as np

from tidal.solver._setup import build_rhs_evaluator
from tidal.solver.fields import FieldSet
from tidal.solver.operators import BCSpec, apply_operator
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.solver.progress import SimulationProgress
    from tidal.solver.rhs import RHSEvaluator
    from tidal.symbolic.json_loader import EquationSystem


def compute_force(  # noqa: PLR0913, PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None,
    y: np.ndarray,
    t: float = 0.0,
    rhs_eval: RHSEvaluator | None = None,
    out: np.ndarray | None = None,
    fieldset: FieldSet | None = None,
) -> np.ndarray:
    """Compute dv/dt = E-L RHS for all velocity slots.

    Parameters
    ----------
    out : np.ndarray, optional
        Pre-allocated output buffer of length ``layout.total_size``.
        If provided, filled in-place (avoids allocation). If ``None``,
        a fresh array is allocated.
    fieldset : FieldSet, optional
        Reusable FieldSet.  If provided, rebound to *y* in-place
        (avoids per-call object allocation).
    """
    shape = grid.shape

    if fieldset is not None:
        fieldset.rebind(y)
    else:
        fieldset = FieldSet.from_flat(layout, shape, y)
    eq_map = spec.equation_map

    if out is not None:
        force = out
        force.fill(0.0)
    else:
        force = np.zeros(layout.total_size)
    for _slot_idx, s, field_name in layout.velocity_slot_groups:
        eq_idx = eq_map.get(field_name)
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


def compute_velocity(
    layout: StateLayout,
    y: np.ndarray,
    out: np.ndarray | None = None,
) -> np.ndarray:
    """Read dq/dt = v directly from velocity slots in the state vector.

    Parameters
    ----------
    out : np.ndarray, optional
        Pre-allocated output buffer. If provided, filled in-place.
    """
    if out is not None:
        velocity = out
        velocity.fill(0.0)
    else:
        velocity = np.zeros(layout.total_size)

    for _slot_idx, s, vel_slot in layout.dynamical_field_slot_groups:
        velocity[s] = y[layout.slot_slice(vel_slot)]

    return velocity


def _half_kick(
    y: np.ndarray,
    force: np.ndarray,
    dt: float,
    layout: StateLayout,
) -> None:
    """Apply half-kick: v += (dt/2) F(q), in-place."""
    for _slot_idx, s, _field_name in layout.velocity_slot_groups:
        y[s] += 0.5 * dt * force[s]


def solve_leapfrog(  # noqa: C901, PLR0912, PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    dt: float,
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    snapshot_interval: float | None = None,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    progress: SimulationProgress | None = None,
) -> SolverResult:
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
        Initial state vector (flat: [q_0, v_0, q_1, v_1, ...]).
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

    # Validate: leapfrog supports second-order (wave) + frozen constraints.
    # First-order (diffusion/transport) equations require IDA.
    constraint_fields: list[str] = []
    for slot in layout.slots:
        if slot.kind == "velocity":
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

    # Build RHSEvaluator if parameters provided
    rhs_eval: RHSEvaluator | None = None
    if parameters is not None:
        rhs_eval = build_rhs_evaluator(spec, grid, parameters, bc)

    # Initialize state
    y = y0.copy()
    t = t_span[0]
    t_end = t_span[1]

    # Snapshot collection — use integer-based indexing to avoid FP drift.
    # Each snapshot trigger recomputes the target time from the integer
    # index, preventing accumulated floating-point error in long runs.
    if snapshot_interval is None:
        snapshot_interval = t_end - t
    snapshot_idx = 0
    times: list[float] = []
    snapshots: list[np.ndarray] = []

    def _save(t_now: float) -> None:
        times.append(t_now)
        snapshots.append(y.copy())
        if snapshot_callback is not None:
            snapshot_callback(t_now, y)

    _save(t)

    # Leapfrog loop — use ceil to guarantee reaching t_end.
    # int() truncates (e.g. int(153.8) = 153), causing t_final < t_end and
    # missing the last snapshot.  ceil() ensures we always reach or slightly
    # overshoot t_end.  The -1e-10 handles exact division (10.0/0.1 = 100.0)
    # without adding an unnecessary extra step.
    n_steps = max(1, math.ceil((t_end - t) / dt - 1e-10))
    force_buf = np.zeros(layout.total_size)
    fieldset_buf = FieldSet.zeros(layout, grid.shape)
    drift_pairs = layout.drift_slot_pairs

    # KDK force caching: compute F(q_0) once before the loop.  Between
    # steps only velocity changes (half-kicks), but force depends only on
    # position, so F(q_{n+1}) from step n is reused at the start of step
    # n+1.  This halves force evaluations (2N → N+1).
    compute_force(
        spec,
        layout,
        grid,
        bc,
        y,
        t,
        rhs_eval,
        out=force_buf,
        fieldset=fieldset_buf,
    )

    for _step in range(n_steps):
        # Half-kick with cached F(q_n)
        _half_kick(y, force_buf, dt, layout)

        # Drift: q += dt * v (zero-copy, reads velocity directly from y)
        for field_slice, vel_slice in drift_pairs:
            y[field_slice] += dt * y[vel_slice]

        # Force at new position F(q_{n+1}) — cached for next step
        compute_force(
            spec,
            layout,
            grid,
            bc,
            y,
            t + dt,
            rhs_eval,
            out=force_buf,
            fieldset=fieldset_buf,
        )

        # Half-kick with F(q_{n+1})
        _half_kick(y, force_buf, dt, layout)

        t += dt

        if progress is not None:
            progress.update(t)

        # Snapshot check (integer-based to avoid FP accumulation)
        if t >= (snapshot_idx + 1) * snapshot_interval - dt * 0.01:
            _save(t)
            snapshot_idx += 1

    # Ensure final state is saved
    if not times or abs(times[-1] - t) > dt * 0.01:
        _save(t)

    if progress is not None:
        progress.finish()

    return {
        "t": np.asarray(times, dtype=np.float64),
        "y": np.asarray(snapshots, dtype=np.float64),
        "success": True,
        "message": "Leapfrog integration completed",
    }
