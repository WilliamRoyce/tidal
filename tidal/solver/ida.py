"""SUNDIALS/IDA integration for TIDAL — DAE solver for mixed systems.

Builds IDA-compatible residual functions from TIDAL equation specs.
Handles arbitrary mixes of:
- Second-order (wave) equations via E-L velocity form
- First-order (diffusion/transport) equations
- Algebraic (constraint) equations

Euler-Lagrange velocity form:
- Velocity slot: dv/dt = E-L RHS (second-order equation)
- Field slot:    dq/dt = v      (trivial kinematic)

Reference: Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers", ACM TOMS 31(3), 2005.
scikit-sundae: NREL, BSD-3 license.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL, SECOND_ORDER
from tidal.solver._setup import configure_linear_solver
from tidal.solver._sksundae import SundialsResult, call_ida
from tidal.solver.fields import FieldSet
from tidal.solver.operators import BCSpec, apply_operator, is_periodic_bc
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.solver.rhs import RHSEvaluator
    from tidal.solver.state import SlotInfo
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
    )


class _ResidualCtx:
    """Bundles pre-computed data and per-call arrays for IDA residual evaluation.

    The ``y``, ``yp``, ``res``, and ``fieldset`` attributes are set per-call
    via ``set_arrays()`` and cleared after each residual evaluation.
    """

    def __init__(
        self,
        spec: EquationSystem,
        layout: StateLayout,
        grid: GridInfo,
        bc: BCSpec | None,
        rhs_eval: RHSEvaluator | None = None,
    ) -> None:
        self.spec = spec
        self.layout = layout
        self.grid = grid
        self.bc = bc
        self.n = grid.num_points
        self.shape = grid.shape
        self.rhs_eval = rhs_eval
        self.eq_map: dict[str, int] = spec.equation_map
        # Detect constraints with no self-referencing terms — the field
        # doesn't appear in its own equation (e.g. momentum constraints
        # from gauge DOF).  These are frozen at zero.
        self._no_self_term_fields = self._detect_no_self_term_fields()
        # Detect constraints that need gauge fixing (singular self-operator
        # with periodic BCs, e.g. pure Laplacian → null space = constants).
        # Standard approach: pin mean(field)=0 for one row (FEniCS/Firedrake).
        self._gauge_fix_fields = self._detect_gauge_fix_fields()
        # Per-call state (set via set_arrays)
        self.y: np.ndarray = np.empty(0)
        self.yp: np.ndarray = np.empty(0)
        self.res: np.ndarray = np.empty(0)
        self.fieldset: FieldSet | None = None
        # Dict for legacy constant-coefficient path in compute_rhs (lazy)
        self.fields: dict[str, np.ndarray] | None = None

    def set_arrays(
        self,
        t: float,
        y: np.ndarray,
        yp: np.ndarray,
        res: np.ndarray,
    ) -> None:
        """Bind per-call arrays and unpack fields from y."""
        self.t = t
        self.y = y
        self.yp = yp
        self.res = res
        self.fieldset = FieldSet.from_flat(self.layout, self.shape, y)

        # Inject constraint velocities from yp so that velocity-dependent
        # operators (first_derivative_t, gradient_x of velocity slots)
        # resolve correctly in the RHSEvaluator.
        for eq in self.spec.equations:
            if eq.time_derivative_order == 0:
                slot_idx = self.layout.field_slot_map[eq.field_name]
                start = slot_idx * self.n
                vel = yp[start : start + self.n].reshape(self.shape)
                self.fieldset.set_aux(f"v_{eq.field_name}", vel)

        self.fields = None  # Reset lazy cache

        # Notify coefficient evaluator of new timestep
        if self.rhs_eval is not None:
            self.rhs_eval.begin_timestep(t)

    def compute_rhs(self, eq_idx: int) -> np.ndarray:
        """Sum operator terms for a single equation."""
        if self.rhs_eval is not None and self.fieldset is not None:
            return self.rhs_eval.evaluate(eq_idx, self.fieldset, self.t)

        # Legacy path: constant coefficients only
        if self.fields is None:
            assert self.fieldset is not None
            self.fields = self.fieldset.as_dict()
        eq = self.spec.equations[eq_idx]
        result = np.zeros(self.shape)
        for term in eq.rhs_terms:
            target_data = self.fields.get(term.field, np.zeros(self.shape))
            operated = apply_operator(term.operator, target_data, self.grid, self.bc)
            result += term.coefficient * operated
        return result

    def _detect_no_self_term_fields(self) -> set[str]:
        """Detect constraint equations where the field has no self-referencing terms.

        When a constraint equation's RHS references only *other* fields (not the
        constraint field itself), the Jacobian ∂F/∂field = 0, making IDA's
        Newton solver singular.  These typically arise from gauge degrees of
        freedom (e.g. lapse/shift in linearized gravity) whose EOM constrain
        momenta rather than the field value.

        For such fields, ``handle_constraint`` freezes them at zero
        (``res = y[field]``), making the Jacobian non-singular (identity block).
        The original equation becomes an initial-data consistency condition
        that should be verified by ``check_no_self_term_ic``.

        **Extensibility:** Fields with ``constraint_solver.enabled = True``
        are excluded — the constraint pre-solve path handles them (and
        raises ``ValueError`` if they truly have no self-terms, which is
        the correct fail-fast for gauge-fixed specs that *should* have
        self-terms).  When Phase B gauge-fixing adds self-terms to these
        equations, they will naturally exit this set.

        Emits ``UserWarning`` for each detected field.
        """
        result: set[str] = set()
        for eq in self.spec.equations:
            if eq.time_derivative_order != 0:
                continue
            # Skip fields with constraint_solver enabled — they go through
            # the pre_solve path which has its own validation.
            if eq.constraint_solver.enabled:
                continue
            has_self = any(t.field == eq.field_name for t in eq.rhs_terms)
            if not has_self:
                result.add(eq.field_name)

        # Describe what each frozen equation originally constrains
        for name in sorted(result):
            eq_idx = self.eq_map[name]
            eq = self.spec.equations[eq_idx]
            other_fields = sorted({t.field for t in eq.rhs_terms})
            constraint_desc = ", ".join(other_fields) if other_fields else "none"

            warnings.warn(
                f"Constraint equation for '{name}' has no self-referencing "
                f"terms (field does not appear in its own RHS). Freezing "
                f"'{name}' at zero. The original equation (involving "
                f"[{constraint_desc}]) becomes an initial-data consistency "
                f"condition — ensure the IC satisfies it. Applying a gauge "
                f"condition (e.g. via [gauge] in the TOML) may add "
                f"self-terms that resolve this automatically.",
                UserWarning,
                stacklevel=2,
            )

        return result

    def _detect_gauge_fix_fields(self) -> set[str]:
        """Detect constraints needing gauge fixing (singular operator + periodic BCs).

        A constraint with pure Laplacian self-terms (no identity/mass) and
        periodic BCs has a singular operator (null space = constants).  IDA's
        Newton solver needs a full-rank Jacobian, so we pin one DOF to zero
        — the standard approach used in FEniCS, Firedrake, and PETSc.

        Emits ``UserWarning`` for each gauge-regularized field so the user
        knows this numerical choice is being made.  If the JSON spec already
        carries gauge metadata mentioning the field, an additional conflict
        warning is emitted.
        """
        if not self._all_periodic_bcs():
            return set()

        result = {
            eq.field_name
            for eq in self.spec.equations
            if eq.time_derivative_order == 0
            and eq.constraint_solver.enabled
            and self._is_pure_laplacian(eq)
        }

        gauge_str = self.spec.metadata.get("gauge", "none")

        for name in sorted(result):
            warnings.warn(
                f"Numerical gauge regularization: pinning {name}[0] = 0 to "
                f"resolve singular operator (pure Laplacian + periodic BCs "
                f"\u2192 null space contains constants). The solution is "
                f"unique only up to an additive constant; this pins that "
                f"constant to zero. To disable, set "
                f"constraint_solver.enabled = false for '{name}' in the "
                f"JSON spec.",
                UserWarning,
                stacklevel=2,
            )
            # Check for potential conflict with explicit gauge fixing
            if gauge_str != "none" and name.split("_")[0] in gauge_str:
                warnings.warn(
                    f"JSON spec has gauge metadata '{gauge_str}' which may "
                    f"already constrain '{name}'. The automatic "
                    f"regularization (pin {name}[0]=0) could conflict with "
                    f"the applied gauge. If the gauge was intended to remove "
                    f"the null space, check that it actually modifies the "
                    f"constraint equation structure (e.g. adds an "
                    f"identity/mass term). To suppress, set "
                    f"constraint_solver.enabled = false for '{name}'.",
                    UserWarning,
                    stacklevel=2,
                )

        return result

    def _all_periodic_bcs(self) -> bool:
        """Check if all BCs are periodic."""
        if self.bc is not None:
            bcs = (
                (self.bc,) * self.grid.ndim
                if isinstance(self.bc, str)
                else tuple(self.bc)
            )
            return all(is_periodic_bc(b) for b in bcs)
        return all(self.grid.periodic)

    @staticmethod
    def _is_pure_laplacian(eq: ComponentEquation) -> bool:
        """Check if eq has Laplacian self-terms but no identity/mass self-term."""
        laplacian_ops = {"laplacian", "laplacian_x", "laplacian_y", "laplacian_z"}
        has_lap = False
        for term in eq.rhs_terms:
            if term.field != eq.field_name:
                continue
            if term.operator in laplacian_ops:
                has_lap = True
            elif term.operator == "identity":
                return False  # Has mass term → not singular
        return has_lap

    def handle_constraint(self, slot_idx: int, slot: SlotInfo) -> None:
        """Algebraic constraint: RHS = 0, with special handling for edge cases.

        Three cases:

        1. **No self-terms** (field absent from its own RHS): freeze at zero.
           Residual ``res = y[field]`` → Jacobian = I (non-singular).
        2. **Gauge regularization** (pure Laplacian + periodic BCs): pin
           ``field[0] = 0`` to remove null-space ambiguity.
        3. **Normal**: ``res = RHS(y, t)`` (algebraic equation).
        """
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        eq_idx = self.eq_map.get(slot.field_name)
        if eq_idx is None:
            self.res[s] = 0.0
            return

        # Case 1: no self-terms — freeze field at zero.
        if slot.field_name in self._no_self_term_fields:
            field_slot = self.layout.field_slot_map[slot.field_name]
            fs = slice(field_slot * self.n, (field_slot + 1) * self.n)
            self.res[s] = self.y[fs]
            return

        rhs = self.compute_rhs(eq_idx).ravel()
        self.res[s] = rhs

        # Case 2: gauge fixing — replace first equation with field[0] = 0.
        # For Poisson constraints (pure Laplacian + periodic BCs), the
        # solution is unique up to a constant.  Pinning one point to zero
        # selects the unique solution.  Uses a single diagonal Jacobian
        # entry, compatible with all linear solvers (dense, sparse, GMRES).
        if slot.field_name in self._gauge_fix_fields:
            field_slot = self.layout.field_slot_map[slot.field_name]
            self.res[slot_idx * self.n] = self.y[field_slot * self.n]

    def handle_velocity(self, slot_idx: int, slot: SlotInfo) -> None:
        """E-L equation: dv/dt = RHS."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        eq_idx = self.eq_map[slot.field_name]
        vel_rhs = self.compute_rhs(eq_idx)
        self.res[s] = self.yp[s] - vel_rhs.ravel()

    def handle_dynamical_field(self, slot_idx: int, slot: SlotInfo) -> None:
        """Trivial kinematic: dq/dt = v."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        n = self.n
        vel_slot = self.layout.velocity_slot_map[slot.field_name]
        v = self.y[vel_slot * n : (vel_slot + 1) * n]
        self.res[s] = self.yp[s] - v

    def handle_first_order(self, slot_idx: int, slot: SlotInfo) -> None:
        """First-order: dy/dt = RHS."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        eq_idx = self.eq_map.get(slot.field_name)
        if eq_idx is not None:
            rhs = self.compute_rhs(eq_idx).ravel()
            self.res[s] = self.yp[s] - rhs
        else:
            self.res[s] = self.yp[s]


def build_residual_fn(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    bc: BCSpec | None = None,
    *,
    parameters: dict[str, float] | None = None,
) -> Callable[[float, np.ndarray, np.ndarray, np.ndarray], None]:
    """Build an IDA-compatible residual function from a TIDAL equation spec.

    The returned function has signature ``resfn(t, y, yp, res)`` where
    ``res`` is written in-place.

    For each slot in the state vector:

    - **Constraint** (time_order=0): ``res = RHS(y, t)`` (algebraic, = 0)
    - **First-order field** (time_order=1): ``res = yp - RHS(y, t)``
    - **Second-order field slot**: ``res = yp - v`` (trivial kinematic)
    - **Second-order velocity slot**: ``res = yp - RHS(y, t)`` (E-L equation)

    Parameters
    ----------
    spec : EquationSystem
        Parsed JSON equation specification.
    layout : StateLayout
        State vector layout descriptor.
    grid : GridInfo
        Spatial grid.
    bc : str or tuple of str, optional
        Boundary conditions for spatial operators.
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients. When
        provided, enables position-dependent and time-dependent
        coefficient evaluation via CoefficientEvaluator.

    Returns
    -------
    Callable
        IDA residual function with signature ``(t, y, yp, res) -> None``.
    """
    # Build RHSEvaluator if parameters provided
    rhs_eval: RHSEvaluator | None = None
    if parameters is not None:
        from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
        from tidal.solver.rhs import RHSEvaluator as _RHSEvaluator  # noqa: PLC0415

        coeff_eval = CoefficientEvaluator(spec, grid, parameters)
        rhs_eval = _RHSEvaluator(spec, grid, coeff_eval, bc=bc)

    ctx = _ResidualCtx(
        spec=spec,
        layout=layout,
        grid=grid,
        bc=bc,
        rhs_eval=rhs_eval,
    )

    def residual(
        t: float,
        y: np.ndarray,
        yp: np.ndarray,
        res: np.ndarray,
    ) -> None:
        """IDA residual: F(t, y, y') = 0."""
        ctx.set_arrays(t, y, yp, res)

        for slot_idx, slot in enumerate(layout.slots):
            if slot.time_order == 0:
                ctx.handle_constraint(slot_idx, slot)
            elif slot.kind == "velocity":
                ctx.handle_velocity(slot_idx, slot)
            elif slot.time_order >= SECOND_ORDER and slot.kind == "field":
                ctx.handle_dynamical_field(slot_idx, slot)
            elif slot.time_order == 1:
                ctx.handle_first_order(slot_idx, slot)

    return residual


def solve_ida(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    num_snapshots: int = 101,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    max_steps: int = 50000,
    snapshot_callback: Callable[..., None] | None = None,
    calc_initcond: str | None = None,
    allow_inconsistent_ic: bool = False,
) -> SolverResult:
    """Solve a TIDAL equation system using SUNDIALS/IDA.

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
    num_snapshots : int
        Number of output time points.
    rtol, atol : float
        Relative and absolute tolerances.
    max_steps : int
        Maximum solver steps.
    snapshot_callback : callable, optional
        Called as ``callback(t, y)`` at each output time.
    calc_initcond : str, optional
        IDA initial condition calculation mode. ``"yp0"`` (default for mixed
        DAE) corrects derivatives given y0. ``"y0"`` corrects algebraic
        variables given yp0 — use this for constraint solving where the
        algebraic field values are unknown.
    allow_inconsistent_ic : bool
        If False (default), raise ValueError when constraint equations are
        violated and cannot be solved. If True, issue a warning instead.

    Returns
    -------
    dict
        Result dictionary with keys: ``t``, ``y``, ``success``, ``message``.

    Warns
    -----
    UserWarning
        When constraint pre-solve encounters singular modes (FFT path), or
        when a constraint field is detected as needing gauge regularization
        (pure Laplacian + periodic BCs → ``field[0] = 0`` pinning).  These
        are numerical choices to resolve the null space of the operator.
        To disable for a specific field, set
        ``constraint_solver.enabled = false`` in the JSON spec.
    """
    layout = StateLayout.from_spec(spec, grid.num_points)

    # Unified constraint IC solver: handles both standard constraints
    # (enabled=True, solve for constraint field) and subsidiary constraints
    # (no-self-term equations, solve for free dynamical fields).
    has_constraints = any(
        eq.time_derivative_order == 0 for eq in spec.equations
    )
    if has_constraints:
        from tidal.solver.constraint_solve import ensure_consistent_ic  # noqa: PLC0415

        y0 = ensure_consistent_ic(
            spec, grid, y0,
            bc=bc, parameters=parameters, t=t_span[0],
            strict=not allow_inconsistent_ic,
        )

    resfn = build_residual_fn(spec, layout, grid, bc, parameters=parameters)

    # Initial yp0 — estimate from residual (IDA will correct via calc_initcond)
    yp0 = np.zeros_like(y0)

    # Identify algebraic variables
    alg_idx = layout.algebraic_indices

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    # Configure IDA solver
    options: dict[str, Any] = {
        "rtol": rtol,
        "atol": atol,
        "max_num_steps": max_steps,
    }

    if alg_idx:
        options["algebraic_idx"] = np.array(alg_idx)

    # Compute consistent initial conditions.  IDA needs yp0 that satisfies
    # F(t0, y0, yp0)=0.  The constraint pre-solve above handles nontrivially
    # violated constraints (e.g. Chern-Simons), producing consistent y0.
    # "yp0" mode then fixes y0 and corrects yp0.
    options["calc_initcond"] = calc_initcond or "yp0"
    options["calc_init_dt"] = float(t_eval[1] - t_eval[0])

    configure_linear_solver(options, layout, spec, grid, bc)

    result: SundialsResult = call_ida(resfn, t_eval, y0, yp0, **options)

    # Call snapshot callback at each output time.
    # IDA provides yp (time-derivative vector) which includes constraint
    # velocities — passed to callback for disk storage.
    if snapshot_callback is not None and result.success:
        yp = result.yp
        for i in range(len(result.t)):
            snapshot_callback(
                result.t[i], result.y[i], yp[i] if yp is not None else None
            )

    return {
        "t": result.t,
        "y": result.y,
        "success": result.success,
        "message": result.message,
    }
