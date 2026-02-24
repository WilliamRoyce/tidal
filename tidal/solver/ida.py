"""SUNDIALS/IDA integration for TIDAL — DAE solver for mixed systems.

Builds IDA-compatible residual functions from TIDAL equation specs.
Handles arbitrary mixes of:
- Second-order (wave) equations via Hamiltonian splitting
- First-order (diffusion/transport) equations
- Algebraic (constraint) equations

The kinetic matrix K is used directly in residual form:
    K_{ij} * dq_j/dt - (pi_i - S_i) = 0
No K^{-1} inversion needed — IDA's Newton iteration handles it.

Reference: Hindmarsh et al., "SUNDIALS: Suite of Nonlinear and
Differential/Algebraic Equation Solvers", ACM TOMS 31(3), 2005.
scikit-sundae: NREL, BSD-3 license.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._sksundae import SundialsResult, call_ida
from tidal.solver._types import DENSE_THRESHOLD, SPARSE_THRESHOLD, SolverResult
from tidal.solver.fields import FieldSet
from tidal.solver.operators import BCSpec, apply_operator, is_periodic_bc
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver.grid import GridInfo
    from tidal.solver.rhs import RHSEvaluator
    from tidal.solver.state import SlotInfo
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        OperatorTerm,
    )

# Time-derivative order threshold for dynamical (wave) equations
_SECOND_ORDER = 2


class _ResidualCtx:
    """Bundles pre-computed data and per-call arrays for IDA residual evaluation.

    The ``y``, ``yp``, ``res``, and ``fieldset`` attributes are set per-call
    via ``set_arrays()`` and cleared after each residual evaluation.
    """

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        spec: EquationSystem,
        layout: StateLayout,
        grid: GridInfo,
        bc: BCSpec | None,
        kinetic: np.ndarray | None,
        spatial_momenta: dict[str, tuple[OperatorTerm, ...]] | None,
        rhs_eval: RHSEvaluator | None = None,
    ) -> None:
        self.spec = spec
        self.layout = layout
        self.grid = grid
        self.bc = bc
        self.n = grid.num_points
        self.shape = grid.shape
        self.kinetic = kinetic
        self.spatial_momenta = spatial_momenta
        self.rhs_eval = rhs_eval
        self.eq_map: dict[str, int] = {
            eq.field_name: i for i, eq in enumerate(spec.equations)
        }
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
        # Legacy dict for kinetic/field_rates handlers that still use raw arrays
        self.fields: dict[str, np.ndarray] = {}
        # Track which fallback warnings have been issued (fire once per field)
        self._warned_field_rates: set[str] = set()
        self._warned_identity_k: set[str] = set()

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

        # Inject constraint velocities from yp so that
        # first_derivative_t(constraint) and gradient_x(pi_constraint)
        # resolve correctly in the RHSEvaluator.
        for eq in self.spec.equations:
            if eq.time_derivative_order == 0:
                slot_idx = self.layout.field_slot_map[eq.field_name]
                start = slot_idx * self.n
                vel = yp[start : start + self.n].reshape(self.shape)
                self.fieldset.set_aux(f"pi_{eq.field_name}", vel)

        self.fields = self.fieldset.as_dict()

        # Notify coefficient evaluator of new timestep
        if self.rhs_eval is not None:
            self.rhs_eval.begin_timestep(t)

    def compute_rhs(self, eq_idx: int) -> np.ndarray:
        """Sum operator terms for a single equation."""
        if self.rhs_eval is not None and self.fieldset is not None:
            return self.rhs_eval.evaluate(eq_idx, self.fieldset, self.t)

        # Legacy path: constant coefficients only
        eq = self.spec.equations[eq_idx]
        result = np.zeros(self.shape)
        for term in eq.rhs_terms:
            target_data = self.fields.get(term.field, np.zeros(self.shape))
            operated = apply_operator(term.operator, target_data, self.grid, self.bc)
            result += term.coefficient * operated
        return result

    def compute_spatial_mom(self, field_name: str) -> np.ndarray:
        """Compute S_i for dynamical field i from spatial_momenta terms."""
        if self.rhs_eval is not None and self.fieldset is not None:
            return self.rhs_eval.evaluate_spatial_momentum(
                field_name, self.fieldset, self.t
            )

        # Legacy path
        if self.spatial_momenta is None or field_name not in self.spatial_momenta:
            return np.zeros(self.n)
        result = np.zeros(self.shape)
        for term in self.spatial_momenta[field_name]:
            target_data = self.fields.get(term.field, np.zeros(self.shape))
            operated = apply_operator(term.operator, target_data, self.grid, self.bc)
            result += term.coefficient * operated
        return result.ravel()

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

    def handle_momentum(self, slot_idx: int, slot: SlotInfo) -> None:
        """Hamilton's 2nd: dpi/dt = RHS."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        eq_idx = self.eq_map[slot.field_name]
        mom_rhs = self.compute_rhs(eq_idx)
        self.res[s] = self.yp[s] - mom_rhs.ravel()

    def handle_dynamical_field(self, slot_idx: int, slot: SlotInfo) -> None:
        """Hamilton's 1st: K*dq/dt = pi - S (or fallback)."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        dyn_i = slot.dynamical_index
        canonical = self.spec.canonical

        if dyn_i is not None and self.kinetic is not None:
            self._handle_kinetic(s, dyn_i, slot.field_name)
        elif canonical and slot.field_name in canonical.field_rates:
            self._handle_field_rates(s, slot.field_name)
        else:
            self._handle_identity_k(s, slot.field_name)

    def _handle_kinetic(self, s: slice, dyn_i: int, field_name: str) -> None:
        """K_{ij} * yp_j - (pi_i - S_i) = 0."""
        assert self.kinetic is not None
        n = self.n
        k_yp = np.zeros(n)
        for j, dyn_field in enumerate(self.layout.dynamical_fields):
            if self.kinetic[dyn_i, j] != 0:
                fs = self.layout.field_slot_map[dyn_field]
                k_yp += self.kinetic[dyn_i, j] * self.yp[fs * n : (fs + 1) * n]

        mom_slot = self.layout.momentum_slot_map[field_name]
        pi = self.y[mom_slot * n : (mom_slot + 1) * n]
        s_i = self.compute_spatial_mom(field_name)
        self.res[s] = k_yp - (pi - s_i)

    def _handle_field_rates(self, s: slice, field_name: str) -> None:
        """Fallback: use field_rates when K not available.

        This path fires when the JSON spec has ``field_rates`` (old format
        with K^{-1} embedded) but no ``kinetic_matrix``.  Results are
        correct, but the spec should be regenerated with the current
        pipeline to get the preferred ``kinetic_matrix`` format.
        """
        if field_name not in self._warned_field_rates:
            self._warned_field_rates.add(field_name)
            warnings.warn(
                f"Using field_rates fallback for '{field_name}' because "
                f"kinetic_matrix is missing from the JSON spec. This uses "
                f"the old K^{{-1}} format. Regenerate JSON with current "
                f"pipeline ('tidal derive') to use the preferred "
                f"kinetic_matrix format.",
                UserWarning,
                stacklevel=2,
            )
        canonical = self.spec.canonical
        assert canonical is not None
        rate = np.zeros(self.shape)
        for term in canonical.field_rates[field_name]:
            target_data = self.fields.get(term.field, np.zeros(self.shape))
            operated = apply_operator(term.operator, target_data, self.grid, self.bc)
            rate += term.coefficient * operated
        self.res[s] = self.yp[s] - rate.ravel()

    def _handle_identity_k(self, s: slice, field_name: str) -> None:
        """Identity K: dq/dt = pi.

        This path fires when neither ``kinetic_matrix`` nor ``field_rates``
        is available.  For single-field systems this is always correct
        (K = 1).  For multi-field systems, K might be non-diagonal, and
        assuming K = I would silently give wrong dynamics.
        """
        if field_name not in self._warned_identity_k:
            self._warned_identity_k.add(field_name)
            n_dyn = len(self.layout.dynamical_fields)
            if n_dyn > 1:
                warnings.warn(
                    f"Multi-field system ({n_dyn} dynamical fields) "
                    f"missing both kinetic_matrix and field_rates for "
                    f"'{field_name}'; assuming K = I (identity). If the "
                    f"kinetic matrix is non-diagonal, results will be "
                    f"incorrect. Regenerate JSON with current pipeline "
                    f"('tidal derive').",
                    UserWarning,
                    stacklevel=2,
                )
        n = self.n
        mom_slot = self.layout.momentum_slot_map[field_name]
        pi = self.y[mom_slot * n : (mom_slot + 1) * n]
        self.res[s] = self.yp[s] - pi

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
    - **Second-order field slot**: ``res = K_{ij} * yp_j - (pi_i - S_i)``
    - **Second-order momentum slot**: ``res = yp - RHS(y, t)``

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
    canonical = spec.canonical

    # Pre-extract kinetic matrix as dense numpy array
    if canonical and canonical.kinetic_matrix:
        kinetic = np.array(canonical.kinetic_matrix.to_dense())
    else:
        kinetic = None

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
        kinetic=kinetic,
        spatial_momenta=canonical.spatial_momenta if canonical else None,
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
            elif slot.kind == "momentum":
                ctx.handle_momentum(slot_idx, slot)
            elif slot.time_order >= _SECOND_ORDER and slot.kind == "field":
                ctx.handle_dynamical_field(slot_idx, slot)
            elif slot.time_order == 1:
                ctx.handle_first_order(slot_idx, slot)

    return residual


_IC_RESIDUAL_TOL = 1e-10
"""Tolerance for initial-data subsidiary constraint residual check."""


def _check_no_self_term_ic(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    y0: np.ndarray,
    bc: BCSpec | None,
) -> None:
    """Verify initial data satisfies no-self-term constraint equations.

    For each constraint equation where the field has no self-referencing
    terms (frozen at zero by IDA), evaluates the original constraint RHS
    at t=0 and warns if the residual is non-negligible.

    These equations are subsidiary conditions (e.g. momentum constraints,
    transversality) that must be **jointly** satisfied by the initial data
    for physical consistency.  When multiple constraints are violated, a
    single summary warning lists all of them so the user can see the full
    set of conditions that need to be met.
    """
    n = grid.num_points

    # Build field dict from y0
    fields: dict[str, np.ndarray] = {}
    for name, slot_idx in layout.field_slot_map.items():
        fields[name] = y0[slot_idx * n : (slot_idx + 1) * n].reshape(grid.shape)
    for name, slot_idx in layout.momentum_slot_map.items():
        fields[f"pi_{name}"] = y0[slot_idx * n : (slot_idx + 1) * n].reshape(grid.shape)

    violations: list[tuple[str, float, list[str]]] = []
    for eq in spec.equations:
        if eq.time_derivative_order != 0:
            continue
        if eq.constraint_solver.enabled:
            continue
        has_self = any(t.field == eq.field_name for t in eq.rhs_terms)
        if has_self:
            continue

        rhs = _eval_constraint_rhs(eq, fields, grid, bc)
        max_res = float(np.max(np.abs(rhs)))
        if max_res > _IC_RESIDUAL_TOL:
            other_fields = sorted({t.field for t in eq.rhs_terms})
            violations.append((eq.field_name, max_res, other_fields))

    if not violations:
        return

    # Build a comprehensive summary of all violated subsidiary constraints
    lines = [
        "Initial data does not satisfy subsidiary constraint(s) "
        "(not enforced at runtime — frozen fields):"
    ]
    for field_name, max_res, involved in violations:
        lines.append(
            f"  {field_name}: max|residual| = {max_res:.2e} "
            f"(involves [{', '.join(involved)}])"
        )
    lines.append(
        "For physical consistency, choose initial conditions that "
        "jointly satisfy all subsidiary constraints."
    )
    warnings.warn("\n".join(lines), UserWarning, stacklevel=2)


def _eval_constraint_rhs(
    eq: ComponentEquation,
    fields: dict[str, np.ndarray],
    grid: GridInfo,
    bc: BCSpec | None,
) -> np.ndarray:
    """Evaluate a constraint equation's RHS using constant coefficients."""
    rhs = np.zeros(grid.shape)
    for term in eq.rhs_terms:
        target = fields.get(term.field, np.zeros(grid.shape))
        operated = apply_operator(term.operator, target, grid, bc)
        rhs += term.coefficient * operated
    return rhs


def solve_ida(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    num_snapshots: int = 101,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_steps: int = 50000,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    calc_initcond: str | None = None,
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

    # Pre-solve constraints to produce consistent y0.
    # Without this, IDA fails on systems like Chern-Simons where
    # the constraint has a nontrivial source at t=0.
    has_enabled_constraints = any(
        eq.time_derivative_order == 0 and eq.constraint_solver.enabled
        for eq in spec.equations
    )
    if has_enabled_constraints:
        from tidal.solver.constraint_solve import pre_solve_constraints  # noqa: PLC0415

        y0 = pre_solve_constraints(
            spec, grid, y0, bc=bc, parameters=parameters, t=t_span[0]
        )

    # Check that IC satisfies subsidiary constraints (no-self-term equations).
    # These are frozen at zero by IDA, so the original equation becomes a
    # consistency condition on the initial data.
    _check_no_self_term_ic(spec, layout, grid, y0, bc)

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

    # Choose linear solver based on system size:
    # - Dense (LU): fast for small systems, O(N^3) / O(N^2) memory
    # - Sparse (SuperLU_MT): analytical sparsity pattern, O(nnz) memory
    # - GMRES: matrix-free, O(krylov_dim * N) memory, no preconditioner
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

    result: SundialsResult = call_ida(resfn, t_eval, y0, yp0, **options)

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
