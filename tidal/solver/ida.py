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

from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver.fields import FieldSet
from tidal.solver.operators import apply_operator
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from tidal.solver.grid import GridInfo
    from tidal.solver.rhs import RHSEvaluator
    from tidal.solver.state import SlotInfo
    from tidal.symbolic.json_loader import EquationSystem

# Time-derivative order threshold for dynamical (wave) equations
_SECOND_ORDER = 2

# System size threshold for switching from dense to iterative linear solver
_DENSE_THRESHOLD = 100_000


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
        bc: str | tuple[str, ...] | None,
        kinetic: np.ndarray | None,
        spatial_momenta: dict | None,
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
        # Per-call state (set via set_arrays)
        self.y: np.ndarray = np.empty(0)
        self.yp: np.ndarray = np.empty(0)
        self.res: np.ndarray = np.empty(0)
        self.fieldset: FieldSet | None = None
        # Legacy dict for kinetic/field_rates handlers that still use raw arrays
        self.fields: dict[str, np.ndarray] = {}

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

    def handle_constraint(self, slot_idx: int, slot: SlotInfo) -> None:
        """Algebraic constraint: RHS = 0."""
        s = slice(slot_idx * self.n, (slot_idx + 1) * self.n)
        eq_idx = self.eq_map.get(slot.field_name)
        if eq_idx is None:
            self.res[s] = 0.0
        else:
            self.res[s] = self.compute_rhs(eq_idx).ravel()

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
        """Fallback: use field_rates when K not available."""
        canonical = self.spec.canonical
        assert canonical is not None
        rate = np.zeros(self.shape)
        for term in canonical.field_rates[field_name]:
            target_data = self.fields.get(term.field, np.zeros(self.shape))
            operated = apply_operator(term.operator, target_data, self.grid, self.bc)
            rate += term.coefficient * operated
        self.res[s] = self.yp[s] - rate.ravel()

    def _handle_identity_k(self, s: slice, field_name: str) -> None:
        """Identity K: dq/dt = pi."""
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
    bc: str | tuple[str, ...] | None = None,
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


def solve_ida(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: str | tuple[str, ...] | None = None,
    parameters: dict[str, float] | None = None,
    num_snapshots: int = 101,
    rtol: float = 1e-8,
    atol: float = 1e-10,
    max_steps: int = 50000,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    calc_initcond: str | None = None,
) -> dict[str, Any]:
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
    """
    from sksundae.ida import IDA  # noqa: PLC0415

    layout = StateLayout.from_spec(spec, grid.num_points)
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
        options["calc_initcond"] = calc_initcond or "yp0"
        options["calc_init_dt"] = float(t_eval[1] - t_eval[0])

    # Choose linear solver based on system size
    if layout.total_size <= _DENSE_THRESHOLD:
        options["linsolver"] = "dense"
    else:
        options["linsolver"] = "gmres"

    solver = IDA(resfn, **options)
    result = solver.solve(t_eval, y0, yp0)

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
