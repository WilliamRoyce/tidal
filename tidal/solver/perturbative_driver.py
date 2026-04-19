"""Iterative perturbative driver (v6 plan, Stage 5).

Given an :class:`~tidal.symbolic.json_loader.EquationSystem` whose terms
carry an ``order_in_eps`` annotation (emitted by ExportJSON.wl when the
theory has a ``[perturbation]`` section), this module composes the
modal solver's Pass 0 base evolution with one or more Pass 1 Duhamel
corrections into a single :class:`PerturbativeResult`.

See :doc:`/home/vscode/.claude/plans/flickering-gathering-orbit.md`
for the full architecture. Key guarantees:

- Ghost-free by construction: the LHS operator at every order is the
  base 2nd-order operator; higher-derivative terms appear only as
  sources on the previous order's solution.
- Theory-agnostic: driven purely by ``order_in_eps`` tags and
  operator strings. No per-theory classification.
- Reuses Pass 0's eigendecomposition for every Pass 1+ call, so the
  added cost is ~1.5x a single Pass 0 at machine precision.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, cast

import numpy as np

from tidal.solver.modal import solve_modal, solve_modal_pass1

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem

# Validity-monitor thresholds. See v6 plan §"TIDAL's parameter space":
#   validity_param = max( |eps| · omega_max^2 · t_end )
#   < 0.1  → "ok"
#   > 0.1  → "warn" (truncation error may exceed 10%)
#   > 1.0  → "error" (EFT regime broken; results not trustworthy)
_VALIDITY_WARN_THRESHOLD = 0.1
_VALIDITY_ERROR_THRESHOLD = 1.0


@dataclass
class PerturbativeResult:
    """Output of :meth:`PerturbativeSolver.solve`.

    Attributes
    ----------
    orders : list[dict]
        Per-order ``SolverResult``-shaped dicts. ``orders[0]`` is Pass 0
        (base), ``orders[n]`` is the Pass-n correction for ``n ≥ 1``.
    total : dict
        Combined solution ``Σ orders[n].y`` in physical space, aligned
        with ``orders[0].t``.
    validity : dict
        Diagnostic metrics from :func:`_compute_validity`:

        * ``"validity_param"`` — max_k( ε · ω² · t_end ) across all
          small parameters and simulated modes. < 0.1 is safe, 1 is
          the EFT-breakdown boundary.
        * ``"omega_max"`` — maximum simulated frequency.
        * ``"warn_level"`` — ``"ok"``, ``"warn"`` (> 0.1), or ``"error"``
          (> 1.0).
    """

    orders: list[dict[str, Any]] = field(default_factory=list)
    total: dict[str, Any] = field(default_factory=dict)
    validity: dict[str, Any] = field(default_factory=dict)


def _compute_validity(
    eigendata: dict[str, Any],
    small_parameters: list[str],
    parameters: dict[str, float] | None,
    t_end: float,
) -> dict[str, Any]:
    """Compute the iterative-expansion validity parameter.

    See plan §"Validity condition in detail". The dominant error of the
    O(ε¹) truncation accumulates linearly in time for TIDAL's linear
    regime, yielding the criterion ``max(ε · ω² · t_end) ≪ 1``.

    The maximum simulation frequency ``ω_max`` is estimated from the
    largest imaginary eigenvalue across all modes and blocks — this
    captures both the dispersion relation and any cross-coupling
    mass shifts present in the base theory.
    """
    omega_max = 0.0
    for block in eigendata.get("blocks", []):
        lam = block["D_diag"]  # (n_modes, bs) complex
        # Use |λ| so damped/tachyonic modes (non-zero Re(λ)) register
        # at their full magnitude.  Previously |Im(λ)| underestimated
        # ω_max for theories with dissipation or with near-zero-mass
        # fields where Re(λ) carries significant amplitude.
        omega_block = float(np.max(np.abs(lam)))
        omega_max = max(omega_max, omega_block)

    params = dict(parameters or {})
    eps_values: dict[str, float] = {}
    for name in small_parameters:
        if name in params:
            eps_values[name] = float(params[name])

    validity_param = 0.0
    dominant: str | None = None
    for name, val in eps_values.items():
        score = abs(val) * omega_max**2 * float(t_end)
        if score > validity_param:
            validity_param = score
            dominant = name

    # Warning bands. Thresholds match the plan's §"TIDAL's parameter
    # space" table: WARN_THRESHOLD = 10% accumulated error,
    # ERROR_THRESHOLD = full EFT breakdown.
    if validity_param > _VALIDITY_ERROR_THRESHOLD:
        warn_level = "error"
    elif validity_param > _VALIDITY_WARN_THRESHOLD:
        warn_level = "warn"
    else:
        warn_level = "ok"

    return {
        "omega_max": omega_max,
        "eps_values": eps_values,
        "validity_param": validity_param,
        "dominant_parameter": dominant,
        "warn_level": warn_level,
    }


def _assemble_full_state_pass_n(
    pass_n_result: dict[str, Any],
    eigendata: dict[str, Any],
    full_state_size: int,
    full_layout: Any,
    grid: GridInfo,
) -> NDArray[np.float64]:
    """Assemble the Pass-n full-layout state from the reduced-layout output.

    Pass n evolves only the dynamical (Schur-reduced) subspace. This
    function reconstructs the full-layout state in physical space:

    * Dynamical slots are mapped from the reduced layout back to their
      original slot indices using ``eigendata["schur_ops"]
      ["orig_to_reduced"]``.
    * Constraint slots are populated via Schur recovery on the
      *Fourier-space* dynamical Pass n output: ``c_hat = recovery_matrix
      @ y_hat_dyn``. This is the O(ε¹) Schur-base recovery. The
      O(ε²) augmentation from the correction's contribution to the
      constraint equation itself (``S_cc^{-1} · source_c(y⁰)``) is
      dropped at order 1.

    Parameters
    ----------
    pass_n_result : dict
        Output of :func:`tidal.solver.modal.solve_modal_pass1`. Must
        contain ``"y"`` (reduced-layout physical state, shape
        ``(n_snap, n_dyn · n_points)``) and ``"y_hat_dyn"``
        (Fourier-space dynamical state, shape
        ``(n_snap, n_dyn, n_modes)``).
    eigendata : dict
        Pass 0 eigendata. Must include ``"schur_ops"`` when the base
        spec had constraints.
    full_state_size : int
        Target flat size ``full_layout.num_slots * full_layout.num_points``.
    full_layout : StateLayout
        The Pass 0 full state layout (not the reduced dynamical one).
        Used to look up the slot index of each constraint field.
    grid : GridInfo
        Spatial grid (for inverse FFT shape).

    Returns
    -------
    pass_n_full : NDArray[float64]
        Shape ``(n_snap, full_state_size)``. Dynamical and constraint
        slots populated in the full layout.

    Raises
    ------
    ValueError
        If the shapes are inconsistent and no ``schur_ops`` is present
        for expansion.
    """
    pass_n_reduced = pass_n_result["y"]
    schur_ops = eigendata.get("schur_ops")

    if schur_ops is None:
        # No constraints: reduced layout == full layout.
        if pass_n_reduced.shape[1] == full_state_size:
            return np.asarray(pass_n_reduced)
        msg = (
            "Pass n output has shape that does not match the expected full "
            "state size and no schur_ops are available for expansion"
        )
        raise ValueError(msg)

    num_points = grid.num_points
    orig_to_reduced: dict[int, int] = schur_ops["orig_to_reduced"]
    c_names: tuple[str, ...] = schur_ops["constraint_field_names"]
    recovery_matrix = schur_ops["recovery_matrix"]  # (n_modes, n_c, n_dyn)

    n_snap = pass_n_reduced.shape[0]
    pass_n_full = np.zeros((n_snap, full_state_size))

    # --- dynamical slots: direct transfer from reduced layout ---
    for orig_si, red_pos in orig_to_reduced.items():
        src = pass_n_reduced[:, red_pos * num_points : (red_pos + 1) * num_points]
        pass_n_full[:, orig_si * num_points : (orig_si + 1) * num_points] = src

    # --- constraint slots: Schur-base recovery from y_hat_dyn ---
    y_hat_dyn = pass_n_result.get("y_hat_dyn")
    if y_hat_dyn is None or not c_names:
        # Nothing to recover; constraint slots stay at zero. This is
        # only a valid state when the correction spec leaves all
        # constraint fields untouched (rare).
        return pass_n_full

    # rfft output shape
    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    for ti in range(n_snap):
        # y_hat_dyn[ti] is (n_dyn, n_modes)
        # recovery_matrix is (n_modes, n_c, n_dyn)
        # c_hat[c_idx, m] = sum_j recovery_matrix[m, c_idx, j] * y_hat_dyn[j, m]
        c_hat = np.einsum("mcj,jm->cm", recovery_matrix, y_hat_dyn[ti])
        for ci, c_name in enumerate(c_names):
            c_slot = full_layout.field_slot_map[c_name]
            c_phys = np.fft.irfftn(
                c_hat[ci].reshape(rfft_shape),
                s=grid.shape,
                axes=list(range(len(grid.shape))),
            ).ravel()
            pass_n_full[ti, c_slot * num_points : (c_slot + 1) * num_points] = np.real(
                c_phys
            )

    return pass_n_full


class PerturbativeSolver:
    """Orchestrate Pass 0 + Pass n corrections for a tagged EquationSystem.

    Usage
    -----
    >>> solver = PerturbativeSolver(spec)  # doctest: +SKIP
    >>> res = solver.solve(y0, grid, t_span, order=1, parameters={"b5": 0.01})  # doctest: +SKIP
    >>> q_total = res.total["y"]  # combined physical-space trajectory

    The solver is a thin orchestrator over :func:`solve_modal` with
    ``return_eigendata=True`` and :func:`solve_modal_pass1`. It assumes
    the modal backend; non-modal backends are out of scope for the v6
    plan (tracked in a follow-up issue).
    """

    def __init__(self, spec: EquationSystem) -> None:
        self.full_spec = spec
        # Use base_spec(...) so LHS kinetic coefficients that vanish at
        # eps=0 trigger demotion to algebraic constraint (v6 Gap B).
        # Small parameter names come from the metadata.perturbation block
        # emitted by ExportJSON.wl when [perturbation] is configured.
        pert_meta = spec.metadata.get("perturbation", {}) or {}
        small_parameters = list(pert_meta.get("small_parameters", []))
        self.base_spec = spec.base_spec(small_parameters)
        self._max_order = spec.max_order()

    @property
    def max_order(self) -> int:
        """Maximum ``order_in_eps`` available in the spec."""
        return self._max_order

    def has_corrections(self) -> bool:
        """Return True when the spec contains any ``order_in_eps > 0`` terms."""
        return self.full_spec.has_corrections()

    def solve(  # noqa: PLR0913, PLR0914
        self,
        y0: NDArray[np.float64],
        grid: GridInfo,
        t_span: tuple[float, float],
        *,
        order: int = 1,
        parameters: dict[str, float] | None = None,
        num_snapshots: int = 101,
        small_parameters: list[str] | None = None,
    ) -> PerturbativeResult:
        """Solve iteratively to the requested order.

        Parameters
        ----------
        y0 : ndarray
            Flat initial state for Pass 0. Pass 1+ starts at zero.
        grid : GridInfo
            Spatial grid (must satisfy modal-solver eligibility).
        t_span : (float, float)
            ``(t_start, t_end)``.
        order : int, default 1
            Highest order of correction to compute. ``0`` returns just
            the base solution.
        parameters : dict, optional
            Runtime parameter overrides. Small parameters (those
            tracked in ``order_in_eps``) should appear here with their
            physical values so the correction coefficients resolve.
        num_snapshots : int
            Number of output time points.
        small_parameters : list[str], optional
            Names of the small parameters for the validity monitor.
            When absent, read from
            ``spec.metadata.get("perturbation", {}).get("small_parameters", [])``.

        Returns
        -------
        PerturbativeResult
            Per-order SolverResults + the combined trajectory + the
            validity diagnostics.

        Raises
        ------
        ValueError
            If ``order`` exceeds the maximum ``order_in_eps`` carried by
            the spec.
        """
        if order > self._max_order:
            msg = (
                f"Requested order={order} but the spec carries terms only "
                f"up to order={self._max_order}. Re-derive the theory with "
                f"a larger truncation or lower --perturbative-order."
            )
            raise ValueError(msg)

        # Pass 0 — base evolution + eigendata capture for Pass 1.
        pass0 = cast(
            "dict[str, Any]",
            solve_modal(
                self.base_spec,
                grid,
                y0,
                t_span=t_span,
                parameters=parameters,
                num_snapshots=num_snapshots,
                return_eigendata=(order >= 1),
            ),
        )

        orders: list[dict[str, Any]] = [pass0]
        total_y = pass0["y"].copy()

        if order >= 1 and self.has_corrections():
            from tidal.solver.state import StateLayout  # noqa: PLC0415

            eigendata = pass0["eigendata"]
            full_state_size = pass0["y"].shape[1]
            # Pass 0's output is in the layout of base_spec (post-Gap-B
            # demotion). Build that layout once so Pass n's constraint
            # recovery can look up constraint slot indices in the same
            # layout the combined solution lives in.
            full_layout = StateLayout.from_spec(self.base_spec, grid.num_points)

            for n in range(1, order + 1):
                correction_spec = self.full_spec.filter_by_order(n)
                # filter_by_order(n) keeps terms whose order_in_eps is
                # exactly n; if any survive, they are corrections by
                # definition, so has_corrections() on the filtered spec
                # is the right gate.
                if not correction_spec.has_corrections():
                    # No terms at this order; append a zero contribution
                    # for shape consistency.
                    zero_result = {
                        "t": pass0["t"].copy(),
                        "y": np.zeros_like(pass0["y"]),
                        "success": True,
                        "message": f"Pass {n}: no order-{n} terms",
                    }
                    orders.append(zero_result)
                    continue

                pass_n = solve_modal_pass1(
                    eigendata,
                    correction_spec,
                    grid,
                    pass0["t"],
                    parameters=parameters,
                )
                # Assemble the full-layout state by mapping dynamical
                # slots directly and recovering constraint slots via the
                # Pass 0 Schur operator applied to y_hat_dyn.
                pass_n_full = _assemble_full_state_pass_n(
                    pass_n,
                    eigendata,
                    full_state_size=full_state_size,
                    full_layout=full_layout,
                    grid=grid,
                )
                pass_n["y"] = pass_n_full
                orders.append(pass_n)
                total_y += pass_n_full

        # Validity monitor (only meaningful when corrections are present).
        if small_parameters is None:
            pert_meta = self.full_spec.metadata.get("perturbation", {}) or {}
            small_parameters = list(pert_meta.get("small_parameters", []))

        t_end = float(t_span[1] - t_span[0])
        if order >= 1 and self.has_corrections() and "eigendata" in pass0:
            validity = _compute_validity(
                pass0["eigendata"], small_parameters, parameters, t_end
            )
        else:
            validity = {
                "omega_max": 0.0,
                "eps_values": {},
                "validity_param": 0.0,
                "dominant_parameter": None,
                "warn_level": "ok",
            }

        total: dict[str, Any] = {
            "t": pass0["t"],
            "y": total_y,
            "success": True,
            "message": (
                f"Perturbative solver completed (order={order}, "
                f"{len(orders) - 1} corrections, validity "
                f"{validity['warn_level']})"
            ),
        }
        return PerturbativeResult(orders=orders, total=total, validity=validity)
