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

        * ``"validity_param"`` — ``max_k( ε · ω · t_end )`` across all
          small parameters and simulated modes. The dimensionless
          secular-phase-error bound. < 0.1 is safe, > 1 breaks the EFT.
          R4.1 / #284 corrected this from the previous (over-counting)
          ``ε·ω²·t`` formula.
        * ``"omega_max"`` — maximum simulated frequency ``max|λ|`` (so
          damped/tachyonic modes register at their full magnitude).
        * ``"warn_level"`` — ``"ok"``, ``"warn"``, or ``"error"``.
          The overall level is the worse of ``correction_level`` and
          ``base_level``.
        * ``"correction_level"`` — warn band for the secular-error
          side only (R4.1 thresholds: 0.1, 1.0).
        * ``"base_stability_param"`` — ``max Re(λ) · t_end``. A
          strictly-positive value means the base theory has an
          unstable mode (tachyon, Jeans, ghost). > 30 → overflow;
          > 1e-10 → warn-level R4.2 / #285.
        * ``"base_level"`` — warn band for the base-stability side.
        * ``"max_real_lambda"`` — ``max Re(λ)`` across all blocks.
        * ``"dominant_tachyon_field"`` — field/velocity slot owning
          the most-positive ``Re(λ)`` (or None if base is stable).
        * ``"correction_drops"`` — list of diagnostic records for any
          correction terms or equations that ``_build_source_matrix_k``
          could not route (e.g. target field has no dynamical slot,
          or equation is for a demoted constraint field — see #272).
          Empty list on a correct pipeline for every shipped JSON at
          shipped parameters; any non-empty list is a signal to
          investigate before trusting Pass 1 magnitudes.
    spec : EquationSystem | None
        The **base spec** (post-Gap-B demotion) whose layout matches
        ``total["y"]`` and each ``orders[n]["y"]``. Callers that hand
        ``total`` to :meth:`SimulationData.from_solver_result` must use
        this spec, NOT the unreduced full spec — the demoted layout
        has a different slot ordering and a reduced component set.
        Pre-R2.2 the CLI silently rebound ``spec = solver.base_spec``
        after the solve; the explicit attribute removes that workaround
        (#276).
    full_spec : EquationSystem | None
        The unreduced spec as loaded from JSON (pre-Gap-B). Kept for
        diagnostic purposes — e.g. inspecting original ``order_in_eps``
        tags or reconstructing the promoted field's kinetic
        coefficient.
    """

    orders: list[dict[str, Any]] = field(default_factory=list)
    total: dict[str, Any] = field(default_factory=dict)
    validity: dict[str, Any] = field(default_factory=dict)
    spec: EquationSystem | None = None
    full_spec: EquationSystem | None = None


def _compute_validity(
    eigendata: dict[str, Any],
    small_parameters: list[str],
    parameters: dict[str, float] | None,
    t_end: float,
) -> dict[str, Any]:
    """Compute the iterative-expansion validity parameter.

    Secular phase accumulation bound (R4.1 / #284). For a mass-shift
    correction ``-ε·φ`` on a base oscillator ``ω₀² = m² + k²``:

    * Fractional frequency shift: ``δω/ω ≈ ε/(2ω²)``.
    * Accumulated phase error over N oscillations (N ≈ ω·t_end/2π):
      ``Δφ ≈ δω · t_end = ε·t_end/(2ω)``.
    * Scaled to a dimensionless "fraction of period" for threshold
      comparison: the relevant bound is ``ε · ω · t_end``, NOT
      ``ε · ω² · t_end``. Pre-R4.1 the driver used the wrong formula,
      which was conservative (over-triggers) for ``ω > 1`` and
      dangerous (under-triggers) for ``ω < 1`` — low-frequency modes
      could see real ``ε·t`` ≈ 1 while the reported score stayed
      below 0.1.

    ``ω_max`` is estimated from the largest absolute eigenvalue across
    all modes and blocks (so damped/tachyonic modes contribute their
    full magnitude via ``|λ|`` rather than ``|Im(λ)|``). R4.2 / #285
    also computes a base-stability diagnostic from ``max Re(λ) · t_end``.
    """
    omega_max = 0.0
    max_real_lambda = 0.0
    dominant_tachyon_block: int | None = None
    dominant_tachyon_slot: int | None = None
    for block_idx, block in enumerate(eigendata.get("blocks", [])):
        lam = block["D_diag"]  # (n_modes, bs) complex
        # |λ| so damped/tachyonic modes (non-zero Re(λ)) register at
        # their full magnitude.
        omega_block = float(np.max(np.abs(lam)))
        omega_max = max(omega_max, omega_block)

        # R4.2 / #285: track max Re(λ) per block for base-theory
        # stability. A positive Re(λ) means the base Pass 0 solution
        # grows exponentially (tachyon / unstable mode); no Pass 1
        # correction can be meaningful on top of it.
        real_lam = np.real(lam)
        block_max_real = float(np.max(real_lam))
        if block_max_real > max_real_lambda:
            max_real_lambda = block_max_real
            # Locate the (mode, slot) where the growth lives.
            flat_idx = int(np.argmax(real_lam))
            n_modes_block = real_lam.shape[0]
            bs = real_lam.shape[1]
            dominant_tachyon_block = block_idx
            dominant_tachyon_slot = int(
                block["slot_indices"][flat_idx % bs]  # slot within this block
            )
            _ = n_modes_block  # silence unused-var lint; retained for clarity

    params = dict(parameters or {})
    eps_values: dict[str, float] = {}
    for name in small_parameters:
        if name in params:
            eps_values[name] = float(params[name])

    validity_param = 0.0
    dominant: str | None = None
    for name, val in eps_values.items():
        # R4.1 / #284: ε·ω·t_end (secular phase-error bound), not
        # ε·ω²·t_end.
        score = abs(val) * omega_max * float(t_end)
        if score > validity_param:
            validity_param = score
            dominant = name

    # R4.2 / #285: base-theory stability parameter. exp(30) ≈ 1e13 so
    # any base_growth > 30 means Pass 0 has already overflowed double
    # precision within the simulation window.
    base_growth = max_real_lambda * float(t_end)
    # Warning bands on the secular-error side.
    if validity_param > _VALIDITY_ERROR_THRESHOLD:
        correction_level = "error"
    elif validity_param > _VALIDITY_WARN_THRESHOLD:
        correction_level = "warn"
    else:
        correction_level = "ok"

    # Base-theory stability side. Thresholds: 30 = double-precision
    # overflow; 1e-10 = numerical noise (floating-point residue of a
    # nominally zero eigenvalue). Anything in between is a real but
    # possibly physical instability (e.g. Jeans, tachyonic torsion)
    # that needs user review before trusting Pass 1.
    if base_growth > 30.0:
        base_level = "error"
    elif max_real_lambda > 1e-10:
        base_level = "warn"
    else:
        base_level = "ok"

    # Overall level = worst of the two diagnostics.
    level_rank = {"ok": 0, "warn": 1, "error": 2}
    warn_level = max((correction_level, base_level), key=lambda s: level_rank[s])

    # Locate the tachyonic field name if one was found.
    dominant_tachyon_field: str | None = None
    if dominant_tachyon_slot is not None:
        layout = eigendata.get("state_layout")
        if layout is not None:
            slot_to_name: dict[int, str] = {}
            for name, slot in layout.field_slot_map.items():
                slot_to_name[slot] = f"{name} (field slot)"
            for name, slot in layout.velocity_slot_map.items():
                slot_to_name[slot] = f"v_{name} (velocity slot)"
            dominant_tachyon_field = slot_to_name.get(
                dominant_tachyon_slot, f"slot {dominant_tachyon_slot}"
            )

    return {
        "omega_max": omega_max,
        "eps_values": eps_values,
        "validity_param": validity_param,
        "dominant_parameter": dominant,
        "warn_level": warn_level,
        "base_stability_param": base_growth,
        "max_real_lambda": max_real_lambda,
        "dominant_tachyon_field": dominant_tachyon_field,
        "dominant_tachyon_block": dominant_tachyon_block,
        "correction_level": correction_level,
        "base_level": base_level,
    }


def _assemble_full_state_pass_n(
    pass_n_result: dict[str, Any],
    eigendata: dict[str, Any],
    full_state_size: int,
    full_layout: Any,
    grid: GridInfo,
    constraint_source_hat: NDArray[np.complex128] | None = None,
) -> NDArray[np.float64]:
    r"""Assemble the Pass-n full-layout state from the reduced-layout output.

    Pass n evolves only the dynamical (Schur-reduced) subspace. This
    function reconstructs the full-layout state in physical space:

    * Dynamical slots are mapped from the reduced layout back to their
      original slot indices using ``eigendata["schur_ops"]
      ["orig_to_reduced"]``.
    * Constraint slots are populated via the **augmented** Schur
      recovery (v6 R8 / #290):

      .. math::

        \\hat h_c^{(1)} = \\mathrm{recovery} \\cdot \\hat y_{dyn}^{(1)}
                     + S_{cc}^{-1} \\cdot \\mathrm{source\\_hat}_c

      where ``source_hat_c = K·d^n_t(h_c⁰) − corr_1(y⁰)`` (see
      :func:`_compute_constraint_source_hat`). Pre-R8, only the
      first term was applied — see #290 for the ~14-orders-of-
      magnitude discrepancy this caused on the flagship R̃² PGT
      theory.

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
    # R8 / #290: S_cc_inv is required when augmentation is applied.
    scc_inv: NDArray[np.complex128] | None = schur_ops.get("S_cc_inv")

    n_snap = pass_n_reduced.shape[0]
    pass_n_full = np.zeros((n_snap, full_state_size))

    # --- dynamical slots: direct transfer from reduced layout ---
    for orig_si, red_pos in orig_to_reduced.items():
        src = pass_n_reduced[:, red_pos * num_points : (red_pos + 1) * num_points]
        pass_n_full[:, orig_si * num_points : (orig_si + 1) * num_points] = src

    # --- constraint slots: Schur-base recovery from y_hat_dyn ---
    # y_hat_dyn is a required field on PerturbativePass1Result (#279);
    # its absence here is a programming error, not a runtime condition.
    assert "y_hat_dyn" in pass_n_result, (
        "PerturbativePass1Result missing required 'y_hat_dyn' key. "
        "See #279: this is a contract violation by solve_modal_pass1, "
        "not a runtime fallback."
    )
    y_hat_dyn = pass_n_result["y_hat_dyn"]
    if not c_names:
        # Legitimate case: the base spec has no constraint fields, so
        # there is nothing to Schur-recover. Constraint slots stay at
        # zero (there aren't any in the full layout either).
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

        # R8 / #290: augmented Schur recovery. Adds
        # S_cc_inv · source_hat(ti) where source_hat encodes the
        # LHS-feedback K·d^n_t(h_c⁰) and the order-1 RHS corrections
        # on demoted constraint rows. Without this, Pass 1 for any
        # theory where b5 corrections live on constraint rows
        # (e.g. R̃² PGT) is identically zero on the physics-relevant
        # fields — see #290 for the full derivation.
        if constraint_source_hat is not None and scc_inv is not None:
            # source_hat[ti]: (n_c, n_modes)
            # scc_inv: (n_modes, n_c, n_c)
            # augmented_hat[c_idx, m] = sum_c' S_cc_inv[m, c_idx, c'] · source_hat[ti, c', m]
            augmented_hat = np.einsum("mij,jm->im", scc_inv, constraint_source_hat[ti])
            c_hat += augmented_hat

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
        self._small_parameters = small_parameters
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
        NotImplementedError
            If ``order >= 2`` — Pass 2 against ``q¹(t)`` is not yet
            implemented. See #273.
        """
        if order > self._max_order:
            msg = (
                f"Requested order={order} but the spec carries terms only "
                f"up to order={self._max_order}. Re-derive the theory with "
                f"a larger truncation or lower --perturbative-order."
            )
            raise ValueError(msg)
        if order >= 2:
            msg = (
                "PerturbativeSolver.solve() currently supports order=0 "
                "and order=1 only. order >= 2 requires Pass 2 Duhamel "
                "against q¹(t), which needs a new IC-amplitude contract "
                "and is not yet implemented. Tracked in #273."
            )
            raise NotImplementedError(msg)

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
        correction_drops: list[dict[str, Any]] = []

        if order >= 1 and self.has_corrections():
            from tidal.solver.state import StateLayout  # noqa: PLC0415

            eigendata = pass0["eigendata"]
            full_state_size = pass0["y"].shape[1]
            full_layout = StateLayout.from_spec(self.base_spec, grid.num_points)

            # R8 / #290: compute once the constraint-row source
            # (K · d^n_t(h_c⁰) − corr_1(y⁰)) that augments the Pass 1
            # Schur recovery. Depends only on Pass 0 output + spec, so
            # it's shared across all correction orders.
            pre_demote = _pre_demote_info(
                self.full_spec, self.base_spec, self._small_parameters
            )
            constraint_source_hat = _compute_constraint_source_hat(
                full_spec=self.full_spec,
                eigendata=eigendata,
                pre_demote_info=pre_demote,
                parameters=parameters,
                pass0_t=pass0["t"],
                grid=grid,
            )

            for n in range(1, order + 1):
                correction_spec = self.full_spec.filter_by_order(n)
                if not correction_spec.has_corrections():
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
                correction_drops.extend(
                    {**d, "pass": n} for d in pass_n["correction_drops"]
                )
                # Pass 1 (n=1) uses the augmented Schur recovery. Higher
                # orders don't have an analogous augmentation derived yet
                # (tracked in #273); current code gates order >= 2.
                pass_n_full = _assemble_full_state_pass_n(
                    pass_n,
                    eigendata,
                    full_state_size=full_state_size,
                    full_layout=full_layout,
                    grid=grid,
                    constraint_source_hat=constraint_source_hat if n == 1 else None,
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
        # v6 R1.1 (#272): surface all skipped correction terms so the
        # CLI / caller can see when Pass 1 silently misses contributions.
        validity["correction_drops"] = correction_drops

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
        return PerturbativeResult(
            orders=orders,
            total=total,
            validity=validity,
            spec=self.base_spec,
            full_spec=self.full_spec,
        )


def _pre_demote_info(
    full_spec: EquationSystem,
    base_spec: EquationSystem,
    small_parameters: list[str],
) -> dict[str, tuple[str, int]]:
    """Return per-demoted-field pre-demotion metadata.

    v6 R8.2 / #290: the augmented Pass 1 constraint recovery needs the
    symbolic kinetic coefficient and the original LHS time-derivative
    order for every demoted field whose pre-demote kinetic is
    non-zero (e.g. R̃² PGT's ``h_4`` with kinetic ``"2*b5"``).

    Symbolic form (not pre-numeric) is retained because the correct
    runtime contribution is ``evaluate(kinetic, runtime_params)``,
    which naturally handles both linear (``"2*b5"`` → 2·ε) and
    higher-order (``"b5**2"`` → ε², which is O(ε²) and therefore
    vanishes at Pass 1) kinetics. The driver resolves at solve-time
    using :func:`evaluate_with_substitutions`.

    Returns
    -------
    dict[str, (kinetic_symbolic, n)]
        Maps demoted field name → (kinetic coefficient symbolic
        expression, pre-demote time-derivative order). Validation
        (that the kinetic evaluates to 0 at small_params=0 i.e. is
        actually small-parameter-dependent) happens at the call
        site.
    """
    from tidal.symbolic._kinetic_eval import lhs_collapses_to_zero  # noqa: PLC0415

    result: dict[str, tuple[str, int]] = {}
    demoted_names = {
        eq.field_name for eq in base_spec.equations if eq.time_derivative_order == 0
    }
    for eq in full_spec.equations:
        if eq.field_name not in demoted_names:
            continue
        if eq.kinetic_coefficient_symbolic is None:
            # Pre-existing algebraic constraint; no LHS feedback.
            continue
        if eq.time_derivative_order <= 0:
            continue
        # Gate: only fields whose kinetic actually collapses at
        # small_params=0 contribute LHS feedback at O(ε¹). Non-
        # collapsing kinetics stayed at time_order > 0 in base_spec
        # anyway (not demoted) — so they wouldn't be in demoted_names.
        if not lhs_collapses_to_zero(eq.kinetic_coefficient_symbolic, small_parameters):
            continue
        result[eq.field_name] = (
            eq.kinetic_coefficient_symbolic,
            int(eq.time_derivative_order),
        )
    return result


def _compute_constraint_source_hat(  # noqa: C901, PLR0912, PLR0914, PLR0915
    *,
    full_spec: EquationSystem,
    eigendata: dict[str, Any],
    pre_demote_info: dict[str, tuple[str, int]],
    parameters: dict[str, float] | None,
    pass0_t: NDArray[np.float64],
    grid: GridInfo,
) -> NDArray[np.complex128]:
    """Compute the augmented constraint-source term for Pass 1 recovery.

    v6 R8.3 / #290. The order-1 equation for a demoted constraint field
    h_c (derived from ``ε·K·d^n_t(h_c) = S_cc·h_c + S_cd·y_dyn +
    ε·corr_1(y)``) yields::

        h_c¹ = recovery · y_dyn¹
             + S_cc⁻¹ · [K · d^n_t(h_c⁰) − corr_1(y⁰)]

    This helper produces the bracketed part (``source_hat`` below) per
    snapshot, per constraint row, per Fourier mode. The caller applies
    ``S_cc_inv`` and adds ``recovery · y_hat_dyn¹`` to get h_c¹.

    Parameters
    ----------
    full_spec
        Unreduced :class:`EquationSystem` (pre-Gap-B demotion). Carries
        the original kinetic coefficients + the order-1 RHS corrections
        on demoted rows.
    eigendata
        Output of :func:`solve_modal` with ``return_eigendata=True`` +
        ``schur_ops`` populated. Must include ``recovery_matrix``.
    pre_demote_info
        Output of :func:`_pre_demote_info`. Maps demoted-field name →
        ``(K_numeric, n)``.
    parameters
        Runtime parameter overrides for resolving coefficient symbols.
    pass0_t
        Snapshot times (shape ``(n_snapshots,)``).
    grid
        Spatial grid; used for Fourier-multiplier evaluation.

    Returns
    -------
    source_hat : complex ndarray
        Shape ``(n_snapshots, n_constraints, n_modes)``. Entry
        ``[ti, c, m]`` carries ``K·d^n_t(h_c⁰)_m − corr_1(y⁰)_m`` at
        time ``pass0_t[ti]``.
    """
    # pylint: disable=import-outside-toplevel
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
    from tidal.solver.modal import (  # noqa: PLC0415
        _OPERATOR_DECOMP,
        _build_k_axes,
        _build_k_grid,
        _resolve_constant_coeff,
    )

    schur_ops = eigendata.get("schur_ops")
    if schur_ops is None or "recovery_matrix" not in schur_ops:
        # No constraints → nothing to source. Caller should handle this
        # earlier; return an empty-shape array as a defensive default.
        n_snap = pass0_t.shape[0]
        rfft_last = grid.shape[-1] // 2 + 1
        n_modes = int(np.prod([*grid.shape[:-1], rfft_last]))
        return np.zeros((n_snap, 0, n_modes), dtype=np.complex128)

    c_names: tuple[str, ...] = schur_ops["constraint_field_names"]
    recovery_matrix: NDArray[np.complex128] = schur_ops["recovery_matrix"]
    # (n_modes, n_c, n_dyn)
    n_modes = recovery_matrix.shape[0]
    n_c = recovery_matrix.shape[1]
    n_dyn = recovery_matrix.shape[2]
    n_snap = pass0_t.shape[0]

    c_idx_map: dict[str, int] = {name: i for i, name in enumerate(c_names)}

    # Reduced dynamical layout (from eigendata).
    dyn_layout = eigendata["state_layout"]

    # --- Fourier k-grid and spatial-multiplier cache --------------------
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    spatial_cache: dict[str, NDArray[np.complex128]] = {}

    def _spatial_multiplier(op: str) -> NDArray[np.complex128]:
        if op in spatial_cache:
            return spatial_cache[op]
        decomp = _OPERATOR_DECOMP.get(op)
        if decomp is None:
            msg = (
                f"_compute_constraint_source_hat: operator {op!r} has no "
                f"spatial multiplier entry in _OPERATOR_DECOMP"
            )
            raise ValueError(msg)
        mult_val = decomp.spatial_fn(k_grid)
        mult_full = np.broadcast_to(mult_val, rfft_shape)
        out = mult_full.ravel().astype(np.complex128)
        spatial_cache[op] = out
        return out

    # --- d^n_t of Pass 0 dynamical state in Fourier space ---------------
    # Collect the set of time-derivative orders that will actually be
    # queried. This INCLUDES the per-demoted-field LHS-feedback orders
    # (from pre_demote_info) AND every order-1 RHS term's operator
    # time_order on demoted rows.
    needed_orders: set[int] = {n for (_kin, n) in pre_demote_info.values()}
    for eq in full_spec.equations:
        if eq.field_name not in c_idx_map:
            continue
        for term in eq.rhs_terms:
            if term.order_in_eps != 1:
                continue
            decomp = _OPERATOR_DECOMP.get(term.operator)
            if decomp is None:
                continue
            needed_orders.add(int(decomp.time_order))

    # Precompute d^n_t(y_dyn⁰)(snap, slot, mode) for every n needed.
    # Shape: (n_snap, n_dyn, n_modes), per n.
    y_hat_dyn_dnt: dict[int, NDArray[np.complex128]] = {}
    for n in needed_orders:
        arr = np.zeros((n_snap, n_dyn, n_modes), dtype=np.complex128)
        for block in eigendata["blocks"]:
            slots = list(block["slot_indices"])
            V = block["V"]  # (n_modes, bs, bs)
            lam = block["D_diag"]  # (n_modes, bs)
            alpha = block["alpha"]  # (n_modes, bs)
            # (n_modes, bs) — λ^n per mode, per eigenmode.
            lam_pow_n = lam**n if n else np.ones_like(lam)
            for ti, t in enumerate(pass0_t):
                # exp_lt[m, i] = exp(λ_i · t) per mode.
                exp_lt = np.exp(lam * float(t))
                # coeff_i = λ_i^n · exp(λ_i·t) · α_i     (n_modes, bs)
                coeff = lam_pow_n * exp_lt * alpha
                # y_hat per slot = V · coeff   shape (n_modes, bs)
                y_hat = np.einsum("mij,mj->mi", V, coeff)
                # Scatter into the (n_dyn,)-indexed array.
                for slot_local, slot_idx in enumerate(slots):
                    arr[ti, slot_idx, :] += y_hat[:, slot_local]
        y_hat_dyn_dnt[n] = arr

    # Helper: d^n_t(h_c⁰)(ti, m) for a constraint target at time ti.
    def _constraint_dnt_hat(n: int, c_idx: int, ti: int) -> NDArray[np.complex128]:
        # h_c⁰_dnt = recovery[m, c_idx, :] · y_hat_dyn_dnt[n][ti, :, m]
        # einsum: "mj, jm -> m"
        return np.einsum("mj,jm->m", recovery_matrix[:, c_idx, :], y_hat_dyn_dnt[n][ti])

    # --- Coefficient evaluator for symbolic corrections ------------------
    coeff_eval = CoefficientEvaluator(full_spec, grid, parameters or {})

    # --- Build source_hat ------------------------------------------------
    source_hat = np.zeros((n_snap, n_c, n_modes), dtype=np.complex128)

    # (a) LHS-feedback: source += +K_eff · d^n_t(h_c⁰) for each demoted
    # field. K_eff is the kinetic evaluated at the RUNTIME parameter
    # values (NOT at small_params=1), so the ε factor is correctly
    # baked in for linear-in-ε kinetics and any higher-order kinetic
    # naturally vanishes at O(ε¹).
    from tidal.symbolic._kinetic_eval import (  # noqa: PLC0415
        evaluate_with_substitutions,
    )

    param_values = {k: float(v) for k, v in (parameters or {}).items()}
    for c_name, (kinetic_symbolic, n) in pre_demote_info.items():
        c_idx = c_idx_map.get(c_name)
        if c_idx is None:
            continue
        K_eff = evaluate_with_substitutions(kinetic_symbolic, param_values)
        if K_eff is None or K_eff == 0.0:
            # Missing small-parameter value at runtime or it's zero —
            # no LHS feedback at this parameter point.
            continue
        for ti in range(n_snap):
            source_hat[ti, c_idx, :] += float(K_eff) * _constraint_dnt_hat(n, c_idx, ti)

    # (b) Explicit order-1 RHS corrections on demoted constraint rows.
    # Sign: the derivation puts ``-corr_1(y⁰)`` on the source side. So
    # we SUBTRACT each evaluated term from source_hat.
    for eq_idx, eq in enumerate(full_spec.equations):
        c_idx = c_idx_map.get(eq.field_name)
        if c_idx is None:
            continue  # row is dynamical, not demoted
        for term_idx, term in enumerate(eq.rhs_terms):
            if term.order_in_eps != 1:
                continue
            # Resolve numeric coefficient via modal's existing helper.
            coeff = _resolve_constant_coeff(
                term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx
            )
            mult = _spatial_multiplier(term.operator)  # (n_modes,)
            decomp = _OPERATOR_DECOMP[term.operator]
            n_op = int(decomp.time_order)

            # Resolve target kind.
            target_field = term.field
            # Velocity reference (``v_<field>``): strip prefix for lookup.
            is_velocity_ref = target_field.startswith("v_")
            base_target = target_field[2:] if is_velocity_ref else target_field

            if base_target in dyn_layout.field_slot_map and not is_velocity_ref:
                slot_idx = dyn_layout.field_slot_map[base_target]
                for ti in range(n_snap):
                    source_hat[ti, c_idx, :] -= (
                        coeff * mult * y_hat_dyn_dnt[n_op][ti, slot_idx, :]
                    )
            elif is_velocity_ref and base_target in dyn_layout.velocity_slot_map:
                slot_idx = dyn_layout.velocity_slot_map[base_target]
                for ti in range(n_snap):
                    source_hat[ti, c_idx, :] -= (
                        coeff * mult * y_hat_dyn_dnt[n_op][ti, slot_idx, :]
                    )
            elif base_target in c_idx_map:
                # Constraint target: expand via base Schur.
                target_c_idx = c_idx_map[base_target]
                for ti in range(n_snap):
                    source_hat[ti, c_idx, :] -= (
                        coeff * mult * _constraint_dnt_hat(n_op, target_c_idx, ti)
                    )
            # else: unresolvable target — the driver's _build_source_matrix_k
            # already logs this as a column-missing drop. Skip here; the
            # existing diagnostic surfaces the issue.

    return source_hat
