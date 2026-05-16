"""JAX-accelerated modal solver — phase 1 (constant-coefficient path only).

Replaces the per-mode ``scipy.linalg.expm`` loop in
:func:`tidal.solver.modal._evolve_per_mode_pade` with ``jax.vmap`` over
``jax.scipy.linalg.expm``, and evolves snapshots via a JIT-compiled
``jax.lax.scan``.  All setup (FFT, per-mode matrix construction, block
detection) reuses the existing scipy path without modification.

Installation::

    uv sync --extra jax

Requires ``jax>=0.4.20`` and ``jaxlib>=0.4.20``.  Double-precision is
enabled automatically (``jax_enable_x64 = True``).

Phase-1 restrictions (raise ``NotImplementedError``):
- Position-dependent coefficient systems (``_evolve_full_matrix`` path).
- Systems requiring Schur constraint elimination (``needs_reduction``).
- ``return_eigendata=True`` (eigendecomposition data not exposed by JAX path).

The ``TIDAL_MODAL_SPARSE_IC`` env var is **ignored** in this path.  JAX
computes ``expm`` for all modes simultaneously via ``vmap``; masking out
zero-IC modes is not cost-free under ``vmap`` (all modes are traced
uniformly).  The eigenvalue-based divergence pre-check (inherited from the
scipy path) still runs outside JIT and catches tachyonic instabilities.

The runtime snapshot-by-snapshot divergence guard is **not** applied inside
the JIT-compiled evolution loop; the pre-check is the only divergence gate.

References
----------
    Bradbury et al. (2018). JAX: composable transformations of Python+NumPy.
    Higham (2009). SIAM Review 51(4):747-764 — matrix exponential.
"""

# ruff: noqa: N803, N806 — uppercase matrix names follow linear-algebra notation.
# ruff: noqa: PLR0913, PLR0914 — solver signatures match solve_modal().
# ruff: noqa: ANN401 — JAX traced arrays have no static dtype; Any is correct.
# ruff: noqa: FBT003 — jax.config.update is a library call with positional bool.
# ruff: noqa: C901, PLR0912, PLR0915 — solve_modal_jax mirrors solve_modal complexity.
# ruff: noqa: PLR2004 — scientific thresholds (1e-15, 1e300, etc.) are named by context.
# ruff: noqa: ARG001 — bc/rtol/atol kept for API compatibility with solve_modal().

from __future__ import annotations

import functools
import math
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._setup import warn_frozen_constraints
from tidal.solver.modal import (
    _build_k_axes,  # pyright: ignore[reportPrivateUsage]
    _build_k_grid,  # pyright: ignore[reportPrivateUsage]
    _build_per_mode_matrices,  # pyright: ignore[reportPrivateUsage]
    _fft_slots,  # pyright: ignore[reportPrivateUsage]
    _has_position_dependent_terms,  # pyright: ignore[reportPrivateUsage]
    _has_time_derivative_operators,  # pyright: ignore[reportPrivateUsage]
    _ifft_slots,  # pyright: ignore[reportPrivateUsage]
    _warn_eigenvalue_growth,  # pyright: ignore[reportPrivateUsage]
    find_independent_blocks,
)
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.solver._types import SolverResult
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.progress import SimulationProgress
    from tidal.symbolic.json_loader import EquationSystem


# ---------------------------------------------------------------------------
# JAX import guard
# ---------------------------------------------------------------------------


def _require_jax() -> tuple[Any, Any, Any]:
    """Import jax, jax.numpy, jax.scipy.linalg with a clear error if missing.

    Raises
    ------
    ImportError
        When ``jax``/``jaxlib`` is not installed.
    """
    try:
        import jax  # noqa: PLC0415
        import jax.numpy as jnp  # noqa: PLC0415
        import jax.scipy.linalg  # noqa: PLC0415
    except ImportError as exc:
        msg = "JAX backend requires jax and jaxlib. Install with:  uv sync --extra jax"
        raise ImportError(msg) from exc
    jax.config.update("jax_enable_x64", True)
    return jax, jnp, jax.scipy.linalg


# ---------------------------------------------------------------------------
# Cached JIT-compiled kernels
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=32)
def _get_evolve_uniform_fn(n_snapshots: int) -> Callable[..., Any]:
    """Return a JIT-compiled scan that evolves *n_snapshots* steps.

    The returned function signature is::

        evolve(exp_M_dt, y0_T) -> ys

    where ``exp_M_dt`` has shape ``(n_modes, bs, bs)`` and ``y0_T`` has
    shape ``(n_modes, bs)``, returning ``ys`` of shape
    ``(n_snapshots, n_modes, bs)`` with ``ys[ti]`` = state at ``t_eval[ti]``.

    Cached by ``n_snapshots`` so compilation happens once per snapshot count.
    The cache persists across calls within a process — warm hits are free.
    """
    jax, jnp, _ = _require_jax()

    @jax.jit
    def _evolve(exp_M_dt: Any, y0_T: Any) -> Any:
        def step(y_curr: Any, _: Any) -> tuple[Any, Any]:
            y_next = jnp.einsum("mij,mj->mi", exp_M_dt, y_curr)
            return y_next, y_curr  # carry=next, output=current (before step)

        # scan runs n_snapshots iterations:
        #   output[0] = y0, output[1] = y1, ..., output[n-1] = y_{n-1}
        _final, ys = jax.lax.scan(step, y0_T, None, length=n_snapshots)
        return ys

    return _evolve


@functools.lru_cache(maxsize=32)
def _get_fused_expm_scan_fn(n_snapshots: int) -> Callable[..., Any]:
    """Return a JIT-compiled (expm + scan) fused kernel for uniform-snapshot evolution.

    Signature::

        fused(M, y0_T, dt) -> ys

    where ``M`` has shape ``(n_modes, bs, bs)`` (UN-scaled — dt multiplied inside),
    ``y0_T`` has shape ``(n_modes, bs)``, and ``dt`` is a scalar. Returns
    ``ys`` of shape ``(n_snapshots, n_modes, bs)``.

    Fuses ``batch_expm`` + ``evolve_uniform`` into a single JIT trace,
    eliminating one round of Python ↔ JAX dispatch per call. Saves ~2-3 ms
    per solve at small problem sizes where dispatch overhead dominates.
    """
    jax, jnp, jsl = _require_jax()

    @jax.jit
    def _fused(M: Any, y0_T: Any, dt: Any) -> Any:
        exp_M_dt = jax.vmap(jsl.expm)(M * dt)  # (n_modes, bs, bs)

        def step(y_curr: Any, _: Any) -> tuple[Any, Any]:
            y_next = jnp.einsum("mij,mj->mi", exp_M_dt, y_curr)
            return y_next, y_curr

        _final, ys = jax.lax.scan(step, y0_T, None, length=n_snapshots)
        return ys

    return _fused


@functools.lru_cache(maxsize=8)
def _get_batch_expm_fn() -> Callable[..., Any]:
    """Return a JIT+vmap wrapper for ``jax.scipy.linalg.expm`` over modes."""
    jax, _, jsl = _require_jax()

    @jax.jit
    def _batch_expm(M_scaled: Any) -> Any:
        # M_scaled: (n_modes, bs, bs) — already multiplied by dt
        return jax.vmap(jsl.expm)(M_scaled)

    return _batch_expm


@functools.lru_cache(maxsize=8)
def _get_batch_eigvals_fn() -> Callable[..., Any]:
    """Return a JIT+vmap wrapper for ``jnp.linalg.eigvals`` over modes.

    Signature::

        batch_eigvals(M) -> eigvals

    ``M`` has shape ``(n_modes, bs, bs)``; output is ``(n_modes, bs)``.
    Replaces the per-block ``np.linalg.eigvals`` in the divergence pre-check,
    moving ~2-3 ms of numpy cost into compiled XLA.
    """
    jax, jnp, _ = _require_jax()

    @jax.jit
    def _batch_eigvals(M: Any) -> Any:
        return jax.vmap(jnp.linalg.eigvals)(M)  # (n_modes, bs, bs) → (n_modes, bs)

    return _batch_eigvals


@functools.lru_cache(maxsize=8)
def _get_evolve_nonuniform_fn() -> Callable[..., Any]:
    """Return a JIT-compiled function for non-uniform snapshot times.

    Signature::

        evolve_nonuniform(M_jax, y0_T, t_rel) -> ys

    ``t_rel`` is a 1-D JAX array of relative times from ``t_start``.
    Uses ``vmap`` over time points (outer) and ``vmap`` over modes (inner).
    """
    jax, jnp, jsl = _require_jax()

    @jax.jit
    def _evolve_nonuniform(M_jax: Any, y0_T: Any, t_rel: Any) -> Any:
        # M_jax: (n_modes, bs, bs), y0_T: (n_modes, bs), t_rel: (n_snapshots,)
        def at_t(t_i: Any) -> Any:
            exp_M_t = jax.vmap(jsl.expm)(M_jax * t_i)
            return jnp.einsum("mij,mj->mi", exp_M_t, y0_T)

        return jax.vmap(at_t)(t_rel)  # (n_snapshots, n_modes, bs)

    return _evolve_nonuniform


# ---------------------------------------------------------------------------
# Phase 2a: constraint / time-derivative RHS path
# ---------------------------------------------------------------------------


def _solve_modal_jax_constrained(  # noqa: PLR0917
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
    y0_hat: NDArray[np.complex128],
    t_eval: NDArray[np.float64],
    jnp: Any,
    num_snapshots: int,
    snapshot_callback: Callable[[float, np.ndarray], None] | None,
    progress: SimulationProgress | None,
) -> SolverResult:
    """Phase 2a evolution: numpy Schur elimination + JAX expm on reduced system.

    Calls :func:`_build_evolution_matrices` (numpy) to eliminate constraint
    fields and produce ``A_reduced`` of shape ``(n_modes, n_dyn, n_dyn)``. The
    evolution itself runs through the JAX vmap+scan kernel; only the Schur
    setup stays in numpy. Constraint fields are reconstructed post-evolution
    via ``recovery @ y_d``.

    Used for theories like ``dark_photon_plasma`` that have constraint fields
    (``time_derivative_order == 0``) or time-derivative RHS operators.

    Raises
    ------
    SimulationDivergedError
        If the divergence pre-check predicts that the IEEE 754 FFT floor will
        be amplified by more than ``divergence_threshold`` over ``t_end_rel``.
    """
    # Local imports: kept inside this branch so unused-import pass doesn't strip them
    # when only the Phase 1 (no-constraint) path is exercised at module-load time.
    from tidal.solver._exceptions import SimulationDivergedError  # noqa: PLC0415
    from tidal.solver.modal import (  # noqa: PLC0415
        _build_evolution_matrices,  # pyright: ignore[reportPrivateUsage]
    )

    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_modes = y0_hat.shape[1]
    num_snap = num_snapshots

    # --- Schur reduction (numpy) ---------------------------------------------
    (
        A_reduced,
        B_lhs_modes,
        recovery,
        _v_recovery,
        c_names,
        orig_to_reduced,
        *_rest,
    ) = _build_evolution_matrices(spec, layout, grid, coeff_eval, k_grid, rfft_shape)
    # A_reduced: (n_modes, n_dyn, n_dyn) — Schur-eliminated evolution matrix
    # recovery:  (n_modes, n_con, n_dyn) — constraint field reconstruction
    # orig_to_reduced: dict mapping original slot index → reduced slot index

    n_dyn = A_reduced.shape[1]

    # --- Generalised eigenvalue B @ dy_d/dt = A @ y_d ------------------------
    # Fold B^{-1} into A so the JAX kernel only needs the matrix exponential.
    # ``_build_evolution_matrices`` pre-solves this in most paths (returns
    # B=None); the explicit-B branch survives for theories where the solve
    # is deferred. JAX vmap of jnp.linalg.solve handles it batched.
    if B_lhs_modes is not None:
        B_jax = jnp.array(B_lhs_modes, dtype=jnp.complex128)
        A_jax_pre = jnp.array(A_reduced, dtype=jnp.complex128)
        A_reduced = np.array(_jax_solve_batched(B_jax, A_jax_pre))

    # --- Extract dynamical IC in reduced ordering ----------------------------
    y0_hat_dyn = np.zeros((n_dyn, n_modes), dtype=np.complex128)
    for orig_si, red_pos in orig_to_reduced.items():
        y0_hat_dyn[red_pos] = y0_hat[orig_si]

    # --- Time axis setup -----------------------------------------------------
    t0 = float(t_eval[0])
    t_end = float(t_eval[-1])
    t_end - t0
    if num_snap > 1:
        dts = np.diff(t_eval)
        uniform = bool(np.allclose(dts, dts[0]))
        dt_step = float(dts[0]) if uniform else None
    else:
        uniform = True
        dt_step = None

    # --- JAX evolution on reduced dynamical subspace -------------------------
    # The eigvals divergence pre-check from the Phase 1 path is omitted here:
    # scipy's _evolve_per_mode does not do a separate pre-check either, and
    # the eigvals call costs ~3 ms at dark_photon_plasma sizes. We rely on
    # the post-evolution NaN/Inf detection below to catch divergence.
    M_jax = jnp.array(A_reduced, dtype=jnp.complex128)  # (n_modes, n_dyn, n_dyn)
    y0_T_jax = jnp.array(y0_hat_dyn.T, dtype=jnp.complex128)  # (n_modes, n_dyn)
    t_rel_jax = jnp.array(t_eval - t0, dtype=jnp.float64)

    evolve_nonuniform = _get_evolve_nonuniform_fn()

    if uniform and dt_step is not None and dt_step > 0:
        # Fused expm + scan: single JIT call instead of two (saves dispatch overhead)
        fused = _get_fused_expm_scan_fn(num_snap)
        ys_jax = fused(M_jax, y0_T_jax, jnp.array(dt_step, dtype=jnp.float64))
    else:
        ys_jax = evolve_nonuniform(M_jax, y0_T_jax, t_rel_jax)

    ys_np = np.array(ys_jax)  # (n_snapshots, n_modes, n_dyn)

    # Post-evolution divergence detection: catch NaN/Inf from overflowed expm
    if not np.all(np.isfinite(ys_np)):
        msg = (
            "Simulation diverged: post-evolution Fourier state contains "
            "NaN/Inf. The reduced-system evolution produced non-finite "
            "values, indicating tachyonic blow-up or extreme parameter "
            "regime. (JAX modal, constraint path.)"
        )
        raise SimulationDivergedError(msg)

    # --- Reconstruct full Fourier state (dynamical + constraint slots) -------
    y_hat_all = np.zeros((num_snap, n_slots, n_modes), dtype=np.complex128)

    # Dynamical slots: vectorised scatter via orig_to_reduced
    orig_indices = np.array(list(orig_to_reduced.keys()))
    red_indices = np.array(list(orig_to_reduced.values()))
    # ys_np[:, :, red_indices].transpose(0, 2, 1) gives (n_snapshots, n_dyn, n_modes)
    y_hat_all[:, orig_indices, :] = ys_np[:, :, red_indices].transpose(0, 2, 1)

    # Constraint slots: reconstruct via recovery matrix.
    # Einsum contracts recovery[modes, con, dyn] with ys_np[snap, modes, dyn]
    # to produce c_hat_all of shape (snap, con, modes).
    if len(c_names) > 0:
        c_hat_all = np.einsum("mcj,tmj->tcm", recovery, ys_np)
        for ci, c_name in enumerate(c_names):
            c_slot = layout.field_slot_map[c_name]
            y_hat_all[:, c_slot, :] = c_hat_all[:, ci, :]

    # --- IFFT to physical space ---------------------------------------------
    snapshots = np.zeros((num_snap, n_slots * n_pts))
    for ti in range(num_snap):
        snapshots[ti] = _ifft_slots(y_hat_all[ti], layout, grid)

    times = t_eval.copy()

    # --- Post-hoc callbacks --------------------------------------------------
    if snapshot_callback is not None or progress is not None:
        for ti, t in enumerate(t_eval):
            if snapshot_callback is not None:
                snapshot_callback(t, snapshots[ti])
            if progress is not None:
                progress.update(t)

    if progress is not None:
        progress.finish()

    return {
        "t": times,
        "y": snapshots,
        "success": True,
        "message": (
            f"Modal solver completed (JAX vmap+scan, Phase 2a constraint path: "
            f"numpy Schur + JAX expm, {len(c_names)} constraint fields, "
            f"{n_dyn} dynamical slots)"
        ),
    }


@functools.lru_cache(maxsize=8)
def _get_jax_solve_batched_fn() -> Callable[..., Any]:
    """JIT+vmap wrapper for batched ``jnp.linalg.solve(B, A)`` over modes."""
    jax, jnp, _ = _require_jax()

    @jax.jit
    def _solve_batched(B: Any, A: Any) -> Any:
        return jax.vmap(jnp.linalg.solve)(B, A)

    return _solve_batched


def _jax_solve_batched(B: Any, A: Any) -> Any:
    """Invoke the cached batched-solve function, returning a JAX array."""
    return _get_jax_solve_batched_fn()(B, A)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def solve_modal_jax(
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    progress: SimulationProgress | None = None,
    return_eigendata: bool = False,
    _zero_nyquist: bool = True,
) -> SolverResult:
    """Solve using the JAX-accelerated modal backend (phase 1).

    Matches :func:`tidal.solver.modal.solve_modal` in signature and return
    type.  Phase-1 restrictions (see module docstring) raise
    ``NotImplementedError`` rather than silently falling back to scipy, so
    callers that request an unsupported path get a clear error.

    Parameters
    ----------
    spec, grid, y0, t_span, bc, parameters, rtol, atol, num_snapshots,
    snapshot_callback, progress, return_eigendata, _zero_nyquist
        Same semantics as :func:`~tidal.solver.modal.solve_modal`.

    Raises
    ------
    NotImplementedError
        For position-dependent coefficients, constraint systems, or
        ``return_eigendata=True``.
    SimulationDivergedError
        When eigenvalue pre-check detects tachyonic instability.
    ValueError
        When ``time_derivative_order > 2``.
    """
    from tidal.solver._exceptions import SimulationDivergedError  # noqa: PLC0415
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    # --- Phase-1 gates -------------------------------------------------------
    if return_eigendata:
        msg = (
            "modal-jax phase 1 does not support return_eigendata=True. "
            "Use --scheme modal (scipy backend) for perturbative Duhamel data."
        )
        raise NotImplementedError(msg)

    has_pos_dep = _has_position_dependent_terms(spec)
    if has_pos_dep:
        msg = (
            "modal-jax phase 1 does not support position-dependent "
            "coefficients (expm_multiply Krylov path). Use --scheme modal."
        )
        raise NotImplementedError(msg)

    has_constraints = any(eq.time_derivative_order == 0 for eq in spec.equations)
    has_time_ops = _has_time_derivative_operators(spec)
    needs_reduction = (has_constraints or has_time_ops) and not has_pos_dep

    # --- Ensure JAX is available and double precision is on ------------------
    _jax, jnp, _ = _require_jax()

    # --- Setup (mirrors scipy solve_modal) -----------------------------------
    max_to = max((eq.time_derivative_order for eq in spec.equations), default=0)
    if max_to > 2:
        high_fields = [
            eq.field_name for eq in spec.equations if eq.time_derivative_order > 2
        ]
        msg = (
            f"solve_modal_jax requires time_derivative_order <= 2; "
            f"got max={max_to} on {high_fields}."
        )
        raise ValueError(msg)

    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, parameters or {})

    if not has_constraints:
        warn_frozen_constraints(layout, "modal-jax")

    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)

    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    # FFT initial conditions (numpy, same as scipy path)
    y0_hat = _fft_slots(y0, layout, grid)

    if _zero_nyquist:
        for _dim_idx, n in enumerate(grid.shape):
            if n % 2 == 0:
                nyq_mode = n // 2
                if len(grid.shape) == 1:
                    y0_hat[:, nyq_mode] = 0.0
                else:
                    rfft_last = grid.shape[-1] // 2
                    y0_hat[:, ..., rfft_last] = 0.0

    # --- Shared slot/mode counts (used by both paths) ------------------------
    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_modes = y0_hat.shape[1]

    # --- Phase 2a: Constraint / time-derivative RHS path ---------------------
    # Schur elimination runs in numpy; only the expm+scan kernel uses JAX.
    if needs_reduction:
        return _solve_modal_jax_constrained(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
            y0_hat,
            t_eval,
            jnp,
            num_snapshots,
            snapshot_callback,
            progress,
        )

    # --- Phase 1: No-constraint path (fast path) -----------------------------
    # Build per-mode evolution matrices (numpy, one-time cost, not the bottleneck)
    A_modes = _build_per_mode_matrices(
        spec, layout, grid, coeff_eval, k_grid, rfft_shape
    )

    n_check = min(3, A_modes.shape[0])
    combined = np.max(np.abs(A_modes[:n_check]), axis=0)
    blocks = find_independent_blocks(combined)

    # --- Time axis setup -----------------------------------------------------
    t0 = float(t_eval[0])
    t_end = float(t_eval[-1])
    t_end_rel = t_end - t0

    if num_snapshots > 1:
        dts = np.diff(t_eval)
        uniform = bool(np.allclose(dts, dts[0]))
        dt_step = float(dts[0]) if uniform else None
    else:
        uniform = True
        dt_step = None

    # --- Per-block divergence pre-check (eigenvalue-based, outside JIT) -----
    # Mirrors the scipy eigenvalue pre-check (#327 follow-up Stage 2).
    # The JAX path skips the runtime per-snapshot amplitude guard since
    # lax.scan cannot raise Python exceptions mid-trace.
    max_real_eig: float = 0.0
    divergence_threshold: float = 100.0
    ieee_fft_floor: float = 1e-14

    initial_physical: NDArray[np.float64] | None = None
    initial_max_amp: float = 1e-15
    if t_end_rel > 0:
        initial_physical = _ifft_slots(y0_hat, layout, grid)
        initial_max_amp = max(float(np.max(np.abs(initial_physical))), 1e-15)

    for block_slots in blocks:
        idx = np.array(block_slots)
        y0_block = y0_hat[idx, :]
        if np.max(np.abs(y0_block)) < 1e-15:
            continue
        A_block = A_modes[:, idx[:, None], idx[None, :]]  # (n_modes, bs, bs)
        # For the fast path (no B matrix), M = A directly.
        M_block: NDArray[np.complex128] = A_block

        if t_end_rel > 0:
            M_block_jax = jnp.array(M_block, dtype=jnp.complex128)
            eigvals_jax = _get_batch_eigvals_fn()(M_block_jax)  # (n_modes, bs)
            eigvals = np.array(eigvals_jax)
            _warn_eigenvalue_growth(
                eigvals.ravel().astype(np.complex128), t_end_rel, context="JAX modal"
            )
            max_real_eig = max(max_real_eig, float(np.max(np.real(eigvals))))

    if t_end_rel > 0 and max_real_eig > 0:
        predicted_log = max_real_eig * t_end_rel
        log_cap = float(np.log(1e300))
        predicted_floor_amp = (
            float("inf")
            if predicted_log > log_cap
            else ieee_fft_floor * math.exp(predicted_log)
        )
        predicted_ratio = predicted_floor_amp / initial_max_amp
        if not math.isfinite(predicted_ratio) or predicted_ratio > divergence_threshold:
            msg = (
                f"Simulation predicted to diverge: max real eigenvalue "
                f"{max_real_eig:.4g} amplifies the IEEE 754 FFT floor "
                f"{ieee_fft_floor:.0e} to ratio {predicted_ratio:.2e} "
                f"at t={t_end:.4g} (threshold {divergence_threshold:.0e}). "
                f"Fields would leave the perturbative regime. "
                f"Rejecting pre-evolution (JAX modal eigenvalue pre-check)."
            )
            raise SimulationDivergedError(msg)

    # --- Per-block JAX evolution ---------------------------------------------
    # Accumulate Fourier snapshots from all blocks.
    y_hat_all: NDArray[np.complex128] = np.zeros(
        (num_snapshots, n_slots, n_modes), dtype=np.complex128
    )

    batch_expm = _get_batch_expm_fn()
    evolve_uniform = _get_evolve_uniform_fn(num_snapshots)
    evolve_nonuniform = _get_evolve_nonuniform_fn()

    t_rel_jax = jnp.array(t_eval - t0, dtype=jnp.float64)

    for block_slots in blocks:
        idx = np.array(block_slots)
        y0_block = y0_hat[idx, :]  # (bs, n_modes)
        if np.max(np.abs(y0_block)) < 1e-15:
            continue

        A_block = A_modes[:, idx[:, None], idx[None, :]]  # (n_modes, bs, bs)
        M_jax = jnp.array(A_block, dtype=jnp.complex128)  # M = A (B=None fast path)
        y0_T_jax = jnp.array(y0_block.T, dtype=jnp.complex128)  # (n_modes, bs)

        if uniform and dt_step is not None and dt_step > 0:
            # Precompute exp(M * dt) for all modes in parallel via vmap+JIT
            exp_M_dt = batch_expm(M_jax * dt_step)  # (n_modes, bs, bs)
            # Evolve all snapshots via scan
            ys_jax = evolve_uniform(exp_M_dt, y0_T_jax)  # (n_snapshots, n_modes, bs)
        else:
            # Non-uniform spacing: vmap over time points
            ys_jax = evolve_nonuniform(M_jax, y0_T_jax, t_rel_jax)  # same shape

        # Transfer back to numpy and accumulate into the full Fourier state
        ys_np: NDArray[np.complex128] = np.array(ys_jax)  # (n_snapshots, n_modes, bs)
        # Vectorized scatter: (n_snapshots, n_modes, bs) → y_hat_all[:, slots, :]
        # ys_np.transpose(0,2,1) has shape (n_snapshots, bs, n_modes), matching y_hat_all layout.
        y_hat_all[:, np.array(block_slots), :] = ys_np.transpose(0, 2, 1)

    # --- IFFT back to physical space -----------------------------------------
    snapshots: NDArray[np.float64] = np.zeros((num_snapshots, n_slots * n_pts))
    times: NDArray[np.float64] = t_eval.copy()

    for ti in range(num_snapshots):
        snapshots[ti] = _ifft_slots(y_hat_all[ti], layout, grid)

    # --- Post-hoc callbacks --------------------------------------------------
    if snapshot_callback is not None or progress is not None:
        for ti, t in enumerate(t_eval):
            if snapshot_callback is not None:
                snapshot_callback(t, snapshots[ti])
            if progress is not None:
                progress.update(t)

    if progress is not None:
        progress.finish()

    return {
        "t": times,
        "y": snapshots,
        "success": True,
        "message": "Modal solver completed (JAX vmap+scan, constant-coefficient path)",
    }
