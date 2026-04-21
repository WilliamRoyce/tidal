"""Pre-simulation tachyonic onset detection for conversion channels.

Builds the Schur complement of the constraint-eliminated system and checks
eigenvalues at the IC wavenumber for positive-real-part excess over the GR
baseline. This catches tachyonic instabilities BEFORE any simulation runs,
preventing artifacts from being misidentified as amplification (#238).

Usage::

    from tidal.measurement._stability import check_conversion_stability

    result = check_conversion_stability(
        spec, grid, params,
        source="h_5", target="a_1",
        baseline_overrides={"delta1": 0.0},
    )
    if not result.stable:
        print(f"TACHYONIC: growth rate {result.max_excess:.4f}")
"""

from __future__ import annotations

import dataclasses
import math
from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


@dataclasses.dataclass(frozen=True)
class ConversionStabilityResult:
    """Result of a pre-simulation tachyonic onset check."""

    stable: bool
    """True if no tachyonic modes detected in the conversion channel."""

    max_excess: float
    """Maximum excess growth rate over GR baseline (0 if stable)."""

    k_tachyonic: float | None
    """Wavenumber where the worst tachyonic mode was found (None if stable)."""

    n_tachyonic_modes: int
    """Number of k modes with tachyonic excess above threshold."""

    message: str
    """Human-readable diagnostic message."""


def check_conversion_stability(  # noqa: C901, PLR0912, PLR0913, PLR0914, PLR0915
    spec: EquationSystem,
    grid: GridInfo,
    parameters: dict[str, float],
    *,
    source: str = "h_5",
    target: str = "a_1",  # noqa: ARG001
    baseline_overrides: dict[str, float] | None = None,  # noqa: ARG001
    ic_wavevector: float | None = None,
    threshold: float = 0.3,
    n_extra_k: int = 4,  # noqa: ARG001
) -> ConversionStabilityResult:
    """Check for tachyonic modes in the block containing the source field.

    Builds the constraint-eliminated first-order system at k=0 and near
    the IC wavenumber, finds the independent block containing the source
    field, and checks its maximum real eigenvalue against an absolute
    threshold.  Uses the generalized eigenvalue problem (A, B) as in the
    modal solver, giving the same eigenvalues the solver would use.

    The check always includes k=0 (DC mode) because plane-wave ICs on a
    finite grid have a non-trivial DC component (sinc leakage when k₀ is
    not an exact grid frequency), which can excite k=0 tachyonic modes
    and cause divergence even for ghost-free parameter points.

    The full source-containing block (not just the 4x4 {source, target}
    sub-matrix) is checked because CDT/PGT torsion models have non-trace
    torsion components that share a block with the source field and carry
    the instability; the 4x4 sub-matrix is always stable regardless.

    Tachyonic modes with negligible IC coupling to the source field are
    excluded from the tachyonic count.  These modes cannot be excited by
    a source-field IC and are suppressed by the solver's
    ``_suppress_tachyonic_noise``.  The coupling is measured via the
    inverse right-eigenvector matrix V⁻¹; a relative coupling
    ``|V⁻¹[i, src_slot]| / max|V⁻¹[:, src_slot]| < 1e-10`` is treated
    as IC-decoupled (matching the solver's ``coupling_threshold=1e-12``,
    with a 100x safety margin).  Values between 1e-10 and 1 will cause
    genuine exponential growth that the solver will not suppress.

    Parameters
    ----------
    spec : EquationSystem
        The equation system (must be kinetic-normalized).
    grid : GridInfo
        Grid configuration.
    parameters : dict
        Full parameter set for the simulation.
    source, target : str
        Source and target field names for the conversion channel.
    baseline_overrides : dict, optional
        Unused; kept for API compatibility.  The check now uses an
        absolute threshold on Re(λ) rather than an excess over a baseline.
    ic_wavevector : float, optional
        IC wavenumber (used only for documentation; all k modes are
        checked regardless).  Default: 2π/L (fundamental mode).
    threshold : float
        Maximum Re(λ) in the source-containing block before the run is
        classified as tachyonic (default: 0.3).  Physical oscillatory
        modes have Re(λ) ≈ 0; values above ~0.3 indicate genuine
        exponential growth.
    n_extra_k : int
        Deprecated.  Previously controlled the k-neighborhood scan
        width; now all k modes are checked and this parameter is
        ignored.

    Returns
    -------
    ConversionStabilityResult
        Contains stable flag, growth rate, and diagnostic message.
    """
    import scipy.linalg as sla  # type: ignore[import-untyped]

    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.modal import (
        _build_evolution_matrices,  # type: ignore[reportPrivateUsage]
        find_independent_blocks,
    )
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid.num_points)

    # Build Fourier wavenumber grid
    N = grid.shape[0]
    dx = grid.dx[0]
    k_vals = np.asarray(2 * math.pi * np.fft.rfftfreq(N, dx), dtype=np.float64)
    k_grid: list[np.ndarray[tuple[int], np.dtype[np.float64]]] = [k_vals]
    rfft_shape = (N // 2 + 1,)

    # Determine IC wavenumber index
    if ic_wavevector is None:
        L = grid.bounds[0][1] - grid.bounds[0][0]
        ic_k = 2 * math.pi / L
    else:
        ic_k = ic_wavevector
    int(np.argmin(np.abs(k_vals - ic_k)))

    # Check ALL k modes.
    #
    # CDT/PGT tachyonic modes can peak at any wavenumber depending on the
    # parameter combination — restricting to a small-k band + IC neighborhood
    # leaves a gap (e.g. k ~ 0.57-1.70 for N=256, L=100) where instabilities
    # are missed, causing runs to slip through the guard and diverge at runtime.
    # With N=256 there are 129 modes; each eigenvalue problem is <=20x20 so the
    # full scan costs only microseconds per sweep point.
    k_indices = list(range(len(k_vals)))

    # Build constraint-eliminated system (with B for generalized eig)
    ce = CoefficientEvaluator(spec, grid, parameters)
    A_test, B_test, _, _, _, _, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        ce,
        k_grid,
        rfft_shape,
    )

    # Find independent blocks; use low-k modes for coupling detection
    combined = np.max(np.abs(A_test[:3]), axis=0)
    blocks = find_independent_blocks(combined, threshold=1e-14)

    # Identify the block containing the source field
    try:
        src_slot = layout.field_slot_map[source]
    except KeyError:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=f"Source field '{source}' not found in reduced system.",
        )

    src_block: list[int] | None = None
    for block_slots in blocks:
        if src_slot in block_slots:
            src_block = list(block_slots)
            break

    if src_block is None:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=f"Source field '{source}' not in any coupled block.",
        )

    idx = np.array(src_block)

    # Check the full source-containing block at each k.
    # For each k where max Re(λ) > threshold, also verify that at least one
    # tachyonic mode has significant IC coupling to the source field.  Modes
    # that are decoupled from the source IC are suppressed by the solver's
    # _suppress_tachyonic_noise and never cause divergence.
    max_excess = 0.0
    worst_k = None
    n_tachyonic = 0
    src_slot_in_block = list(idx).index(src_slot)

    for ki in k_indices:
        Ak = A_test[ki][idx[:, None], idx[None, :]]
        if B_test is not None:
            Bk = B_test[ki][idx[:, None], idx[None, :]]
            # Null-space projection for rank-deficient B (issue #264).
            # Must match the fix in _evolve_per_mode: project A and B
            # onto range(B) before QZ to avoid finite spurious eigenvalues
            # from kinetic-null DOF (e.g. non-trace CDT torsion components).
            _u_bk, s_bk, Vt_bk = cast(
                "tuple[NDArray[np.complex128], NDArray[np.float64], NDArray[np.complex128]]",
                np.linalg.svd(Bk),
            )
            null_thresh: float = float(s_bk[0]) * 1e-10 if s_bk[0] > 0 else 1e-14
            rank_bk = int(np.sum(s_bk > null_thresh))
            null_dim = len(s_bk) - rank_bk
            if null_dim > 0:
                Vphys: NDArray[np.complex128] = np.asarray(
                    Vt_bk[:rank_bk].T,
                    dtype=np.complex128,
                )
                Vnull: NDArray[np.complex128] = np.asarray(
                    Vt_bk[rank_bk:].T,
                    dtype=np.complex128,
                )
                eig_r = sla.eig(  # type: ignore[reportUnknownVariableType]
                    Vphys.T @ Ak @ Vphys,
                    Vphys.T @ Bk @ Vphys,
                )
                ev_red = np.asarray(eig_r[0], dtype=np.complex128)  # type: ignore[reportUnknownArgumentType]
                vr_red = np.asarray(eig_r[1], dtype=np.complex128)  # type: ignore[reportUnknownArgumentType]
                ev = np.concatenate([ev_red, np.zeros(null_dim, dtype=np.complex128)])
                vr_mat: NDArray[np.complex128] = np.hstack([Vphys @ vr_red, Vnull])
            else:
                ev, vr_mat = cast(
                    "tuple[NDArray[np.complexfloating], NDArray[np.complexfloating]]",
                    sla.eig(Ak, Bk),  # type: ignore[arg-type]
                )
            # Zero out gauge/infinite eigenvalues (unphysical DOF)
            gauge = ~np.isfinite(ev) | (np.abs(ev) > 1e12)  # noqa: PLR2004
            ev = ev.copy()
            ev[gauge] = 0.0
        else:
            ev, vr_mat = cast(
                "tuple[NDArray[np.complexfloating], NDArray[np.complexfloating]]",
                np.linalg.eig(Ak),  # type: ignore[assignment]
            )

        max_re = float(np.max(np.real(ev)))

        if max_re > threshold:
            # IC coupling check: of the tachyonic modes (Re(λ) > threshold),
            # does any one project significantly onto the source field?
            # V_inv[i, j] = how much mode i is excited by a unit IC in slot j.
            tach_mask = np.real(ev) > threshold
            V_inv = np.linalg.pinv(vr_mat)
            src_col = np.abs(V_inv[:, src_slot_in_block])
            max_col = float(np.max(src_col))
            if max_col > 1e-20:  # noqa: PLR2004
                rel_coupling = float(np.max(src_col[tach_mask])) / max_col
            else:
                rel_coupling = 0.0

            # Condition-aware coupling threshold.  When vr_mat is
            # ill-conditioned (e.g. CDT non-trace torsion DOF at large ξ),
            # pinv(vr_mat) entries are unreliable: numerical noise of order
            # 1/cond(vr_mat) contaminates the coupling column, making truly
            # IC-decoupled tachyonic modes appear coupled.  The solver's
            # _suppress_tachyonic_noise uses the actual IC vector and
            # correctly finds ~1e-16 for these modes; the guard must
            # account for the conditioning to avoid false rejections.
            # See issue #266.
            cond_vr_mat = float(np.linalg.cond(vr_mat))
            # Floor: entries of pinv(vr_mat) have noise ~cond(vr_mat)*eps, so
            # rel_coupling has noise ~cond(vr_mat)*eps / max_col.  With
            # max_col ~ O(1), noise floor ~ cond(vr_mat) * eps.  Use a
            # conservative 100x margin on machine epsilon (1e-16).
            coupling_floor = max(1e-10, min(1.0, cond_vr_mat * 1e-14))
            if rel_coupling < coupling_floor:
                continue

            n_tachyonic += 1
            if max_re > max_excess:
                max_excess = max_re
                worst_k = float(k_vals[ki])

    if n_tachyonic > 0:
        msg = (
            f"Tachyonic mode in source-block: {n_tachyonic}/{len(k_indices)} "
            f"k-modes have max Re(λ) > {threshold}. "
            f"Worst: max_Re={max_excess:.4f} at k={worst_k:.4f} "
            f"(block size {len(src_block)} fields). "
            f"System will diverge exponentially; any P_max is an artefact."
        )
        return ConversionStabilityResult(
            stable=False,
            max_excess=max_excess,
            k_tachyonic=worst_k,
            n_tachyonic_modes=n_tachyonic,
            message=msg,
        )

    return ConversionStabilityResult(
        stable=True,
        max_excess=max_excess,
        k_tachyonic=None,
        n_tachyonic_modes=0,
        message="Stable: max Re(λ) in source-block within tolerance.",
    )


@dataclasses.dataclass(frozen=True)
class FullStabilityResult:
    """Result of a full-system stability check across all field blocks."""

    all_stable: bool
    """True if ALL blocks are stable (no tachyonic modes anywhere)."""

    n_tachyonic_blocks: int
    """Number of independent blocks with tachyonic excess."""

    max_growth_rate: float
    """Maximum excess growth rate across all blocks (0 if all stable)."""

    block_results: tuple[tuple[tuple[str, ...], bool, float], ...]
    """Per-block results: ((field_names,...), stable, excess) for each block."""

    message: str
    """Human-readable diagnostic."""


def check_full_stability(  # noqa: PLR0913, PLR0914
    spec: EquationSystem,
    grid: GridInfo,
    parameters: dict[str, float],
    *,
    baseline_overrides: dict[str, float] | None = None,
    threshold: float = 0.01,
    n_k_check: int = 5,
) -> FullStabilityResult:
    """Check ALL independent field blocks for tachyonic modes.

    Unlike ``check_conversion_stability`` which only examines the source→target
    block, this function examines every independent block in the reduced system
    to provide a complete stability map of the parameter point.

    Parameters
    ----------
    spec : EquationSystem
        The equation system (must be kinetic-normalized).
    grid : GridInfo
        Grid configuration.
    parameters : dict
        Full parameter set.
    baseline_overrides : dict, optional
        Parameters to override for GR baseline. Default: sets all non-κ/B₀
        coupling params to zero (detected automatically).
    threshold : float
        Minimum excess growth rate to flag as tachyonic (default: 0.01).
    n_k_check : int
        Number of k modes to check per block (default: 5, evenly spaced).

    Returns
    -------
    FullStabilityResult
    """
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.modal import (
        _build_evolution_matrices,  # type: ignore[reportPrivateUsage]
        find_independent_blocks,
    )
    from tidal.solver.state import StateLayout

    if baseline_overrides is None:
        # Default: zero out all params except kappa and B0
        baseline_overrides = {
            k: 0.0
            for k in parameters
            if k not in {"kappa", "B0", "grid_shape", "t_end"}
        }

    layout = StateLayout.from_spec(spec, grid.num_points)
    N = grid.shape[0]
    k_vals = np.asarray(2 * math.pi * np.fft.rfftfreq(N, grid.dx[0]), dtype=np.float64)
    k_grid: list[np.ndarray[tuple[int], np.dtype[np.float64]]] = [k_vals]
    rfft_shape = (N // 2 + 1,)
    n_modes = len(k_vals)

    # Build systems
    ce = CoefficientEvaluator(spec, grid, parameters)
    A_test, _, _, _, _, mapping, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        ce,
        k_grid,
        rfft_shape,
    )

    baseline_params = {**parameters, **baseline_overrides}
    ce_bl = CoefficientEvaluator(spec, grid, baseline_params)
    A_bl, _, _, _, _, _, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        ce_bl,
        k_grid,
        rfft_shape,
    )

    # Find independent blocks
    combined = np.max(np.abs(A_test[:3]), axis=0)
    blocks = find_independent_blocks(combined, threshold=1e-14)

    inv_map = {v: k for k, v in mapping.items()}

    # Select k indices to check (evenly spaced)
    k_indices = [
        int(i * (n_modes - 1) / max(n_k_check - 1, 1)) for i in range(n_k_check)
    ]
    k_indices = sorted(set(k_indices))

    block_results: list[tuple[tuple[str, ...], bool, float]] = []
    n_tachyonic_blocks = 0
    max_growth = 0.0

    for block_slots in blocks:
        idx = np.array(block_slots)
        field_names = tuple(
            layout.slots[inv_map[s]].name for s in block_slots if s in inv_map
        )

        # Use RELATIVE excess over baseline for ALL blocks.
        # Even torsion blocks have real eigenvalue pairs ±λ (wave
        # propagation characteristics) in the first-order system,
        # analogous to ±k for graviton/photon. The excess measures
        # how much the test parameters INCREASE the max growth rate
        # compared to the decoupled baseline (all coupling = 0).
        block_excess = 0.0
        for ki in k_indices:
            block_test = A_test[ki][idx[:, None], idx[None, :]]
            block_bl = A_bl[ki][idx[:, None], idx[None, :]]

            eigs_test = cast(
                "NDArray[np.complexfloating]",
                np.linalg.eigvals(block_test),
            )
            eigs_bl = cast("NDArray[np.complexfloating]", np.linalg.eigvals(block_bl))

            max_re_test = float(np.max(np.real(eigs_test)))
            max_re_bl = float(np.max(np.real(eigs_bl)))
            excess = max_re_test - max_re_bl
            block_excess = max(block_excess, excess)

        block_stable = block_excess <= threshold
        if not block_stable:
            n_tachyonic_blocks += 1
        max_growth = max(max_growth, block_excess)
        block_results.append((field_names, block_stable, block_excess))

    all_stable = n_tachyonic_blocks == 0
    if all_stable:
        msg = f"All {len(blocks)} blocks stable (max excess {max_growth:.4f})."
    else:
        unstable_names = [
            ",".join(names) for names, stable, _ in block_results if not stable
        ]
        msg = (
            f"{n_tachyonic_blocks}/{len(blocks)} blocks tachyonic "
            f"(max gamma={max_growth:.4f}). "
            f"Unstable: {'; '.join(unstable_names)}"
        )

    return FullStabilityResult(
        all_stable=all_stable,
        n_tachyonic_blocks=n_tachyonic_blocks,
        max_growth_rate=max_growth,
        block_results=tuple(block_results),
        message=msg,
    )
