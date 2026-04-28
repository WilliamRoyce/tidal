"""Pre-simulation tachyonic onset detection for conversion channels.

Probes the constraint-eliminated source-containing block by evolving a unit
IC vector at the source slot through the same Padé matrix exponential that
the modal solver's Pass 0 uses (`scipy.linalg.expm`, Higham 2009 scaling-
and-squaring).  Catches tachyonic instabilities BEFORE any simulation runs,
preventing artifacts from being misidentified as amplification (#238).

Usage::

    from tidal.measurement._stability import check_conversion_stability

    result = check_conversion_stability(
        spec, grid, params,
        source="h_5", target="a_1",
        t_test=10.0,    # match the simulation t_end
    )
    if not result.stable:
        print(f"TACHYONIC: effective growth rate {result.max_excess:.4f}/s")

Architecture (post-#322 refactor, 2026-04-26):
==============================================

The earlier eigenvalue+``pinv(V)`` implementation became unreliable on
models with intrinsically high ``cond(V)`` (e.g. T4 Ricci-EM had cond(V)
~ 1e14--1e17, making the IC-coupling filter fire on every parameter
point — see issue #322).  The fix mirrors the modal solver's path-D
evolution:

* Build ``M = B⁻¹·A`` per Fourier mode using the same
  ``_build_m_with_null_projection`` helper as Pass 0 evolution.
* Construct a **unit IC at the source field's slot** (the actual h_5
  plane-wave IC reduces to this in the source-containing block).
* Compute ``y(t) = expm(M·t_test) @ y0`` via Padé scaling-and-squaring.
* Effective growth rate ``gamma_eff = log(‖y(t)‖ / ‖y₀‖) / t_test`` is the
  multiplicative-growth-aware analogue of ``max Re(λ)`` from the old
  approach, but **cond(V)-independent** because no eigendecomposition or
  pseudoinverse is involved.

Key correctness property: IC-decoupled growing modes correctly contribute
**zero growth** to ``y(t)`` because they don't appear in the source IC's
projection onto the eigenbasis.  The previous eigenvalue+``pinv(V)`` path
attempted to recover this via a coupling filter on ``pinv(V)[:, src_slot]``,
which collapsed at high ``cond(V)``.  The Padé probe sidesteps the
problem entirely by operating on the actual IC vector.

References
----------
- Higham, N.J. (2009). "The scaling and squaring method for the matrix
  exponential revisited." SIAM Review 51(4):747-764.
- Higham, N.J. (2008). Functions of Matrices, SIAM, §2.3 (cond(V)
  abandonment threshold for diagonalization).
- Issue #322 (root cause + resolution discussion).
- ``docs/tex/modal_solver.tex`` §"Robust Matrix-Exponential Evolution".
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
    threshold: float = 0.1,
    t_test: float = 10.0,
    n_extra_k: int = 4,  # noqa: ARG001
    conservative: bool = False,  # noqa: ARG001
) -> ConversionStabilityResult:
    """Check whether the source-IC excites a growing mode in any Fourier block.

    Mirrors the modal solver's Pass 0 path-D evolution (`_evolve_per_mode_pade`,
    `tidal/solver/modal.py:1842`).  For each Fourier mode of the source-
    containing constraint-eliminated block:

    1. Build ``M = B⁻¹·A`` via ``_build_m_with_null_projection`` — same
       construction as the solver evolves with.
    2. Construct a unit IC ``y₀`` at the source field's slot in the block
       (a plane-wave IC on `source` reduces to this in the per-block
       reduced system).
    3. Compute ``y(t_test) = expm(M·t_test) @ y₀`` via Padé scaling-and-
       squaring (`scipy.linalg.expm`, Higham 2009).
    4. Effective growth rate ``gamma_eff = log(‖y(t_test)‖ / ‖y₀‖) / t_test``.
       If ``gamma_eff > threshold`` for any k mode, classify the parameter
       point as tachyonic.

    Why this is correct, and why it replaces the previous eigenvalue-and-
    pseudoinverse approach
    --------------------------------------------------------------------

    The earlier implementation computed the spectrum ``Re(λ_i)`` and the
    eigenvector matrix ``V``, then estimated the IC's projection onto
    each eigenmode via ``pinv(V)[:, src_slot]``.  At ``cond(V) > 1e13``
    (Higham 2008 §2.3 diagonalization-abandonment threshold) ``pinv(V)``
    is numerically meaningless — the IC-coupling filter cannot
    distinguish "eigenmode genuinely uncoupled to source" from
    "eigenvector noise at floor of ``cond(V)·ε``".  T4 (Ricci-EM) hits
    ``cond(V) ~ 1e14--1e17`` across most parameter points, causing 100%
    rejection on inputs that ``tidal simulate`` evolves cleanly to
    ``t_end``.  See issue #322 root-cause discussion.

    The Padé probe sidesteps the entire problem by operating on the
    *actual* IC vector, exactly as the modal solver does.  IC-decoupled
    growing modes contribute zero to ``y(t_test)`` because they're not
    excited by ``y₀`` — the same reason ``tidal simulate`` doesn't
    diverge.  No eigendecomposition, no pseudoinverse, no
    ``cond(V)``-sensitive bookkeeping.

    Threshold semantics are unchanged: ``threshold`` is still the
    rate above which growth is judged non-perturbative.  The reported
    ``max_excess`` is the worst ``gamma_eff`` across k modes, so existing
    callers comparing it against rate cutoffs (e.g. 0.3/s) get the same
    interpretation.

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
        Unused; kept for API compatibility.
    ic_wavevector : float, optional
        IC wavenumber (used only for documentation; all k modes are
        checked regardless).  Default: ``2π/L`` (fundamental mode).
    threshold : float
        Maximum effective growth rate ``gamma_eff`` in the source-block
        before the run is classified as tachyonic (default: 0.1).
        Physical oscillatory modes excited by the IC give
        ``gamma_eff ≈ 0``; values above ~0.1 indicate exponential growth
        that approaches Hwang-Noh non-perturbative regime within the
        simulation window.

        **Threshold history** (#323): the original 0.3 was chosen as
        "exp(0.3·t_end=10)≈20× growth — well past linear regime", but
        the Stage D2.0 (Bahamonde-PGT) smoke test (commit 712336d era)
        found a parameter point with γ_eff=0.275 that the probe
        accepted but where ``tidal simulate`` reached P_max=1.06 by
        t=10 (well into non-perturbative regime).  Tightened to 0.1
        on 2026-04-27 to eliminate this false-negative regime.
        exp(0.1·10) = 2.7× growth allowed — closer to Hwang-Noh's
        P<<1 linear-regime requirement.
    t_test : float
        Probe time for the Padé evolution (default: 10.0).  Ideally
        match the simulation's ``t_end`` so the probe reflects what the
        simulation will see.  Has only second-order effect on
        ``gamma_eff = log(‖y(t)‖/‖y₀‖)/t``: in the linear-eigenvalue limit
        ``gamma_eff → max(Re(λ_i))`` (subject to IC coupling), independent
        of ``t_test``.
    n_extra_k : int
        Deprecated.  All k modes are checked.
    conservative : bool
        Deprecated.  No longer used — the Padé probe is robust at any
        ``cond(V)``, so the conservative-fallback path that used to
        skip the unreliable ``pinv(V)`` filter is no longer needed.
        Kept for backward-compatible API; ignored.

    Returns
    -------
    ConversionStabilityResult
        ``stable=True`` iff no k-mode produces ``gamma_eff > threshold``
        when the source-IC is evolved through ``expm(M·t_test)``.
        ``max_excess`` is the worst observed ``gamma_eff``.
    """
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.modal import (
        _build_evolution_matrices,  # type: ignore[reportPrivateUsage]
        _build_m_with_null_projection,  # type: ignore[reportPrivateUsage]
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

    # Determine IC wavenumber index (documentation only — all k checked)
    if ic_wavevector is None:
        L = grid.bounds[0][1] - grid.bounds[0][0]
        ic_k = 2 * math.pi / L
    else:
        ic_k = ic_wavevector
    int(np.argmin(np.abs(k_vals - ic_k)))

    # Check ALL k modes.  CDT/PGT tachyonic modes can peak at any wavenumber.
    k_indices = list(range(len(k_vals)))

    # Build constraint-eliminated system.  ``mapping`` is
    # ``orig_to_reduced``: it translates full-layout slot indices to
    # the reduced-system indices that A_test/B_test/blocks all use.
    ce = CoefficientEvaluator(spec, grid, parameters)
    A_test, B_test, _, _, _, mapping, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        ce,
        k_grid,
        rfft_shape,
    )

    # Find independent blocks using low-k coupling structure.  block_slots
    # returned here are REDUCED-system indices, not full-layout indices.
    combined = np.max(np.abs(A_test[:3]), axis=0)
    blocks = find_independent_blocks(combined, threshold=1e-14)

    # Identify the block containing the source field.  Translate the
    # source field's full-layout slot to its reduced-system index before
    # comparing against block_slots.  (Pre-refactor code happened to work
    # when no constraint fields preceded the source slot, but failed
    # silently otherwise — see issue #322 root-cause analysis.)
    try:
        src_slot_full = layout.field_slot_map[source]
    except KeyError:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=f"Source field '{source}' not found in reduced system.",
        )

    if src_slot_full not in mapping:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=(
                f"Source field '{source}' (full slot {src_slot_full}) is "
                f"a constraint field eliminated from the reduced system."
            ),
        )
    src_reduced = mapping[src_slot_full]

    src_block: list[int] | None = None
    for block_slots in blocks:
        if src_reduced in block_slots:
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
    src_slot_in_block = list(idx).index(src_reduced)
    block_size = len(src_block)

    # Build the per-mode evolution generator M = B⁻¹·A for the source block,
    # using the same null-space-aware construction as the modal solver
    # (handles rank-deficient B via SVD projection onto range(B); null
    # directions get M=0 so they stay at IC, matching solver semantics).
    A_block = A_test[:, idx[:, None], idx[None, :]].astype(np.complex128)
    if B_test is not None:
        B_block = B_test[:, idx[:, None], idx[None, :]].astype(np.complex128)
    else:
        B_block = None
    M_block = _build_m_with_null_projection(A_block, B_block)

    # Probe each k mode with a unit IC at the source slot.  In the per-block
    # reduced system, a plane-wave IC on `source` *is* a unit vector at
    # src_slot_in_block at the IC's k mode (and zero at other k); but each k
    # mode probes the local M at that k, so the unit-IC probe is the right
    # universal stress test for "would a source-IC at this k get amplified?".
    #
    # Performance: empirically ~13-22 ms per call on the campaign workloads
    # (T1 17 fields x 17 k modes; T4 10 fields x 17 k modes).  Per inference
    # run (~10000 likelihood evals on 76 MPI ranks) the guard overhead is
    # ~3 s -- negligible compared to the ~30 min total wall time.  The
    # cheap-skip below for spectral-norm-bounded growth lets stable points
    # return almost immediately; the expm cost is paid only when the
    # spectral radius bound can't rule out growth.
    import scipy.linalg as sla  # type: ignore[import-untyped]

    y0 = np.zeros(block_size, dtype=np.complex128)
    y0[src_slot_in_block] = 1.0
    y0_norm = float(np.linalg.norm(y0))  # = 1.0 by construction; explicit for clarity

    max_excess = 0.0
    worst_k: float | None = None
    n_tachyonic = 0

    # Spectral-radius cutoff (Henrici 1962, gershgorin bound):
    # norm(expm(M*t), 2) <= exp(sigma_max(M*t)) <= exp(norm(M*t, 2)).
    # If norm(M*t, 2) <= threshold*t_test, then gamma_eff <= threshold,
    # i.e. cannot exceed the threshold.  Skipping the expm call in that
    # regime saves the bulk of the per-call cost on stable parameter
    # points.
    cutoff_norm_M_t = threshold * t_test

    for ki in k_indices:
        M_k = M_block[ki]
        # Cheap operator-norm prefilter: skip Padé when M·t is small enough
        # that no growth above the threshold is geometrically possible.
        m_norm_t = float(np.linalg.norm(M_k, ord=2)) * t_test
        if m_norm_t <= cutoff_norm_M_t:
            continue
        # Padé scaling-and-squaring: robust for arbitrary cond(V).
        try:
            exp_M = cast(
                "NDArray[np.complex128]",
                sla.expm(M_k * t_test),  # type: ignore[no-untyped-call]
            )
        except (ValueError, np.linalg.LinAlgError):
            # If expm itself fails (rare; only on pathological inputs),
            # treat as unstable to be safe.
            n_tachyonic += 1
            max_excess = max(max_excess, float("inf"))
            worst_k = float(k_vals[ki])
            continue
        y_t = exp_M @ y0
        y_t_norm = float(np.linalg.norm(y_t))
        if not math.isfinite(y_t_norm) or y_t_norm <= 0.0:
            # Overflow or numerical pathology → declare unstable.
            n_tachyonic += 1
            if not math.isfinite(y_t_norm):
                max_excess = float("inf")
            worst_k = float(k_vals[ki])
            continue
        # Effective growth rate over t_test.  In the eigenvalue limit this
        # converges to max(Re(λ_i)) restricted to modes the IC excites.
        gamma_eff = math.log(y_t_norm / y0_norm) / t_test
        if gamma_eff > threshold:
            n_tachyonic += 1
            if gamma_eff > max_excess:
                max_excess = gamma_eff
                worst_k = float(k_vals[ki])

    if n_tachyonic > 0:
        msg = (
            f"Source-IC excites a growing mode: {n_tachyonic}/{len(k_indices)} "
            f"k-modes have effective growth rate gamma_eff > {threshold} "
            f"(Padé probe at t_test={t_test}). "
            f"Worst: gamma_eff={max_excess:.4f}/s at k={worst_k:.4f} "
            f"(block size {block_size} fields). "
            f"Source IC will grow non-perturbatively; any P_max is an artefact."
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
        message=(
            f"Stable: source-IC effective growth rate gamma_eff ≤ {threshold} "
            f"across all k modes (Padé probe at t_test={t_test})."
        ),
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
