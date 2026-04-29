"""Pre-simulation tachyonic onset detection for conversion channels.

Canonical probe: a **unit IC at the source field's slot** in the
constraint-eliminated reduced system, evolved at **every** Fourier
mode through the same Padé matrix exponential the modal solver's
Pass 0 uses (``scipy.linalg.expm``, Higham 2009).  Catches tachyonic
instabilities BEFORE any simulation runs, preventing artifacts from
being misidentified as amplification (#238).

Usage::

    from tidal.measurement._stability import check_conversion_stability

    result = check_conversion_stability(
        spec, grid, params,
        source="h_5", target="a_1",
        # t_test default is 20.0; can be overridden to match a sim t_end.
    )
    if not result.stable:
        print(f"TACHYONIC: effective growth rate {result.max_excess:.4f}/s")

Architecture (post-#323 Stage C refactor, 2026-04-29)
=====================================================

* Build ``M = B⁻¹·A`` per Fourier mode using the same
  ``_build_m_with_null_projection`` helper as Pass 0 evolution.
* Construct a **unit IC at the source field's slot** in the reduced
  per-block system; evolve at **every** k mode (no IC-coupling skip).
* Compute ``y(t) = expm(M·t_test) @ y₀`` via Padé scaling-and-squaring.
* Effective growth rate ``γ_eff = log(‖y(t)‖ / ‖y₀‖) / t_test``;
  reject if ``γ_eff > threshold`` at any k.

Why this design
---------------

The Stage C truth-table investigation (``stage_c_truth_table.csv``,
joining ``hpc_results/28520217/samples.profile_replay.csv`` with the
audit) showed three probe variants on the D1 amplify chain's 57
audited contamination samples:

* ``v1`` (unit IC at source slot, all k, ``threshold=0.3``,
  ``t_test=10``): **56/57** caught.
* ``v2`` (consistent-IC FFT'd into k-space, only IC-coupled k bins
  checked): **0/57**.
* ``v3`` candidate (consistent-ic-solve via
  ``ensure_consistent_ic`` + per-block FFT): **0/57**.

The shared blind spot of v2 / v3 is that the linearised IC has zero
overlap with non-fundamental k bins, so they never check those bins;
but the actual numerical simulation reaches those bins anyway via
off-grid plane-wave snap residual, snapshot rounding, and finite-
precision arithmetic.  v1's all-k unit-IC scan is therefore the
empirically correct probe.  The modal solver's Schur complement
already eliminates algebraic constraints from the per-block ``M``
matrix, so the probe operates on the same reduced system the solver
evolves; an additional ``ensure_consistent_ic`` step is redundant
*and* catches nothing not already caught by the all-k scan.

The default ``t_test`` is **20.0** (was 10.0 in v1).  The increase
catches additional slow-growth Hwang–Noh cases inline at probe time.
Probe cost is dominated by the spectral-radius prefilter (committed
6a7d638) which skips ``expm`` whenever ``‖M·t_test‖ ≤ threshold·t_test``;
the longer ``t_test`` adds one extra Padé squaring on the modes that
do trigger the call (≈5 % overall cost).  Sample 391's Hwang–Noh
slow-growth case has γ_eff(t=20)=0.272 — still below the 0.3
threshold; this last residual is caught by the inline ``P_max > 0.5``
gate in ``tidal/inference/_likelihood.py`` after the simulation runs.

References
----------
- Higham, N.J. (2009). "The scaling and squaring method for the matrix
  exponential revisited." SIAM Review 51(4):747-764.
- Higham, N.J. (2008). Functions of Matrices, SIAM, §2.3 (cond(V)
  abandonment threshold for diagonalization).
- Issue #322 (Padé refactor; replaced eigenvalue+pinv approach).
- Issue #323 (Stage C investigation; canonical probe selection).
- ``docs/tex/stability_probe.tex`` (architecture + Stage C evidence).
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


# Stable provenance string written to chain CSVs (``stability_profile``
# column) so the chain-time probe used by a given chain can be
# identified after the fact.  Bumping the probe is a public-API change:
# pick a new name and update ``docs/tex/stability_probe.tex``.
PROBE_PROFILE_NAME: str = "unit-ic-all-k-0.3"

# Lower bound of the borderline strip used by the post-hoc audit's
# sample selection.  Empirically the false-negative regime is the top
# ~30 % of the ``[0, threshold]`` interval (#323).
_BORDERLINE_LOW_FRACTION: float = 0.7


def _borderline_low(threshold: float) -> float:
    """Lower bound of the borderline strip for the given threshold."""
    return _BORDERLINE_LOW_FRACTION * threshold


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

    borderline: bool = False
    """True if the worst gamma_eff falls in the borderline strip
    (``[0.7·threshold, threshold]``) for any IC-coupled k mode.

    The borderline strip is the empirical false-negative regime
    documented in #323: probe verdict ``stable=True`` but the simulation
    can still evolve into the non-perturbative regime within ``t_end``.
    Post-hoc audit prioritises borderline samples for re-simulation."""

    profile_name: str = PROBE_PROFILE_NAME
    """Stable probe-version identifier persisted into chain CSV
    metadata (``stability_profile`` column).  Always
    :data:`PROBE_PROFILE_NAME` for results produced by this module."""


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
    t_test: float = 20.0,
    n_extra_k: int = 4,  # noqa: ARG001
    conservative: bool = False,  # noqa: ARG001
) -> ConversionStabilityResult:
    """Check whether the source-IC excites a growing mode in any Fourier block.

    Mirrors the modal solver's Pass 0 path-D evolution
    (``_evolve_per_mode_pade``, ``tidal/solver/modal.py:1842``).  For
    each Fourier mode of the source-containing constraint-eliminated
    block:

    1. Build ``M = B⁻¹·A`` via ``_build_m_with_null_projection`` — same
       construction as the solver evolves with.
    2. Construct a **unit IC** at the source field's slot in the
       reduced per-block system.
    3. Compute ``y(t_test) = expm(M·t_test) @ y₀`` via Padé scaling-
       and-squaring (``scipy.linalg.expm``, Higham 2009).
    4. ``γ_eff = log(‖y(t_test)‖ / ‖y₀‖) / t_test``.  If
       ``γ_eff > threshold`` at any k mode, classify the point
       tachyonic.

    Architecture justification (#323 Stage C)
    -----------------------------------------

    The unit-IC all-k design is empirically the only probe that
    catches the contamination samples on the D1 amplify chain (see
    module docstring + ``docs/tex/stability_probe.tex``).  The
    consistent-IC variants (v2 FFT, v3 ``ensure_consistent_ic``) place
    zero IC weight on non-fundamental k bins and miss tachyons that
    those bins carry — bins the actual numerical simulation reaches
    via off-grid plane-wave snap residual, snapshot rounding, and
    finite-precision arithmetic.  The modal solver's Schur complement
    already eliminates algebraic constraints from the reduced ``M``
    matrix, so the probe operates on the same reduced system the
    solver evolves; no separate ``ensure_consistent_ic`` step is
    needed.

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
        ``target`` is currently unused (kept for API stability).
    baseline_overrides : dict, optional
        Unused; kept for API compatibility.
    ic_wavevector : float, optional
        IC wavenumber (used only for diagnostic messages; all k modes
        are checked regardless).  Default: ``2π/L`` (fundamental mode).
    threshold : float
        Maximum effective growth rate ``γ_eff`` before the run is
        classified as tachyonic (default: 0.3).  ``exp(0.3·t_end) ≈
        20×`` source-norm growth at ``t_end=10`` — typically the
        boundary where the simulation's own divergence pre-check
        (#298) catches things, matching the simulation's own
        conservatism.
    t_test : float
        Probe time for the Padé evolution (default: 20.0; was 10.0 in
        the v1 probe).  Larger ``t_test`` catches additional slow-
        growth Hwang–Noh cases by averaging the gain over a longer
        evolution window.  Probe cost increases ≈5 % (one extra Padé
        squaring on the modes that survive the spectral-radius
        prefilter).  Sample 391's slow-growth case has
        ``γ_eff(t=20)=0.272 < 0.3``; that final residual is caught
        by the inline ``P_max > 0.5`` gate in the likelihood after
        the simulation runs.
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
        ``stable=True`` iff no k-mode produces ``γ_eff > threshold``
        when the unit-IC is evolved through ``expm(M·t_test)``.
        ``max_excess`` is the worst observed ``γ_eff``.  ``borderline``
        flags points whose worst ``γ_eff`` lies in
        ``[0.7·threshold, threshold]`` — the empirical false-negative
        strip; surface these to the post-hoc audit.
    """
    from tidal.measurement._stability_cache import get_structural_context
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.modal import (
        _build_evolution_matrices,  # type: ignore[reportPrivateUsage]
        _build_m_with_null_projection,  # type: ignore[reportPrivateUsage]
        find_independent_blocks,
    )

    # Cached parameter-independent setup: layout, k_vals, rfft_shape,
    # ic_k_default, source-slot lookup.  See _stability_cache.py for
    # the cache-key safety analysis.  This shaves the StateLayout +
    # FFT-grid + source-lookup overhead off every call.
    static = get_structural_context(spec, grid, source)
    layout = static.layout
    k_vals = static.k_vals
    rfft_shape = static.rfft_shape
    k_grid: list[np.ndarray[tuple[int], np.dtype[np.float64]]] = [k_vals]
    int(grid.shape[0])

    # Source field absent from the layout → trivially stable; the
    # simulation cannot exercise this conversion channel.
    if static.early_return_reason is not None:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=static.early_return_reason,
            profile_name=PROBE_PROFILE_NAME,
        )

    # IC wavenumber is used only for the diagnostic message.  All k
    # modes are checked regardless because the actual numerical
    # simulation reaches non-fundamental bins via off-grid plane-wave
    # snap residual, snapshot rounding, and finite-precision arithmetic
    # (#323 Stage C investigation).
    ic_k = ic_wavevector if ic_wavevector is not None else static.ic_k_default

    # All k modes are checked.  ``y₀`` is a unit vector at the source
    # slot, identical at every k — see the architecture-justification
    # block in the docstring for why this is the empirically correct
    # design.  Cost is dominated by the spectral-radius prefilter below.
    k_indices = list(range(len(k_vals)))

    # Build constraint-eliminated system.  ``mapping`` is
    # ``orig_to_reduced``: it translates full-layout slot indices to
    # the reduced-system indices that A_test/B_test/blocks all use.
    # Some parameter points produce singular constraint matrices
    # (np.linalg.inv on S_cc_reg fails); those are physically degenerate
    # and the simulation cannot run there either, so we mark them
    # tachyonic with infinite excess.
    ce = CoefficientEvaluator(spec, grid, parameters)
    try:
        A_test, B_test, _, _, _, mapping, _, _ = _build_evolution_matrices(
            spec,
            layout,
            grid,
            ce,
            k_grid,
            rfft_shape,
        )
    except (np.linalg.LinAlgError, ValueError) as exc:
        return ConversionStabilityResult(
            stable=False,
            max_excess=float("inf"),
            k_tachyonic=None,
            n_tachyonic_modes=1,
            message=(
                f"Constraint elimination failed at this parameter point "
                f"({type(exc).__name__}: {exc}); simulation cannot run, "
                f"treating as tachyonic."
            ),
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
    # The "field not found" path is handled upstream via the cached
    # ProbeStructuralContext; ``static.src_slot_full`` is guaranteed
    # non-None at this point.
    src_slot_full = static.src_slot_full
    assert src_slot_full is not None

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

    # Build the per-mode evolution generator M = B⁻¹·A for the source
    # block, using the same null-space-aware construction as the modal
    # solver (handles rank-deficient B via SVD projection onto range(B);
    # null directions get M=0 so they stay at IC, matching solver
    # semantics).
    A_block = A_test[:, idx[:, None], idx[None, :]].astype(np.complex128)
    if B_test is not None:
        B_block = B_test[:, idx[:, None], idx[None, :]].astype(np.complex128)
    else:
        B_block = None
    M_block = _build_m_with_null_projection(A_block, B_block)

    # Unit IC at the source slot, replicated at every k mode.  Same
    # vector at every k by construction — no need to allocate a full
    # ``(block_size, n_modes)`` matrix; the IC norm is exactly 1 at
    # every k.
    y0 = np.zeros(block_size, dtype=np.complex128)
    y0[src_slot_in_block] = 1.0

    import scipy.linalg as sla  # type: ignore[import-untyped]

    max_excess = 0.0
    worst_k: float | None = None
    worst_excess = 0.0
    n_tachyonic = 0
    n_modes_checked = 0

    # Spectral-radius cutoff (Henrici 1962, Gershgorin bound):
    # ``norm(expm(M·t), 2) ≤ exp(σ_max(M·t)) ≤ exp(norm(M·t, 2))``.
    # If ``norm(M·t, 2) ≤ threshold·t_test``, then ``γ_eff ≤
    # threshold`` and the mode cannot exceed the threshold.  Skipping
    # the ``expm`` call in that regime saves the bulk of the per-call
    # cost on stable parameter points (committed 6a7d638).
    cutoff_norm_M_t = threshold * t_test

    for ki in k_indices:
        n_modes_checked += 1
        M_k = M_block[ki]
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
            # If ``expm`` itself fails (rare; only on pathological
            # inputs), treat as unstable to be safe (committed 0581fdb).
            n_tachyonic += 1
            max_excess = max(max_excess, float("inf"))
            worst_k = float(k_vals[ki])
            continue
        y_t = exp_M @ y0
        y_t_norm = float(np.linalg.norm(y_t))
        if not math.isfinite(y_t_norm) or y_t_norm <= 0.0:
            n_tachyonic += 1
            if not math.isfinite(y_t_norm):
                max_excess = float("inf")
            worst_k = float(k_vals[ki])
            continue
        # ``‖y₀‖ = 1`` by construction.
        gamma_eff = math.log(y_t_norm) / t_test
        worst_excess = max(worst_excess, gamma_eff)
        if gamma_eff > threshold:
            n_tachyonic += 1
            if gamma_eff > max_excess:
                max_excess = gamma_eff
                worst_k = float(k_vals[ki])

    borderline = _borderline_low(threshold) <= worst_excess <= threshold

    if n_tachyonic > 0:
        msg = (
            f"Source-IC excites a growing mode: {n_tachyonic}/{n_modes_checked} "
            f"k-modes have effective growth rate gamma_eff > {threshold} "
            f"(Padé probe at t_test={t_test}, profile={PROBE_PROFILE_NAME!r}, "
            f"ic_k={ic_k:.4f}). "
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
            borderline=False,  # rejected — borderline only meaningful for stable
            profile_name=PROBE_PROFILE_NAME,
        )

    return ConversionStabilityResult(
        stable=True,
        max_excess=worst_excess,
        k_tachyonic=None,
        n_tachyonic_modes=0,
        message=(
            f"Stable: source-IC effective growth rate gamma_eff ≤ {threshold} "
            f"across {n_modes_checked} k modes "
            f"(worst gamma_eff={worst_excess:.4f}/s; "
            f"profile={PROBE_PROFILE_NAME!r}, t_test={t_test}, "
            f"ic_k={ic_k:.4f}). "
            + ("BORDERLINE — flagged for post-hoc audit." if borderline else "")
        ),
        borderline=borderline,
        profile_name=PROBE_PROFILE_NAME,
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
