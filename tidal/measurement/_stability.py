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
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
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


def check_conversion_stability(
    spec: EquationSystem,
    grid: GridInfo,
    parameters: dict[str, float],
    *,
    source: str = "h_5",
    target: str = "a_1",
    baseline_overrides: dict[str, float] | None = None,
    ic_wavevector: float | None = None,
    threshold: float = 0.01,
    n_extra_k: int = 4,
) -> ConversionStabilityResult:
    """Check for tachyonic modes in a source→target conversion channel.

    Builds the Schur complement at the IC wavenumber and a few surrounding
    k values, extracts the source+target field block eigenvalues, and
    compares with the GR baseline to detect tachyonic onset.

    Parameters
    ----------
    spec : EquationSystem
        The equation system (must be kinetic-normalised).
    grid : GridInfo
        Grid configuration.
    parameters : dict
        Full parameter set for the simulation.
    source, target : str
        Source and target field names for the conversion channel.
    baseline_overrides : dict, optional
        Parameters to override for the GR baseline comparison.
        Default: {"delta1": 0.0} (standard for torsion models).
    ic_wavevector : float, optional
        IC wavenumber to check. Default: 2π/L (fundamental mode).
    threshold : float
        Minimum excess growth rate to flag as tachyonic (default: 0.01).
    n_extra_k : int
        Number of extra k values to check around the IC mode (default: 4).

    Returns
    -------
    ConversionStabilityResult
        Contains stable flag, growth rate, and diagnostic message.
    """
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.modal import _build_constraint_eliminated_matrices
    from tidal.solver.state import StateLayout

    if baseline_overrides is None:
        baseline_overrides = {"delta1": 0.0}

    layout = StateLayout.from_spec(spec, grid.num_points)

    # Build Fourier wavenumber grid
    N = grid.shape[0]
    dx = grid.dx[0]
    k_vals = 2 * math.pi * np.fft.rfftfreq(N, dx)
    k_grid = [k_vals]
    rfft_shape = (N // 2 + 1,)

    # Determine IC wavenumber index
    if ic_wavevector is None:
        L = grid.bounds[0][1] - grid.bounds[0][0]
        ic_k = 2 * math.pi / L
    else:
        ic_k = ic_wavevector
    ic_k_idx = int(np.argmin(np.abs(k_vals - ic_k)))

    # Build k indices to check: IC mode + surrounding modes
    k_indices = [ic_k_idx]
    for i in range(1, n_extra_k + 1):
        if ic_k_idx + i < len(k_vals):
            k_indices.append(ic_k_idx + i)
        if ic_k_idx - i >= 0:
            k_indices.append(ic_k_idx - i)
    k_indices = sorted(set(k_indices))

    # Build constraint-eliminated system for both test and baseline
    ce = CoefficientEvaluator(spec, grid, parameters)
    A_test, _, _, _, mapping = _build_constraint_eliminated_matrices(
        spec, layout, grid, ce, k_grid, rfft_shape
    )

    baseline_params = {**parameters, **baseline_overrides}
    ce_bl = CoefficientEvaluator(spec, grid, baseline_params)
    A_baseline, _, _, _, _ = _build_constraint_eliminated_matrices(
        spec, layout, grid, ce_bl, k_grid, rfft_shape
    )

    # Extract source+target block indices
    try:
        src_r = mapping[layout.field_slot_map[source]]
        vsrc_r = mapping[layout.velocity_slot_map[source]]
        tgt_r = mapping[layout.field_slot_map[target]]
        vtgt_r = mapping[layout.velocity_slot_map[target]]
    except KeyError as exc:
        return ConversionStabilityResult(
            stable=True,
            max_excess=0.0,
            k_tachyonic=None,
            n_tachyonic_modes=0,
            message=f"Could not find source/target fields in reduced system: {exc}",
        )

    idx = np.array([src_r, vsrc_r, tgt_r, vtgt_r])

    # Check eigenvalues at each k
    max_excess = 0.0
    worst_k = None
    n_tachyonic = 0

    for ki in k_indices:
        block_test = A_test[ki][idx[:, None], idx[None, :]]
        block_bl = A_baseline[ki][idx[:, None], idx[None, :]]

        eigs_test = np.linalg.eigvals(block_test)
        eigs_bl = np.linalg.eigvals(block_bl)

        max_re_test = float(np.max(np.real(eigs_test)))
        max_re_bl = float(np.max(np.real(eigs_bl)))

        excess = max_re_test - max_re_bl
        if excess > threshold:
            n_tachyonic += 1
            if excess > max_excess:
                max_excess = excess
                worst_k = float(k_vals[ki])

    if n_tachyonic > 0:
        msg = (
            f"Tachyonic mode detected: {n_tachyonic}/{len(k_indices)} k-modes "
            f"have excess growth rate > {threshold}. "
            f"Worst: γ={max_excess:.4f} at k={worst_k:.4f}. "
            f"Conversion P will grow as exp(2γt), not oscillate. "
            f"Any measured 'amplification' A >> 1 is instability, not physics."
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
        message="Stable: no tachyonic modes in conversion channel.",
    )
