"""Resonance analysis for coupled-field systems.

Identifies resonant modes where ``omega_source(k) ~ omega_target(k)`` and
computes conversion bandwidth, peak conversion wavenumber, and the
adiabaticity parameter.

Requires flat + spatially homogeneous systems (inherits from dispersion).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._dispersion import DispersionResult, compute_dispersion

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class ResonanceResult:
    """Resonance analysis for coupled fields.

    Attributes
    ----------
    wavenumbers : ndarray, shape ``(n_shared,)``
        Shared wavenumber bins where both fields have active modes.
    omega_source : ndarray, shape ``(n_shared,)``
        Source field frequencies at shared k-bins.
    omega_target : ndarray, shape ``(n_shared,)``
        Target field frequencies at shared k-bins.
    resonance_mismatch : ndarray, shape ``(n_shared,)``
        ``|omega_s - omega_t| / omega_avg`` per shared mode.
    resonant_modes : ndarray, shape ``(n_shared,)``
        Boolean mask: which modes are near-resonant.
    n_resonant_modes : int
        Count of resonant modes.
    conversion_bandwidth : float
        FWHM of the resonance mismatch curve in k-space (0.0 if < 2 resonant).
    peak_conversion_k : float
        Wavenumber with minimum mismatch (closest to exact resonance).
    source_field : str
        Source field name.
    target_field : str
        Target field name.
    """

    wavenumbers: NDArray[np.float64]
    omega_source: NDArray[np.float64]
    omega_target: NDArray[np.float64]
    resonance_mismatch: NDArray[np.float64]
    resonant_modes: NDArray[np.bool_]
    n_resonant_modes: int
    conversion_bandwidth: float
    peak_conversion_k: float
    source_field: str
    target_field: str


def _compute_bandwidth(
    wavenumbers: NDArray[np.float64],
    mismatch: NDArray[np.float64],
) -> float:
    """Compute FWHM of the resonance window in k-space.

    The "resonance window" is defined by inverting the mismatch: regions of
    low mismatch are the resonance peak. We find the width at half-minimum.
    """
    if len(wavenumbers) < 2:  # noqa: PLR2004
        return 0.0

    # Invert: peak where mismatch is minimum
    min_mismatch = mismatch.min()
    max_mismatch = mismatch.max()
    if max_mismatch <= min_mismatch:
        return 0.0

    half_level = (min_mismatch + max_mismatch) / 2.0

    # Find where mismatch crosses half_level (going from below to above)
    below_half = mismatch < half_level
    if not np.any(below_half):
        return 0.0

    k_below = wavenumbers[below_half]
    return float(k_below[-1] - k_below[0])


def _find_shared_modes(
    disp_src: DispersionResult,
    disp_tgt: DispersionResult,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Find wavenumber bins where both fields have active modes.

    Returns
    -------
    tuple
        (shared_k, omega_source, omega_target) arrays.

    Raises
    ------
    ValueError
        If no shared active modes found.
    """
    active_src = disp_src.peak_frequencies > 0.0
    active_tgt = disp_tgt.peak_frequencies > 0.0

    idx_src: list[int] = []
    idx_tgt: list[int] = []
    for i, k_s in enumerate(disp_src.wavenumbers):
        if not active_src[i]:
            continue
        matches = np.where(
            np.isclose(disp_tgt.wavenumbers, k_s, rtol=1e-6) & active_tgt
        )[0]
        if len(matches) > 0:
            idx_src.append(i)
            idx_tgt.append(matches[0])

    if not idx_src:
        msg = "No shared active wavenumber modes between source and target fields"
        raise ValueError(msg)

    return (
        disp_src.wavenumbers[idx_src],
        disp_src.peak_frequencies[idx_src],
        disp_tgt.peak_frequencies[idx_tgt],
    )


def compute_resonance_analysis(
    data: SimulationData,
    source_field: str | Sequence[str],
    target_field: str | Sequence[str],
    *,
    resonance_threshold: float = 0.1,
    min_amplitude: float = 1e-12,
) -> ResonanceResult:
    """Identify resonant modes and compute conversion bandwidth.

    A mode k is resonant when
    ``|omega_source(k) - omega_target(k)| / omega_avg(k) < threshold``.

    Parameters
    ----------
    data : SimulationData
        Simulation output.
    source_field : str or sequence of str
        Source field name(s).
    target_field : str or sequence of str
        Target field name(s).
    resonance_threshold : float, optional
        Relative frequency mismatch threshold for resonance.
    min_amplitude : float, optional
        Minimum Fourier amplitude for active modes.

    Returns
    -------
    ResonanceResult
    """
    disp_src = compute_dispersion(data, source_field, min_amplitude=min_amplitude)
    disp_tgt = compute_dispersion(data, target_field, min_amplitude=min_amplitude)

    shared_k, omega_s, omega_t = _find_shared_modes(disp_src, disp_tgt)

    # Relative frequency mismatch
    omega_avg = (omega_s + omega_t) / 2.0
    safe_avg = np.where(omega_avg > 0, omega_avg, 1.0)
    mismatch = np.abs(omega_s - omega_t) / safe_avg

    resonant = mismatch < resonance_threshold
    peak_k = float(shared_k[int(np.argmin(mismatch))])
    bandwidth = _compute_bandwidth(shared_k, mismatch)

    src_name = (
        source_field if isinstance(source_field, str) else ", ".join(source_field)
    )
    tgt_name = (
        target_field if isinstance(target_field, str) else ", ".join(target_field)
    )

    return ResonanceResult(
        wavenumbers=shared_k,
        omega_source=omega_s,
        omega_target=omega_t,
        resonance_mismatch=mismatch,
        resonant_modes=resonant,
        n_resonant_modes=int(resonant.sum()),
        conversion_bandwidth=bandwidth,
        peak_conversion_k=peak_k,
        source_field=src_name,
        target_field=tgt_name,
    )
