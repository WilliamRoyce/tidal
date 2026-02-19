"""Dispersion relation extraction from simulation data.

Extracts the dispersion relation omega(k) by performing a temporal FFT on
the spatial Fourier coefficients phi_hat(k, t) for each spatial mode.  The
resulting 2D power spectrum S(k, omega) reveals how energy is distributed
across wavenumber and frequency, with peak detection extracting the dominant
frequency per k-bin to form the dispersion curve.

The key insight is to use complex spatial FFT coefficients (not amplitudes
|phi_hat|) to avoid frequency doubling artifacts.  The temporal FFT of the
complex coefficient phi_hat(k, t) correctly identifies the oscillation
frequency omega(k).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._spectral import (
    build_k_grid,
    radial_bin,
)
from tidal.measurement._utils import (
    check_no_position_dependent_terms,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


_MIN_SNAPSHOTS = 3


@dataclass(frozen=True)
class DispersionResult:
    """Dispersion relation extracted from simulation output.

    Attributes
    ----------
    wavenumbers : ndarray, shape ``(n_modes,)``
        Radially binned wavenumber magnitudes ``|k|``.
    frequencies : ndarray, shape ``(n_freq,)``
        Angular frequencies ``omega`` (rad/time), excluding DC.
    power : ndarray, shape ``(n_modes, n_freq)``
        Spectral power ``S(k, omega) = |A_hat(k, omega)|^2``.
    peak_frequencies : ndarray, shape ``(n_modes,)``
        Dominant angular frequency at each k-bin (0.0 for inactive modes).
    peak_powers : ndarray, shape ``(n_modes,)``
        Spectral power at the dominant frequency per k-bin.
    field_name : str
        Which field this dispersion was computed for.
    rayleigh_resolution : float
        ``2*pi/T`` -- minimum resolvable frequency difference.
    """

    wavenumbers: NDArray[np.float64]
    frequencies: NDArray[np.float64]
    power: NDArray[np.float64]
    peak_frequencies: NDArray[np.float64]
    peak_powers: NDArray[np.float64]
    field_name: str
    rayleigh_resolution: float


def _validate_timestep(data: SimulationData) -> None:
    """Check snapshot count and uniform timestep.

    Raises
    ------
    ValueError
        If fewer than 3 snapshots or non-uniform timestep.
    """
    if data.n_snapshots < _MIN_SNAPSHOTS:
        msg = (
            f"Need at least {_MIN_SNAPSHOTS} snapshots for dispersion "
            f"extraction, got {data.n_snapshots}"
        )
        raise ValueError(msg)

    times = data.times
    dt = float(times[1] - times[0])
    diffs = np.diff(times)
    if not np.allclose(diffs, dt, rtol=1e-6):
        msg = (
            "Non-uniform timestep -- temporal FFT requires uniform sampling "
            f"(dt range: {float(diffs.min()):.6g} to {float(diffs.max()):.6g})"
        )
        raise ValueError(msg)


def _spacetime_fft(
    field_snapshots: NDArray[np.float64],
    dt: float,
) -> tuple[NDArray[np.float64], NDArray[np.complex128], NDArray[np.float64]]:
    """Spatial + temporal FFT pipeline.

    Returns (angular_freqs, spatial_fft, spacetime_power) where
    spacetime_power has shape ``(n_freq, *rfft_shape)`` and
    spatial_fft has shape ``(n_snapshots, *rfft_shape)``.
    """
    n_snapshots = field_snapshots.shape[0]

    # Spatial rfftn per snapshot -> complex coefficients
    spatial_fft = np.stack([
        np.fft.rfftn(field_snapshots[t]) for t in range(n_snapshots)
    ])

    # Temporal fft (complex input -> must use fft, not rfft).
    # Keep only strictly positive frequencies (skip DC and negative Nyquist).
    full_fft = np.fft.fft(spatial_fft, axis=0)
    n_pos = n_snapshots // 2
    spacetime_power = np.abs(full_fft[1:n_pos]) ** 2

    # Angular frequencies for the positive temporal bins
    raw_t_freqs = np.fft.fftfreq(n_snapshots, d=dt)
    angular_freqs = np.asarray(
        2.0 * np.pi * raw_t_freqs[1:n_pos], dtype=np.float64,
    )

    return angular_freqs, spatial_fft, spacetime_power


def _bin_and_detect(  # noqa: PLR0913, PLR0917
    angular_freqs: NDArray[np.float64],
    spatial_fft: NDArray[np.complex128],
    spacetime_power: NDArray[np.float64],
    grid_shape: tuple[int, ...],
    grid_spacing: tuple[float, ...],
    min_amplitude: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Radially bin S(k, omega) and detect peaks.

    Returns (wavenumbers, power, peak_frequencies, peak_powers).
    """
    _k_grid, k_mag = build_k_grid(grid_shape, grid_spacing)

    # Establish bin structure
    wn_ref, _ = radial_bin(k_mag, np.zeros_like(k_mag), grid_spacing, grid_shape)
    n_modes = len(wn_ref)
    n_freq = len(angular_freqs)
    power = np.zeros((n_modes, n_freq), dtype=np.float64)

    for fi in range(n_freq):
        _, binned = radial_bin(k_mag, spacetime_power[fi], grid_spacing, grid_shape)
        power[:len(binned), fi] = binned

    # Max spatial amplitude per k-bin for activity detection
    max_amp = np.zeros(n_modes, dtype=np.float64)
    for t in range(spatial_fft.shape[0]):
        _, binned_amp = radial_bin(k_mag, np.abs(spatial_fft[t]), grid_spacing, grid_shape)
        max_amp = np.maximum(max_amp, binned_amp[:n_modes])

    # Peak detection per k-bin
    peak_frequencies = np.zeros(n_modes, dtype=np.float64)
    peak_powers = np.zeros(n_modes, dtype=np.float64)

    for ki in range(n_modes):
        if max_amp[ki] < min_amplitude:
            continue
        mode_power = power[ki, :]
        if np.max(mode_power) > 0.0:
            peak_idx = int(np.argmax(mode_power))
            peak_frequencies[ki] = angular_freqs[peak_idx]
            peak_powers[ki] = mode_power[peak_idx]

    return wn_ref, power, peak_frequencies, peak_powers


def compute_dispersion(
    data: SimulationData,
    field_name: str,
    *,
    min_amplitude: float = 1e-12,
) -> DispersionResult:
    """Extract dispersion relation omega(k) from simulation output.

    Algorithm
    ---------
    1. For each snapshot, compute spatial ``rfftn`` to get complex
       Fourier coefficients ``phi_hat(k, t)``.
    2. For each spatial mode ``k``, temporal FFT of the **complex**
       coefficient gives ``S(k, omega)``.
    3. Radially bin ``S(k, omega)`` into ``|k|`` shells.
    4. Peak detection: ``argmax(S)`` per k-bin extracts ``omega(k)``.
    5. Modes with max power below threshold are marked inactive.

    Using complex coefficients (not ``|phi_hat|``) avoids frequency
    doubling artifacts from taking the absolute value before FFT.

    Parameters
    ----------
    data : SimulationData
        Simulation output with time-resolved field snapshots.
    field_name : str
        Which field to extract the dispersion relation for.
    min_amplitude : float, optional
        Modes with max ``|phi_hat(k, t)|`` below this threshold are
        treated as inactive.  Default ``1e-12``.

    Returns
    -------
    DispersionResult

    Raises
    ------
    ValueError
        If *field_name* is unknown, fewer than 3 snapshots, the
        timestep is non-uniform, or any equation term is
        position-dependent (uniform medium required).
    """
    check_no_position_dependent_terms(data, "Dispersion relation omega(k)")

    if field_name not in data.spec.component_names:
        msg = f"Field '{field_name}' not in spec fields: {data.spec.component_names}"
        raise ValueError(msg)

    _validate_timestep(data)

    field_snapshots = data.fields[field_name]
    dt = float(data.times[1] - data.times[0])
    rayleigh = 2.0 * np.pi / float(data.times[-1] - data.times[0])

    angular_freqs, spatial_fft, spacetime_power = _spacetime_fft(field_snapshots, dt)
    wn, power, peak_freqs, peak_pow = _bin_and_detect(
        angular_freqs, spatial_fft, spacetime_power,
        field_snapshots.shape[1:], data.grid_spacing, min_amplitude,
    )

    return DispersionResult(
        wavenumbers=wn,
        frequencies=angular_freqs,
        power=power,
        peak_frequencies=peak_freqs,
        peak_powers=peak_pow,
        field_name=field_name,
        rayleigh_resolution=rayleigh,
    )
