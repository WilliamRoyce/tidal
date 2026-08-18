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

For a group of fields the group spectral power is summed:

    S_group(k, omega) = sum_i |FFT(phi_i)(k, omega)|^2

This is rotationally covariant within the field group and reveals the actual
physical mode structure without depending on a particular component choice.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._spectral import (
    _build_k_grid,  # pyright: ignore[reportPrivateUsage]
    _radial_bin,  # pyright: ignore[reportPrivateUsage]
)
from tidal.measurement._utils import (
    _check_no_position_dependent_terms,  # pyright: ignore[reportPrivateUsage]
)

if TYPE_CHECKING:
    from collections.abc import Sequence

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
        Spectral power ``S(k, omega) = sum_i |A_hat_i(k, omega)|^2``
        summed over all fields in the group.
    peak_frequencies : ndarray, shape ``(n_modes,)``
        Dominant angular frequency at each k-bin (0.0 for inactive modes).
    peak_powers : ndarray, shape ``(n_modes,)``
        Spectral power at the dominant frequency per k-bin.
    field_name : str
        Which field (or comma-joined group) this dispersion was computed for.
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
    spatial_fft = np.stack(
        [np.fft.rfftn(field_snapshots[t]) for t in range(n_snapshots)],
    )

    # Temporal fft (complex input -> must use fft, not rfft).
    # Keep only strictly positive frequencies (skip DC and negative Nyquist).
    full_fft = np.fft.fft(spatial_fft, axis=0)
    n_pos = n_snapshots // 2
    spacetime_power = np.abs(full_fft[1:n_pos]) ** 2

    # Angular frequencies for the positive temporal bins
    raw_t_freqs = np.fft.fftfreq(n_snapshots, d=dt)
    angular_freqs = np.asarray(
        2.0 * np.pi * raw_t_freqs[1:n_pos],
        dtype=np.float64,
    )

    return angular_freqs, spatial_fft, spacetime_power


def _bin_and_detect(  # noqa: PLR0913, PLR0917
    angular_freqs: NDArray[np.float64],
    spatial_amp: NDArray[np.float64],
    spacetime_power: NDArray[np.float64],
    grid_shape: tuple[int, ...],
    grid_spacing: tuple[float, ...],
    min_amplitude: float,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Radially bin S(k, omega) and detect peaks.

    Parameters
    ----------
    spatial_amp : ndarray, shape ``(n_snapshots, *rfft_shape)``
        Spatial amplitude array (absolute values, already real).  For a
        single field this is ``|spatial_fft|``; for a group it is the
        element-wise sum of ``|spatial_fft_i|`` over all fields.

    Returns (wavenumbers, power, peak_frequencies, peak_powers).
    """
    _k_grid, k_mag = _build_k_grid(grid_shape, grid_spacing)

    # Establish bin structure
    wn_ref, _ = _radial_bin(k_mag, np.zeros_like(k_mag), grid_spacing, grid_shape)
    n_modes = len(wn_ref)
    n_freq = len(angular_freqs)
    power = np.zeros((n_modes, n_freq), dtype=np.float64)

    for fi in range(n_freq):
        _, binned = _radial_bin(k_mag, spacetime_power[fi], grid_spacing, grid_shape)
        power[: len(binned), fi] = binned

    # Max spatial amplitude per k-bin for activity detection
    max_amp = np.zeros(n_modes, dtype=np.float64)
    for t in range(spatial_amp.shape[0]):
        _, binned_amp = _radial_bin(k_mag, spatial_amp[t], grid_spacing, grid_shape)
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


def compute_dispersion(  # noqa: PLR0914
    data: SimulationData,
    field_names: str | Sequence[str],
    *,
    min_amplitude: float = 1e-12,
) -> DispersionResult:
    """Extract dispersion relation omega(k) from simulation output.

    Algorithm
    ---------
    1. For each field in *field_names*, compute spatial ``rfftn`` per
       snapshot to get complex Fourier coefficients ``phi_hat_i(k, t)``.
    2. For each spatial mode ``k``, temporal FFT of the **complex**
       coefficient gives ``S_i(k, omega)``.
    3. Sum spectral power across all fields:
       ``S_group(k, omega) = sum_i S_i(k, omega)``.
    4. Radially bin ``S_group(k, omega)`` into ``|k|`` shells.
    5. Peak detection: ``argmax(S_group)`` per k-bin extracts ``omega(k)``.
    6. Modes with combined max amplitude below threshold are inactive.

    Using complex coefficients (not ``|phi_hat|``) avoids frequency
    doubling artifacts from taking the absolute value before FFT.
    Summing power over a field group makes the measurement rotationally
    covariant within the group (no dependence on which single component
    is selected).

    Parameters
    ----------
    data : SimulationData
        Simulation output with time-resolved field snapshots.
    field_names : str or sequence of str
        Field or group of fields to extract the dispersion relation for.
        All fields must be dynamical (``time_derivative_order >= 2``);
        constraint fields raise ``ValueError``.
    min_amplitude : float, optional
        Modes with combined max ``|phi_hat(k, t)|`` below this threshold
        are treated as inactive.  Default ``1e-12``.

    Returns
    -------
    DispersionResult

    Raises
    ------
    ValueError
        If any field is unknown, is a constraint field
        (``time_derivative_order < 2``), fewer than 3 snapshots, the
        timestep is non-uniform, or any equation term is
        position-dependent (uniform medium required).
    """
    _check_no_position_dependent_terms(data, "Dispersion relation omega(k)")

    # Normalize to list
    if isinstance(field_names, str):
        names: list[str] = [field_names]
    else:
        names = list(field_names)

    if not names:
        msg = "field_names must not be empty"
        raise ValueError(msg)

    # Validate all field names and check they are dynamical
    for fname in names:
        if fname not in data.spec.component_names:
            msg = f"Field '{fname}' not in spec fields: {data.spec.component_names}"
            raise ValueError(msg)
        eq = next(e for e in data.spec.equations if e.field_name == fname)
        if eq.time_derivative_order < 2:  # noqa: PLR2004
            msg = (
                f"Field '{fname}' is a constraint (time_derivative_order="
                f"{eq.time_derivative_order}). Dispersion requires a dynamical "
                f"field (time_derivative_order >= 2)."
            )
            raise ValueError(msg)

    _validate_timestep(data)

    dt = float(data.times[1] - data.times[0])
    rayleigh = 2.0 * np.pi / float(data.times[-1] - data.times[0])

    # First field: seed the accumulators
    first_snapshots = data.fields[names[0]]
    angular_freqs, first_sfft, spacetime_power_total = _spacetime_fft(
        first_snapshots,
        dt,
    )
    spatial_amp_sum = np.abs(first_sfft)
    grid_shape = first_snapshots.shape[1:]

    # Remaining fields: add their spectral power and amplitude
    for fname in names[1:]:
        field_snapshots = data.fields[fname]
        _, sfft, power = _spacetime_fft(field_snapshots, dt)
        spacetime_power_total += power
        spatial_amp_sum += np.abs(sfft)

    wn, power_binned, peak_freqs, peak_pow = _bin_and_detect(
        angular_freqs,
        spatial_amp_sum,
        spacetime_power_total,
        grid_shape,
        data.grid_spacing,
        min_amplitude,
    )

    return DispersionResult(
        wavenumbers=wn,
        frequencies=angular_freqs,
        power=power_binned,
        peak_frequencies=peak_freqs,
        peak_powers=peak_pow,
        field_name=", ".join(names),
        rayleigh_resolution=rayleigh,
    )
