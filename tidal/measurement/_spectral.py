"""Fourier spectral decomposition of field data.

Decomposes fields into wavenumber (k-space) components for tracking
which Fourier modes participate in wave conversion.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._energy import (
    _validate_array,  # pyright: ignore[reportPrivateUsage]
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class SpectralSnapshot:
    """Spectral decomposition of a field at one time.

    Attributes
    ----------
    wavenumbers : ndarray
        Radially binned wavenumber magnitudes ``|k|``.
    power_spectrum : ndarray
        ``|φ̂(k)|²`` averaged over shells of constant ``|k|``.
    """

    wavenumbers: NDArray[np.float64]
    power_spectrum: NDArray[np.float64]


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _apply_hann_window(
    field_data: NDArray[np.float64],
    periodic: tuple[bool, ...],
) -> NDArray[np.float64]:
    """Apply Hann window on non-periodic axes; emit warning if needed."""
    data = field_data.copy()
    any_non_periodic = False
    for axis in range(len(periodic)):
        if not periodic[axis]:
            any_non_periodic = True
            n = data.shape[axis]
            window = np.hanning(n)
            shape = [1] * data.ndim
            shape[axis] = n
            data *= window.reshape(shape)
    if any_non_periodic:
        warnings.warn(
            "Applying Hann window on non-periodic axes; spectral amplitudes "
            "are approximate (reduced by windowing factor).",
            UserWarning,
            stacklevel=3,
        )
    return data


def _build_k_grid(
    field_shape: tuple[int, ...],
    grid_spacing: tuple[float, ...],
) -> tuple[list[NDArray[np.float64]], NDArray[np.float64]]:
    """Build wavenumber arrays and ``|k|`` magnitude grid."""
    k_arrays: list[NDArray[np.float64]] = []
    ndim = len(grid_spacing)
    for axis in range(ndim):
        n = field_shape[axis]
        dx = grid_spacing[axis]
        if axis == ndim - 1:
            freq = np.asarray(np.fft.rfftfreq(n, d=dx) * (2.0 * np.pi), dtype=np.float64)
        else:
            freq = np.asarray(np.fft.fftfreq(n, d=dx) * (2.0 * np.pi), dtype=np.float64)
        k_arrays.append(freq)
    k_grid: list[NDArray[np.float64]] = [
        np.asarray(g, dtype=np.float64) for g in np.meshgrid(*k_arrays, indexing="ij")
    ]
    k_mag = np.sqrt(np.asarray(sum(ki**2 for ki in k_grid), dtype=np.float64))
    return k_grid, k_mag


def _radial_bin(
    k_mag: NDArray[np.float64],
    values: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    field_shape: tuple[int, ...],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Radially bin *values* by ``|k|`` magnitude."""
    k_max = float(np.max(k_mag))
    if k_max == 0.0:
        return np.array([0.0]), np.array([float(values.sum())])

    dk = min(
        2.0 * np.pi / (field_shape[ax] * grid_spacing[ax])
        for ax in range(len(grid_spacing))
    )
    n_bins = max(1, int(np.ceil(k_max / dk)))
    bin_edges = np.linspace(0.0, k_max + dk, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    binned = np.zeros(n_bins, dtype=np.float64)
    bin_indices = np.clip(np.digitize(k_mag.ravel(), bin_edges) - 1, 0, n_bins - 1)
    v_flat = values.ravel()
    for b_idx in range(n_bins):
        mask = bin_indices == b_idx
        binned[b_idx] = float(v_flat[mask].sum())
    return bin_centers, binned


# ------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------


def compute_spectrum(
    field_data: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
) -> SpectralSnapshot:
    """Compute the radially-averaged power spectrum of a field snapshot.

    Parameters
    ----------
    field_data : ndarray, shape ``(*grid_shape)``
        Real-valued field on the spatial grid.
    grid_spacing : tuple[float, ...]
        Cell sizes per axis.
    periodic : tuple[bool, ...]
        Per-axis periodicity.

    Returns
    -------
    SpectralSnapshot
    """
    _validate_array(field_data, "field_data")

    data_for_fft = _apply_hann_window(field_data, periodic)
    fhat = np.fft.rfftn(data_for_fft)
    power_full = np.abs(fhat) ** 2

    _k_grid, k_mag = _build_k_grid(field_data.shape, grid_spacing)
    wavenumbers, binned = _radial_bin(k_mag, power_full, grid_spacing, field_data.shape)
    return SpectralSnapshot(wavenumbers=wavenumbers, power_spectrum=binned)


def compute_spectral_energy(
    field_data: NDArray[np.float64],
    momentum_data: NDArray[np.float64] | None,
    mass_squared: float,
    grid_spacing: tuple[float, ...],
    _periodic: tuple[bool, ...],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute per-mode energy ``E(k) = 0.5 * [|π̂_k|² + (k² + m²)|φ̂_k|²]``.

    Returns the radially-averaged spectral energy density.

    Parameters
    ----------
    field_data : ndarray
        Field snapshot.
    momentum_data : ndarray or None
        Momentum snapshot (None for constraint fields).
    mass_squared : float
        Mass-squared term.
    grid_spacing : tuple[float, ...]
        Cell sizes per axis.
    _periodic : tuple[bool, ...]
        Per-axis periodicity (reserved for future windowing).

    Returns
    -------
    wavenumbers : ndarray
        Radially binned ``|k|`` values.
    spectral_energy : ndarray
        Energy per wavenumber bin.
    """
    _validate_array(field_data, "field_data")

    phi_hat = np.fft.rfftn(field_data)
    k_grid, k_mag = _build_k_grid(field_data.shape, grid_spacing)
    k_sq = sum(ki**2 for ki in k_grid)

    n_total = float(np.array(field_data.shape).prod())
    energy_field = 0.5 * (k_sq + mass_squared) * np.abs(phi_hat) ** 2 / n_total

    if momentum_data is not None:
        _validate_array(momentum_data, "momentum_data")
        pi_hat = np.fft.rfftn(momentum_data)
        energy_total = energy_field + 0.5 * np.abs(pi_hat) ** 2 / n_total
    else:
        energy_total = energy_field

    return _radial_bin(k_mag, energy_total, grid_spacing, field_data.shape)


def compute_mode_amplitudes(
    data: SimulationData,
    field_name: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Track mode amplitudes ``|φ̂(k)|`` over time.

    Parameters
    ----------
    data : SimulationData
    field_name : str
        Which field to decompose.

    Returns
    -------
    times : ndarray, shape ``(n_snapshots,)``
    wavenumbers : ndarray, shape ``(n_modes,)``
    amplitudes : ndarray, shape ``(n_snapshots, n_modes)``
        ``|FFT coefficient|`` at each ``(time, |k|)`` bin.

    Raises
    ------
    ValueError
        If *field_name* is not in the spec.
    """
    if field_name not in data.spec.component_names:
        msg = f"Field '{field_name}' not in spec fields: {data.spec.component_names}"
        raise ValueError(msg)

    # Compute spectrum at first snapshot to get wavenumber bins
    first = compute_spectrum(
        data.fields[field_name][0], data.grid_spacing, data.periodic,
    )
    n_modes = len(first.wavenumbers)

    amplitudes = np.zeros((data.n_snapshots, n_modes), dtype=np.float64)
    amplitudes[0] = np.sqrt(first.power_spectrum)

    for t_idx in range(1, data.n_snapshots):
        snap = compute_spectrum(
            data.fields[field_name][t_idx], data.grid_spacing, data.periodic,
        )
        amplitudes[t_idx] = np.sqrt(snap.power_spectrum)

    return data.times.copy(), first.wavenumbers.copy(), amplitudes
