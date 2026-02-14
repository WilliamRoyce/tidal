"""Mixing length and mixing spectrum extraction from conversion probability.

Model-independent measurement of the characteristic energy-exchange timescale
in coupled field systems.  The mixing length ``L_mix`` is derived from the
dominant peak of the temporal FFT of ``P(t)`` — the half-period ``π/ω`` of
the strongest oscillation frequency.

This spectral approach correctly identifies the physically meaningful
oscillation timescale even in multi-scale systems where rapid noise
oscillations sit on top of a slower mixing envelope.  Uncertainty is
estimated from the full width at half maximum (FWHM) of the spectral peak.

The mixing spectrum is the temporal FFT of ``P(t)``, showing which oscillation
frequencies participate in the energy exchange.  Users interpret the spectrum
through their own theory (Gertsenshtein, axion-photon, Rabi, etc.).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._conversion import ConversionResult


@dataclass(frozen=True)
class SpectralPeak:
    """A detected peak in the mixing power spectrum.

    Each peak represents a frequency at which ``P(t)`` oscillates.
    The mixing length ``π/ω`` gives the half-period of energy exchange
    at that frequency.  FWHM measures how sharply defined the
    oscillation is — a narrow peak means a coherent, well-defined
    oscillation; a broad peak means the frequency is less certain.

    Attributes
    ----------
    frequency : float
        Angular frequency ``ω`` (rad/time) of the spectral peak.
    power : float
        ``|P̂(ω)|²`` — spectral power at this frequency.
    mixing_length : float
        ``π/ω`` — half-period of energy exchange at this frequency.
    fwhm : float
        Full width at half maximum of the peak (rad/time).
    mixing_length_uncertainty : float
        Propagated uncertainty: ``(pi/omega**2) * fwhm``.
    """

    frequency: float
    power: float
    mixing_length: float
    fwhm: float
    mixing_length_uncertainty: float


@dataclass(frozen=True)
class MixingResult:
    """Result of spectral mixing length extraction.

    The mixing length is derived from the dominant peak of the temporal
    FFT of ``P(t)``, rather than from time-domain peak detection.  This
    correctly identifies the physically meaningful oscillation timescale
    even in multi-scale systems where rapid noise oscillations sit on
    top of a slower mixing envelope.

    The uncertainty comes from the FWHM of the spectral peak:
    ``dL = (pi/omega**2) * d_omega``, where ``d_omega`` is the FWHM.  A sharp spectral
    peak (small FWHM) gives a precise mixing length; a broad peak
    indicates the oscillation frequency is less well-defined.

    Attributes
    ----------
    mixing_length : float
        ``π/ω_dom`` — half-period of the dominant oscillation frequency.
    mixing_length_uncertainty : float
        Propagated from FWHM of the dominant spectral peak.
    dominant_frequency : float
        ``ω_dom`` — angular frequency of the strongest spectral peak.
    frequency_fwhm : float
        FWHM of the dominant peak (rad/time).
    max_conversion : float
        ``max(P(t))`` — peak conversion probability over the full timeseries.
    peaks : tuple of SpectralPeak
        All detected spectral peaks, sorted by power descending.
        ``peaks[0]`` is the dominant peak.
    """

    mixing_length: float
    mixing_length_uncertainty: float
    dominant_frequency: float
    frequency_fwhm: float
    max_conversion: float
    peaks: tuple[SpectralPeak, ...]


@dataclass(frozen=True)
class MixingSpectrum:
    """Frequency decomposition of the conversion probability timeseries.

    Reports angular frequencies at which ``P(t)`` oscillates — no
    theory-specific conversion.  For each frequency ``ω``, the half-period
    ``π/ω`` gives the mixing timescale at that frequency.

    Attributes
    ----------
    frequencies : ndarray
        Angular frequencies ``ω`` (rad/time), excluding DC.
    power : ndarray
        ``|P̂(ω)|²`` at each frequency.
    dominant_frequency : float
        ``ω`` of the strongest oscillation peak.
    dominant_mixing_length : float
        ``π / dominant_frequency`` — half-period at the dominant frequency.
    """

    frequencies: NDArray[np.float64]
    power: NDArray[np.float64]
    dominant_frequency: float
    dominant_mixing_length: float


_MIN_POINTS = 3  # Minimum time points for peak detection / FFT

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _compute_fwhm(
    freqs: NDArray[np.float64],
    power: NDArray[np.float64],
    peak_idx: int,
) -> float:
    """Full width at half maximum of a spectral peak.

    Uses linear interpolation between bins for sub-bin accuracy.
    If the peak is at a spectrum boundary, uses the available one-sided
    width doubled as an estimate.
    """
    half_max = float(power[peak_idx]) / 2.0
    n = len(freqs)

    # --- Scan left ---
    left_freq = float(freqs[0])  # default: left edge
    for j in range(peak_idx - 1, -1, -1):
        if power[j] < half_max:
            # Linear interpolation between j and j+1
            denom = float(power[j + 1] - power[j])
            if abs(denom) > 0.0:
                frac = (half_max - float(power[j])) / denom
                left_freq = float(freqs[j]) + frac * (
                    float(freqs[j + 1]) - float(freqs[j])
                )
            else:
                left_freq = float(freqs[j])
            break

    # --- Scan right ---
    right_freq = float(freqs[-1])  # default: right edge
    for j in range(peak_idx + 1, n):
        if power[j] < half_max:
            # Linear interpolation between j-1 and j
            denom = float(power[j - 1] - power[j])
            if abs(denom) > 0.0:
                frac = (half_max - float(power[j])) / denom
                right_freq = float(freqs[j]) - frac * (
                    float(freqs[j]) - float(freqs[j - 1])
                )
            else:
                right_freq = float(freqs[j])
            break

    fwhm = right_freq - left_freq

    # If only one side was found (peak at boundary), double the half-width
    if peak_idx == 0:
        fwhm = 2.0 * (right_freq - float(freqs[peak_idx]))
    elif peak_idx == n - 1:
        fwhm = 2.0 * (float(freqs[peak_idx]) - left_freq)

    # Floor: at least the frequency resolution
    min_bins_for_resolution = 2
    if n >= min_bins_for_resolution:
        df = float(freqs[1] - freqs[0])
        fwhm = max(fwhm, df)

    return fwhm


def _find_spectral_peaks(
    freqs: NDArray[np.float64],
    power: NDArray[np.float64],
    min_prominence: float,
) -> tuple[SpectralPeak, ...]:
    """Find peaks in the power spectrum above the prominence threshold.

    Returns peaks sorted by power descending.
    """
    if len(power) < _MIN_POINTS:
        return ()

    max_power = float(np.max(power))
    if max_power <= 0.0:
        return ()

    threshold = min_prominence * max_power

    # Find local maxima
    peak_indices: list[int] = [
        i
        for i in range(1, len(power) - 1)
        if power[i] > power[i - 1]
        and power[i] > power[i + 1]
        and float(power[i]) > threshold
    ]

    # Build SpectralPeak for each
    peaks: list[SpectralPeak] = []
    for idx in peak_indices:
        omega = float(freqs[idx])
        p = float(power[idx])
        fwhm = _compute_fwhm(freqs, power, idx)
        ml = np.pi / omega
        ml_unc = (np.pi / (omega * omega)) * fwhm
        peaks.append(
            SpectralPeak(
                frequency=omega,
                power=p,
                mixing_length=ml,
                fwhm=fwhm,
                mixing_length_uncertainty=ml_unc,
            )
        )

    # Sort by power descending
    peaks.sort(key=lambda pk: pk.power, reverse=True)
    return tuple(peaks)


# ------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------


def compute_mixing_length(
    conversion: ConversionResult,
    *,
    min_prominence: float = 0.01,
) -> MixingResult:
    """Extract the characteristic mixing length from the power spectrum of ``P(t)``.

    Computes the temporal FFT of the conversion probability ``P(t)``, finds
    spectral peaks, and derives the mixing length from the dominant peak.
    This spectral approach correctly identifies the physically meaningful
    oscillation timescale even in multi-scale systems where rapid
    oscillations sit on top of a slower mixing envelope.

    The uncertainty comes from the full width at half maximum (FWHM) of
    the spectral peak.  Error propagation: ``L = pi/omega``,
    ``dL = (pi/omega**2) * d_omega``.

    Parameters
    ----------
    conversion : ConversionResult
        Output of :func:`compute_conversion_probability` or
        :func:`compute_group_conversion`.
    min_prominence : float, optional
        Minimum spectral peak power as a fraction of the maximum power.
        Peaks below this threshold are treated as noise.  Default ``0.01``
        (1% of dominant peak power).

    Returns
    -------
    MixingResult

    Raises
    ------
    ValueError
        If fewer than 3 time points, if timestep is non-uniform, if
        ``min_prominence`` is not in ``(0, 1)``, or if no spectral peaks
        are found above the prominence threshold.
    """
    if not (0.0 < min_prominence < 1.0):
        msg = f"min_prominence must be in (0, 1), got {min_prominence}"
        raise ValueError(msg)

    # Compute the spectrum (validates time points + uniform timestep)
    spectrum = compute_mixing_spectrum(conversion)

    # Find spectral peaks
    peaks = _find_spectral_peaks(spectrum.frequencies, spectrum.power, min_prominence)

    if not peaks:
        msg = (
            "No spectral peaks found above prominence threshold — "
            "P(t) may be monotonic, too short, or coupling too weak "
            "for observable oscillation"
        )
        raise ValueError(msg)

    dominant = peaks[0]

    return MixingResult(
        mixing_length=dominant.mixing_length,
        mixing_length_uncertainty=dominant.mixing_length_uncertainty,
        dominant_frequency=dominant.frequency,
        frequency_fwhm=dominant.fwhm,
        max_conversion=float(np.max(conversion.probability)),
        peaks=peaks,
    )


def compute_mixing_spectrum(
    conversion: ConversionResult,
) -> MixingSpectrum:
    """Compute the frequency decomposition of ``P(t)``.

    Returns the power spectrum of the conversion probability timeseries,
    showing which oscillation frequencies participate in the energy exchange.
    This is the temporal analog of spatial spectral decomposition — answering
    "at what timescales does energy transfer occur?"

    Angular frequencies are reported directly (consistent with
    :func:`compute_spectrum` which reports spatial wavenumbers as angular
    frequencies).  No theory-specific conversion is applied.  For each
    frequency ``ω``, the half-period ``π/ω`` gives the corresponding
    mixing timescale.

    Parameters
    ----------
    conversion : ConversionResult
        Output of :func:`compute_conversion_probability` or
        :func:`compute_group_conversion`.

    Returns
    -------
    MixingSpectrum

    Raises
    ------
    ValueError
        If the timestep is non-uniform or if fewer than 3 time points.
    """
    times = conversion.times
    prob = conversion.probability

    if len(times) < _MIN_POINTS:
        msg = "Need at least 3 time points for spectral analysis"
        raise ValueError(msg)

    # Validate uniform timestep
    dt = float(times[1] - times[0])
    diffs = np.diff(times)
    if not np.allclose(diffs, dt, rtol=1e-6):
        msg = (
            "Non-uniform timestep — FFT requires uniform sampling "
            f"(dt range: {float(diffs.min()):.6g} to {float(diffs.max()):.6g})"
        )
        raise ValueError(msg)

    # Subtract mean to remove DC component
    prob_centered = prob - np.mean(prob)

    # FFT → power spectrum
    fft_vals = np.fft.rfft(prob_centered)
    raw_freqs = np.fft.rfftfreq(len(prob_centered), d=dt)

    # Skip DC bin (index 0), convert to angular frequency
    angular_freqs = np.asarray(2.0 * np.pi * raw_freqs[1:], dtype=np.float64)
    power = np.asarray(np.abs(fft_vals[1:]) ** 2, dtype=np.float64)

    if len(angular_freqs) == 0:
        msg = "Not enough frequency bins — need more time points"
        raise ValueError(msg)

    # Dominant peak
    peak_idx = int(np.argmax(power))
    omega_dom = float(angular_freqs[peak_idx])

    return MixingSpectrum(
        frequencies=angular_freqs,
        power=power,
        dominant_frequency=omega_dom,
        dominant_mixing_length=np.pi / omega_dom,
    )
