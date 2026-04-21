"""Mixing length and mixing spectrum extraction from conversion probability.

Model-independent measurement of the characteristic energy-exchange timescale
in coupled field systems.  The mixing length ``L_mix`` is derived from the
dominant peak of the temporal FFT of ``P(t)`` — the half-period ``pi/omega`` of
the strongest oscillation frequency.

This spectral approach correctly identifies the physically meaningful
oscillation timescale even in multi-scale systems where rapid noise
oscillations sit on top of a slower mixing envelope.  Uncertainty is
estimated from the half-width at half-maximum (HWHM = FWHM/2) of the
spectral peak — this represents the distance from the peak center to
the half-power point.

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
    The mixing length ``pi/omega`` gives the half-period of energy exchange
    at that frequency.  FWHM measures how sharply defined the
    oscillation is — a narrow peak means a coherent, well-defined
    oscillation; a broad peak means the frequency is less certain.

    Attributes
    ----------
    frequency : float
        Angular frequency ``omega`` (rad/time) of the spectral peak.
    power : float
        ``|P_hat(omega)|**2`` — spectral power at this frequency.
    mixing_length : float
        ``pi/omega`` — half-period of energy exchange at this frequency.
    fwhm : float
        Full width at half maximum of the peak (rad/time).
    mixing_length_uncertainty : float
        Propagated from half-width: ``(pi/omega**2) * HWHM``
        where ``HWHM = FWHM / 2``.
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

    The uncertainty comes from the half-width at half-maximum (HWHM)
    of the spectral peak: ``dL = (pi/omega**2) * HWHM`` where
    ``HWHM = FWHM / 2``.  A sharp spectral peak (small FWHM) gives a
    precise mixing length; a broad peak indicates the oscillation
    frequency is less well-defined.

    Attributes
    ----------
    mixing_length : float
        ``pi/omega_dom`` — half-period of the dominant oscillation frequency.
    mixing_length_uncertainty : float
        Propagated from HWHM of the dominant spectral peak.
    dominant_frequency : float
        ``omega_dom`` — angular frequency of the strongest spectral peak.
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
    theory-specific conversion.  For each frequency ``omega``, the half-period
    ``pi/omega`` gives the mixing timescale at that frequency.

    The frequency resolution is ``rayleigh_resolution = 2*pi/T`` where ``T``
    is the observation window duration.

    Attributes
    ----------
    frequencies : ndarray
        Angular frequencies ``omega`` (rad/time), excluding DC.
    power : ndarray
        ``|P_hat(omega)|**2`` at each frequency.
    dominant_frequency : float
        ``omega`` of the strongest oscillation peak.
    dominant_mixing_length : float
        ``pi / dominant_frequency`` — half-period at the dominant frequency.
    rayleigh_resolution : float
        ``2*pi/T`` — fundamental frequency resolution from the observation
        window duration.  This is the minimum resolvable frequency difference.
    """

    frequencies: NDArray[np.float64]
    power: NDArray[np.float64]
    dominant_frequency: float
    dominant_mixing_length: float
    rayleigh_resolution: float


_MIN_POINTS = 3  # Minimum time points for peak detection / FFT

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _compute_fwhm(
    freqs: NDArray[np.float64],
    power: NDArray[np.float64],
    peak_idx: int,
    rayleigh_resolution: float,
) -> float:
    """Full width at half maximum of a spectral peak.

    Uses linear interpolation between bins for sub-bin accuracy.
    If one side of the peak reaches the spectrum boundary without
    crossing half-max, uses the available one-sided width doubled
    as an estimate.

    The FWHM is floored at the Rayleigh resolution (``2*pi/T``),
    which is the minimum resolvable frequency width.
    """
    half_max = float(power[peak_idx]) / 2.0
    n = len(freqs)
    peak_freq = float(freqs[peak_idx])

    # --- Scan left ---
    left_freq = float(freqs[0])  # default: left edge
    left_found = False
    for j in range(peak_idx - 1, -1, -1):
        if power[j] < half_max:
            left_found = True
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
    right_found = False
    for j in range(peak_idx + 1, n):
        if power[j] < half_max:
            right_found = True
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

    # Combine: use one-sided doubling when a scan reached the edge
    if left_found and right_found:
        fwhm = right_freq - left_freq
    elif right_found and not left_found:
        fwhm = 2.0 * (right_freq - peak_freq)
    elif left_found and not right_found:
        fwhm = 2.0 * (peak_freq - left_freq)
    else:
        # Both sides hit boundary — use full span (best effort)
        fwhm = right_freq - left_freq

    # Floor: at least the Rayleigh resolution (2*pi/T)
    return max(fwhm, rayleigh_resolution)


def _find_spectral_peaks(
    freqs: NDArray[np.float64],
    power: NDArray[np.float64],
    min_prominence: float,
    rayleigh_resolution: float,
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
        if omega <= 0.0:
            continue  # skip non-positive frequencies (avoid division by zero)
        p = float(power[idx])
        fwhm = _compute_fwhm(freqs, power, idx, rayleigh_resolution)
        ml = np.pi / omega
        # Uncertainty from HWHM (half-width at half-maximum = FWHM/2)
        hwhm = fwhm / 2.0
        ml_unc = (np.pi / (omega * omega)) * hwhm
        peaks.append(
            SpectralPeak(
                frequency=omega,
                power=p,
                mixing_length=ml,
                fwhm=fwhm,
                mixing_length_uncertainty=ml_unc,
            ),
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

    The uncertainty comes from the half-width at half-maximum (HWHM) of
    the spectral peak.  Error propagation: ``L = pi/omega``,
    ``dL = (pi/omega**2) * HWHM`` where ``HWHM = FWHM/2``.

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
    peaks = _find_spectral_peaks(
        spectrum.frequencies,
        spectrum.power,
        min_prominence,
        spectrum.rayleigh_resolution,
    )

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
    frequency ``omega``, the half-period ``pi/omega`` gives the corresponding
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

    # Rayleigh resolution: fundamental frequency resolution from observation window
    rayleigh_resolution = 2.0 * np.pi / float(times[-1] - times[0])

    # FFT (subtract mean to remove DC component)
    prob_centered = prob - np.mean(prob)
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
    if omega_dom <= 0.0:
        msg = "Dominant frequency is zero — cannot compute mixing length"
        raise ValueError(msg)

    return MixingSpectrum(
        frequencies=angular_freqs,
        power=power,
        dominant_frequency=omega_dom,
        dominant_mixing_length=np.pi / omega_dom,
        rayleigh_resolution=rayleigh_resolution,
    )
