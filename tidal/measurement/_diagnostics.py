"""Energy conservation and summary diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.measurement._energy import (
    ENERGY_FLOOR,
    compute_energy_timeseries,
)

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class EnergyDiagnostics:
    """Energy conservation diagnostic result.

    Attributes
    ----------
    times : ndarray, shape ``(n_snapshots,)``
    total_energy : ndarray, shape ``(n_snapshots,)``
        Spatially-averaged energy density ⟨ε⟩ at each snapshot.
    relative_error : ndarray, shape ``(n_snapshots,)``
        ``(⟨ε⟩(t) - ⟨ε⟩(0)) / ⟨ε⟩(0)``.
    max_relative_error : float
        Peak relative energy density drift.
    is_conserved : bool
        Whether ``max_relative_error < threshold``.
    """

    times: NDArray[np.float64]
    total_energy: NDArray[np.float64]
    relative_error: NDArray[np.float64]
    max_relative_error: float
    is_conserved: bool


def check_energy_conservation(
    data: SimulationData,
    threshold: float = 1e-3,
) -> EnergyDiagnostics:
    """Check whether energy density is conserved over the simulation.

    Parameters
    ----------
    data : SimulationData
    threshold : float
        Maximum allowed ``|ΔE/E₀|``.  Default ``1e-3`` (0.1%).

    Returns
    -------
    EnergyDiagnostics

    Raises
    ------
    ValueError
        If *threshold* is not positive.
    """
    if threshold <= 0:
        msg = f"threshold must be positive, got {threshold}"
        raise ValueError(msg)

    times, _per_field, _interaction, total = compute_energy_timeseries(data)

    e0 = total[0]
    relative_error = (total - e0) / e0 if e0 >= ENERGY_FLOOR else np.zeros_like(total)

    max_err = float(np.max(np.abs(relative_error)))

    return EnergyDiagnostics(
        times=times,
        total_energy=total,
        relative_error=relative_error,
        max_relative_error=max_err,
        is_conserved=max_err < threshold,
    )


def summarize(data: SimulationData) -> dict[str, Any]:
    """Compute a measurement summary of the simulation.

    Returns
    -------
    dict with keys:
        - ``per_field_energy``: ``dict[str, list[float]]`` time series
        - ``interaction_energy``: ``list[float]``
        - ``total_energy``: ``list[float]``
        - ``energy_conservation``: :class:`EnergyDiagnostics`
        - ``field_peaks``: ``dict[str, tuple[float, float]]`` (initial, final peak amplitude)
    """
    times, per_field, interaction, total = compute_energy_timeseries(data)

    # Peak amplitudes
    field_peaks: dict[str, tuple[float, float]] = {}
    for name in data.fields:
        initial_peak = float(np.max(np.abs(data.fields[name][0])))
        final_peak = float(np.max(np.abs(data.fields[name][-1])))
        field_peaks[name] = (initial_peak, final_peak)

    return {
        "times": times.tolist(),
        "per_field_energy": {k: v.tolist() for k, v in per_field.items()},
        "interaction_energy": interaction.tolist(),
        "total_energy": total.tolist(),
        "energy_conservation": check_energy_conservation(data),
        "field_peaks": field_peaks,
    }
