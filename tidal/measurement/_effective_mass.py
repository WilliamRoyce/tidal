"""Effective mass extraction from the dispersion relation.

Computes the Lorentz-invariant effective mass squared:

    m²_eff = ω² - k²

from the dominant frequency ω(k) at each active wavenumber k.  This is
the 4-momentum norm and is a **Lorentz scalar** — frame-independent by
construction.

The effective mass is extracted from the dispersion relation computed by
:func:`compute_dispersion`.  For a free Klein-Gordon field, m²_eff
should match the parameter mass.  Deviations indicate interaction or
medium effects (mass renormalization from coupling).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._dispersion import compute_dispersion

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class EffectiveMassResult:
    """Result of effective mass extraction.

    Attributes
    ----------
    m2_eff : float
        Median effective mass squared across active modes.
    m2_eff_std : float
        Standard deviation of m² across active modes (spread indicator).
    m2_eff_values : ndarray
        Per-mode m² = ω² - k² for all active modes.
    wavenumbers : ndarray
        Wavenumber |k| for each active mode.
    frequencies : ndarray
        Dominant ω(k) for each active mode.
    n_active_modes : int
        Number of active k-modes used in the estimate.
    field_name : str
        Field group name (comma-joined).
    """

    m2_eff: float
    m2_eff_std: float
    m2_eff_values: NDArray[np.float64]
    wavenumbers: NDArray[np.float64]
    frequencies: NDArray[np.float64]
    n_active_modes: int
    field_name: str


def compute_effective_mass(
    data: SimulationData,
    field_names: str | Sequence[str],
    *,
    min_amplitude: float = 1e-12,
) -> EffectiveMassResult:
    """Extract effective mass from the dispersion relation.

    Parameters
    ----------
    data : SimulationData
        Simulation output.
    field_names : str or sequence of str
        Field name(s) to analyze.  Multiple fields are summed
        (rotationally covariant within the group).
    min_amplitude : float
        Minimum Fourier amplitude for a mode to be considered active.

    Returns
    -------
    EffectiveMassResult
        Effective mass and per-mode breakdown.

    Raises
    ------
    ValueError
        If no active modes are found or dispersion cannot be computed.
    """
    disp = compute_dispersion(data, field_names, min_amplitude=min_amplitude)

    # Extract active modes
    active = disp.peak_frequencies > 0.0
    if not np.any(active):
        msg = f"No active modes found for {disp.field_name} — cannot extract effective mass"
        raise ValueError(msg)

    k_active = disp.wavenumbers[active]
    omega_active = disp.peak_frequencies[active]
    m2_values = omega_active**2 - k_active**2

    return EffectiveMassResult(
        m2_eff=float(np.median(m2_values)),
        m2_eff_std=float(np.std(m2_values)),
        m2_eff_values=m2_values,
        wavenumbers=k_active,
        frequencies=omega_active,
        n_active_modes=int(active.sum()),
        field_name=disp.field_name,
    )
