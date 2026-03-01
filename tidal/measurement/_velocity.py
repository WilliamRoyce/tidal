"""Group and phase velocity extraction from dispersion relation.

Computes:

- **Phase velocity:** ``v_p(k) = omega(k) / k``
- **Group velocity:** ``v_g(k) = d omega / dk`` (numerical derivative)
- **Velocity mismatch:** ``|v_g_source(k) - v_g_target(k)|`` for two-field systems

Requires flat + spatially homogeneous systems (inherits from dispersion).
For curved spacetimes, wave-packet tracking would be needed instead.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._dispersion import compute_dispersion
from tidal.measurement._mode_utils import find_shared_modes

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class VelocityResult:
    """Velocity analysis from dispersion relation.

    Attributes
    ----------
    wavenumbers : ndarray, shape ``(n_active,)``
        Wavenumber |k| for active modes.
    group_velocity : ndarray, shape ``(n_active,)``
        Group velocity ``d omega / dk`` per active mode.
    phase_velocity : ndarray, shape ``(n_active,)``
        Phase velocity ``omega / k`` per active mode.
    group_velocity_mean : float
        Amplitude-weighted mean group velocity.
    phase_velocity_mean : float
        Amplitude-weighted mean phase velocity.
    n_active_modes : int
        Number of active modes used.
    field_name : str
        Field group name (comma-joined).
    """

    wavenumbers: NDArray[np.float64]
    group_velocity: NDArray[np.float64]
    phase_velocity: NDArray[np.float64]
    group_velocity_mean: float
    phase_velocity_mean: float
    n_active_modes: int
    field_name: str


def _compute_group_velocity(
    wavenumbers: NDArray[np.float64],
    frequencies: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute group velocity d omega/dk via finite differences.

    Uses central differences where possible, forward/backward at edges.
    """
    n = len(wavenumbers)
    if n < 2:  # noqa: PLR2004
        return np.zeros(n, dtype=np.float64)

    vg = np.zeros(n, dtype=np.float64)

    # Central differences for interior points
    if n > 2:  # noqa: PLR2004
        dk = wavenumbers[2:] - wavenumbers[:-2]
        domega = frequencies[2:] - frequencies[:-2]
        nonzero = dk > 0
        vg[1:-1] = np.where(nonzero, domega / dk, 0.0)

    # Forward difference at left edge
    dk0 = wavenumbers[1] - wavenumbers[0]
    if dk0 > 0:
        vg[0] = (frequencies[1] - frequencies[0]) / dk0

    # Backward difference at right edge
    dkn = wavenumbers[-1] - wavenumbers[-2]
    if dkn > 0:
        vg[-1] = (frequencies[-1] - frequencies[-2]) / dkn

    return vg


def compute_velocities(
    data: SimulationData,
    field_names: str | Sequence[str],
    *,
    min_amplitude: float = 1e-12,
) -> VelocityResult:
    """Extract group and phase velocities from the dispersion relation.

    Parameters
    ----------
    data : SimulationData
        Simulation output with time-resolved field snapshots.
    field_names : str or sequence of str
        Field or group of fields to analyze.
    min_amplitude : float, optional
        Minimum Fourier amplitude for a mode to be considered active.

    Returns
    -------
    VelocityResult

    Raises
    ------
    ValueError
        If no active modes found or dispersion cannot be computed.
    """
    disp = compute_dispersion(data, field_names, min_amplitude=min_amplitude)

    active = disp.peak_frequencies > 0.0
    if not np.any(active):
        msg = f"No active modes found for {disp.field_name} — cannot extract velocities"
        raise ValueError(msg)

    k_active = disp.wavenumbers[active]
    omega_active = disp.peak_frequencies[active]
    power_active = disp.peak_powers[active]

    # Phase velocity: v_p = omega / k (skip k=0 to avoid division by zero)
    phase_vel = np.zeros_like(k_active)
    nonzero_k = k_active > 0
    phase_vel[nonzero_k] = omega_active[nonzero_k] / k_active[nonzero_k]

    # Group velocity: v_g = d omega / dk
    group_vel = _compute_group_velocity(k_active, omega_active)

    # Amplitude-weighted means (use spectral power as weight)
    weights = (
        power_active / power_active.sum()
        if power_active.sum() > 0
        else np.ones_like(power_active) / len(power_active)
    )
    vg_mean = float(np.dot(group_vel, weights))
    vp_mean = float(np.dot(phase_vel, weights))

    return VelocityResult(
        wavenumbers=k_active,
        group_velocity=group_vel,
        phase_velocity=phase_vel,
        group_velocity_mean=vg_mean,
        phase_velocity_mean=vp_mean,
        n_active_modes=int(active.sum()),
        field_name=disp.field_name,
    )


@dataclass(frozen=True)
class VelocityMismatchResult:
    """Velocity mismatch between two field groups.

    Attributes
    ----------
    source_velocity : VelocityResult
        Velocity analysis for the source field.
    target_velocity : VelocityResult
        Velocity analysis for the target field.
    mismatch : ndarray
        ``|v_g_source(k) - v_g_target(k)|`` at shared wavenumber bins.
    shared_wavenumbers : ndarray
        Wavenumbers at which both fields have active modes.
    max_mismatch : float
        Maximum velocity mismatch across all shared modes.
    mean_mismatch : float
        Mean velocity mismatch across shared modes.
    """

    source_velocity: VelocityResult
    target_velocity: VelocityResult
    mismatch: NDArray[np.float64]
    shared_wavenumbers: NDArray[np.float64]
    max_mismatch: float
    mean_mismatch: float


def compute_velocity_mismatch(
    data: SimulationData,
    source_field: str | Sequence[str],
    target_field: str | Sequence[str],
    *,
    min_amplitude: float = 1e-12,
) -> VelocityMismatchResult:
    """Compute group velocity mismatch between two field groups.

    Parameters
    ----------
    data : SimulationData
        Simulation output.
    source_field : str or sequence of str
        Source field name(s).
    target_field : str or sequence of str
        Target field name(s).
    min_amplitude : float, optional
        Minimum Fourier amplitude for active modes.

    Returns
    -------
    VelocityMismatchResult

    Raises
    ------
    ValueError
        If no shared active modes between source and target.
    """
    vel_src = compute_velocities(data, source_field, min_amplitude=min_amplitude)
    vel_tgt = compute_velocities(data, target_field, min_amplitude=min_amplitude)

    # Find shared wavenumber bins (matching by value with tolerance)
    idx_src, idx_tgt = find_shared_modes(vel_src.wavenumbers, vel_tgt.wavenumbers)

    if not idx_src:
        msg = "No shared active wavenumber modes between source and target fields"
        raise ValueError(msg)

    shared_k_arr = vel_src.wavenumbers[idx_src]
    mismatch = np.abs(vel_src.group_velocity[idx_src] - vel_tgt.group_velocity[idx_tgt])

    return VelocityMismatchResult(
        source_velocity=vel_src,
        target_velocity=vel_tgt,
        mismatch=mismatch,
        shared_wavenumbers=shared_k_arr,
        max_mismatch=float(mismatch.max()),
        mean_mismatch=float(mismatch.mean()),
    )
