"""Asymptotic conversion measurement — frame-independent scattering observables.

Computes:

- **P_final**: ``E_target(t_final) / E_source(t=0)`` — total conversion at end
  of simulation.  Uses field-group energies (summed over all components),
  making this a Lorentz scalar for any field type.

- **P_forward / P_reflected**: Directional decomposition of the target-field
  energy at the final snapshot.  The forward/backward split is defined by
  the **source field's initial propagation direction** — the spectral centroid
  wavevector ``⟨k⟩`` of the source at ``t=0``.

  - Forward (transmitted): target modes with ``k · k̂_source > 0``
  - Reflected: target modes with ``k · k̂_source < 0``

  This definition is:
  - **Coordinate-independent**: the reference direction comes from the physics
    (the initial wave), not from the choice of spatial axes.
  - **Dimension-independent**: uses the full wavevector dot product, working
    identically in 1+1D, 2+1D, and 3+1D.
  - **Rotation-invariant**: ``k · k̂`` is a scalar under spatial rotations.

  It is NOT Lorentz-invariant (wavevectors transform under boosts), but for
  simulations on a fixed spatial grid, the grid frame is the natural frame
  and the definition is unambiguous.

For scalar fields, φ² is a Lorentz scalar so P_final is automatically
frame-independent.  For vector fields, the group sum over all components
preserves rotational covariance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._energy import (
    _resolve_mass_squared,  # pyright: ignore[reportPrivateUsage]
    compute_field_energy,
)
from tidal.measurement._utils import (
    _normalize_group,  # pyright: ignore[reportPrivateUsage]
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData

# Floor for energy comparisons (avoid division by zero)
_ENERGY_FLOOR = 1e-30
# Floor for wavevector magnitude (detect standing waves)
_WAVEVECTOR_FLOOR = 1e-12


@dataclass(frozen=True)
class AsymptoticConversionResult:
    """Result of asymptotic conversion measurement.

    Attributes
    ----------
    P_final : float
        Total conversion at final snapshot: ``E_target(t_final) / E_source(0)``.
    P_reflected : float
        Fraction of total conversion in backward-propagating target modes.
    P_forward : float
        Fraction of total conversion in forward-propagating target modes.
    E_source_initial : float
        Source group energy at ``t=0``.
    E_target_final : float
        Target group energy at final snapshot.
    source_field : str
        Source field group name.
    target_field : str
        Target field group name.
    source_wavevector : tuple[float, ...]
        Spectral centroid wavevector of the source at ``t=0`` (defines the
        forward direction).
    """

    P_final: float
    P_reflected: float
    P_forward: float
    E_source_initial: float
    E_target_final: float
    source_field: str
    target_field: str
    source_wavevector: tuple[float, ...]


def _group_energy_at_snapshot(
    data: SimulationData,
    field_names: list[str],
    t_idx: int,
) -> float:
    """Compute total energy for a group of fields at snapshot *t_idx*."""
    total = 0.0
    for fname in field_names:
        field_arr = data.fields[fname][t_idx]
        vel_all = data.velocities.get(fname)
        vel_arr = vel_all[t_idx] if vel_all is not None else None
        m2 = _resolve_mass_squared(data, fname)
        fe = compute_field_energy(
            field_arr,
            vel_arr,
            m2,
            data.grid_spacing,
            data.periodic,
            volume_weight=getattr(data, "volume_weight", 1.0),
        )
        total += fe.total
    return total


def _build_k_grids(
    field_shape: tuple[int, ...],
    grid_spacing: tuple[float, ...],
) -> list[NDArray[np.float64]]:
    """Build per-axis wavenumber grids (full FFT, not rFFT)."""
    ndim = len(grid_spacing)
    k_arrays: list[NDArray[np.float64]] = []
    for ax in range(ndim):
        n = field_shape[ax]
        k_ax = np.fft.fftfreq(n, d=grid_spacing[ax]) * (2.0 * np.pi)
        k_arrays.append(k_ax)
    return [
        np.asarray(g, dtype=np.float64)
        for g in np.meshgrid(*k_arrays, indexing="ij")
    ]


def _source_wavevector(
    data: SimulationData,
    source_fields: list[str],
) -> NDArray[np.float64]:
    """Compute the spectral centroid wavevector of the source at t=0.

    The spectral centroid is the energy-weighted average wavevector:

        ⟨k⟩ = Σ_k (k · |φ̂(k)|²) / Σ_k |φ̂(k)|²

    summed over all source fields (group-covariant).  This defines the
    propagation direction of the initial wave packet.

    For a plane wave with wavevector k₀, this returns k₀ exactly.
    For a wave packet, it returns the central wavevector.
    """
    ndim = len(data.grid_spacing)
    k_centroid = np.zeros(ndim)
    total_power = 0.0

    for fname in source_fields:
        field_0 = np.asarray(data.fields[fname][0], dtype=np.float64)
        fhat = np.fft.fftn(field_0)
        power = np.abs(fhat) ** 2

        k_grids = _build_k_grids(field_0.shape, data.grid_spacing)
        for ax in range(ndim):
            k_centroid[ax] += float(np.sum(k_grids[ax] * power))
        total_power += float(np.sum(power))

    if total_power > 0:
        k_centroid /= total_power

    return k_centroid


def _field_spectral_energy(
    data: SimulationData,
    fname: str,
    t_idx: int,
) -> NDArray[np.float64]:
    """Spectral energy density for one field: |φ̂(k)|² + |v̂(k)|²."""
    field_arr = np.asarray(data.fields[fname][t_idx], dtype=np.float64)
    power = np.abs(np.fft.fftn(field_arr)) ** 2

    vel_all = data.velocities.get(fname)
    if vel_all is not None:
        vel_arr = np.asarray(vel_all[t_idx], dtype=np.float64)
        power += np.abs(np.fft.fftn(vel_arr)) ** 2

    return power


def _directional_split(
    data: SimulationData,
    target_fields: list[str],
    t_idx: int,
    k_hat: NDArray[np.float64],
) -> tuple[float, float]:
    """Split target field spectral energy into forward and reflected fractions.

    Forward = modes where ``k · k̂_source > 0``.
    Reflected = modes where ``k · k̂_source < 0``.

    Parameters
    ----------
    data : SimulationData
        Simulation output.
    target_fields : list[str]
        Target field names.
    t_idx : int
        Snapshot index.
    k_hat : ndarray, shape ``(ndim,)``
        Unit vector in the source propagation direction.

    Returns
    -------
    forward_frac, reflected_frac : float
        Fractions of total spectral energy (each in [0, 1], sum ≤ 1).
    """
    field_shape = data.fields[target_fields[0]][t_idx].shape
    k_grids = _build_k_grids(field_shape, data.grid_spacing)

    # k_dot = k · k̂_source — scalar field on the k-grid
    k_dot = sum(k_grids[ax] * k_hat[ax] for ax in range(len(data.grid_spacing)))

    forward_mask = k_dot > 0
    reflected_mask = k_dot < 0

    total_forward = 0.0
    total_reflected = 0.0

    for fname in target_fields:
        spectral_energy = _field_spectral_energy(data, fname, t_idx)
        total_forward += float(np.sum(spectral_energy[forward_mask]))
        total_reflected += float(np.sum(spectral_energy[reflected_mask]))

    total = total_forward + total_reflected
    if total < _ENERGY_FLOOR:
        return 0.0, 0.0
    return total_forward / total, total_reflected / total


def compute_asymptotic_conversion(
    data: SimulationData,
    source_fields: str | Sequence[str],
    target_fields: str | Sequence[str] | None = None,
) -> AsymptoticConversionResult:
    """Compute asymptotic scattering observables.

    The forward/backward directional split is defined by the source
    field's initial propagation direction (spectral centroid wavevector
    at ``t=0``).  See module docstring for details on frame independence.

    Parameters
    ----------
    data : SimulationData
        Simulation output.
    source_fields : str or sequence of str
        Source field name(s).
    target_fields : str, sequence of str, or None
        Target field name(s).  ``None`` → all dynamical fields not in source.

    Returns
    -------
    AsymptoticConversionResult

    Raises
    ------
    ValueError
        If source/target fields are invalid or overlap, or if the source
        has zero initial energy.
    """
    source_list = _normalize_group(source_fields)
    if target_fields is None:
        all_dyn = list(data.dynamical_fields)
        target_list = [f for f in all_dyn if f not in source_list]
        if not target_list:
            msg = "No target fields: all dynamical fields are in the source set"
            raise ValueError(msg)
    else:
        target_list = _normalize_group(target_fields)

    overlap = set(source_list) & set(target_list)
    if overlap:
        msg = f"Source and target overlap: {overlap}"
        raise ValueError(msg)

    # Energies
    e_source_0 = _group_energy_at_snapshot(data, source_list, 0)
    e_target_final = _group_energy_at_snapshot(data, target_list, -1)

    if e_source_0 < _ENERGY_FLOOR:
        msg = "Source initial energy is zero — cannot compute conversion"
        raise ValueError(msg)

    p_final = e_target_final / e_source_0

    # Source propagation direction (spectral centroid at t=0)
    k_source = _source_wavevector(data, source_list)
    k_mag = float(np.linalg.norm(k_source))

    if k_mag < _WAVEVECTOR_FLOOR:
        # Source has no net propagation direction (e.g., standing wave).
        # Directional decomposition is undefined — report equal split.
        p_forward = p_final / 2.0
        p_reflected = p_final / 2.0
    else:
        k_hat = k_source / k_mag
        fwd_frac, ref_frac = _directional_split(data, target_list, -1, k_hat)
        p_forward = p_final * fwd_frac
        p_reflected = p_final * ref_frac

    return AsymptoticConversionResult(
        P_final=p_final,
        P_reflected=p_reflected,
        P_forward=p_forward,
        E_source_initial=e_source_0,
        E_target_final=e_target_final,
        source_field=",".join(source_list),
        target_field=",".join(target_list),
        source_wavevector=tuple(float(x) for x in k_source),
    )
