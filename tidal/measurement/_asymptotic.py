"""Asymptotic conversion measurement — frame-independent scattering observables.

Computes:

- **P_final**: ``E_target(t_final) / E_source(t=0)`` — total conversion at end
  of simulation.  Uses field-group energies (summed over all components),
  making this a Lorentz scalar for any field type.

- **P_transmitted / P_reflected**: Directional decomposition of the
  target-field energy using **flux-based mode classification**.  For each
  Fourier mode k, the directional flux
  ``(k · k̂_source) · Im[v̂*(k) · φ̂(k)]`` determines propagation direction.
  The mode's spectral energy is attributed accordingly.

  This correctly handles real-valued fields where |φ̂(k)|² = |φ̂(-k)|²:
  the cross-spectrum ``Im[v̂* · φ̂]`` is antisymmetric, so for a rightward
  wave both the +k₀ and -k₀ modes have positive directional flux and are
  both attributed to forward.

  Properties:
  - **Coordinate-independent**: reference direction from source physics.
  - **Dimension-independent**: full wavevector dot product, works in 1-3D.
  - **Rotation-invariant**: ``k · k̂`` is a scalar under spatial rotations.
  - **Correctly handles real fields**: flux breaks the power-spectrum symmetry.

For scalar fields, φ² is a Lorentz scalar so P_final is automatically
frame-independent.  For vector fields, the group sum over all components
preserves rotational covariance.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._energy import (
    _ENERGY_FLOOR,  # pyright: ignore[reportPrivateUsage]
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
    P_transmitted : float
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
    P_transmitted: float
    E_source_initial: float
    E_target_final: float
    source_field: str
    target_field: str
    source_wavevector: tuple[float, ...]


def _group_energy_at_snapshot(
    data: SimulationData,
    field_names: Sequence[str],
    t_idx: int,
) -> float:
    """Compute total energy for a group of fields at snapshot *t_idx*."""
    names = [eq.field_name for eq in data.spec.equations]
    total = 0.0
    for fname in field_names:
        field_arr = data.fields[fname][t_idx]
        vel_all = data.velocities.get(fname)
        vel_arr = vel_all[t_idx] if vel_all is not None else None
        m2 = _resolve_mass_squared(data, names.index(fname))
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
        k_ax: NDArray[np.float64] = np.asarray(
            np.fft.fftfreq(n, d=grid_spacing[ax]) * (2.0 * np.pi), dtype=np.float64
        )
        k_arrays.append(k_ax)
    return [
        np.asarray(g, dtype=np.float64) for g in np.meshgrid(*k_arrays, indexing="ij")
    ]


def _source_wavevector(
    data: SimulationData,
    source_fields: Sequence[str],
) -> NDArray[np.float64]:
    """Compute the propagation direction of the source at t=0.

    Uses the **energy flux** (momentum density) in Fourier space to determine
    the net propagation direction.  For real fields, the power spectrum
    ``|φ̂(k)|²`` is symmetric (``|φ̂(k)| = |φ̂(-k)|``), so a power-weighted
    centroid is always zero.  The energy flux breaks this symmetry:

        ⟨k⟩_flux = Σ_k k · Im[v̂*(k) · φ̂(k)]

    This is proportional to the momentum density ``T^{0i} = -v · ∂_i φ`` and
    correctly identifies the propagation direction for travelling waves while
    returning zero for standing waves.

    Summed over all source fields for group covariance.

    For a plane wave with wavevector k₀, this returns a vector parallel to k₀.
    The magnitude is proportional to energy flux, not |k₀| itself.
    """
    ndim = len(data.grid_spacing)
    k_flux = np.zeros(ndim)
    total_flux = 0.0

    for fname in source_fields:
        field_0 = np.asarray(data.fields[fname][0], dtype=np.float64)
        fhat = np.fft.fftn(field_0)

        vel_all = data.velocities.get(fname)
        if vel_all is None:
            continue
        vel_0 = np.asarray(vel_all[0], dtype=np.float64)
        vhat = np.fft.fftn(vel_0)

        # Cross-spectrum: Im[v̂*(k) · φ̂(k)] gives directional weight
        cross = np.imag(np.conj(vhat) * fhat)

        k_grids = _build_k_grids(field_0.shape, data.grid_spacing)
        for ax in range(ndim):
            flux_ax = float(np.sum(k_grids[ax] * cross))
            k_flux[ax] += flux_ax
        total_flux += float(np.sum(np.abs(cross)))

    # Normalize to get a unit-like direction vector scaled by net flux
    if total_flux > 0:
        k_flux /= total_flux

    return k_flux


def _directional_split(  # noqa: PLR0914
    data: SimulationData,
    target_fields: Sequence[str],
    t_idx: int,
    k_hat: NDArray[np.float64],
) -> tuple[float, float]:
    """Split target field energy into forward (transmitted) and reflected.

    Uses **flux-based mode classification**: for each Fourier mode k, the
    directional energy flux ``(k · k̂_source) · Im[v̂*(k) · φ̂(k)]``
    determines whether the mode propagates forward or backward.  The mode's
    spectral energy ``|φ̂(k)|² + |v̂(k)|²`` is then attributed to the
    corresponding direction.

    For real-valued fields, a naive spectral power split by sign of k always
    gives ~50/50 because |φ̂(k)|² = |φ̂(-k)|².  The flux-based approach
    correctly handles real fields because the directional flux uses the
    cross-spectrum ``Im[v̂* · φ̂]`` which is antisymmetric: for a rightward
    travelling wave, both the +k and -k modes have positive directional
    flux and are both attributed to the forward direction.

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
        Fractions of total spectral energy (each in [0, 1], sum ≈ 1).
    """
    total_forward = 0.0
    total_reflected = 0.0
    total_zero = 0.0

    # Build k-grids once (invariant across fields — depends only on grid geometry)
    first_field = np.asarray(data.fields[target_fields[0]][t_idx], dtype=np.float64)
    k_grids = _build_k_grids(first_field.shape, data.grid_spacing)

    for fname in target_fields:
        field_arr = np.asarray(data.fields[fname][t_idx], dtype=np.float64)
        fhat = np.fft.fftn(field_arr)
        spectral_energy = np.abs(fhat) ** 2

        vel_all = data.velocities.get(fname)
        if vel_all is not None:
            vel_arr = np.asarray(vel_all[t_idx], dtype=np.float64)
            vhat = np.fft.fftn(vel_arr)
            spectral_energy += np.abs(vhat) ** 2

            # Directional flux: (k · k̂_source) · Im[v̂*(k) · φ̂(k)]
            k_dot_source = sum(
                k_grids[ax] * k_hat[ax] for ax in range(len(data.grid_spacing))
            )
            cross = np.imag(np.conj(vhat) * fhat)
            flux = k_dot_source * cross

            flux_scale = np.max(np.abs(flux))
            eps = 1e-12 * flux_scale if flux_scale > 0 else 0.0
            fwd_mask = flux > eps
            ref_mask = flux < -eps
            zero_mask = ~fwd_mask & ~ref_mask

            total_forward += float(np.sum(spectral_energy[fwd_mask]))
            total_reflected += float(np.sum(spectral_energy[ref_mask]))
            total_zero += float(np.sum(spectral_energy[zero_mask]))
        else:
            # No velocity → cannot determine direction; split equally
            total_zero += float(np.sum(spectral_energy))

    # Split zero-flux modes equally between forward and reflected
    total_forward += total_zero / 2.0
    total_reflected += total_zero / 2.0

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
        p_transmitted = p_final / 2.0
        p_reflected = p_final / 2.0
    else:
        k_hat = k_source / k_mag
        fwd_frac, ref_frac = _directional_split(data, target_list, -1, k_hat)
        p_transmitted = p_final * fwd_frac
        p_reflected = p_final * ref_frac

    return AsymptoticConversionResult(
        P_final=p_final,
        P_reflected=p_reflected,
        P_transmitted=p_transmitted,
        E_source_initial=e_source_0,
        E_target_final=e_target_final,
        source_field=",".join(source_list),
        target_field=",".join(target_list),
        source_wavevector=tuple(float(x) for x in k_source),
    )
