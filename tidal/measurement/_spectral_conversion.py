"""Per-mode spectral conversion P(k,t) for coupled field systems.

Decomposes the scalar conversion probability P(t) = E_target(t) / E_source(0)
into a 2D array P(k,t) = E_target(k,t) / E_source(k, t=0), revealing which
Fourier modes participate in energy conversion.

For systems with derivative coupling (gradient/cross-derivative cross-field
terms), the mixing angle theta(k) is k-dependent, so different modes have
different oscillation amplitudes --- not just different frequencies.

    P(k,t) = E_target(k,t) / E_source(k, t=0)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from tidal.measurement._energy import (
    _ENERGY_FLOOR,  # pyright: ignore[reportPrivateUsage]
    _resolve_mass_squared,  # pyright: ignore[reportPrivateUsage]
)
from tidal.measurement._spectral import compute_spectral_energy
from tidal.measurement._utils import (
    _check_no_position_dependent_terms,  # pyright: ignore[reportPrivateUsage]
    _normalize_group,  # pyright: ignore[reportPrivateUsage]
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData


@dataclass(frozen=True)
class SpectralConversion:
    """Per-mode spectral conversion probability P(k,t).

    Tracks which Fourier modes participate in energy conversion between
    source and target fields.  For each wavenumber bin *k*,
    ``P(k,t) = E_target(k,t) / E_source(k, t=0)``.

    Attributes
    ----------
    times : ndarray, shape ``(n_snapshots,)``
        Snapshot times.
    wavenumbers : ndarray, shape ``(n_modes,)``
        Radially binned wavenumber magnitudes ``|k|``.
    probability : ndarray, shape ``(n_snapshots, n_modes)``
        Per-mode conversion probability.  Zero for inactive modes.
    source_spectral_energy : ndarray, shape ``(n_snapshots, n_modes)``
        Source spectral energy ``E_source(k,t)`` at every snapshot.
    target_spectral_energy : ndarray, shape ``(n_snapshots, n_modes)``
        Target spectral energy ``E_target(k,t)`` at every snapshot.
    active_modes : ndarray, shape ``(n_modes,)``, dtype bool
        ``True`` for k-bins where ``E_source(k, t=0) >= energy_floor``.
        Inactive bins have ``probability[:, k] = 0``.
    source_field : str
        Name of the source field (or comma-joined group).
    target_field : str
        Name of the target field (or comma-joined group).
    """

    times: NDArray[np.float64]
    wavenumbers: NDArray[np.float64]
    probability: NDArray[np.float64]
    source_spectral_energy: NDArray[np.float64]
    target_spectral_energy: NDArray[np.float64]
    active_modes: NDArray[np.bool_]
    source_field: str
    target_field: str


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _field_spectral_energy_series(
    data: SimulationData,
    field_name: str,
    mass_squared: float | NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Compute per-mode spectral energy E(k,t) for one field at all snapshots.

    Returns
    -------
    wavenumbers : ndarray, shape ``(n_modes,)``
        Radially binned ``|k|`` values.
    spectral_energy : ndarray, shape ``(n_snapshots, n_modes)``
        Energy per wavenumber bin at each snapshot.
    """
    # First snapshot — establishes the k-binning
    mom_all = data.momenta.get(field_name)
    mom_0 = mom_all[0] if mom_all is not None else None
    wavenumbers, e0 = compute_spectral_energy(
        data.fields[field_name][0],
        mom_0,
        mass_squared,
        data.grid_spacing,
        data.periodic,
    )
    n_modes = len(wavenumbers)
    energy = np.zeros((data.n_snapshots, n_modes), dtype=np.float64)
    energy[0] = e0

    for t_idx in range(1, data.n_snapshots):
        mom_t = mom_all[t_idx] if mom_all is not None else None
        _wn, e_t = compute_spectral_energy(
            data.fields[field_name][t_idx],
            mom_t,
            mass_squared,
            data.grid_spacing,
            data.periodic,
        )
        energy[t_idx] = e_t

    return wavenumbers, energy


# ------------------------------------------------------------------
# Public API
# ------------------------------------------------------------------


def compute_spectral_conversion(
    data: SimulationData,
    source_field: str,
    target_field: str,
    *,
    energy_floor: float = _ENERGY_FLOOR,
) -> SpectralConversion:
    """Compute per-mode spectral conversion ``P(k,t)``.

    ``P(k,t) = E_target(k,t) / E_source(k, t=0)`` for each wavenumber
    bin, revealing which Fourier modes participate in conversion.

    Parameters
    ----------
    data : SimulationData
        Full simulation output.
    source_field : str
        Name of the source field (e.g. ``"phi_0"``).
    target_field : str
        Name of the target field (e.g. ``"chi_0"``).
    energy_floor : float
        Per-mode energy floor.  k-bins where ``E_source(k,0) < floor``
        are masked (``active_modes=False``, ``P=0``).

    Returns
    -------
    SpectralConversion

    Raises
    ------
    ValueError
        If field names are invalid, the same, no modes are above the
        energy floor, or any equation term is position-dependent
        (translation invariance required for Fourier decomposition).
    """
    _check_no_position_dependent_terms(data, "Spectral conversion P(k,t)")

    names = data.spec.component_names
    if source_field not in names:
        msg = f"Source field '{source_field}' not in spec fields: {names}"
        raise ValueError(msg)
    if target_field not in names:
        msg = f"Target field '{target_field}' not in spec fields: {names}"
        raise ValueError(msg)
    if source_field == target_field:
        msg = f"Source and target must be different fields, got '{source_field}'"
        raise ValueError(msg)

    src_m2 = _resolve_mass_squared(data, names.index(source_field))
    tgt_m2 = _resolve_mass_squared(data, names.index(target_field))

    wavenumbers, src_energy = _field_spectral_energy_series(data, source_field, src_m2)
    _wn_tgt, tgt_energy = _field_spectral_energy_series(data, target_field, tgt_m2)

    # Determine active modes from initial source energy
    active = src_energy[0] >= energy_floor
    if not np.any(active):
        msg = (
            f"Source field '{source_field}' has no modes above energy floor "
            f"({energy_floor:.1e}) — cannot compute spectral conversion"
        )
        raise ValueError(msg)

    # P(k,t) = E_target(k,t) / E_source(k, t=0) for active modes
    probability = np.zeros_like(tgt_energy)
    denom = np.maximum(src_energy[0, active], energy_floor)
    probability[:, active] = tgt_energy[:, active] / denom

    return SpectralConversion(
        times=data.times.copy(),
        wavenumbers=wavenumbers,
        probability=probability,
        source_spectral_energy=src_energy,
        target_spectral_energy=tgt_energy,
        active_modes=active,
        source_field=source_field,
        target_field=target_field,
    )


def compute_group_spectral_conversion(
    data: SimulationData,
    source_fields: str | Sequence[str],
    target_fields: str | Sequence[str] | None = None,
    *,
    energy_floor: float = _ENERGY_FLOOR,
) -> SpectralConversion:
    """Compute per-mode spectral conversion for field groups.

    Sums spectral energies across group members before computing
    ``P(k,t) = E_targets(k,t) / E_sources(k, t=0)``.

    Parameters
    ----------
    data : SimulationData
        Full simulation output.
    source_fields : str or sequence of str
        Source field name(s).
    target_fields : str, sequence of str, or None
        Target field name(s).  ``None`` means all other dynamical fields.
    energy_floor : float
        Per-mode energy floor.

    Returns
    -------
    SpectralConversion
        ``source_field`` / ``target_field`` are comma-joined group names.

    Raises
    ------
    ValueError
        If any field name is invalid, groups overlap, target is empty,
        no modes are above the energy floor, or any equation term is
        position-dependent (translation invariance required).
    """
    _check_no_position_dependent_terms(data, "Spectral conversion P(k,t)")

    names = data.spec.component_names

    src = _normalize_group(source_fields)
    if target_fields is None:
        src_set = set(src)
        tgt = tuple(f for f in data.dynamical_fields if f not in src_set)
    else:
        tgt = _normalize_group(target_fields)

    # Validate field names
    for name in (*src, *tgt):
        if name not in names:
            msg = f"Field '{name}' not in spec fields: {names}"
            raise ValueError(msg)

    overlap = set(src) & set(tgt)
    if overlap:
        msg = f"Source and target groups overlap on: {sorted(overlap)}"
        raise ValueError(msg)
    if not tgt:
        msg = (
            "Target group is empty — no dynamical fields remain after excluding source"
        )
        raise ValueError(msg)

    # Sum spectral energies across group members
    def _group_spectral(
        group: tuple[str, ...],
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        wn, total = _field_spectral_energy_series(
            data,
            group[0],
            _resolve_mass_squared(data, names.index(group[0])),
        )
        for f in group[1:]:
            _wn, e = _field_spectral_energy_series(
                data,
                f,
                _resolve_mass_squared(data, names.index(f)),
            )
            total += e
        return wn, total

    wavenumbers, src_energy = _group_spectral(src)
    _wn_tgt, tgt_energy = _group_spectral(tgt)

    active = src_energy[0] >= energy_floor
    if not np.any(active):
        msg = (
            f"Source group {src} has no modes above energy floor "
            f"({energy_floor:.1e}) — cannot compute spectral conversion"
        )
        raise ValueError(msg)

    probability = np.zeros_like(tgt_energy)
    probability[:, active] = tgt_energy[:, active] / src_energy[0, active]

    return SpectralConversion(
        times=data.times.copy(),
        wavenumbers=wavenumbers,
        probability=probability,
        source_spectral_energy=src_energy,
        target_spectral_energy=tgt_energy,
        active_modes=active,
        source_field=",".join(src),
        target_field=",".join(tgt),
    )
