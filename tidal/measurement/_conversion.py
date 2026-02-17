"""Conversion probability measurement for coupled field systems.

The primary measurement for Gertsenshtein-type physics: how much energy
transfers from a source field to a target field over time.

    P(t) = E_target(t) / E_source(0)
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


@dataclass(frozen=True)
class ConversionResult:
    """Result of a conversion probability measurement.

    Attributes
    ----------
    times : ndarray, shape ``(n_snapshots,)``
        Snapshot times.
    probability : ndarray, shape ``(n_snapshots,)``
        ``P(t) = E_target(t) / E_source(0)``.
    source_energy : ndarray, shape ``(n_snapshots,)``
        Source field energy over time.
    target_energy : ndarray, shape ``(n_snapshots,)``
        Target field energy over time.
    total_energy : ndarray, shape ``(n_snapshots,)``
        Total system energy (source + target) over time.
    relative_energy_error : ndarray, shape ``(n_snapshots,)``
        ``(E_total(t) - E_total(0)) / E_total(0)``.
    source_field : str
        Name of the source field.
    target_field : str
        Name of the target field.
    """

    times: NDArray[np.float64]
    probability: NDArray[np.float64]
    source_energy: NDArray[np.float64]
    target_energy: NDArray[np.float64]
    total_energy: NDArray[np.float64]
    relative_energy_error: NDArray[np.float64]
    source_field: str
    target_field: str


def _field_energy_series(
    data: SimulationData,
    field_name: str,
    mass_squared: float | NDArray[np.float64],
) -> NDArray[np.float64]:
    """Compute total energy of *field_name* at every snapshot."""
    energies: list[float] = []
    for t_idx in range(data.n_snapshots):
        field_arr = data.fields[field_name][t_idx]
        mom_all = data.momenta.get(field_name)
        mom_arr = mom_all[t_idx] if mom_all is not None else None
        fe = compute_field_energy(
            field_arr, mom_arr, mass_squared, data.grid_spacing, data.periodic,
        )
        energies.append(fe.total)
    return np.array(energies, dtype=np.float64)


def compute_conversion_probability(
    data: SimulationData,
    source_field: str,
    target_field: str,
) -> ConversionResult:
    """Compute wave conversion probability ``P(t) = E_target(t) / E_source(0)``.

    This is the primary measurement for the Gertsenshtein effect.  The source
    field is excited with some initial energy; the target field starts at zero.
    Coupling terms transfer energy between them over time.

    Parameters
    ----------
    data : SimulationData
        Full simulation output.
    source_field : str
        Name of the source field (e.g. ``"phi_0"``).
    target_field : str
        Name of the target field (e.g. ``"chi_0"``).

    Returns
    -------
    ConversionResult

    Raises
    ------
    ValueError
        If *source_field* or *target_field* is not a valid field name,
        if they are the same field, or if the source has zero initial energy.
    """
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

    source_arr = _field_energy_series(
        data, source_field, _resolve_mass_squared(data, names.index(source_field)),
    )
    target_arr = _field_energy_series(
        data, target_field, _resolve_mass_squared(data, names.index(target_field)),
    )
    total_arr = source_arr + target_arr

    e_source_0 = source_arr[0]
    if e_source_0 < _ENERGY_FLOOR:
        msg = (
            f"Source field '{source_field}' has zero initial energy — "
            f"cannot compute conversion probability"
        )
        raise ValueError(msg)

    probability = target_arr / e_source_0
    e_total_0 = total_arr[0]
    relative_error = (total_arr - e_total_0) / e_total_0 if e_total_0 >= _ENERGY_FLOOR else np.zeros_like(total_arr)

    return ConversionResult(
        times=data.times.copy(),
        probability=probability,
        source_energy=source_arr,
        target_energy=target_arr,
        total_energy=total_arr,
        relative_energy_error=relative_error,
        source_field=source_field,
        target_field=target_field,
    )


def _group_energy_series(
    data: SimulationData,
    field_group: tuple[str, ...],
) -> NDArray[np.float64]:
    """Sum per-field energy timeseries across a group of fields."""
    names = data.spec.component_names
    total: NDArray[np.float64] = np.zeros(data.n_snapshots, dtype=np.float64)
    for f in field_group:
        total += _field_energy_series(data, f, _resolve_mass_squared(data, names.index(f)))
    return total


def compute_group_conversion(
    data: SimulationData,
    source_fields: str | Sequence[str],
    target_fields: str | Sequence[str] | None = None,
) -> ConversionResult:
    """Measure energy conversion between field groups.

    Computes ``P(t) = E_targets(t) / E_sources(0)`` where each energy
    is the sum of per-field canonical energies across the group.  This is
    the natural measurement for the Gertsenshtein effect where source and
    target are multi-component tensor/vector field groups.

    Parameters
    ----------
    data : SimulationData
        Full simulation output.
    source_fields : str or sequence of str
        Source field name(s).  A single string is treated as a one-element
        group.
    target_fields : str, sequence of str, or None
        Target field name(s).  ``None`` means all other dynamical fields
        (``time_derivative_order >= 2``) not in *source_fields*.

    Returns
    -------
    ConversionResult
        ``source_field`` / ``target_field`` are comma-joined group names.

    Raises
    ------
    ValueError
        If any field name is invalid, if source and target overlap, if
        the target group is empty, or if the source group has zero
        initial energy.
    """
    names = data.spec.component_names

    # Normalize to tuples
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

    # Validate no overlap and non-empty target
    overlap = set(src) & set(tgt)
    if overlap:
        msg = f"Source and target groups overlap on: {sorted(overlap)}"
        raise ValueError(msg)
    if not tgt:
        msg = "Target group is empty — no dynamical fields remain after excluding source"
        raise ValueError(msg)

    # Sum energies across group members
    source_arr = _group_energy_series(data, src)
    target_arr = _group_energy_series(data, tgt)
    total_arr = source_arr + target_arr

    e_source_0 = float(source_arr[0])
    if e_source_0 < _ENERGY_FLOOR:
        msg = (
            f"Source group {src} has zero initial energy — "
            f"cannot compute conversion probability"
        )
        raise ValueError(msg)

    probability = target_arr / e_source_0
    e_total_0 = float(total_arr[0])
    relative_error: NDArray[np.float64] = (
        (total_arr - e_total_0) / e_total_0
        if e_total_0 >= _ENERGY_FLOOR
        else np.zeros_like(total_arr)
    )

    return ConversionResult(
        times=data.times.copy(),
        probability=probability,
        source_energy=source_arr,
        target_energy=target_arr,
        total_energy=total_arr,
        relative_energy_error=relative_error,
        source_field=",".join(src),
        target_field=",".join(tgt),
    )
