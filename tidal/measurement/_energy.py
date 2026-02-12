"""Per-field canonical energy and interaction energy computation.

Computes the Hamiltonian energy density for each field in a multi-field
coupled system:

    E_i = ∫ dV [ ½ π_i² + ½ |∇φ_i|² + ½ m_ii² φ_i² ]

and the interaction energy between fields:

    E_int = ½ Σ_{i≠j} C_ij ∫ φ_i φ_j dV
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData

# Threshold below which energy is treated as zero.  Prevents division by
# near-zero values from floating-point integration noise.
_ENERGY_FLOOR: float = 1e-12


@dataclass(frozen=True)
class FieldEnergy:
    """Energy decomposition for a single field at one snapshot.

    Attributes
    ----------
    kinetic : float
        ``0.5 * ∫ π² dV``
    gradient : float
        ``0.5 * ∫ |∇φ|² dV``
    mass : float
        ``0.5 * m² * ∫ φ² dV``
    total : float
        Sum of kinetic + gradient + mass.
    """

    kinetic: float
    gradient: float
    mass: float
    total: float


@dataclass(frozen=True)
class SystemEnergy:
    """Energy for the full coupled system at one snapshot.

    Attributes
    ----------
    per_field : dict[str, FieldEnergy]
        Energy breakdown per field.
    interaction : float
        ``0.5 * Σ_{i≠j} C_ij * ∫ φ_i φ_j dV``
    total : float
        Sum of all per-field energies + interaction.
    """

    per_field: dict[str, FieldEnergy]
    interaction: float
    total: float


# ------------------------------------------------------------------
# Gradient computation
# ------------------------------------------------------------------


def _gradient_energy_density(
    field: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
) -> NDArray[np.float64]:
    """Compute ``|∇φ|²`` on the grid.

    For periodic axes, uses spectral (FFT) gradient for accuracy.
    For non-periodic axes, uses 2nd-order central differences.
    """
    grad_sq: NDArray[np.float64] = np.zeros_like(field)

    for axis in range(len(grid_spacing)):
        dx = grid_spacing[axis]
        if periodic[axis]:
            # FFT-based gradient: exact for periodic fields
            n = field.shape[axis]
            freq = np.fft.fftfreq(n, d=dx) * (2.0 * np.pi)
            # Build shape for broadcasting: (1, ..., n, ..., 1)
            shape = [1] * field.ndim
            shape[axis] = n
            ik = 1j * freq.reshape(shape)
            fhat = np.fft.fft(field, axis=axis)
            grad_axis = np.real(np.fft.ifft(ik * fhat, axis=axis))
        else:
            # Finite difference gradient
            grad_axis = np.gradient(field, dx, axis=axis)

        grad_sq += grad_axis**2

    return grad_sq


# ------------------------------------------------------------------
# Single-field energy
# ------------------------------------------------------------------


def compute_field_energy(
    field_data: NDArray[np.float64],
    momentum_data: NDArray[np.float64] | None,
    mass_squared: float,
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
) -> FieldEnergy:
    """Compute canonical energy for a single field at one snapshot.

    Parameters
    ----------
    field_data : ndarray, shape ``(*grid_shape)``
        Field values on the spatial grid.
    momentum_data : ndarray or None
        Conjugate momentum ``π``.  ``None`` for constraint fields.
    mass_squared : float
        Diagonal mass matrix entry ``m²``.
    grid_spacing : tuple[float, ...]
        Cell size per spatial axis.
    periodic : tuple[bool, ...]
        Per-axis periodicity.

    Returns
    -------
    FieldEnergy
    """
    _validate_array(field_data, "field_data")
    if momentum_data is not None:
        _validate_array(momentum_data, "momentum_data")

    dv = float(np.prod(np.array(grid_spacing)))

    # Kinetic energy: 0.5 * ∫ π² dV
    if momentum_data is not None:
        kinetic = 0.5 * float(np.sum(momentum_data**2)) * dv
    else:
        kinetic = 0.0

    # Gradient energy: 0.5 * ∫ |∇φ|² dV
    grad_sq = _gradient_energy_density(field_data, grid_spacing, periodic)
    gradient = 0.5 * float(np.sum(grad_sq)) * dv

    # Mass energy: 0.5 * m² * ∫ φ² dV
    mass_energy = 0.5 * mass_squared * float(np.sum(field_data**2)) * dv

    total = kinetic + gradient + mass_energy
    return FieldEnergy(kinetic=kinetic, gradient=gradient, mass=mass_energy, total=total)


# ------------------------------------------------------------------
# System energy (all fields + interaction)
# ------------------------------------------------------------------


def _check_position_independent(
    data: SimulationData, eq_idx: int, target_field: str, label: str,
) -> None:
    """Raise if the identity term for *target_field* in equation *eq_idx* is position-dependent.

    Raises
    ------
    ValueError
        If the identity operator for *target_field* has position-dependent coefficients.
    """
    eq = data.spec.equations[eq_idx]
    for term in eq.rhs_terms:
        if term.operator == "identity" and term.field == target_field and term.position_dependent:
            msg = (
                f"Position-dependent {label} term in equation '{eq.field_name}' "
                f"targeting '{target_field}' — energy computation requires "
                f"constant coefficients (coordinate_dependent={term.coordinate_dependent})"
            )
            raise ValueError(msg)


def _resolve_mass_squared(
    data: SimulationData, field_idx: int,
) -> float:
    """Get the diagonal mass matrix entry for a field.

    Tries symbolic resolution first (using ``data.parameters``), then
    falls back to the pre-computed numeric matrix.
    """
    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
    )

    field_name = data.spec.component_names[field_idx]
    _check_position_independent(data, field_idx, field_name, "mass")

    sym_row = data.spec.mass_matrix_symbolic[field_idx] if data.spec.mass_matrix_symbolic else ()
    if field_idx < len(sym_row):
        sym_val = sym_row[field_idx]
        if sym_val is not None:
            resolved = _resolve_symbolic_coeff(sym_val, data.parameters)
            if resolved is not None:
                return float(resolved)

    return float(data.spec.mass_matrix[field_idx][field_idx])


def _resolve_coupling(
    data: SimulationData, i: int, j: int,
) -> float:
    """Get the off-diagonal coupling matrix entry C_ij."""
    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
    )

    target_field = data.spec.component_names[j]
    _check_position_independent(data, i, target_field, "coupling")

    sym_row = data.spec.coupling_matrix_symbolic[i] if data.spec.coupling_matrix_symbolic else ()
    if j < len(sym_row):
        sym_val = sym_row[j]
        if sym_val is not None:
            resolved = _resolve_symbolic_coeff(sym_val, data.parameters)
            if resolved is not None:
                return float(resolved)

    return float(data.spec.coupling_matrix[i][j])


def compute_system_energy(
    data: SimulationData,
    t_idx: int,
) -> SystemEnergy:
    """Compute total system energy at snapshot *t_idx*.

    Skips constraint fields (``time_derivative_order == 0``) as they are
    not dynamical degrees of freedom.

    Parameters
    ----------
    data : SimulationData
    t_idx : int
        Snapshot index.

    Raises
    ------
    ValueError
        If *t_idx* is out of range.
    """
    if t_idx < 0 or t_idx >= data.n_snapshots:
        msg = f"t_idx={t_idx} out of range [0, {data.n_snapshots})"
        raise ValueError(msg)

    names = data.spec.component_names
    per_field: dict[str, FieldEnergy] = {}

    for field_idx, eq in enumerate(data.spec.equations):
        name = eq.field_name
        if eq.time_derivative_order == 0:
            # Constraint field — no dynamical energy
            continue

        field_snapshot = data.fields[name][t_idx]
        mom_snapshot = data.momenta.get(name)
        mom_arr = mom_snapshot[t_idx] if mom_snapshot is not None else None

        m2 = _resolve_mass_squared(data, field_idx)

        per_field[name] = compute_field_energy(
            field_snapshot,
            mom_arr,
            m2,
            data.grid_spacing,
            data.periodic,
        )

    # Interaction energy: 0.5 * Σ_{i≠j} C_ij * ∫ φ_i φ_j dV
    dv = data.volume_element
    interaction = 0.0
    for i, eq_i in enumerate(data.spec.equations):
        if eq_i.time_derivative_order == 0:
            continue
        for j, eq_j in enumerate(data.spec.equations):
            if i == j or eq_j.time_derivative_order == 0:
                continue
            c_ij = _resolve_coupling(data, i, j)
            if c_ij == 0.0:
                continue
            phi_i = data.fields[names[i]][t_idx]
            phi_j = data.fields[names[j]][t_idx]
            interaction += 0.5 * c_ij * float(np.sum(phi_i * phi_j)) * dv

    total = sum(fe.total for fe in per_field.values()) + interaction
    return SystemEnergy(per_field=per_field, interaction=interaction, total=total)


def compute_energy_timeseries(
    data: SimulationData,
) -> tuple[
    NDArray[np.float64],
    dict[str, NDArray[np.float64]],
    NDArray[np.float64],
    NDArray[np.float64],
]:
    """Compute energy for every snapshot in the simulation.

    Returns
    -------
    times : ndarray, shape ``(n_snapshots,)``
    per_field : dict[str, ndarray]
        Each value is shape ``(n_snapshots,)`` — total energy of that field.
    interaction : ndarray, shape ``(n_snapshots,)``
    total : ndarray, shape ``(n_snapshots,)``
    """
    n = data.n_snapshots
    per_field_arrays: dict[str, list[float]] = {}
    interaction_list: list[float] = []
    total_list: list[float] = []

    for t_idx in range(n):
        se = compute_system_energy(data, t_idx)
        for name, fe in se.per_field.items():
            per_field_arrays.setdefault(name, []).append(fe.total)
        interaction_list.append(se.interaction)
        total_list.append(se.total)

    per_field_np: dict[str, NDArray[np.float64]] = {
        name: np.array(vals, dtype=np.float64)
        for name, vals in per_field_arrays.items()
    }

    return (
        data.times.copy(),
        per_field_np,
        np.array(interaction_list, dtype=np.float64),
        np.array(total_list, dtype=np.float64),
    )


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------


def _validate_array(arr: NDArray[np.float64], label: str) -> None:
    """Check array for non-finite values.

    Raises
    ------
    ValueError
        If *arr* contains NaN or Inf.
    """
    if not np.all(np.isfinite(arr)):
        msg = f"{label} contains NaN or Inf values"
        raise ValueError(msg)
