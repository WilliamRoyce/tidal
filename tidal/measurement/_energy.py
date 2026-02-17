"""Hamiltonian energy computation from Lagrangian-derived equations.

Computes the complete Hamiltonian H = T + V for any quadratic Lagrangian
by reconstructing it from the Euler-Lagrange equations in the JSON spec:

    H = ½ Σ_{dyn} ∫ π_sim² dV     (kinetic, using simulation momenta)
      + V_virial                    (from dynamical fields' spatial RHS terms)
      + V_constraint_self           (constraint field gradient + mass, sign-flipped)

The virial potential uses Euler's homogeneous function theorem for degree-2
functionals: V = -½ Σ ∫ φ_i · RHS_i^{spatial} dV.

Constraint fields (temporal gauge components) have NEGATIVE self-energy
due to the Minkowski metric g^{00} = -1.  This sign flip is automatic.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData
    from tidal.symbolic.json_loader import ComponentEquation, OperatorTerm

# Threshold below which energy is treated as zero.  Prevents division by
# near-zero values from floating-point integration noise.
_ENERGY_FLOOR: float = 1e-12

# Axis letter → numpy axis index (spatial axes only).
_AXIS_MAP: dict[str, int] = {"x": 0, "y": 1, "z": 2}

# Pattern for momentum field references: pi_0, pi_1, pi0, pi1, etc.
_MOMENTUM_RE = re.compile(r"^pi_?(\d+)$")


@dataclass(frozen=True)
class FieldEnergy:
    """Energy decomposition for a single field at one snapshot.

    Attributes
    ----------
    kinetic : float
        ``0.5 * ∫ π² dV``
    gradient : float
        ``0.5 * ∫ |∇_self φ|² dV`` — gradient energy over self-laplacian
        axes only.  For scalar fields (full ``laplacian``), this equals
        ``0.5 * ∫ |∇φ|² dV``.  For vector components with directional
        laplacians (e.g. ``laplacian_y``), only the relevant axes.
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
        Energy breakdown per field (operator-aware gradient).
    interaction : float
        Cross-field coupling energy: total potential minus per-field
        self-potentials.  Uses operator-aware gradient axes, so this is
        zero when fields are uncoupled and only one field is excited.
    total : float
        Complete Hamiltonian: kinetic + virial potential + constraint self-energy.
    """

    per_field: dict[str, FieldEnergy]
    interaction: float
    total: float


# ------------------------------------------------------------------
# Derivative helpers
# ------------------------------------------------------------------


def _pad_dirichlet(
    field: NDArray[np.float64],
    axis: int,
) -> NDArray[np.float64]:
    """Pad *field* with one anti-symmetric ghost cell per side on *axis*.

    For Dirichlet BCs (``field = 0`` at the boundary), the ghost cell
    value is the negative of the adjacent interior cell.  This is
    equivalent to ``np.pad(..., mode='reflect', reflect_type='odd')``.
    """
    pad_width = [(0, 0)] * field.ndim
    pad_width[axis] = (1, 1)
    return np.pad(field, pad_width, mode="reflect", reflect_type="odd")


def _first_derivative(
    field: NDArray[np.float64],
    axis: int,
    dx: float,
    *,
    is_periodic: bool,
) -> NDArray[np.float64]:
    """Single-axis first derivative.

    Uses FFT for periodic axes, central differences with Dirichlet
    ghost cells for non-periodic axes.
    """
    if is_periodic:
        n = field.shape[axis]
        freq = np.fft.fftfreq(n, d=dx) * (2.0 * np.pi)
        shape = [1] * field.ndim
        shape[axis] = n
        ik = 1j * freq.reshape(shape)
        fhat = np.fft.fft(field, axis=axis)
        return np.real(np.fft.ifft(ik * fhat, axis=axis))

    # Central difference with Dirichlet ghost cells: (f[i+1] - f[i-1]) / (2dx)
    padded = _pad_dirichlet(field, axis)
    slc_plus: list[slice] = [slice(None)] * field.ndim
    slc_minus: list[slice] = [slice(None)] * field.ndim
    slc_plus[axis] = slice(2, None)
    slc_minus[axis] = slice(None, -2)
    return (padded[tuple(slc_plus)] - padded[tuple(slc_minus)]) / (2.0 * dx)


def _second_derivative(
    field: NDArray[np.float64],
    axis: int,
    dx: float,
    *,
    is_periodic: bool,
) -> NDArray[np.float64]:
    """Single-axis second derivative.

    Uses ``-k²`` in FFT for periodic axes, central differences with
    Dirichlet ghost cells for non-periodic axes.
    """
    if is_periodic:
        n = field.shape[axis]
        freq = np.fft.fftfreq(n, d=dx) * (2.0 * np.pi)
        shape = [1] * field.ndim
        shape[axis] = n
        neg_k2 = -(freq**2).reshape(shape)
        fhat = np.fft.fft(field, axis=axis)
        return np.real(np.fft.ifft(neg_k2 * fhat, axis=axis))

    # Standard 3-point stencil with Dirichlet ghost cells:
    # (f[i+1] - 2f[i] + f[i-1]) / dx²
    padded = _pad_dirichlet(field, axis)
    slc_center: list[slice] = [slice(None)] * field.ndim
    slc_plus: list[slice] = [slice(None)] * field.ndim
    slc_minus: list[slice] = [slice(None)] * field.ndim
    slc_center[axis] = slice(1, -1)
    slc_plus[axis] = slice(2, None)
    slc_minus[axis] = slice(None, -2)
    return (
        padded[tuple(slc_plus)]
        - 2.0 * padded[tuple(slc_center)]
        + padded[tuple(slc_minus)]
    ) / (dx * dx)


# ------------------------------------------------------------------
# Operator-aware gradient axes
# ------------------------------------------------------------------


def _self_gradient_axes(eq: ComponentEquation) -> list[int] | None:
    """Return spatial axes that have self-laplacian operators in *eq*.

    Inspects ``eq.rhs_terms`` for laplacian-type operators acting on the
    equation's own field.  Returns:

    - ``None`` if the equation contains a full ``laplacian`` (all axes).
    - A sorted list of axis indices for directional laplacians
      (e.g. ``[1]`` for ``laplacian_y``).

    This is used to compute operator-aware per-field gradient energy:
    only the gradient axes that appear as self-laplacian terms contribute
    to per-field energy; the remaining axes contribute to interaction.
    """
    axes: set[int] = set()
    for term in eq.rhs_terms:
        if term.field != eq.field_name:
            continue
        if term.operator == "laplacian":
            return None  # full laplacian → use all axes
        if term.operator.startswith("laplacian_"):
            letter = term.operator[len("laplacian_"):]
            if letter in _AXIS_MAP:
                axes.add(_AXIS_MAP[letter])
    return sorted(axes)


# ------------------------------------------------------------------
# Gradient energy density (used by compute_field_energy)
# ------------------------------------------------------------------


def _gradient_energy_density(
    field: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    axes: list[int] | None = None,
) -> NDArray[np.float64]:
    """Compute ``|∇φ|²`` on the grid.

    Parameters
    ----------
    axes : list[int] | None
        If ``None``, sum over all spatial axes (isotropic gradient).
        Otherwise, sum only over the specified axis indices.
    """
    grad_sq: NDArray[np.float64] = np.zeros_like(field)
    iter_axes = range(len(grid_spacing)) if axes is None else axes
    for axis in iter_axes:
        grad_axis = _first_derivative(field, axis, grid_spacing[axis], is_periodic=periodic[axis])
        grad_sq += grad_axis**2
    return grad_sq


# ------------------------------------------------------------------
# Spatial operator application (for virial potential)
# ------------------------------------------------------------------


def _apply_spatial_operator(
    operator: str,
    field: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
) -> NDArray[np.float64]:
    """Apply a named spatial operator to a field array.

    Raises
    ------
    ValueError
        If the operator is unknown.
    """
    if operator == "identity":
        return field.copy()

    # gradient_{x,y,z}
    if operator.startswith("gradient_"):
        axis_letter = operator[len("gradient_"):]
        if axis_letter in _AXIS_MAP:
            ax = _AXIS_MAP[axis_letter]
            return _first_derivative(field, ax, grid_spacing[ax], is_periodic=periodic[ax])

    # laplacian_{x,y,z}
    if operator.startswith("laplacian_"):
        axis_letter = operator[len("laplacian_"):]
        if axis_letter in _AXIS_MAP:
            ax = _AXIS_MAP[axis_letter]
            return _second_derivative(field, ax, grid_spacing[ax], is_periodic=periodic[ax])

    # laplacian (isotropic sum)
    if operator == "laplacian":
        result: NDArray[np.float64] = np.zeros_like(field)
        for ax in range(len(grid_spacing)):
            result += _second_derivative(field, ax, grid_spacing[ax], is_periodic=periodic[ax])
        return result

    # cross_derivative_{xy,xz,yz}
    if operator.startswith("cross_derivative_"):
        axes_str = operator[len("cross_derivative_"):]
        if len(axes_str) == 2 and axes_str[0] in _AXIS_MAP and axes_str[1] in _AXIS_MAP:  # noqa: PLR2004
            ax0 = _AXIS_MAP[axes_str[0]]
            ax1 = _AXIS_MAP[axes_str[1]]
            tmp = _first_derivative(field, ax0, grid_spacing[ax0], is_periodic=periodic[ax0])
            return _first_derivative(tmp, ax1, grid_spacing[ax1], is_periodic=periodic[ax1])

    msg = f"Unknown spatial operator for energy measurement: '{operator}'"
    raise ValueError(msg)


# ------------------------------------------------------------------
# Term resolution helpers
# ------------------------------------------------------------------


def _is_momentum_field(field_name: str) -> bool:
    """Check if a field name is a momentum reference (pi_N / piN)."""
    return _MOMENTUM_RE.match(field_name) is not None


def _resolve_term_coefficient(
    term: OperatorTerm,
    parameters: dict[str, float],
) -> float:
    """Resolve a term's numeric coefficient.

    Tries symbolic resolution first, falls back to numeric.
    No negation — RHS coefficients are already the correct sign.
    """
    if term.coefficient_symbolic is not None:
        from tidal.symbolic.json_loader import (  # noqa: PLC0415
            _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
        )

        resolved = _resolve_symbolic_coeff(term.coefficient_symbolic, parameters)
        if resolved is not None:
            return float(resolved)

    return float(term.coefficient)


def _build_coord_arrays(
    data: SimulationData,
) -> dict[str, NDArray[np.float64]]:
    """Build spatial coordinate arrays from ``SimulationData`` grid info.

    Constructs cell-centered coordinate meshgrids matching the field
    data shape (``*grid_shape``).
    """
    spatial_coords = data.spec.spatial_coordinates
    axes: list[NDArray[np.float64]] = []
    for i, name in enumerate(spatial_coords):
        lo, hi = data.grid_bounds[i]
        dx = data.grid_spacing[i]
        n = round((hi - lo) / dx)
        # Cell-centered: offset by dx/2
        axes.append(np.linspace(lo + dx / 2, hi - dx / 2, n))
        _ = name  # used below

    grids = np.meshgrid(*axes, indexing="ij")
    return {
        name: np.asarray(g, dtype=np.float64)
        for name, g in zip(spatial_coords, grids, strict=True)
    }


def _resolve_coefficient_on_grid(
    term: OperatorTerm,
    data: SimulationData,
    coord_arrays: dict[str, NDArray[np.float64]],
) -> float | NDArray[np.float64]:
    """Resolve a term's coefficient, returning a grid array if position-dependent.

    For constant coefficients, returns a scalar float (same as
    ``_resolve_term_coefficient``).  For position-dependent coefficients,
    evaluates the symbolic expression on the grid using
    :func:`~tidal.symbolic._eval_utils.evaluate_coefficient`.
    """
    if not term.position_dependent:
        return _resolve_term_coefficient(term, data.parameters)

    sym = term.coefficient_symbolic
    if sym is None:
        return float(term.coefficient)

    from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

    return evaluate_coefficient(
        sym,
        data.parameters,
        data.spec.coordinates,
        coord_arrays=coord_arrays,
        t=0.0,
    )


def _resolve_term_target(
    data: SimulationData,
    field_name: str,
    t_idx: int,
) -> NDArray[np.float64] | None:
    """Resolve the field data a term acts on.

    Returns
    -------
    NDArray or None
        The field/momentum snapshot, or ``None`` if the target is a
        zero-momentum constraint field (expected case).

    Raises
    ------
    ValueError
        If *field_name* cannot be resolved to any known field or momentum.
    """
    # Direct field reference
    if field_name in data.fields:
        return data.fields[field_name][t_idx]

    # Momentum reference: pi_N or piN
    m = _MOMENTUM_RE.match(field_name)
    if m is not None:
        idx = int(m.group(1))
        names = data.spec.component_names
        if idx >= len(names):
            msg = (
                f"Momentum reference '{field_name}' resolves to index {idx}, "
                f"but spec only has {len(names)} fields: {names}"
            )
            raise ValueError(msg)
        target_name = names[idx]
        # Constraint field → zero momentum (expected None)
        eq = data.spec.equations[idx]
        if eq.time_derivative_order == 0:
            return None
        mom = data.momenta.get(target_name)
        if mom is not None:
            return mom[t_idx]
        msg = (
            f"Momentum reference '{field_name}' resolves to field "
            f"'{target_name}', but no momentum data found"
        )
        raise ValueError(msg)

    msg = (
        f"Unresolvable field reference '{field_name}' — "
        f"not a known field ({list(data.fields.keys())}) "
        f"or momentum pattern (pi_N)"
    )
    raise ValueError(msg)


# ------------------------------------------------------------------
# Single-field energy
# ------------------------------------------------------------------


def compute_field_energy(  # noqa: PLR0913
    field_data: NDArray[np.float64],
    momentum_data: NDArray[np.float64] | None,
    mass_squared: float | NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    *,
    gradient_axes: list[int] | None = None,
) -> FieldEnergy:
    """Compute canonical energy for a single field at one snapshot.

    Parameters
    ----------
    field_data : ndarray, shape ``(*grid_shape)``
        Field values on the spatial grid.
    momentum_data : ndarray or None
        Conjugate momentum ``π``.  ``None`` for constraint fields.
    mass_squared : float | ndarray
        Diagonal mass matrix entry ``m²``.  May be a scalar (constant mass)
        or a grid-shaped ndarray (position-dependent mass).
    grid_spacing : tuple[float, ...]
        Cell size per spatial axis.
    periodic : tuple[bool, ...]
        Per-axis periodicity.
    gradient_axes : list[int] | None
        Spatial axes to include in gradient energy.  ``None`` uses all
        axes (isotropic gradient).  Pass a subset for operator-aware
        gradient (e.g. ``[1]`` when the PDE has only ``laplacian_y``).

    Returns
    -------
    FieldEnergy
    """
    _validate_array(field_data, "field_data")
    if momentum_data is not None:
        _validate_array(momentum_data, "momentum_data")

    dv = float(np.array(grid_spacing).prod())

    # Kinetic energy: 0.5 * ∫ π² dV
    if momentum_data is not None:
        kinetic = 0.5 * float((momentum_data**2).sum()) * dv
    else:
        kinetic = 0.0

    # Gradient energy: 0.5 * ∫ |∇φ|² dV (over specified axes)
    grad_sq = _gradient_energy_density(field_data, grid_spacing, periodic, axes=gradient_axes)
    gradient = 0.5 * float(grad_sq.sum()) * dv

    # Mass energy: 0.5 * m² * ∫ φ² dV (m² may be scalar or ndarray)
    mass_energy = 0.5 * float((mass_squared * field_data**2).sum()) * dv

    total = kinetic + gradient + mass_energy
    return FieldEnergy(kinetic=kinetic, gradient=gradient, mass=mass_energy, total=total)


# ------------------------------------------------------------------
# System energy (virial + constraint)
# ------------------------------------------------------------------


def _resolve_mass_squared(
    data: SimulationData,
    field_idx: int,
    coord_arrays: dict[str, NDArray[np.float64]] | None = None,
) -> float | NDArray[np.float64]:
    """Get the diagonal mass matrix entry for a field.

    Tries symbolic resolution first (using ``data.parameters``), then
    falls back to the pre-computed numeric matrix.  For position-dependent
    mass terms, evaluates the symbolic expression on the grid and returns
    an ndarray (negated, per mass matrix convention).

    Parameters
    ----------
    data : SimulationData
        Simulation data.
    field_idx : int
        Field equation index.
    coord_arrays : dict[str, NDArray] | None
        Pre-built coordinate arrays (from ``_build_coord_arrays``).
        Required if the mass term is position-dependent.
    """
    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
    )

    # Check for position-dependent mass term → evaluate on grid
    eq = data.spec.equations[field_idx]
    field_name = eq.field_name
    for term in eq.rhs_terms:
        if term.operator == "identity" and term.field == field_name and term.position_dependent:
            sym = term.coefficient_symbolic
            if sym is None:
                return float(term.coefficient)
            if coord_arrays is None:
                coord_arrays = _build_coord_arrays(data)
            from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

            coeff = evaluate_coefficient(
                sym, data.parameters, data.spec.coordinates,
                coord_arrays=coord_arrays,
            )
            # Convention: mass_matrix[i][i] = -(coefficient of identity(field_i))
            if isinstance(coeff, np.ndarray):
                return -coeff
            return -float(coeff)

    sym_row = data.spec.mass_matrix_symbolic[field_idx] if data.spec.mass_matrix_symbolic else ()
    if field_idx < len(sym_row):
        sym_val = sym_row[field_idx]
        if sym_val is not None:
            resolved = _resolve_symbolic_coeff(sym_val, data.parameters)
            if resolved is not None:
                # Negate: symbolic matrix stores raw coefficient_symbolic,
                # but convention is matrix[i][j] = -(coefficient).
                return float(-resolved)

    return float(data.spec.mass_matrix[field_idx][field_idx])


def _compute_virial_potential(
    data: SimulationData,
    t_idx: int,
) -> float:
    """Virial potential from dynamical fields' spatial RHS terms.

    ``V_virial = -½ Σ_{i: dynamical} ∫ φ_i · RHS_i^{spatial} dV``

    Excludes ``first_derivative_t`` (gyroscopic, do no work) and
    ``pi_N`` momentum references (velocity-dependent forces).

    Supports position-dependent coefficients by evaluating them on the
    grid and performing elementwise integration.
    """
    dv = data.volume_element
    potential = 0.0

    # Build coordinate arrays once (lazy, only if needed)
    coord_arrays: dict[str, NDArray[np.float64]] | None = None
    has_posdep = any(
        term.position_dependent
        for eq in data.spec.equations
        for term in eq.rhs_terms
    )
    if has_posdep:
        coord_arrays = _build_coord_arrays(data)

    for eq in data.spec.equations:
        if eq.time_derivative_order < 2:  # noqa: PLR2004
            continue  # skip constraints and first-order

        phi_i = data.fields[eq.field_name][t_idx]

        for term in eq.rhs_terms:
            # Skip time-derivative terms (gyroscopic forces)
            if term.operator == "first_derivative_t":
                continue

            # Skip momentum-field references (velocity-dependent)
            if _is_momentum_field(term.field):
                continue

            target = _resolve_term_target(data, term.field, t_idx)
            if target is None:
                continue

            coeff = _resolve_coefficient_on_grid(term, data, coord_arrays or {})
            operated = _apply_spatial_operator(
                term.operator, target, data.grid_spacing, data.periodic,
            )
            # coeff may be scalar or ndarray — numpy handles both
            potential += float((coeff * phi_i * operated).sum()) * dv

    return -0.5 * potential


def _compute_constraint_self_energy(
    data: SimulationData,
    t_idx: int,
) -> float:
    """Constraint field self-energy with sign flip (g^{00} = -1).

    ``V_constraint = Σ_{j: constraint} [ -½ |∇C_j|² - ½ m_j² C_j² ] dV``

    Temporal gauge components have NEGATIVE gradient and mass self-energy
    relative to spatial/scalar fields, due to the Minkowski metric.
    """
    energy = 0.0

    for field_idx, eq in enumerate(data.spec.equations):
        if eq.time_derivative_order != 0:
            continue  # only constraint fields
        name = eq.field_name
        if name not in data.fields:
            continue  # constraint field not stored in data

        c_field = data.fields[name][t_idx]
        dv = data.volume_element

        # Gradient: -½ ∫ |∇C|² dV  (NEGATIVE)
        grad_sq = _gradient_energy_density(c_field, data.grid_spacing, data.periodic)
        energy -= 0.5 * float(grad_sq.sum()) * dv

        # Mass: -½ m² ∫ C² dV  (NEGATIVE, m² may be scalar or ndarray)
        m2 = _resolve_mass_squared(data, field_idx)
        energy -= 0.5 * float((m2 * c_field**2).sum()) * dv

    return energy


def compute_system_energy(  # noqa: PLR0914
    data: SimulationData,
    t_idx: int,
) -> SystemEnergy:
    """Compute total Hamiltonian energy at snapshot *t_idx*.

    Uses the complete formula derived from the Lagrangian:

        H = ½ Σ π_sim² + V_virial + V_constraint_self

    where V_virial is computed from dynamical equations' spatial RHS terms
    (Euler's theorem for quadratic functionals), and V_constraint_self
    accounts for temporal gauge components' negative self-energy.

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

    # Per-field canonical energy (kinetic + gradient + mass) — dynamical only.
    # Gradient uses operator-aware axes: only the spatial axes that appear
    # as self-laplacian operators in each field's equation.  For scalar
    # fields with a full ``laplacian``, this is all axes (unchanged).
    # For vector components with directional laplacians (e.g. laplacian_y),
    # only the corresponding axes contribute to per-field gradient energy.

    # Pre-build coordinate arrays once if any field has position-dependent mass
    coord_arrays: dict[str, NDArray[np.float64]] | None = None
    has_posdep_mass = any(
        term.operator == "identity"
        and term.field == eq.field_name
        and term.position_dependent
        for eq in data.spec.equations
        for term in eq.rhs_terms
    )
    if has_posdep_mass:
        coord_arrays = _build_coord_arrays(data)

    per_field: dict[str, FieldEnergy] = {}
    for field_idx, eq in enumerate(data.spec.equations):
        name = eq.field_name
        if eq.time_derivative_order == 0:
            continue

        field_snapshot = data.fields[name][t_idx]
        mom_snapshot = data.momenta.get(name)
        mom_arr = mom_snapshot[t_idx] if mom_snapshot is not None else None

        m2 = _resolve_mass_squared(data, field_idx, coord_arrays=coord_arrays)
        axes = _self_gradient_axes(eq)
        per_field[name] = compute_field_energy(
            field_snapshot, mom_arr, m2, data.grid_spacing, data.periodic,
            gradient_axes=axes,
        )

    # Virial potential from dynamical fields' spatial RHS terms
    v_virial = _compute_virial_potential(data, t_idx)

    # Constraint field self-energy (negative, from g^{00} = -1)
    v_constraint = _compute_constraint_self_energy(data, t_idx)

    # Total potential = virial + constraint
    v_total = v_virial + v_constraint

    # Interaction = total potential minus per-field self-potentials
    self_potential = sum(fe.gradient + fe.mass for fe in per_field.values())
    interaction = v_total - self_potential

    # Total energy = kinetic + total potential
    total_kinetic = sum(fe.kinetic for fe in per_field.values())
    total = total_kinetic + v_total

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
    if not np.isfinite(arr).all():
        msg = f"{label} contains NaN or Inf values"
        raise ValueError(msg)
