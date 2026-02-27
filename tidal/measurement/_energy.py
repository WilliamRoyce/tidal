"""Hamiltonian energy density computation from Lagrangian-derived equations.

Computes the spatially-averaged Hamiltonian energy density ⟨ε⟩ = H / V_domain
for any quadratic Lagrangian by reconstructing it from the Euler-Lagrange
equations in the JSON spec:

    ⟨ε⟩ = ½ Σ_{dyn} ⟨v²⟩             (kinetic, using simulation velocities)
         + ⟨v_virial⟩                 (from dynamical fields' spatial RHS terms)
         + ⟨v_constraint_self⟩        (constraint field gradient + mass, sign-flipped)
         + ⟨v_constraint_cross⟩       (cross-constraint identity coupling)

The virial potential density uses Euler's homogeneous function theorem for
degree-2 functionals: ⟨v⟩ = -½ Σ ⟨φ_i · RHS_i^{spatial}⟩.

Constraint fields (temporal gauge components) have NEGATIVE self-energy
due to the Minkowski metric g^{00} = -1.  This sign flip is automatic.

All energy values are **average energy densities** (intensive quantities),
independent of the domain size.  Ratios (conservation, conversion probability)
are invariant under this normalization.
"""

from __future__ import annotations

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

# Gradient operator → axis index.  Used by the integration-by-parts
# Hamiltonian evaluation to convert gradient-product terms into the
# equivalent second-order operators (laplacian / cross_derivative).
_GRADIENT_AXES: dict[str, int] = {
    "gradient_x": 0,
    "gradient_y": 1,
    "gradient_z": 2,
}


@dataclass(frozen=True)
class FieldEnergy:
    """Energy density decomposition for a single field at one snapshot.

    All values are spatially-averaged energy densities ⟨ε⟩ = E / V_domain.

    Attributes
    ----------
    kinetic : float
        ``0.5 * ⟨π²⟩``
    gradient : float
        ``0.5 * ⟨|∇_self φ|²⟩`` — gradient energy density over
        self-laplacian axes only.  For scalar fields (full ``laplacian``),
        this equals ``0.5 * ⟨|∇φ|²⟩``.  For vector components with
        directional laplacians (e.g. ``laplacian_y``), only the relevant axes.
    mass : float
        ``0.5 * m² * ⟨φ²⟩``
    total : float
        Sum of kinetic + gradient + mass.
    """

    kinetic: float
    gradient: float
    mass: float
    total: float


@dataclass(frozen=True)
class SystemEnergy:
    """Energy density for the full coupled system at one snapshot.

    All values are spatially-averaged energy densities ⟨ε⟩ = E / V_domain.

    Attributes
    ----------
    per_field : dict[str, FieldEnergy]
        Energy density breakdown per field (operator-aware gradient).
    interaction : float
        Cross-field coupling energy density: total potential density minus
        per-field self-potential densities.  Uses operator-aware gradient
        axes, so this is zero when fields are uncoupled and only one field
        is excited.
    total : float
        Complete Hamiltonian density: kinetic + virial + constraint self-energy.
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


def _pad_neumann(
    field: NDArray[np.float64],
    axis: int,
) -> NDArray[np.float64]:
    """Pad *field* with one symmetric ghost cell per side on *axis*.

    For Neumann BCs (``∂field/∂n = 0`` at the boundary), the ghost cell
    value equals the adjacent interior cell.  This is equivalent to
    ``np.pad(..., mode='reflect', reflect_type='even')``.
    """
    pad_width = [(0, 0)] * field.ndim
    pad_width[axis] = (1, 1)
    return np.pad(field, pad_width, mode="reflect", reflect_type="even")


def _first_derivative(
    field: NDArray[np.float64],
    axis: int,
    dx: float,
    *,
    is_periodic: bool = False,
    bc_type: str | None = None,
) -> NDArray[np.float64]:
    """Single-axis first derivative via central differences.

    Uses ``np.roll`` for periodic wrapping, ghost-cell padding for
    non-periodic axes.  Both paths use the same 2-point central stencil
    ``(f[i+1] - f[i-1]) / (2 dx)``, matching py-pde's finite-difference
    operators so the virial energy tracks the PDE solver's conserved
    Hamiltonian exactly.

    Parameters
    ----------
    bc_type : str or None
        Boundary condition type: ``"periodic"``, ``"neumann"``, or
        ``"dirichlet"``.  When provided, overrides *is_periodic*.
    is_periodic : bool
        Legacy parameter — used only when *bc_type* is ``None``.
    """
    effective_bc = (
        bc_type if bc_type is not None else ("periodic" if is_periodic else "dirichlet")
    )

    if effective_bc == "periodic":
        return (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)) / (
            2.0 * dx
        )

    # Central difference with ghost cells: (f[i+1] - f[i-1]) / (2dx)
    padded = (
        _pad_neumann(field, axis)
        if effective_bc == "neumann"
        else _pad_dirichlet(field, axis)
    )
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
    is_periodic: bool = False,
    bc_type: str | None = None,
) -> NDArray[np.float64]:
    """Single-axis second derivative via 3-point central stencil.

    Uses ``np.roll`` for periodic wrapping, ghost-cell padding for
    non-periodic axes.  Both paths use ``(f[i+1] - 2f[i] + f[i-1]) / dx²``,
    matching py-pde's finite-difference operators so the virial energy
    tracks the PDE solver's conserved Hamiltonian exactly.

    Parameters
    ----------
    bc_type : str or None
        Boundary condition type: ``"periodic"``, ``"neumann"``, or
        ``"dirichlet"``.  When provided, overrides *is_periodic*.
    is_periodic : bool
        Legacy parameter — used only when *bc_type* is ``None``.
    """
    effective_bc = (
        bc_type if bc_type is not None else ("periodic" if is_periodic else "dirichlet")
    )

    if effective_bc == "periodic":
        return (
            np.roll(field, -1, axis=axis) - 2.0 * field + np.roll(field, 1, axis=axis)
        ) / (dx * dx)

    # Standard 3-point stencil with ghost cells:
    # (f[i+1] - 2f[i] + f[i-1]) / dx²
    padded = (
        _pad_neumann(field, axis)
        if effective_bc == "neumann"
        else _pad_dirichlet(field, axis)
    )
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
            letter = term.operator[len("laplacian_") :]
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
    bc_types: tuple[str, ...] | None = None,
) -> NDArray[np.float64]:
    """Gradient energy density consistent with the virial potential formula.

    For **periodic** axes, uses the Laplacian-based identity
    ``-φ · ∂²φ/∂x²`` so that ``0.5 * Σ * dV`` exactly matches the
    virial potential's self-laplacian contribution (discrete integration
    by parts is exact for periodic wrapping).

    For **non-periodic** axes, uses ``(∂φ/∂x)²`` (central-difference
    gradient squared), because the discrete IBP has nonzero boundary
    terms that would make the Laplacian form inconsistent.

    Parameters
    ----------
    axes : list[int] | None
        If ``None``, sum over all spatial axes (isotropic gradient).
        Otherwise, sum only over the specified axis indices.
    bc_types : tuple[str, ...] | None
        Per-axis BC type (``"periodic"``, ``"neumann"``, ``"dirichlet"``).
        When ``None``, falls back to ``periodic`` booleans (legacy behavior).
    """
    result: NDArray[np.float64] = np.zeros_like(field)
    iter_axes = range(len(grid_spacing)) if axes is None else axes
    for axis in iter_axes:
        dx = grid_spacing[axis]
        bc = (
            bc_types[axis]
            if bc_types is not None
            else ("periodic" if periodic[axis] else "dirichlet")
        )
        if bc == "periodic":
            # Virial-consistent: -φ · ∂²φ/∂x² (exact discrete IBP)
            result -= field * _second_derivative(field, axis, dx, bc_type="periodic")
        else:
            # Gradient squared: (∂φ/∂x)² (correct for non-periodic BCs)
            grad = _first_derivative(field, axis, dx, bc_type=bc)
            result += grad**2
    return result


# ------------------------------------------------------------------
# Spatial operator application (for virial potential)
# ------------------------------------------------------------------


def _effective_bc(
    axis: int,
    periodic: tuple[bool, ...],
    bc_types: tuple[str, ...] | None,
) -> str:
    """Resolve effective BC type for a single axis."""
    if bc_types is not None:
        return bc_types[axis]
    return "periodic" if periodic[axis] else "dirichlet"


def _apply_spatial_operator(
    operator: str,
    field: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    bc_types: tuple[str, ...] | None = None,
) -> NDArray[np.float64]:
    """Apply a named spatial operator to a field array.

    Parameters
    ----------
    bc_types : tuple[str, ...] | None
        Per-axis BC type.  When ``None``, falls back to ``periodic`` booleans.

    Raises
    ------
    ValueError
        If the operator is unknown.
    """
    if operator == "identity":
        return field.copy()

    # gradient_{x,y,z}
    if operator.startswith("gradient_"):
        axis_letter = operator[len("gradient_") :]
        if axis_letter in _AXIS_MAP:
            ax = _AXIS_MAP[axis_letter]
            return _first_derivative(
                field,
                ax,
                grid_spacing[ax],
                bc_type=_effective_bc(ax, periodic, bc_types),
            )

    # laplacian_{x,y,z}
    if operator.startswith("laplacian_"):
        axis_letter = operator[len("laplacian_") :]
        if axis_letter in _AXIS_MAP:
            ax = _AXIS_MAP[axis_letter]
            return _second_derivative(
                field,
                ax,
                grid_spacing[ax],
                bc_type=_effective_bc(ax, periodic, bc_types),
            )

    # laplacian (isotropic sum)
    if operator == "laplacian":
        result: NDArray[np.float64] = np.zeros_like(field)
        for ax in range(len(grid_spacing)):
            result += _second_derivative(
                field,
                ax,
                grid_spacing[ax],
                bc_type=_effective_bc(ax, periodic, bc_types),
            )
        return result

    # cross_derivative_{xy,xz,yz}
    if operator.startswith("cross_derivative_"):
        axes_str = operator[len("cross_derivative_") :]
        if len(axes_str) == 2 and axes_str[0] in _AXIS_MAP and axes_str[1] in _AXIS_MAP:  # noqa: PLR2004
            ax0 = _AXIS_MAP[axes_str[0]]
            ax1 = _AXIS_MAP[axes_str[1]]
            tmp = _first_derivative(
                field,
                ax0,
                grid_spacing[ax0],
                bc_type=_effective_bc(ax0, periodic, bc_types),
            )
            return _first_derivative(
                tmp,
                ax1,
                grid_spacing[ax1],
                bc_type=_effective_bc(ax1, periodic, bc_types),
            )

    msg = f"Unknown spatial operator for energy measurement: '{operator}'"
    raise ValueError(msg)


# ------------------------------------------------------------------
# Term resolution helpers
# ------------------------------------------------------------------


def _is_velocity_field(field_name: str) -> bool:
    """Check if a field name is a velocity reference (v_field_name, e.g. v_A_1)."""
    return field_name.startswith("v_") and len(field_name) > len("v_")


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
        del name  # used only in the dict comprehension below

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
        data.spec.effective_coordinates,
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
        The field/velocity snapshot, or ``None`` if the target is a
        zero-velocity constraint field (expected case).

    Raises
    ------
    ValueError
        If *field_name* cannot be resolved to any known field or velocity.
    """
    # Direct field reference
    if field_name in data.fields:
        return data.fields[field_name][t_idx]

    # Velocity reference: v_field_name (e.g. "v_A_1")
    if field_name.startswith("v_") and len(field_name) > len("v_"):
        suffix = field_name[2:]
        names = data.spec.component_names

        if suffix not in names:
            msg = (
                f"Velocity reference '{field_name}': suffix '{suffix}' "
                f"is not a known field ({names})"
            )
            raise ValueError(msg)

        target_name = suffix
        eq_idx = names.index(target_name)

        # Constraint field → zero velocity (expected None)
        eq = data.spec.equations[eq_idx]
        if eq.time_derivative_order == 0:
            return None
        vel = data.velocities.get(target_name)
        if vel is not None:
            return vel[t_idx]
        msg = (
            f"Velocity reference '{field_name}' resolves to field "
            f"'{target_name}', but no velocity data found"
        )
        raise ValueError(msg)

    msg = (
        f"Unresolvable field reference '{field_name}' — "
        f"not a known field ({list(data.fields.keys())}) "
        f"or velocity pattern (v_field_name)"
    )
    raise ValueError(msg)


# ------------------------------------------------------------------
# Single-field energy
# ------------------------------------------------------------------


def compute_field_energy(  # noqa: PLR0913
    field_data: NDArray[np.float64],
    velocity_data: NDArray[np.float64] | None,
    mass_squared: float | NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    *,
    gradient_axes: list[int] | None = None,
    bc_types: tuple[str, ...] | None = None,
    volume_weight: float | NDArray[np.float64] = 1.0,
) -> FieldEnergy:
    """Compute canonical energy density for a single field at one snapshot.

    Returns spatially-averaged energy density ⟨ε⟩ = E / V_domain.

    Parameters
    ----------
    field_data : ndarray, shape ``(*grid_shape)``
        Field values on the spatial grid.
    velocity_data : ndarray or None
        Velocity ``v = dq/dt``.  ``None`` for constraint fields.
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
    bc_types : tuple[str, ...] | None
        Per-axis BC type.  When ``None``, falls back to ``periodic``.
    volume_weight : float or ndarray
        Spatial volume element ``sqrt|g_spatial|`` for curved coordinates.
        Defaults to 1.0 (flat/Cartesian).

    Returns
    -------
    FieldEnergy
    """
    _validate_array(field_data, "field_data")
    if velocity_data is not None:
        _validate_array(velocity_data, "velocity_data")

    # Kinetic energy density: 0.5 * ⟨v² * sqrt|g|⟩
    if velocity_data is not None:
        kinetic = 0.5 * float((velocity_data**2 * volume_weight).mean())
    else:
        kinetic = 0.0

    # Gradient energy density: 0.5 * ⟨|∇φ|² * sqrt|g|⟩ (over specified axes)
    grad_sq = _gradient_energy_density(
        field_data,
        grid_spacing,
        periodic,
        axes=gradient_axes,
        bc_types=bc_types,
    )
    gradient = 0.5 * float((grad_sq * volume_weight).mean())

    # Mass energy density: 0.5 * ⟨m² φ² * sqrt|g|⟩ (m² may be scalar or ndarray)
    mass_energy = 0.5 * float((mass_squared * field_data**2 * volume_weight).mean())

    total = kinetic + gradient + mass_energy
    return FieldEnergy(
        kinetic=kinetic, gradient=gradient, mass=mass_energy, total=total
    )


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
        if (
            term.operator == "identity"
            and term.field == field_name
            and term.position_dependent
        ):
            sym = term.coefficient_symbolic
            if sym is None:
                return float(term.coefficient)
            if coord_arrays is None:
                coord_arrays = _build_coord_arrays(data)
            from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

            coeff = evaluate_coefficient(
                sym,
                data.parameters,
                data.spec.effective_coordinates,
                coord_arrays=coord_arrays,
            )
            # Convention: mass_matrix[i][i] = -(coefficient of identity(field_i))
            if isinstance(coeff, np.ndarray):
                return -coeff
            return -float(coeff)

    sym_row = (
        data.spec.mass_matrix_symbolic[field_idx]
        if data.spec.mass_matrix_symbolic
        else ()
    )
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
    """Virial potential density from dynamical fields' spatial RHS terms.

    ``⟨v_virial⟩ = -½ Σ_{i: dynamical} ⟨φ_i · RHS_i^{spatial}⟩``

    Excludes ``first_derivative_t`` (gyroscopic, do no work) and
    ``v_N`` velocity references (velocity-dependent forces).

    Supports position-dependent coefficients by evaluating them on the
    grid and performing elementwise averaging.
    """
    potential = 0.0

    # Build coordinate arrays once (lazy, only if needed)
    coord_arrays: dict[str, NDArray[np.float64]] | None = None
    has_posdep = any(
        term.position_dependent for eq in data.spec.equations for term in eq.rhs_terms
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

            # Skip velocity-field references (velocity-dependent)
            if _is_velocity_field(term.field):
                continue

            target = _resolve_term_target(data, term.field, t_idx)
            if target is None:
                continue

            coeff = _resolve_coefficient_on_grid(term, data, coord_arrays or {})
            operated = _apply_spatial_operator(
                term.operator,
                target,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )
            # coeff may be scalar or ndarray — numpy handles both
            potential += float((coeff * phi_i * operated).mean())

    return -0.5 * potential


def _compute_constraint_self_energy(
    data: SimulationData,
    t_idx: int,
) -> float:
    """Constraint field self-energy density with sign flip (g^{00} = -1).

    ``⟨v_constraint⟩ = Σ_{j: constraint} [ -½ ⟨|∇C_j|²⟩ - ½ m_j² ⟨C_j²⟩ ]``

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

        # Gradient: -½ ⟨|∇C|²⟩  (NEGATIVE)
        grad_sq = _gradient_energy_density(
            c_field,
            data.grid_spacing,
            data.periodic,
            bc_types=data.bc_types,
        )
        energy -= 0.5 * float(grad_sq.mean())

        # Mass: -½ ⟨m² C²⟩  (NEGATIVE, m² may be scalar or ndarray)
        m2 = _resolve_mass_squared(data, field_idx)
        energy -= 0.5 * float((m2 * c_field**2).mean())

    return energy


def _compute_constraint_coupling_energy(
    data: SimulationData,
    t_idx: int,
) -> float:
    """Cross-constraint coupling energy density from RHS terms between constraints.

    For each constraint equation i with a term ``c * op(C_j)`` referencing
    another constraint field C_j, accumulates:

        ``⟨v_cross⟩ += +½ * c * ⟨C_i * op(C_j)⟩``

    The ``+½`` sign (opposite of the dynamical virial ``-½``) arises because
    constraint fields have ``π = 0``, so their Hamiltonian contribution is
    ``H = -L``.  The RHS coefficients in the JSON already embed the Minkowski
    sign, so the formula reproduces the correct Hamiltonian coupling.

    For symmetric coupling (c_ij in eq_i AND c_ji in eq_j), the two halves
    sum to the full coupling energy.

    This handles any spatial operator (identity, gradient, laplacian, etc.)
    between constraint fields, not just identity coupling.
    """
    # Identify constraint field names for fast lookup
    constraint_names: set[str] = set()
    for eq in data.spec.equations:
        if eq.time_derivative_order == 0:
            constraint_names.add(eq.field_name)

    if len(constraint_names) < 2:  # noqa: PLR2004
        return 0.0  # need at least 2 constraints for cross terms

    # Build coordinate arrays lazily (only if position-dependent)
    coord_arrays: dict[str, NDArray[np.float64]] | None = None
    has_posdep = any(
        term.position_dependent
        for eq in data.spec.equations
        if eq.time_derivative_order == 0
        for term in eq.rhs_terms
        if term.field in constraint_names and term.field != eq.field_name
    )
    if has_posdep:
        coord_arrays = _build_coord_arrays(data)

    return _accumulate_cross_constraint_terms(
        data,
        t_idx,
        constraint_names,
        coord_arrays,
    )


def _accumulate_cross_constraint_terms(
    data: SimulationData,
    t_idx: int,
    constraint_names: set[str],
    coord_arrays: dict[str, NDArray[np.float64]] | None,
) -> float:
    """Sum cross-constraint term densities: +½ c_ij ⟨C_i·op(C_j)⟩."""
    energy = 0.0

    for eq in data.spec.equations:
        if eq.time_derivative_order != 0 or eq.field_name not in data.fields:
            continue
        c_i = data.fields[eq.field_name][t_idx]

        for term in eq.rhs_terms:
            if term.field == eq.field_name or term.field not in constraint_names:
                continue
            if _is_velocity_field(term.field):
                continue

            target = _resolve_term_target(data, term.field, t_idx)
            if target is None:
                continue

            coeff = _resolve_coefficient_on_grid(term, data, coord_arrays or {})
            operated = _apply_spatial_operator(
                term.operator,
                target,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )
            energy += 0.5 * float((coeff * c_i * operated).mean())

    return energy


def _evaluate_hamiltonian_factor(
    factor_field: str,
    factor_operator: str,
    data: SimulationData,
    t_idx: int,
) -> NDArray[np.float64] | None:
    """Evaluate a single Hamiltonian factor on the grid.

    For ``time_derivative`` operator, reads the velocity directly from
    ``data.velocities`` (which stores velocities v = dq/dt in the E-L form).

    For spatial operators, applies the operator to the field data.
    For ``identity``, returns the field data directly.

    Returns None if the factor cannot be evaluated (e.g., constraint
    field without stored velocity for time_derivative).
    """
    if factor_operator == "time_derivative":
        vel = data.velocities.get(factor_field)
        if vel is not None:
            return vel[t_idx]
        return None

    # Get the field data
    field_arr = _resolve_term_target(data, factor_field, t_idx)
    if field_arr is None:
        return None

    if factor_operator == "identity":
        return field_arr

    # Apply spatial operator
    return _apply_spatial_operator(
        factor_operator,
        field_arr,
        data.grid_spacing,
        data.periodic,
        bc_types=data.bc_types,
    )


def _gradient_pair_to_second_order(op_a: str, op_b: str) -> str:
    """Map a pair of gradient operators to the equivalent 2nd-order operator.

    Uses integration-by-parts identity (exact for periodic BCs):
    ``⟨∂_a u, ∂_b v⟩ = -⟨u, ∂²_ab v⟩`` where ``∂²_ab`` is laplacian
    (same axis) or cross_derivative (different axes).
    """
    ax_a = _GRADIENT_AXES[op_a]
    ax_b = _GRADIENT_AXES[op_b]
    if ax_a == ax_b:
        return f"laplacian_{'xyz'[ax_a]}"
    lo, hi = sorted([ax_a, ax_b])
    return f"cross_derivative_{'xyz'[lo]}{'xyz'[hi]}"


def _gradient_product_density(  # noqa: PLR0913, PLR0917
    op_a: str,
    field_a: NDArray[np.float64],
    op_b: str,
    field_b: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    bc_types: tuple[str, ...] | None = None,
) -> NDArray[np.float64]:
    """Pointwise density for the gradient inner product ⟨∂_a f, ∂_b g⟩.

    Returns a grid array whose spatial mean gives the gradient inner product.
    Returning an array (not a scalar) allows callers to weight by
    position-dependent coefficients before taking ``.mean()``.

    For **periodic** axes: uses integration by parts (IBP)
        ``density = -f · ∂²_ab(g)``
    matching the solver's 3-point laplacian stencil (exact discrete IBP).

    For **non-periodic** axes: uses direct central-difference
        ``density = ∂_a(f) · ∂_b(g)``
    because discrete IBP has boundary contributions
    (cf. :func:`_gradient_energy_density`).

    This is the **single source of truth** for gradient-product evaluation.
    Both the standalone gradient x gradient Hamiltonian path and the kinetic
    bilinear expansion dispatch here, guaranteeing stencil consistency for
    terms that must cancel (e.g. ``½(∂_x A_0)²`` from kinetic ``- ½(∂_x A_0)²``
    standalone in Proca/CS theories).
    """
    ax_a = _GRADIENT_AXES[op_a]
    bc_a = _effective_bc(ax_a, periodic, bc_types)

    if bc_a == "periodic":
        # IBP: ⟨∂_a f, ∂_b g⟩ = mean(-f · ∂²_ab g)
        second_op = _gradient_pair_to_second_order(op_a, op_b)
        operated = _apply_spatial_operator(
            second_op,
            field_b,
            grid_spacing,
            periodic,
            bc_types=bc_types,
        )
        return -(field_a * operated)

    # Non-periodic: direct gradient product
    grad_a = _apply_spatial_operator(
        op_a,
        field_a,
        grid_spacing,
        periodic,
        bc_types=bc_types,
    )
    grad_b = _apply_spatial_operator(
        op_b,
        field_b,
        grid_spacing,
        periodic,
        bc_types=bc_types,
    )
    return grad_a * grad_b


def _merge_parameters(data: SimulationData) -> dict[str, float]:
    """Merge spec metadata parameters with runtime parameters.

    Spec metadata provides defaults (from theory.toml); runtime parameters
    (from ``--param`` CLI flags) override them.
    """
    import contextlib  # noqa: PLC0415

    raw_meta = data.spec.metadata.get("parameters", {})
    params: dict[str, float] = {}
    if isinstance(raw_meta, dict):
        for k, v in raw_meta.items():  # type: ignore[union-attr]
            with contextlib.suppress(ValueError, TypeError):
                params[str(k)] = float(v)  # type: ignore[arg-type]
    params.update(data.parameters)
    return params


def _compute_hamiltonian_from_canonical(  # noqa: C901, PLR0912, PLR0914
    data: SimulationData,
    t_idx: int,
) -> float:
    """Evaluate the symbolic Hamiltonian from canonical structure.

    For spatial gradient terms, uses :func:`_gradient_product_density` —
    the single source of truth for gradient inner products.  BC-aware:
    periodic → IBP, non-periodic → central-difference.

    For kinetic (``time_derivative x time_derivative``) terms, reads
    velocities directly from ``data.velocities`` (which stores v = dq/dt
    in the E-L velocity form).  No field_rates expansion needed.

    This ensures the measured Hamiltonian uses the **same** finite-difference
    stencils as the solver (3-point laplacian, cascaded-gradient cross
    derivative), making it exactly the conserved quantity of the discrete
    system.  Without IBP, the central-difference gradient squared (5-point
    "wide" stencil) differs from the solver's 3-point laplacian by O(dx²),
    producing spurious time-varying energy drift.

    Parameters
    ----------
    data : SimulationData
    t_idx : int
        Snapshot index.

    Returns
    -------
    float
        Spatially-averaged Hamiltonian energy density.

    Raises
    ------
    ValueError
        If ``data.spec.canonical`` is None.
    """
    canonical = data.spec.canonical
    if canonical is None:
        msg = "_compute_hamiltonian_from_canonical called without canonical structure"
        raise ValueError(msg)

    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
    )

    params = _merge_parameters(data)
    coord_arrays: dict[str, NDArray[np.float64]] | None = None  # lazy-initialized

    # Volume element: sqrt|g_spatial| for curved coordinates.
    # None -> flat spacetime, volume_weight stays 1.0 (no grid allocation).
    volume_weight: float | NDArray[np.float64] = 1.0
    if canonical.volume_element is not None:
        if coord_arrays is None:
            coord_arrays = _build_coord_arrays(data)
        from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

        volume_weight = evaluate_coefficient(
            canonical.volume_element,
            params,
            data.spec.effective_coordinates,
            coord_arrays=coord_arrays,
            t=0.0,
        )

    total = 0.0
    for term in canonical.hamiltonian_terms:
        # --- Coefficient resolution ---
        # Position-dependent coefficients (e.g. Gaussian coupling, Csc[y]/x² in
        # spherical coordinates) must be evaluated on the spatial grid.  The scalar
        # path via _resolve_symbolic_coeff() cannot handle expressions that contain
        # coordinate calls like x[] or y[].
        if term.position_dependent:
            if coord_arrays is None:
                coord_arrays = _build_coord_arrays(data)
            from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

            assert term.coefficient_symbolic is not None  # guaranteed when position_dependent
            coeff: float | NDArray[np.float64] = evaluate_coefficient(
                term.coefficient_symbolic,
                params,
                data.spec.effective_coordinates,
                coord_arrays=coord_arrays,
                t=0.0,
            )
        else:
            coeff = float(term.coefficient)
            if term.coefficient_symbolic is not None and params:
                resolved = _resolve_symbolic_coeff(
                    term.coefficient_symbolic,
                    params,
                )
                if resolved is not None:
                    coeff = float(resolved)

        op_a = term.factor_a.operator
        op_b = term.factor_b.operator

        # Gradient x gradient path: use _gradient_product_density (single source
        # of truth).  BC-aware: periodic→IBP, non-periodic→CD.
        # The helper returns a pointwise density array; coeff (possibly
        # position-dependent NDArray) is multiplied in before .mean().
        if op_a in _GRADIENT_AXES and op_b in _GRADIENT_AXES:
            field_a = _resolve_term_target(data, term.factor_a.field, t_idx)
            field_b = _resolve_term_target(data, term.factor_b.field, t_idx)
            if field_a is None or field_b is None:
                continue
            density = _gradient_product_density(
                op_a,
                field_a,
                op_b,
                field_b,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )
            total += float((coeff * density * volume_weight).mean())
            continue

        # Kinetic: time_derivative x time_derivative — direct velocity lookup.
        # In E-L velocity form, data.velocities stores v = dq/dt directly.
        # No field_rates expansion needed: vel_A = data.velocities[field_A].
        if op_a == "time_derivative" and op_b == "time_derivative":
            fname_a = term.factor_a.field
            fname_b = term.factor_b.field
            vel_a = data.velocities.get(fname_a)
            vel_b = data.velocities.get(fname_b)
            if vel_a is not None and vel_b is not None:
                total += float(
                    (coeff * vel_a[t_idx] * vel_b[t_idx] * volume_weight).mean()
                )
            continue

        # All other terms: identity, mixed operator x identity, etc.
        fa = _evaluate_hamiltonian_factor(
            term.factor_a.field,
            term.factor_a.operator,
            data,
            t_idx,
        )
        fb = _evaluate_hamiltonian_factor(
            term.factor_b.field,
            term.factor_b.operator,
            data,
            t_idx,
        )
        if fa is None or fb is None:
            continue

        total += float((coeff * fa * fb * volume_weight).mean())

    return total


def compute_system_energy(  # noqa: PLR0914
    data: SimulationData,
    t_idx: int,
) -> SystemEnergy:
    """Compute Hamiltonian energy density at snapshot *t_idx*.

    When the spec includes canonical structure (from ``tidal derive``),
    evaluates the Legendre-transform Hamiltonian directly from structured
    quadratic terms. Otherwise, falls back to the virial-based formula.

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
    # or we need the volume element for curved coordinates.
    coord_arrays: dict[str, NDArray[np.float64]] | None = None
    has_posdep_mass = any(
        term.operator == "identity"
        and term.field == eq.field_name
        and term.position_dependent
        for eq in data.spec.equations
        for term in eq.rhs_terms
    )
    needs_coord_arrays = has_posdep_mass or (
        data.spec.canonical is not None
        and data.spec.canonical.volume_element is not None
    )
    if needs_coord_arrays:
        coord_arrays = _build_coord_arrays(data)

    # Volume element: sqrt|g_spatial| for curved coordinates.
    # 1.0 for flat spacetimes (no grid allocation, scalar multiply is no-op).
    volume_weight: float | NDArray[np.float64] = 1.0
    if (
        data.spec.canonical is not None
        and data.spec.canonical.volume_element is not None
        and coord_arrays is not None
    ):
        from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

        volume_weight = evaluate_coefficient(
            data.spec.canonical.volume_element,
            _merge_parameters(data),
            data.spec.effective_coordinates,
            coord_arrays=coord_arrays,
            t=0.0,
        )

    per_field: dict[str, FieldEnergy] = {}
    for field_idx, eq in enumerate(data.spec.equations):
        name = eq.field_name
        if eq.time_derivative_order == 0:
            continue

        field_snapshot = data.fields[name][t_idx]
        vel_snapshot = data.velocities.get(name)
        vel_arr = vel_snapshot[t_idx] if vel_snapshot is not None else None

        m2 = _resolve_mass_squared(data, field_idx, coord_arrays=coord_arrays)
        axes = _self_gradient_axes(eq)
        per_field[name] = compute_field_energy(
            field_snapshot,
            vel_arr,
            m2,
            data.grid_spacing,
            data.periodic,
            gradient_axes=axes,
            bc_types=data.bc_types,
            volume_weight=volume_weight,
        )

    # Use canonical Hamiltonian when available (Phase K: Legendre transform)
    if data.spec.canonical is not None:
        total = _compute_hamiltonian_from_canonical(data, t_idx)
        self_sum = sum(fe.total for fe in per_field.values())
        interaction = total - self_sum
        return SystemEnergy(per_field=per_field, interaction=interaction, total=total)

    # Fallback: virial-based formula for legacy specs without canonical structure.
    v_virial = _compute_virial_potential(data, t_idx)
    v_constraint = _compute_constraint_self_energy(data, t_idx)
    v_constraint_cross = _compute_constraint_coupling_energy(data, t_idx)
    v_total = v_virial + v_constraint + v_constraint_cross
    self_potential = sum(fe.gradient + fe.mass for fe in per_field.values())
    interaction = v_total - self_potential
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
    """Compute energy density for every snapshot in the simulation.

    Returns
    -------
    times : ndarray, shape ``(n_snapshots,)``
    per_field : dict[str, ndarray]
        Each value is shape ``(n_snapshots,)`` — energy density of that field.
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
