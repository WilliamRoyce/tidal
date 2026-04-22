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
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.measurement._io import SimulationData
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        HamiltonianTerm,
        OperatorTerm,
    )

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
# Derivative helpers — delegated to operators.py
# ------------------------------------------------------------------
#
# Energy measurement MUST use the same finite-difference stencils as the
# PDE solver (operators.py) to track the discrete conserved Hamiltonian
# exactly.  Previous versions reimplemented stencils here independently;
# this caused inconsistencies when the solver FD order changed.
#
# Now all spatial operators are called through the solver module, which
# automatically inherits the current FD order from set_fd_order().
#
# Reference: the stencil-matching requirement is discussed in
# Hairer, Lubich & Wanner, "Geometric Numerical Integration", Springer,
# 2006, Chapter VI — virial energy tracks the discrete (shadow)
# Hamiltonian only when the measurement stencils match the solver's.


def _make_grid_info(
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    shape: tuple[int, ...],
    bc_types: tuple[str, ...] | None = None,
) -> object:
    """Build a minimal ``GridInfo`` from measurement parameters.

    Lazily imports ``GridInfo`` to avoid circular imports at module load.
    The ``bounds`` are reconstructed from ``shape * dx`` per axis.
    """
    from tidal.solver.grid import GridInfo  # noqa: PLC0415

    bounds = tuple(
        (0.0, float(n) * dx) for n, dx in zip(shape, grid_spacing, strict=True)
    )
    bc: tuple[str, ...] | None = bc_types
    return GridInfo(
        bounds=bounds,
        shape=shape,
        periodic=periodic,
        bc=bc,
    )


def _resolve_bc_spec(
    bc_types: tuple[str, ...] | None,
    periodic: tuple[bool, ...],
) -> str | tuple[str, ...]:
    """Build a BC spec for operators.py from energy module parameters."""
    if bc_types is not None:
        return bc_types
    return tuple("periodic" if p else "dirichlet" for p in periodic)


# ------------------------------------------------------------------
# Operator-aware gradient axes
# ------------------------------------------------------------------


def _self_gradient_axes(eq: ComponentEquation) -> list[int] | None:  # pyright: ignore[reportUnusedFunction]
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
# Spatial operator application (for Hamiltonian evaluation)
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


def _apply_spatial_operator(  # noqa: PLR0913, PLR0917
    operator: str,
    field: NDArray[np.float64],
    grid_spacing: tuple[float, ...],
    periodic: tuple[bool, ...],
    bc_types: tuple[str, ...] | None = None,
    grid: object | None = None,
) -> NDArray[np.float64]:
    """Apply a named spatial operator to a field array.

    Delegates to ``operators.apply_operator()`` from the solver module,
    ensuring the same FD stencils are used for energy measurement and
    PDE integration.  This is critical for virial energy tracking the
    discrete (shadow) Hamiltonian exactly.

    Parameters
    ----------
    bc_types : tuple[str, ...] | None
        Per-axis BC type.  When ``None``, falls back to ``periodic`` booleans.
    grid : GridInfo | None
        Pre-built grid info.  When provided, skips ``_make_grid_info()``
        reconstruction (avoids redundant allocation per call).
    """
    from tidal.solver.operators import apply_operator as _solver_apply  # noqa: PLC0415

    if operator == "identity":
        return field

    bc_spec = _resolve_bc_spec(bc_types, periodic)
    if grid is None:
        grid = _make_grid_info(grid_spacing, periodic, field.shape, bc_types)
    return _solver_apply(operator, field, grid, bc_spec)  # type: ignore[arg-type]


# ------------------------------------------------------------------
# Term resolution helpers
# ------------------------------------------------------------------


def _is_velocity_field(field_name: str) -> bool:
    """Check if a field name is a velocity reference (v_field_name, e.g. v_A_1)."""
    return field_name.startswith("v_") and len(field_name) > len("v_")


# ------------------------------------------------------------------
# Mixed time-space operator support
# ------------------------------------------------------------------
#
# Mixed operators arise from Lagrangians containing second-order time
# derivatives (e.g. linearized Einstein-Hilbert: L^(2) ~ h * d²_t h).
# The Legendre transform doesn't eliminate these, so the Hamiltonian
# retains operators like mixed_1_0_0_1 (d_t d_z) or mixed_2_0_0_0 (d²_t).
#
# Strategy: decompose into time part (velocity/EOM) + spatial part.


def _parse_mixed_operator(operator: str) -> tuple[int, tuple[int, ...]]:
    """Parse ``mixed_T_S1_S2_...`` → ``(time_order, spatial_orders)``.

    Examples
    --------
    >>> _parse_mixed_operator("mixed_1_0_0_1")
    (1, (0, 0, 1))
    >>> _parse_mixed_operator("mixed_2_0_0")
    (2, (0, 0))
    """
    parts = operator.split("_")[1:]  # strip "mixed" prefix
    return int(parts[0]), tuple(int(p) for p in parts[1:])


def _compute_acceleration_from_eom(
    data: SimulationData,
    field: str,
    t_idx: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Compute d²_t(field) by evaluating the E-L equation RHS.

    For dynamical fields whose equation has the form ``d²_t f = Σ c_k op_k(g_k)``,
    evaluates the RHS using stored field/velocity snapshots and spatial operators.
    Returns ``None`` if no matching dynamical equation exists.
    """
    for eq in data.spec.equations:
        if eq.field_name != field or eq.time_derivative_order != 2:  # noqa: PLR2004
            continue

        shape: tuple[int, ...] = next(iter(data.fields.values()))[t_idx].shape
        result: NDArray[np.float64] = np.zeros(shape)

        for term in eq.rhs_terms:
            if term.operator == "first_derivative_t":
                vel = data.velocities.get(term.field)
                if vel is None:
                    continue
                operated: NDArray[np.float64] = vel[t_idx]
            else:
                target = _resolve_term_target(data, term.field, t_idx)
                if target is None:
                    continue
                operated = _apply_spatial_operator(
                    term.operator,
                    target,
                    data.grid_spacing,
                    data.periodic,
                    bc_types=data.bc_types,
                )

            coeff = _resolve_term_coefficient(term, params)
            result += coeff * operated

        return result

    return None


def _differentiate_constraint(
    data: SimulationData,
    field: str,
    t_idx: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Compute d_t of a constraint field by differentiating its equation.

    For a constraint ``f = Σ c_k op_k(source_k)``, the time derivative is
    ``d_t(f) = Σ c_k op_k(d_t(source_k))`` where ``d_t(source)`` is the
    velocity (for field sources) or acceleration from EOM (for velocity sources).
    """
    for eq in data.spec.equations:
        if eq.field_name != field or eq.time_derivative_order != 0:
            continue

        shape: tuple[int, ...] = next(iter(data.fields.values()))[t_idx].shape
        result: NDArray[np.float64] = np.zeros(shape)

        for term in eq.rhs_terms:
            # Skip self-referential terms (e.g. trace constraint h_9 = h_4 + h_7 + h_9)
            if term.field == field:
                continue

            if _is_velocity_field(term.field):
                # Source is velocity v_X → d_t(v_X) = acceleration from EOM
                base_field = term.field[2:]  # strip "v_" prefix
                dt_source = _compute_acceleration_from_eom(
                    data,
                    base_field,
                    t_idx,
                    params,
                )
            else:
                # Source is field X → d_t(X) = velocity
                vel = data.velocities.get(term.field)
                dt_source = vel[t_idx] if vel is not None else None

            if dt_source is None:
                continue

            operated: NDArray[np.float64] = _apply_spatial_operator(
                term.operator,
                dt_source,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )

            coeff = _resolve_term_coefficient(term, params)
            result += coeff * operated

        return result

    return None


def _differentiate_constraint_twice(
    data: SimulationData,
    field: str,
    t_idx: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Compute d²_t of a constraint field.

    Handles constraints whose sources are plain fields (not velocities):
    ``d²_t(field_source) = acceleration`` from EOM.  Constraints with
    velocity sources would need the jerk (d³_t), which is not supported.
    """
    for eq in data.spec.equations:
        if eq.field_name != field or eq.time_derivative_order != 0:
            continue

        shape: tuple[int, ...] = next(iter(data.fields.values()))[t_idx].shape
        result: NDArray[np.float64] = np.zeros(shape)

        for term in eq.rhs_terms:
            if term.field == field:
                continue

            if _is_velocity_field(term.field):
                # d²_t(velocity) = jerk — not supported
                return None

            accel = _compute_acceleration_from_eom(data, term.field, t_idx, params)
            if accel is None:
                return None

            operated: NDArray[np.float64] = _apply_spatial_operator(
                term.operator,
                accel,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )

            coeff = _resolve_term_coefficient(term, params)
            result += coeff * operated

        return result

    return None


def _resolve_time_derivative(
    data: SimulationData,
    field: str,
    t_idx: int,
    time_order: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Resolve the Nth time derivative of a field from stored data + EOM.

    ========== ===============================================================
    time_order Resolution
    ========== ===============================================================
    0          Field value (identity)
    1          Velocity (stored), or differentiated constraint equation
    2          Acceleration from EOM RHS, or double-differentiated constraint
    ≥ 3        Not supported (returns None)
    ========== ===============================================================
    """
    if time_order == 0:
        return _resolve_term_target(data, field, t_idx)

    if time_order == 1:
        vel = data.velocities.get(field)
        if vel is not None:
            return vel[t_idx]
        return _differentiate_constraint(data, field, t_idx, params)

    if time_order == 2:  # noqa: PLR2004
        accel = _compute_acceleration_from_eom(data, field, t_idx, params)
        if accel is not None:
            return accel
        return _differentiate_constraint_twice(data, field, t_idx, params)

    # time_order >= 3: not currently supported
    return None


_SPATIAL_AXIS_LETTERS: tuple[str, ...] = ("x", "y", "z")


def _evaluate_mixed_factor(
    factor_field: str,
    factor_operator: str,
    data: SimulationData,
    t_idx: int,
    params: dict[str, float],
) -> NDArray[np.float64] | None:
    """Evaluate a ``mixed_T_S1_S2_...`` operator on a field.

    Decomposition: apply T time derivatives (velocity/EOM), then apply
    spatial derivatives sequentially.
    """
    time_order, spatial_orders = _parse_mixed_operator(factor_operator)

    # Time part
    base = _resolve_time_derivative(data, factor_field, t_idx, time_order, params)
    if base is None:
        return None

    # Spatial part: apply gradient operators for each nonzero spatial order.
    # Handle dimension mismatch from plane-wave reduction: when the mixed
    # operator has more spatial indices than grid dimensions, collapse all
    # derivative orders onto the available axes (axis 0 for 1D).
    n_spatial = len(data.grid_spacing)
    if len(spatial_orders) > n_spatial:
        total_order = sum(spatial_orders)
        effective_orders: tuple[int, ...] = (total_order,) + (0,) * (n_spatial - 1)
    else:
        effective_orders = spatial_orders

    for ax_idx, order in enumerate(effective_orders):
        for _ in range(order):
            ax_letter = _SPATIAL_AXIS_LETTERS[ax_idx]
            base = _apply_spatial_operator(
                f"gradient_{ax_letter}",
                base,
                data.grid_spacing,
                data.periodic,
                bc_types=data.bc_types,
            )

    return base


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


def _resolve_coefficient_on_grid(  # pyright: ignore[reportUnusedFunction]
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
# Mass resolution
# ------------------------------------------------------------------


def _resolve_mass_squared(  # pyright: ignore[reportUnusedFunction]
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


# ------------------------------------------------------------------
# Hamiltonian evaluation helpers
# ------------------------------------------------------------------


def _evaluate_hamiltonian_factor(
    factor_field: str,
    factor_operator: str,
    data: SimulationData,
    t_idx: int,
    grid: object | None = None,
) -> NDArray[np.float64] | None:
    """Evaluate a single Hamiltonian factor on the grid.

    For ``time_derivative`` operator, reads the velocity directly from
    ``data.velocities`` (which stores velocities v = dq/dt in the E-L form).

    For ``mixed_*`` operators (time-space cross derivatives), decomposes into
    time derivative resolution (velocity/EOM) + spatial operator application.

    For spatial operators, applies the operator to the field data.
    For ``identity``, returns the field data directly.

    Returns None if the factor cannot be evaluated (e.g., constraint
    field without stored velocity for time_derivative).

    Parameters
    ----------
    grid : GridInfo | None
        Pre-built grid info from context (avoids per-call reconstruction).
    """
    if factor_operator == "time_derivative":
        vel = data.velocities.get(factor_field)
        if vel is not None:
            return vel[t_idx]
        return None

    # Mixed time-space operators (e.g. mixed_1_0_0_1 = d_t d_z)
    if factor_operator.startswith("mixed_"):
        params = _merge_parameters(data)
        return _evaluate_mixed_factor(
            factor_field,
            factor_operator,
            data,
            t_idx,
            params,
        )

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
        grid=grid,
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
    grid: object | None = None,
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
    because discrete IBP requires a constant integration weight and breaks
    down when the energy integral includes a position-dependent volume
    element (e.g. ``r²`` in spherical coordinates).

    This is the **single source of truth** for gradient-product evaluation.
    Both the standalone gradient x gradient Hamiltonian path and the kinetic
    bilinear expansion dispatch here, guaranteeing stencil consistency for
    terms that must cancel (e.g. ``½(∂_x A_0)²`` from kinetic ``- ½(∂_x A_0)²``
    standalone in Proca/CS theories).
    """
    ax_a = _GRADIENT_AXES[op_a]
    bc_a = _effective_bc(ax_a, periodic, bc_types)

    from tidal.solver.operators import get_spectral as _get_spectral  # noqa: PLC0415

    if bc_a == "periodic" and not _get_spectral():
        # FD periodic: IBP is more accurate than two separate FD gradients
        # because it uses a single laplacian stencil (fewer FD applications).
        # IBP: ⟨∂_a f, ∂_b g⟩ = mean(-f · ∂²_ab g)
        second_op = _gradient_pair_to_second_order(op_a, op_b)
        operated = _apply_spatial_operator(
            second_op,
            field_b,
            grid_spacing,
            periodic,
            bc_types=bc_types,
            grid=grid,
        )
        return -(field_a * operated)

    # Spectral or non-periodic: direct gradient product.
    # For spectral operators, direct gradient*gradient is more accurate
    # than IBP (laplacian) because rfft Nyquist-mode handling introduces
    # O(1/N) discrepancy between mean(|∂f|²) and -mean(f·∂²f).
    grad_a = _apply_spatial_operator(
        op_a,
        field_a,
        grid_spacing,
        periodic,
        bc_types=bc_types,
        grid=grid,
    )
    grad_b = _apply_spatial_operator(
        op_b,
        field_b,
        grid_spacing,
        periodic,
        bc_types=bc_types,
        grid=grid,
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


@dataclass
class _HamiltonianContext:
    """Pre-computed snapshot-invariant data for Hamiltonian evaluation.

    Hoisting these computations out of the per-snapshot loop avoids
    redundant parameter merging, coordinate array construction, volume
    element evaluation, per-term coefficient resolution, and GridInfo
    reconstruction.
    """

    params: dict[str, float]
    coord_arrays: dict[str, NDArray[np.float64]] | None
    volume_weight: float | NDArray[np.float64]
    term_coeffs: list[float | NDArray[np.float64]]
    resolve_symbolic: Callable[..., float | None]
    grid_info: object | None = None  # GridInfo, typed as object to avoid import


def _prepare_hamiltonian_context(data: SimulationData) -> _HamiltonianContext:
    """Build a reusable context for repeated Hamiltonian evaluations."""
    canonical = data.spec.canonical
    assert canonical is not None

    from tidal.symbolic.json_loader import (  # noqa: PLC0415
        _resolve_symbolic_coeff,  # pyright: ignore[reportPrivateUsage]
    )

    params = _merge_parameters(data)
    coord_arrays: dict[str, NDArray[np.float64]] | None = None

    # Volume element
    volume_weight: float | NDArray[np.float64] = 1.0
    if canonical.volume_element is not None:
        coord_arrays = _build_coord_arrays(data)
        from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

        volume_weight = evaluate_coefficient(
            canonical.volume_element,
            params,
            data.spec.effective_coordinates,
            coord_arrays=coord_arrays,
            t=0.0,
        )

    # Pre-resolve all term coefficients
    term_coeffs: list[float | NDArray[np.float64]] = []
    for term in canonical.hamiltonian_terms:
        if term.position_dependent:
            if coord_arrays is None:
                coord_arrays = _build_coord_arrays(data)
            from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

            assert term.coefficient_symbolic is not None
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
        term_coeffs.append(coeff)

    # Pre-build GridInfo once (avoids reconstruction per operator call)
    grid_info = None
    if data.fields:
        first_field = next(iter(data.fields.values()))
        field_shape = first_field[0].shape  # spatial shape from first snapshot
        grid_info = _make_grid_info(
            data.grid_spacing,
            data.periodic,
            field_shape,
            data.bc_types,
        )

    return _HamiltonianContext(
        params=params,
        coord_arrays=coord_arrays,
        volume_weight=volume_weight,
        term_coeffs=term_coeffs,
        resolve_symbolic=_resolve_symbolic_coeff,
        grid_info=grid_info,
    )


def _compute_hamiltonian_from_canonical(  # pyright: ignore[reportUnusedFunction]
    data: SimulationData,
    t_idx: int,
    ctx: _HamiltonianContext | None = None,
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

    # Use pre-computed context when available (timeseries path),
    # otherwise compute inline (single-snapshot callers).
    if ctx is None:
        ctx = _prepare_hamiltonian_context(data)

    volume_weight = ctx.volume_weight

    total = 0.0
    grid = ctx.grid_info
    for term_idx, term in enumerate(canonical.hamiltonian_terms):
        contrib = _evaluate_single_hamiltonian_term(
            term,
            ctx.term_coeffs[term_idx],
            volume_weight,
            data,
            t_idx,
            grid=grid,
        )
        total += contrib

    return total


def _evaluate_single_hamiltonian_term(  # noqa: PLR0913, PLR0917
    term: HamiltonianTerm,
    coeff: float | NDArray[np.float64],
    volume_weight: float | NDArray[np.float64],
    data: SimulationData,
    t_idx: int,
    grid: object | None = None,
) -> float:
    """Evaluate a single Hamiltonian term's contribution to energy density.

    Returns the spatially-averaged contribution ``⟨coeff * f_a * f_b * √|g|⟩``.

    Parameters
    ----------
    grid : GridInfo | None
        Pre-built grid info from context (avoids per-call reconstruction).
    """
    op_a = term.factor_a.operator
    op_b = term.factor_b.operator

    # Gradient x gradient path: use _gradient_product_density (single source
    # of truth).  BC-aware: periodic→IBP, non-periodic→CD.
    if op_a in _GRADIENT_AXES and op_b in _GRADIENT_AXES:
        field_a = _resolve_term_target(data, term.factor_a.field, t_idx)
        field_b = _resolve_term_target(data, term.factor_b.field, t_idx)
        if field_a is None or field_b is None:
            return 0.0
        density = _gradient_product_density(
            op_a,
            field_a,
            op_b,
            field_b,
            data.grid_spacing,
            data.periodic,
            bc_types=data.bc_types,
            grid=grid,
        )
        return float((coeff * density * volume_weight).mean())

    # Kinetic: time_derivative x time_derivative — direct velocity lookup.
    if op_a == "time_derivative" and op_b == "time_derivative":
        fname_a = term.factor_a.field
        fname_b = term.factor_b.field
        vel_a = data.velocities.get(fname_a)
        vel_b = data.velocities.get(fname_b)
        if vel_a is not None and vel_b is not None:
            return float((coeff * vel_a[t_idx] * vel_b[t_idx] * volume_weight).mean())
        return 0.0

    # All other terms: identity, mixed operator x identity, etc.
    fa = _evaluate_hamiltonian_factor(
        term.factor_a.field,
        term.factor_a.operator,
        data,
        t_idx,
        grid=grid,
    )
    fb = _evaluate_hamiltonian_factor(
        term.factor_b.field,
        term.factor_b.operator,
        data,
        t_idx,
        grid=grid,
    )
    if fa is None or fb is None:
        return 0.0

    return float((coeff * fa * fb * volume_weight).mean())


def _compute_hamiltonian_per_field(
    data: SimulationData,
    t_idx: int,
    ctx: _HamiltonianContext | None = None,
    order_filter: int | None = None,
) -> tuple[dict[str, float], float]:
    """Decompose Hamiltonian into per-field self-energy and interaction.

    Self-energy: terms where both factors reference the same base field.
    Interaction: terms where factors reference different fields.

    Parameters
    ----------
    data : SimulationData
    t_idx : int
        Snapshot index.
    ctx : _HamiltonianContext or None
        Pre-computed context (for timeseries path).
    order_filter : int or None
        When set, only include Hamiltonian terms whose ``order_in_eps``
        equals this value. ``None`` includes all terms (default).

    Returns
    -------
    per_field : dict[str, float]
        Per-field self-energy density (base field name → energy).
    interaction : float
        Total interaction energy density (cross-field terms).

    Raises
    ------
    ValueError
        If the spec lacks hamiltonian_terms.
    """
    canonical = data.spec.canonical
    if canonical is None or not canonical.hamiltonian_terms:
        msg = (
            "No hamiltonian_terms in JSON spec. Re-derive the theory to populate "
            "the canonical Hamiltonian (required for correct energy measurement). "
            "Run: uv run tidal derive <theory.toml>"
        )
        raise ValueError(msg)

    if ctx is None:
        ctx = _prepare_hamiltonian_context(data)

    volume_weight = ctx.volume_weight
    grid = ctx.grid_info
    per_field: dict[str, float] = {}
    interaction = 0.0

    for term_idx, term in enumerate(canonical.hamiltonian_terms):
        if order_filter is not None and term.order_in_eps != order_filter:
            continue
        contrib = _evaluate_single_hamiltonian_term(
            term,
            ctx.term_coeffs[term_idx],
            volume_weight,
            data,
            t_idx,
            grid=grid,
        )
        if term.is_self_energy:
            base = term.base_field_a
            per_field[base] = per_field.get(base, 0.0) + contrib
        else:
            interaction += contrib

    return per_field, interaction


def compute_system_energy(
    data: SimulationData,
    t_idx: int,
    _ctx: _HamiltonianContext | None = None,
    order: int | None = None,
) -> SystemEnergy:
    """Compute Hamiltonian energy density at snapshot *t_idx*.

    Uses the Legendre-transform Hamiltonian from canonical structure,
    decomposed into per-field self-energy and cross-field interaction.
    This is the only physically correct energy computation — the
    Hamiltonian coefficients include all metric/kinetic prefactors
    from the covariant Lagrangian, making it coordinate-invariant.

    Parameters
    ----------
    data : SimulationData
    t_idx : int
        Snapshot index.
    order : int or None
        When set, restrict to Hamiltonian terms with ``order_in_eps``
        equal to this value. ``0`` gives the leading-order Maxwell part;
        ``1`` gives only the first EH/perturbative correction. ``None``
        (default) sums all terms.

    Raises
    ------
    ValueError
        If *t_idx* is out of range, or if the spec lacks hamiltonian_terms
        (re-derive with ``tidal derive`` to populate them).
    """
    if t_idx < 0 or t_idx >= data.n_snapshots:
        msg = f"t_idx={t_idx} out of range [0, {data.n_snapshots})"
        raise ValueError(msg)

    if data.spec.canonical is None or not data.spec.canonical.hamiltonian_terms:
        msg = (
            "No hamiltonian_terms in JSON spec. Re-derive the theory to "
            "populate the canonical Hamiltonian (required for correct energy "
            "measurement). Run: uv run tidal derive <theory.toml>"
        )
        raise ValueError(msg)

    if data.spec.has_constraint_velocity_terms:
        import warnings  # noqa: PLC0415

        warnings.warn(
            "Hamiltonian contains kinetic coupling between constraint fields "
            "(time_derivative_order=0). The naive Legendre transform "
            "H = sum(pi_i * v_i) - L is not the correct conserved quantity "
            "for theories with non-trivial constraints (Dirac-Bergmann "
            "theory). Energy measurements may show systematic drift. "
            "See: https://github.com/WilliamRoyce/torsion-gertsenshtein/issues/178",
            stacklevel=2,
        )

    per_field_totals, interaction = _compute_hamiltonian_per_field(
        data,
        t_idx,
        ctx=_ctx,
        order_filter=order,
    )
    per_field = {
        name: FieldEnergy(kinetic=0.0, gradient=0.0, mass=0.0, total=total)
        for name, total in per_field_totals.items()
    }
    total = sum(per_field_totals.values()) + interaction
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

    # Pre-compute Hamiltonian context once (avoids re-merging parameters,
    # rebuilding coord arrays, and re-resolving coefficients per snapshot).
    ctx: _HamiltonianContext | None = None
    if data.spec.canonical is not None:
        ctx = _prepare_hamiltonian_context(data)

    # Compute first snapshot to discover field names, then pre-allocate.
    se0 = compute_system_energy(data, 0, _ctx=ctx)
    field_names = list(se0.per_field.keys())
    per_field_np: dict[str, NDArray[np.float64]] = {
        name: np.empty(n, dtype=np.float64) for name in field_names
    }
    interaction = np.empty(n, dtype=np.float64)
    total = np.empty(n, dtype=np.float64)

    # Fill snapshot 0
    for name, fe in se0.per_field.items():
        per_field_np[name][0] = fe.total
    interaction[0] = se0.interaction
    total[0] = se0.total

    # Fill remaining snapshots
    for t_idx in range(1, n):
        se = compute_system_energy(data, t_idx, _ctx=ctx)
        for name, fe in se.per_field.items():
            per_field_np[name][t_idx] = fe.total
        interaction[t_idx] = se.interaction
        total[t_idx] = se.total

    return (
        data.times.copy(),
        per_field_np,
        interaction,
        total,
    )


# ------------------------------------------------------------------
# Validation
# ------------------------------------------------------------------


def _validate_array(arr: NDArray[np.float64], label: str) -> None:  # pyright: ignore[reportUnusedFunction]
    """Check array for non-finite values.

    Raises
    ------
    ValueError
        If *arr* contains NaN or Inf.
    """
    if not np.isfinite(arr).all():
        msg = f"{label} contains NaN or Inf values"
        raise ValueError(msg)
