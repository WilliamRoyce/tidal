"""Solver-agnostic validation for TIDAL equation specifications.

Extracted from ``pde_builder.py`` so that IDA, leapfrog, and any future
solver can share the same validation logic without importing py-pde.

All functions are module-level (no shared state) and raise ``ValueError``
on invalid specs per the fail-fast-and-loud convention.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from tidal.solver.operators import operator_min_dim

if TYPE_CHECKING:
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


def validate_operator_dimensions(spec: EquationSystem) -> None:
    """Check that all operators are compatible with the spec's spatial dimension.

    Raises
    ------
    ValueError
        If any operator requires more dimensions than ``spec.spatial_dimension``.
    """
    spatial_dim = spec.spatial_dimension
    for eq in spec.equations:
        for term in eq.rhs_terms:
            try:
                min_dim = operator_min_dim(term.operator)
            except ValueError:
                continue  # Unknown ops will fail at apply_operator time
            if min_dim > spatial_dim:
                msg = (
                    f"Operator '{term.operator}' in equation for "
                    f"'{eq.field_name}' requires at least {min_dim}D "
                    f"spatial grid, but the spec has "
                    f"spatial_dimension={spatial_dim} "
                    f"(from {spec.dimension}D spacetime)."
                )
                raise ValueError(msg)


def validate_field_references(spec: EquationSystem) -> None:
    """Check that all term field references point to valid fields.

    Also checks canonical field_rates and spatial_momenta references.

    Raises
    ------
    ValueError
        If a field reference is invalid.
    """
    valid_fields = set(spec.component_names)
    # Also accept momentum names (pi_N format)
    valid_fields.update(f"pi_{i}" for i in range(spec.n_components))

    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field not in valid_fields:
                msg = (
                    f"Unknown field reference '{term.field}' "
                    f"in equation for '{eq.field_name}'. "
                    f"Valid fields: {sorted(valid_fields)}."
                )
                raise ValueError(msg)

    # Check canonical references
    if spec.canonical is not None:
        _validate_canonical_refs(spec, valid_fields)


def _validate_canonical_refs(spec: EquationSystem, valid_fields: set[str]) -> None:
    """Check canonical field_rates and spatial_momenta references.

    Raises
    ------
    ValueError
        If a canonical reference is invalid.
    """
    canonical = spec.canonical
    assert canonical is not None

    for field_name, terms in canonical.field_rates.items():
        for term in terms:
            if term.field not in valid_fields and not term.field.startswith("pi_"):
                msg = (
                    f"Unknown field reference '{term.field}' "
                    f"in field_rates for '{field_name}'."
                )
                raise ValueError(msg)

    if canonical.spatial_momenta is not None:
        for field_name, terms in canonical.spatial_momenta.items():
            for term in terms:
                if term.field not in valid_fields:
                    msg = (
                        f"Unknown field reference '{term.field}' "
                        f"in spatial_momenta for '{field_name}'."
                    )
                    raise ValueError(msg)


def check_cfl_stability(
    spec: EquationSystem,
    grid: GridInfo,
    dt: float,
) -> list[str]:
    """Check CFL stability condition for explicit time-steppers.

    Returns a list of warning strings (empty if all clear).
    The CFL condition for the wave equation is dt <= dx / c where
    c is the maximum wave speed (estimated from the laplacian coefficient).
    """
    warnings: list[str] = []
    dx_min = min(grid.dx)

    for eq in spec.equations:
        if eq.time_derivative_order < 2:  # noqa: PLR2004
            continue

        # Find the largest laplacian coefficient
        max_lap_coeff = 0.0
        for term in eq.rhs_terms:
            if "laplacian" in term.operator:
                max_lap_coeff = max(max_lap_coeff, abs(term.coefficient))

        if max_lap_coeff > 0:
            import math  # noqa: PLC0415

            c_est = math.sqrt(max_lap_coeff)
            cfl_dt = dx_min / c_est
            if dt > cfl_dt:
                warnings.append(
                    f"CFL violation for '{eq.field_name}': "
                    f"dt={dt:.4g} > CFL limit={cfl_dt:.4g} "
                    f"(c_est={c_est:.4g}, dx_min={dx_min:.4g}). "
                    f"Consider reducing dt or increasing grid resolution."
                )
    return warnings


def check_mass_sign(
    coeff_eval: CoefficientEvaluator,
    spec: EquationSystem,
) -> list[str]:
    """Check for sign-changing position-dependent mass terms.

    Returns a list of warning strings for tachyonic diagnostics.
    """
    import numpy as np  # noqa: PLC0415

    warnings: list[str] = []
    for eq_idx, eq in enumerate(spec.equations):
        for term_idx, term in enumerate(eq.rhs_terms):
            if (
                term.operator != "identity"
                or term.field != eq.field_name
                or term.coefficient_symbolic is None
                or not term.position_dependent
                or term.time_dependent
            ):
                continue

            result = coeff_eval.resolve(term, t=0.0, eq_idx=eq_idx, term_idx=term_idx)
            if (
                isinstance(result, np.ndarray)
                and float(result.min()) * float(result.max()) < 0
            ):
                warnings.append(
                    f"Position-dependent mass term "
                    f"'{term.coefficient_symbolic}' for field "
                    f"'{eq.field_name}' changes sign across "
                    f"the grid (min={float(result.min()):.4g}, "
                    f"max={float(result.max()):.4g})."
                )
    return warnings


def check_robin_stability(grid: GridInfo) -> list[str]:
    """Check Robin BC ghost-cell formula stability.

    The ghost-cell formula denominator is ``gamma * dx + 2``.  When
    ``gamma * dx >= 2`` the mirror factor ``(2 - gamma*dx)/(gamma*dx + 2)``
    becomes non-positive, which can destabilize the scheme.

    Returns a list of warning strings (empty if all clear).
    """
    warnings: list[str] = []
    if grid.axis_bcs is None:
        return warnings

    for i, abc in enumerate(grid.axis_bcs):
        if abc.periodic:
            continue
        dx = grid.dx[i]
        for side_label, side in [("low", abc.low), ("high", abc.high)]:
            if side is None or side.kind != "robin":
                continue
            if side.gamma * dx >= 2.0:  # noqa: PLR2004
                warnings.append(
                    f"Robin BC on axis {i} ({side_label}): "
                    f"gamma*dx = {side.gamma * dx:.4g} >= 2. "
                    f"Ghost-cell formula becomes unstable. "
                    f"Increase grid resolution or decrease gamma."
                )
    return warnings
