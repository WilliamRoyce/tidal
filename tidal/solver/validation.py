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
    import numpy as np
    from numpy.typing import NDArray

    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


# ---------------------------------------------------------------------------
# Dimension checks
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Field reference checks
# ---------------------------------------------------------------------------


def validate_field_references(spec: EquationSystem) -> None:
    """Check that all term field references point to valid fields.

    Raises
    ------
    ValueError
        If a field reference is invalid.
    """
    valid_fields = set(spec.component_names)
    # Accept velocity names in v_field_name format (e.g. v_A_1)
    valid_fields.update(f"v_{eq.field_name}" for eq in spec.equations)

    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.field not in valid_fields:
                msg = (
                    f"Unknown field reference '{term.field}' "
                    f"in equation for '{eq.field_name}'. "
                    f"Valid fields: {sorted(valid_fields)}."
                )
                raise ValueError(msg)


# ---------------------------------------------------------------------------
# CFL / mass diagnostics
# ---------------------------------------------------------------------------


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
                    f"Consider reducing dt or increasing grid resolution.",
                )
    return warnings


def _summed_self_mass_profile(
    coeff_eval: CoefficientEvaluator,
    spec: EquationSystem,
    eq_idx: int,
) -> tuple[NDArray[np.float64] | None, str]:
    """Return an equation's total position-dependent self-mass profile.

    A component can carry several ``identity`` terms acting on its own field —
    in the dual-Gaussian backgrounds each carries three — and the physical mass
    is their **sum**.  Testing them one at a time can warn about a term that
    crosses zero while the sum never does, or stay silent when each term is
    single-signed but the sum is not.

    Returns ``(None, "")`` when the equation has no position-dependent
    self-mass, otherwise the summed profile and a label naming the terms.
    """
    import numpy as np  # noqa: PLC0415

    eq = spec.equations[eq_idx]
    total: NDArray[np.float64] | None = None
    labels: list[str] = []

    for term_idx, term in enumerate(eq.rhs_terms):
        if (
            term.operator != "identity"
            or term.field != eq.field_name
            or term.coefficient_symbolic is None
            or not term.position_dependent
            or term.time_dependent
        ):
            continue
        resolved = coeff_eval.resolve(term, t=0.0, eq_idx=eq_idx, term_idx=term_idx)
        if not isinstance(resolved, np.ndarray):
            continue
        total = resolved.copy() if total is None else total + resolved
        labels.append(term.coefficient_symbolic)

    if total is None:
        return None, ""
    return total, " + ".join(labels)


def check_mass_sign(
    coeff_eval: CoefficientEvaluator,
    spec: EquationSystem,
) -> list[str]:
    """Check for sign-changing position-dependent mass terms.

    The check is on the **summed** self-mass of each equation, not on individual
    terms: several ``identity`` terms may act on the same field, and only their
    sum is the physical mass.

    Returns a list of warning strings for tachyonic diagnostics.
    """
    warnings: list[str] = []
    for eq_idx, eq in enumerate(spec.equations):
        profile, label = _summed_self_mass_profile(coeff_eval, spec, eq_idx)
        if profile is None:
            continue
        low, high = float(profile.min()), float(profile.max())
        if low * high < 0:
            warnings.append(
                f"Position-dependent mass term "
                f"'{label}' for field "
                f"'{eq.field_name}' changes sign across "
                f"the grid (min={low:.4g}, max={high:.4g}).",
            )
    return warnings


# ---------------------------------------------------------------------------
# Stability analysis
# ---------------------------------------------------------------------------


class StabilityResult:
    """Result from :func:`check_pointwise_mass_stability`.

    Separates fatal stability *errors* (negative eigenvalues) from
    informational *notes* (e.g. asymmetric matrix detected).
    """

    __slots__ = ("errors", "notes")

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.notes: list[str] = []


def check_pointwise_mass_stability(  # noqa: C901, PLR0912, PLR0914, PLR0915
    coeff_eval: CoefficientEvaluator,
    spec: EquationSystem,
    grid: GridInfo,
) -> StabilityResult:
    """Check eigenvalues of the pointwise mass/coupling matrix M[i,j](x,y).

    Builds M[i,j](x,y) from identity-operator terms in **dynamical**
    equations (``time_derivative_order > 0``), then verifies that all
    eigenvalues are positive at every grid point.  A negative eigenvalue
    indicates an exponentially growing (tachyonic) mode.

    Constraint equations (``time_derivative_order == 0``) are excluded
    because their algebraic form ``0 = +m²φ + ...`` produces mass signs
    opposite to the dynamical form ``d²t φ = -m²φ + ...``.  Including
    both in one matrix creates a block-indefinite matrix with spurious
    negative eigenvalues (false positives for vector field systems).

    This check runs once pre-simulation using the pre-computed spatial cache
    in CoefficientEvaluator — zero runtime cost during the actual simulation.

    Returns a :class:`StabilityResult` with ``errors`` (instability) and
    ``notes`` (informational diagnostics like asymmetry detection).

    Each row is divided by its equation's ``kinetic_coefficient_symbolic``
    before the eigenvalue check, so the eigenvalues are those of ``K⁻¹M``
    rather than of ``M`` alone (GH #237).  Without that scaling every
    non-unit-kinetic theory is judged on the wrong matrix — all torsion models
    carry ``kin = -xi`` and the gravitons ``-kappa^(-2)`` — which is why the
    propagating-torsion model was reported unstable when it is not.

    .. warning::

       This establishes positive-semidefiniteness of ``K⁻¹M``; it is **not** a
       complete stability analysis.  In particular a *negative* kinetic
       coefficient signals a ghost, which inverts how the eigenvalue sign
       should be read, and that case is not diagnosed here.  Full spectral
       health — the generalized eigenproblem and ghost identification — is
       GH #360 (PSALTer integration).

    Returns a :class:`StabilityResult` with ``errors`` (instability) and
    ``notes`` (informational diagnostics like asymmetry detection).

    Notes
    -----
    The stability condition is that the potential matrix M (defined as the
    *negative* of the identity-operator coefficient matrix) is
    positive-semidefinite at every grid point.  For coupled scalars with
    Gaussian coupling ``G(x,y) = g0*exp(-r^2/2R^2)``, this reduces to
    ``mPhi2 * mChi2 > G(x,y)^2`` everywhere, i.e. ``mPhi2 * mChi2 > g0^2``
    at the coupling peak.
    """
    import numpy as np  # noqa: PLC0415

    from tidal.solver._kinetic import build_inverse_kinetic_diag  # noqa: PLC0415

    result = StabilityResult()
    grid_shape = grid.shape

    # Only include dynamical equations (time_derivative_order > 0).
    # Constraint equations (time_derivative_order == 0) have algebraically
    # inverted mass signs that produce false-positive instability signals.
    dyn_eqs = [
        (eq_idx, eq)
        for eq_idx, eq in enumerate(spec.equations)
        if eq.time_derivative_order > 0
    ]
    dyn_names = [eq.field_name for _, eq in dyn_eqs]
    if not dyn_names:
        return result  # No dynamical fields → nothing to check

    n = len(dyn_names)

    # Build pot[i,j](x,y) as ndarray of shape (n, n, *grid_shape).
    # Convention: pot[i,j] = -(coefficient of identity(field_j) in equation_i)
    pot = np.zeros((n, n, *grid_shape))

    for eq_idx, eq in dyn_eqs:
        i = dyn_names.index(eq.field_name)
        for term_idx, term in enumerate(eq.rhs_terms):
            if term.operator != "identity":
                continue
            if term.field not in dyn_names:
                continue  # Skip references to constraint fields
            j = dyn_names.index(term.field)
            coeff = coeff_eval.resolve(term, t=0.0, eq_idx=eq_idx, term_idx=term_idx)
            if isinstance(coeff, np.ndarray):
                pot[i, j] -= coeff  # position-dependent: subtract broadcast array
            else:
                pot[i, j] -= float(coeff)  # constant: subtract scalar

    # Divide each row by its equation's kinetic coefficient, so the spectrum
    # is that of K^-1 M rather than of M alone (GH #237). Rows whose kinetic
    # coefficient is absent or unity are left untouched.
    inv_kinetic = build_inverse_kinetic_diag(
        spec,
        dict(coeff_eval.parameters),
        grid,
    )
    if inv_kinetic:
        ghosts: list[str] = []
        for i, name in enumerate(dyn_names):
            scale = inv_kinetic.get(name)
            if scale is None:
                continue
            pot[i] *= scale
            if float(np.min(scale)) < 0:
                ghosts.append(name)
        if ghosts:
            # A negative kinetic coefficient means the field enters with the
            # wrong-sign kinetic term. The eigenvalue test below still runs,
            # but its sign no longer reads as "stable"; say so rather than
            # let a caller over-trust the verdict. Full analysis is GH #360.
            result.notes.append(
                f"Fields {sorted(ghosts)} have a negative kinetic coefficient "
                "(wrong-sign kinetic term). Eigenvalue signs below are those "
                "of K^-1 M and do not by themselves establish stability for "
                "these fields; see GH #360 for ghost analysis.",
            )

    # Vectorized eigenvalue check: reshape to (n_grid, n, n) batch.
    pot_flat = pot.reshape(n, n, -1).transpose(2, 0, 1)  # (n_grid, n, n)

    # Check symmetry with *relative* tolerance: scale by matrix norm so that
    # large-amplitude systems (O(1e6)) don't trigger false asymmetry warnings
    # from floating-point roundoff.
    sym_diff = float(np.abs(pot_flat - pot_flat.transpose(0, 2, 1)).max())
    mat_scale = max(float(np.abs(pot_flat).max()), 1.0)
    if n > 1 and sym_diff > 1e-12 * mat_scale:
        result.notes.append(
            f"Mass/coupling matrix is asymmetric (max |M-M^T| = {sym_diff:.2e}). "
            f"Using general eigenvalues; stability check may be less precise.",
        )
        eigenvalues = np.linalg.eigvals(pot_flat).real  # general case
    else:
        eigenvalues = np.linalg.eigvalsh(pot_flat)  # faster, guaranteed real

    min_per_point = eigenvalues.min(axis=1)  # (n_grid,) minimum eigenvalue per point
    max_per_point = eigenvalues.max(axis=1)
    global_min = float(min_per_point.min())
    global_max = float(max_per_point.max())

    tolerance = 1e-10

    # Metric-convention-aware stability check.
    # The sign of the mass matrix that indicates stability depends on the
    # metric signature: for (-,+,+,+) the EOM is d²φ/dt² = -M·φ, so M
    # must be positive-definite.  For (+,-,-,-), M must be negative-definite.
    g_tt = spec.signature[0] if spec.signature else -1  # default (-,+,+,+)

    if g_tt < 0:
        # (-,+,+,+): stability requires M positive-definite (min eigenvalue ≥ 0)
        if global_min >= -tolerance:
            return result
    # (+,-,-,-): stability requires M negative-definite (max eigenvalue ≤ 0)
    elif global_max <= tolerance:
        return result

    # Find the worst grid point for a diagnostic message
    if g_tt < 0:
        worst_flat = int(min_per_point.argmin())
        worst_val = global_min
        requirement = "positive-definite"
    else:
        worst_flat = int(max_per_point.argmax())
        worst_val = global_max
        requirement = "negative-definite"

    worst_idx = np.unravel_index(worst_flat, grid_shape)
    # spatial coordinates = effective_coordinates[1:] (skip the time coordinate)
    spatial_coords = spec.effective_coordinates[1:]
    coord_strs = [
        f"{spatial_coords[d]}={grid.axes_coords(d)[worst_idx[d]]:.4g}"
        for d in range(len(grid_shape))
    ]
    sig_str = (
        ",".join(f"{s:+d}" for s in spec.signature) if spec.signature else "-,+,..."
    )
    result.errors.append(
        f"Coupled mass matrix has {'minimum' if g_tt < 0 else 'maximum'} "
        f"eigenvalue {worst_val:.4g} at ({', '.join(coord_strs)}). "
        f"The system has exponentially growing modes -- it will be unstable. "
        f"With signature ({sig_str}), mass matrix must be {requirement} "
        f"at all grid points.",
    )
    return result


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
                    f"Increase grid resolution or decrease gamma.",
                )
    return warnings
