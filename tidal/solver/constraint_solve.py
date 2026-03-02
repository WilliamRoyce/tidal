"""Pre-solve algebraic constraints to produce consistent initial conditions.

Before IDA starts, constraint equations (time_order=0) must be satisfied.
For example, Gauss's law ``laplacian(A_0) = -div(pi)`` requires A_0 to be
non-trivially solved when the source (pi fields) is nonzero.

Three-tier solver hierarchy (automatically selected):

    Tier 1 — **FFT** (O(N log N)): All axes periodic, all self-coefficients
    constant, all operators have known Fourier multipliers.  Uses modified
    wavenumbers ``k_mod = (2/dx)*sin(k*dx/2)`` for exact FD-consistency.

    Tier 2 — **Operator probing → sparse matrix** (O(N²) build, O(N) solve):
    Universal fallback.  Applies ``apply_operator()`` to unit vectors e_j to
    build the operator matrix column by column.  Automatically handles
    position-dependent coefficients, arbitrary BCs, and unknown/future operators.

    Tier 3 — **Iterative** (future): CG/GMRES for very large grids where
    direct factorization is impractical.

Design choices:

- Tier 2 probing reuses the **exact same** ``apply_operator()`` code path as
  the simulation, guaranteeing mathematical consistency.
- Position-dependent coefficients (NDArrays from ``CoefficientEvaluator``)
  multiply element-wise in probing, yielding correct spatially-varying matrix
  entries.  This is critical for curved-spacetime or background-field
  constraints where translational symmetry is broken.
- FFT eligibility is checked conservatively: any position-dependent self-coeff,
  non-periodic axis, or unknown multiplier → fall back to Tier 2.
- Zero-mode handling for singular operators (pure Poisson/periodic): enforce
  zero mean on the solution by setting ``u_hat[0,...,0] = 0``, with a
  compatibility check that the source has zero mean.

References
----------
- Modified wavenumbers for finite differences: Lele, J. Comp. Phys., 1992.
- Spectral constraint projection: Dedalus (Burns et al., PRR 2020).
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from tidal.solver._defaults import SECOND_ORDER
from tidal.solver._scipy_types import SparseMatrix, lil_matrix, sparse_solve
from tidal.solver.operators import BCSpec, apply_operator, is_periodic_bc

# Numerical tolerance thresholds
_SINGULAR_TOL = 1e-14  # Below this, a Fourier multiplier is treated as singular
_COMPAT_TOL = 1e-10  # Source projection tolerance for compatibility check

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.fields import FieldSet
    from tidal.solver.grid import GridInfo
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import (
        ConstraintSolverConfig,
        EquationSystem,
        OperatorTerm,
    )


# ---------------------------------------------------------------------------
# Term classification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _ConstraintTerms:
    """Classified terms for one constraint equation.

    ``self_terms`` act on the constraint field itself (form the linear
    operator L).  ``source_terms`` act on other fields (form the RHS b).

    The equation is: L[u] + b = 0  →  L[u] = -b  →  solve for u.
    """

    field_name: str
    self_terms: list[tuple[float | NDArray[np.float64], str]]
    source_terms: list[tuple[float | NDArray[np.float64], str, str]]
    eq_idx: int
    has_position_dependent_self: bool
    config: ConstraintSolverConfig


def _classify_terms(  # noqa: PLR0913, PLR0917
    eq_idx: int,
    rhs_terms: tuple[OperatorTerm, ...],
    constraint_field: str,
    coeff_eval: CoefficientEvaluator,
    t: float,
    config: ConstraintSolverConfig,
) -> _ConstraintTerms:
    """Separate RHS terms into self-operator and source terms."""
    self_terms: list[tuple[float | NDArray[np.float64], str]] = []
    source_terms: list[tuple[float | NDArray[np.float64], str, str]] = []
    has_pos_dep = False

    for term_idx, term in enumerate(rhs_terms):
        coeff = coeff_eval.resolve(term, t, eq_idx=eq_idx, term_idx=term_idx)
        # first_derivative_t(X) = dX/dt = velocity of X.
        # Resolve to identity(v_X) since this is not a spatial operator.
        if term.operator == "first_derivative_t":
            source_terms.append((coeff, "identity", f"v_{term.field}"))
        elif term.field == constraint_field:
            self_terms.append((coeff, term.operator))
            if isinstance(coeff, np.ndarray):
                has_pos_dep = True
        else:
            source_terms.append((coeff, term.operator, term.field))

    return _ConstraintTerms(
        field_name=constraint_field,
        self_terms=self_terms,
        source_terms=source_terms,
        eq_idx=eq_idx,
        has_position_dependent_self=has_pos_dep,
        config=config,
    )


# ---------------------------------------------------------------------------
# Method selection
# ---------------------------------------------------------------------------


def _lap_axis(k: NDArray[np.float64], h: float) -> NDArray[np.float64]:
    """Modified-wavenumber Laplacian along one axis: -(2/h²)(1 - cos(k*h))."""
    return -(2.0 / (h * h)) * (1.0 - np.cos(k * h))


def _grad_axis(k: NDArray[np.float64], h: float) -> NDArray[np.complex128]:
    """Modified-wavenumber gradient along one axis: i*sin(k*h)/h."""
    return 1j * np.sin(k * h) / h


# Type for Fourier multiplier functions: (kvecs, dx) -> NDArray
_MultiplierFn = Callable[
    [list[np.ndarray], tuple[float, ...]],
    np.ndarray,
]


def _biharmonic_mult(
    kvecs: list[np.ndarray],
    dx: tuple[float, ...],
) -> np.ndarray:
    """Fourier multiplier for biharmonic = (laplacian)^2."""
    lap = sum(_lap_axis(kvecs[i], dx[i]) for i in range(len(kvecs)))
    return lap * lap  # type: ignore[return-value]


def _cross_deriv(
    kvecs: list[np.ndarray],
    dx: tuple[float, ...],
    a: int,
    b: int,
) -> np.ndarray:
    """Fourier multiplier for cross_derivative: d^2/(dx_a dx_b).

    FD stencil uses central differences along each axis, giving
    multiplier = (i*sin(k_a*h_a)/h_a) * (i*sin(k_b*h_b)/h_b)
               = -sin(k_a*h_a)*sin(k_b*h_b)/(h_a*h_b).
    """
    return _grad_axis(kvecs[a], dx[a]) * _grad_axis(kvecs[b], dx[b])


_MULTIPLIERS: dict[str, _MultiplierFn] = {
    "identity": lambda kv, _dx: np.ones_like(kv[0]),
    "laplacian": lambda kv, dx: sum(  # type: ignore[return-value,arg-type]
        _lap_axis(kv[i], dx[i]) for i in range(len(kv))
    ),
    "laplacian_x": lambda kv, dx: _lap_axis(kv[0], dx[0]),
    "laplacian_y": lambda kv, dx: _lap_axis(kv[1], dx[1]),
    "laplacian_z": lambda kv, dx: _lap_axis(kv[2], dx[2]),
    "gradient_x": lambda kv, dx: _grad_axis(kv[0], dx[0]),
    "gradient_y": lambda kv, dx: _grad_axis(kv[1], dx[1]),
    "gradient_z": lambda kv, dx: _grad_axis(kv[2], dx[2]),
    "cross_derivative_xy": lambda kv, dx: _cross_deriv(kv, dx, 0, 1),
    "cross_derivative_xz": lambda kv, dx: _cross_deriv(kv, dx, 0, 2),
    "cross_derivative_yz": lambda kv, dx: _cross_deriv(kv, dx, 1, 2),
    "biharmonic": _biharmonic_mult,
}


def _select_method(
    terms: _ConstraintTerms,
    grid: GridInfo,
    bc: BCSpec | None,
) -> str:
    """Determine solver method: 'fft' or 'matrix'.

    User override via config.method takes precedence.  Otherwise, checks
    FFT eligibility (all-periodic, constant coefficients, known multipliers)
    and falls back to matrix (operator probing).
    """
    config = terms.config
    if config.method not in {"auto", "poisson"}:
        return config.method

    # Check all-periodic
    all_periodic = all(grid.periodic)
    if bc is not None:
        bcs = (bc,) * grid.ndim if isinstance(bc, str) else tuple(bc)
        all_periodic = all(is_periodic_bc(b) for b in bcs)

    if not all_periodic:
        return "matrix"

    if terms.has_position_dependent_self:
        return "matrix"

    for _, op_name in terms.self_terms:
        if op_name not in _MULTIPLIERS:
            return "matrix"

    return "fft"


# ---------------------------------------------------------------------------
# Field name resolution
# ---------------------------------------------------------------------------


def _build_name_map(spec: EquationSystem) -> dict[str, str]:
    """Build a map from JSON field references to FieldSet slot names.

    The canonical naming convention is ``v_{field_name}`` (e.g.
    ``v_A_1``), matching ``StateLayout`` slot names.
    """
    name_map: dict[str, str] = {}
    for eq in spec.equations:
        # Field names map to themselves
        name_map[eq.field_name] = eq.field_name
        # Velocity references: canonical v_field_name
        if eq.time_derivative_order >= SECOND_ORDER:
            vel_slot = f"v_{eq.field_name}"
            name_map[vel_slot] = vel_slot
    return name_map


# ---------------------------------------------------------------------------
# Source evaluation
# ---------------------------------------------------------------------------


def _evaluate_source(
    source_terms: list[tuple[float | NDArray[np.float64], str, str]],
    fields: FieldSet,
    grid: GridInfo,
    bc: BCSpec | None,
    name_map: dict[str, str] | None = None,
) -> NDArray[np.float64]:
    """Compute source RHS: b = Σ coeff_i * operator_i(field_i)."""
    result: NDArray[np.float64] = np.zeros(grid.shape)
    for coeff, op_name, field_name in source_terms:
        resolved = field_name
        if name_map and field_name in name_map:
            resolved = name_map[field_name]
        data = fields[resolved] if resolved in fields else np.zeros(grid.shape)
        operated = apply_operator(op_name, data, grid, bc)
        result += coeff * operated
    return result


# ---------------------------------------------------------------------------
# Tier 1: FFT solver
# ---------------------------------------------------------------------------


def _build_wavenumber_grids(
    grid: GridInfo,
) -> tuple[list[Any], tuple[float, ...]]:
    """Build wavenumber grids for FFT-based solving."""
    k_1d = [
        2.0 * np.pi * np.fft.fftfreq(grid.shape[ax], d=grid.dx[ax])
        for ax in range(grid.ndim)
    ]
    kvecs = list(np.meshgrid(*k_1d, indexing="ij")) if grid.ndim > 1 else [k_1d[0]]
    return kvecs, grid.dx


def _fft_solve_single(
    terms: _ConstraintTerms,
    grid: GridInfo,
    source_rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Solve a single constraint via FFT.

    Computes: u = F^{-1}[ -F[source] / multiplier ]

    Raises
    ------
    ValueError
        If the operator is singular and the source has nonzero projection
        onto the corresponding mode.
    """
    kvecs, dx = _build_wavenumber_grids(grid)

    # Build combined multiplier from self-terms
    multiplier: NDArray[np.complex128] = np.zeros(grid.shape, dtype=np.complex128)
    for coeff, op_name in terms.self_terms:
        fn = _MULTIPLIERS[op_name]
        multiplier += coeff * fn(kvecs, dx)

    # Transform source
    source_hat = np.fft.fftn(source_rhs)

    # Detect singular modes
    is_singular = np.abs(multiplier) < _SINGULAR_TOL

    if np.any(is_singular):
        source_at_singular = np.abs(source_hat[is_singular])
        max_source = float(np.max(np.abs(source_hat))) + 1e-30
        max_incompatible = (
            float(np.max(source_at_singular)) if source_at_singular.size > 0 else 0.0
        )
        if max_incompatible > _COMPAT_TOL * max_source:
            msg = (
                f"Constraint for '{terms.field_name}' is incompatible: "
                f"source has nonzero projection (max={max_incompatible:.4g}) "
                f"onto the null space of the self-operator. "
                f"Check that source terms have zero mean for "
                f"Poisson-type constraints."
            )
            raise ValueError(msg)
        n_singular = int(np.sum(is_singular))
        warnings.warn(
            f"FFT constraint pre-solve: field '{terms.field_name}' has "
            f"{n_singular} singular mode(s) in Fourier space (null space "
            f"of operator). Setting u_hat = 0 at these modes (zero-mean "
            f"gauge). Solution is unique up to these modes.",
            UserWarning,
            stacklevel=3,
        )
        safe_mult = np.where(is_singular, 1.0, multiplier)
        u_hat = np.where(is_singular, 0.0, -source_hat / safe_mult)
    else:
        u_hat = -source_hat / multiplier

    return np.fft.ifftn(u_hat).real.astype(np.float64)


def _fft_solve_coupled(  # noqa: PLR0914
    groups: list[_ConstraintTerms],
    grid: GridInfo,
    fields: FieldSet,
    bc: BCSpec | None,
    name_map: dict[str, str] | None = None,
) -> dict[str, NDArray[np.float64]]:
    """Solve coupled constraints via block FFT.

    At each wavenumber k, assembles and solves the n-by-n system.

    Raises
    ------
    ValueError
        If the coupled system is singular and the source has nonzero mean.
    """
    n_c = len(groups)
    kvecs, dx = _build_wavenumber_grids(grid)

    constraint_names = [g.field_name for g in groups]
    name_to_idx = {name: i for i, name in enumerate(constraint_names)}

    # Evaluate sources from non-constraint fields
    sources: list[NDArray[np.float64]] = []
    for group in groups:
        non_constraint = [
            (c, op, f) for c, op, f in group.source_terms if f not in name_to_idx
        ]
        sources.append(_evaluate_source(non_constraint, fields, grid, bc, name_map))

    # Build multiplier matrix and RHS in Fourier space
    m_hat = np.zeros((*grid.shape, n_c, n_c), dtype=np.complex128)
    rhs_hat = np.zeros((*grid.shape, n_c), dtype=np.complex128)

    for i, group in enumerate(groups):
        # Self-terms → diagonal
        for coeff, op_name in group.self_terms:
            fn = _MULTIPLIERS[op_name]
            m_hat[..., i, i] += coeff * fn(kvecs, dx)

        # Cross-constraint terms → off-diagonal
        for coeff, op_name, field_name in group.source_terms:
            if field_name in name_to_idx:
                j = name_to_idx[field_name]
                fn = _MULTIPLIERS[op_name]
                m_hat[..., i, j] += coeff * fn(kvecs, dx)

        rhs_hat[..., i] = -np.fft.fftn(sources[i])

    # Handle singular zero-mode
    zero_idx = tuple(0 for _ in range(grid.ndim))
    m_zero = m_hat[zero_idx]
    if abs(np.linalg.det(m_zero)) < _SINGULAR_TOL:
        rhs_zero = rhs_hat[zero_idx]
        if float(np.max(np.abs(rhs_zero))) > _COMPAT_TOL:
            msg = (
                "Coupled constraint system is singular at zero wavenumber "
                "and source has nonzero mean. Check compatibility."
            )
            raise ValueError(msg)
        coupled_names = ", ".join(constraint_names)
        warnings.warn(
            f"FFT coupled constraint pre-solve: system "
            f"[{coupled_names}] is singular at zero wavenumber (null "
            f"space of operator). Setting zero-mode to identity/zero "
            f"(zero-mean gauge). Solution is unique up to constants.",
            UserWarning,
            stacklevel=3,
        )
        m_hat[zero_idx] = np.eye(n_c)
        rhs_hat[zero_idx] = 0.0

    # np.linalg.solve dispatches to the matrix path (m,m),(m,n)->(m,n)
    # when ndim >= 2, but rhs_hat is (..., n_c) — a vector per wavenumber.
    # Add a trailing dim so numpy sees (..., n_c, 1), then squeeze it out.
    u_hat = np.linalg.solve(m_hat, rhs_hat[..., np.newaxis])[..., 0]

    results: dict[str, NDArray[np.float64]] = {}
    for i, name in enumerate(constraint_names):
        results[name] = np.fft.ifftn(u_hat[..., i]).real.astype(np.float64)

    return results


# ---------------------------------------------------------------------------
# Tier 2: Operator probing → sparse matrix
# ---------------------------------------------------------------------------


def _probe_operator_matrix(
    self_terms: list[tuple[float | NDArray[np.float64], str]],
    grid: GridInfo,
    bc: BCSpec | None,
) -> SparseMatrix:
    """Build sparse matrix by probing apply_operator() with unit vectors.

    For single constant-coefficient terms, delegates to the shared
    ``build_operator_matrix()`` utility.  For multiple or position-dependent
    terms, builds the matrix by probing the combined operator.
    """
    # Fast path: single term with constant coefficient
    if len(self_terms) == 1:
        coeff, op_name = self_terms[0]
        if np.ndim(coeff) == 0:
            from tidal.solver.analytical_jacobian import (  # noqa: PLC0415
                build_operator_matrix,
            )

            mat = build_operator_matrix(op_name, grid, bc)
            return float(coeff) * mat if float(coeff) != 1.0 else mat

    # General path: multiple terms and/or position-dependent coefficients
    n = grid.num_points
    mat = lil_matrix((n, n))

    for j in range(n):
        e_j: NDArray[np.float64] = np.zeros(grid.shape)
        e_j.flat[j] = 1.0

        col: NDArray[np.float64] = np.zeros(grid.shape)
        for coeff, op_name in self_terms:
            col += coeff * apply_operator(op_name, e_j, grid, bc)

        col_flat = col.ravel()
        nz = np.nonzero(col_flat)[0]
        for row in nz:
            mat[row, j] = col_flat[row]

    return mat.tocsc()


def _matrix_solve(
    op_matrix: SparseMatrix,
    source_rhs: NDArray[np.float64],
    grid_shape: tuple[int, ...],
) -> NDArray[np.float64]:
    """Solve op_matrix @ u = -source_rhs via sparse direct factorization."""
    rhs = -source_rhs.ravel()
    return sparse_solve(op_matrix, rhs).reshape(grid_shape)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def pre_solve_constraints(  # noqa: PLR0913
    spec: EquationSystem,
    grid: GridInfo,
    y0: NDArray[np.float64],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    t: float = 0.0,
) -> NDArray[np.float64]:
    """Solve algebraic constraints to produce consistent initial conditions.

    Only processes equations where ``time_derivative_order == 0`` AND
    ``constraint_solver.enabled == True``.  Returns a copy of y0 with
    constraint field slots overwritten by the solution.

    Parameters
    ----------
    spec : EquationSystem
        Parsed JSON equation specification.
    grid : GridInfo
        Spatial grid.
    y0 : NDArray[np.float64]
        Initial state vector (flat, from StateLayout).
    bc : str or tuple of str, optional
        Boundary conditions for spatial operators.
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients.
    t : float
        Time at which to evaluate coefficients (usually t_span[0]).

    Returns
    -------
    NDArray[np.float64]
        Updated y0 with constraint fields solved.

    Raises
    ------
    ValueError
        If a constraint is incompatible (singular operator with nonzero
        projection in null space), or has no self-terms.

    Warns
    -----
    UserWarning
        When the FFT solver encounters singular modes (null space of the
        operator) and regularizes by setting ``u_hat = 0`` at those modes.
        This is a numerical gauge choice (zero-mean).  To disable
        automatic constraint solving for a field, set
        ``constraint_solver.enabled = false`` in the JSON spec.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
    from tidal.solver.fields import FieldSet  # noqa: PLC0415
    from tidal.solver.state import StateLayout  # noqa: PLC0415

    # Find enabled constraint equations
    constraint_eqs = [
        (i, eq)
        for i, eq in enumerate(spec.equations)
        if eq.time_derivative_order == 0 and eq.constraint_solver.enabled
    ]

    if not constraint_eqs:
        return y0

    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, parameters)
    coeff_eval.begin_timestep(t)

    # Copy y0 BEFORE creating FieldSet — FieldSet wraps its backing array
    # with zero-copy views, so mutations via fields[name]=... propagate back.
    y0_out = y0.copy()
    fields = FieldSet.from_flat(layout, grid.shape, y0_out)

    # Build name map: JSON field references (e.g. "v_A_1") → FieldSet slot names
    name_map = _build_name_map(spec)

    # Classify terms for each constraint
    groups: list[_ConstraintTerms] = []
    for eq_idx, eq in constraint_eqs:
        terms = _classify_terms(
            eq_idx,
            eq.rhs_terms,
            eq.field_name,
            coeff_eval,
            t,
            eq.constraint_solver,
        )
        if not terms.self_terms:
            msg = (
                f"Constraint equation for '{eq.field_name}' has no "
                f"self-referencing terms — cannot solve for the field."
            )
            raise ValueError(msg)
        groups.append(terms)

    # Detect coupled constraints (check both original and resolved names)
    constraint_names = {g.field_name for g in groups}
    is_coupled = any(
        any(
            f in constraint_names or name_map.get(f, f) in constraint_names
            for _, _, f in g.source_terms
        )
        for g in groups
    )

    if is_coupled:
        _solve_coupled(groups, grid, bc, fields, layout, y0_out, name_map)
    else:
        _solve_independent(groups, grid, bc, fields, layout, y0_out, name_map)

    return y0_out


def _solve_independent(  # noqa: PLR0913, PLR0917
    groups: list[_ConstraintTerms],
    grid: GridInfo,
    bc: BCSpec | None,
    fields: FieldSet,
    layout: StateLayout,
    y0: NDArray[np.float64],
    name_map: dict[str, str] | None = None,
) -> None:
    """Solve independent (non-coupled) constraints one at a time."""
    n = grid.num_points

    for terms in groups:
        source = _evaluate_source(terms.source_terms, fields, grid, bc, name_map)
        method = _select_method(terms, grid, bc)

        if method == "fft":
            solution = _fft_solve_single(terms, grid, source)
        else:
            op_mat = _probe_operator_matrix(terms.self_terms, grid, bc)
            solution = _matrix_solve(op_mat, source, grid.shape)

        slot_idx = layout.field_slot_map[terms.field_name]
        y0[slot_idx * n : (slot_idx + 1) * n] = solution.ravel()
        fields[terms.field_name] = solution


def _solve_coupled(  # noqa: PLR0912, PLR0913, PLR0914, PLR0917, C901
    groups: list[_ConstraintTerms],
    grid: GridInfo,
    bc: BCSpec | None,
    fields: FieldSet,
    layout: StateLayout,
    y0: NDArray[np.float64],
    name_map: dict[str, str] | None = None,
) -> None:
    """Solve coupled constraints (e.g., coupled Proca A_0 ↔ B_0)."""
    n = grid.num_points

    # Check if FFT path is available for ALL coupled constraints.
    # Must also verify that cross-constraint source operators have FFT
    # multipliers, since _select_method only checks self-terms.
    constraint_names = {g.field_name for g in groups}
    all_fft = all(_select_method(g, grid, bc) == "fft" for g in groups)
    if all_fft:
        for g in groups:
            for _, op_name, field_name in g.source_terms:
                if field_name in constraint_names and op_name not in _MULTIPLIERS:
                    all_fft = False
                    break

    if all_fft:
        solutions = _fft_solve_coupled(groups, grid, fields, bc, name_map)
        for name, sol in solutions.items():
            slot_idx = layout.field_slot_map[name]
            y0[slot_idx * n : (slot_idx + 1) * n] = sol.ravel()
            fields[name] = sol
    else:
        # Gauss-Seidel iteration
        max_iter = max(g.config.max_iterations for g in groups)
        tol = min(g.config.tolerance for g in groups)

        # Pre-compute operator matrices and methods (loop-invariant)
        methods: dict[str, str] = {}
        op_matrices: dict[str, np.ndarray] = {}
        for terms in groups:
            method = _select_method(terms, grid, bc)
            methods[terms.field_name] = method
            if method != "fft":
                op_matrices[terms.field_name] = _probe_operator_matrix(
                    terms.self_terms, grid, bc
                )

        for _iteration in range(max_iter):
            max_change = 0.0
            for terms in groups:
                source = _evaluate_source(
                    terms.source_terms, fields, grid, bc, name_map
                )
                method = methods[terms.field_name]

                if method == "fft":
                    solution = _fft_solve_single(terms, grid, source)
                else:
                    op_mat = op_matrices[terms.field_name]
                    solution = _matrix_solve(op_mat, source, grid.shape)

                old = (
                    fields[terms.field_name]
                    if terms.field_name in fields
                    else np.zeros(grid.shape)
                )
                change = float(np.max(np.abs(solution - old)))
                max_change = max(max_change, change)

                slot_idx = layout.field_slot_map[terms.field_name]
                y0[slot_idx * n : (slot_idx + 1) * n] = solution.ravel()
                fields[terms.field_name] = solution

            if max_change < tol:
                break


# ---------------------------------------------------------------------------
# Unified constraint IC solver
# ---------------------------------------------------------------------------

# Tolerance for considering a field "zero" (uninitialized)
_FIELD_ZERO_TOL = 1e-30

# Tolerance for constraint verification residual
_IC_RESIDUAL_TOL = 1e-10

# Maximum number of propagation iterations
_MAX_PROPAGATION_ITERATIONS = 20


def _is_field_zero(data: np.ndarray) -> bool:
    """Check if a field array is effectively zero (uninitialized)."""
    return float(np.max(np.abs(data))) < _FIELD_ZERO_TOL


def _solve_for_target(
    terms: _ConstraintTerms,
    grid: GridInfo,
    bc: BCSpec | None,
    source_rhs: np.ndarray,
) -> np.ndarray:
    """Solve L[target] = -source for a single target field.

    Reuses the existing FFT/matrix solver infrastructure.
    """
    method = _select_method(terms, grid, bc)
    if method == "fft":
        return _fft_solve_single(terms, grid, source_rhs)
    op_mat = _probe_operator_matrix(terms.self_terms, grid, bc)
    return _matrix_solve(op_mat, source_rhs, grid.shape)


def _find_target_field(  # noqa: PLR0913, PLR0917
    eq_idx: int,
    rhs_terms: tuple[OperatorTerm, ...],
    eq_field_name: str,  # noqa: ARG001 — kept for caller consistency
    determined: set[str],
    coeff_eval: CoefficientEvaluator,
    t: float,
    config: ConstraintSolverConfig,
) -> tuple[str | None, _ConstraintTerms | None]:
    """Identify which field to solve for in a constraint equation.

    Returns ``(target_field_name, classified_terms)`` or ``(None, None)``
    if the constraint cannot be solved in this iteration.

    Strategy: find exactly ONE free (undetermined) field among all
    referenced fields.  That field becomes the target.  If the target
    has self-referencing operator terms, it can be solved.

    - 0 free fields → all determined, verification only (return None).
    - 1 free field → solvable (the free field is the target).
    - 2+ free fields → underdetermined, skip this iteration.

    The single-free-field check naturally handles both:

    - **Standard constraints** (constraint field is the free one):
      e.g., ``laplacian(A_0) + source = 0`` when A_0 is free.
    - **Subsidiary constraints** (dynamical field is the free one):
      e.g., ``gradient(h_4) + gradient(h_7) = 0`` when only h_7 is free.
    - **Cascade ordering**: ``identity(h_4) + identity(h_7) + identity(h_9) = 0``
      waits until h_7 is determined before solving for h_9.
    """
    # Collect all field names referenced by this equation
    all_referenced: set[str] = {term.field for term in rhs_terms}

    # Find free (undetermined) fields
    free_fields: list[str] = [
        f for f in all_referenced if f not in determined
    ]

    if len(free_fields) != 1:
        return None, None

    target = free_fields[0]
    terms = _classify_terms(
        eq_idx, rhs_terms, target, coeff_eval, t, config,
    )
    if terms.self_terms:
        return target, terms

    return None, None


def ensure_consistent_ic(  # noqa: PLR0912, PLR0913, PLR0914, PLR0915, C901
    spec: EquationSystem,
    grid: GridInfo,
    y0: NDArray[np.float64],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    t: float = 0.0,
    strict: bool = True,
) -> NDArray[np.float64]:
    """Unified constraint IC solver.

    Given user-supplied initial conditions, iteratively solves ALL
    constraint equations (``time_derivative_order == 0``) to produce
    consistent initial data.  Handles three cases uniformly:

    1. **Standard constraints** (field has self-terms): solve for the
       constraint field.  Applies to all constraint equations, not just
       ``constraint_solver.enabled=True``.

    2. **Subsidiary constraints** (no self-terms, one free source field):
       solve for the free dynamical field.  E.g., transverse gauge
       constraint ``gradient_x(h_4) + gradient_x(h_7) = 0`` determines
       ``h_7 = -h_4`` when ``h_4`` is user-initialized.

    3. **Verification** (all referenced fields determined): check that
       the constraint is satisfied; error or warn if not.

    The algorithm iterates until no more fields can be determined.
    Each iteration may unlock new solvable constraints (cascade).

    Parameters
    ----------
    spec : EquationSystem
        Parsed JSON equation specification.
    grid : GridInfo
        Spatial grid.
    y0 : NDArray[np.float64]
        Initial state vector (flat, from StateLayout).
    bc : str or tuple of str, optional
        Boundary conditions for spatial operators.
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients.
    t : float
        Time at which to evaluate coefficients.
    strict : bool
        If True (default), raise ``ValueError`` when a constraint is
        violated and cannot be solved.  If False, issue a warning.

    Returns
    -------
    NDArray[np.float64]
        Updated y0 with all solvable constraint fields determined.

    Raises
    ------
    ValueError
        If ``strict=True`` and any constraint is violated after solving.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415
    from tidal.solver.fields import FieldSet  # noqa: PLC0415
    from tidal.solver.state import StateLayout  # noqa: PLC0415

    # Collect ALL constraint equations
    constraint_eqs = [
        (i, eq)
        for i, eq in enumerate(spec.equations)
        if eq.time_derivative_order == 0
    ]

    if not constraint_eqs:
        return y0

    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, parameters)
    coeff_eval.begin_timestep(t)
    n = grid.num_points

    y0_out = y0.copy()
    fields = FieldSet.from_flat(layout, grid.shape, y0_out)
    name_map = _build_name_map(spec)

    # Track which fields are "determined" (non-zero in y0).
    determined: set[str] = set()
    for name, slot_idx in layout.field_slot_map.items():
        data = y0_out[slot_idx * n : (slot_idx + 1) * n]
        if not _is_field_zero(data):
            determined.add(name)
    for name, slot_idx in layout.velocity_slot_map.items():
        data = y0_out[slot_idx * n : (slot_idx + 1) * n]
        if not _is_field_zero(data):
            determined.add(f"v_{name}")

    # Phase 1: Handle enabled constraints with coupled detection
    # (preserves existing coupled-constraint behavior for EM, Chern-Simons)
    enabled_eqs = [
        (i, eq) for i, eq in constraint_eqs if eq.constraint_solver.enabled
    ]
    if enabled_eqs:
        enabled_groups: list[_ConstraintTerms] = []
        for eq_idx, eq in enabled_eqs:
            terms = _classify_terms(
                eq_idx, eq.rhs_terms, eq.field_name,
                coeff_eval, t, eq.constraint_solver,
            )
            if not terms.self_terms:
                msg = (
                    f"Constraint equation for '{eq.field_name}' has "
                    f"constraint_solver.enabled=True but no self-referencing "
                    f"terms — cannot solve for the field."
                )
                raise ValueError(msg)
            enabled_groups.append(terms)

        # Detect coupled constraints among enabled group
        enabled_names = {g.field_name for g in enabled_groups}
        is_coupled = any(
            any(
                f in enabled_names or name_map.get(f, f) in enabled_names
                for _, _, f in g.source_terms
            )
            for g in enabled_groups
        )

        if is_coupled:
            _solve_coupled(
                enabled_groups, grid, bc, fields, layout, y0_out, name_map,
            )
        else:
            _solve_independent(
                enabled_groups, grid, bc, fields, layout, y0_out, name_map,
            )

        determined.update(g.field_name for g in enabled_groups)

    # Phase 2: Iterative propagation for remaining constraints
    remaining_eqs = [
        (i, eq) for i, eq in constraint_eqs
        if not eq.constraint_solver.enabled
    ]

    for _iteration in range(_MAX_PROPAGATION_ITERATIONS):
        progress = False

        for eq_idx, eq in remaining_eqs:
            target, terms = _find_target_field(
                eq_idx, eq.rhs_terms, eq.field_name,
                determined, coeff_eval, t, eq.constraint_solver,
            )

            if target is None or terms is None:
                continue

            source = _evaluate_source(
                terms.source_terms, fields, grid, bc, name_map,
            )

            solution = _solve_for_target(terms, grid, bc, source)

            # Write solution into y0 and FieldSet
            if target in layout.field_slot_map:
                slot_idx = layout.field_slot_map[target]
                y0_out[slot_idx * n : (slot_idx + 1) * n] = solution.ravel()
                fields[target] = solution
            elif target.startswith("v_"):
                vel_field = target.removeprefix("v_")
                if vel_field in layout.velocity_slot_map:
                    slot_idx = layout.velocity_slot_map[vel_field]
                    y0_out[slot_idx * n : (slot_idx + 1) * n] = (
                        solution.ravel()
                    )
                    fields[target] = solution

            determined.add(target)
            progress = True

            eq_desc = eq.field_name
            if target != eq.field_name:
                eq_desc = f"{eq.field_name} (subsidiary)"
            warnings.warn(
                f"Constraint IC propagation: solved '{target}' "
                f"from {eq_desc} constraint equation.",
                UserWarning,
                stacklevel=2,
            )

        if not progress:
            break

    # Phase 3: Final verification of ALL constraint equations
    violations: list[tuple[str, float, list[str]]] = []
    for eq_idx, eq in constraint_eqs:
        rhs = np.zeros(grid.shape)
        for term_idx, term in enumerate(eq.rhs_terms):
            coeff = coeff_eval.resolve(
                term, t, eq_idx=eq_idx, term_idx=term_idx,
            )
            # first_derivative_t(X) → identity(v_X)
            op_name = term.operator
            field_ref = term.field
            if op_name == "first_derivative_t":
                op_name = "identity"
                field_ref = f"v_{field_ref}"
            if name_map and field_ref in name_map:
                field_ref = name_map[field_ref]
            data = (
                fields[field_ref]
                if field_ref in fields
                else np.zeros(grid.shape)
            )
            operated = apply_operator(op_name, data, grid, bc)
            rhs += coeff * operated

        max_res = float(np.max(np.abs(rhs)))
        if max_res > _IC_RESIDUAL_TOL:
            involved = sorted({term.field for term in eq.rhs_terms})
            violations.append((eq.field_name, max_res, involved))

    if violations:
        lines = [
            "Initial data does not satisfy constraint equation(s):"
        ]
        for field_name, max_res, involved in violations:
            lines.append(
                f"  {field_name}: max|residual| = {max_res:.2e} "
                f"(involves [{', '.join(involved)}])"
            )
        lines.append(
            "For physical consistency, choose initial conditions that "
            "jointly satisfy all constraint equations, or use "
            "--allow-inconsistent-ic to proceed with a warning."
        )
        msg = "\n".join(lines)
        if strict:
            raise ValueError(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)

    return y0_out
