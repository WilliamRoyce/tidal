"""Fourier modal solver — exact spectral time evolution for linear PDEs.

Transforms the spatial grid to Fourier space, builds a per-mode evolution
matrix, and eigendecomposes to obtain the exact solution y(t) = exp(A·t)·y₀.
Eliminates spatial discretization error entirely and provides machine-precision
solutions for time-independent linear systems.

Applicable to any linear PDE system with:
- Flat (Minkowski) metric
- All-periodic boundary conditions
- Time-independent coefficients (position-dependent OK via convolution)
- Operators with known exact Fourier multipliers

Two algorithm paths are used depending on coefficient structure:
- Constant coefficients: per-mode eigendecomposition with block-aware independent
  blocks (machine-precision, ~14x faster).
- Position-dependent coefficients: Krylov matrix exponential (expm_multiply) which
  is backward-stable for non-normal convolution matrices where eigendecomposition
  gives incorrect results due to pseudospectral overflow.

References
----------
    Moler & Van Loan (2003), SIAM Review 45(1):3-49 — matrix exponential.
    Hairer, Lubich & Wanner (2006), Geometric Numerical Integration, §4.
    Burns et al. (2020), Phys. Rev. Research 2:023068 — pseudo-spectral.
    Raffelt & Stodolsky (1988), PRD 37:1237 — mixing-matrix formalism.
"""

# ruff: noqa: N803, N806 — uppercase names for matrices (A, V, T, Z) follow
#   standard linear-algebra notation.
# ruff: noqa: PLR0913, PLR0917, PLR0914, PLR0912, PLR0911, PLR0915, PLR2004
#   — numerical code inherently requires many arguments, local variables,
#   return statements, statements, and literal comparisons.
# ruff: noqa: C901 — complexity and Unicode math symbols.
# ruff: noqa: ERA001, ARG001 — commented-out code serves as documentation;
#   unused args (bc, grid) kept for interface consistency with other solvers.
# ruff: noqa: B903, PLR1702 — _OperatorDecomp uses __slots__ for memory efficiency;
#   nested block depth is inherent to multi-field modal algebra.

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._setup import warn_frozen_constraints
from tidal.solver.operators import get_wavenumbers, is_periodic_bc
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from tidal.solver._types import SolverResult
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.progress import SimulationProgress
    from tidal.symbolic.json_loader import (
        ComponentEquation,
        EquationSystem,
        OperatorTerm,
    )

# ---------------------------------------------------------------------------
# Exact Fourier multipliers (angular wavenumber convention: k = 2π·rfftfreq)
# ---------------------------------------------------------------------------
# These use the EXACT wavenumber k, consistent with operators.py spectral
# mode (gradient → ik, laplacian → -k²).  NOT the modified-wavenumber
# convention from constraint_solve.py which matches FD stencils.

_ExactMultFn = Callable[[list[NDArray[np.float64]]], NDArray[Any] | int]

# ---------------------------------------------------------------------------
# Operator decomposition: (spatial Fourier multiplier, time derivative order)
# ---------------------------------------------------------------------------
# Every operator in flat spacetime decomposes as spatial_multiplier(k) x ∂ⁿ_t.
# Derivatives commute in Minkowski: ∂²_t ∂_x = ∂_x ∂²_t.
#
# Time order classification:
#   0 = position operator   → stiffness matrix K
#   1 = velocity operator   → damping matrix D
#   2 = acceleration operator → mass matrix M (off-diagonal, implicit coupling)
#   3 = jerk operator       → eliminated via EOM substitution
#
# References:
#   Golub & Van Loan (2013), Matrix Computations, 4th ed. §7.7
#   Hairer & Lubich (2003), ZAMM 83(1) — mass matrices in structural dynamics


class _OperatorDecomp:
    """Operator decomposition into spatial multiplier and time derivative order."""

    __slots__ = ("spatial_fn", "time_order")

    def __init__(self, spatial_fn: _ExactMultFn, time_order: int) -> None:
        self.spatial_fn = spatial_fn
        self.time_order = time_order


_OPERATOR_DECOMP: dict[str, _OperatorDecomp] = {
    # --- Pure spatial operators (time_order=0) ---
    "identity": _OperatorDecomp(lambda k_axes: np.ones_like(k_axes[0]), 0),
    "laplacian": _OperatorDecomp(
        lambda k_axes: -sum(ki**2 for ki in k_axes),
        0,
    ),
    "laplacian_x": _OperatorDecomp(lambda k_axes: -(k_axes[0] ** 2), 0),
    "laplacian_y": _OperatorDecomp(lambda k_axes: -(k_axes[1] ** 2), 0),
    "laplacian_z": _OperatorDecomp(lambda k_axes: -(k_axes[2] ** 2), 0),
    "gradient_x": _OperatorDecomp(lambda k_axes: 1j * k_axes[0], 0),
    "gradient_y": _OperatorDecomp(lambda k_axes: 1j * k_axes[1], 0),
    "gradient_z": _OperatorDecomp(lambda k_axes: 1j * k_axes[2], 0),
    "cross_derivative_xy": _OperatorDecomp(lambda k_axes: -(k_axes[0] * k_axes[1]), 0),
    "cross_derivative_xz": _OperatorDecomp(lambda k_axes: -(k_axes[0] * k_axes[2]), 0),
    "cross_derivative_yz": _OperatorDecomp(lambda k_axes: -(k_axes[1] * k_axes[2]), 0),
    "biharmonic": _OperatorDecomp(
        lambda k_axes: sum(ki**2 for ki in k_axes) ** 2,
        0,
    ),
    "derivative_3_x": _OperatorDecomp(lambda k_axes: -1j * k_axes[0] ** 3, 0),
    "derivative_3_y": _OperatorDecomp(lambda k_axes: -1j * k_axes[1] ** 3, 0),
    "derivative_3_z": _OperatorDecomp(lambda k_axes: -1j * k_axes[2] ** 3, 0),
    # --- Velocity operators (time_order=1) ---
    "first_derivative_t": _OperatorDecomp(lambda k_axes: np.ones_like(k_axes[0]), 1),
    "mixed_T1_S1x": _OperatorDecomp(lambda k_axes: 1j * k_axes[0], 1),
    "mixed_T1_S1y": _OperatorDecomp(lambda k_axes: 1j * k_axes[1], 1),
    "mixed_T1_S1z": _OperatorDecomp(lambda k_axes: 1j * k_axes[2], 1),
    # --- Acceleration operators (time_order=2) ---
    "d2_t": _OperatorDecomp(lambda k_axes: np.ones_like(k_axes[0]), 2),
    "mixed_T2_S1x": _OperatorDecomp(lambda k_axes: 1j * k_axes[0], 2),
    "mixed_T2_S1y": _OperatorDecomp(lambda k_axes: 1j * k_axes[1], 2),
    "mixed_T2_S1z": _OperatorDecomp(lambda k_axes: 1j * k_axes[2], 2),
    "mixed_T2_S2x": _OperatorDecomp(lambda k_axes: -(k_axes[0] ** 2), 2),
    "mixed_T2_S2y": _OperatorDecomp(lambda k_axes: -(k_axes[1] ** 2), 2),
    "mixed_T2_S2z": _OperatorDecomp(lambda k_axes: -(k_axes[2] ** 2), 2),
    # --- Jerk operators (time_order=3, eliminated via EOM substitution) ---
    "d3_t": _OperatorDecomp(lambda k_axes: np.ones_like(k_axes[0]), 3),
    "mixed_T3_S1x": _OperatorDecomp(lambda k_axes: 1j * k_axes[0], 3),
    "mixed_T3_S1y": _OperatorDecomp(lambda k_axes: 1j * k_axes[1], 3),
    "mixed_T3_S1z": _OperatorDecomp(lambda k_axes: 1j * k_axes[2], 3),
}

# Backward-compatible mapping: operator name → spatial multiplier function.
# Used by existing code paths that only need the spatial part.
_EXACT_MULTIPLIERS: dict[str, _ExactMultFn] = {
    name: dec.spatial_fn for name, dec in _OPERATOR_DECOMP.items()
}


# ---------------------------------------------------------------------------
# Eligibility check
# ---------------------------------------------------------------------------


def can_use_modal(
    spec: EquationSystem,
    grid: GridInfo,
    bc: BCSpec | None,
) -> bool:
    """Check whether the modal solver is applicable to this system.

    Requirements (checked in order):
    1. Flat metric (volume_element is None)
    2. Constraints (time_order=0) must be Fourier-eliminable via Schur complement
    3. All boundary conditions periodic
    4. All RHS operators have exact Fourier multipliers
    5. No time-dependent coefficients
    """
    # 1. Flat metric — curved metrics have non-None volume_element
    if spec.canonical is not None and spec.canonical.volume_element is not None:
        return False
    # Also reject if canonical is None but any term is position-dependent
    # with non-Cartesian coordinate names (heuristic for curved metrics
    # without canonical structure)
    if spec.canonical is None:
        for eq in spec.equations:
            for term in eq.rhs_terms:
                if term.position_dependent:
                    coords = set(term.coordinate_dependent)
                    # Non-Cartesian coordinate references suggest curved metric
                    cartesian = {"x", "y", "z"}
                    if coords - cartesian:
                        return False

    # 2. Constraints — allow if Fourier-eliminable via Schur complement
    constraint_eqs = [eq for eq in spec.equations if eq.time_derivative_order == 0]
    if constraint_eqs and not _constraints_fourier_eliminable(spec, constraint_eqs):
        return False

    # 3. All-periodic BCs
    if not all(grid.periodic):
        return False
    if bc is not None:
        if isinstance(bc, str):
            if not is_periodic_bc(bc):
                return False
        else:
            for b in bc:
                if not is_periodic_bc(b):
                    return False

    # 4. All operators supported (spatial or time-derivative decomposable)
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.operator not in _OPERATOR_DECOMP:
                return False

    # 5. No time-dependent coefficients
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.time_dependent:
                return False

    return True


# ---------------------------------------------------------------------------
# FFT state transforms
# ---------------------------------------------------------------------------


def _fft_slots(
    y: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
) -> NDArray[np.complex128]:
    """Transform each slot from physical space to Fourier space (rfft).

    Returns a complex array of shape (n_slots, n_modes) where n_modes is
    the rfft output length for the 1D case (n//2+1), or the product of
    rfft output lengths for multi-D.
    """
    n_slots = layout.num_slots
    n_pts = layout.num_points
    shape = grid.shape

    # For 1D: rfft output length is shape[0]//2 + 1
    # For nD: rfftn produces shape[:-1] + (shape[-1]//2+1,)
    # Compute analytically instead of probing with a zero FFT.
    rfft_shape = list(shape)
    rfft_shape[-1] = shape[-1] // 2 + 1
    n_modes = int(np.prod(rfft_shape))

    y_hat = np.zeros((n_slots, n_modes), dtype=np.complex128)
    for slot_idx in range(n_slots):
        start = slot_idx * n_pts
        end = start + n_pts
        field_data = y[start:end].reshape(shape)
        y_hat[slot_idx] = np.fft.rfftn(field_data).ravel()

    return y_hat


def _ifft_slots(
    y_hat: NDArray[np.complex128],
    layout: StateLayout,
    grid: GridInfo,
) -> NDArray[np.float64]:
    """Transform each slot from Fourier space back to physical space (irfft).

    Returns a real flat array of shape (n_slots * n_pts,).
    """
    n_slots = layout.num_slots
    n_pts = layout.num_points
    shape = grid.shape

    # Determine rfftn output shape for irfftn reconstruction
    rfft_shape = list(shape)
    rfft_shape[-1] = shape[-1] // 2 + 1
    rfft_shape_tuple = tuple(rfft_shape)

    y_out = np.zeros(n_slots * n_pts)
    for slot_idx in range(n_slots):
        hat_data = y_hat[slot_idx].reshape(rfft_shape_tuple)
        physical = np.fft.irfftn(hat_data, s=shape, axes=list(range(len(shape))))
        y_out[slot_idx * n_pts : (slot_idx + 1) * n_pts] = physical.ravel()

    return y_out


# ---------------------------------------------------------------------------
# Wavenumber grid construction
# ---------------------------------------------------------------------------


def _build_k_axes(grid: GridInfo) -> list[NDArray[np.float64]]:
    """Build wavenumber arrays for each spatial axis.

    Uses the same convention as operators.get_wavenumbers:
    k = 2π · rfftfreq(N, d=dx) for the last axis (rfft),
    k = 2π · fftfreq(N, d=dx) for all other axes (full fft).
    """
    k_axes: list[NDArray[np.float64]] = []
    ndim = grid.ndim

    for ax in range(ndim):
        n = grid.shape[ax]
        dx = grid.dx[ax]
        if ax == ndim - 1:
            # Last axis uses rfft (half-complex)
            k = get_wavenumbers(n, dx)
        else:
            # Other axes use full fft
            k = np.asarray(2.0 * np.pi * np.fft.fftfreq(n, d=dx), dtype=np.float64)
        k_axes.append(k)

    return k_axes


def _build_k_grid(
    k_axes: list[NDArray[np.float64]],
) -> list[NDArray[np.float64]]:
    """Build broadcasted k-grid arrays from per-axis wavenumbers.

    Returns a list of arrays, one per axis, each broadcastable to the
    full rfft output shape.
    """
    ndim = len(k_axes)
    k_grid: list[NDArray[np.float64]] = []

    for ax in range(ndim):
        shape = [1] * ndim
        shape[ax] = len(k_axes[ax])
        k_grid.append(k_axes[ax].reshape(shape))

    return k_grid


# ---------------------------------------------------------------------------
# Constraint elimination (Fourier Schur complement)
# ---------------------------------------------------------------------------

# Ref: Hairer & Wanner (1996), Solving ODEs II, Ch. VII — DAE reduction.
# Ref: Ascher & Petzold (1998), Computer Methods for ODEs/DAEs, §10.2.


def _constraints_fourier_eliminable(
    spec: EquationSystem,
    constraint_eqs: Sequence[ComponentEquation],
) -> bool:
    """Check if all constraint equations can be eliminated in Fourier space.

    Requirements:
    - Each constraint operator must be decomposable (spatial x time)
    - No time-dependent coefficients in constraints

    Constraints may contain acceleration operators (mixed_T2_S1x, d2_t)
    which are handled by substituting the dynamical equations of motion
    before Schur elimination.
    """
    for eq in constraint_eqs:
        for term in eq.rhs_terms:
            if term.operator not in _OPERATOR_DECOMP:
                return False
            if term.time_dependent:
                return False
    return True


# ---------------------------------------------------------------------------
# Generalized mass-matrix evolution (M·ẍ = K·x + D·ẋ + J·x⃛)
# ---------------------------------------------------------------------------
# For systems with implicit acceleration coupling (d2_t, mixed_T2_S*) and
# jerk coupling (d3_t, mixed_T3_S*).  The mass matrix M may be singular,
# creating hidden constraints analogous to time_order=0 fields.
#
# Algorithm:
#   1. Build M, D, K, J matrices from operator decomposition
#   2. Eigendecompose M per mode — zero eigenvalues → constraints
#   3. Schur-eliminate mass-matrix constraints (same as constraint fields)
#   4. Substitute jerk terms using equations of motion
#   5. Build first-order evolution matrix A = [[0,I],[M⁻¹K, M⁻¹D]]
#   6. Combine with existing constraint field Schur elimination
#
# References:
#   Golub & Van Loan (2013), Matrix Computations §7.7 (generalized eigenvalue)
#   Hairer & Lubich (2003), ZAMM 83(1) (mass matrices in dynamics)
#   Ostrogradsky (1850), Mem. Acad. St. Petersbourg VI 4, 385


def _has_time_derivative_operators(spec: EquationSystem) -> bool:
    """Check whether any equation has time-derivative operators on its RHS."""
    for eq in spec.equations:
        for term in eq.rhs_terms:
            decomp = _OPERATOR_DECOMP.get(term.operator)
            if decomp is not None and decomp.time_order > 0:
                return True
    return False


def _build_evolution_matrices(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: object,  # CoefficientEvaluator
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
) -> tuple[
    NDArray[np.complex128],  # A_rhs (n_modes, n_dyn_slots, n_dyn_slots)
    NDArray[np.complex128] | None,  # B_lhs (n_modes, n_dyn, n_dyn) or None
    NDArray[np.complex128],  # recovery (n_modes, n_total_constraints, n_dyn_slots)
    NDArray[np.complex128] | None,  # v_recovery (n_modes, n_c, n_dyn) or None
    list[str],  # all constraint field names
    dict[int, int],  # orig_to_reduced slot mapping
]:
    """Build per-mode matrices for systems with mass-matrix coupling.

    Returns A_rhs and optionally B_lhs for the generalized eigenvalue
    problem B·d' = A·d, where B = I - vel_coupling may be singular
    (gauge freedom from circular constraint velocity dependencies).

    When B_lhs is not None, the caller should use scipy.linalg.eig(A, B)
    (QZ decomposition) instead of np.linalg.eig(A) for eigendecomposition.
    Infinite eigenvalues correspond to gauge-constrained directions.

    Handles the generalized second-order system:

        M(k)·ẍ = K(k)·x + D(k)·ẋ + J(k)·x⃛

    where M may be singular (creating hidden algebraic constraints) and J
    encodes jerk coupling from d3_t/mixed_T3_S* operators.

    The algorithm:
    1. Separates constraint (time_order=0) and dynamical fields
    2. Builds M, D, K matrices for dynamical fields from operator decomposition
    3. SVD of M to detect singular directions (rank-deficient mass matrix)
    4. Schur-eliminates mass-null directions: solves K_cc · z_c = -K_cd · z_d
       for determinable constraint modes.  When K_cc is itself rank-deficient
       (e.g. CDT dark photon at critical kinetic-mixing point), uses
       ``np.linalg.pinv(K_cc)`` to give the minimum-norm solution — determinable
       modes are solved exactly, undetermined modes are frozen at zero.
    5. Back-projects to original field basis via V_eff = V_d + V_c · recovery,
       where V_d, V_c are right singular vectors of M and recovery is the Schur
       constraint-recovery matrix.  The effective operators are
       K_orig = V_eff · E · V_eff⁺ (pseudoinverse), not V_d · E · V_dᴴ (#260).
    6. Substitutes jerk terms using the (now-invertible) dynamical equations
    7. Combines both constraint levels and builds the first-order evolution matrix
    8. Pre-solves ``B_lhs · A_final = A_rhs`` via ``np.linalg.solve`` + ``lstsq``
       fallback, returning a single first-order evolution matrix (``B_lhs=None``).
       This replaces the old generalized-eigenvalue output after #256.

    Returns
    -------
    A_final : ndarray, shape (n_modes, n_dyn, n_dyn)
        First-order evolution matrix (pre-solved, no separate ``B_lhs``).
    B_lhs : None
        Always ``None`` after #256; retained in the tuple for API compat
        with callers that destructure all six return values.
    recovery : ndarray
        Schur-recovery matrix for algebraic constraint fields.
    v_recovery : ndarray or None
        Exact constraint-velocity recovery (``recovery @ A_final``).
    constraint_field_names : list[str]
        Names of the algebraic constraint fields.
    orig_to_reduced : dict[int, int]
        Original → reduced slot index mapping.
    """
    import logging  # noqa: PLC0415

    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    assert isinstance(coeff_eval, CoefficientEvaluator)
    logger = logging.getLogger(__name__)

    n_modes = int(np.prod(rfft_shape))

    # ---- Identify constraint and dynamical fields ----
    constraint_field_names: list[str] = [
        eq.field_name for eq in spec.equations if eq.time_derivative_order == 0
    ]
    c_idx_map: dict[str, int] = {
        name: i for i, name in enumerate(constraint_field_names)
    }
    n_c = len(constraint_field_names)

    # Build dynamical-only slot mapping (excluding constraint field slots)
    orig_to_reduced: dict[int, int] = {}
    red_idx = 0
    for si, slot in enumerate(layout.slots):
        if slot.kind == "constraint":
            continue
        orig_to_reduced[si] = red_idx
        red_idx += 1
    n_dyn_slots = red_idx

    # Map field/velocity names → reduced slot indices
    dyn_slot_map: dict[str, int] = {}
    for si, slot in enumerate(layout.slots):
        if si in orig_to_reduced:
            dyn_slot_map[slot.name] = orig_to_reduced[si]
    for fname, si in layout.velocity_slot_map.items():
        v_name = f"v_{fname}"
        if si in orig_to_reduced:
            dyn_slot_map[v_name] = orig_to_reduced[si]

    # ---- Evaluate spatial Fourier multipliers ----
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op not in multiplier_cache:
                decomp = _OPERATOR_DECOMP[op]
                mult_val = decomp.spatial_fn(k_grid)
                mult_full = np.broadcast_to(mult_val, rfft_shape)
                multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    # ---- Identify dynamical fields and their indices ----
    # Map: dynamical field name → index in the n_f dynamical field array
    dyn_field_names: list[str] = []
    dyn_field_idx: dict[str, int] = {}
    for eq in spec.equations:
        if eq.time_derivative_order > 0:
            dyn_field_idx[eq.field_name] = len(dyn_field_names)
            dyn_field_names.append(eq.field_name)
    n_f = len(dyn_field_names)  # number of dynamical FIELDS (not slots)

    # ---- Build M, D, K matrices for dynamical fields (n_f x n_f) ----
    # These are the FIELD-level matrices, not slot-level.
    # M·ẍ = K·x + D·ẋ where x is the vector of field values.
    M_mat = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)
    D_mat = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)
    K_mat = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)
    J_mat = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)

    # Diagonal of M: each 2nd-order equation has kinetic_coefficient · ẍ_i
    # on the LHS.  Previously hardcoded to 1.0, which was wrong for equations
    # whose Wolfram-exported kinetic coefficient is symbolic (e.g. -kappa^(-2)
    # for gravitons, -xi for massive vector torsion).  The CLI workaround
    # `normalize_kinetic_coefficients` divided each RHS by its own kinetic
    # coefficient so the assumption `diagonal=1` held, but that introduced
    # asymmetric off-diagonals under any Schur elimination that assumes
    # symmetric M (np.linalg.eigh).  Root fix: read the kinetic coefficient
    # directly into the diagonal, and handle asymmetry via SVD below.
    from tidal.symbolic._eval_utils import (  # noqa: PLC0415
        evaluate_coefficient,
    )

    eq_for_field: dict[str, object] = {
        eq.field_name: eq for eq in spec.equations if eq.time_derivative_order > 0
    }
    for fi, fname in enumerate(dyn_field_names):
        eq_f = eq_for_field[fname]
        kin_sym = getattr(eq_f, "kinetic_coefficient_symbolic", None)
        if kin_sym is None:
            M_mat[:, fi, fi] = 1.0
        else:
            try:
                kin_val = evaluate_coefficient(
                    kin_sym,
                    coeff_eval._parameters,  # noqa: SLF001
                    spec.effective_coordinates,
                )
            except Exception:  # noqa: BLE001
                # Fall back to the old hardcoded value if resolution fails.
                # Matches the behaviour of normalize_kinetic_coefficients
                # which skips normalisation in this case.
                M_mat[:, fi, fi] = 1.0
                continue
            if isinstance(kin_val, np.ndarray):
                # Position-dependent kinetic — broadcast to per-mode shape.
                # (Rare; most theories have constant kin coefficients.)
                M_mat[:, fi, fi] = complex(kin_val.ravel()[0])
            else:
                M_mat[:, fi, fi] = complex(float(kin_val))

    # Constraint matrices — built in two phases:
    #   Phase 1: collect terms from constraint equations
    #   Phase 2: substitute acceleration/velocity terms after M inversion
    S_cd = np.zeros((n_modes, n_c, n_dyn_slots), dtype=np.complex128)
    S_cc = np.zeros((n_modes, n_c, n_c), dtype=np.complex128)
    # A_dc: constraint → dynamical (field + velocity references)
    A_dc_field = np.zeros((n_modes, n_dyn_slots, n_c), dtype=np.complex128)
    A_dc_vel = np.zeros((n_modes, n_dyn_slots, n_c), dtype=np.complex128)

    # Deferred constraint terms with time_order > 0 on dynamical fields.
    # These need acceleration/velocity substitution after M inversion.
    # Each entry: (ci, coeff, spatial_mult, time_order, target_field_idx)
    deferred_constraint_terms: list[
        tuple[int, complex, NDArray[np.complex128], int, int]
    ] = []

    # ---- Populate matrices from equations ----
    for eq_idx, eq in enumerate(spec.equations):
        is_constraint = eq.time_derivative_order == 0

        if is_constraint:
            ci = c_idx_map[eq.field_name]
            for term_idx, term in enumerate(eq.rhs_terms):
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx
                )
                mult = multiplier_cache[term.operator]
                decomp = _OPERATOR_DECOMP[term.operator]
                t_order = decomp.time_order

                if term.field in c_idx_map:
                    cj = c_idx_map[term.field]
                    S_cc[:, ci, cj] += coeff * mult
                elif t_order == 0:
                    # Pure spatial operator on dynamical field/velocity
                    if term.field in dyn_slot_map:
                        dj = dyn_slot_map[term.field]
                        S_cd[:, ci, dj] += coeff * mult
                elif t_order == 1 and term.field in dyn_field_idx:
                    # Velocity operator on dynamical field (e.g. mixed_T1_S1x)
                    # This references ẋ_field → use velocity slot
                    fj = dyn_field_idx[term.field]
                    vel_j = orig_to_reduced[layout.velocity_slot_map[term.field]]
                    S_cd[:, ci, vel_j] += coeff * mult
                elif t_order >= 2 and term.field in dyn_field_idx:
                    # Acceleration/jerk on dynamical field — defer until M inverted
                    fj = dyn_field_idx[term.field]
                    deferred_constraint_terms.append(
                        (ci, complex(coeff), mult, t_order, fj)
                    )
                elif term.field in dyn_slot_map:
                    # Fallback: direct slot reference
                    dj = dyn_slot_map[term.field]
                    S_cd[:, ci, dj] += coeff * mult
            continue

        # Dynamical equation
        fi = dyn_field_idx[eq.field_name]
        field_slot = orig_to_reduced[layout.field_slot_map[eq.field_name]]
        vel_slot = orig_to_reduced[layout.velocity_slot_map[eq.field_name]]

        for term_idx, term in enumerate(eq.rhs_terms):
            coeff = _resolve_constant_coeff(
                term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx
            )
            mult = multiplier_cache[term.operator]
            decomp = _OPERATOR_DECOMP[term.operator]
            t_order = decomp.time_order

            # Determine which field this term targets
            target_field = term.field
            # Strip v_ prefix to get base field name for velocity references
            is_vel_ref = target_field.startswith("v_")
            base_field = target_field[2:] if is_vel_ref else target_field

            if target_field in c_idx_map:
                # Direct reference to constraint field
                cj = c_idx_map[target_field]
                A_dc_field[:, vel_slot, cj] += coeff * mult
            elif is_vel_ref and base_field in c_idx_map:
                # Velocity of constraint field
                cj = c_idx_map[base_field]
                A_dc_vel[:, vel_slot, cj] += coeff * mult
            elif base_field in dyn_field_idx:
                fj = dyn_field_idx[base_field]
                if t_order == 0:
                    if is_vel_ref:
                        # Velocity reference with spatial operator
                        # → damping matrix D[fi, fj]
                        D_mat[:, fi, fj] += coeff * mult
                    else:
                        # Position reference with spatial operator
                        # → stiffness matrix K[fi, fj]
                        K_mat[:, fi, fj] += coeff * mult
                elif t_order == 1:
                    if is_vel_ref:
                        # first_derivative_t(v_X) = ẍ_X → acceleration
                        # This should be rare; treat as M coupling
                        M_mat[:, fi, fj] -= coeff * mult
                    else:
                        # first_derivative_t(X) = ẋ_X → velocity
                        D_mat[:, fi, fj] += coeff * mult
                elif t_order == 2:
                    # d2_t or mixed_T2: acceleration coupling → mass matrix
                    # RHS has coeff·ẍ_j, move to LHS: M[fi,fj] -= coeff·mult
                    M_mat[:, fi, fj] -= coeff * mult
                elif t_order == 3:
                    # d3_t or mixed_T3: jerk coupling → substitute later
                    J_mat[:, fi, fj] += coeff * mult
            elif target_field in dyn_slot_map:
                # Direct slot reference (velocity name like v_h_3)
                dj = dyn_slot_map[target_field]
                # Just put it in the A matrix directly later
                # For now, track separately if needed

    # ---- Mass-matrix constraint elimination ----
    # Eigendecompose M per mode to find singular directions.
    # Zero eigenvalues → hidden constraints; nonzero → dynamical.
    #
    # For each mode k:
    #   M(k) = Q(k) · Λ(k) · Q(k)ᵀ
    #   Rotate: K̃ = QᵀKQ, D̃ = QᵀDQ
    #   Singular rows (λ=0) → constraint: 0 = K̃_c·z + D̃_c·ż
    #   Dynamical rows (λ≠0) → ODE: Λ_d·z̈ = K̃_d·z + D̃_d·ż

    # Check if any mode has singular M
    dets = np.linalg.det(M_mat)
    has_singular_M = np.any(np.abs(dets) < 1e-12)
    m_k_independent = False
    con_mask = np.zeros(n_f, dtype=bool)  # will be updated if singular

    if has_singular_M:
        logger.info(
            "Generalized mass matrix: singular M detected — "
            "applying mass-matrix Schur elimination"
        )

    # Build the first-order evolution matrix A in the FULL dynamical slot space.
    # A has shape (n_modes, n_dyn_slots, n_dyn_slots).
    # For 2nd-order fields: rows for field_slot get dq/dt = v (kinematic),
    # rows for vel_slot get dv/dt = M⁻¹(K·x + D·v).
    A_dd = np.zeros((n_modes, n_dyn_slots, n_dyn_slots), dtype=np.complex128)

    # Kinematic equations: dq/dt = v
    for fname in dyn_field_names:
        field_slot = orig_to_reduced[layout.field_slot_map[fname]]
        vel_slot = orig_to_reduced[layout.velocity_slot_map[fname]]
        A_dd[:, field_slot, vel_slot] = 1.0

    if has_singular_M:
        # Use eigendecomposition to handle singular M
        # We work with each mode separately for modes where M is singular,
        # and batch-process modes where M is invertible.

        # For simplicity and correctness, process per-mode where needed.
        # M is typically k-independent for d2_t coupling (spatial_mult=1),
        # so use the k=0 mode's eigenstructure as representative.
        # For k-dependent M (from mixed_T2_S*), process per mode.

        # Check if M is k-independent
        M_spread = np.max(np.abs(M_mat - M_mat[0:1, :, :]))
        m_k_independent = M_spread < 1e-14

        if m_k_independent:
            # M is the same for all modes — single SVD-based Schur.
            #
            # The system is M·ẍ = K·x + D·ẋ + J·x⃛.  When M is singular
            # (rank-deficient), the naive Moore-Penrose pseudoinverse
            # M⁺ · (K·x + D·ẋ) gives the minimum-norm ẍ but misses the
            # **algebraic constraints** encoded in the left null space
            # of M.  For the torsion dark-photon model, the rank-3
            # tensor produces a rank-deficient M where multiple equations
            # share the same mass-matrix row up to sign; these rows
            # impose algebraic constraints on x (and ẋ) via their
            # non-proportional K and D rows.
            #
            # Full Schur elimination via SVD:
            #   M = U·Σ·Vᵀ  (asymmetric-safe, unlike eigh)
            #   Left null space = columns of U with s=0 -- row constraints
            #   Right null space = columns of V with s=0 -- gauge modes
            #   Rotate K, D, J as K~ = Ut*K*V, etc.
            #   Partition by dyn/con masks: dynamical rows (s>0)
            #     evolve, constraint rows (s=0) give algebraic
            #     constraints on the gauge column variables.
            #   Solve K̃_cc·z_c = -K̃_cd·z_d for the gauge modes.
            #   Substitute back into the dynamical block to get
            #     K_eff, D_eff, then diagonal inversion by Σ_d.
            #   Back-project via V_d to the original basis.
            #
            # This replaces the previous np.linalg.eigh-based approach,
            # which assumed M was real-symmetric.  The derived-JSON M
            # is NOT symmetric in general (each equation has its own
            # LHS kinetic coefficient, and the off-diagonal d2_t
            # cross-couplings between distinct fields don't
            # antisymmetrise), so eigh silently symmetrized via
            # (M + Mᵀ)/2 and gave a wrong eigenbasis for rank-deficient
            # cases (leading to the torsion dark-photon null result).
            #
            # Ref: Golub & Van Loan (2013), Matrix Computations §5.5
            # (generalized Schur decomposition for DAE systems).
            M0 = M_mat[0]
            U_svd, S_svd, Vh_svd = np.linalg.svd(M0)
            tol = 1e-10 * max(1.0, float(np.max(S_svd)) if S_svd.size else 1.0)
            dyn_mask = S_svd > tol
            con_mask = ~dyn_mask
            n_mass_con = int(np.sum(con_mask))
            n_mass_dyn = int(np.sum(dyn_mask))

            logger.info(
                "Mass matrix singular values: %s (dynamical: %d, constrained: %d)",
                S_svd,
                n_mass_dyn,
                n_mass_con,
            )

            if n_mass_con > 0:
                # U = left singular vectors (rows of M -> row basis)
                # V = right singular vectors (columns of M -> col basis)
                # Σ_d = diag of positive singular values
                V_full = Vh_svd.conj().T  # (n_f, n_f)
                U_full = U_svd  # (n_f, n_f)
                V_d = V_full[:, dyn_mask]  # (n_f, n_mass_dyn)
                V_c = V_full[:, con_mask]  # (n_f, n_mass_con)
                Sigma_d_inv = np.diag(1.0 / S_svd[dyn_mask])  # (n_mass_dyn, n_mass_dyn)

                # Rotate K, D, J: K̃ = Uᵀ · K · V
                K_rot = np.einsum("ij,mjk,kl->mil", U_full.conj().T, K_mat, V_full)
                D_rot = np.einsum("ij,mjk,kl->mil", U_full.conj().T, D_mat, V_full)
                J_rot = np.einsum("ij,mjk,kl->mil", U_full.conj().T, J_mat, V_full)

                # Partition into dyn/con blocks
                d_idx = np.where(dyn_mask)[0]
                c_idx = np.where(con_mask)[0]
                K_dd = K_rot[:, np.ix_(d_idx, d_idx)[0], np.ix_(d_idx, d_idx)[1]]
                K_dc = K_rot[:, np.ix_(d_idx, c_idx)[0], np.ix_(d_idx, c_idx)[1]]
                K_cd = K_rot[:, np.ix_(c_idx, d_idx)[0], np.ix_(c_idx, d_idx)[1]]
                K_cc = K_rot[:, np.ix_(c_idx, c_idx)[0], np.ix_(c_idx, c_idx)[1]]
                D_dd = D_rot[:, np.ix_(d_idx, d_idx)[0], np.ix_(d_idx, d_idx)[1]]
                D_dc = D_rot[:, np.ix_(d_idx, c_idx)[0], np.ix_(d_idx, c_idx)[1]]
                J_dd = J_rot[:, np.ix_(d_idx, d_idx)[0], np.ix_(d_idx, d_idx)[1]]
                J_dc = J_rot[:, np.ix_(d_idx, c_idx)[0], np.ix_(d_idx, c_idx)[1]]

                # Constraint rows (s=0): 0 = K_cd*z_d + K_cc*z_c + D_cd*dz_d/dt + ...
                # Solve K_cc · z_c = -K_cd · z_d (position-only constraints).
                #
                # K_cc may itself be rank-deficient (e.g. CDT dark photon at
                # dm = √ξ/2 where extra mass-matrix null directions open).
                # Use pseudoinverse instead of regularized inverse: pinv gives
                # the minimum-norm solution for determinable constraint modes
                # and zeroes undetermined ones.  For full-rank K_cc, pinv ≡ inv.
                K_cc_norms = np.linalg.norm(K_cc, axis=(1, 2))
                has_k_con = np.any(K_cc_norms > 1e-14)

                if has_k_con:
                    K_cc_pinv = np.linalg.pinv(K_cc)
                    mass_recovery = -np.einsum("mij,mjk->mik", K_cc_pinv, K_cd)
                    K_eff = K_dd + np.einsum("mij,mjk->mik", K_dc, mass_recovery)
                    D_eff = D_dd + np.einsum("mij,mjk->mik", D_dc, mass_recovery)
                    J_eff = J_dd + np.einsum("mij,mjk->mik", J_dc, mass_recovery)
                else:
                    # Null-space decouples trivially from the rest —
                    # keep only the dyn/dyn block.
                    K_eff = K_dd
                    D_eff = D_dd
                    J_eff = J_dd

                # Invert Σ_d (diagonal): E = Σ_d⁻¹·K_eff, F = Σ_d⁻¹·D_eff
                E = np.einsum("ij,mjk->mik", Sigma_d_inv, K_eff)
                F = np.einsum("ij,mjk->mik", Sigma_d_inv, D_eff)

                # Jerk substitution: when J ≠ 0, the d3_t(x_j) term on
                # the RHS couples third-order time derivatives. Reduce
                # to second-order by iterating the ẍ expression once.
                J_eff_inv = np.einsum("ij,mjk->mik", Sigma_d_inv, J_eff)
                has_jerk = float(np.max(np.abs(J_eff_inv))) > 1e-15
                if has_jerk:
                    logger.info("Jerk substitution: applying d3_t elimination")
                    FE = np.einsum("mij,mjk->mik", F, E)
                    K_jerk = np.einsum("mij,mjk->mik", J_eff_inv, FE)
                    FF = np.einsum("mij,mjk->mik", F, F)
                    D_jerk = np.einsum("mij,mjk->mik", J_eff_inv, E + FF)
                    E_final = E + K_jerk
                    F_final = F + D_jerk
                else:
                    E_final = E
                    F_final = F

                # Back-project to original basis.  After Schur
                # elimination, z_c = recovery · z_d, so the original
                # coordinates are:
                #   x = V_d · z_d + V_c · z_c
                #     = (V_d + V_c · recovery) · z_d  =  V_eff · z_d
                # The effective operators in original field space are:
                #   K_orig = V_eff · E_final · V_eff⁺
                #   D_orig = V_eff · F_final · V_eff⁺
                # where V_eff⁺ = (V_effᴴ V_eff)⁻¹ V_effᴴ is the
                # left pseudoinverse.  When recovery = 0 (trivially
                # decoupled constraints), V_eff = V_d and V_eff⁺ = V_dᴴ,
                # recovering the previous formula.
                if has_k_con:
                    # mass_recovery: (n_modes, n_mass_con, n_mass_dyn)
                    V_eff = V_d[np.newaxis, :, :] + np.einsum(
                        "ic,mcj->mij", V_c, mass_recovery
                    )  # (n_modes, n_f, n_mass_dyn)
                    # V_eff^+ = pinv(V_eff) handles rank-deficient cases
                    # (e.g. mA2 near zero where mass modes become degenerate).
                    # Using pinv(V_eff) directly is more numerically stable than
                    # the two-step (V_eff^H V_eff)^{-1} V_eff^H formula, which
                    # requires inv(V_eff^H V_eff) and fails when V_eff has
                    # linearly-dependent columns.
                    V_eff_pinv = np.linalg.pinv(V_eff)  # (n_modes, n_mass_dyn, n_f)
                    K_orig = np.einsum("mia,mab,mbj->mij", V_eff, E_final, V_eff_pinv)
                    D_orig = np.einsum("mia,mab,mbj->mij", V_eff, F_final, V_eff_pinv)
                else:
                    # Trivially decoupled: V_eff = V_d, V_eff⁺ = V_dᴴ
                    K_orig = np.einsum("ia,mab,jb->mij", V_d, E_final, V_d.conj())
                    D_orig = np.einsum("ia,mab,jb->mij", V_d, F_final, V_d.conj())

                # Fill A_dd velocity rows
                for i, fname_i in enumerate(dyn_field_names):
                    vel_i = orig_to_reduced[layout.velocity_slot_map[fname_i]]
                    for j, fname_j in enumerate(dyn_field_names):
                        field_j = orig_to_reduced[layout.field_slot_map[fname_j]]
                        vel_j = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                        A_dd[:, vel_i, field_j] += K_orig[:, i, j]
                        A_dd[:, vel_i, vel_j] += D_orig[:, i, j]
            else:
                # No singular directions — M is invertible
                m_inv = np.linalg.inv(M_mat)
                eff_k = np.einsum("mij,mjk->mik", m_inv, K_mat)
                eff_d = np.einsum("mij,mjk->mik", m_inv, D_mat)

                # Jerk substitution
                j_inv = np.einsum("mij,mjk->mik", m_inv, J_mat)
                has_jerk = np.max(np.abs(j_inv)) > 1e-15
                if has_jerk:
                    fd_k = np.einsum("mij,mjk->mik", eff_d, eff_k)
                    k_jerk = np.einsum("mij,mjk->mik", j_inv, fd_k)
                    fd_d = np.einsum("mij,mjk->mik", eff_d, eff_d)
                    d_jerk = np.einsum("mij,mjk->mik", j_inv, eff_k + fd_d)
                    eff_k += k_jerk
                    eff_d += d_jerk

                for i, fname_i in enumerate(dyn_field_names):
                    vel_i = orig_to_reduced[layout.velocity_slot_map[fname_i]]
                    for j, fname_j in enumerate(dyn_field_names):
                        field_j = orig_to_reduced[layout.field_slot_map[fname_j]]
                        vel_j = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                        A_dd[:, vel_i, field_j] += eff_k[:, i, j]
                        A_dd[:, vel_i, vel_j] += eff_d[:, i, j]
        else:
            # M is k-dependent — process per mode
            # For now, treat each mode independently
            for m in range(n_modes):
                M_m = M_mat[m]
                # Use SVD (not eigh) for rank detection since M may be
                # asymmetric per the non-uniform kinetic-coefficient
                # convention; eigh would silently symmetrise it.
                sv_m = np.linalg.svd(M_m, compute_uv=False)
                tol = 1e-10 * max(1.0, float(np.max(sv_m)) if sv_m.size else 1.0)
                dyn_m = sv_m > tol
                if np.all(dyn_m):
                    # Invertible for this mode
                    m_inv_m = np.linalg.inv(M_m)  # pyright: ignore[reportUnknownVariableType]
                    ek_m = m_inv_m @ K_mat[m]  # pyright: ignore[reportUnknownVariableType]
                    ed_m = m_inv_m @ D_mat[m]  # pyright: ignore[reportUnknownVariableType]
                    j_inv_m = m_inv_m @ J_mat[m]  # pyright: ignore[reportUnknownVariableType]
                    if np.max(np.abs(j_inv_m)) > 1e-15:  # pyright: ignore[reportUnknownArgumentType]
                        fd_k_m = ed_m @ ek_m  # pyright: ignore[reportUnknownVariableType]
                        ek_m += j_inv_m @ fd_k_m  # pyright: ignore[reportUnknownVariableType]
                        ed_m += j_inv_m @ (ek_m + ed_m @ ed_m)  # pyright: ignore[reportUnknownVariableType]
                    for i, fname_i in enumerate(dyn_field_names):
                        vi = orig_to_reduced[layout.velocity_slot_map[fname_i]]
                        for j, fname_j in enumerate(dyn_field_names):
                            fj = orig_to_reduced[layout.field_slot_map[fname_j]]
                            vj = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                            A_dd[m, vi, fj] += ek_m[i, j]
                            A_dd[m, vi, vj] += ed_m[i, j]
                else:
                    # Singular mode — would need per-mode Schur elimination
                    # This is rare for k-dependent M; log and use pseudoinverse
                    m_pinv = np.linalg.pinv(M_m)  # pyright: ignore[reportUnknownVariableType]
                    ek_m2 = m_pinv @ K_mat[m]  # pyright: ignore[reportUnknownVariableType]
                    ed_m2 = m_pinv @ D_mat[m]  # pyright: ignore[reportUnknownVariableType]
                    for i, fname_i in enumerate(dyn_field_names):
                        vi = orig_to_reduced[layout.velocity_slot_map[fname_i]]
                        for j, fname_j in enumerate(dyn_field_names):
                            fj = orig_to_reduced[layout.field_slot_map[fname_j]]
                            vj = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                            A_dd[m, vi, fj] += ek_m2[i, j]
                            A_dd[m, vi, vj] += ed_m2[i, j]
    else:
        # M is invertible for all modes — standard path
        m_inv = np.linalg.inv(M_mat)
        eff_k = np.einsum("mij,mjk->mik", m_inv, K_mat)
        eff_d = np.einsum("mij,mjk->mik", m_inv, D_mat)

        # Jerk substitution
        j_inv = np.einsum("mij,mjk->mik", m_inv, J_mat)
        has_jerk = np.max(np.abs(j_inv)) > 1e-15
        if has_jerk:
            logger.info("Jerk substitution: applying d3_t elimination")
            fd_k = np.einsum("mij,mjk->mik", eff_d, eff_k)
            k_jerk = np.einsum("mij,mjk->mik", j_inv, fd_k)
            fd_d = np.einsum("mij,mjk->mik", eff_d, eff_d)
            d_jerk = np.einsum("mij,mjk->mik", j_inv, eff_k + fd_d)
            eff_k += k_jerk
            eff_d += d_jerk

        for i, fname_i in enumerate(dyn_field_names):
            vel_i = orig_to_reduced[layout.velocity_slot_map[fname_i]]
            for j, fname_j in enumerate(dyn_field_names):
                field_j = orig_to_reduced[layout.field_slot_map[fname_j]]
                vel_j = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                A_dd[:, vel_i, field_j] += eff_k[:, i, j]
                A_dd[:, vel_i, vel_j] += eff_d[:, i, j]

    # ---- Substitute deferred constraint acceleration/velocity terms ----
    # Constraints may contain time_order>=2 operators on dynamical fields
    # (e.g., mixed_T2_S1x(t_3) = ik_x x ẍ_{t_3}).  After mass-matrix
    # inversion, ẍ_j = Σ_k E[j,k]·field_k + F[j,k]·vel_k.  Substitute
    # this into the constraint's S_cd matrix.
    if deferred_constraint_terms:
        # Extract effective acceleration matrices from A_dd.
        # A_dd[m, vel_i, field_j] = K_eff[i,j] (position → acceleration)
        # A_dd[m, vel_i, vel_j] = D_eff[i,j] (velocity → acceleration)
        K_eff = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)
        D_eff = np.zeros((n_modes, n_f, n_f), dtype=np.complex128)
        for i, fname_i in enumerate(dyn_field_names):
            vel_i = orig_to_reduced[layout.velocity_slot_map[fname_i]]
            for j, fname_j in enumerate(dyn_field_names):
                field_j = orig_to_reduced[layout.field_slot_map[fname_j]]
                vel_j = orig_to_reduced[layout.velocity_slot_map[fname_j]]
                K_eff[:, i, j] = A_dd[:, vel_i, field_j]
                D_eff[:, i, j] = A_dd[:, vel_i, vel_j]

        for ci, coeff_val, spatial_mult, t_order, fj in deferred_constraint_terms:
            if t_order == 2:
                # ẍ_fj = Σ_k K_eff[fj,k]·field_k + D_eff[fj,k]·vel_k
                for k, fname_k in enumerate(dyn_field_names):
                    fk_slot = orig_to_reduced[layout.field_slot_map[fname_k]]
                    vk_slot = orig_to_reduced[layout.velocity_slot_map[fname_k]]
                    # Position contribution: coeff x spatial x K_eff[fj, k]
                    S_cd[:, ci, fk_slot] += coeff_val * spatial_mult * K_eff[:, fj, k]
                    # Velocity contribution: coeff x spatial x D_eff[fj, k]
                    S_cd[:, ci, vk_slot] += coeff_val * spatial_mult * D_eff[:, fj, k]
            # time_order=3 in constraints is very rare; log and skip
            elif t_order >= 3:
                logger.warning(
                    "Constraint has time_order=%d operator — not yet handled",
                    t_order,
                )

    # ---- Constraint field Schur elimination ----
    if n_c > 0:
        # Batch-invert S_cc
        cc_dets = np.linalg.det(S_cc) if n_c > 0 else np.ones(n_modes)
        singular_mask = np.abs(cc_dets) < 1e-14
        S_cc_reg = S_cc.copy()
        if np.any(singular_mask):
            S_cc_reg[singular_mask] += 1e-14 * np.eye(n_c, dtype=np.complex128)
        S_cc_inv = np.linalg.inv(S_cc_reg)

        # Recovery: c = -S_cc⁻¹ · S_cd · d
        recovery = -np.einsum("mij,mjk->mik", S_cc_inv, S_cd)

        # Field correction: A_dc_field · recovery
        field_correction = np.einsum("mij,mjk->mik", A_dc_field, recovery)

        # Velocity coupling: A_dc_vel · recovery
        vel_coupling = np.einsum("mij,mjk->mik", A_dc_vel, recovery)
        has_vel = np.max(np.abs(vel_coupling)) > 1e-15

        A_rhs = A_dd + field_correction

        if has_vel:
            eye = np.broadcast_to(
                np.eye(n_dyn_slots, dtype=np.complex128),
                (n_modes, n_dyn_slots, n_dyn_slots),
            ).copy()
            B_lhs: NDArray[np.complex128] | None = eye - vel_coupling
        else:
            B_lhs = None
    else:
        recovery = np.zeros((n_modes, 0, n_dyn_slots), dtype=np.complex128)
        A_rhs = A_dd
        B_lhs = None

    n_mass_con_total = int(np.sum(con_mask))
    logger.info(
        "Generalized evolution: %d constraint fields, %d mass-matrix constraints, "
        "%d dynamical slots, jerk=%s, vel_coupling=%s",
        n_c,
        n_mass_con_total,
        n_dyn_slots,
        "yes" if np.max(np.abs(J_mat)) > 1e-15 else "no",
        "generalized_eig" if B_lhs is not None else "none",
    )

    # Velocity recovery for constraint fields: v_c = recovery · d' where
    # d' = B_lhs⁻¹ · A_rhs · d (computed per mode with lstsq fallback for
    # near-singular B).  This is needed only for the constraint-velocity
    # measurement; the main evolution uses generalized eig(A, B) in
    # _evolve_per_mode when B_lhs is not None, which handles rank-deficient
    # or near-singular B correctly via QZ decomposition.
    #
    # Why NOT pre-solve A_final = B⁻¹·A here and return a single first-order
    # matrix: for systems with rank-deficient mass matrix (d2_t cross-
    # couplings), the pre-solved A_final has nearly-degenerate eigenvalues
    # which make the eigendecomposition V ill-conditioned (cond(V) > 1e9).
    # Passing (A, B) to scipy.linalg.eig uses QZ decomposition, which
    # handles rank-deficiency cleanly via the Schur form without requiring
    # V to be well-conditioned.  See #256 for the empirical test case.
    if B_lhs is not None and recovery.size > 0:
        A_eff = np.zeros_like(A_rhs)
        fallback_count = 0
        for m in range(n_modes):
            try:
                A_eff[m] = np.linalg.solve(B_lhs[m], A_rhs[m])
            except np.linalg.LinAlgError:
                A_eff[m] = np.asarray(
                    np.linalg.lstsq(B_lhs[m], A_rhs[m], rcond=None)[0],  # type: ignore[reportUnknownMemberType]
                    dtype=np.complex128,
                )
                fallback_count += 1
        if fallback_count > 0:
            logger.info(
                "Constraint-velocity recovery: %d/%d modes used lstsq fallback",
                fallback_count,
                n_modes,
            )
        v_recovery = np.einsum("mci,mij->mcj", recovery, A_eff)
    elif recovery.size > 0:
        v_recovery = np.einsum("mci,mij->mcj", recovery, A_rhs)
    else:
        v_recovery = None

    return A_rhs, B_lhs, recovery, v_recovery, constraint_field_names, orig_to_reduced


# ---------------------------------------------------------------------------
# Evolution matrix construction
# ---------------------------------------------------------------------------


def _build_per_mode_matrices(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
) -> NDArray[np.complex128]:
    """Build evolution matrices for the all-constant-coefficient case.

    Returns array of shape (n_modes, n_state_slots, n_state_slots) where
    each [m, :, :] is the evolution matrix for mode m.

    The matrix has block structure:
        For second-order fields: [0, 1; L(k), 0]  (velocity coupling)
        For first-order fields:  [L(k)]            (direct evolution)
    """
    n_slots = layout.num_slots
    n_modes = int(np.prod(rfft_shape))

    def _resolve_target_slot(field_ref: str) -> int | None:
        """Resolve a field/velocity reference to a slot index."""
        if field_ref in layout.field_slot_map:
            return layout.field_slot_map[field_ref]
        if field_ref.startswith("v_") and field_ref[2:] in layout.velocity_slot_map:
            return layout.velocity_slot_map[field_ref[2:]]
        return None

    # Evaluate Fourier multipliers on the k-grid
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op not in multiplier_cache:
                mult_fn = _EXACT_MULTIPLIERS[op]
                mult_val = mult_fn(k_grid)
                # Broadcast to full rfft shape and flatten
                mult_full = np.broadcast_to(mult_val, rfft_shape)
                multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    # Build matrices: A[m, i, j] for each mode m
    A = np.zeros((n_modes, n_slots, n_slots), dtype=np.complex128)

    for _eq_idx, eq in enumerate(spec.equations):
        field_name = eq.field_name
        is_second_order = eq.time_derivative_order >= 2

        if is_second_order:
            # Field slot and velocity slot
            field_slot = layout.field_slot_map[field_name]
            vel_slot = layout.velocity_slot_map[field_name]

            # dq/dt = v  →  A[field_slot, vel_slot] = 1
            A[:, field_slot, vel_slot] = 1.0

            # dv/dt = Σ coeff * operator(target_field)
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = _resolve_target_slot(term.field)
                if target_slot is None:
                    continue
                coeff = _resolve_constant_coeff(
                    term,
                    coeff_eval,
                    eq_idx=_eq_idx,
                    term_idx=_term_idx,
                )
                mult = multiplier_cache[term.operator]
                A[:, vel_slot, target_slot] += coeff * mult

        else:
            # First-order: du/dt = Σ coeff * operator(target_field)
            this_slot = layout.field_slot_map[field_name]
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = _resolve_target_slot(term.field)
                if target_slot is None:
                    continue
                coeff = _resolve_constant_coeff(
                    term,
                    coeff_eval,
                    eq_idx=_eq_idx,
                    term_idx=_term_idx,
                )
                mult = multiplier_cache[term.operator]
                A[:, this_slot, target_slot] += coeff * mult

    return A


def _resolve_constant_coeff(
    term: OperatorTerm,
    coeff_eval: CoefficientEvaluator,
    *,
    eq_idx: int = -1,
    term_idx: int = -1,
) -> complex:
    """Resolve a constant (non-position-dependent) coefficient to a scalar.

    Uses CoefficientEvaluator.resolve() which returns a float for constant
    coefficients or an ndarray for position-dependent ones (the latter should
    not reach this function).
    """
    resolved = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
    if isinstance(resolved, np.ndarray):
        # Position-dependent — should not happen for constant-coeff path
        return complex(resolved.ravel()[0])
    return complex(resolved)


def _build_convolution_matrix(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
) -> NDArray[np.complex128]:
    """Build full evolution matrix for position-dependent coefficient case.

    Position-dependent coefficients c(x) create convolution coupling in
    k-space: FFT[c(x)·u(x)] = ĉ * û (convolution).  This couples
    different k-modes, producing a full (n_total x n_total) matrix where
    n_total = n_slots x n_modes.

    For localized c(x) (e.g. Gaussian B₀), the convolution kernel ĉ(q)
    decays exponentially, making the matrix effectively banded.  The
    downstream ``_evolve_full_matrix`` exploits this by thresholding small
    entries and converting to sparse CSC format for faster expm_multiply.

    Reference: Burns et al. (2020), Phys. Rev. Research 2:023068.
    """
    n_slots = layout.num_slots
    n_modes = int(np.prod(rfft_shape))
    n_total = n_slots * n_modes

    def _resolve_target_slot(field_ref: str) -> int | None:
        """Resolve a field/velocity reference to a slot index."""
        if field_ref in layout.field_slot_map:
            return layout.field_slot_map[field_ref]
        if field_ref.startswith("v_") and field_ref[2:] in layout.velocity_slot_map:
            return layout.velocity_slot_map[field_ref[2:]]
        return None

    # Evaluate Fourier multipliers on the k-grid (for constant terms)
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op not in multiplier_cache:
                mult_fn = _EXACT_MULTIPLIERS[op]
                mult_val = mult_fn(k_grid)
                mult_full = np.broadcast_to(mult_val, rfft_shape)
                multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    A = np.zeros((n_total, n_total), dtype=np.complex128)

    for _eq_idx, eq in enumerate(spec.equations):
        field_name = eq.field_name
        is_second_order = eq.time_derivative_order >= 2

        if is_second_order:
            field_slot = layout.field_slot_map[field_name]
            vel_slot = layout.velocity_slot_map[field_name]

            # dq/dt = v → diagonal identity coupling between field and velocity
            for m in range(n_modes):
                row = field_slot * n_modes + m
                col = vel_slot * n_modes + m
                A[row, col] = 1.0

            # dv/dt = Σ coeff(x) * operator(target_field)
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = _resolve_target_slot(term.field)
                if target_slot is None:
                    continue
                mult = multiplier_cache[term.operator]

                if not term.position_dependent:
                    # Constant coefficient: diagonal in mode space
                    coeff = _resolve_constant_coeff(
                        term,
                        coeff_eval,
                        eq_idx=_eq_idx,
                        term_idx=_term_idx,
                    )
                    for m in range(n_modes):
                        row = vel_slot * n_modes + m
                        col = target_slot * n_modes + m
                        A[row, col] += coeff * mult[m]
                else:
                    # Position-dependent: convolution coupling
                    _add_convolution_coupling(
                        A,
                        vel_slot,
                        target_slot,
                        term,
                        coeff_eval,
                        mult,
                        grid,
                        rfft_shape,
                        n_modes,
                        eq_idx=_eq_idx,
                        term_idx=_term_idx,
                    )
        else:
            # First-order
            this_slot = layout.field_slot_map[field_name]
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = _resolve_target_slot(term.field)
                if target_slot is None:
                    continue
                mult = multiplier_cache[term.operator]

                if not term.position_dependent:
                    coeff = _resolve_constant_coeff(
                        term,
                        coeff_eval,
                        eq_idx=_eq_idx,
                        term_idx=_term_idx,
                    )
                    for m in range(n_modes):
                        row = this_slot * n_modes + m
                        col = target_slot * n_modes + m
                        A[row, col] += coeff * mult[m]
                else:
                    _add_convolution_coupling(
                        A,
                        this_slot,
                        target_slot,
                        term,
                        coeff_eval,
                        mult,
                        grid,
                        rfft_shape,
                        n_modes,
                        eq_idx=_eq_idx,
                        term_idx=_term_idx,
                    )

    return A


def _add_convolution_coupling(
    A: NDArray[np.complex128],
    row_slot: int,
    col_slot: int,
    term: OperatorTerm,
    coeff_eval: CoefficientEvaluator,
    operator_mult: NDArray[np.complex128],
    grid: GridInfo,
    rfft_shape: tuple[int, ...],
    n_modes: int,
    *,
    eq_idx: int = -1,
    term_idx: int = -1,
) -> None:
    """Add convolution coupling from a position-dependent coefficient.

    The product c(x)·op(u(x)) in k-space becomes a convolution:
    FFT[c·op(u)]_k = Σ_k' ĉ(k-k') · mult(k') · û(k')

    This creates off-diagonal entries in the evolution matrix coupling
    different k-modes.
    """
    # Get the coefficient array on the spatial grid
    coeff_array = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
    if isinstance(coeff_array, (int, float)):
        coeff_array = np.full(grid.shape, float(coeff_array))

    # For each pair of output mode m and input mode m',
    # the coupling is (1/N) * ĉ(m-m') * mult(m')
    # This is a Toeplitz-like structure in 1D.
    #
    # Build via outer product approach for efficiency:
    # We compute the full convolution matrix using FFT properties.
    #
    # For rfftn: the convolution of real functions in rfft space requires
    # care with the Hermitian symmetry. We use the identity:
    # FFT[c·u]_k = (1/N) Σ_{k'} ĉ_{k-k'} · û_{k'}
    #
    # Build the convolution matrix C where C[k, k'] = (1/N) * ĉ_{k-k'}
    # using probe vectors (unit impulse per mode).
    for m_prime in range(n_modes):
        # Probe: unit impulse at mode m_prime
        probe_hat = np.zeros(n_modes, dtype=np.complex128)
        probe_hat[m_prime] = 1.0

        # Reconstruct to physical space, multiply by coefficient, FFT back
        probe_physical = np.fft.irfftn(
            probe_hat.reshape(rfft_shape),
            s=grid.shape,
            axes=list(range(len(grid.shape))),
        )
        product = coeff_array * probe_physical
        result_hat = np.fft.rfftn(product).ravel()

        # result_hat[m] = Σ_{k'} (1/N) ĉ_{m-k'} δ_{k',m'} = (1/N) ĉ_{m-m'}
        # multiplied by operator multiplier at m'
        row_start = row_slot * n_modes
        col = col_slot * n_modes + m_prime
        A[row_start : row_start + n_modes, col] += result_hat * operator_mult[m_prime]


# ---------------------------------------------------------------------------
# Block decomposition
# ---------------------------------------------------------------------------


def find_independent_blocks(
    A: NDArray[np.complex128],
    threshold: float = 1e-14,
) -> list[list[int]]:
    """Find independent (decoupled) blocks in an evolution matrix.

    Analyzes the sparsity pattern of A: slots i and j are coupled if
    |A[i,j]| > threshold or |A[j,i]| > threshold.  Returns a list of
    slot-index groups (connected components).

    This prevents degenerate-eigenvalue mixing when ``np.linalg.eig``
    processes a block-diagonal matrix with repeated eigenvalues across
    independent blocks — a common situation in symmetric multi-field
    theories (e.g. Gertsenshtein h₅↔a₁ + h₇↔a₂).
    """
    n = A.shape[0]
    # Union-find (path compression + union by rank)
    parent = list(range(n))
    rank = [0] * n

    def _find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x

    def _union(x: int, y: int) -> None:
        rx, ry = _find(x), _find(y)
        if rx == ry:
            return
        if rank[rx] < rank[ry]:
            rx, ry = ry, rx
        parent[ry] = rx
        if rank[rx] == rank[ry]:
            rank[rx] += 1

    # Build coupling graph from matrix entries
    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i, j]) > threshold or abs(A[j, i]) > threshold:
                _union(i, j)

    # Group by root
    groups: dict[int, list[int]] = {}
    for i in range(n):
        root = _find(i)
        groups.setdefault(root, []).append(i)
    return list(groups.values())


# ---------------------------------------------------------------------------
# Eigendecomposition and time evolution
# ---------------------------------------------------------------------------


def _has_position_dependent_terms(spec: EquationSystem) -> bool:
    """Check if any RHS term has a position-dependent coefficient."""
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.position_dependent:
                return True
    return False


def _suppress_tachyonic_noise(
    eig_vals: NDArray[np.complex128],
    y0_eigen: NDArray[np.complex128],
    *,
    growth_threshold: float = 1e-8,
    coupling_threshold: float = 1e-12,
) -> tuple[NDArray[np.complex128], int]:
    """Suppress tachyonic eigenvectors with zero physical coupling.

    In theories like PGT with R̃, the mass matrix can have tachyonic eigenvalues
    (positive real part) corresponding to tensor/axial torsion sectors.  These
    modes may have **zero physical coupling** — the source terms project entirely
    onto the stable (trace) eigenvector.  However, machine-precision noise
    (~10⁻¹⁶) in these uncoupled modes grows as exp(Re(λ)·t), overwhelming the
    physical solution.

    This function identifies modes that are:

    1. Tachyonic: ``Re(λ) > growth_threshold`` (exponentially growing)
    2. Uncoupled: ``|IC projection| < coupling_threshold * max(|IC|)``
       (zero physical coupling within numerical precision)

    For such modes, the eigenvalue is set to zero (freezing the mode at its
    initial noise level) instead of allowing exponential amplification.

    This is a **numerical noise fix**, not a physics change: if the IC genuinely
    projects onto a tachyonic eigenvector, the mode is kept.  The threshold
    scales with IC magnitude to be robust across different amplitudes.

    Follows the existing gauge-mode suppression pattern (infinite eigenvalues
    set to zero for gauge degrees of freedom).

    Parameters
    ----------
    eig_vals : ndarray, shape (n_modes, block_size)
        Eigenvalues per k-mode and eigenvector.
    y0_eigen : ndarray, shape (n_modes, block_size)
        IC amplitude projected onto each eigenvector.
    growth_threshold : float
        Minimum positive real part to classify as tachyonic.
    coupling_threshold : float
        Relative threshold for "zero coupling" (fraction of max IC).

    Returns
    -------
    eig_vals : ndarray
        Modified eigenvalues (tachyonic uncoupled modes set to 0).
    n_suppressed : int
        Number of suppressed mode-components.
    """
    max_ic = np.max(np.abs(y0_eigen))
    if max_ic == 0:
        return eig_vals, 0

    noise_floor = coupling_threshold * max_ic

    tachyonic = np.real(eig_vals) > growth_threshold
    uncoupled = np.abs(y0_eigen) < noise_floor
    suppress_mask = tachyonic & uncoupled

    n_suppressed = int(np.sum(suppress_mask))
    if n_suppressed > 0:
        eig_vals = eig_vals.copy()
        eig_vals[suppress_mask] = 0.0

    return eig_vals, n_suppressed


def _warn_eigenvalue_growth(
    eigenvalues: NDArray[np.complex128],
    dt_total: float,
    context: str = "",
) -> None:
    """Warn if eigenvalues have positive real parts that could cause overflow."""
    import warnings  # noqa: PLC0415

    max_real = float(np.max(np.real(eigenvalues)))
    if max_real > 1e-10:
        max_growth = max_real * dt_total
        if max_growth > 30:  # exp(30) ≈ 1e13
            ctx = f" ({context})" if context else ""
            warnings.warn(
                f"Modal solver{ctx}: eigenvalues with positive real parts "
                f"(max Re(λ)={max_real:.3e}). Growth factor exp({max_growth:.1f}) "
                f"over Δt={dt_total:.1f} may cause overflow. "
                f"Consider --scheme cvode for numerical stability.",
                stacklevel=3,
            )


def _evolve_per_mode(
    A_modes: NDArray[np.complex128],
    y0_hat: NDArray[np.complex128],
    t_eval: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
    snapshot_callback: Callable[[float, NDArray[np.float64]], None] | None,
    progress: SimulationProgress | None,
    *,
    return_fourier: bool = False,
    return_derivative_fourier: bool = False,
    B_modes: NDArray[np.complex128] | None = None,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.complex128] | None,
    NDArray[np.complex128] | None,
]:
    """Evolve system with per-mode independent matrices (constant coefficients).

    A_modes has shape (n_modes, n_slots, n_slots).
    y0_hat has shape (n_slots, n_modes).

    Uses block-aware eigendecomposition: independent field blocks are detected
    and eigendecomposed separately to prevent degenerate-eigenvalue mixing.
    Blocks with all-zero initial conditions are skipped entirely.

    Ref: Golub & Van Loan (1996), Matrix Computations, §4.8.

    Raises
    ------
    SimulationDivergedError
        If field amplitudes grow beyond 100x the initial maximum, which
        indicates the simulation has left the perturbative regime where
        the linearized equations are valid.
    """
    from tidal.solver._exceptions import SimulationDivergedError  # noqa: PLC0415

    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_snapshots = len(t_eval)
    t0 = t_eval[0]
    dt_total = float(t_eval[-1] - t0)

    # Detect independent blocks from the first mode's matrix.
    # Block structure is k-independent for constant coefficients, so we only
    # need to analyze one representative mode (use max across a few modes for
    # robustness against accidental zeros at specific k).
    n_check = min(3, A_modes.shape[0])
    combined = np.max(np.abs(A_modes[:n_check]), axis=0)
    blocks = find_independent_blocks(combined)

    # Pre-compute eigendecomposition for each active block
    block_data: list[
        tuple[
            list[int],  # slot indices
            NDArray[np.complex128],  # eigenvalues (n_modes, block_size)
            NDArray[np.complex128],  # V (n_modes, block_size, block_size)
            NDArray[np.complex128],  # y0_eigen (n_modes, block_size)
        ]
    ] = []

    for block_slots in blocks:
        # Extract IC for this block
        y0_block = y0_hat[block_slots, :]  # (block_size, n_modes)

        # Skip blocks with all-zero IC — output stays at zero
        if np.max(np.abs(y0_block)) < 1e-15:
            continue

        # Extract block sub-matrices: (n_modes, block_size, block_size)
        idx = np.array(block_slots)
        A_block = A_modes[:, idx[:, None], idx[None, :]]

        if B_modes is not None:
            # Generalized eigenvalue problem: B · d' = A · d
            # Uses QZ decomposition via scipy.linalg.eig(A, B).
            # Infinite eigenvalues (gauge DOF) are zeroed — they don't evolve.
            # Essential for rank-deficient M: pre-solving via np.linalg.solve
            # produces ill-conditioned eigenvectors (cond(V) > 1e9), whereas
            # QZ handles rank-deficiency via the Schur form without requiring
            # V to be well-conditioned.
            # Ref: Golub & Van Loan (2013), Matrix Computations §7.7.6
            import scipy.linalg as sla  # type: ignore[import-untyped]  # noqa: PLC0415

            B_block = B_modes[:, idx[:, None], idx[None, :]]
            bs = len(block_slots)
            n_block_modes = A_block.shape[0]
            eig_vals = np.zeros((n_block_modes, bs), dtype=np.complex128)
            v_mat = np.zeros((n_block_modes, bs, bs), dtype=np.complex128)
            n_gauge_total = 0
            for m in range(A_block.shape[0]):
                A_m = A_block[m]
                B_m = B_block[m]
                # Null-space projection before QZ: rank-deficient B (e.g. CDT
                # rank-3 PGT torsion where only the trace subspace has kinetic
                # terms) causes scipy.linalg.eig to return FINITE spurious
                # eigenvalues not caught by the |λ|>1e12 filter.  Project A
                # and B onto range(B) first; null directions get eigenvalue 0
                # (frozen at IC).  See issue #257.
                _, s_b, Vt_b = np.linalg.svd(B_m)
                null_thresh_b = s_b[0] * 1e-10 if s_b[0] > 0 else 1e-14
                rank_b = int(np.sum(s_b > null_thresh_b))
                null_dim_b = bs - rank_b
                if null_dim_b > 0:
                    Vphys = Vt_b[:rank_b].T  # (bs, rank_b) — physical subspace
                    Vnull = Vt_b[rank_b:].T  # (bs, null_dim_b) — null(B)
                    eig_r = sla.eig(  # pyright: ignore[reportUnknownVariableType]
                        Vphys.T @ A_m @ Vphys,
                        Vphys.T @ B_m @ Vphys,
                        right=True,
                    )
                    ev_red = np.asarray(eig_r[0], dtype=np.complex128)  # pyright: ignore[reportUnknownArgumentType]
                    vr_red = np.asarray(eig_r[1], dtype=np.complex128)  # pyright: ignore[reportUnknownArgumentType]
                    # Lift eigenvectors to full space; null modes frozen at IC
                    ev_m: NDArray[np.complex128] = np.concatenate(
                        [ev_red, np.zeros(null_dim_b, dtype=np.complex128)]
                    )
                    # V_full = [Vphys @ vr_red | Vnull] is always invertible:
                    # physical cols live in range(Vphys), null cols in null(B),
                    # and SVD guarantees these two subspaces are orthogonal.
                    vr_m: NDArray[np.complex128] = np.hstack([Vphys @ vr_red, Vnull])
                    n_gauge_total += null_dim_b
                else:
                    eig_r2 = sla.eig(A_m, B_m, right=True)  # pyright: ignore[reportUnknownVariableType]
                    ev_m = np.asarray(eig_r2[0], dtype=np.complex128)  # pyright: ignore[reportUnknownArgumentType]
                    vr_m = np.asarray(eig_r2[1], dtype=np.complex128)  # pyright: ignore[reportUnknownArgumentType]
                # Filter any remaining infinite/very-large eigenvalues
                # (gauge modes not caught by null-space projection).
                gauge = ~np.isfinite(ev_m) | (np.abs(ev_m) > 1e12)
                ev_m[gauge] = 0.0  # gauge modes frozen at IC
                n_gauge_total += int(np.sum(gauge))
                eig_vals[m] = ev_m
                v_mat[m] = vr_m
            if n_gauge_total > 0:
                import logging as _log  # noqa: PLC0415

                _log.getLogger(__name__).info(
                    "Generalized eigenvalue: %d gauge modes zeroed across %d modes",
                    n_gauge_total,
                    A_block.shape[0],
                )
            v_inv = np.linalg.inv(v_mat)
        else:
            # Regular eigendecomposition for non-generalized systems.
            eig_vals, v_mat = np.linalg.eig(A_block)
            v_inv = np.linalg.inv(v_mat)

        # Warn about potential overflow
        _warn_eigenvalue_growth(eig_vals, dt_total, context="per-mode")

        # Transform IC to eigenbasis
        y0_eigen = np.einsum("mij,mj->mi", v_inv, y0_block.T)

        # Suppress tachyonic modes with zero physical coupling.
        # These are numerical noise amplifiers, not physics — see #222.
        eig_vals, n_suppressed = _suppress_tachyonic_noise(eig_vals, y0_eigen)
        if n_suppressed > 0:
            import logging as _log_tach  # noqa: PLC0415

            _log_tach.getLogger(__name__).info(
                "Suppressed %d tachyonic modes with zero IC coupling "
                "(numerical noise prevention)",
                n_suppressed,
            )

        block_data.append((block_slots, eig_vals, v_mat, y0_eigen))

    # Evolve at each time point.
    # Pre-multiply V @ diag(y0_eigen) for each block so the inner loop only
    # needs element-wise exp + matrix-vector product, not a full einsum.
    block_evolved: list[
        tuple[
            list[int],  # slot indices
            NDArray[np.complex128],  # V_y0: V * y0_eigen, (n_modes, bs, bs)
            NDArray[np.complex128] | None,  # V_y0_deriv (n_modes, bs, bs) or None
            NDArray[np.complex128],  # eigenvalues (n_modes, bs)
        ]
    ] = []
    for block_slots, eig_vals, v_mat, y0_eigen in block_data:
        # V_y0[m, i, j] = v_mat[m, i, j] * y0_eigen[m, j]
        # so y(t) = V_y0 @ exp(λ*dt) is just a matvec
        V_y0 = v_mat * y0_eigen[:, np.newaxis, :]  # (n_modes, bs, bs)
        # V_y0_deriv[m, i, j] = V_y0[m, i, j] * λ[m, j]
        # so y'(t) = V_y0_deriv @ exp(λ*dt) gives exact time derivative
        V_y0_deriv = (
            V_y0 * eig_vals[:, np.newaxis, :] if return_derivative_fourier else None
        )
        block_evolved.append((block_slots, V_y0, V_y0_deriv, eig_vals))

    snapshots = np.zeros((n_snapshots, n_slots * n_pts))
    times = np.zeros(n_snapshots)
    n_modes = y0_hat.shape[1]

    # Optionally collect Fourier-space snapshots (avoids re-FFT in constraint
    # recovery — the Fourier data is already computed here).
    fourier_snaps: NDArray[np.complex128] | None = None
    if return_fourier:
        fourier_snaps = np.zeros(
            (n_snapshots, n_slots, n_modes),
            dtype=np.complex128,
        )

    # Optionally collect Fourier-space TIME DERIVATIVE snapshots.
    # d'(t) = V · diag(λ · exp(λt)) · y0_eigen — exact, no numerical diff.
    # Used for machine-precision constraint velocity: v_c = recovery · d'.
    deriv_fourier_snaps: NDArray[np.complex128] | None = None
    dy_hat_t: NDArray[np.complex128] | None = None
    if return_derivative_fourier:
        deriv_fourier_snaps = np.zeros(
            (n_snapshots, n_slots, n_modes),
            dtype=np.complex128,
        )
        dy_hat_t = np.zeros((n_slots, n_modes), dtype=np.complex128)

    # Pre-allocate buffer reused each timestep (avoids n_snapshots allocations)
    y_hat_t = np.zeros((n_slots, n_modes), dtype=np.complex128)

    # Divergence guard: perturbation-theory validity check.
    #
    # Physics reasoning: the solver operates on LINEARIZED equations,
    # derived by expanding to first order in a small perturbation
    # parameter epsilon ~ IC amplitude.  Higher-order terms (F^4, R*F^2,
    # torsion * F^2, ...) were dropped.  The linearized evolution is only
    # physical while all fields remain O(epsilon) — once any field grows
    # to O(100 * epsilon), the dropped terms become comparable to the
    # retained ones and the physics is unreliable regardless of numerical
    # stability.
    #
    # We therefore enforce: max(|y(t)|) / max(|y(0)|) <= divergence_threshold.
    # This is a ratio check (not absolute) so it scales with IC amplitude.
    # The threshold 100 puts us at the boundary of the perturbative regime
    # (well before overflow).
    #
    # Discrimination properties:
    #   - Normal evolution: ratio ~ O(1), passes.
    #   - Legitimate amplification A up to ~10^5: target field grows to
    #     sqrt(A * P_GR) * epsilon ~ a few * epsilon, ratio < 100, passes.
    #   - Tachyonic instability: source field itself grows exponentially,
    #     ratio rapidly exceeds 100, rejected.
    #   - Partial blow-up (target blows up, source stays bounded): target
    #     reaches ~100 * epsilon meaning perturbation theory has broken,
    #     ratio exceeds 100, rejected.
    #
    # The threshold is chosen for physics (perturbation-theory validity),
    # not numerics — it errs on the side of being strict about "is this
    # still a linearized simulation" rather than "will this overflow".
    initial_max_amp: float = 0.0
    divergence_threshold: float = 100.0

    # NOTE: The eigenvalue pre-check guard (commit 2b94172) was removed.
    # It detected modes with Re(λ) > 0 and physical IC projection, but for
    # PGT models with constrained torsion fields, the constraint-eliminated
    # first-order system generically has eigenmodes with large positive real
    # parts that represent wave propagation characteristics (not physical
    # instabilities).  The _suppress_tachyonic_noise function handles modes
    # with zero physical coupling, and the per-snapshot divergence guard
    # (below, in the time loop) catches genuinely diverging simulations.
    # The pre-check was overly aggressive: it blocked the GR baseline of
    # the nonminimal PGT model, which produces correct oscillatory P(t).

    for ti, t in enumerate(t_eval):
        dt = t - t0
        y_hat_t[:] = 0.0
        if dy_hat_t is not None:
            dy_hat_t[:] = 0.0

        for block_slots, V_y0, V_y0_deriv, eig_vals in block_evolved:
            # exp_lambda shape: (n_modes, block_size)
            exp_lambda = np.exp(eig_vals * dt)
            # y_evolved[m, i] = Σ_j V_y0[m, i, j] * exp(λ_j * dt)
            y_evolved = np.einsum("mij,mj->mi", V_y0, exp_lambda)
            y_hat_t[block_slots, :] = y_evolved.T

            # Exact time derivative: dy[m,i] = Σ_j V_y0_deriv[m,i,j] * exp(λ_j*dt)
            if V_y0_deriv is not None and dy_hat_t is not None:
                dy_evolved = np.einsum("mij,mj->mi", V_y0_deriv, exp_lambda)
                dy_hat_t[block_slots, :] = dy_evolved.T

        if fourier_snaps is not None:
            fourier_snaps[ti] = y_hat_t
        if deriv_fourier_snaps is not None and dy_hat_t is not None:
            deriv_fourier_snaps[ti] = dy_hat_t

        y_physical = _ifft_slots(y_hat_t, layout, grid)
        snapshots[ti] = y_physical
        times[ti] = t

        # Divergence guard: enforce perturbation-theory validity via
        # global amplitude ratio.  See declaration comment for physics.
        max_amp = float(np.max(np.abs(y_physical)))
        if ti == 0:
            initial_max_amp = max(max_amp, 1e-15)
        elif (
            not np.isfinite(max_amp) or max_amp / initial_max_amp > divergence_threshold
        ):
            msg = (
                f"Simulation diverged at t={t:.4g}: amplitude ratio "
                f"{max_amp / initial_max_amp:.2e} exceeds threshold "
                f"{divergence_threshold:.0e}. Fields have left the "
                f"perturbative regime (linearized approximation invalid)."
            )
            # Fill remaining snapshots so partial results are usable
            for tj in range(ti + 1, n_snapshots):
                snapshots[tj] = y_physical
                times[tj] = t_eval[tj]
            raise SimulationDivergedError(msg)

        if snapshot_callback is not None:
            snapshot_callback(t, y_physical)

        if progress is not None:
            progress.update(t)

    return times, snapshots, fourier_snaps, deriv_fourier_snaps


def _evolve_full_matrix(
    A_full: NDArray[np.complex128],
    y0_hat: NDArray[np.complex128],
    t_eval: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
    snapshot_callback: Callable[[float, NDArray[np.float64]], None] | None,
    progress: SimulationProgress | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evolve system with full coupled matrix (position-dependent coefficients).

    A_full has shape (n_total, n_total) where n_total = n_slots x n_modes.
    y0_hat has shape (n_slots, n_modes).

    Uses ``scipy.sparse.linalg.expm_multiply`` to compute exp(A·t)·y₀ at each
    output time without eigendecomposition.  This is backward-stable for
    non-normal matrices (position-dependent gradient coupling creates non-normal
    convolution matrices whose eigenvalues have large real parts, but the true
    dynamics are bounded).  The algorithm uses scaling + truncated Taylor series
    in matrix-vector products, avoiding individual exp(λ·t) overflow.

    **Why not eigendecomposition?**  The original full-matrix eigendecomposition
    gave incorrect physics for localized Gertsenshtein (P=0.477 vs correct
    P=0.3437) because non-normal convolution matrices have eigenvalues with
    significant positive real parts despite conservative physics — individual
    exp(λ·t) overflow while exp(A·t)·y₀ is bounded (pseudospectral phenomenon;
    Trefethen & Embree 2005, Ch. 14).

    **Sparse optimization:** For localized coefficients (e.g. Gaussian B₀), the
    convolution kernel ĉ(q) decays exponentially, making the matrix effectively
    banded.  Entries below a relative threshold (1e-14 x max|A|) are zeroed, and
    if density < 30% the matrix is converted to sparse CSC format.  This
    accelerates expm_multiply's internal matrix-vector products.

    Ref: Al-Mohy & Higham (2011), "Computing the Action of the Matrix
    Exponential", SIAM J. Sci. Comput. 33(2):488-511.
    """
    import scipy.sparse  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
    from scipy.sparse.linalg import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        expm_multiply,  # pyright: ignore[reportUnknownVariableType]
    )

    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_modes = y0_hat.shape[1]
    n_snapshots = len(t_eval)

    # Flatten y0_hat to (n_total,) — slot-major order
    y0_flat = y0_hat.ravel()

    # --- Sparse matrix optimization ---
    # Position-dependent convolution matrices are effectively banded: for
    # Gaussian B₀(x) the kernel ĉ(q) decays exponentially, so most off-
    # diagonal entries are negligibly small.  Thresholding and converting to
    # sparse CSC format accelerates expm_multiply's internal matrix-vector
    # products (the dominant cost) without affecting accuracy.
    abs_max = float(np.max(np.abs(A_full)))
    if abs_max > 0:
        threshold = abs_max * 1e-14
        A_work = A_full.copy()
        A_work[np.abs(A_work) < threshold] = 0.0
        density = np.count_nonzero(A_work) / A_work.size
        if density < 0.3:
            A_op: NDArray[np.complex128] | scipy.sparse.csc_array = (
                scipy.sparse.csc_array(A_work)
            )
        else:
            A_op = A_work
    else:
        A_op = A_full
        density = 1.0

    snapshots = np.zeros((n_snapshots, n_slots * n_pts))
    times = np.zeros(n_snapshots)

    # Compute exp(A·t)·y₀ at all requested times using Krylov/Taylor method.
    # expm_multiply handles the scaling internally — no manual dt stepping.
    t0 = float(t_eval[0])
    t_end = float(t_eval[-1])

    if n_snapshots > 1 and t_end > t0:
        # Use expm_multiply's built-in multi-point evaluation
        y_all: NDArray[np.complex128] = np.asarray(
            expm_multiply(
                A_op,
                y0_flat,
                start=t0,
                stop=t_end,
                num=n_snapshots,
            ),
            dtype=np.complex128,
        )
        # y_all has shape (n_snapshots, n_total)
        for ti in range(n_snapshots):
            t = float(t_eval[ti])
            y_hat_t = y_all[ti].reshape(n_slots, n_modes)
            y_physical = _ifft_slots(y_hat_t, layout, grid)
            snapshots[ti] = y_physical
            times[ti] = t

            if snapshot_callback is not None:
                snapshot_callback(t, y_physical)
            if progress is not None:
                progress.update(t)
    else:
        # Single time point or t0 == t_end
        for ti, t in enumerate(t_eval):
            if t == t0:
                y_evolved = y0_flat.copy()
            else:
                y_evolved = np.asarray(
                    expm_multiply(A_op, y0_flat, start=t0, stop=float(t), num=2)[-1],
                    dtype=np.complex128,
                )
            y_hat_t = y_evolved.reshape(n_slots, n_modes)
            y_physical = _ifft_slots(y_hat_t, layout, grid)
            snapshots[ti] = y_physical
            times[ti] = t

            if snapshot_callback is not None:
                snapshot_callback(t, y_physical)
            if progress is not None:
                progress.update(t)

    return times, snapshots


# ---------------------------------------------------------------------------
# Public solver entry point
# ---------------------------------------------------------------------------


def solve_modal(
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    *,
    bc: BCSpec | None = None,
    parameters: dict[str, float] | None = None,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    num_snapshots: int = 101,
    snapshot_callback: Callable[[float, np.ndarray], None] | None = None,
    progress: SimulationProgress | None = None,
) -> SolverResult:
    """Solve a TIDAL equation system using Fourier modal decomposition.

    Transforms the spatial grid to Fourier space, builds the per-mode or
    full evolution matrix, eigendecomposes, and evaluates the exact solution
    at each output time.

    Parameters
    ----------
    spec : EquationSystem
        Parsed equation specification (from JSON).
    grid : GridInfo
        Spatial grid (must be all-periodic).
    y0 : np.ndarray
        Initial state vector (flat).
    t_span : tuple[float, float]
        (t_start, t_end).
    bc : str or tuple, optional
        Boundary conditions (must be all-periodic).
    parameters : dict[str, float], optional
        Runtime parameter overrides for symbolic coefficients.
    rtol, atol : float
        Tolerances (unused for eigendecomposition; reserved for solve_ivp
        fallback with time-dependent coefficients).
    num_snapshots : int
        Number of output time points.
    snapshot_callback : callable, optional
        Called as ``callback(t, y)`` at each output time.
    progress : SimulationProgress, optional
        Progress tracker for tqdm display.

    Returns
    -------
    SolverResult
        Dict with keys: ``t``, ``y``, ``success``, ``message``.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, parameters or {})

    # Detect constraint fields
    has_constraints = any(eq.time_derivative_order == 0 for eq in spec.equations)
    if not has_constraints:
        warn_frozen_constraints(layout, "modal")

    # Build time evaluation points
    t_eval = np.linspace(t_span[0], t_span[1], num_snapshots)

    # Build k-grid
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)

    # Compute rfft output shape
    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    # FFT initial conditions
    y0_hat = _fft_slots(y0, layout, grid)

    # Zero the Nyquist mode(s) in the IC.  The rfft Nyquist bin (last mode
    # in each dimension) must be real for real-valued fields.  The modal
    # evolution matrix has complex entries (from gradient coupling ik),
    # which creates imaginary components at the Nyquist bin.  irfft
    # silently drops these, causing energy non-conservation proportional
    # to the Nyquist power.  Zeroing the Nyquist IC prevents this entirely.
    # This is standard practice in spectral methods — the Nyquist mode
    # aliases with its conjugate and cannot represent physical content.
    # Ref: Boyd (2001), Chebyshev & Fourier Spectral Methods, §11.5.
    for _dim_idx, n in enumerate(grid.shape):
        if n % 2 == 0:  # Nyquist mode exists only for even N
            nyq_mode = n // 2  # last rfft bin
            if len(grid.shape) == 1:
                y0_hat[:, nyq_mode] = 0.0
            else:
                # Multi-D: zero along the last-axis Nyquist slice
                rfft_last = grid.shape[-1] // 2
                y0_hat[:, ..., rfft_last] = 0.0

    has_pos_dep = _has_position_dependent_terms(spec)
    has_time_ops = _has_time_derivative_operators(spec)

    # Single evolution-matrix builder handles:
    #   - algebraic constraint fields (time_derivative_order=0) via Schur
    #   - rank-deficient mass matrices (cross d2_t couplings) via mass-
    #     eigendecomposition Schur
    #   - near-singular velocity-coupling (#166/#175) via np.linalg.solve
    #     with lstsq fallback (LU+partial-pivot stable up to cond~1e15,
    #     unlike the generalized-eig QZ path that used to produce spurious
    #     large finite eigenvalues)
    #
    # The builder returns (A_rhs, B_lhs).  If vel_coupling is present,
    # B_lhs is non-None and the downstream eigendecomposition uses
    # scipy.linalg.eig(A, B) (QZ decomposition), which handles
    # rank-deficient M and near-singular (I - vel_coupling) correctly
    # via the Schur form — pre-solving via np.linalg.solve produces
    # ill-conditioned V (cond > 1e9) for rank-deficient M and was rejected.
    needs_reduction = (has_constraints or has_time_ops) and not has_pos_dep
    constraint_vel_arrays: dict[str, NDArray[np.float64]] = {}  # populated below

    if needs_reduction:
        (
            A_reduced,
            B_lhs_modes,
            recovery_matrix,
            _v_recovery_gen,  # unused — constraint vel from eigendata
            c_names,
            orig_to_reduced,
        ) = _build_evolution_matrices(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
        )

        n_dyn = A_reduced.shape[1]
        n_modes = y0_hat.shape[1]
        n_pts = layout.num_points

        # Extract dynamical IC in reduced ordering
        y0_hat_dyn = np.zeros((n_dyn, n_modes), dtype=np.complex128)
        for orig_si, red_pos in orig_to_reduced.items():
            y0_hat_dyn[red_pos] = y0_hat[orig_si]

        # Build a reduced StateLayout for eigendecomposition
        sorted_orig = sorted(orig_to_reduced.keys())
        red_slots = tuple(layout.slots[si] for si in sorted_orig)
        red_field_map: dict[str, int] = {}
        red_vel_map: dict[str, int] = {}
        for new_i, si in enumerate(sorted_orig):
            s = layout.slots[si]
            if s.kind == "field":
                red_field_map[s.field_name] = new_i
            elif s.kind == "velocity":
                red_vel_map[s.field_name] = new_i
        dyn_layout = StateLayout(
            slots=red_slots,
            num_points=n_pts,
            field_slot_map=red_field_map,
            velocity_slot_map=red_vel_map,
            dynamical_fields=layout.dynamical_fields,
        )

        # Evolve dynamical fields (return Fourier + derivative data)
        times, dyn_snapshots, dyn_fourier, dyn_deriv_fourier = _evolve_per_mode(
            A_reduced,
            y0_hat_dyn,
            t_eval,
            dyn_layout,
            grid,
            None,
            progress,  # callback handled below with full state
            return_fourier=True,
            return_derivative_fourier=True,  # for exact constraint velocities
            B_modes=B_lhs_modes,  # generalized eig(A, B) when vel_coupling present
        )

        # Reconstruct full state (including constraints) at each snapshot
        n_full = layout.num_slots * n_pts
        snapshots = np.zeros((len(t_eval), n_full))
        assert dyn_fourier is not None  # guaranteed by return_fourier=True

        # Populate constraint velocity arrays: exact ∂_t(c) from eigendata.
        # "Constraint" is a solver concept (algebraic evolution), not a physics
        # statement — these fields have physically meaningful velocities.
        # v_c(t) = recovery · d'(t), where d'(t) is computed from eigendata
        # inside _evolve_per_mode (V·diag(λ·exp(λt))·y0_eigen — exact).
        for c_name in c_names:
            constraint_vel_arrays[c_name] = np.zeros((len(t_eval), *grid.shape))

        for ti in range(len(t_eval)):
            dyn_phys = dyn_snapshots[ti]
            # Use Fourier data directly (already computed in _evolve_per_mode)
            y_hat_dyn_t = dyn_fourier[ti]  # (n_dyn, n_modes)

            # Recover constraint fields: c_hat = recovery @ d_hat
            c_hat = np.einsum("mcj,jm->cm", recovery_matrix, y_hat_dyn_t)

            # Recover constraint velocities: v_c_hat = recovery @ d'_hat
            # d'_hat comes from eigendata — exact, no numerical differentiation.
            assert dyn_deriv_fourier is not None
            dy_hat_dyn_t = dyn_deriv_fourier[ti]  # (n_dyn, n_modes)
            v_c_hat = np.einsum("mcj,jm->cm", recovery_matrix, dy_hat_dyn_t)

            # Assemble full physical state
            full_state = np.zeros(n_full)
            for orig_si, red_pos in orig_to_reduced.items():
                full_state[orig_si * n_pts : (orig_si + 1) * n_pts] = dyn_phys[
                    red_pos * n_pts : (red_pos + 1) * n_pts
                ]
            for ci, c_name in enumerate(c_names):
                c_slot = layout.field_slot_map[c_name]
                c_phys = np.fft.irfftn(
                    c_hat[ci].reshape(rfft_shape),
                    s=grid.shape,
                    axes=list(range(len(grid.shape))),
                ).ravel()
                full_state[c_slot * n_pts : (c_slot + 1) * n_pts] = np.real(
                    c_phys,
                )
                # Store exact constraint velocity (from eigendata d')
                v_c_phys = np.fft.irfftn(
                    v_c_hat[ci].reshape(rfft_shape),
                    s=grid.shape,
                    axes=list(range(len(grid.shape))),
                )
                constraint_vel_arrays[c_name][ti] = np.real(v_c_phys)

            snapshots[ti] = full_state
            if snapshot_callback is not None:
                snapshot_callback(t_eval[ti], full_state)

        n_c = len(c_names)
        method_desc = (
            f"per-mode eigendecomposition with unified Schur elimination "
            f"({n_c} constraint fields, {n_dyn} dynamical slots, mass-matrix)"
        )

    elif not has_pos_dep:
        # All-constant coefficients: per-mode independent evolution
        A_modes = _build_per_mode_matrices(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
        )
        times, snapshots, _, _ = _evolve_per_mode(
            A_modes,
            y0_hat,
            t_eval,
            layout,
            grid,
            snapshot_callback,
            progress,
        )
        method_desc = "per-mode eigendecomposition (constant coefficients)"
    else:
        # Position-dependent coefficients: full convolution matrix
        A_full = _build_convolution_matrix(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
        )
        times, snapshots = _evolve_full_matrix(
            A_full,
            y0_hat,
            t_eval,
            layout,
            grid,
            snapshot_callback,
            progress,
        )
        n_total = A_full.shape[0]
        method_desc = f"expm_multiply ({n_total}x{n_total}, position-dependent)"

    if progress is not None:
        progress.finish()

    result: SolverResult = {
        "t": times,
        "y": snapshots,
        "success": True,
        "message": f"Modal solver completed ({method_desc})",
    }
    # Attach constraint velocity arrays (exact ∂_t for constraint fields).
    # These are computed from v_recovery = recovery @ A_reduced inside the
    # Schur elimination path. For generalized eigenvalue or non-constraint
    # systems, constraint_vel_arrays is empty.
    if constraint_vel_arrays:
        result["constraint_velocities"] = constraint_vel_arrays  # type: ignore[typeddict-unknown-key]
    return result
