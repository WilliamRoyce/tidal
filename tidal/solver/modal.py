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

For second-order wave equations the Hamiltonian structure A = [[0, I], [L, 0]]
is exploited: only L is eigendecomposed (halves the eigenproblem), and the
solution uses cos/sin matrix functions.

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
# ruff: noqa: C901, RUF001, RUF002 — complexity and Unicode math symbols.
# ruff: noqa: ERA001, ARG001 — commented-out code serves as documentation;
#   unused args (bc, grid) kept for interface consistency with other solvers.

from __future__ import annotations

from collections.abc import Callable
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._setup import warn_frozen_constraints
from tidal.solver.operators import _get_wavenumbers, is_periodic_bc
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from tidal.solver._types import SolverResult
    from tidal.solver.coefficients import CoefficientEvaluator
    from tidal.solver.grid import GridInfo
    from tidal.solver.operators import BCSpec
    from tidal.solver.progress import SimulationProgress
    from tidal.symbolic.json_loader import EquationSystem, OperatorTerm

# ---------------------------------------------------------------------------
# Exact Fourier multipliers (angular wavenumber convention: k = 2π·rfftfreq)
# ---------------------------------------------------------------------------
# These use the EXACT wavenumber k, consistent with operators.py spectral
# mode (gradient → ik, laplacian → -k²).  NOT the modified-wavenumber
# convention from constraint_solve.py which matches FD stencils.

_ExactMultFn = Callable[[list[NDArray[np.float64]]], NDArray[np.complex128]]

_EXACT_MULTIPLIERS: dict[str, _ExactMultFn] = {
    "identity": lambda k_axes: np.ones_like(k_axes[0]),
    "laplacian": lambda k_axes: -sum(ki**2 for ki in k_axes),  # type: ignore[return-value]
    "laplacian_x": lambda k_axes: -(k_axes[0] ** 2),
    "laplacian_y": lambda k_axes: -(k_axes[1] ** 2),
    "laplacian_z": lambda k_axes: -(k_axes[2] ** 2),
    "gradient_x": lambda k_axes: 1j * k_axes[0],
    "gradient_y": lambda k_axes: 1j * k_axes[1],
    "gradient_z": lambda k_axes: 1j * k_axes[2],
    "cross_derivative_xy": lambda k_axes: -(k_axes[0] * k_axes[1]),
    "cross_derivative_xz": lambda k_axes: -(k_axes[0] * k_axes[2]),
    "cross_derivative_yz": lambda k_axes: -(k_axes[1] * k_axes[2]),
    "biharmonic": lambda k_axes: sum(ki**2 for ki in k_axes) ** 2,  # type: ignore[return-value]
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

    # 4. All operators supported
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.operator not in _EXACT_MULTIPLIERS:
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
    # For nD: use rfftn which produces shape[:-1] + (shape[-1]//2+1,)
    # We flatten the modal output for uniform handling.
    sample_data = np.zeros(shape)
    sample_hat = np.fft.rfftn(sample_data)
    n_modes = sample_hat.size

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

    Uses the same convention as operators._get_wavenumbers:
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
            k = _get_wavenumbers(n, dx)
        else:
            # Other axes use full fft
            k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
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
    constraint_eqs: list[object],
) -> bool:
    """Check if all constraint equations can be eliminated in Fourier space.

    Requirements:
    - Each constraint self-operator must have exact Fourier multipliers
    - Constraint source terms must only reference fields/velocities with
      exact multipliers
    """
    for eq in constraint_eqs:
        for term in eq.rhs_terms:
            if term.operator not in _EXACT_MULTIPLIERS:
                return False
            if term.time_dependent:
                return False
    return True


def _build_constraint_eliminated_matrices(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: object,  # CoefficientEvaluator
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
) -> tuple[
    NDArray[np.complex128],  # A_reduced (n_modes, n_dyn, n_dyn)
    NDArray[np.complex128],  # recovery (n_modes, n_constraints, n_dyn)
    list[str],               # constraint_field_names
    dict[int, int],          # orig_to_reduced slot mapping
]:
    """Build reduced per-mode matrices with constraints algebraically eliminated.

    For a mixed system with dynamical (d) and constraint (c) fields:

        d/dt[d] = A_dd·d + A_dc·c     (dynamical equations)
        0       = S_cd·d + S_cc·c      (constraint: solve for c)

    The constraint gives c = -S_cc⁻¹·S_cd·d. Substituting:

        d/dt[d] = (A_dd - A_dc·S_cc⁻¹·S_cd)·d

    This handles v_A₀ references in dynamical equations by recognizing that
    v_A₀ = dA₀/dt = d/dt[-S_cc⁻¹·S_cd·d] = -S_cc⁻¹·S_cd·d', creating an
    implicit equation (I - A_dc_v·S_cc⁻¹·S_cd)·d' = (A_dd + A_dc_f·f)·d
    which is resolved by matrix inversion of the LHS factor.

    All operations are purely numeric (CoefficientEvaluator returns floats).
    S_cc⁻¹ in Fourier space is diagonal per mode — just 1/(m²+k²).

    Returns
    -------
    A_reduced : ndarray
        Per-mode matrices (n_modes, n_dyn, n_dyn) for dynamical fields only.
    recovery : ndarray
        Per-mode recovery (n_modes, n_constraints, n_dyn) for reconstructing
        constraint fields from dynamical state.
    constraint_field_names : list[str]
        Names of eliminated constraint fields.
    orig_to_reduced : dict
        Mapping from original layout slot index to reduced slot index.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    assert isinstance(coeff_eval, CoefficientEvaluator)

    n_modes = int(np.prod(rfft_shape))

    # Identify constraint and dynamical fields
    constraint_field_names: list[str] = []
    constraint_eq_map: dict[str, int] = {}  # field_name → eq_idx
    for eq_idx, eq in enumerate(spec.equations):
        if eq.time_derivative_order == 0:
            constraint_field_names.append(eq.field_name)
            constraint_eq_map[eq.field_name] = eq_idx
    n_c = len(constraint_field_names)

    # Build dynamical-only slot mapping (excluding constraint field slots)
    orig_to_reduced: dict[int, int] = {}
    red_idx = 0
    for si, slot in enumerate(layout.slots):
        if slot.kind == "constraint":
            continue
        orig_to_reduced[si] = red_idx
        red_idx += 1
    n_dyn = red_idx

    # Map field names to slot indices in the REDUCED layout
    dyn_slot_map: dict[str, int] = {}
    for si, slot in enumerate(layout.slots):
        if si in orig_to_reduced:
            dyn_slot_map[slot.name] = orig_to_reduced[si]
    # Also map velocity names v_X for dynamical fields
    for fname, si in layout.velocity_slot_map.items():
        v_name = f"v_{fname}"
        if si in orig_to_reduced:
            dyn_slot_map[v_name] = orig_to_reduced[si]

    # Constraint slot map (constraint field names → constraint index 0..n_c-1)
    c_idx_map: dict[str, int] = {
        name: i for i, name in enumerate(constraint_field_names)
    }

    # Evaluate Fourier multipliers
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op not in multiplier_cache:
                mult_fn = _EXACT_MULTIPLIERS[op]
                mult_val = mult_fn(k_grid)
                mult_full = np.broadcast_to(mult_val, rfft_shape)
                multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    # --- Build the four coupling matrices per mode ---

    # A_dd: dynamical → dynamical (n_modes, n_dyn, n_dyn)
    A_dd = np.zeros((n_modes, n_dyn, n_dyn), dtype=np.complex128)
    # A_dc: constraint → dynamical via FIELD references (n_modes, n_dyn, n_c)
    A_dc_field = np.zeros((n_modes, n_dyn, n_c), dtype=np.complex128)
    # A_dc_vel: constraint VELOCITY → dynamical (n_modes, n_dyn, n_c)
    A_dc_vel = np.zeros((n_modes, n_dyn, n_c), dtype=np.complex128)
    # S_cd: dynamical → constraint source (n_modes, n_c, n_dyn)
    S_cd = np.zeros((n_modes, n_c, n_dyn), dtype=np.complex128)
    # S_cc: constraint self-coupling (n_modes, n_c, n_c)
    S_cc = np.zeros((n_modes, n_c, n_c), dtype=np.complex128)

    for eq_idx, eq in enumerate(spec.equations):
        is_constraint = eq.time_derivative_order == 0
        is_second_order = eq.time_derivative_order >= 2

        if is_constraint:
            # Constraint equation: 0 = Σ coeff·op(target)
            ci = c_idx_map[eq.field_name]
            for term_idx, term in enumerate(eq.rhs_terms):
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx,
                )
                mult = multiplier_cache[term.operator]

                if term.field in c_idx_map:
                    # Self/cross constraint coupling
                    cj = c_idx_map[term.field]
                    S_cc[:, ci, cj] += coeff * mult
                elif term.field in dyn_slot_map:
                    # Source coupling to dynamical state
                    dj = dyn_slot_map[term.field]
                    S_cd[:, ci, dj] += coeff * mult

        elif is_second_order:
            field_slot = orig_to_reduced[layout.field_slot_map[eq.field_name]]
            vel_slot = orig_to_reduced[layout.velocity_slot_map[eq.field_name]]

            # Kinematic: dq/dt = v
            A_dd[:, field_slot, vel_slot] = 1.0

            # RHS terms: dv/dt = Σ coeff·op(target)
            for term_idx, term in enumerate(eq.rhs_terms):
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx,
                )
                mult = multiplier_cache[term.operator]

                if term.field in c_idx_map:
                    # References constraint field directly
                    cj = c_idx_map[term.field]
                    A_dc_field[:, vel_slot, cj] += coeff * mult
                elif term.field.startswith("v_") and term.field[2:] in c_idx_map:
                    # References constraint velocity v_A₀
                    cj = c_idx_map[term.field[2:]]
                    A_dc_vel[:, vel_slot, cj] += coeff * mult
                elif term.field in dyn_slot_map:
                    # Normal dynamical reference
                    dj = dyn_slot_map[term.field]
                    A_dd[:, vel_slot, dj] += coeff * mult

        else:
            # First-order: du/dt = Σ coeff·op(target)
            this_slot = orig_to_reduced[layout.field_slot_map[eq.field_name]]
            for term_idx, term in enumerate(eq.rhs_terms):
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=eq_idx, term_idx=term_idx,
                )
                mult = multiplier_cache[term.operator]

                if term.field in c_idx_map:
                    cj = c_idx_map[term.field]
                    A_dc_field[:, this_slot, cj] += coeff * mult
                elif term.field in dyn_slot_map:
                    dj = dyn_slot_map[term.field]
                    A_dd[:, this_slot, dj] += coeff * mult

    # --- Compute Schur complement ---

    # Invert S_cc per mode (small matrix, typically 1x1 or 2x2)
    S_cc_inv = np.zeros_like(S_cc)
    for m in range(n_modes):
        det = np.linalg.det(S_cc[m]) if n_c > 0 else 1.0
        if abs(det) < 1e-14:
            # Singular at k=0 (gauge freedom) — regularize
            S_cc_inv[m] = np.linalg.inv(
                S_cc[m] + 1e-14 * np.eye(n_c, dtype=np.complex128),
            )
        else:
            S_cc_inv[m] = np.linalg.inv(S_cc[m])

    # Recovery: c = -S_cc⁻¹ · S_cd · d
    # recovery[m, ci, dj] = -Σ_cj S_cc_inv[m,ci,cj] · S_cd[m,cj,dj]
    recovery = -np.einsum("mij,mjk->mik", S_cc_inv, S_cd)

    # Substitution: A_dc_field · c = A_dc_field · recovery · d
    # field_correction[m] = A_dc_field[m] @ recovery[m]
    field_correction = np.einsum("mij,mjk->mik", A_dc_field, recovery)

    # For constraint velocity: v_c = d/dt[c] = recovery · d'
    # where d' = A_reduced · d. So A_dc_vel · v_c = A_dc_vel · recovery · d'.
    # This creates implicit coupling:
    #   d' = A_dd · d + field_correction · d + A_dc_vel · recovery · d'
    #   (I - A_dc_vel · recovery) · d' = (A_dd + field_correction) · d
    #   d' = (I - A_dc_vel · recovery)⁻¹ · (A_dd + field_correction) · d
    vel_coupling = np.einsum("mij,mjk->mik", A_dc_vel, recovery)

    # Check if vel_coupling is nonzero (constraint velocity referenced)
    has_vel_coupling = np.max(np.abs(vel_coupling)) > 1e-15

    A_rhs = A_dd + field_correction

    if has_vel_coupling:
        # Implicit solve: (I - vel_coupling) · d' = A_rhs · d
        eye = np.broadcast_to(
            np.eye(n_dyn, dtype=np.complex128), (n_modes, n_dyn, n_dyn),
        ).copy()
        lhs = eye - vel_coupling
        # A_reduced = lhs⁻¹ · A_rhs (per mode)
        A_reduced = np.zeros((n_modes, n_dyn, n_dyn), dtype=np.complex128)
        for m in range(n_modes):
            A_reduced[m] = np.linalg.solve(lhs[m], A_rhs[m])
    else:
        A_reduced = A_rhs

    return A_reduced, recovery, constraint_field_names, orig_to_reduced


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
                target_slot = layout.field_slot_map[term.field]
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=_eq_idx, term_idx=_term_idx,
                )
                mult = multiplier_cache[term.operator]
                A[:, vel_slot, target_slot] += coeff * mult

        else:
            # First-order: du/dt = Σ coeff * operator(target_field)
            this_slot = layout.field_slot_map[field_name]
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = layout.field_slot_map[term.field]
                coeff = _resolve_constant_coeff(
                    term, coeff_eval, eq_idx=_eq_idx, term_idx=_term_idx,
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
    different k-modes, producing a full (n_total × n_total) matrix where
    n_total = n_slots × n_modes.

    For Gaussian c(x), the convolution kernel ĉ(q) decays exponentially,
    making the matrix effectively banded.

    Reference: Burns et al. (2020), Phys. Rev. Research 2:023068.
    """
    n_slots = layout.num_slots
    n_modes = int(np.prod(rfft_shape))
    n_total = n_slots * n_modes

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
                target_slot = layout.field_slot_map[term.field]
                mult = multiplier_cache[term.operator]

                if not term.position_dependent:
                    # Constant coefficient: diagonal in mode space
                    coeff = _resolve_constant_coeff(
                        term, coeff_eval, eq_idx=_eq_idx, term_idx=_term_idx,
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
                target_slot = layout.field_slot_map[term.field]
                mult = multiplier_cache[term.operator]

                if not term.position_dependent:
                    coeff = _resolve_constant_coeff(
                        term, coeff_eval, eq_idx=_eq_idx, term_idx=_term_idx,
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
            probe_hat.reshape(rfft_shape), s=grid.shape,
            axes=list(range(len(grid.shape))),
        )
        product = coeff_array * probe_physical
        result_hat = np.fft.rfftn(product).ravel()

        # result_hat[m] = Σ_{k'} (1/N) ĉ_{m-k'} δ_{k',m'} = (1/N) ĉ_{m-m'}
        # multiplied by operator multiplier at m'
        for m in range(n_modes):
            row = row_slot * n_modes + m
            col = col_slot * n_modes + m_prime
            A[row, col] += result_hat[m] * operator_mult[m_prime]


# ---------------------------------------------------------------------------
# Block decomposition
# ---------------------------------------------------------------------------


def _find_independent_blocks(
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


def _is_all_second_order(spec: EquationSystem) -> bool:
    """Check if all equations are second-order (wave-type)."""
    return all(eq.time_derivative_order >= 2 for eq in spec.equations)


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
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evolve system with per-mode independent matrices (constant coefficients).

    A_modes has shape (n_modes, n_slots, n_slots).
    y0_hat has shape (n_slots, n_modes).

    Uses block-aware eigendecomposition: independent field blocks are detected
    and eigendecomposed separately to prevent degenerate-eigenvalue mixing.
    Blocks with all-zero initial conditions are skipped entirely.

    Ref: Golub & Van Loan (1996), Matrix Computations, §4.8.
    """
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
    blocks = _find_independent_blocks(combined)

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

        # Batch eigendecomposition for this block
        eig_vals, V = np.linalg.eig(A_block)
        V_inv = np.linalg.inv(V)

        # Warn about potential overflow
        _warn_eigenvalue_growth(eig_vals, dt_total, context="per-mode")

        # Transform IC to eigenbasis
        y0_eigen = np.einsum("mij,mj->mi", V_inv, y0_block.T)

        block_data.append((block_slots, eig_vals, V, y0_eigen))

    # Evolve at each time point.
    # Pre-multiply V @ diag(y0_eigen) for each block so the inner loop only
    # needs element-wise exp + matrix-vector product, not a full einsum.
    block_evolved: list[
        tuple[
            list[int],  # slot indices
            NDArray[np.complex128],  # V_y0: V * y0_eigen, (n_modes, bs, bs)
            NDArray[np.complex128],  # eigenvalues (n_modes, bs)
        ]
    ] = []
    for block_slots, eig_vals, V, y0_eigen in block_data:
        # V_y0[m, i, j] = V[m, i, j] * y0_eigen[m, j]
        # so y(t) = V_y0 @ exp(λ*dt) is just a matvec
        V_y0 = V * y0_eigen[:, np.newaxis, :]  # (n_modes, bs, bs)
        block_evolved.append((block_slots, V_y0, eig_vals))

    snapshots = np.zeros((n_snapshots, n_slots * n_pts))
    times = np.zeros(n_snapshots)
    n_modes = y0_hat.shape[1]

    for ti, t in enumerate(t_eval):
        dt = t - t0
        y_hat_t = np.zeros((n_slots, n_modes), dtype=np.complex128)

        for block_slots, V_y0, eig_vals in block_evolved:
            # exp_lambda shape: (n_modes, block_size)
            exp_lambda = np.exp(eig_vals * dt)
            # y_evolved[m, i] = Σ_j V_y0[m, i, j] * exp(λ_j * dt)
            y_evolved = np.einsum("mij,mj->mi", V_y0, exp_lambda)
            y_hat_t[block_slots, :] = y_evolved.T

        y_physical = _ifft_slots(y_hat_t, layout, grid)
        snapshots[ti] = y_physical
        times[ti] = t

        if snapshot_callback is not None:
            snapshot_callback(t, y_physical)

        if progress is not None:
            progress.update(t)

    return times, snapshots


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

    A_full has shape (n_total, n_total) where n_total = n_slots × n_modes.
    y0_hat has shape (n_slots, n_modes).

    Uses eigendecomposition of the full matrix for exact evolution.
    For large matrices (>2000), uses Schur decomposition for stability.
    Ref: Moler & Van Loan (2003), SIAM Review 45(1):3-49.
    """
    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_modes = y0_hat.shape[1]
    n_total = A_full.shape[0]
    n_snapshots = len(t_eval)

    # Flatten y0_hat to (n_total,) — interleave as [slot0_mode0, slot0_mode1, ...]
    y0_flat = y0_hat.ravel()  # (n_slots * n_modes,) — slot-major order

    # Choose algorithm based on matrix size
    # Ref: Moler & Van Loan (2003), Table 1
    SCHUR_THRESHOLD = 2000

    if n_total <= SCHUR_THRESHOLD:
        # Dense eigendecomposition
        eigenvalues, V = np.linalg.eig(A_full)
        V_inv = np.linalg.inv(V)
        y0_eigen = V_inv @ y0_flat
    else:
        # Schur decomposition for numerical stability on large matrices
        from scipy.linalg import schur  # noqa: PLC0415

        T, Z = schur(A_full, output="complex")
        eigenvalues = np.diag(T)
        V = Z
        V_inv = Z.conj().T  # Unitary for Schur
        y0_eigen = V_inv @ y0_flat

    # Clean up eigencoefficient roundoff: entries far below the dominant
    # IC amplitude are numerical noise that can seed exponential growth
    # in systems with positive-real-part eigenvalues.
    y0_max = np.max(np.abs(y0_eigen))
    if y0_max > 0:
        cleanup_threshold = y0_max * 1e-14
        y0_eigen[np.abs(y0_eigen) < cleanup_threshold] = 0.0

    # Warn about potential overflow from growing eigenvalues
    dt_total = float(t_eval[-1] - t_eval[0])
    _warn_eigenvalue_growth(eigenvalues, dt_total, context="full-matrix")

    snapshots = np.zeros((n_snapshots, n_slots * n_pts))
    times = np.zeros(n_snapshots)

    t0 = t_eval[0]
    for ti, t in enumerate(t_eval):
        dt = t - t0  # elapsed time from IC
        exp_lambda = np.exp(eigenvalues * dt)
        y_evolved = V @ (exp_lambda * y0_eigen)

        # Reshape from (n_total,) = (n_slots * n_modes,) back to (n_slots, n_modes)
        y_hat_t = y_evolved.reshape(n_slots, n_modes)

        # IFFT back to physical space
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
    has_constraints = any(
        eq.time_derivative_order == 0 for eq in spec.equations
    )
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

    has_pos_dep = _has_position_dependent_terms(spec)

    if has_constraints and not has_pos_dep:
        # Constraint elimination via Fourier Schur complement
        # Ref: Hairer & Wanner (1996), Solving ODEs II, Ch. VII
        (
            A_reduced, recovery_matrix, c_names, orig_to_reduced,
        ) = _build_constraint_eliminated_matrices(
            spec, layout, grid, coeff_eval, k_grid, rfft_shape,
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

        # Evolve dynamical fields
        times, dyn_snapshots = _evolve_per_mode(
            A_reduced, y0_hat_dyn, t_eval, dyn_layout, grid,
            None, progress,  # callback handled below with full state
        )

        # Reconstruct full state (including constraints) at each snapshot
        n_full = layout.num_slots * n_pts
        snapshots = np.zeros((len(t_eval), n_full))

        for ti in range(len(t_eval)):
            dyn_phys = dyn_snapshots[ti]
            # Re-FFT dynamical fields for constraint recovery
            y_hat_dyn_t = np.zeros((n_dyn, n_modes), dtype=np.complex128)
            for di in range(n_dyn):
                y_hat_dyn_t[di] = np.fft.rfftn(
                    dyn_phys[di * n_pts : (di + 1) * n_pts].reshape(
                        grid.shape,
                    ),
                ).ravel()

            # Recover constraint fields: c_hat = recovery @ d_hat
            c_hat = np.einsum("mcj,jm->cm", recovery_matrix, y_hat_dyn_t)

            # Assemble full physical state
            full_state = np.zeros(n_full)
            for orig_si, red_pos in orig_to_reduced.items():
                full_state[orig_si * n_pts : (orig_si + 1) * n_pts] = (
                    dyn_phys[red_pos * n_pts : (red_pos + 1) * n_pts]
                )
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

            snapshots[ti] = full_state
            if snapshot_callback is not None:
                snapshot_callback(t_eval[ti], full_state)

        n_c = len(c_names)
        method_desc = (
            f"per-mode eigendecomposition with Schur constraint elimination "
            f"({n_c} constraints, {n_dyn} dynamical slots)"
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
        times, snapshots = _evolve_per_mode(
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
        method_desc = (
            f"full eigendecomposition ({n_total}×{n_total}, position-dependent)"
        )

    if progress is not None:
        progress.finish()

    return {
        "t": times,
        "y": snapshots,
        "success": True,
        "message": f"Modal solver completed ({method_desc})",
    }
