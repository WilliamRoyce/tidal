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

Kinetic-coefficient (mass-matrix) handling — four primary dispatches
--------------------------------------------------------------------
Equations with a non-trivial ``kinetic_coefficient_symbolic`` (``M ẍ = K x``
rather than ``ẍ = K x``) reach the solver via **four distinct code paths**
plus the Pass 1 correction-matrix builder. All five must consume ``M``
identically, or cross-path regressions appear silently. Two such regressions
to date:

- **GH #367** (fixed v0.42.0, 2026-05-19): path 3 (pos-dep convolution)
  silently dropped ``M⁻¹`` for theories with non-trivial
  ``kinetic_coefficient_symbolic``.
- **GH #379** (fixed 2026-05-24): path 4 (the pos-dep + constraint case)
  did not exist; the dispatch silently routed those theories to path 3
  which has no constraint handling, causing 10²⁵² divergence at any
  nonzero coupling. Now handled by
  :func:`_build_convolution_matrix_with_constraints` (convolution analog
  of path 2's Schur elimination).

1. **Fast path** — :func:`_build_per_mode_matrices` (called when
   ``needs_reduction = False`` *and* no position-dependent coefficients).
   Builds the first-order evolution matrix directly as
   ``A[velocity_slot, target_slot] = M⁻¹(field) · coeff · multiplier``.
   M⁻¹ is pre-computed once from
   :func:`tidal.solver._kinetic.build_inverse_kinetic_diag` and folded in
   via the shared :func:`tidal.solver._kinetic.velocity_row_scale` helper.

2. **Generalized-eig path** — :func:`_build_evolution_matrices` (called
   when the spec has constraints or time-derivative RHS operators, and no
   position-dependent coefficients). Populates ``M_mat`` from
   ``kinetic_coefficient_symbolic`` on the diagonal and solves
   ``(K − λM) v = 0`` via ``scipy.linalg.eig(A, B)`` with QZ decomposition.
   Handles rank-deficient M via null-space projection and Schur elimination
   for hidden algebraic constraints. *Structurally distinct* from paths 1,
   3, and 4: writes ``M`` to the generalized eigenvalue problem rather than
   post-multiplying a velocity-row contribution, so it intentionally does
   NOT call :func:`velocity_row_scale`.

3. **Position-dependent convolution path** — :func:`_build_convolution_matrix`
   (called when the spec has position-dependent coefficients — RHS *or*
   kinetic, GH #427 — and *no* algebraic constraints). Builds
   ``A[m, m'] = ĉ(m − m') · op_mult(k_{m'})`` over the full Fourier-mode
   block, threading ``M⁻¹`` through both the constant-coefficient and the
   convolution-coupling paths via :func:`velocity_row_scale`. A
   position-dependent kinetic yields an ndarray scale — a REAL-SPACE
   ``M⁻¹(x)`` profile folded into the coefficient before the FFT
   (:func:`_conv_block_with_kinetic`); it must never multiply a k-space
   block. Rejects (via assertion) systems with
   ``time_derivative_order == 0`` equations to make the GH #379 bug class
   a static error.

4. **Position-dependent + constraint convolution path** —
   :func:`_build_convolution_matrix_with_constraints` (called when both
   pos-dep coefficients *and* algebraic constraints are present). Builds
   convolution sub-matrices ``A_dd``, ``A_dc_field/vel``, ``K_cd``, ``K_cc``
   and Schur-eliminates: ``recovery = -K_cc⁻¹ · K_cd``, ``A_reduced = A_dd
   + A_dc_field · recovery`` (plus B_lhs pre-solve for v_constraint and
   Newton correction for ẍ_constraint references in dyn RHS). Mirrors
   path 2's structure but with full 2D convolution matrices instead of
   per-mode 3-tensors. Applies M⁻¹ on dynamical velocity rows via
   :func:`velocity_row_scale`, same contract as paths 1 and 3.

Plus the Pass 1 correction-matrix builder
:func:`_build_pass1_source_matrices` (the perturbative Duhamel correction
path; see :class:`tidal.solver.perturbative_driver.PerturbativeSolver`),
which also emits velocity-row entries into a per-mode source matrix and
likewise consumes ``M⁻¹`` through :func:`velocity_row_scale`.

**Cross-builder contract** (enforced by the regression test
:class:`tests.test_solver_kinetic_consistency.TestAllModalPathsRespectKinetic`):
any matrix builder that populates ``A[velocity_slot, target_slot]`` entries
for a second-order equation **MUST** apply
:func:`tidal.solver._kinetic.velocity_row_scale(field_name, m_inv)` to the
contribution. Path 2 (generalized-eig) is the documented exception —
it writes ``M`` to a separate matrix that participates in the eigenvalue
problem. Adding a new matrix builder without satisfying this contract was
exactly the GH #367 regression; please don't repeat it. Path 4 (pos-dep +
constraint) follows the same contract.

The choice between dispatches is transparent to callers of :func:`solve_modal`
but load-bearing for reviewers: all five sites must be updated together
when the kinetic handling changes. The time-domain backends (cvode/ida/
leapfrog/scipy) use :func:`tidal.solver._kinetic.build_inverse_kinetic_diag`
as a shared entry point, keeping the six backends (modal-fast,
modal-genEig, modal-conv, modal-conv-Schur, cvode, ida, leapfrog/scipy)
numerically consistent per :mod:`tests.test_solver_kinetic_consistency`.

Algorithm paths for coefficient structure
-----------------------------------------
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
# ruff: noqa: DOC501, DOC502 — _evolve_per_mode_pade raises via nested helper;
#   _evolve_per_mode wrapper documents the exception from the delegated call.
# ruff: noqa: PLW2901 — y_curr intentionally overwritten in snapshot loop.

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING, Any, cast

import numpy as np
from numpy.typing import NDArray

from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL
from tidal.solver._exceptions import KineticEvaluationError
from tidal.solver._setup import warn_frozen_constraints
from tidal.solver.operators import get_wavenumbers, is_periodic_bc
from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from tidal.solver._types import PerturbativePass1Result, SolverResult
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
    # --- Fourth-time-derivative operators (only valid as correction
    # sources in the v6 perturbative pipeline — never in the base
    # system, which base_spec demotes to time_order <= 2).  Spatial
    # Fourier multiplier is 1 for the pure-time d4_t and the usual
    # ik/−k² factors for the mixed variants. ---
    "d4_t": _OperatorDecomp(lambda k_axes: np.ones_like(k_axes[0]), 4),
    "mixed_T4_S1x": _OperatorDecomp(lambda k_axes: 1j * k_axes[0], 4),
    "mixed_T4_S1y": _OperatorDecomp(lambda k_axes: 1j * k_axes[1], 4),
    "mixed_T4_S1z": _OperatorDecomp(lambda k_axes: 1j * k_axes[2], 4),
    "mixed_T4_S2x": _OperatorDecomp(lambda k_axes: -(k_axes[0] ** 2), 4),
    "mixed_T4_S2y": _OperatorDecomp(lambda k_axes: -(k_axes[1] ** 2), 4),
    "mixed_T4_S2z": _OperatorDecomp(lambda k_axes: -(k_axes[2] ** 2), 4),
}

# Backward-compatible mapping: operator name → spatial multiplier function.
# Used by existing code paths that only need the spatial part.
_EXACT_MULTIPLIERS: dict[str, _ExactMultFn] = {
    name: dec.spatial_fn for name, dec in _OPERATOR_DECOMP.items()
}


# ---------------------------------------------------------------------------
# Closed-form Duhamel kernel for Pass 1 of the v6 perturbative solver.
# ---------------------------------------------------------------------------
# Given the inhomogeneous linear ODE
#     dy/dt = A·y + M_src·y0(t),   y(0) = 0,
# where y0(t) is the Pass 0 solution (a sum of eigenmode exponentials)
# and A is the same base operator from Pass 0, the solution in the
# eigenbasis of A (A = V·D·V^-1) is
#     z_i(t) = Σ_j β_ij · α_j · G(λ_i, λ_j; t),
# with β = V^-1·M_src·V, α = V^-1·y0(0), and
#     G(λ, μ; t) = [exp(μt) − exp(λt)] / (μ − λ)          (μ ≠ λ)
#     G(λ, λ; t) = t·exp(λt)                              (degenerate)
# The naive μ ≠ λ formula suffers catastrophic cancellation when
# |μ − λ|·t is small. Follow Al-Mohy & Higham (2011), SIAM J. Sci.
# Comput. 33, §3 (Table 3.1): switch to a Taylor expansion of the
# entire function φ₁(z) = (exp(z) − 1)/z = Σ_{k≥0} z^k/(k+1)! whenever
# |(μ−λ)·t| is smaller than a conservative threshold. The crossover
# point for double precision is ~ √ε_mach ≈ 1.5e-8; setting the
# threshold an order of magnitude higher (1e-5) leaves both branches
# accurate while avoiding edge-case switching noise.

_PHI1_DEGENERACY_THRESHOLD = 1e-5
"""Switch-over |(μ−λ)·t| for Duhamel kernel Taylor fallback.

Cites Al-Mohy & Higham 2011 §3 Table 3.1. Set one order of magnitude
above √ε_mach ≈ 1.5e-8 so both the direct and Taylor branches evaluate
in their numerically stable regime.
"""


def _duhamel_kernel(  # pyright: ignore[reportUnusedFunction]
    lam: complex | NDArray[np.complex128],
    mu: complex | NDArray[np.complex128],
    t: float,
) -> complex | NDArray[np.complex128]:
    """Closed-form Duhamel kernel G(λ, μ; t) for linear inhomogeneous ODEs.

    Evaluates ``G(λ, μ; t) = [exp(μ t) − exp(λ t)] / (μ − λ)`` with a
    12-term Taylor fallback on the φ₁ entire function near μ = λ to
    prevent catastrophic cancellation. Accepts scalar or
    broadcast-compatible complex array inputs; returns the element-wise
    kernel.

    Parameters
    ----------
    lam, mu : complex or complex ndarray
        Eigenvalues λ (row index in the source formula) and μ
        (column index). Typically arrays of shape ``(bs,)`` per block
        or ``(n_modes, bs)``.
    t : float
        Evaluation time.

    Returns
    -------
    G : complex or complex ndarray
        Same shape as ``broadcast(lam, mu)``.

    References
    ----------
    Al-Mohy & Higham 2011, "Computing the action of the matrix
    exponential", SIAM J. Sci. Comput. 33, §3 Table 3.1.
    """
    lam_arr = np.asarray(lam, dtype=np.complex128)
    mu_arr = np.asarray(mu, dtype=np.complex128)
    z = (mu_arr - lam_arr) * t
    e_lam_t = np.exp(lam_arr * t)

    # Direct formula: G = exp(λt) · (exp(z) − 1) / (μ − λ) = exp(λt) · t ·
    # expm1(z) / z. Using expm1 preserves accuracy for z near 0, but the
    # 1/z blows up near degeneracy — hence the Taylor branch below.
    direct_denom = np.where(np.abs(z) > 0, z, np.ones_like(z))
    with np.errstate(divide="ignore", invalid="ignore"):
        direct = e_lam_t * t * np.expm1(z) / direct_denom

    # Taylor expansion of φ₁(z) = Σ_{k≥0} z^k / (k+1)! with G = t·exp(λt)·φ₁(z).
    # Twelve terms saturate double precision for |z| ≲ 1; this branch
    # runs only for |z| < _PHI1_DEGENERACY_THRESHOLD so truncation error
    # is well below 1e-13 in practice.
    phi1 = np.ones_like(z)
    z_power = np.ones_like(z)
    for k in range(1, 12):
        z_power = z_power * z / (k + 1)
        phi1 += z_power
    taylor = t * e_lam_t * phi1

    mask = np.abs(z) > _PHI1_DEGENERACY_THRESHOLD
    result = np.where(mask, direct, taylor)
    if np.ndim(lam) == 0 and np.ndim(mu) == 0:
        return complex(result)
    return result


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
    6. No position-dependent kinetic coefficients (GH #421) — position-
       dependent RHS coefficients are fine (handled via the ĉ(k−k′)
       convolution paths), but a position-dependent M(x) needs the
       mass-side convolution M̂(k−k′), which is not implemented (GH #427).
       Auto-selection falls through to the time-domain backends, which
       evaluate such kinetics on the grid (GH #382).
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

    # Position-dependent kinetic coefficients are supported (GH #427): the
    # kinetics-aware routing predicate sends such specs to the convolution
    # builders, which fold M⁻¹(x) into the real-space coefficients. The
    # former requirement 6 (GH #421 refusal) is retired.
    return True


def _kinetic_eval_error(
    field_name: str,
    kin_sym: str,
    coeff_eval: object,
    exc: Exception,
) -> KineticEvaluationError:
    """Build the GH #447 refusal for an unresolvable kinetic coefficient.

    Names the field and the expression, and — for the common case of an
    unbound symbol — the missing name and the parameters that WERE
    supplied, since "forgot a ``--param``" is the failure this most often
    represents. ``evaluate_coefficient`` wraps the original exception, so
    the unbound-name case is recovered from the ``__cause__`` chain.
    """
    cause: BaseException | None = exc
    while cause is not None and not isinstance(cause, NameError):
        cause = cause.__cause__
    hint = ""
    if isinstance(cause, NameError):
        params = getattr(coeff_eval, "_parameters", {}) or {}
        supplied = ", ".join(sorted(params)) or "(none)"
        hint = (
            f" The expression references a symbol that has no value: "
            f"{cause}. Parameters supplied: {supplied}. Pass the missing "
            f"one with --param NAME=VALUE (or add it to the theory's "
            f"[parameters] section)."
        )
    msg = (
        f"Cannot evaluate the kinetic coefficient for field {field_name!r} "
        f"(expression: {kin_sym!r}): {exc}.{hint} The modal mass matrix "
        f"needs a concrete M; continuing with M = 1 would silently solve a "
        f"different theory, so this is refused (GH #447)."
    )
    return KineticEvaluationError(msg)


def _raise_position_dependent_kinetic(spec: EquationSystem, where: str) -> None:
    """Raise the per-mode-regime refusal for a position-dependent kinetic coefficient.

    Post-GH #427, plain :func:`solve_modal` supports M(x) — the routing
    predicate sends such specs to the convolution builders, which fold the
    per-grid-point M⁻¹(x) into the real-space coefficients. This guard
    remains only on the strictly PER-MODE engines, where a position-
    dependent M has no per-k representation:

    * :func:`_build_evolution_matrices` — the genEig/Schur engine writes M
      onto the per-mode generalized-eigenvalue problem ``(K − λM)v = 0``,
      which ceases to exist when M(x) couples modes. Reached directly by
      the stability probe (:mod:`tidal.measurement._stability`, every
      gated sweep point / likelihood evaluation) and :mod:`modal_jax`.
    * :func:`solve_modal_pass1` — the per-mode Duhamel correction path
      (no eigendata exists on the convolution route).

    Failing loudly here replaces the pre-GH #421 behaviors: a grid-less
    ``evaluate_coefficient`` ValueError on some paths, and a silent
    ``M = 1`` fallback on the genEig path.
    """
    posdep_fields = [
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order > 0 and eq.kinetic_position_dependent
    ]
    msg = (
        f"{where} operates per k-mode and cannot represent a position-"
        f"dependent kinetic coefficient (field(s) {posdep_fields}): "
        "M(x) couples Fourier modes, so no per-mode mass matrix exists "
        "(GH #421 is this guard; GH #441 tracks a mode-coupling-aware "
        "stability gate). Alternatives: (a) plain solve_modal / "
        "auto-selection — the convolution modal path supports M(x) "
        "(GH #427); (b) --scheme cvode or --scheme ida — the time-domain "
        "solvers evaluate position-dependent kinetics on the grid "
        "(GH #382); (c) for perturbative flows, declare the position-"
        "dependence's small parameter in [perturbation].small_parameters "
        "so kinetic canonicalization keeps the base spec constant-"
        "coefficient (GH #380, the Euler-Heisenberg flow)."
    )
    raise NotImplementedError(msg)


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


def _fft_slots_full(
    y: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
) -> NDArray[np.complex128]:
    """Transform each slot to Fourier space over the FULL complex spectrum.

    GH #445: the position-dependent convolution paths operate in the full
    fftn basis (all N frequencies per axis, ``n_modes = prod(grid.shape)``)
    because pointwise multiplication by c(x) is exactly C-linear there,
    whereas the rfft half-spectrum basis has no imaginary degree of
    freedom at the DC/Nyquist bins and misrepresents sin-phase content that the
    coefficient's harmonics fold onto them. The constant-coefficient
    per-mode paths keep the cheaper rfft basis (:func:`_fft_slots`) — a
    constant coefficient is diagonal in k, so no folding occurs.
    """
    n_slots = layout.num_slots
    n_pts = layout.num_points
    shape = grid.shape
    n_modes = int(np.prod(shape))

    y_hat = np.zeros((n_slots, n_modes), dtype=np.complex128)
    for slot_idx in range(n_slots):
        start = slot_idx * n_pts
        field_data = y[start : start + n_pts].reshape(shape)
        y_hat[slot_idx] = np.fft.fftn(field_data).ravel()
    return y_hat


def _ifft_slots_full(
    y_hat: NDArray[np.complex128],
    layout: StateLayout,
    grid: GridInfo,
) -> NDArray[np.float64]:
    """Inverse of :func:`_fft_slots_full`; discards residual imaginary parts.

    Real coefficients and real/conjugate-symmetric operators keep the
    evolution conjugate-symmetric, so the imaginary residue is roundoff.
    """
    n_slots = layout.num_slots
    n_pts = layout.num_points
    shape = grid.shape

    y_out = np.zeros(n_slots * n_pts)
    for slot_idx in range(n_slots):
        hat_data = y_hat[slot_idx].reshape(shape)
        physical = np.fft.ifftn(hat_data)
        y_out[slot_idx * n_pts : (slot_idx + 1) * n_pts] = physical.real.ravel()
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


def _zero_nyquist_full(
    y_hat: NDArray[np.complex128],
    grid: GridInfo,
) -> None:
    """Zero the Nyquist bin(s) of full-fftn slot spectra, in place.

    Same rationale as the rfft-basis zeroing in :func:`solve_modal`
    (Boyd 2001 §11.5): the Nyquist mode aliases with its conjugate and
    cannot represent physical content; complex evolution entries would
    otherwise create imaginary Nyquist components that the real-part
    projection silently drops, losing energy. In the full basis the
    Nyquist bin is the single ``fftfreq`` index ``n // 2`` per even axis.
    """
    shape = grid.shape
    spectral = y_hat.reshape(y_hat.shape[0], *shape)
    for ax, n in enumerate(shape):
        if n % 2 == 0:
            idx: list[Any] = [slice(None)] * (len(shape) + 1)
            idx[ax + 1] = n // 2
            spectral[tuple(idx)] = 0.0


def _build_k_axes_full(grid: GridInfo) -> list[NDArray[np.float64]]:
    """Wavenumber arrays with the FULL fft convention on every axis.

    GH #445: companion to :func:`_fft_slots_full` for the position-
    dependent convolution paths — ``k = 2π · fftfreq(N, d=dx)`` on all
    axes, matching ``np.fft.fftn`` bin ordering.
    """
    return [
        np.asarray(
            2.0 * np.pi * np.fft.fftfreq(grid.shape[ax], d=grid.dx[ax]),
            dtype=np.float64,
        )
        for ax in range(grid.ndim)
    ]


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

    **Pos-dep coefficients in constraints**: supported as of GH #379. The
    modal solver routes such theories to
    :func:`_build_convolution_matrix_with_constraints`, which performs
    Schur elimination on full convolution matrices. Do not add a pos-dep
    rejection here without a fresh bug — the previous regime where this
    function effectively gated pos-dep + constraints out of the modal
    path was exactly the #379 bug class (dispatch silently routed to the
    wrong builder rather than rejecting the theory).
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
    diagnostics: PencilDiagnostics | None = None,
) -> tuple[
    NDArray[np.complex128],  # A_rhs (n_modes, n_dyn_slots, n_dyn_slots)
    NDArray[np.complex128] | None,  # B_lhs (n_modes, n_dyn, n_dyn) or None
    NDArray[np.complex128],  # recovery (n_modes, n_total_constraints, n_dyn_slots)
    NDArray[np.complex128] | None,  # v_recovery (n_modes, n_c, n_dyn) or None
    list[str],  # all constraint field names
    dict[int, int],  # orig_to_reduced slot mapping
    NDArray[np.complex128] | None,  # S_cc_inv (n_modes, n_c, n_c) or None
    NDArray[np.bool_] | None,  # singular_mask (n_modes,) or None
    NDArray[np.complex128] | None,  # manifold projector (n_modes, n, n) or None
]:
    """Build the per-mode evolution generator for constraint-coupled systems.

    Handles the generalized second-order system

        M(k)·ẍ = K(k)·x + D(k)·ẋ + J(k)·x⃛

    plus residual algebraic rows, uniformly (GH #457 rework):

    1. The promoted second-order sector (``spec.second_order_sector``) —
       order-0 rows carrying inter-constraint time derivatives — joins the
       dynamical field set with pure off-diagonal mass rows; only RESIDUAL
       rows are Schur-eliminated (``recovery = −S_cc⁻¹·S_cd``).
    2. The dynamical sector is assembled as a first-order PENCIL
       ``B·ẏ = A·y`` with all mass content on the B side and raw K/D on
       the A side; residual couplings fold in as
       ``A += A_dc_field·recovery``, ``B −= A_dc_vel·recovery`` — scale-
       correct by construction (GH #458).
    3. Composition: batched ``B⁻¹A`` when B is invertible on every mode;
       otherwise the ordered-QZ deflating-subspace engine
       (:func:`_pencil_deflate`), which handles singular mass matrices,
       k-dependent rank deficiency, and constraint chains of any index,
       and raises :class:`SingularPencilError` on true gauge freedom.

    Returns
    -------
    A_rhs : ndarray, shape (n_modes, n_dyn_slots, n_dyn_slots)
        Fully composed first-order evolution generator.
    B_lhs : None
        Composition is complete; kept for tuple-shape compatibility.
    recovery : ndarray, shape (n_modes, n_c, n_dyn_slots)
        Residual-constraint recovery.
    v_recovery : ndarray or None
        Exact constraint-velocity recovery (``recovery @ A_rhs``).
    constraint_field_names : list[str]
        Names of the RESIDUAL algebraic constraint fields.
    orig_to_reduced : dict[int, int]
        Original → reduced slot index mapping.
    S_cc_inv, singular_mask :
        Residual Schur internals for the Pass-1 augmented recovery (#290).
    manifold_proj : ndarray or None
        Per-mode orthogonal projector onto the constraint manifold when
        the QZ engine ran (project the IC with it); None on the batched
        fast path.
    """
    import logging  # noqa: PLC0415

    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    assert isinstance(coeff_eval, CoefficientEvaluator)
    logger = logging.getLogger(__name__)

    # GH #421 guard AT THE BUILDER, not only at solve_modal: this builder
    # has five production call sites (solve_modal, modal_jax, and three in
    # tidal/measurement/_stability.py — the stability probe that gates
    # every likelihood evaluation in an inference chain).  Without this
    # check the inline kinetic evaluation below hits `except Exception`
    # and silently proceeds with M = 1 — wrong physics with no error.
    # NotImplementedError deliberately does NOT match the probe's
    # (LinAlgError, ValueError) handler, so it propagates instead of
    # being mislabeled "tachyonic".
    if spec.has_position_dependent_kinetic():
        _raise_position_dependent_kinetic(spec, "_build_evolution_matrices")

    n_modes = int(np.prod(rfft_shape))

    # ---- Identify constraint and dynamical fields ----
    # GH #457: order-0 rows carrying inter-constraint time derivatives
    # (the promoted second-order sector, per the ONE classification in
    # EquationSystem.second_order_sector) are NOT algebraic — they join
    # the dynamical sector below. Only the residual rows are Schur-
    # eliminated here.
    promoted_fields = spec.second_order_sector.promoted
    constraint_field_names: list[str] = [
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order == 0 and eq.field_name not in promoted_fields
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
    # Map: dynamical field name → index in the n_f dynamical field array.
    # Promoted order-0 fields (GH #457) are dynamical-sector members: they
    # have field + velocity slots from the layout and mass-matrix ROWS
    # (pure off-diagonal, from their d2_t cross terms) below.
    dyn_field_names: list[str] = []
    dyn_field_idx: dict[str, int] = {}
    for eq in spec.equations:
        if eq.time_derivative_order > 0 or eq.field_name in promoted_fields:
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
        eq.field_name: eq
        for eq in spec.equations
        if eq.time_derivative_order > 0 or eq.field_name in promoted_fields
    }
    for fi, fname in enumerate(dyn_field_names):
        eq_f = eq_for_field[fname]
        if fname in promoted_fields:
            # GH #457: a promoted row has an ALGEBRAIC LHS — there is no
            # ẍ_self term, so its mass-matrix diagonal is exactly 0 (its
            # mass content is the off-diagonal d2_t cross terms collected
            # in the term loop below). The M = 1 JSON convention applies
            # only to true order > 0 equations.
            continue
        kin_sym = getattr(eq_f, "kinetic_coefficient_symbolic", None)
        if kin_sym is None:
            # No kinetic declared — M = 1 by JSON convention (not a fallback).
            M_mat[:, fi, fi] = 1.0
        else:
            # GH #447: no silent fallback. This builder previously caught
            # every exception and proceeded with M = 1, so a kinetic
            # referencing a parameter absent from `--param` produced wrong
            # amplitudes with exit code 0 — and did so inside the stability
            # probe that gates every sweep point / likelihood evaluation.
            # There is no correct fallback value: M = 1 is a different
            # theory. Raise KineticEvaluationError (a RuntimeError, so the
            # probe's (LinAlgError, ValueError) handler cannot relabel it
            # "tachyonic") with the field, the expression and the cause.
            try:
                kin_val = evaluate_coefficient(
                    kin_sym,
                    coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
                    spec.effective_coordinates,
                )
            except (ValueError, TypeError, ArithmeticError) as exc:
                raise _kinetic_eval_error(fname, kin_sym, coeff_eval, exc) from exc
            if isinstance(kin_val, np.ndarray):
                # Position-dependent kinetics are refused by the guard at
                # the top of this function, so an array here means the
                # coefficient varies over some other axis (e.g. time) that
                # the per-mode mass matrix cannot represent either.
                # Collapsing to element 0 would silently pick one point
                # (the GH #438 defect class).
                msg = (
                    f"Kinetic coefficient for field {fname!r} evaluated to an "
                    f"array rather than a scalar (expression: {kin_sym!r}). "
                    f"The generalized-eigenvalue builder writes M onto a "
                    f"per-mode mass matrix and has no representation for a "
                    f"varying M. Use plain solve_modal / auto-selection "
                    f"(position-dependent kinetics run on the convolution "
                    f"paths, GH #427) or a time-domain scheme "
                    f"(--scheme cvode / --scheme ida)."
                )
                raise KineticEvaluationError(msg)
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
        is_constraint = (
            eq.time_derivative_order == 0 and eq.field_name not in promoted_fields
        )

        if is_constraint:
            ci = c_idx_map[eq.field_name]
            for term_idx, term in enumerate(eq.rhs_terms):
                coeff = _resolve_constant_coeff(
                    term,
                    coeff_eval,
                    eq_idx=eq_idx,
                    term_idx=term_idx,
                )
                mult = multiplier_cache[term.operator]
                decomp = _OPERATOR_DECOMP[term.operator]
                t_order = decomp.time_order

                # Classification invariant (GH #457): a residual row can
                # carry NO time derivative of another residual field —
                # such an edge would have promoted both endpoints. The
                # old code folded d2_t(C) into S_cc as identity and
                # silently dropped v_C refs here; both are now impossible
                # by construction, and this assert keeps them so.
                res_base = term.field.removeprefix("v_")
                if res_base in c_idx_map and (
                    t_order >= 1 or term.field.startswith("v_")
                ):
                    msg = (
                        f"classification invariant violated: residual row "
                        f"{eq.field_name!r} carries a time reference "
                        f"{term.operator}({term.field}) to residual field "
                        f"{res_base!r} — second_order_sector should have "
                        f"promoted both (GH #457)"
                    )
                    raise AssertionError(msg)

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
                        (ci, complex(coeff), mult, t_order, fj),
                    )
                elif term.field in dyn_slot_map:
                    # Fallback: direct slot reference
                    dj = dyn_slot_map[term.field]
                    S_cd[:, ci, dj] += coeff * mult
            continue

        # Dynamical equation
        fi = dyn_field_idx[eq.field_name]
        orig_to_reduced[layout.field_slot_map[eq.field_name]]
        vel_slot = orig_to_reduced[layout.velocity_slot_map[eq.field_name]]

        for term_idx, term in enumerate(eq.rhs_terms):
            coeff = _resolve_constant_coeff(
                term,
                coeff_eval,
                eq_idx=eq_idx,
                term_idx=term_idx,
            )
            mult = multiplier_cache[term.operator]
            decomp = _OPERATOR_DECOMP[term.operator]
            t_order = decomp.time_order

            # Determine which field this term targets
            target_field = term.field
            # Strip v_ prefix to get base field name for velocity references
            is_vel_ref = target_field.startswith("v_")
            base_field = target_field[2:] if is_vel_ref else target_field

            if base_field in c_idx_map:
                # Reference to a RESIDUAL algebraic field: dispatch by the
                # TOTAL time order (operator time order + 1 for a v_ ref).
                # The old dispatch matched `target_field in c_idx_map`
                # first, so first_derivative_t(C) folded ḣ_c into
                # A_dc_field as h_c — GH #458 mechanism 2 (measured 0.98
                # on a_3, which has M = 1 and no other coupling).
                cj = c_idx_map[base_field]
                total_t = t_order + (1 if is_vel_ref else 0)
                if total_t == 0:
                    A_dc_field[:, vel_slot, cj] += coeff * mult
                elif total_t == 1:
                    A_dc_vel[:, vel_slot, cj] += coeff * mult
                else:
                    # Classification invariant (GH #457): an order-0
                    # field whose acceleration any row references is
                    # promoted by second_order_sector (algebraic recovery
                    # cannot supply q̈), so it cannot be in c_idx_map
                    # here. The OLD dispatch folded q̈_c as q_c
                    # (t_order-blind — GH #458).
                    msg = (
                        f"classification invariant violated: dynamical row "
                        f"{eq.field_name!r} carries {term.operator}"
                        f"({term.field}) (total time order {total_t}) of "
                        f"residual field {base_field!r} — "
                        f"second_order_sector should have promoted it "
                        f"(GH #457)"
                    )
                    raise AssertionError(msg)
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

    # ---- Second-order-sector composition: the per-mode pencil (GH #457) ----
    # The dynamical-sector rows form a first-order pencil  B·ẏ = A·y  over
    # the reduced slots, with ALL mass content on the B side: the LHS
    # kinetic on the velocity-row diagonal, RHS d2_t/accel cross terms
    # off-diagonal (moved to the LHS), and promoted rows (GH #457)
    # carrying no kinetic at all. K/D content sits RAW on the A side — the
    # row scale lives in B, which makes the residual-constraint couplings
    # (A_dc, folded in below) scale-correct by construction: the old
    # post-hoc field_correction/vel_coupling path added them onto
    # M⁻¹-scaled rows without the mass factor (GH #458, measured at
    # exactly (M−1)/M on the synthetic and a sign flip on κ⁻² rows).
    #
    # Composition: batched B⁻¹A when B is invertible on every mode (the
    # common case — mathematically identical to the old m_inv path);
    # otherwise the ordered-QZ deflating-subspace engine
    # (_pencil_deflate), which subsumes the old SVD mass-Schur branch,
    # its k-dependent bare-pinv gap (GH #460), and the velocity-aware
    # null-row cases uniformly, for ANY constraint index. Genuinely
    # singular pencils (true gauge) raise SingularPencilError.
    dets = np.linalg.det(M_mat)
    has_singular_M = bool(np.any(np.abs(dets) < 1e-12))
    use_engine = has_singular_M or bool(promoted_fields)
    if use_engine:
        logger.info(
            "Pencil composition: QZ deflating-subspace engine "
            "(singular mass: %s; promoted sector: %s)",
            has_singular_M,
            sorted(promoted_fields) if promoted_fields else "none",
        )

    # Jerk (d3_t) substitution requires an invertible mass matrix; apply
    # it at the field level so the pencil below stays first-order.
    if float(np.max(np.abs(J_mat))) > 1e-15:
        if use_engine:
            msg = (
                "Jerk (d3_t) terms combined with a singular mass matrix or "
                "a promoted second-order sector are not supported — the "
                "jerk substitution needs M⁻¹. GH #457."
            )
            raise NotImplementedError(msg)
        logger.info("Jerk substitution: applying d3_t elimination")
        m_inv = np.linalg.inv(M_mat)
        eff_k = m_inv @ K_mat
        eff_d = m_inv @ D_mat
        j_inv = m_inv @ J_mat
        fd_k = eff_d @ eff_k
        k_jerk = j_inv @ fd_k
        fd_d = eff_d @ eff_d
        d_jerk = j_inv @ (eff_k + fd_d)
        K_mat = M_mat @ (eff_k + k_jerk)
        D_mat = M_mat @ (eff_d + d_jerk)

    # ---- Assemble the pencil over reduced slots ----
    A0 = np.zeros((n_modes, n_dyn_slots, n_dyn_slots), dtype=np.complex128)
    B0 = np.zeros((n_modes, n_dyn_slots, n_dyn_slots), dtype=np.complex128)
    for fname in dyn_field_names:
        f_slot = orig_to_reduced[layout.field_slot_map[fname]]
        v_slot = orig_to_reduced[layout.velocity_slot_map[fname]]
        # Kinematic rows: dq/dt = v.
        B0[:, f_slot, f_slot] = 1.0
        A0[:, f_slot, v_slot] = 1.0
    for i, fname_i in enumerate(dyn_field_names):
        vel_i = orig_to_reduced[layout.velocity_slot_map[fname_i]]
        for j, fname_j in enumerate(dyn_field_names):
            field_j = orig_to_reduced[layout.field_slot_map[fname_j]]
            vel_j = orig_to_reduced[layout.velocity_slot_map[fname_j]]
            B0[:, vel_i, vel_j] += M_mat[:, i, j]
            A0[:, vel_i, field_j] += K_mat[:, i, j]
            A0[:, vel_i, vel_j] += D_mat[:, i, j]

    def _compose(
        A_p: NDArray[np.complex128],
        B_p: NDArray[np.complex128],
        what: str,
    ) -> tuple[NDArray[np.complex128], NDArray[np.complex128] | None]:
        """Per-mode generator (and manifold projector when reduced)."""
        if not use_engine:
            try:
                return (
                    np.asarray(np.linalg.solve(B_p, A_p), dtype=np.complex128),
                    None,
                )
            except np.linalg.LinAlgError:
                logger.info(
                    "Pencil composition (%s): batched solve hit a singular "
                    "mode — switching to the QZ engine",
                    what,
                )
        A_out = np.zeros_like(A_p)
        proj_out = np.zeros_like(A_p)
        for m in range(n_modes):
            A_out[m], proj_out[m] = _pencil_deflate(
                A_p[m],
                B_p[m],
                context=f"{what}, mode {m}",
                diagnostics=diagnostics,
                tag=(what, m),
            )
        return A_out, proj_out

    # ---- Deferred residual-row acceleration substitution ----
    # Residual algebraic rows may reference q̈ of dynamical-sector fields
    # (time_order>=2 operators). Substitute using the PRE-FOLD generator's
    # velocity rows: q̈_j = G0[vel_j, :]·y. Documented approximation (the
    # A_dc corrections folded below are not yet included — same semantics
    # as the old A_dd-based extraction; measured 8e-10 on E.cal).
    if deferred_constraint_terms:
        G0, _proj0 = _compose(A0, B0, "deferred-substitution")
        for ci, coeff_val, spatial_mult, t_order, fj in deferred_constraint_terms:
            if t_order == 2:
                vel_j = orig_to_reduced[layout.velocity_slot_map[dyn_field_names[fj]]]
                # S_cd[:, ci, :] += coeff·mult·(q̈_j row over all slots)
                S_cd[:, ci, :] += (
                    coeff_val * spatial_mult[:, np.newaxis] * G0[:, vel_j, :]
                )
            elif t_order >= 3:
                logger.warning(
                    "Constraint has time_order=%d operator — not yet handled",
                    t_order,
                )

    # ---- Residual-constraint Schur elimination ----
    if n_c > 0:
        # GH #459 (fixed here): the old |det| < 1e-14 gate was
        # scale-sensitive and its 1e-14 Tikhonov shift amplified
        # off-range roundoff by 1/ε into the recovery — measured as
        # 1e-4..1e-2 closure noise on dark_photon_plasma's rank-1
        # repeated constraint rows, and as a spurious ~2e5 abscissa once
        # the recovery fed the pencil fold. Rank-revealing pseudoinverse
        # instead: consistent (possibly redundant) systems are solved
        # exactly with the minimum-norm convention — the same
        # "undetermined modes frozen at zero" choice the mass machinery
        # has always documented — and no ε enters anywhere.
        sv_cc = np.linalg.svd(S_cc, compute_uv=False)
        sv_scale = np.maximum(sv_cc[:, 0], 1e-300)
        singular_mask = (sv_cc[:, -1] / sv_scale) < 1e-10
        S_cc_inv = np.asarray(
            np.linalg.pinv(S_cc, rcond=1e-10),
            dtype=np.complex128,
        )

        # See #290: expose S_cc_inv + singular_mask so the augmented
        # Pass 1 constraint recovery can use them.
        Scc_inv_out: NDArray[np.complex128] | None = np.asarray(
            S_cc_inv,
            dtype=np.complex128,
        )
        Scc_singular_mask_out: NDArray[np.bool_] | None = singular_mask

        # Recovery: c = -S_cc⁻¹ · S_cd · d.
        recovery = -(S_cc_inv @ S_cd)

        # ---- Fold residual couplings into the pencil (GH #458 fix) ----
        # identity/spatial(C) content enters the A side; v_C /
        # first_derivative_t(C) content is d/dt(recovery·y) = recovery·ẏ,
        # i.e. B-side content: B1 = B0 − A_dc_vel·recovery. Because B
        # carries the row mass, no separate M⁻¹ scaling exists to forget.
        a_fold = A0 + A_dc_field @ recovery
        b_fold = B0 - A_dc_vel @ recovery
    else:
        recovery = np.zeros((n_modes, 0, n_dyn_slots), dtype=np.complex128)
        a_fold = A0
        b_fold = B0
        Scc_inv_out = None
        Scc_singular_mask_out = None

    # ---- Final composition ----
    A_rhs, manifold_proj = _compose(a_fold, b_fold, "evolution")
    if diagnostics is not None:
        _record_pin_overlap(
            diagnostics,
            layout,
            orig_to_reduced,
            n_modes,
            n_dyn_slots,
            "evolution",
            per_mode=True,
        )

    logger.info(
        "Generalized evolution: %d residual constraint fields, %d dynamical "
        "slots (%d promoted), composition=%s",
        n_c,
        n_dyn_slots,
        len(promoted_fields),
        "qz-engine" if use_engine else "batched-solve",
    )

    # Constraint-velocity recovery: v_c = recovery · ẏ = recovery · A_rhs·y.
    v_recovery = recovery @ A_rhs if recovery.size > 0 else None

    return (
        A_rhs,
        None,  # B_lhs: composition is complete — no pencil is returned
        recovery,
        v_recovery,
        constraint_field_names,
        orig_to_reduced,
        Scc_inv_out,
        Scc_singular_mask_out,
        manifold_proj,
    )


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

    Fast path for systems without constraints and without time-derivative RHS
    operators — directly produces the pre-solved first-order evolution matrix
    without needing a generalized eigenvalue decomposition.

    Kinetic-coefficient (mass matrix) handling: when any dynamical equation
    carries ``kinetic_coefficient_symbolic`` (``M ẍ = K x`` form), ``M⁻¹`` is
    pre-computed per field via
    :func:`tidal.solver._kinetic.build_inverse_kinetic_diag` and folded into
    the velocity-row coefficients so the emitted matrix represents
    ``dv/dt = M⁻¹ · K(q)``. The companion generalized-eig path
    (:func:`_build_evolution_matrices`) reads the same field into ``M_mat``
    instead; the two paths must stay consistent (see module docstring).

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

    # #301 / #302: apply M⁻¹ to velocity-row entries so dv/dt = M⁻¹ K(q) for
    # theories with non-trivial kinetic_coefficient_symbolic. The generalized
    # eig path (_build_evolution_matrices) reads kinetic into M_mat directly;
    # this fast path instead folds M⁻¹ into the pre-solved evolution matrix.
    # build_inverse_kinetic_diag returns None when every dyn M ≈ 1 (fast path).
    # Shared `velocity_row_scale` helper enforces the cross-builder contract
    # (see modal.py module docstring; regression guard:
    # tests.test_solver_kinetic_consistency::TestAllModalPathsRespectKinetic).
    from tidal.solver._kinetic import (  # noqa: PLC0415
        build_inverse_kinetic_diag,
        velocity_row_scale,
    )

    m_inv = build_inverse_kinetic_diag(
        spec,
        coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
    )

    for _eq_idx, eq in enumerate(spec.equations):
        field_name = eq.field_name
        is_second_order = eq.time_derivative_order >= 2
        scale = velocity_row_scale(field_name, m_inv)

        if is_second_order:
            # Field slot and velocity slot
            field_slot = layout.field_slot_map[field_name]
            vel_slot = layout.velocity_slot_map[field_name]

            # dq/dt = v  →  A[field_slot, vel_slot] = 1
            A[:, field_slot, vel_slot] = 1.0

            # dv/dt = M⁻¹ · Σ coeff * operator(target_field)
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
                A[:, vel_slot, target_slot] += scale * coeff * mult

        else:
            # First-order: du/dt = Σ coeff * operator(target_field)
            # First-order kinetic scaling is not current-scope; leave unscaled
            # so M⁻¹ only affects 2nd-order EL velocity updates.
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
    coefficients or an ndarray for position-dependent ones. An ndarray
    reaching this function means a position-dependent coefficient reached a
    per-mode (constant-coefficient) builder — refuse loudly (GH #438).
    Before this guard the array was silently collapsed to
    ``resolved.ravel()[0]``, i.e. the coefficient's value at the FIRST grid
    point: the stability probe evaluated localized backgrounds at the
    domain edge (near-vacuum) on every gated sweep point / likelihood
    evaluation, and Pass-1 source matrices corner-collapsed O(ε)
    correction coefficients. NotImplementedError deliberately does not
    match the probe's ``(LinAlgError, ValueError)`` handler — same choice
    as the GH #421 kinetic guard — so the refusal propagates instead of
    being converted into a fabricated "tachyonic" verdict.
    """
    resolved = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
    if isinstance(resolved, np.ndarray):
        msg = (
            f"Position-dependent coefficient (field {term.field!r}, "
            f"operator {term.operator!r}) reached a per-mode builder that "
            "requires constant coefficients. Per-mode engines (genEig/"
            "Schur, per-mode Duhamel, the stability probe) have no per-k "
            "representation of a position-dependent coefficient — refusing "
            "instead of silently evaluating it at the first grid point "
            "(GH #438). Alternatives: plain solve_modal / auto-selection "
            "route position-dependent specs to the convolution builders "
            "(GH #427); time-domain schemes (--scheme cvode/ida) evaluate "
            "coefficients on the grid; a mode-coupling-aware stability "
            "gate is tracked in GH #441."
        )
        raise NotImplementedError(msg)
    return complex(resolved)


def _build_convolution_matrix(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
) -> NDArray[np.complex128]:
    """Build full evolution matrix for position-dependent, constraint-free systems.

    Position-dependent coefficients c(x) create convolution coupling in
    k-space: FFT[c(x)·u(x)] = ĉ * û (convolution).  This couples
    different k-modes, producing a full (n_total x n_total) matrix where
    n_total = n_slots x n_modes.

    For localized c(x) (e.g. Gaussian B₀), the convolution kernel ĉ(q)
    decays exponentially, making the matrix effectively banded.  The
    downstream ``_evolve_full_matrix`` exploits this by thresholding small
    entries and converting to sparse CSC format for faster expm_multiply.

    **Constraint handling:** This function handles dynamical (``time_order≥2``)
    and pure first-order (``time_order==1``) equations only. Systems with
    algebraic constraints (``time_order==0``) must go through
    :func:`_build_convolution_matrix_with_constraints`, which performs
    Schur elimination of constraint fields before returning the reduced
    matrix; the dispatch in :func:`solve_modal` routes accordingly. GH #379
    fixed a silent regression where pos-dep + constraint theories were
    routed here and the constraint rows were treated as first-order ODEs,
    producing divergent dynamics at any nonzero coupling.

    Reference: Burns et al. (2020), Phys. Rev. Research 2:023068.
    """
    if any(eq.time_derivative_order == 0 for eq in spec.equations):
        msg = (
            "_build_convolution_matrix called on a spec with algebraic "
            "constraints (time_derivative_order==0). The dispatch in "
            "solve_modal should route this case to "
            "_build_convolution_matrix_with_constraints (GH #379)."
        )
        raise AssertionError(msg)
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

    # #301 / #302 / GH #367 fix: apply M⁻¹ to velocity-row contributions so
    # `dv/dt = M⁻¹·K(q)` for theories with non-trivial `kinetic_coefficient_symbolic`
    # (e.g. `-1/kappa²` on graviton h-modes in Gertsenshtein). The per-mode path
    # `_build_per_mode_matrices` already does this; the convolution path previously
    # did not, producing wrong-sign EOM that drove `expm_multiply` to ~10⁹× error
    # — diagnosed 2026-05-19 as the actual root cause of GH #367 (originally
    # believed to be a non-normal Fourier convolution artifact). Shared
    # `velocity_row_scale` helper enforces the cross-builder contract (see
    # modal.py module docstring; regression guard:
    # tests.test_solver_kinetic_consistency::TestAllModalPathsRespectKinetic).
    from tidal.solver._kinetic import (  # noqa: PLC0415
        build_inverse_kinetic_diag,
        velocity_row_scale,
    )

    # GH #427: pass `grid` so a position-dependent kinetic evaluates to a
    # per-grid-point M⁻¹(x) profile (GH #382 machinery) instead of raising.
    m_inv = build_inverse_kinetic_diag(
        spec,
        coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        grid,
    )

    for _eq_idx, eq in enumerate(spec.equations):
        field_name = eq.field_name
        is_second_order = eq.time_derivative_order >= 2
        scale = velocity_row_scale(field_name, m_inv)

        if is_second_order:
            field_slot = layout.field_slot_map[field_name]
            vel_slot = layout.velocity_slot_map[field_name]

            # dq/dt = v → diagonal identity coupling between field and velocity
            for m in range(n_modes):
                row = field_slot * n_modes + m
                col = vel_slot * n_modes + m
                A[row, col] = 1.0

            # dv/dt = M⁻¹ · Σ coeff(x) * operator(target_field)
            for _term_idx, term in enumerate(eq.rhs_terms):
                target_slot = _resolve_target_slot(term.field)
                if target_slot is None:
                    continue
                mult = multiplier_cache[term.operator]

                if not term.position_dependent and not isinstance(scale, np.ndarray):
                    # Constant coefficient AND constant kinetic: diagonal in
                    # mode space. An ndarray scale means M⁻¹(x) makes even a
                    # constant RHS coefficient effectively position-dependent
                    # (GH #427) — routed to the convolution branch below.
                    coeff = _resolve_constant_coeff(
                        term,
                        coeff_eval,
                        eq_idx=_eq_idx,
                        term_idx=_term_idx,
                    )
                    for m in range(n_modes):
                        row = vel_slot * n_modes + m
                        col = target_slot * n_modes + m
                        A[row, col] += scale * coeff * mult[m]
                else:
                    # Position-dependent coefficient and/or kinetic:
                    # convolution coupling
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
                        scale=scale,
                    )
        else:
            # First-order ODE (time_derivative_order == 1): direct evolution,
            # dq/dt = Σ coeff·op(target). M⁻¹ scaling is a 2nd-order concept
            # and does not apply here. Algebraic constraints (time_order==0)
            # are rejected above; this branch is now reachable only by genuine
            # 1st-order fields (Bundled fix #2 of GH #379, locks the invariant).
            assert eq.time_derivative_order == 1, (
                f"Unexpected time_derivative_order={eq.time_derivative_order} "
                f"for field {field_name!r} in convolution path"
            )
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
    scale: complex | NDArray[np.float64] = 1.0,
) -> None:
    """Add convolution coupling from a position-dependent coefficient.

    The product c(x)·op(u(x)) in k-space becomes a convolution:
    FFT[c·op(u)]_k = Σ_k' ĉ(k-k') · mult(k') · û(k')

    This creates off-diagonal entries in the evolution matrix coupling
    different k-modes.

    ``scale`` multiplies the convolution contribution — used by the caller
    to apply the inverse-kinetic factor `M⁻¹` for 2nd-order equations with
    `kinetic_coefficient_symbolic != 1` (Gertsenshtein h-modes, etc.). See
    `_build_convolution_matrix` (modal.py) for the rationale (GH #367 fix).
    A scalar ``scale`` multiplies the k-space block; an ndarray ``scale``
    is a REAL-SPACE ``M⁻¹(x)`` profile (GH #427) and is folded into the
    physical-space coefficient before the FFT — it must never multiply a
    k-space block (wrong semantics, and a shape mismatch: n_points vs
    n_modes).

    GH #384 Phase A′: for BSM-separable terms the block of size
    (n_modes × n_modes) is independent of the BSM scalar and reusable
    across PolyChord likelihood calls. The cache layer in
    ``tidal/solver/_conv_block_cache.py`` memoizes it.
    """
    row_start = row_slot * n_modes
    col_start = col_slot * n_modes
    if isinstance(scale, np.ndarray):
        block = _conv_block_with_kinetic(
            term,
            operator_mult,
            scale,
            coeff_eval,
            grid,
            rfft_shape,
            n_modes,
            eq_idx=eq_idx,
            term_idx=term_idx,
        )
        A[row_start : row_start + n_modes, col_start : col_start + n_modes] += block
        return
    block = _compute_conv_block_cached_or_fresh(
        term,
        operator_mult,
        coeff_eval,
        grid,
        rfft_shape,
        n_modes,
        eq_idx=eq_idx,
        term_idx=term_idx,
    )
    # Accumulate: A[row_block, col_block] += scale * block
    A[row_start : row_start + n_modes, col_start : col_start + n_modes] += scale * block


def _compute_conv_block_uncached(
    coeff_array: NDArray[np.float64] | float,
    operator_mult: NDArray[np.complex128],
    grid_shape: tuple[int, ...],
    spectral_shape: tuple[int, ...],
    n_modes: int,
) -> NDArray[np.complex128]:
    """Build the (n_modes × n_modes) convolution-matrix block at scale=1.

    Pure function of ``(coeff_array, operator_mult, grid_shape)``. Used by
    both ``_add_convolution_coupling`` and ``_term_conv_block`` (the latter
    nested inside ``_build_convolution_matrix_with_constraints``).

    For constant coefficients, returns a diagonal block ``diag(c · mult)``.
    For position-dependent coefficients (GH #445), the block is the exact
    circulant of pointwise multiplication in the FULL fftn basis::

        block[k, k'] = ĉ[(k − k') mod N] / N_tot · mult[k']

    with ``ĉ = fftn(c)``, representing ``fftn(c(x) · op(y))`` for
    ``op(y) = ifftn(mult · ŷ)``. This replaces the former probe-vector
    loop over the rfft HALF-spectrum, which was only R-linear-correct:
    the DC/Nyquist bins carry no imaginary degree of freedom there, so
    sin-phase content folding onto them was misrepresented (7–25% per-mode
    action errors — see GH #445 for the repro). The circulant needs one
    FFT of ``c`` instead of ``n_modes`` FFT pairs, and is exact for every
    mode. ``spectral_shape`` must equal ``grid_shape`` on this path (the
    pos-dep dispatch passes the full-basis shape).

    GH #384 Phase A′: this function is the cache target. Its output for
    a BSM-separable term is identical across PolyChord likelihood calls.
    """
    if isinstance(coeff_array, (int, float)):
        block = np.zeros((n_modes, n_modes), dtype=np.complex128)
        scaled = float(coeff_array)
        for m in range(n_modes):
            block[m, m] = scaled * operator_mult[m]
        return block

    if spectral_shape != grid_shape:
        msg = (
            "position-dependent convolution blocks require the full fftn "
            f"basis (spectral_shape {spectral_shape} != grid_shape "
            f"{grid_shape}); the rfft half-spectrum basis misrepresents "
            "DC/Nyquist-fold coupling (GH #445)"
        )
        raise ValueError(msg)

    n_tot = int(np.prod(grid_shape))
    chat = np.fft.fftn(np.asarray(coeff_array, dtype=np.float64).reshape(grid_shape))
    # Per-axis index differences (k − k') mod N_axis, raveled over modes.
    axis_indices = np.unravel_index(np.arange(n_tot), grid_shape)
    diff = tuple(
        (idx[:, None] - idx[None, :]) % grid_shape[a]
        for a, idx in enumerate(axis_indices)
    )
    block = (chat[diff] / float(n_tot)).astype(np.complex128)
    block *= operator_mult[None, :]
    return block


def _conv_block_with_kinetic(
    term: OperatorTerm,
    operator_mult: NDArray[np.complex128],
    m_inv_profile: NDArray[np.float64],
    coeff_eval: CoefficientEvaluator,
    grid: GridInfo,
    rfft_shape: tuple[int, ...],
    n_modes: int,
    *,
    eq_idx: int = -1,
    term_idx: int = -1,
) -> NDArray[np.complex128]:
    """Build a term's convolution block with a position-dependent ``M⁻¹(x)`` folded in.

    GH #427: a position-dependent kinetic coefficient makes the velocity-row
    scale a real-space profile. Dividing pointwise in real space BEFORE the
    forward FFT is mathematically identical to applying the ``M̂⁻¹(k−k′)``
    convolution in k-space, and the resulting effective coefficient
    ``M⁻¹(x)·c(x)`` is exactly what the existing ĉ(k−k′) probe-vector
    machinery handles. This applies equally when the RHS coefficient is
    constant — ``M⁻¹(x)·c`` is position-dependent regardless, which is why
    the constant-coefficient diagonal branches must route here too.

    Bypasses the GH #384 BSM block cache: the cached block is keyed on the
    RHS geometry alone and would not include the kinetic profile. Extending
    the cache key with a kinetic identity is a possible future optimization;
    the inference-relevant canonicalized flows have constant-kinetic base
    specs and never reach this path, so nothing hot is de-cached.
    """
    coeff = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
    profile = np.asarray(m_inv_profile, dtype=np.float64).reshape(grid.shape)
    # `coeff` is a float for constant terms and a grid-shaped array for
    # position-dependent ones; either way the product is the grid-shaped
    # effective coefficient.
    coeff_eff = profile * coeff
    return _compute_conv_block_uncached(
        coeff_eff, operator_mult, grid.shape, rfft_shape, n_modes
    )


def _compute_conv_block_cached_or_fresh(
    term: OperatorTerm,
    operator_mult: NDArray[np.complex128],
    coeff_eval: CoefficientEvaluator,
    grid: GridInfo,
    rfft_shape: tuple[int, ...],
    n_modes: int,
    *,
    eq_idx: int = -1,
    term_idx: int = -1,
) -> NDArray[np.complex128]:
    """Compute the convolution-matrix block for ``term``, using cache when possible.

    GH #384 Phase A′ entry point. For BSM-separable terms (identified by
    ``coeff_eval.bsm_separable_factors``), the block is computed once per
    (geometry, term) and stored in the module-level
    ``_conv_block_cache``; subsequent PolyChord likelihood calls multiply
    by the new BSM scalar. For non-separable or non-BSM-tagged terms,
    every call computes fresh (same as pre-#384).
    """
    from tidal.solver._conv_block_cache import (  # noqa: PLC0415
        get_or_compute,
        make_geometry_hash,
        make_key,
    )
    from tidal.symbolic._eval_utils import evaluate_coefficient  # noqa: PLC0415

    factors = coeff_eval.bsm_separable_factors(eq_idx, term_idx)
    if factors is None:
        # No BSM separability info → compute fresh (no cache benefit, but
        # correctness preserved). Path taken for non-separable terms and
        # non-PolyChord callers that don't tag sampled params.
        coeff_array = coeff_eval.resolve(term, 0.0, eq_idx=eq_idx, term_idx=term_idx)
        return _compute_conv_block_uncached(
            coeff_array, operator_mult, grid.shape, rfft_shape, n_modes
        )

    bsm_str, geom_str = factors
    # Build the cache key. geometry_hash depends only on non-BSM params,
    # so PolyChord calls at varying BSM sample the same cached block.
    sampled: tuple[str, ...] = coeff_eval._spec.metadata.get(  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        "_inference_sampled_params", ()
    )
    bsm_set = set(sampled)
    geom_hash = make_geometry_hash(
        coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        bsm_set,
    )
    key = make_key(geom_str, term.operator, tuple(grid.shape), geom_hash)

    def _compute() -> NDArray[np.complex128]:
        # Evaluate the BSM-stripped geometric expression on the grid. The
        # result is an N-element array (1D) for position-dependent terms.
        coeff_array = evaluate_coefficient(
            geom_str,
            coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
            coeff_eval._coordinates,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
            coeff_eval._coord_arrays,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
            0.0,
        )
        return _compute_conv_block_uncached(
            coeff_array, operator_mult, grid.shape, rfft_shape, n_modes
        )

    block_unit = get_or_compute(key, _compute)

    # Multiply by current BSM scalar to recover the full block.
    if bsm_str and bsm_str != "1":
        bsm_val = evaluate_coefficient(
            bsm_str,
            coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
            coeff_eval._coordinates,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
            None,  # bsm_str has no coord deps by construction
            0.0,
        )
        return complex(bsm_val) * block_unit
    return block_unit


# ---------------------------------------------------------------------------
# Position-dependent + constraint convolution builder (GH #379)
# ---------------------------------------------------------------------------


def _build_convolution_matrix_with_constraints(
    spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
    diagnostics: PencilDiagnostics | None = None,
) -> tuple[
    NDArray[np.complex128],  # A_reduced (n_dyn_tot x n_dyn_tot)
    NDArray[np.complex128],  # recovery (n_c_tot x n_dyn_tot)
    list[str],  # RESIDUAL constraint field names
    dict[int, int],  # orig_to_reduced slot mapping
    NDArray[np.complex128] | None,  # manifold projector (n_dyn_tot²) or None
]:
    """Build reduced evolution matrix for pos-dep theories with constraints.

    GH #457: promoted order-0 rows (``spec.second_order_sector``) join
    the dynamical sector as pencil rows with mass content on the B side;
    the composition runs the deflating-subspace engine when B is
    singular, returning the manifold projector for on-manifold ICs.
    Only RESIDUAL rows are Schur-eliminated.

    This is the position-dependent analog of :func:`_build_evolution_matrices`
    (the constant-coefficient generalized-eigenvalue path). It builds
    convolution sub-matrices for the four blocks of the first-order system:

    * ``A_dd``: dynamical → dynamical (size ``n_dyn_slots·n_modes``)
    * ``A_dc_field``, ``A_dc_vel``: dynamical → constraint (field / velocity)
    * ``K_cd``: constraint → dynamical (field/velocity, plus deferred ẍ)
    * ``K_cc``: constraint → constraint

    Then Schur-eliminates the constraint fields. The recovery matrix
    ``recovery = -K_cc⁻¹ · K_cd`` gives ``q_c(t) = recovery · y_d(t)``;
    the reduced matrix is ``A_reduced = A_dd + A_dc_field · recovery``.

    **JSON convention:** the constraint equation `q_c = RHS` is encoded
    with the LHS already absorbed (the RHS terms include any self-coupling
    needed to make the canonical form ``K_cc·q_c + K_cd·q_d = 0``); see
    the synthetic test in ``tests/test_solver_constraint_posdep.py``.

    **dyn-RHS time-derivative-on-constraint handling**: a dynamical RHS
    term referencing ``v_constraint``, ``first_derivative_t(constraint)``,
    or ``d2_t(constraint)`` requires substituting derivatives of q_c.
    These are folded in via composition with ``recovery`` and a
    ``B_lhs`` pre-solve when needed (mirrors genEig's vel_coupling
    machinery at modal.py:1304-1322).

    Returns
    -------
    A_reduced : ndarray
        First-order evolution matrix on dynamical slots only,
        shape ``(n_dyn_slots·n_modes, n_dyn_slots·n_modes)``.
    recovery : ndarray
        Constraint-recovery matrix; ``q_c = recovery @ y_d``.
        Shape ``(n_c·n_modes, n_dyn_slots·n_modes)``.
    constraint_field_names : list[str]
        Names of algebraic constraint fields (input order).
    orig_to_reduced : dict[int, int]
        Mapping from original slot indices to reduced (dyn-only) indices.

    Notes
    -----
    GH #379: this function did not exist; pos-dep + constraint theories
    were silently routed to :func:`_build_convolution_matrix` which had
    no constraint handling, producing divergent dynamics at any nonzero
    coupling. The bug manifested as ~10²⁵² blowup with a misleading
    "GH #367 small-grid expm_multiply" error message.
    """
    import logging  # noqa: PLC0415
    import warnings  # noqa: PLC0415

    from tidal.solver._kinetic import (  # noqa: PLC0415
        build_inverse_kinetic_diag,
        velocity_row_scale,
    )

    logger = logging.getLogger(__name__)

    n_modes = int(np.prod(rfft_shape))

    # ---- Classify equations ----
    # GH #457: promoted order-0 rows (second_order_sector) join the
    # dynamical sector; only RESIDUAL rows are Schur-eliminated.
    promoted_fields = spec.second_order_sector.promoted
    dyn_field_names: list[str] = [
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order >= 2 or eq.field_name in promoted_fields
    ]
    dyn_field_idx: dict[str, int] = {n: i for i, n in enumerate(dyn_field_names)}
    constraint_field_names: list[str] = [
        eq.field_name
        for eq in spec.equations
        if eq.time_derivative_order == 0 and eq.field_name not in promoted_fields
    ]
    c_idx_map: dict[str, int] = {n: i for i, n in enumerate(constraint_field_names)}
    n_c = len(constraint_field_names)

    # GH #468 (constraints-as-slots): when a promoted second-order sector
    # is present, the residual constraints are NOT pre-eliminated —
    # Schur recovery on the localized residual K_cc injects ‖·‖ ~ 1/B₀⁴
    # entries into the pencil, collapsing the deflation's tolerance
    # headroom (measured: genuine smin ~6e-12 on E.cal's 528-dim
    # pencil). Instead the residual constraint fields become pencil
    # SLOTS and their algebraic rows pure-A pencil rows: every entry
    # stays at natural equation scale, the deferred d2_t(dyn)
    # substitution becomes EXACT B-side content, and the deflation
    # handles elimination + promotion + IC projection in one object.
    # Non-promoted specs keep the classic Schur route bit-for-bit.
    constraints_as_slots = bool(promoted_fields)

    # ---- Reduced slot mapping ----
    orig_to_reduced: dict[int, int] = {}
    for si, slot in enumerate(layout.slots):
        if slot.kind == "constraint" and not constraints_as_slots:
            continue
        orig_to_reduced[si] = len(orig_to_reduced)
    n_dyn_slots = len(orig_to_reduced)

    # ---- M⁻¹ for kinetic-coefficient handling on dynamical rows ----
    # GH #427: pass `grid` so a position-dependent kinetic evaluates to a
    # per-grid-point M⁻¹(x) profile (GH #382 machinery) instead of raising.
    m_inv = build_inverse_kinetic_diag(
        spec,
        coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        grid,
    )

    # ---- Multiplier cache ----
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op not in multiplier_cache:
                mult_val = _EXACT_MULTIPLIERS[op](k_grid)
                mult_full = np.broadcast_to(mult_val, rfft_shape)
                multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    # ---- Helpers --------------------------------------------------------
    def _emit_term(
        out_matrix: NDArray[np.complex128],
        row_block: int,
        col_block: int,
        term: OperatorTerm,
        mult: NDArray[np.complex128],
        scale: complex | NDArray[np.float64] = 1.0,
        eq_idx: int = -1,
        term_idx: int = -1,
    ) -> None:
        """Emit a term's contribution to a (row_block, col_block) of out_matrix."""
        if not term.position_dependent and not isinstance(scale, np.ndarray):
            coeff = _resolve_constant_coeff(
                term,
                coeff_eval,
                eq_idx=eq_idx,
                term_idx=term_idx,
            )
            scaled = scale * coeff
            for m in range(n_modes):
                out_matrix[
                    row_block * n_modes + m,
                    col_block * n_modes + m,
                ] += scaled * mult[m]
        else:
            # Position-dependent coefficient and/or M⁻¹(x) kinetic (GH
            # #427) — _add_convolution_coupling folds an ndarray scale
            # into the real-space coefficient before the FFT.
            _add_convolution_coupling(
                out_matrix,
                row_block,
                col_block,
                term,
                coeff_eval,
                mult,
                grid,
                rfft_shape,
                n_modes,
                eq_idx=eq_idx,
                term_idx=term_idx,
                scale=scale,
            )

    def _term_conv_block(
        term: OperatorTerm,
        mult: NDArray[np.complex128],
        scale: complex | NDArray[np.float64] = 1.0,
        eq_idx: int = -1,
        term_idx: int = -1,
    ) -> NDArray[np.complex128]:
        """Build the (n_modes × n_modes) convolution block for a single term.

        Delegates to the shared cache-aware helper
        ``_compute_conv_block_cached_or_fresh`` (GH #384 Phase A′).
        For BSM-separable terms the block lives in
        ``tidal.solver._conv_block_cache`` keyed on the geometry and is
        reused across PolyChord likelihood calls.
        """
        if isinstance(scale, np.ndarray):
            # M⁻¹(x) kinetic (GH #427): fold the real-space profile into
            # the coefficient before the FFT — never scale a k-space block
            # by a real-space array. Covers constant RHS coefficients too
            # (M⁻¹(x)·c is position-dependent regardless).
            return _conv_block_with_kinetic(
                term,
                mult,
                scale,
                coeff_eval,
                grid,
                rfft_shape,
                n_modes,
                eq_idx=eq_idx,
                term_idx=term_idx,
            )
        if not term.position_dependent:
            block = np.zeros((n_modes, n_modes), dtype=np.complex128)
            coeff = _resolve_constant_coeff(
                term,
                coeff_eval,
                eq_idx=eq_idx,
                term_idx=term_idx,
            )
            scaled = scale * coeff
            for m in range(n_modes):
                block[m, m] = scaled * mult[m]
            return block
        block_unit = _compute_conv_block_cached_or_fresh(
            term,
            mult,
            coeff_eval,
            grid,
            rfft_shape,
            n_modes,
            eq_idx=eq_idx,
            term_idx=term_idx,
        )
        return scale * block_unit

    # ---- Allocate matrices ----------------------------------------------
    n_dyn_tot = n_dyn_slots * n_modes
    # GH #468: in constraints-as-slots mode the residual Schur machinery
    # (K_cc/K_cd/A_dc/recovery) is bypassed — allocate it empty so every
    # downstream matmul no-ops naturally.
    n_c_tot = (0 if constraints_as_slots else n_c) * n_modes
    A_dd = np.zeros((n_dyn_tot, n_dyn_tot), dtype=np.complex128)
    # GH #457: B-side mass content of the pencil B·ẏ = A·y. Collected as
    # POSITIVE contributions of RHS q̈-target terms (and promoted-row
    # mass); composed below as B = I_adj − B_mass − vel_coupling, where
    # I_adj has 0 on promoted velocity-row diagonals (their LHS is
    # algebraic — no v̇_self term).
    B_mass = np.zeros((n_dyn_tot, n_dyn_tot), dtype=np.complex128)
    has_b_mass = bool(promoted_fields)
    A_dc_field = np.zeros((n_dyn_tot, n_c_tot), dtype=np.complex128)
    A_dc_vel = np.zeros((n_dyn_tot, n_c_tot), dtype=np.complex128)
    K_cd = np.zeros((n_c_tot, n_dyn_tot), dtype=np.complex128)
    K_cc = np.zeros((n_c_tot, n_c_tot), dtype=np.complex128)

    # deferred_terms_constraint: constraint RHS references d2_t/first_dt(v) of dyn → ẍ_fj
    # entries: (ci, eq_idx, term_idx, term, mult, fj_dyn_idx)
    deferred_terms_constraint: list[
        tuple[int, int, int, OperatorTerm, NDArray[np.complex128], int]
    ] = []
    # deferred_terms_dyn_velc: dyn RHS references v_constraint with t_order=0 OR
    # first_derivative_t(constraint) — both require v_c substitution after recovery.
    # entries: (vel_red, eq_idx, term_idx, term, mult, cj_constraint_idx, scale)
    # GH #444: `scale` is the emitting equation's velocity-row M⁻¹
    # (velocity_row_scale). The deferral previously dropped it and the
    # application sites used scale=1.0, so every deferred contribution on a
    # row with a non-unit kinetic (−κ⁻², −ξ, …: all 8 dual-Gaussian roster
    # specs) was off by a factor of M — a kinetic-contract violation
    # (module docstring; TestAllModalPathsRespectKinetic).
    deferred_terms_dyn_velc: list[
        tuple[
            int,
            int,
            int,
            OperatorTerm,
            NDArray[np.complex128],
            int,
            complex | NDArray[np.float64],
        ]
    ] = []

    # ---- A_dd kinematic rows: dq/dt = v ---------------------------------
    for fname in dyn_field_names:
        fs_red = orig_to_reduced[layout.field_slot_map[fname]]
        vs_red = orig_to_reduced[layout.velocity_slot_map[fname]]
        for m in range(n_modes):
            A_dd[fs_red * n_modes + m, vs_red * n_modes + m] = 1.0

    # ---- Iterate equations ---------------------------------------------
    warned_unhandled: set[tuple[str, str, str]] = set()

    for eq_idx, eq in enumerate(spec.equations):
        is_promoted_row = eq.field_name in promoted_fields
        if eq.time_derivative_order >= 2 or is_promoted_row:
            # ---------------- Dynamical-sector equation ----------------
            # Promoted rows (GH #457) have an algebraic LHS: no kinetic
            # to divide by (scale = 1) and no v̇_self — their mass
            # content is the off-diagonal q̈ cross terms collected into
            # B_mass below.
            scale = 1.0 if is_promoted_row else velocity_row_scale(eq.field_name, m_inv)
            vel_red = orig_to_reduced[layout.velocity_slot_map[eq.field_name]]
            for term_idx, term in enumerate(eq.rhs_terms):
                target_field = term.field
                is_vel_ref = target_field.startswith("v_")
                base_field = target_field[2:] if is_vel_ref else target_field
                decomp = _OPERATOR_DECOMP[term.operator]
                mult = multiplier_cache[term.operator]
                t_order = decomp.time_order

                if base_field in c_idx_map:
                    # ----- target is a RESIDUAL algebraic field -----
                    total_t = t_order + (1 if is_vel_ref else 0)
                    if constraints_as_slots:
                        # GH #468: the residual field is a pencil SLOT.
                        c_col = orig_to_reduced[layout.field_slot_map[base_field]]
                        if total_t == 0:
                            _emit_term(
                                A_dd,
                                vel_red,
                                c_col,
                                term,
                                mult,
                                scale=scale,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                        elif total_t == 1:
                            # c·Ċ — exact B-side content (replaces the
                            # velc deferral for this mode).
                            _emit_term(
                                B_mass,
                                vel_red,
                                c_col,
                                term,
                                mult,
                                scale=scale,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                            has_b_mass = True
                        else:
                            msg = (
                                f"classification invariant violated: row "
                                f"{eq.field_name!r} carries {term.operator}"
                                f"({term.field}) of residual field "
                                f"{base_field!r} (GH #457)"
                            )
                            raise AssertionError(msg)
                        continue
                    cj = c_idx_map[base_field]
                    if total_t >= 2:
                        # Classification invariant (GH #457): an order-0
                        # field whose acceleration any row references is
                        # PROMOTED — it cannot be residual here.
                        msg = (
                            f"classification invariant violated: row "
                            f"{eq.field_name!r} carries {term.operator}"
                            f"({term.field}) of residual field "
                            f"{base_field!r} (GH #457)"
                        )
                        raise AssertionError(msg)
                    if total_t == 1:
                        # coeff · v_c (or first_dt(q_c)) — needs
                        # v_c = d/dt(recovery·y_d); defer until recovery.
                        deferred_terms_dyn_velc.append(
                            (vel_red, eq_idx, term_idx, term, mult, cj, scale),
                        )
                    else:
                        # coeff · op · q_c → A_dc_field
                        _emit_term(
                            A_dc_field,
                            vel_red,
                            cj,
                            term,
                            mult,
                            scale=scale,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                elif base_field in dyn_field_idx:
                    # ----- target is a DYNAMICAL field -----
                    if is_vel_ref:
                        if t_order == 0:
                            target_red = orig_to_reduced[
                                layout.velocity_slot_map[base_field]
                            ]
                            _emit_term(
                                A_dd,
                                vel_red,
                                target_red,
                                term,
                                mult,
                                scale=scale,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                        else:
                            # first_dt(v_X) = ẍ_X: mass-side content of
                            # the pencil (GH #457) — collect into B_mass.
                            target_red = orig_to_reduced[
                                layout.velocity_slot_map[base_field]
                            ]
                            _emit_term(
                                B_mass,
                                vel_red,
                                target_red,
                                term,
                                mult,
                                scale=scale,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                            has_b_mass = True
                    elif t_order == 0:
                        target_red = orig_to_reduced[layout.field_slot_map[base_field]]
                        _emit_term(
                            A_dd,
                            vel_red,
                            target_red,
                            term,
                            mult,
                            scale=scale,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                    elif t_order == 1:
                        # ẋ_dyn = v_dyn — go to velocity slot
                        target_red = orig_to_reduced[
                            layout.velocity_slot_map[base_field]
                        ]
                        _emit_term(
                            A_dd,
                            vel_red,
                            target_red,
                            term,
                            mult,
                            scale=scale,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                    else:
                        # d2_t(dyn_field): q̈ cross coupling = mass-side
                        # content of the pencil B·ẏ = A·y (GH #457) —
                        # this includes the promoted rows' M_cc terms.
                        target_red = orig_to_reduced[
                            layout.velocity_slot_map[base_field]
                        ]
                        _emit_term(
                            B_mass,
                            vel_red,
                            target_red,
                            term,
                            mult,
                            scale=scale,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                        has_b_mass = True
            continue

        if eq.time_derivative_order == 0:
            # ---------------- Algebraic constraint equation ----------------
            if constraints_as_slots:
                # GH #468: the residual row is a pure-A pencil row at the
                # constraint's own slot (zero B-row — algebraic LHS).
                # d2_t(dyn)/first_dt(v_dyn) content is EXACT B-side mass
                # (q̈ = v̇), replacing the old A_dd-based deferred
                # substitution approximation.
                c_row = orig_to_reduced[layout.field_slot_map[eq.field_name]]
                for term_idx, term in enumerate(eq.rhs_terms):
                    target_field = term.field
                    is_vel_ref = target_field.startswith("v_")
                    base_field = target_field[2:] if is_vel_ref else target_field
                    decomp = _OPERATOR_DECOMP[term.operator]
                    mult = multiplier_cache[term.operator]
                    t_order = decomp.time_order
                    total_t = t_order + (1 if is_vel_ref else 0)

                    if base_field in c_idx_map:
                        if total_t > 0:
                            msg = (
                                f"classification invariant violated: "
                                f"residual row {eq.field_name!r} carries "
                                f"{term.operator}({term.field}) of residual "
                                f"field {base_field!r} (GH #457)"
                            )
                            raise AssertionError(msg)
                        tgt = orig_to_reduced[layout.field_slot_map[base_field]]
                        _emit_term(
                            A_dd,
                            c_row,
                            tgt,
                            term,
                            mult,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                    elif base_field in dyn_field_idx:
                        if total_t == 0:
                            tgt = orig_to_reduced[layout.field_slot_map[base_field]]
                            _emit_term(
                                A_dd,
                                c_row,
                                tgt,
                                term,
                                mult,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                        elif total_t == 1:
                            tgt = orig_to_reduced[layout.velocity_slot_map[base_field]]
                            _emit_term(
                                A_dd,
                                c_row,
                                tgt,
                                term,
                                mult,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                        elif total_t == 2:
                            tgt = orig_to_reduced[layout.velocity_slot_map[base_field]]
                            _emit_term(
                                B_mass,
                                c_row,
                                tgt,
                                term,
                                mult,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                            has_b_mass = True
                        else:
                            msg = (
                                f"residual row {eq.field_name!r}: "
                                f"{term.operator}({term.field}) is "
                                f"jerk-level content — not supported "
                                f"(GH #469)."
                            )
                            raise NotImplementedError(msg)
                    else:
                        warned_unhandled.add(
                            (eq.field_name, term.operator, target_field),
                        )
                continue

            ci = c_idx_map[eq.field_name]
            for term_idx, term in enumerate(eq.rhs_terms):
                target_field = term.field
                is_vel_ref = target_field.startswith("v_")
                base_field = target_field[2:] if is_vel_ref else target_field
                decomp = _OPERATOR_DECOMP[term.operator]
                mult = multiplier_cache[term.operator]
                t_order = decomp.time_order

                if base_field in c_idx_map:
                    # residual → residual (K_cc). Time content here is
                    # impossible: any inter-constraint time reference
                    # promotes BOTH endpoints (GH #457) — the old code
                    # silently dropped these (the h_3-class defect).
                    cj = c_idx_map[base_field]
                    if is_vel_ref or t_order > 0:
                        msg = (
                            f"classification invariant violated: residual "
                            f"row {eq.field_name!r} carries "
                            f"{term.operator}({term.field}) of residual "
                            f"field {base_field!r} (GH #457)"
                        )
                        raise AssertionError(msg)
                    _emit_term(
                        K_cc,
                        ci,
                        cj,
                        term,
                        mult,
                        eq_idx=eq_idx,
                        term_idx=term_idx,
                    )
                elif base_field in dyn_field_idx:
                    fj = dyn_field_idx[base_field]
                    if is_vel_ref:
                        if t_order == 0:
                            # v_dyn → velocity slot
                            target_red = orig_to_reduced[
                                layout.velocity_slot_map[base_field]
                            ]
                            _emit_term(
                                K_cd,
                                ci,
                                target_red,
                                term,
                                mult,
                                eq_idx=eq_idx,
                                term_idx=term_idx,
                            )
                        elif t_order == 1:
                            # first_dt(v_dyn) = ẍ_dyn → defer
                            deferred_terms_constraint.append(
                                (ci, eq_idx, term_idx, term, mult, fj),
                            )
                        else:
                            warned_unhandled.add(
                                (eq.field_name, term.operator, target_field),
                            )
                    elif t_order == 0:
                        target_red = orig_to_reduced[layout.field_slot_map[base_field]]
                        _emit_term(
                            K_cd,
                            ci,
                            target_red,
                            term,
                            mult,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                    elif t_order == 1:
                        target_red = orig_to_reduced[
                            layout.velocity_slot_map[base_field]
                        ]
                        _emit_term(
                            K_cd,
                            ci,
                            target_red,
                            term,
                            mult,
                            eq_idx=eq_idx,
                            term_idx=term_idx,
                        )
                    elif t_order == 2:
                        # d2_t(dyn_field) → ẍ_dyn → defer
                        deferred_terms_constraint.append(
                            (ci, eq_idx, term_idx, term, mult, fj),
                        )
                    else:
                        warned_unhandled.add(
                            (eq.field_name, term.operator, target_field),
                        )
            continue

        # ---------------- First-order equation ----------------
        assert eq.time_derivative_order == 1
        this_slot_red = orig_to_reduced[layout.field_slot_map[eq.field_name]]
        for term_idx, term in enumerate(eq.rhs_terms):
            target_field = term.field
            is_vel_ref = target_field.startswith("v_")
            base_field = target_field[2:] if is_vel_ref else target_field
            decomp = _OPERATOR_DECOMP[term.operator]
            mult = multiplier_cache[term.operator]
            t_order = decomp.time_order
            if t_order > 0:
                warned_unhandled.add((eq.field_name, term.operator, target_field))
                continue
            if base_field in dyn_field_idx:
                target_slot = (
                    layout.velocity_slot_map[base_field]
                    if is_vel_ref
                    else layout.field_slot_map[base_field]
                )
                target_red = orig_to_reduced[target_slot]
                _emit_term(
                    A_dd,
                    this_slot_red,
                    target_red,
                    term,
                    mult,
                    eq_idx=eq_idx,
                    term_idx=term_idx,
                )
            elif base_field in c_idx_map:
                cj = c_idx_map[base_field]
                out = A_dc_vel if is_vel_ref else A_dc_field
                _emit_term(
                    out,
                    this_slot_red,
                    cj,
                    term,
                    mult,
                    eq_idx=eq_idx,
                    term_idx=term_idx,
                )

    # ---- Substitute deferred d2_t(dyn) / first_dt(v_dyn) in constraint RHS ----
    # ẍ_fj = (M⁻¹·K)_{fj,k}·q_k + (M⁻¹·D)_{fj,k}·v_k
    # which equals A_dd[vel_slot(fj), field_slot(k)] for K and
    #              A_dd[vel_slot(fj), vel_slot(k)] for D, already M⁻¹-scaled.
    # The deferred-term contribution: K_cd[ci_block, :] += conv(term) @ A_dd[vel_slot(fj)_block, :]
    # GH #444 note: unlike the dyn-row velc/ddc lists below, this list
    # correctly uses scale=1.0 — the M⁻¹ it needs is the DYNAMICAL row's,
    # and it arrives through the composed A_dd factor (whose velocity rows
    # were emitted with velocity_row_scale). Constraint rows themselves
    # carry no kinetic coefficient (time_order == 0).
    for ci, eq_idx, term_idx, term, mult, fj in deferred_terms_constraint:
        fj_vel_red = orig_to_reduced[layout.velocity_slot_map[dyn_field_names[fj]]]
        temp = _term_conv_block(term, mult, scale=1.0, eq_idx=eq_idx, term_idx=term_idx)
        vel_block = A_dd[
            fj_vel_red * n_modes : (fj_vel_red + 1) * n_modes,
            :,
        ]
        K_cd[ci * n_modes : (ci + 1) * n_modes, :] += temp @ vel_block

    # ---- Schur eliminate constraint fields ------------------------------
    if n_c_tot == 0:
        recovery = np.zeros((0, n_dyn_tot), dtype=np.complex128)
    else:
        # Solve K_cc · q_c + K_cd · q_d = 0  →  q_c = -K_cc⁻¹ · K_cd · q_d.
        # GH #459 (mirrored from the path-2 S_cc fix): rank-revealing
        # pseudoinverse — a near-singular localized K_cc under plain
        # solve/lstsq produced ~1e12-scale recovery that poisoned the
        # pencil downstream; pinv solves consistent (possibly redundant)
        # systems exactly with the documented min-norm convention.
        recovery = -(
            np.asarray(np.linalg.pinv(K_cc, rcond=1e-10), dtype=np.complex128) @ K_cd
        )

    # ---- Compose A_reduced from dynamical block + constraint substitution ----
    A_reduced = A_dd + A_dc_field @ recovery

    # ---- Handle deferred dyn-RHS terms referencing v_c -------------------
    # v_c = d/dt(q_c) = d/dt(recovery · y_d) = recovery · dy_d/dt = recovery · A_reduced · y_d
    # A dyn-RHS term "coeff · op · v_c" contributes (conv(term) @ recovery_block @ A_reduced) to dv/dt rows.
    # This creates a cyclic dependence in A_reduced; handled via B_lhs pre-solve below.
    # First accumulate "linear-in-A_reduced" contribution:
    has_vel_coupling = False
    vel_coupling_mat = np.zeros_like(
        A_reduced
    )  # term × recovery, will be composed with A_reduced
    for vel_red, eq_idx, term_idx, term, mult, cj, row_scale in deferred_terms_dyn_velc:
        # GH #444: apply the emitting row's M⁻¹ (recorded at deferral time).
        # _term_conv_block handles both scalar and position-dependent
        # (ndarray) scales — the latter folds M⁻¹(x) into the coefficient
        # before the FFT (GH #427), so arbitrary kinetics are covered.
        temp = _term_conv_block(
            term, mult, scale=row_scale, eq_idx=eq_idx, term_idx=term_idx
        )
        # vel_c block in recovery = recovery[cj_block, :]
        recovery_block = recovery[cj * n_modes : (cj + 1) * n_modes, :]
        # Contribution: (vel_red row block) += temp @ recovery_block @ A_reduced
        # We store the prefactor: vel_coupling_mat[vel_red_block, :] += temp @ recovery_block
        vel_coupling_mat[
            vel_red * n_modes : (vel_red + 1) * n_modes,
            :,
        ] += temp @ recovery_block
        has_vel_coupling = True

    # If A_dc_vel is non-trivial (first-order eq RHS → v_c), its contribution is
    # A_dc_vel @ v_c = A_dc_vel @ recovery · A_reduced → folded into vel_coupling.
    # A_dc_vel is zero-size when the spec has no constraint fields (n_c == 0,
    # reachable post-GH #427: pos-dep kinetics + time-derivative RHS ops but
    # no algebraic constraints route here) — np.max on a zero-size array
    # raises, and there is nothing to fold.
    if A_dc_vel.size > 0 and np.max(np.abs(A_dc_vel)) > 1e-15:
        vel_coupling_mat += A_dc_vel @ recovery
        has_vel_coupling = True

    # ---- Final composition of the pencil B·ẏ = A·y ----------------------
    #   A = A_dd + A_dc_field·recovery   (accumulated above)
    #   B = I_adj − B_mass − vel_coupling_mat
    # where I_adj zeroes the promoted velocity-row diagonals (their LHS is
    # algebraic — no v̇_self, GH #457). When B has mass content the pencil
    # is composed by the deflating-subspace engine (singular B: the
    # promoted sector's constraint chains are deflated exactly, with the
    # manifold projector returned for on-manifold ICs); when only velocity
    # coupling is present, B = I − vc is the classic invertible pre-solve.
    manifold_proj_out: NDArray[np.complex128] | None = None
    if has_b_mass:
        eye_adj = np.eye(n_dyn_tot, dtype=np.complex128)
        for fname in promoted_fields:
            vel_red_p = orig_to_reduced[layout.velocity_slot_map[fname]]
            for m in range(n_modes):
                idx = vel_red_p * n_modes + m
                eye_adj[idx, idx] = 0.0
        if constraints_as_slots:
            # Residual constraint rows are algebraic pencil rows too.
            for cname in constraint_field_names:
                c_red = orig_to_reduced[layout.field_slot_map[cname]]
                for m in range(n_modes):
                    idx = c_red * n_modes + m
                    eye_adj[idx, idx] = 0.0
        b_total = eye_adj - B_mass - vel_coupling_mat
        A_reduced, manifold_proj_out = _pencil_deflate(
            A_reduced,
            b_total,
            context="pos-dep evolution pencil",
            diagnostics=diagnostics,
            tag=("evolution", 0),
        )
        if diagnostics is not None:
            _record_pin_overlap(
                diagnostics,
                layout,
                orig_to_reduced,
                n_modes,
                n_dyn_tot // n_modes,
                "evolution",
                per_mode=False,
            )
    elif has_vel_coupling:
        B_lhs = np.eye(A_reduced.shape[0], dtype=np.complex128) - vel_coupling_mat
        try:
            A_reduced = np.linalg.solve(B_lhs, A_reduced)
        except np.linalg.LinAlgError:
            logger.warning(
                "B_lhs singular in pos-dep velocity-coupling solve; falling back to lstsq",
            )
            A_reduced = np.linalg.lstsq(B_lhs, A_reduced, rcond=1e-12)[0]

    if warned_unhandled:
        warnings.warn(
            f"Modal solver: {len(warned_unhandled)} RHS term type(s) not yet "
            f"handled in pos-dep + constraint path (see logger for details).",
            stacklevel=2,
        )
        for fld, op, target in sorted(warned_unhandled):
            logger.warning(
                "Unhandled RHS term: equation %s, operator %s, target %s",
                fld,
                op,
                target,
            )

    return (
        A_reduced.astype(np.complex128, copy=False),
        recovery.astype(np.complex128, copy=False),
        # GH #468 slots mode: constraint values live in the evolved
        # slots (o2r covers them) — no recovery-based reconstruction.
        [] if constraints_as_slots else constraint_field_names,
        orig_to_reduced,
        manifold_proj_out,
    )


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
    """Check if the spec has any position-dependent coefficient — RHS or kinetic.

    This is the ROUTING predicate that selects the ĉ(k−k′) convolution
    builders over the constant-coefficient per-mode ones. Position
    dependence is a REGIME property: whether it lives in an RHS
    coefficient c(x) or in the kinetic coefficient M(x), the k-space
    object is a mode-coupling matrix and only the convolution builders
    represent it. The kinetic case (GH #427) is handled by folding the
    per-grid-point M⁻¹(x) profile into each velocity-row coefficient in
    real space before the FFT (see :func:`_conv_block_with_kinetic`) —
    mathematically identical to the mass-side M̂⁻¹(k−k′) convolution, so
    no separate mass-side machinery is needed.

    Callers outside :func:`solve_modal` use this for eligibility:
    ``modal_jax`` bails to the scipy modal path on ``True`` (the JAX
    per-mode path is constant-coefficient), which is the correct verdict
    for position-dependent kinetics too.
    """
    if spec.has_position_dependent_kinetic():
        return True
    for eq in spec.equations:
        for term in eq.rhs_terms:
            if term.position_dependent:
                return True
    return False


def _suppress_tachyonic_noise(  # pyright: ignore[reportUnusedFunction]
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


class SingularPencilError(RuntimeError):
    """The per-mode pencil λB − A is singular (det ≡ 0): true gauge freedom.

    Raised by :func:`_pencil_deflate` when an (α, β) ≈ (0, 0) pair is
    found — some state direction is genuinely undetermined by the
    equations. There is no correct evolution for such a direction; a
    silent choice (freeze, minimum-norm, Tikhonov) is a wrong theory
    with exit code 0. Gauge-fix the theory (``[[gauge]]`` in the TOML)
    or file the spec on GH #457's follow-up.
    """


class PencilDiagnostics:
    """Pinned (gauge-quotiented) directions recorded during deflation.

    GH #468 "pin + certify": every state direction the gauge quotient
    sets to zero (the explicit min-norm choice) is recorded here, and the
    builders turn the records into ``pin_overlap[mode, slot]`` — the
    squared norm of the pinned subspace's rows on that state slot. A
    measurement whose field support carries no pinned content is
    provably independent of the pin (certified); one that does is
    flagged with the magnitude by ``tidal measure`` instead of being
    reported as if it were gauge-invariant.
    """

    def __init__(self) -> None:
        self.entries: list[tuple[tuple[str, int], NDArray[np.complex128], float]] = []
        self.slot_names: list[str] = []
        self.pin_overlap: NDArray[np.float64] | None = None  # (n_modes, n_slots)
        self.pinned_dims: int = 0
        self.determination_floor: float = 0.0

    def record(
        self, tag: tuple[str, int], pinned: NDArray[np.complex128], tau: float
    ) -> None:
        """Record one deflation's pinned basis (``(n, w_dim)``, may be empty)."""
        self.entries.append((tag, pinned, tau))
        self.determination_floor = max(self.determination_floor, tau)

    def to_dict(self) -> dict[str, Any]:
        """JSON-ready summary for run metadata (``solver_diagnostics``)."""
        overlap = (
            self.pin_overlap
            if self.pin_overlap is not None
            else np.zeros((0, len(self.slot_names)))
        )
        return {
            "slot_names": list(self.slot_names),
            "pin_overlap_max": [
                float(v) for v in (overlap.max(axis=0) if overlap.size else [])
            ],
            "pin_overlap": [[float(v) for v in row] for row in overlap],
            "pinned_dims": int(self.pinned_dims),
            "determination_floor": float(self.determination_floor),
        }


def _reduced_slot_names(
    layout: StateLayout, orig_to_reduced: dict[int, int], n_red: int
) -> list[str]:
    """Names of the reduced state slots (fields and ``v_``-velocities)."""
    names = [""] * n_red
    for name, s in layout.field_slot_map.items():
        if s in orig_to_reduced:
            names[orig_to_reduced[s]] = name
    for name, s in layout.velocity_slot_map.items():
        if s in orig_to_reduced:
            names[orig_to_reduced[s]] = f"v_{name}"
    return names


def _record_pin_overlap(
    diagnostics: PencilDiagnostics,
    layout: StateLayout,
    orig_to_reduced: dict[int, int],
    n_modes: int,
    n_red: int,
    kind: str,
    *,
    per_mode: bool,
) -> None:
    """Turn recorded pinned bases into ``pin_overlap[mode, slot]``.

    ``per_mode=True``: one pencil per mode, rows = reduced slots.
    ``per_mode=False``: one full-FFT pencil, row ``r*n_modes + m``.
    """
    overlap = np.zeros((n_modes, n_red))
    total = 0
    for (tag_kind, m), w, _tau in diagnostics.entries:
        if tag_kind != kind or w.shape[1] == 0:
            continue
        total += w.shape[1]
        rows = np.sum(np.abs(w) ** 2, axis=1)
        if per_mode:
            overlap[m, :] += rows
        else:
            overlap += rows.reshape(n_red, n_modes).T
    diagnostics.slot_names = _reduced_slot_names(layout, orig_to_reduced, n_red)
    diagnostics.pin_overlap = overlap
    diagnostics.pinned_dims = total


def _gauge_quotient_pencil(
    A_n: NDArray[np.complex128],
    B_n: NDArray[np.complex128],
    probe_points: tuple[complex, ...],
    context: str,
    *,
    tau: float,
) -> tuple[
    NDArray[np.complex128],
    NDArray[np.complex128],
    NDArray[np.complex128],
    NDArray[np.complex128],
]:
    """Quotient a singular pencil by its Kronecker column-chain space.

    For a singular pencil, ``λ₀B − A`` is singular at EVERY λ₀, and its
    null vectors at generic probe points are evaluations of the
    polynomial (Kronecker chain) null family — e.g. for a byte-identical
    equation pair (GH #465), ``x(λ) = q_diff + λ·v_diff``. The union of
    those null spaces over a few generic λ₀ spans the gauge chain space
    W: the state combinations the equations genuinely never determine.

    The quotient: gauge-fix W to ZERO (the min-norm convention the
    pre-#457 machinery documented for undetermined modes), restrict the
    pencil to W's orthogonal complement C, and row-compress the
    restricted equations to a square regular pencil. Returns
    ``(A_red, B_red, C, W)`` with the reduced pencil expressed in the
    C-basis and ``W`` the orthonormal basis of the pinned (zeroed)
    directions; the caller composes the generator back through C, uses
    ``C·C^H`` as the (gauge-zeroing) part of the IC projection, and
    records ``W`` for the measurement-time gauge certificate (GH #468).

    KNOWN LIMIT (GH #474): the row compression assumes the restricted
    equations have exactly ``n_red`` independent combinations — true for
    exact chains (GH #465), false for near-null quotients of localized
    backgrounds, where it discards genuine equations; the general fix is
    the Kronecker-like staircase reduction (GH #473).

    Emits a warning — a spec that hits this every run should be
    re-derived without the redundancy (GH #465) or gauge-fixed in the
    theory (``[[gauge]]``). Raises :class:`SingularPencilError` when the
    quotient does not yield a regular pencil (structure beyond gauge
    chains).
    """
    import warnings as _warnings  # noqa: PLC0415

    n = A_n.shape[0]
    tol = tau  # determination floor from the adaptive loop (GH #468)

    # Gauge chain space W: union of null(λ₀B − A) over generic probes.
    null_vecs: list[NDArray[np.complex128]] = []
    for lam0 in probe_points:
        _u, s_p, vh_p = np.linalg.svd(lam0 * B_n - A_n)
        scale = float(s_p[0]) if s_p.size else 1.0
        n_null = int(np.sum(s_p <= tol * max(scale, 1.0)))
        if n_null:
            null_vecs.append(vh_p[-n_null:, :].conj().T)  # (n, n_null)
    if not null_vecs:
        ctx = f" ({context})" if context else ""
        msg = f"gauge quotient called on a regular pencil{ctx}"
        raise SingularPencilError(msg)
    w_stack = np.hstack(null_vecs)
    # Orthonormalize W and take the complement C from a full SVD.
    u_w, s_w, _ = np.linalg.svd(w_stack, full_matrices=True)
    w_dim = int(np.sum(s_w > 1e-10 * max(float(s_w[0]), 1.0)))
    c_basis = u_w[:, w_dim:] if w_dim < n else np.zeros((n, 0), dtype=np.complex128)
    n_red = c_basis.shape[1]
    if n_red == 0:
        ctx = f" ({context})" if context else ""
        msg = f"Singular pencil{ctx}: every state direction is gauge. GH #457."
        raise SingularPencilError(msg)

    # Restrict to the complement and row-compress to a square system.
    a_rest = A_n @ c_basis  # (n, n_red)
    b_rest = B_n @ c_basis
    u_r, s_r, _ = np.linalg.svd(np.hstack([a_rest, b_rest]))
    row_scale = float(s_r[0]) if s_r.size else 1.0
    row_rank = int(np.sum(s_r > tol * max(row_scale, 1.0)))
    if row_rank < n_red:
        ctx = f" ({context})" if context else ""
        msg = (
            f"Singular pencil{ctx}: the gauge quotient leaves "
            f"{n_red - row_rank} dependent equation combination(s) beyond "
            f"the chain space — structure no frozen-gauge choice "
            f"resolves. Gauge-fix the theory ([[gauge]]). GH #457."
        )
        raise SingularPencilError(msg)
    u_keep = u_r[:, :n_red].conj().T  # (n_red, n)
    a_red = u_keep @ a_rest
    b_red = u_keep @ b_rest
    _warnings.warn(
        f"Singular pencil ({context}): {w_dim} undetermined (gauge) state "
        f"combination(s) quotiented out and set to zero at determination "
        f"floor τ = {tau:.1e} — the explicit form of the min-norm choice "
        f"the pre-#457 machinery made silently. For a localized "
        f"background this is the far-field gauge restoration (GH #468); "
        f"for redundant equations re-derive the theory (GH #465) or "
        f"gauge-fix it ([[gauge]] in the TOML).",
        stacklevel=4,
    )
    return a_red, b_red, c_basis, u_w[:, :w_dim]


def _pencil_deflate(
    A_pencil: NDArray[np.complex128],
    B_pencil: NDArray[np.complex128],
    *,
    context: str = "",
    diagnostics: PencilDiagnostics | None = None,
    tag: tuple[str, int] = ("pencil", 0),
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Deflating-subspace reduction of one mode's pencil ``B·ẏ = A·y``.

    The general engine for DAE structure of ANY index (GH #457): ordered
    QZ separates finite eigenvalues (dynamics) from infinite ones
    (algebraic constraint chains). For the HOMOGENEOUS free evolution the
    infinite sector is exactly ``z_i ≡ 0`` (the nilpotent chain
    ``T_ii·ż_i = S_ii·z_i`` with ``S_ii`` invertible forces it), so the
    constraint manifold is the finite deflating subspace with an
    ORTHONORMAL basis ``Z_f`` and reduced generator
    ``G_red = T_ff⁻¹·S_ff`` (triangular solve). Conditioning note: the
    v0.31 eigendecomposition retirement (cond(V) blowup, Higham 2008
    §2.3) does not apply here — Q/Z are unitary and the only inversion
    is triangular.

    Returns
    -------
    A_eff : ndarray (n, n)
        ``Z_f · G_red · Z_f^H`` — the evolution generator on the full
        slot space. Off-manifold components are annihilated by ``A_eff``
        and therefore stay frozen under ``exp(A_eff·t)`` (matching the
        existing null-projection semantics); project the IC with
        ``proj`` for on-manifold initialization.
    proj : ndarray (n, n)
        ``Z_f · Z_f^H`` — orthogonal projector onto the manifold.

    Raises
    ------
    SingularPencilError
        If the pencil is singular ((α, β) ≈ (0, 0) pair): true gauge
        freedom that no evolution choice can silently resolve.
    """
    n = A_pencil.shape[0]
    # ROW-equilibrate the pencil: each row (equation) is scaled to unit
    # max magnitude. This is a left diagonal preconditioner — it changes
    # neither the solution set, the generalized eigenvalues, nor the
    # RIGHT deflating subspaces (Z), so no un-scaling is needed below.
    # Without it, rows carrying small physical scales (e.g. B0²-weight
    # constraint rows) make both the regularity probe and QZ's (α, β)
    # magnitudes unreliable — the same scale-blindness as the det-gated
    # regularizer of GH #459.
    row_scale = np.maximum(
        np.max(np.abs(A_pencil), axis=1), np.max(np.abs(B_pencil), axis=1)
    )
    row_scale[row_scale == 0.0] = 1.0  # all-zero row → the probe rejects it
    A_n = A_pencil / row_scale[:, np.newaxis]
    B_n = B_pencil / row_scale[:, np.newaxis]

    probe_points = (0.5488 + 0.7152j, -1.3113 + 0.4227j)

    def _probe_smin(
        a_mat: NDArray[np.complex128], b_mat: NDArray[np.complex128]
    ) -> float:
        smin = 0.0
        for lam0 in probe_points:
            sv = np.linalg.svd(lam0 * b_mat - a_mat, compute_uv=False)
            smin = max(smin, float(sv[-1]) if sv.size else 0.0)
        return smin

    # ---- Adaptive determination-floor quotient (GH #468) ----
    # Localized backgrounds produce a CONTINUOUS cascade of weakly
    # determined state combinations: the Gaussian B₀(x) tail restores
    # gauge symmetry smoothly toward the far field (measured on the
    # E.cal constraints-as-slots pencil: singular values 2e-13 … 1e-9
    # with no gap), so no single tolerance separates "gauge" from
    # "weakly determined". Instead: quotient at a determination floor
    # τ, verify the retained sector against the deflation contract, and
    # RAISE τ decade by decade until the contract passes. Each step is
    # SELF-VERIFIED, so the choice cannot silently violate the retained
    # equations; the quotiented combinations are gauge-fixed to zero
    # (the documented min-norm convention) and the warning reports the
    # final floor and dimension. If no floor up to τ_max yields a
    # contract-consistent operator, the last refusal propagates
    # (genuine near-index breakdown, GH #467).
    tau0 = n * 1e-12
    # tau_max = 1e-4 is the ACCEPTED determination floor for the
    # localized promoted class (user decision 2026-08-26, GH #468/#470):
    # the far-field gauge-restoration continuum makes machine-precision
    # closure unachievable at float64, and content determined more
    # weakly than ~1e-4 is physically fading (B0 -> 0 there). The
    # contract tolerance is floor-linked, so every acceptance is
    # quantified and warned.
    tau_max = 1e-4
    tau = tau0
    last_err: SingularPencilError | None = None
    while tau <= tau_max:
        pinned = np.zeros((n, 0), dtype=np.complex128)
        try:
            if _probe_smin(A_n, B_n) < tau:
                a_red, b_red, c_basis, pinned = _gauge_quotient_pencil(
                    A_n, B_n, probe_points, context, tau=tau
                )
                a_eff_red, proj_red = _deflate_regular_pencil(
                    a_red, b_red, context, tau
                )
                a_eff_full = (c_basis @ a_eff_red) @ c_basis.conj().T
                proj_full = (c_basis @ proj_red) @ c_basis.conj().T
                result = (
                    np.asarray(a_eff_full, dtype=np.complex128),
                    np.asarray(proj_full, dtype=np.complex128),
                )
            else:
                result = _deflate_regular_pencil(A_n, B_n, context, tau)
        except SingularPencilError as exc:
            last_err = exc
            tau *= 10.0
            continue
        if diagnostics is not None:
            diagnostics.record(tag, pinned, tau)
        if tau > tau0:
            import warnings as _warnings  # noqa: PLC0415

            _warnings.warn(
                f"Pencil deflation ({context}): determination floor raised "
                f"to τ = {tau:.1e} — state content determined more weakly "
                f"than τ is treated as algebraic/gauge (contract-verified "
                f"at each step; near-breakdown or far-field content beyond "
                f"the floor is quotiented). See GH #467/#468.",
                stacklevel=3,
            )
        return result
    assert last_err is not None
    raise last_err


def _deflate_regular_pencil(
    A_n: NDArray[np.complex128],
    B_n: NDArray[np.complex128],
    context: str,
    tau: float,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128]]:
    """Ordered-QZ deflation of a (nominally regular) equilibrated pencil.

    ``tau`` is the determination floor: generalized-eigenvalue pairs
    with ``|β|/(|α|+|β|) ≤ tau``-scale are classified infinite. The
    deflation contract ``B·Z_f·G = A·Z_f`` is verified before returning
    — its violation raises :class:`SingularPencilError` (caught by the
    adaptive floor loop in :func:`_pencil_deflate`, which retries with a
    larger floor; the final failure is the honest near-index-breakdown
    refusal of GH #467).
    """
    from scipy.linalg import ordqz, solve_triangular  # noqa: PLC0415

    n = A_n.shape[0]

    def _is_finite_eig(
        alpha: NDArray[np.complex128], beta: NDArray[np.complex128]
    ) -> NDArray[np.bool_]:
        return np.abs(beta) > tau * (np.abs(alpha) + np.abs(beta) + tau)

    qz_out = ordqz(  # pyright: ignore[reportUnknownVariableType, reportCallIssue]
        A_n,
        B_n,
        sort=_is_finite_eig,  # pyright: ignore[reportArgumentType]
        output="complex",
    )
    s_mat = np.asarray(qz_out[0], dtype=np.complex128)
    t_mat = np.asarray(qz_out[1], dtype=np.complex128)
    alpha = np.asarray(qz_out[2], dtype=np.complex128)
    beta = np.asarray(qz_out[3], dtype=np.complex128)
    z_mat = np.asarray(qz_out[5], dtype=np.complex128)

    finite = _is_finite_eig(alpha, beta)
    n_fin = int(np.sum(finite))
    # ordqz sorts selected (finite) eigenvalues to the leading block.
    z_f = z_mat[:, :n_fin]
    if n_fin == 0:
        zero = np.zeros((n, n), dtype=np.complex128)
        return zero, zero
    t_ff = t_mat[:n_fin, :n_fin]
    s_ff = s_mat[:n_fin, :n_fin]
    # Row equilibration is a LEFT transformation: the equilibrated system
    # has the same solutions, so G_red is the physical generator directly.
    g_red = np.asarray(solve_triangular(t_ff, s_ff), dtype=np.complex128)

    # SELF-VERIFICATION of the deflation contract B·Z_f·G = A·Z_f
    # (measured necessity, GH #467): near constraint-index breakdown a
    # genuine eigenvalue |λ| ~ 1/eps^(1/2..1) appears whose float64
    # representation cannot be both faithful and integrable — keeping it
    # pollutes the composition (eps·|λ| = O(1) on ordinary rows),
    # zeroing it violates the equations. Harmless huge-but-consistent
    # modes (β-rounding of true infinite pairs) pass this check; only a
    # composition that fails its own defining equation refuses.
    tol_contract = max(1e-8, tau)
    contract_res = _contract_residual(A_n, B_n, z_f, g_red)
    if contract_res > tol_contract:
        # ordqz's eigenvalue REORDERING through many ill-conditioned
        # infinite pairs is numerically unstable at large n (measured:
        # contract 3.7e2 → 1.3e-2 on the 600-dim E.cal pencil, never
        # reaching tolerance). Second rung: a reordering-free
        # construction — right eigenvectors of the retained finite
        # spectrum (|λ| ≤ 1/τ), orthonormalized, with G from the
        # contract's own least-squares — judged by the SAME contract.
        z_f, g_red, contract_res = _deflate_eig_subspace(A_n, B_n, tau)
        if contract_res > tol_contract:
            ctx = f" ({context})" if context else ""
            msg = (
                f"Deflation contract violated{ctx}: best residual "
                f"{contract_res:.2e} exceeds the floor-linked tolerance "
                f"{tol_contract:.1e} at determination floor τ = {tau:.1e} "
                f"(both the ordered-QZ and the reordering-free "
                f"constructions). The pencil is within machine precision "
                f"of constraint-index breakdown at these parameters — no "
                f"double-precision operator is both faithful to the "
                f"equations and integrable. This parameter point cannot "
                f"be simulated as posed. GH #467/#468."
            )
            raise SingularPencilError(msg)

    a_eff = (z_f @ g_red) @ z_f.conj().T
    proj = z_f @ z_f.conj().T
    return np.asarray(a_eff, dtype=np.complex128), np.asarray(proj, dtype=np.complex128)


def _contract_residual(
    A_n: NDArray[np.complex128],
    B_n: NDArray[np.complex128],
    z_f: NDArray[np.complex128],
    g_red: NDArray[np.complex128],
) -> float:
    """Relative residual of the deflation contract ``B·Z_f·G = A·Z_f``."""
    contract_lhs = B_n @ (z_f @ g_red)
    contract_rhs = A_n @ z_f
    contract_scale = max(float(np.max(np.abs(contract_rhs))), 1e-300)
    return float(np.max(np.abs(contract_lhs - contract_rhs))) / contract_scale


def _deflate_eig_subspace(
    A_n: NDArray[np.complex128],
    B_n: NDArray[np.complex128],
    tau: float,
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], float]:
    """Reordering-free deflating-subspace construction (GH #468).

    Span the retained finite spectrum (|lambda| <= 1/tau — the
    determination floor's rate horizon) with right eigenvectors,
    orthonormalize, and obtain the reduced generator from the contract's
    own least-squares ``G = argmin ||B·Q·G − A·Q||``. Eigenvector
    conditioning is NOT trusted: the returned contract residual is the
    sole acceptance criterion (checked by the caller against the
    floor-linked tolerance). Avoids ordqz's eigenvalue reordering, whose
    swaps through large ill-conditioned infinite clusters are the
    measured unstable step on big localized pencils.
    """
    from scipy.linalg import eig, lstsq  # noqa: PLC0415

    out = eig(A_n, B_n, homogeneous_eigvals=True)
    alpha, beta = out[0]
    vr = np.asarray(out[1], dtype=np.complex128)
    with np.errstate(divide="ignore", invalid="ignore"):
        lam = np.where(
            np.abs(beta) > 0.0,
            alpha / np.where(np.abs(beta) > 0.0, beta, 1.0),
            np.inf,
        )
    # Retention-horizon scan: the contract residual is non-monotonic in
    # the horizon R (too small excludes genuine modes; too large retains
    # ill-conditioned near-horizon content) — take the best over decades.
    # NOTE (measured 2026-08-26, GH #468/#470): on the localized E.cal
    # class even the best retention carries RESOLUTION-DIVERGENT spurious
    # growth (maxRe +85→+472 for N 12→32), so the floor-linked contract
    # refusal is expected to stand there; this construction serves the
    # cases where a clean subspace exists.
    n = A_n.shape[0]
    best: tuple[NDArray[np.complex128], NDArray[np.complex128], float] | None = None
    for horizon in (1e1, 1e2, 1e3, 1e4):
        keep = np.isfinite(lam) & (np.abs(lam) <= horizon)
        if not np.any(keep):
            continue
        q_basis, _ = np.linalg.qr(vr[:, keep])
        bq = B_n @ q_basis
        aq = A_n @ q_basis
        lstsq_out = lstsq(bq, aq)
        assert lstsq_out is not None  # scipy stubs type the result Optional
        g_red = np.asarray(lstsq_out[0], dtype=np.complex128)
        res = _contract_residual(A_n, B_n, q_basis, g_red)
        if best is None or res < best[2]:
            best = (np.asarray(q_basis, dtype=np.complex128), g_red, res)
    if best is None:
        zero = np.zeros((n, 0), dtype=np.complex128)
        return zero, np.zeros((0, 0), dtype=np.complex128), 0.0
    return best


def _build_m_with_null_projection(
    A_block: NDArray[np.complex128],
    B_block: NDArray[np.complex128] | None,
) -> NDArray[np.complex128]:
    """Compute ``M = B^-1 · A`` per mode, handling rank-deficient ``B``.

    For null directions of ``B`` (gauge / non-propagating components), the
    evolution generator ``M`` is set to zero so those directions stay at IC
    under ``exp(M · t)``. For ``range(B)`` we solve normally.

    Parameters
    ----------
    A_block : ndarray, shape (n_modes, bs, bs)
        RHS matrix per mode.
    B_block : ndarray, shape (n_modes, bs, bs) or None
        Kinetic matrix per mode. ``None`` means ``B = I`` (no LHS coupling),
        in which case ``M = A`` directly.

    Returns
    -------
    M_block : ndarray, shape (n_modes, bs, bs)
        Evolution generator. Null(B) directions are zero (IC-frozen under
        ``exp(M · t)``).
    """
    if B_block is None:
        return A_block.copy()
    n_modes_local, bs_local, _ = A_block.shape
    M_block = np.zeros_like(A_block)
    for m in range(n_modes_local):
        A_m = A_block[m]
        B_m = B_block[m]
        _u, s, Vt = cast(
            "tuple[NDArray[np.complex128], NDArray[np.float64], NDArray[np.complex128]]",
            np.linalg.svd(B_m),
        )
        thresh = s[0] * 1e-10 if s.size > 0 and s[0] > 0 else 1e-14
        rank = int(np.sum(s > thresh))
        if rank == bs_local:
            M_block[m] = np.linalg.solve(B_m, A_m)
        elif rank > 0:
            # Project onto range(B): A_red, B_red of size (rank, rank)
            Vphys = Vt[:rank].T  # (bs, rank)
            A_red = Vphys.conj().T @ A_m @ Vphys
            B_red = Vphys.conj().T @ B_m @ Vphys
            M_red = np.linalg.solve(B_red, A_red)
            # Lift back; null directions stay at zero (IC-frozen).
            M_block[m] = Vphys @ M_red @ Vphys.conj().T
        # else: rank == 0 — entire block is null(B); M stays zero (IC-frozen).
    return M_block


def _evolve_per_mode_pade(
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
    collect_pass0_state: list[dict[str, Any]] | None = None,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.complex128] | None,
    NDArray[np.complex128] | None,
]:
    """Per-mode matrix-exponential evolution via Padé precompute + matvec (path D).

    For each independent field block, computes ``M = B^-1·A`` per mode (with
    null-space projection if ``B`` is rank-deficient), then ``exp(M·dt)`` once
    per mode via ``scipy.linalg.expm`` (Higham 2009 Padé scaling-and-
    squaring). The snapshot loop is pure ``O(bs^2)`` matvec, vectorized across
    modes via einsum --- structurally identical to the eigendecomposition
    snapshot loop, but **robust for arbitrary ``cond(V)``** and wall-time
    competitive (0.49--1.06× across the campaign workload envelope; see
    ``docs/tex/modal_solver.tex`` §"Robust Matrix-Exponential Evolution").

    This is the default Pass 0 path. The eigendecomposition path remains
    available via ``collect_eigendata`` for Pass 1 perturbative-Duhamel
    callers (which need per-eigenvalue kernels) --- see
    :func:`_evolve_per_mode` and ``perturbative_solver.tex`` §Pass~1.

    References
    ----------
    Higham, N.J. (2009). "The scaling and squaring method for the matrix
    exponential revisited." SIAM Review 51(4):747-764.
    Higham, N.J. (2008). Functions of Matrices, SIAM, §2.3 (cond(V)
    abandonment threshold).
    """
    import os  # noqa: PLC0415

    import scipy.linalg as sla  # type: ignore[import-untyped]  # noqa: PLC0415

    from tidal.solver._exceptions import SimulationDivergedError  # noqa: PLC0415

    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_snapshots = len(t_eval)
    n_modes = y0_hat.shape[1]
    t0 = float(t_eval[0])
    t_end_rel = float(t_eval[-1] - t0)

    # Sparse-IC mode-skip kill-switch (#327 follow-up): when env var
    # ``TIDAL_MODAL_SPARSE_IC=0`` is set, fall back to the legacy
    # all-modes precompute and the all-modes expm-based divergence
    # pre-check.  Default (unset / non-zero) enables sparse-IC +
    # eigenvalue-based pre-check.  See
    # ``docs/tex/stability_probe.tex`` §sec:ieee-floor for the
    # noise-floor mechanism the eigenvalue pre-check captures.
    sparse_ic_enabled: bool = os.environ.get("TIDAL_MODAL_SPARSE_IC", "1") != "0"

    # Block detection: same as the eigendecomposition path.
    n_check = min(3, A_modes.shape[0])
    combined = np.max(np.abs(A_modes[:n_check]), axis=0)
    blocks = find_independent_blocks(combined)

    # Per-block precompute: M_block, exp_M_dt (uniform spacing) or M only,
    # plus the active-modes mask used by sparse-IC (#327 follow-up).
    # block_state entries are (slot_indices, M_block, exp_M_dt or None,
    #                          y_curr, y0_block, active_modes).
    # y_curr has shape (n_modes, bs) — the per-mode current Fourier state.
    if n_snapshots > 1:
        dts = np.diff(t_eval)
        uniform = bool(np.allclose(dts, dts[0]))
        dt_step = float(dts[0]) if uniform else None
    else:
        uniform = True
        dt_step = None

    block_state: list[
        tuple[
            list[int],
            NDArray[np.complex128],  # M_block (n_modes, bs, bs)
            NDArray[np.complex128] | None,  # exp_M_dt (n_modes, bs, bs) or None
            NDArray[np.complex128],  # y_curr (n_modes, bs); starts as y0_block.T
            NDArray[np.complex128],  # y0_block (bs, n_modes); kept for non-uniform path
            NDArray[np.intp],  # active_modes — indices to evolve via expm
        ]
    ] = []
    # Worst-case max real eigenvalue across all blocks/modes — drives the
    # eigenvalue-based divergence pre-check (replaces the per-mode expm
    # pre-check; reuses ``eigvals_diag`` already computed for the growth-
    # rate warning).  See plan §Stage 2.
    max_real_eig: float = 0.0

    for block_slots in blocks:
        y0_block = y0_hat[block_slots, :]  # (bs, n_modes)
        if np.max(np.abs(y0_block)) < 1e-15:
            continue
        idx = np.array(block_slots)
        A_block = A_modes[:, idx[:, None], idx[None, :]]
        B_block = (
            B_modes[:, idx[:, None], idx[None, :]] if B_modes is not None else None
        )
        M_block = _build_m_with_null_projection(A_block, B_block)

        # Sparse-IC active-mode mask (#327 follow-up).  Modes whose IC
        # amplitude is at or below the IEEE 754 FFT floor (~ N·ε ≈
        # 10⁻¹⁴ for plane-wave IC) carry no physical signal; their
        # ``y_curr`` evolves from noise and contributes nothing to the
        # ifft'd output.  Skipping the per-mode ``expm`` for those bins
        # cuts the precompute cost by 50–65× on inference workloads
        # (plane-wave IC has ~1 active rfft bin out of N//2+1).  When
        # ``sparse_ic_enabled`` is False the mask covers every mode,
        # recovering legacy behavior bit-exactly.
        n_modes_block: int = int(M_block.shape[0])
        if sparse_ic_enabled:
            y0_norms = np.linalg.norm(y0_block, axis=0)  # (n_modes,)
            norm_floor = 1e-12 * max(1.0, float(np.max(np.abs(y0_block))))
            active_modes: NDArray[np.intp] = np.flatnonzero(
                y0_norms > norm_floor,
            ).astype(np.intp)
        else:
            active_modes = np.arange(n_modes_block, dtype=np.intp)

        # Use ``zeros_like`` (not ``empty_like``) so inactive modes'
        # ``exp_M_dt[m]`` is exactly zero — the snapshot matvec
        # (``np.einsum("mij,mj->mi", exp_M_dt, y_curr)``) then produces
        # 0 for those modes regardless of ``y_curr[m]``'s value, with
        # no risk of NaN/inf garbage propagating from uninitialized
        # memory.
        exp_M_dt: NDArray[np.complex128] | None = None
        if uniform and dt_step is not None and dt_step > 0:
            exp_M_dt = np.zeros_like(M_block)
            for m in active_modes:
                exp_M_dt[m] = sla.expm(M_block[m] * dt_step)  # pyright: ignore[reportUnknownArgumentType]

        # Diagnostic: warn if growth factor over t_end_rel exceeds exp(30).
        # path D doesn't form eigenvectors, but eigenvalues alone are cheap
        # and reliable for this purpose. Mirrors the eigen-path warning so
        # downstream tests / users get the same overflow signal.  Reuse the
        # max real eigenvalue here for the divergence pre-check (#327
        # follow-up Stage 2) — eigvals are computed regardless, so the
        # extra ``np.max(np.real(...))`` is cents of microseconds.
        if t_end_rel > 0:
            eigvals_diag = np.linalg.eigvals(M_block)
            _warn_eigenvalue_growth(
                eigvals_diag,  # pyright: ignore[reportArgumentType]
                t_end_rel,
                context="per-mode (Padé)",
            )
            max_real_eig = max(max_real_eig, float(np.max(np.real(eigvals_diag))))

        y_curr = y0_block.T.astype(np.complex128, copy=True)  # (n_modes, bs)
        block_state.append(
            (list(block_slots), M_block, exp_M_dt, y_curr, y0_block, active_modes),
        )

        # Collect Pass 0 state for the augmented-exp Pass 1 solver. Replaces
        # the legacy {V, V_inv, D_diag, alpha} eigendata schema; consumed by
        # :func:`_evolve_duhamel_per_mode` to build the (2bs × 2bs) augmented
        # operator [[A, S], [0, A]] without any eigendecomposition.
        # See docs/tex/modal_solver.tex §"Robust Matrix-Exponential Evolution".
        if collect_pass0_state is not None:
            collect_pass0_state.append(
                {
                    "slot_indices": list(block_slots),
                    "M_block": M_block.copy(),
                    "y0_block": y0_block.copy(),
                },
            )

    # Outputs
    snapshots = np.zeros((n_snapshots, n_slots * n_pts))
    times = np.zeros(n_snapshots)
    fourier_snaps: NDArray[np.complex128] | None = (
        np.zeros((n_snapshots, n_slots, n_modes), dtype=np.complex128)
        if return_fourier
        else None
    )
    deriv_fourier_snaps: NDArray[np.complex128] | None = (
        np.zeros((n_snapshots, n_slots, n_modes), dtype=np.complex128)
        if return_derivative_fourier
        else None
    )
    y_hat_t = np.zeros((n_slots, n_modes), dtype=np.complex128)
    dy_hat_t: NDArray[np.complex128] | None = (
        np.zeros((n_slots, n_modes), dtype=np.complex128)
        if return_derivative_fourier
        else None
    )

    # Divergence pre-check at t_end (#327 follow-up Stage 2).
    #
    # When ``sparse_ic_enabled`` (default), use the worst per-mode
    # max real eigenvalue we already computed for the growth-rate
    # warning: the IEEE 754 FFT floor at non-fundamental bins is
    # ~1e-14 (see docs/tex/stability_probe.tex §sec:ieee-floor), so a
    # tachyonic mode with eigenvalue γ at t_end_rel amplifies to
    # ~1e-14 · exp(γ · t_end_rel).  Trip the guard if this exceeds
    # ``divergence_threshold × initial_max_amp``.  This is bit-
    # equivalent to the per-mode expm pre-check on tachyonic samples
    # (sample 2860 in the docs verifies the mechanism) but cents-of-
    # microseconds vs ~10–50 ms.
    #
    # When the kill-switch is set, fall back to the legacy per-mode
    # expm pre-check.  Same physics, same numbers — slower.
    initial_max_amp: float = 0.0
    divergence_threshold: float = 100.0
    if t_end_rel > 0:
        initial_physical = _ifft_slots(y0_hat, layout, grid)
        initial_max_amp = max(float(np.max(np.abs(initial_physical))), 1e-15)

        if sparse_ic_enabled:
            # Eigenvalue-based check.  ``ieee_fft_floor`` is the
            # double-precision FFT round-off scaled by the IC amplitude
            # (see stability_probe.tex Eq. for sample 2860); 1e-14 is a
            # conservative upper bound.
            ieee_fft_floor: float = 1e-14
            predicted_log_growth = max_real_eig * t_end_rel
            # Cap the exponent before ``math.exp`` to avoid overflow
            # warnings — anything > log(1e308) is divergence by any
            # measure.
            log_overflow_cap = float(np.log(1e300))
            if predicted_log_growth > log_overflow_cap:
                predicted_floor_amp = float("inf")
            else:
                predicted_floor_amp = ieee_fft_floor * float(
                    np.exp(predicted_log_growth),
                )
            predicted_ratio = predicted_floor_amp / initial_max_amp
            if (
                not np.isfinite(predicted_ratio)
                or predicted_ratio > divergence_threshold
            ):
                msg = (
                    f"Simulation predicted to diverge: max real eigenvalue "
                    f"{max_real_eig:.4g} amplifies the IEEE 754 FFT floor "
                    f"{ieee_fft_floor:.0e} to ratio "
                    f"{predicted_ratio:.2e} at t={t_eval[-1]:.4g} "
                    f"(threshold {divergence_threshold:.0e}). Fields would "
                    f"leave the perturbative regime (linearized approximation "
                    f"invalid).  Rejecting pre-evolution based on the "
                    f"eigenvalue-based pre-check (#327 follow-up Stage 2)."
                )
                raise SimulationDivergedError(msg)
        else:
            # Legacy expm-based pre-check — kept under the kill-switch
            # for bit-exact rollback on a single env var.
            y_hat_predict = np.zeros((n_slots, n_modes), dtype=np.complex128)
            for (
                block_slots,
                M_block,
                _exp_M_dt,
                _y_curr,
                y0_block,
                _active_modes,
            ) in block_state:
                y_pred = np.empty(
                    (M_block.shape[0], M_block.shape[1]),
                    dtype=np.complex128,
                )
                for m in range(M_block.shape[0]):
                    y_pred[m] = (
                        sla.expm(M_block[m] * t_end_rel) @ y0_block[:, m]  # pyright: ignore[reportUnknownArgumentType]
                    )
                y_hat_predict[block_slots, :] = y_pred.T
            predicted_physical = _ifft_slots(y_hat_predict, layout, grid)
            predicted_max = float(np.max(np.abs(predicted_physical)))
            if (
                not np.isfinite(predicted_max)
                or predicted_max / initial_max_amp > divergence_threshold
            ):
                msg = (
                    f"Simulation predicted to diverge: amplitude at t={t_eval[-1]:.4g} "
                    f"would reach ratio {predicted_max / initial_max_amp:.2e} "
                    f"(threshold {divergence_threshold:.0e}). Fields would leave "
                    f"the perturbative regime (linearized approximation invalid). "
                    f"Rejecting pre-evolution based on Padé matrix exponential."
                )
                raise SimulationDivergedError(msg)

    # Time loop. For uniform spacing, advance y_curr by exp_M_dt between
    # snapshots; for non-uniform, recompute exp(M · t_rel) per snapshot per
    # mode (slower but correct).
    for ti, t in enumerate(t_eval):
        dt_from_t0 = float(t - t0)
        y_hat_t[:] = 0.0
        if dy_hat_t is not None:
            dy_hat_t[:] = 0.0

        for k, (
            block_slots,
            M_block,
            exp_M_dt,
            y_curr,
            y0_block,
            active_modes,
        ) in enumerate(block_state):
            if ti == 0:
                # First snapshot: y_curr already equals y0_block.T
                pass
            elif uniform and exp_M_dt is not None:
                # Advance one step: y_curr <- exp_M_dt @ y_curr.
                # Inactive modes have ``exp_M_dt[m] == 0`` (set by
                # ``np.zeros_like`` above) so the einsum produces 0
                # for them — bit-exactly correct since they carried
                # no signal at IC.
                y_curr_new = np.einsum("mij,mj->mi", exp_M_dt, y_curr)
                block_state[k] = (
                    block_slots,
                    M_block,
                    exp_M_dt,
                    y_curr_new,
                    y0_block,
                    active_modes,
                )
                y_curr = y_curr_new
            else:
                # Non-uniform spacing: compute exp(M · dt_from_t0) ·
                # y0_block per mode.  Sparse-IC: only iterate active
                # modes; inactive modes are zeroed in y_curr_new (the
                # ifft sum then drops their contribution).
                y_curr_new = np.zeros_like(y_curr)
                for m in active_modes:
                    y_curr_new[m] = (
                        sla.expm(M_block[m] * dt_from_t0) @ y0_block[:, m]  # pyright: ignore[reportUnknownArgumentType]
                    )
                block_state[k] = (
                    block_slots,
                    M_block,
                    exp_M_dt,
                    y_curr_new,
                    y0_block,
                    active_modes,
                )
                y_curr = y_curr_new

            y_hat_t[block_slots, :] = y_curr.T

            if dy_hat_t is not None:
                # dy/dt = M @ y(t) — exact, no numerical differentiation.
                dy = np.einsum("mij,mj->mi", M_block, y_curr)
                dy_hat_t[block_slots, :] = dy.T

        if fourier_snaps is not None:
            fourier_snaps[ti] = y_hat_t
        if deriv_fourier_snaps is not None and dy_hat_t is not None:
            deriv_fourier_snaps[ti] = dy_hat_t

        y_physical = _ifft_slots(y_hat_t, layout, grid)
        snapshots[ti] = y_physical
        times[ti] = t

        # Runtime divergence guard (same physics threshold as eigen path).
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
            for tj in range(ti + 1, n_snapshots):
                snapshots[tj] = y_physical
                times[tj] = t_eval[tj]
            raise SimulationDivergedError(msg)

        if snapshot_callback is not None:
            snapshot_callback(t, y_physical)
        if progress is not None:
            progress.update(t)

    return times, snapshots, fourier_snaps, deriv_fourier_snaps


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
    collect_eigendata: list[dict[str, Any]] | None = None,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.complex128] | None,
    NDArray[np.complex128] | None,
]:
    """Evolve system with per-mode independent matrices (constant coefficients).

    Thin wrapper over :func:`_evolve_per_mode_pade`. The eigendecomposition
    backend was retired in v0.31+ in favor of unconditional Padé scaling-
    and-squaring (path D), which is robust for arbitrary ``cond(V)`` AND
    wall-time competitive (0.49--1.06× across the campaign workload
    envelope; see ``docs/tex/modal_solver.tex`` §"Robust Matrix-Exponential
    Evolution").

    The ``collect_eigendata`` parameter is kept for backwards compatibility
    with Pass 1 callers but now collects ``{M_block, y0_block, slot_indices}``
    per block instead of the legacy ``{V, V_inv, D_diag, alpha}``. Pass 1
    consumes this via the augmented matrix exponential
    (Al-Mohy & Higham 2011 §5.2) — see :func:`_evolve_duhamel_per_mode`.

    Parameters
    ----------
    A_modes : ndarray, shape (n_modes, n_slots, n_slots)
        Per-mode RHS matrix.
    y0_hat : ndarray, shape (n_slots, n_modes)
        Initial condition in Fourier space.
    B_modes : ndarray, optional, shape (n_modes, n_slots, n_slots)
        Per-mode kinetic matrix; ``M = B^-1 A`` is solved with null-space
        projection if ``B`` is rank-deficient (gauge / constraint).
    collect_eigendata : list of dict, optional
        When provided, populated per block with ``{slot_indices, M_block,
        y0_block}`` for downstream Pass 1 augmented-exp evaluation.

    Raises
    ------
    SimulationDivergedError
        If field amplitudes grow beyond 100x the initial maximum (linearized
        equations are no longer physical).

    References
    ----------
    Higham 2009 (Padé scaling-and-squaring), Higham 2008 §2.3 (eigendecomp
    abandonment threshold), Al-Mohy & Higham 2011 §5.2 (augmented-exp Pass 1).
    """
    return _evolve_per_mode_pade(
        A_modes,
        y0_hat,
        t_eval,
        layout,
        grid,
        snapshot_callback,
        progress,
        return_fourier=return_fourier,
        return_derivative_fourier=return_derivative_fourier,
        B_modes=B_modes,
        collect_pass0_state=collect_eigendata,
    )


def _evolve_full_matrix(
    A_full: NDArray[np.complex128],
    y0_hat: NDArray[np.complex128],
    t_eval: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
    snapshot_callback: Callable[[float, NDArray[np.float64]], None] | None,
    progress: SimulationProgress | None,
    *,
    full_spectrum: bool = False,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evolve system with full coupled matrix (position-dependent coefficients).

    A_full has shape (n_total, n_total) where n_total = n_slots x n_modes.
    y0_hat has shape (n_slots, n_modes).

    ``full_spectrum`` selects the snapshot inverse transform: True for the
    full fftn basis (GH #445, the pos-dep convolution paths), False for
    the rfft half-spectrum basis.

    Uses ``scipy.sparse.linalg.expm_multiply`` (Al-Mohy & Higham 2011) to compute
    exp(A·t)·y₀ at each output time without eigendecomposition.

    **Why not eigendecomposition?**  Non-normal convolution matrices have
    eigenvalues with positive real parts that cause exp(λ·t) overflow while
    exp(A·t)·y₀ remains bounded.  Eigendecomposition was permanently retired
    in v0.31 (GH #320).

    **Why not dense Padé (scipy.linalg.expm)?**  For strongly non-normal matrices
    (GH #367: max real eigenvalue > 0 from pseudospectral Fourier truncation of
    a localized background field), exp(A·t) has norm ~10⁹ while the physically
    correct exp(A·t)·y₀ is ~0.005 — 12 orders of magnitude of catastrophic
    cancellation, impossible at float64 precision.  Dense Padé computes exp(A·t)
    correctly but the subsequent product with y₀ loses all significant digits.
    **Workaround for this case: ``--scheme cvode``** (integrates in physical space,
    avoids the Fourier coupling matrix entirely).

    **Sparse optimization:** For Gaussian B₀(x), the convolution kernel decays
    exponentially → effectively banded matrix.  Entries below 1e-14 × max|A|
    are zeroed; if density < 30% the matrix converts to sparse CSC format.  This
    accelerates expm_multiply's internal matrix-vector products.

    Ref: Al-Mohy & Higham (2011), "Computing the Action of the Matrix
    Exponential", SIAM J. Sci. Comput. 33(2):488-511.
    """
    import scipy.sparse  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
    from scipy.sparse.linalg import (  # noqa: PLC0415  # pyright: ignore[reportMissingTypeStubs]
        expm_multiply,  # pyright: ignore[reportUnknownVariableType]
    )

    from tidal.solver._exceptions import SimulationDivergedError  # noqa: PLC0415

    n_slots = layout.num_slots
    n_pts = layout.num_points
    n_modes = y0_hat.shape[1]
    n_snapshots = len(t_eval)

    # Historical context: GH #367 (fixed v0.42.0) and GH #379 (fixed 2026-05-24)
    # both caused this guard to fire when the convolution matrix had wrong-sign
    # eigenvalues at large k. Both root causes were in the matrix builder (path 3
    # missing M⁻¹; path 4 not existing) rather than in expm_multiply; see the
    # module docstring at modal.py:14-104 for the path-by-path contract. The
    # post-evolution amplitude check below (10⁶ threshold) remains as a defensive
    # net for future builder regressions.
    #
    # We do NOT pre-flight gate on k_max·t_end here because that heuristic would
    # falsely reject tests like test_position_dependent_correctness that mark
    # coefficients as symbolically position-dependent but analytically constant
    # ("D*(1 + 0*x[])" at k_max·t_end ≈ 40 is correct because the convolution
    # matrix is diagonal). Detecting the artifact properly requires computing
    # eigenvalues of A_full (O(n³)) — more expensive than just running the solve
    # and checking the post-evolution amplitude.

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
        inverse_transform = _ifft_slots_full if full_spectrum else _ifft_slots
        for ti in range(n_snapshots):
            t = float(t_eval[ti])
            y_hat_t = y_all[ti].reshape(n_slots, n_modes)
            y_physical = inverse_transform(y_hat_t, layout, grid)
            snapshots[ti] = y_physical
            times[ti] = t

            if snapshot_callback is not None:
                snapshot_callback(t, y_physical)
            if progress is not None:
                progress.update(t)
    else:
        # Single time point or t0 == t_end
        inverse_transform = _ifft_slots_full if full_spectrum else _ifft_slots
        for ti, t in enumerate(t_eval):
            if t == t0:
                y_evolved = y0_flat.copy()
            else:
                y_evolved = np.asarray(
                    expm_multiply(A_op, y0_flat, start=t0, stop=float(t), num=2)[-1],
                    dtype=np.complex128,
                )
            y_hat_t = y_evolved.reshape(n_slots, n_modes)
            y_physical = inverse_transform(y_hat_t, layout, grid)
            snapshots[ti] = y_physical
            times[ti] = t

            if snapshot_callback is not None:
                snapshot_callback(t, y_physical)
            if progress is not None:
                progress.update(t)

    # --- Post-evolution divergence check (GH #367, GH #379) -------------------
    # Defensive net: with both #367 (M⁻¹ in convolution path) and #379 (Schur
    # elimination for pos-dep + constraints) now fixed, reaching this branch
    # indicates either (a) a new builder regression in path 3/4, or (b) genuine
    # tachyonic instability in the physics (run with α=0 to disambiguate;
    # see docs/CLAUDE.md "t_end independence test").
    y0_max = float(np.max(np.abs(y0_flat))) or 1.0
    y_final_max = float(np.max(np.abs(snapshots[-1])))
    if y_final_max > y0_max * 1e6 and t_end > t0:
        # Classify the theory by slot composition to make the message actionable.
        has_constraints = any(s.kind == "constraint" for s in layout.slots)
        if has_constraints:
            hint = (
                "Theory has algebraic constraints; reaching this guard with "
                "constraints present indicates a regression in path 4 "
                "(_build_convolution_matrix_with_constraints) — file a bug "
                "with the theory's JSON and parameter values. CVODE and IDA "
                "typically cannot run constraint-bearing theories with d2_t "
                "RHS operators; rerun with α=0 to check for tachyonic "
                "instability vs. matrix-builder bug."
            )
        else:
            hint = (
                "If physics is expected to be smooth at these parameters, "
                "this indicates a builder regression in path 3 "
                "(_build_convolution_matrix) — file a bug. For genuine "
                "tachyonic instability, reduce coupling or shorten t_end. "
                "Fallback: --scheme cvode (works for theories without "
                "constraint or d2_t-RHS structure)."
            )
        msg = (
            f"Simulation diverged in _evolve_full_matrix: final amplitude "
            f"{y_final_max:.2e} is >{y_final_max / y0_max:.1e}× the initial "
            f"amplitude {y0_max:.2e}.\n{hint}"
        )
        raise SimulationDivergedError(msg)

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
    return_eigendata: bool = False,
    _zero_nyquist: bool = True,
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
    return_eigendata : bool, optional
        When True, attaches a ``"eigendata"`` key to the returned
        ``SolverResult`` containing per-block eigendecomposition data
        (``V``, ``D_diag``, ``V_inv``, ``alpha``), the ``mode_k`` grid,
        the dynamical ``state_layout``, and (when constraints are
        present) a ``schur_ops`` sub-dict with the recovery matrix and
        the original → reduced slot map. Used by the v6 iterative
        :class:`PerturbativeSolver` to drive a Pass 1 closed-form
        Duhamel evaluation that reuses Pass 0's eigendecomposition.
        Not supported for position-dependent coefficient systems (they
        use ``expm_multiply`` without an explicit eigendecomposition);
        attempting to combine the two raises ``NotImplementedError``.

    Returns
    -------
    SolverResult
        Dict with keys: ``t``, ``y``, ``success``, ``message``. When
        ``return_eigendata=True`` also has ``"eigendata"``.

    Raises
    ------
    ValueError
        If ``spec`` contains time-derivative orders > 2, or incompatible
        options are given for the position-dependent coefficient path.
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    # #278: solve_modal's companion-matrix builder assumes
    # time_derivative_order <= 2 on every equation. The `d4_t` and
    # `mixed_T4_*` operators in ``_EXACT_MULTIPLIERS`` exist for
    # _build_source_matrix_k (correction spec only) and must never
    # reach here in the base spec — if they do, the eigendecomposition
    # would silently operate on a truncated system. Fail loudly.
    max_to = max((eq.time_derivative_order for eq in spec.equations), default=0)
    if max_to > 2:
        high_fields = [
            eq.field_name for eq in spec.equations if eq.time_derivative_order > 2
        ]
        msg = (
            f"solve_modal requires time_derivative_order <= 2 on every "
            f"equation; got max={max_to} on field(s) {high_fields}. Use "
            f"PerturbativeSolver.solve() with a [perturbation] section "
            f"in theory.toml so base_spec() demotes higher-derivative "
            f"LHS kinetics to algebraic constraints before they reach "
            f"the modal backend. See "
            f"docs/PERTURBATIVE_REDUCTION_IMPLEMENTATION.md."
        )
        raise ValueError(msg)

    # GH #427: position-dependent kinetics are supported here — the
    # kinetics-aware _has_position_dependent_terms routes such specs to the
    # convolution builders, which fold M⁻¹(x) into the real-space
    # coefficients (see _conv_block_with_kinetic). The former GH #421 entry
    # guard remains only in _build_evolution_matrices (the per-mode genEig
    # engine, reached directly by the stability probe and modal_jax) and in
    # solve_modal_pass1 (per-mode Duhamel).

    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, parameters or {})

    # Detect constraint fields
    has_constraints = any(eq.time_derivative_order == 0 for eq in spec.equations)
    if not has_constraints:
        warn_frozen_constraints(layout, "modal")

    # Pass 1 eigendata collector (v6 Stage 3). Populated per block by
    # _evolve_per_mode when return_eigendata=True. Empty (and ignored) in
    # the position-dependent expm_multiply path where there is no
    # eigendecomposition to expose.
    block_eigendata: list[dict[str, Any]] = []

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
    # _zero_nyquist=False reproduces the pre-fix behavior for benchmarking.
    if _zero_nyquist:
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
    # GH #468 "pin + certify": the builders record every gauge-quotiented
    # (pinned) direction here; the per-slot overlap goes into the run's
    # metadata so `tidal measure` can certify each observable.
    diagnostics = PencilDiagnostics()

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
    #
    # GH #379: when pos-dep AND constraints are both present, we route to
    # _build_convolution_matrix_with_constraints (the convolution analog of
    # the genEig Schur path). Pre-fix, the dispatch had `and not has_pos_dep`
    # which silently skipped Schur elimination for this case and treated
    # algebraic constraints `q_c = RHS` as `dq_c/dt = RHS`, producing
    # divergent dynamics at any nonzero coupling.
    needs_reduction = has_constraints or has_time_ops
    constraint_vel_arrays: dict[str, NDArray[np.float64]] = {}  # populated below

    # Variables assigned inside `if needs_reduction` and used later in the same
    # function — initialized here so pyright can track definite assignment.
    dyn_layout: StateLayout | None = None
    recovery_matrix: NDArray[np.complex128] | None = None
    c_names: list[str] | None = None
    orig_to_reduced: dict[int, int] | None = None
    Scc_inv_modes: NDArray[np.complex128] | None = None
    Scc_singular_mask_modes: NDArray[np.bool_] | None = None

    if needs_reduction and has_pos_dep:
        # GH #379: pos-dep + constraints — convolution Schur path
        if return_eigendata:
            msg = (
                "return_eigendata=True is not supported for position-"
                "dependent coefficient systems with algebraic constraints "
                "(GH #379 path uses expm_multiply Krylov with no explicit "
                "eigendecomposition). v6 Pass 1 Duhamel requires constant "
                "coefficients."
            )
            raise NotImplementedError(msg)

        # GH #445: the convolution paths operate in the FULL fftn basis —
        # pointwise multiplication by c(x) is exactly C-linear there,
        # unlike the rfft half-spectrum (see _fft_slots_full).
        k_grid_full = _build_k_grid(_build_k_axes_full(grid))
        full_shape = tuple(grid.shape)
        y0_hat_full = _fft_slots_full(y0, layout, grid)
        if _zero_nyquist:
            _zero_nyquist_full(y0_hat_full, grid)

        (
            A_reduced_2d,
            recovery_2d,
            c_names_list,
            orig_to_reduced_map,
            conv_manifold_proj,
        ) = _build_convolution_matrix_with_constraints(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid_full,
            full_shape,
            diagnostics=diagnostics,
        )
        recovery_matrix = recovery_2d
        c_names = c_names_list
        orig_to_reduced = orig_to_reduced_map

        n_modes = y0_hat_full.shape[1]
        n_pts = layout.num_points
        n_dyn_slots = len(orig_to_reduced)
        n_c_count = len(c_names)

        # Build reduced StateLayout (dynamical slots only, in reduced order)
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

        # Extract dynamical IC in reduced ordering, as (n_dyn_slots, n_modes)
        y0_hat_dyn = np.zeros((n_dyn_slots, n_modes), dtype=np.complex128)
        for orig_si, red_pos in orig_to_reduced.items():
            y0_hat_dyn[red_pos] = y0_hat_full[orig_si]

        # GH #457: promoted second-order sector — project the IC onto
        # the constraint manifold (off-manifold directions are frozen by
        # the deflated generator; a component left off the manifold
        # would persist as a constant constraint violation).
        if conv_manifold_proj is not None:
            y0_hat_dyn = (conv_manifold_proj @ y0_hat_dyn.ravel()).reshape(
                y0_hat_dyn.shape
            )

        # Evolve the reduced matrix via Krylov expm_multiply
        times, dyn_snapshots = _evolve_full_matrix(
            A_reduced_2d,
            y0_hat_dyn,
            t_eval,
            dyn_layout,
            grid,
            None,  # callback handled below with full state
            progress,
            full_spectrum=True,
        )

        # Reconstruct full state (dyn slots + constraint slots)
        n_full = layout.num_slots * n_pts
        snapshots = np.zeros((len(t_eval), n_full))

        # GH #468 slots mode: c_names is empty and constraint values live
        # in the evolved slots; their velocities come from the
        # generator's constraint-slot rows below.
        slots_mode_constraints: list[tuple[str, int]] = []
        if conv_manifold_proj is not None:
            for _si, _slot in enumerate(layout.slots):
                if _slot.kind == "constraint" and _si in orig_to_reduced:
                    slots_mode_constraints.append(
                        (_slot.field_name, orig_to_reduced[_si])
                    )
        for c_name in c_names:
            constraint_vel_arrays[c_name] = np.zeros((len(t_eval), *grid.shape))
        for c_name, _red in slots_mode_constraints:
            constraint_vel_arrays[c_name] = np.zeros((len(t_eval), *grid.shape))

        for ti in range(len(t_eval)):
            dyn_phys = dyn_snapshots[ti]
            # FFT each dynamical slot back to k-space (full basis, GH #445)
            # to apply recovery
            y_hat_dyn_t = np.zeros((n_dyn_slots, n_modes), dtype=np.complex128)
            for red_si in range(n_dyn_slots):
                field_data = dyn_phys[red_si * n_pts : (red_si + 1) * n_pts].reshape(
                    grid.shape
                )
                y_hat_dyn_t[red_si] = np.fft.fftn(field_data).ravel()

            # Flatten and apply recovery (q_c_hat = recovery @ y_d_hat)
            y_hat_dyn_flat = y_hat_dyn_t.ravel()
            c_hat = (recovery_matrix @ y_hat_dyn_flat).reshape(n_c_count, n_modes)

            # Constraint velocities: v_c_hat = recovery @ (A_reduced @ y_d_hat)
            dy_hat_dyn_flat = A_reduced_2d @ y_hat_dyn_flat
            v_c_hat = (recovery_matrix @ dy_hat_dyn_flat).reshape(n_c_count, n_modes)

            # Assemble full physical state
            full_state = np.zeros(n_full)
            for orig_si, red_pos in orig_to_reduced.items():
                full_state[orig_si * n_pts : (orig_si + 1) * n_pts] = dyn_phys[
                    red_pos * n_pts : (red_pos + 1) * n_pts
                ]
            for ci, c_name in enumerate(c_names):
                c_slot = layout.field_slot_map[c_name]
                c_phys = np.fft.ifftn(c_hat[ci].reshape(grid.shape)).ravel()
                full_state[c_slot * n_pts : (c_slot + 1) * n_pts] = np.real(c_phys)
                v_c_phys = np.fft.ifftn(v_c_hat[ci].reshape(grid.shape))
                constraint_vel_arrays[c_name][ti] = np.real(v_c_phys)
            for c_name, c_red in slots_mode_constraints:
                # v_C = Ċ from the generator's constraint-slot row (exact).
                vc_hat_row = dy_hat_dyn_flat[c_red * n_modes : (c_red + 1) * n_modes]
                v_c_phys = np.fft.ifftn(vc_hat_row.reshape(grid.shape))
                constraint_vel_arrays[c_name][ti] = np.real(v_c_phys)

            snapshots[ti] = full_state
            if snapshot_callback is not None:
                snapshot_callback(t_eval[ti], full_state)

        n_total_red = A_reduced_2d.shape[0]
        method_desc = (
            f"expm_multiply with constraint Schur elimination "
            f"({n_c_count} constraint fields, {n_dyn_slots} dynamical slots, "
            f"{n_total_red}x{n_total_red} reduced matrix, pos-dep)"
        )
    elif needs_reduction:
        (
            A_reduced,
            B_lhs_modes,
            recovery_matrix,
            _v_recovery_gen,  # unused — constraint vel from eigendata
            c_names,
            orig_to_reduced,
            Scc_inv_modes,
            Scc_singular_mask_modes,
            manifold_proj_modes,
        ) = _build_evolution_matrices(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid,
            rfft_shape,
            diagnostics=diagnostics,
        )

        n_dyn = A_reduced.shape[1]
        n_modes = y0_hat.shape[1]
        n_pts = layout.num_points

        # Extract dynamical IC in reduced ordering
        y0_hat_dyn = np.zeros((n_dyn, n_modes), dtype=np.complex128)
        for orig_si, red_pos in orig_to_reduced.items():
            y0_hat_dyn[red_pos] = y0_hat[orig_si]

        # GH #457: when the pencil engine reduced the sector (promoted or
        # singular-mass), project the IC onto the constraint manifold —
        # off-manifold directions are frozen by the generator, so a
        # component left off the manifold would persist as a constant constraint
        # violation for the whole run.
        if manifold_proj_modes is not None:
            y0_hat_dyn = np.einsum("mij,jm->im", manifold_proj_modes, y0_hat_dyn)

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
            collect_eigendata=block_eigendata if return_eigendata else None,
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
            collect_eigendata=block_eigendata if return_eigendata else None,
        )
        method_desc = "per-mode eigendecomposition (constant coefficients)"
    else:
        # Position-dependent coefficients: full convolution matrix in the
        # FULL fftn basis (GH #445; see _fft_slots_full).
        k_grid_full = _build_k_grid(_build_k_axes_full(grid))
        full_shape = tuple(grid.shape)
        y0_hat_full = _fft_slots_full(y0, layout, grid)
        if _zero_nyquist:
            _zero_nyquist_full(y0_hat_full, grid)

        A_full = _build_convolution_matrix(
            spec,
            layout,
            grid,
            coeff_eval,
            k_grid_full,
            full_shape,
        )
        if return_eigendata:
            msg = (
                "return_eigendata=True is not supported for position-"
                "dependent coefficient systems (expm_multiply Krylov path "
                "has no explicit eigendecomposition). v6 Pass 1 Duhamel "
                "requires constant coefficients."
            )
            raise NotImplementedError(msg)
        times, snapshots = _evolve_full_matrix(
            A_full,
            y0_hat_full,
            t_eval,
            layout,
            grid,
            snapshot_callback,
            progress,
            full_spectrum=True,
        )
        n_total = A_full.shape[0]
        method_desc = f"expm_multiply ({n_total}x{n_total}, position-dependent)"

    if progress is not None:
        progress.finish()

    # A modal run that never needed the pencil engine (invertible fast path,
    # plain position-dependent builder) pinned nothing — say so explicitly
    # so `tidal measure` can certify "no pins" rather than stay silent.
    if not diagnostics.slot_names:
        diagnostics.slot_names = [*layout.field_slot_map] + [
            f"v_{name}" for name in layout.velocity_slot_map
        ]
        diagnostics.pin_overlap = np.zeros((1, len(diagnostics.slot_names)))

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
    if diagnostics.slot_names:
        result["diagnostics"] = diagnostics.to_dict()  # type: ignore[typeddict-unknown-key]

    # Stage 3 (v6): expose Pass 0 eigendecomposition for Pass 1 Duhamel.
    if return_eigendata:
        eigendata: dict[str, Any] = {
            "blocks": block_eigendata,
            "mode_k": k_grid,
            "state_layout": dyn_layout if needs_reduction else layout,
        }
        if needs_reduction:
            assert recovery_matrix is not None
            assert c_names is not None
            assert orig_to_reduced is not None
            schur_ops_out: dict[str, Any] = {
                "recovery_matrix": recovery_matrix,
                "constraint_field_names": tuple(c_names),
                "orig_to_reduced": dict(orig_to_reduced),
            }
            # #290: expose S_cc_inv + singular_mask so the
            # augmented Pass 1 constraint recovery can apply
            # ``h_c¹ += S_cc_inv · (K·d^n_t(h_c⁰) − corr_1(y⁰))``.
            if Scc_inv_modes is not None:
                schur_ops_out["S_cc_inv"] = Scc_inv_modes
            if Scc_singular_mask_modes is not None:
                schur_ops_out["S_cc_singular_mask"] = Scc_singular_mask_modes
            eigendata["schur_ops"] = schur_ops_out
        result["eigendata"] = eigendata  # type: ignore[typeddict-unknown-key]
    return result


# ---------------------------------------------------------------------------
# Pass 1: closed-form Duhamel evolution (v6 Stage 4)
# ---------------------------------------------------------------------------
# Solves  dy⁽¹⁾/dt = A(k)·y⁽¹⁾ + M_src(k)·y⁰(t),  y⁽¹⁾(0) = 0
# where y⁰(t) is the Pass 0 solution provided via eigendata. Reuses the
# Stage 3 eigendecomposition of A(k) so the only per-mode cost is a
# dense (bs × bs) β-projection and (n_t) Duhamel-kernel evaluations.


def _build_source_matrix_k(
    correction_spec: EquationSystem,
    layout: StateLayout,
    grid: GridInfo,
    coeff_eval: CoefficientEvaluator,
    k_grid: list[NDArray[np.float64]],
    rfft_shape: tuple[int, ...],
    schur_ops: dict[str, Any] | None = None,
) -> tuple[dict[int, NDArray[np.complex128]], list[dict[str, Any]]]:
    """Build M_src(k) for Pass 1: correction RHS as a linear source on the Pass 0 state.

    Shape ``(n_modes, n_slots, n_slots)``. For each correction RHS term
    ``c·op(field)`` appearing in equation ``eq``:

    - The row is the velocity slot when ``eq.time_derivative_order >= 2``
      (source enters the ``dv/dt`` entry of the companion system) or the
      field slot when ``eq.time_derivative_order == 1``.
    - The column is the target-field slot when ``term.field`` is
      dynamical. When it is a constraint field (eliminated by Schur in
      Pass 0), the row ``recovery_matrix[m, c_idx, :]`` from
      ``schur_ops`` is used to expand the reference into the dynamical
      subspace — ``h_c → Σ_j recovery[m, c_idx, j] · d_j``.

    ``correction_spec`` is typically ``spec.filter_by_order(n)`` for the
    Pass ``n`` correction (v6 plan, Stage 4).

    #293: returns a DICT keyed by operator time-derivative order
    instead of a single matrix. For a correction term with
    ``_OperatorDecomp.time_order == n``, the contribution goes into
    ``M_src[n]``; :func:`_evolve_duhamel_per_mode` pre-multiplies each
    order's source column by ``diag(λⁿ)`` before the Duhamel kernel so
    the ``d^n_t(target)`` operator contributes the correct eigenvalue
    power. The pre-fix code silently dropped ``time_order`` and any
    correction term like ``d²_t(dyn_field)`` produced wrong Pass 1
    output (latent bug, not triggered by shipped theories because
    constraint-row corrections route through R8 augmented recovery).

    Returns
    -------
    M_src_by_order : dict[int, ndarray]
        Maps time-derivative order ``n`` → ``(n_modes, n_slots,
        n_slots)`` complex array holding the source contributions of
        all correction terms whose operator has ``_OperatorDecomp.
        time_order == n``. Non-present keys (no terms at that order)
        are omitted. The empty dict represents "no sources".
    drops : list of diagnostic records
        Each entry describes a correction term or equation whose
        contribution was NOT written into any ``M_src[n]`` because
        either (a) the equation is for a demoted constraint field
        (row skip; routed to R8 augmented recovery), or (b) the
        target field reference could not be resolved to a dynamical
        or Schur-recoverable slot (column skip).

    Raises
    ------
    ValueError
        If an operator string from the correction spec cannot be resolved.
    """
    n_slots = layout.num_slots
    n_modes = int(np.prod(rfft_shape))

    def _resolve_target_slot(field_ref: str) -> int | None:
        if field_ref in layout.field_slot_map:
            return layout.field_slot_map[field_ref]
        if field_ref.startswith("v_") and field_ref[2:] in layout.velocity_slot_map:
            return layout.velocity_slot_map[field_ref[2:]]
        return None

    # Build constraint-field name → recovery row index map
    constraint_idx: dict[str, int] = {}
    recovery_matrix: NDArray[np.complex128] | None = None
    if schur_ops is not None:
        c_names_list = list(schur_ops.get("constraint_field_names", ()))
        constraint_idx.update({cname: ci for ci, cname in enumerate(c_names_list)})
        recovery_matrix = schur_ops.get("recovery_matrix")

    # Cache Fourier multipliers per operator
    multiplier_cache: dict[str, NDArray[np.complex128]] = {}
    for eq in correction_spec.equations:
        for term in eq.rhs_terms:
            op = term.operator
            if op in multiplier_cache:
                continue
            if op not in _EXACT_MULTIPLIERS:
                msg = (
                    f"Operator {op!r} in correction spec has no Fourier "
                    f"multiplier — modal Duhamel requires one of "
                    f"{sorted(_EXACT_MULTIPLIERS)}"
                )
                raise ValueError(msg)
            mult_fn = _EXACT_MULTIPLIERS[op]
            mult_val = mult_fn(k_grid)
            mult_full = np.broadcast_to(mult_val, rfft_shape)
            multiplier_cache[op] = mult_full.ravel().astype(np.complex128)

    # #293: per-time-order matrices so correction terms with
    # d^n_t operators scale by λⁿ in the Duhamel kernel. Key lookup is
    # cheap (small int keys), and we only allocate n_modes×n_slots²
    # arrays for orders that actually appear.
    M_src_by_order: dict[int, NDArray[np.complex128]] = {}

    # #301 / #302: Pass 1 source mirrors the Pass 0 fast-path contract —
    # each velocity-row contribution is scaled by M₀⁻¹ for the
    # corresponding equation's field. Without this, Pass 1 produces a
    # source that omits the M₀⁻¹(K₁ − M₁·M₀⁻¹·K₀) factor from the
    # perturbative identity; the error scales linearly with ε and breaks
    # Phase-3 canonicalized theories (e.g., Euler-Heisenberg with σ in
    # small_parameters — the synthesized corrections come out wrong).
    # Returns None when every M ≈ 1 so the existing fast path continues
    # with zero overhead for theories unaffected by #301.
    from tidal.solver._kinetic import (  # noqa: PLC0415
        build_inverse_kinetic_diag,
        velocity_row_scale,
    )

    # GH #427: `grid` is threaded through for consistency with the
    # convolution builders. The solve_modal_pass1 entry guard still refuses
    # position-dependent correction kinetics (per-mode Duhamel has no
    # convolution route), so array-valued entries cannot legitimately
    # appear here — the eq_scale check below enforces that invariant.
    m_inv_src = build_inverse_kinetic_diag(
        correction_spec,
        coeff_eval._parameters,  # noqa: SLF001  # type: ignore[reportPrivateUsage]
        grid,
    )

    def _get_m_src(n: int) -> NDArray[np.complex128]:
        if n not in M_src_by_order:
            M_src_by_order[n] = np.zeros(
                (n_modes, n_slots, n_slots),
                dtype=np.complex128,
            )
        return M_src_by_order[n]

    drops: list[dict[str, Any]] = []

    for eq_idx, eq in enumerate(correction_spec.equations):
        field_name = eq.field_name
        # Identify row
        if eq.time_derivative_order >= 2:
            if field_name not in layout.velocity_slot_map:
                # Equation for a field not in the reduced dynamical layout.
                # This is the "correction on a demoted constraint field"
                # case. The contribution is routed to
                # ``_compute_constraint_source_hat`` (R8 / #290) and
                # applied via augmented Schur recovery in
                # ``_assemble_full_state_pass_n``; NOT actually lost.
                # Record as "row-routed-to-augmented" so the caller can
                # distinguish this from genuine column-missing bugs.
                drops.append(
                    {
                        "field": field_name,
                        "operator": None,
                        "reason": "row-routed-to-augmented: correction equation is "
                        "on a demoted constraint field; O(ε¹) contribution is "
                        "applied via R8 augmented Schur recovery (#290) in the "
                        "driver, not via the Pass 1 dynamical source path",
                        "eq_idx": eq_idx,
                        "term_idx": None,
                        "n_terms": len(eq.rhs_terms),
                    },
                )
                continue
            row_slot = layout.velocity_slot_map[field_name]
        elif field_name in layout.field_slot_map:
            row_slot = layout.field_slot_map[field_name]
        else:
            # Time-order-1 equation on a demoted constraint field
            # (no velocity slot, no field slot in reduced layout). R8
            # augmented Schur recovery picks this up in the driver
            # (#290).
            if field_name in constraint_idx:
                drops.append(
                    {
                        "field": field_name,
                        "operator": None,
                        "reason": "row-routed-to-augmented: correction equation is "
                        "on a demoted constraint field; O(ε¹) contribution is "
                        "applied via R8 augmented Schur recovery (#290) in the "
                        "driver, not via the Pass 1 dynamical source path",
                        "eq_idx": eq_idx,
                        "term_idx": None,
                        "n_terms": len(eq.rhs_terms),
                    },
                )
            else:
                drops.append(
                    {
                        "field": field_name,
                        "operator": None,
                        "reason": "row-missing: correction equation on field not in "
                        "base-spec dynamical layout AND not a known constraint",
                        "eq_idx": eq_idx,
                        "term_idx": None,
                        "n_terms": len(eq.rhs_terms),
                    },
                )
            continue

        # #301 / #302: scale every contribution by 1/M₀ for this equation.
        # Mirrors the modal fast path's folding of M⁻¹ into the evolution
        # matrix (see _build_per_mode_matrices). Shared `velocity_row_scale`
        # helper enforces the cross-builder contract (see modal.py module
        # docstring; regression guard:
        # tests.test_solver_kinetic_consistency::TestAllModalPathsRespectKinetic).
        eq_scale = velocity_row_scale(field_name, m_inv_src)
        if isinstance(eq_scale, np.ndarray):
            # Invariant enforcement: the solve_modal_pass1 entry guard
            # refuses position-dependent correction kinetics before this
            # builder runs. A real-space M⁻¹(x) profile must never scale
            # the k-diagonal source (wrong semantics, n_points vs n_modes).
            msg = (
                f"Pass 1 source builder received a position-dependent "
                f"kinetic profile for field {field_name!r}; the per-mode "
                "Duhamel path requires constant kinetics (guard in "
                "solve_modal_pass1, GH #421/#427)."
            )
            raise NotImplementedError(msg)
        for term_idx, term in enumerate(eq.rhs_terms):
            coeff = _resolve_constant_coeff(
                term,
                coeff_eval,
                eq_idx=eq_idx,
                term_idx=term_idx,
            )
            mult = multiplier_cache[term.operator]
            # #293: route the contribution to the matrix keyed
            # by this operator's time-derivative order. ``d²_t(...)``
            # has ``_OperatorDecomp.time_order == 2`` which makes
            # ``_evolve_duhamel_per_mode`` scale the source column by
            # ``λ²`` before the Duhamel kernel.
            op_time_order = _OPERATOR_DECOMP[term.operator].time_order

            target_slot = _resolve_target_slot(term.field)
            if target_slot is not None:
                # Dynamical target: direct write into the slot column.
                M_n = _get_m_src(op_time_order)
                M_n[:, row_slot, target_slot] += eq_scale * coeff * mult
            elif term.field in constraint_idx and recovery_matrix is not None:
                # Constraint field: expand via Schur row.
                # recovery_matrix[m, c_idx, j] → contributes to column j.
                c_idx = constraint_idx[term.field]
                # (coeff * mult[m]) * recovery_matrix[m, c_idx, j]
                # → M_src_n[m, row_slot, j] += ...
                term_contrib = (
                    eq_scale * coeff * mult[:, None] * recovery_matrix[:, c_idx, :]
                )
                M_n = _get_m_src(op_time_order)
                M_n[:, row_slot, :] += term_contrib
            else:
                # Target field reference cannot be resolved — either
                # references a demoted field's velocity slot or an
                # unknown name. Record it so the driver surfaces the
                # drop on PerturbativeResult.validity (#272).
                drops.append(
                    {
                        "field": term.field,
                        "operator": term.operator,
                        "reason": "column-missing: target field has no dynamical "
                        "slot and is not a Schur-recoverable constraint (likely "
                        "a demoted field's velocity slot, v_<field>)",
                        "eq_idx": eq_idx,
                        "term_idx": term_idx,
                        "row_field": field_name,
                    },
                )
                continue

    return M_src_by_order, drops


def _evolve_duhamel_per_mode(
    eigendata: dict[str, Any],
    M_src_k_by_order: dict[int, NDArray[np.complex128]],
    t_eval: NDArray[np.float64],
    layout: StateLayout,
    grid: GridInfo,
) -> tuple[
    NDArray[np.float64],
    NDArray[np.float64],
    NDArray[np.complex128],
]:
    r"""Evaluate Pass 1 inhomogeneous solution via augmented matrix exponential.

    Solves ``dy⁽¹⁾/dt = A·y⁽¹⁾ + S·y⁰(t),  y⁽¹⁾(0) = 0``  where ``y⁰(t) =
    exp(t·A)·y⁰(0)`` is the Pass 0 solution and ``S = Σ_n M_src_n · Aⁿ``
    folds the d^n_t time-derivative orders into a single effective source
    matrix (since ``∂ⁿ_t y⁰(t) = Aⁿ · y⁰(t)`` for time-independent ``A``).

    The closed-form solution uses the **augmented matrix exponential trick**
    (Al-Mohy & Higham 2011 §5.2 "Inhomogeneous problems"):

    .. math::

        \\exp\\!\\left(t \\begin{pmatrix} A & S \\\\ 0 & A \\end{pmatrix}\\right)
        \\begin{pmatrix} 0 \\\\ y^0(0) \\end{pmatrix}
        = \\begin{pmatrix} y^{(1)}(t) \\\\ y^0(t) \\end{pmatrix}

    Per mode-block this is one ``scipy.linalg.expm`` call on a ``(2bs × 2bs)``
    matrix --- robust for arbitrary ``cond(A)``, no eigendecomposition, no
    per-eigenvalue Duhamel kernel, no V/V⁻¹ projection. Replaces the legacy
    closed-form Duhamel-kernel implementation that depended on Pass 0
    eigendata; the augmented form has backward error bounded by the function
    condition number, not the eigenvector condition number, matching the
    Pass 0 path-D guarantee (``docs/tex/modal_solver.tex`` §"Robust
    Matrix-Exponential Evolution").

    Returns ``(t, y_phys, y_hat_snap)``. ``y_phys`` is the dynamical Pass 1
    output in physical space; ``y_hat_snap`` is the Fourier-space Pass 1
    output the constraint-recovery step consumes via
    ``recovery_matrix @ y_hat_snap``.

    References
    ----------
    Al-Mohy, A.H. & Higham, N.J. (2011). "Computing the action of the matrix
    exponential." *SIAM J. Sci. Comput.* 33(2):488-511, §5.2.
    Higham, N.J. (2009). "The scaling and squaring method for the matrix
    exponential revisited." *SIAM Review* 51(4):747-764.
    """
    import logging  # noqa: PLC0415

    import scipy.linalg as sla  # type: ignore[import-untyped]  # noqa: PLC0415

    logger = logging.getLogger(__name__)

    n_slots = layout.num_slots
    n_pts = layout.num_points
    # Empty dict = no sources; still need n_modes for allocation.
    if M_src_k_by_order:
        any_mat = next(iter(M_src_k_by_order.values()))
        n_modes = any_mat.shape[0]
    else:
        rfft_last = grid.shape[-1] // 2 + 1
        n_modes = int(np.prod([*grid.shape[:-1], rfft_last]))

    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    n_snapshots = len(t_eval)
    t0 = float(t_eval[0])

    y_hat_snap = np.zeros((n_snapshots, n_slots, n_modes), dtype=np.complex128)

    # Cross-block scale: the guard below compares off-block source entries
    # against a tolerance RELATIVE to this global source magnitude.
    m_scale_total = 0.0
    for _M_n in M_src_k_by_order.values():
        if _M_n.size:
            m_scale_total = max(m_scale_total, float(np.max(np.abs(_M_n))))
    if m_scale_total == 0.0:
        m_scale_total = 1.0

    # Detect uniform spacing for snapshot precompute.
    if n_snapshots > 1:
        dts = np.diff(t_eval)
        uniform = bool(np.allclose(dts, dts[0]))
        dt_step = float(dts[0]) if uniform else None
    else:
        uniform = True
        dt_step = None

    max_order = max(M_src_k_by_order.keys()) if M_src_k_by_order else 0

    for block in eigendata["blocks"]:
        slot_indices = list(block["slot_indices"])
        M_block = block["M_block"]  # (n_modes, bs, bs) — Pass 0 evolution generator
        y0_block = block["y0_block"]  # (bs, n_modes) — Pass 0 IC in Fourier space
        idx = np.array(slot_indices)
        bs = len(slot_indices)
        n_block_modes = M_block.shape[0]

        # Cross-block coupling check: any non-zero in M_src_n's rows of this
        # block but columns outside it would mean Pass 1 mixes Pass 0 sectors
        # that were independent. Pass 0's block decomposition cannot represent
        # that — raise early. Threshold is PURELY scale-relative (#429): the
        # former absolute 1e-14 floor (`max(1e-14, rel)`) masked *physical*
        # coupling whenever max|M_src| < 1e-4 — the EH dual-Gaussian spec at
        # Phase E geometry has max|M_src| ~ 8e-15, so the floor exceeded the
        # entire source by ~10 orders while the intra-block slicing below
        # silently discarded a 3%-of-norm coupling. Schur-tail roundoff
        # (#275) scales WITH max|M_src|, so the relative form still tolerates
        # it (pinned by test_cross_block_guard_tolerates_schur_tail_noise).
        mask = np.ones(n_slots, dtype=bool)
        mask[idx] = False
        out_cols = np.flatnonzero(mask)
        # Max |entry| per out-of-block column, over modes/rows/orders — used
        # both for the verdict and to name the offending slots on refusal.
        col_max = np.zeros(out_cols.size)
        for M_n in M_src_k_by_order.values():
            if M_n.size == 0 or out_cols.size == 0:
                continue
            cross_n = np.abs(M_n[:, idx, :][:, :, out_cols])
            col_max = np.maximum(col_max, cross_n.max(axis=(0, 1)))
        cross_max = float(col_max.max()) if col_max.size else 0.0
        atol = m_scale_total * 1e-10
        if cross_max > atol:
            in_names = [layout.slots[i].name for i in slot_indices]
            offending = [
                layout.slots[int(c)].name
                for c, mag in zip(out_cols, col_max, strict=True)
                if mag > atol
            ]
            msg = (
                f"Pass 1 source matrix couples blocks that Pass 0 evolved "
                f"independently: block {in_names} is sourced by "
                f"out-of-block slot(s) {offending} "
                f"(max|cross|={cross_max:.3e} > {atol:.3e}, relative to "
                f"max|M_src|={m_scale_total:.3e}). Coupling above the "
                f"roundoff tolerance means the correction theory mixes "
                f"sectors Pass 0 evolved separately, which the per-mode "
                f"Duhamel kernel cannot represent (GH #439 tracks "
                f"cross-block support). Alternatives: run the full spec on "
                f"a time-domain scheme (--scheme cvode or --scheme ida), "
                f"or solve the full spec on plain modal with "
                f"--perturbative-order 0 (requires modal support for the "
                f"full spec's coefficients; see GH #427)."
            )
            raise NotImplementedError(msg)
        if cross_max > 0.0:
            # Sub-tolerance entries are roundoff by construction (anything
            # physical raises above); note the discard for auditability
            # since the intra-block slicing below drops them silently.
            logger.info(
                "Pass 1: discarding sub-roundoff cross-block entries for "
                "block %s (max|cross|=%.3e <= atol=%.3e)",
                [layout.slots[i].name for i in slot_indices],
                cross_max,
                atol,
            )

        # Build the effective source S = Σ_n M_src_n · Aⁿ per mode, where
        # A = M_block. For the augmented form, S absorbs the time-derivative
        # orders so we never need to evaluate `dⁿ_t y⁰(t)` separately.
        # Cache M_block^n powers up to max_order.
        if max_order > 0:
            M_block_powers: list[NDArray[np.complex128]] = [
                np.broadcast_to(
                    np.eye(bs, dtype=np.complex128),
                    (n_block_modes, bs, bs),
                ).copy(),
            ]
            for _ in range(max_order):
                # Batched matmul instead of einsum (#327).
                M_block_powers.append(M_block_powers[-1] @ M_block)
        else:
            M_block_powers = [
                np.broadcast_to(
                    np.eye(bs, dtype=np.complex128),
                    (n_block_modes, bs, bs),
                ).copy(),
            ]

        S = np.zeros((n_block_modes, bs, bs), dtype=np.complex128)
        for n, M_n in M_src_k_by_order.items():
            M_sub_n = M_n[:, idx[:, None], idx[None, :]]  # (n_modes, bs, bs)
            if n == 0:
                S += M_sub_n  # pyright: ignore[reportConstantRedefinition]
            else:
                # Batched matmul instead of einsum (#327).
                S += M_sub_n @ M_block_powers[n]  # pyright: ignore[reportConstantRedefinition]

        # Build augmented per-mode generator Aug = [[A, S], [0, A]] of size 2bs×2bs.
        Aug = np.zeros((n_block_modes, 2 * bs, 2 * bs), dtype=np.complex128)
        Aug[:, :bs, :bs] = M_block
        Aug[:, :bs, bs:] = S
        Aug[:, bs:, bs:] = M_block

        # Augmented IC: [0; y⁰(0)] per mode → shape (n_modes, 2bs).
        ic_aug = np.zeros((n_block_modes, 2 * bs), dtype=np.complex128)
        ic_aug[:, bs:] = y0_block.T

        if uniform and dt_step is not None and dt_step > 0:
            # Precompute exp(Aug · dt) once per mode, apply iteratively.
            exp_Aug_dt = np.empty_like(Aug)
            for m in range(n_block_modes):
                exp_Aug_dt[m] = sla.expm(Aug[m] * dt_step)  # pyright: ignore[reportUnknownArgumentType]
            y_aug = ic_aug.copy()
            for ti in range(n_snapshots):
                if ti > 0:
                    y_aug = np.einsum("mij,mj->mi", exp_Aug_dt, y_aug)
                # Top half (y_aug[:, :bs]) is Pass 1; bottom half is Pass 0
                # (validated by the augmented identity but unused here).
                pass1_out = y_aug[:, :bs]  # (n_modes, bs)
                for local_i, slot_idx in enumerate(slot_indices):
                    y_hat_snap[ti, slot_idx] += pass1_out[:, local_i]
        else:
            # Non-uniform spacing: compute exp(Aug · dt) per snapshot per mode.
            # Skip ti=0 (Pass 1 IC = 0 by construction).
            for ti in range(n_snapshots):
                dt = float(t_eval[ti] - t0)
                if dt == 0.0:
                    continue
                y_aug_t = np.empty_like(ic_aug)
                for m in range(n_block_modes):
                    y_aug_t[m] = sla.expm(Aug[m] * dt) @ ic_aug[m]  # pyright: ignore[reportUnknownArgumentType]
                pass1_out = y_aug_t[:, :bs]
                for local_i, slot_idx in enumerate(slot_indices):
                    y_hat_snap[ti, slot_idx] += pass1_out[:, local_i]

    # Inverse FFT each slot at each snapshot back to physical space
    y_phys = np.zeros((n_snapshots, n_slots * n_pts))
    for ti in range(n_snapshots):
        for si in range(n_slots):
            y_phys_block = np.fft.irfftn(
                y_hat_snap[ti, si].reshape(rfft_shape),
                s=grid.shape,
                axes=list(range(len(grid.shape))),
            ).ravel()
            y_phys[ti, si * n_pts : (si + 1) * n_pts] = np.real(y_phys_block)

    return t_eval.astype(np.float64), y_phys, y_hat_snap


def solve_modal_pass1(
    eigendata: dict[str, Any],
    correction_spec: EquationSystem,
    grid: GridInfo,
    t_eval: NDArray[np.float64],
    *,
    parameters: dict[str, float] | None = None,
) -> PerturbativePass1Result:
    """Solve the Pass 1 correction  dy⁽¹⁾/dt = A·y⁽¹⁾ + M_src·y⁰(t).

    Reuses ``eigendata`` from a prior ``solve_modal(..., return_eigendata=True)``
    call on the base (order-0) spec. ``correction_spec`` is the
    ``spec.filter_by_order(n)`` spec containing only the order-``n``
    RHS terms. The closed-form Duhamel kernel is applied per mode in
    the eigenbasis of the base operator — no ODE integration.

    Parameters
    ----------
    eigendata : dict
        Output of Pass 0 with ``return_eigendata=True``. Must include
        ``"blocks"``, ``"mode_k"``, ``"state_layout"``, and (for
        constraint systems) ``"schur_ops"``.
    correction_spec : EquationSystem
        The correction-only equation system (RHS terms filtered by
        ``order_in_eps``). The LHS structure is unused — only the
        right-hand side operator terms matter for the source.
    grid : GridInfo
        Same grid used in Pass 0 (must match the eigendata's
        ``mode_k``).
    t_eval : array of floats
        Evaluation times. Pass 1's IC is zero, so ``t_eval[0]`` can be
        any reference time; output at ``t = t_eval[0]`` is zero.
    parameters : dict, optional
        Runtime parameter overrides for the correction coefficients.

    Returns
    -------
    PerturbativePass1Result
        ``t``, ``y``, ``success``, ``message`` (from SolverResult) plus:

        * ``y_hat_dyn`` — Fourier-space dynamical output with shape
          ``(n_snapshots, n_slots, n_modes)``. Consumed by
          :func:`PerturbativeSolver` to recover constraint fields via
          Schur (v6 Gap C).
        * ``correction_drops`` — diagnostic records for any correction
          terms or equations that could not be routed. Aggregated by
          the driver onto ``PerturbativeResult.validity`` (#272).
    """
    from tidal.solver.coefficients import CoefficientEvaluator  # noqa: PLC0415

    # GH #421: the Pass 1 source builder consumes the correction spec's
    # kinetic coefficients without a grid (`_build_source_matrix_k` has no
    # grid parameter), so a position-dependent kinetic surviving into the
    # corrections must be refused here rather than fail deep inside
    # evaluate_coefficient.  Canonicalized flows (GH #380) never hit this:
    # their correction kinetics are constants.
    if correction_spec.has_position_dependent_kinetic():
        _raise_position_dependent_kinetic(correction_spec, "solve_modal_pass1")

    layout = eigendata["state_layout"]

    # Correction coefficient evaluator — driven by the correction spec
    # so only correction-term coefficients are pre-resolved.
    coeff_eval = CoefficientEvaluator(correction_spec, grid, parameters or {})

    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape_list = list(grid.shape)
    rfft_shape_list[-1] = grid.shape[-1] // 2 + 1
    rfft_shape = tuple(rfft_shape_list)

    M_src_k_by_order, correction_drops = _build_source_matrix_k(
        correction_spec,
        layout,
        grid,
        coeff_eval,
        k_grid,
        rfft_shape,
        schur_ops=eigendata.get("schur_ops"),
    )

    if correction_drops:
        import logging  # noqa: PLC0415

        logger = logging.getLogger(__name__)
        for d in correction_drops:
            logger.warning(
                "Pass 1 correction drop: field=%r operator=%r reason=%s "
                "(eq_idx=%s, term_idx=%s). See #272.",
                d.get("field"),
                d.get("operator"),
                d.get("reason"),
                d.get("eq_idx"),
                d.get("term_idx"),
            )

    times, y_phys, y_hat_dyn = _evolve_duhamel_per_mode(
        eigendata,
        M_src_k_by_order,
        t_eval,
        layout,
        grid,
    )

    result: PerturbativePass1Result = {
        "t": times,
        "y": y_phys,
        "success": True,
        "message": "Pass 1 closed-form Duhamel",
        "y_hat_dyn": y_hat_dyn,
        "correction_drops": correction_drops,
    }
    return result
