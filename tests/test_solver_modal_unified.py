"""Regression tests for the unified modal-solver evolution-matrix builder (#256).

Pins the failure mode and the correct behaviour before and after the
unification of `build_constraint_eliminated_matrices` and
`_build_evolution_matrices` into one path.

The core failure: a theory with algebraic constraint fields **and**
cross-`d2_t` couplings on dynamical equations (e.g. the TorsionCDT-trace
dark photon) routes through `build_constraint_eliminated_matrices` which
silently mis-classifies `d2_t(other_field)` RHS terms as spatial operators
(because it uses `_EXACT_MULTIPLIERS[op]` which returns only the spatial
part, `ones_like(k) = 1` for `d2_t`). The resulting evolution matrix has
an implicit diagonal mass matrix instead of the actual rank-deficient
pattern, producing spurious tachyons.

The correct path is `_build_evolution_matrices`, which
classifies terms by `_OPERATOR_DECOMP[op].time_order` and handles
rank-deficient mass matrices via per-mode eigendecomposition.

These tests use a minimal synthetic 3-field spec that exhibits the
failure mode: three dynamical fields `t_0, t_1, t_2` coupled via
cross-`d2_t` terms in a trace-like pattern whose mass matrix is
`[[1,1,1],[1,1,1],[1,1,1]]` (rank 1, 2 null directions).
"""

# ruff: noqa: RUF002, ERA001 — math symbols in docstrings; commented stubs.

from __future__ import annotations

from typing import Any, cast

import numpy as np

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_evolution_matrices,  # type: ignore[attr-defined]
    _build_k_axes,  # type: ignore[attr-defined]
    _build_k_grid,  # type: ignore[attr-defined]
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Minimal synthetic spec: rank-deficient mass matrix
# ---------------------------------------------------------------------------
#
# 3 dynamical fields (t_0, t_1, t_2), each with kinetic -1 on the LHS
# after normalization, and cross-d2_t RHS couplings of coefficient -1.
# After normalization (dividing each row by its LHS -1):
#     d2_t(t_i) + d2_t(t_{i+1}) + d2_t(t_{i+2}) = k^2 * t_i   (laplacian)
# Mass matrix is thus [[1,1,1],[1,1,1],[1,1,1]] — rank 1, eigenvalues (3,0,0).
# The generalized path must detect and eliminate the 2-dim null subspace.
# The constraint path (if exercised) will silently mis-classify the d2_t
# terms and produce wrong dynamics.

_RANK_DEFICIENT_MASS_SPEC: dict[str, object] = {
    "metadata": {"source": "inline-test"},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "t_0", "index": 0, "is_dynamical": True},
        {"name": "t_1", "index": 1, "is_dynamical": True},
        {"name": "t_2", "index": 2, "is_dynamical": True},
    ],
    "equations": [
        {
            "field": f"t_{i}",
            "lhs": {
                "expression": f"d2_t(t_{i})",
                "order": {"time": 2, "space": 0},
                "kinetic_coefficient_symbolic": "-1",
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    # Cross-d2_t: trace structure
                    *[
                        {
                            "coefficient": 1.0,
                            "operator": "d2_t",
                            "field": f"t_{j}",
                        }
                        for j in range(3)
                        if j != i
                    ],
                    # Laplacian on own field (gives stable propagation on the physical mode)
                    {
                        "coefficient": -1.0,
                        "operator": "laplacian_x",
                        "field": f"t_{i}",
                    },
                ],
            },
        }
        for i in range(3)
    ],
}


def _make_rank_deficient_spec() -> EquationSystem:
    """Construct the synthetic rank-deficient-mass spec."""
    from tidal.symbolic.json_loader import normalize_kinetic_coefficients

    spec = cast(
        "EquationSystem",
        EquationSystem.from_dict(_RANK_DEFICIENT_MASS_SPEC),  # type: ignore[arg-type]
    )
    return normalize_kinetic_coefficients(spec, {})


def _build_eval_context(
    spec: EquationSystem,
) -> tuple[
    StateLayout, GridInfo, CoefficientEvaluator, list[np.ndarray], tuple[int, ...],
]:
    grid = GridInfo(shape=(32,), bounds=((0.0, 10.0),), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, {})
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape = (grid.shape[0] // 2 + 1,)
    return layout, grid, coeff_eval, k_grid, rfft_shape


# ---------------------------------------------------------------------------
# Historical note (test 1 removed): previously this file contained
# `test_constraint_path_misclassifies_d2_t_cross_couplings`, which pinned
# the #256 silent bug in `build_constraint_eliminated_matrices` — namely
# that it used `_EXACT_MULTIPLIERS` (spatial part only) and silently
# reinterpreted `d2_t(t_j)` RHS terms as `identity(field_j)` stiffness
# couplings. That function has now been deleted from the solver, so the
# test is no longer reachable. The unified path tests below replace it.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test 1: unified builder correctly eliminates the mass-null subspace
# ---------------------------------------------------------------------------


def test_unified_builder_handles_rank_deficient_mass() -> None:
    """The generalized path detects rank-deficient M and eliminates null subspace.

    For M = [[1,1,1],[1,1,1],[1,1,1]], eigenvalues are (3, 0, 0). The
    generalized path should:
      1. Detect `has_singular_M` via `det(M) ~ 0`.
      2. Eigendecompose M, find 2 null directions.
      3. Schur-eliminate them into K_eff.
      4. Return an evolution matrix whose eigenvalues are finite and
         physical (specifically: imaginary for the surviving trace
         direction, zero for the frozen null directions).
    """
    spec = _make_rank_deficient_spec()
    layout, grid, coeff_eval, k_grid, rfft_shape = _build_eval_context(spec)

    A_rhs, B_lhs, _, _, _, _, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        coeff_eval,
        k_grid,
        rfft_shape,
    )

    # Check eigenvalues at a non-DC mode
    m_idx = 2
    if B_lhs is not None:
        import scipy.linalg as sla  # type: ignore[import-untyped]

        w = cast(
            "np.ndarray[Any, np.dtype[np.complex128]]",
            sla.eig(A_rhs[m_idx], B_lhs[m_idx], right=False),
        )
    else:
        w = cast(
            "np.ndarray[Any, np.dtype[np.complex128]]",
            np.linalg.eigvals(A_rhs[m_idx]),
        )

    w_finite: np.ndarray[Any, np.dtype[np.complex128]] = w[np.isfinite(w)]
    # After mass-matrix Schur elimination, no eigenvalue should have a
    # large positive real part (spurious tachyon). The physical
    # trace-direction mode should have pure imaginary eigenvalue ±ik,
    # and the null directions should be frozen (eigenvalue 0).
    max_re = float(np.max(np.abs(w_finite.real)))
    assert max_re < 1e-6, (
        f"Generalized path produced a tachyonic eigenvalue "
        f"(max |Re(λ)| = {max_re:.3e}) for a rank-deficient mass "
        f"matrix that should have been Schur-eliminated. "
        f"Eigenvalues: {sorted(w_finite, key=lambda z: -abs(z.imag))[:6]}"
    )


# ---------------------------------------------------------------------------
# Test 3: unified builder produces stable generalized eigenvalues across all modes
# ---------------------------------------------------------------------------


def test_unified_builder_generalized_eig_tachyon_free() -> None:
    """The unified builder's (A_rhs, B_lhs) pair yields tachyon-free generalized eig.

    The builder returns A_rhs and B_lhs separately so downstream can use
    scipy.linalg.eig(A, B) (QZ decomposition), which handles rank-deficient
    mass matrices via the Schur form.  This test verifies that every mode's
    generalized eigenvalue spectrum is tachyon-free (max Re(lambda) ~ 0)
    for the synthetic rank-deficient spec.
    """
    import scipy.linalg as sla  # type: ignore[import-untyped]

    spec = _make_rank_deficient_spec()
    layout, grid, coeff_eval, k_grid, rfft_shape = _build_eval_context(spec)

    A_rhs, B_lhs, _, _, _, _, _, _ = _build_evolution_matrices(
        spec,
        layout,
        grid,
        coeff_eval,
        k_grid,
        rfft_shape,
    )

    # Scan every mode for spurious tachyons (positive Re part) in the
    # generalized-eigenvalue spectrum.  Infinite/large eigenvalues are
    # filtered out — they correspond to gauge/null directions that the
    # downstream modal evolution freezes at IC.
    max_re_across_modes = 0.0
    for m in range(A_rhs.shape[0]):
        if B_lhs is not None:
            w = cast(
                "np.ndarray[Any, Any]",
                sla.eig(A_rhs[m], B_lhs[m], right=False),
            )
        else:
            w = cast(
                "np.ndarray[Any, Any]",
                np.linalg.eigvals(A_rhs[m]),
            )
        # Filter |λ| > 1e12 (gauge modes, zeroed downstream) and infinite
        w_phys = w[np.isfinite(w) & (np.abs(w) < 1e12)]
        if len(w_phys) > 0:
            max_re_across_modes = max(max_re_across_modes, float(np.max(w_phys.real)))

    assert max_re_across_modes < 1e-6, (
        f"Unified builder produced a tachyonic eigenvalue "
        f"(max Re(λ) = {max_re_across_modes:.3e}) across all modes. "
        f"The generalized eig(A, B) dispatch must handle rank-deficient M."
    )


# ---------------------------------------------------------------------------
# Test 3: SVD Schur back-projection continuity through critical point (#260)
# ---------------------------------------------------------------------------


def _make_asymmetric_m_spec(deltam: float, xi: float = 1.0) -> EquationSystem:
    """Build a 2-field spec with asymmetric mass matrix.

    Models the FV dark-photon kinetic-mixing system:
      d2_t(a) + 2*deltam * d2_t(t) = -k^2 * a     (photon)
      -2*deltam * d2_t(a) + (-xi) * d2_t(t) = -k^2 * t   (dark photon)

    After normalizing by LHS kinetic coefficient:
      M = [[1, 2*deltam], [-2*deltam/xi, 1]]
    (pre-normalization: M = [[1, 2*deltam], [+2*deltam, -xi]])

    det(M) = -xi + 4*deltam^2.  Critical point: deltam = sqrt(xi)/2.
    """
    from tidal.symbolic.json_loader import normalize_kinetic_coefficients

    spec_dict: dict[str, object] = {
        "metadata": {"source": "inline-test"},
        "spacetime": {
            "dimension": 2,
            "signature": [-1, 1],
            "coordinates": ["t", "x"],
        },
        "fields": [
            {"name": "a_1", "index": 0, "is_dynamical": True},
            {"name": "t_1", "index": 1, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "a_1",
                "lhs": {
                    "expression": "d2_t(a_1)",
                    "order": {"time": 2, "space": 0},
                    "kinetic_coefficient_symbolic": "-1",
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 2 * deltam,
                            "operator": "d2_t",
                            "field": "t_1",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "laplacian_x",
                            "field": "a_1",
                        },
                    ],
                },
            },
            {
                "field": "t_1",
                "lhs": {
                    "expression": "d2_t(t_1)",
                    "order": {"time": 2, "space": 0},
                    "kinetic_coefficient_symbolic": str(-xi),
                },
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -2 * deltam,
                            "operator": "d2_t",
                            "field": "a_1",
                        },
                        {
                            "coefficient": -1.0,
                            "operator": "laplacian_x",
                            "field": "t_1",
                        },
                    ],
                },
            },
        ],
    }
    spec = cast(
        "EquationSystem",
        EquationSystem.from_dict(spec_dict),  # type: ignore[arg-type]
    )
    return normalize_kinetic_coefficients(spec, {})


def test_svd_schur_continuity_through_critical_point() -> None:
    """Eigenvalues of A are continuous through the mass-matrix critical point.

    At deltam = sqrt(xi)/2, det(M) = 0 and the solver switches from
    M^{-1}*K (invertible) to SVD Schur elimination (rank-deficient).
    The physical eigenvalues must be continuous across this boundary.

    Regression test for #260: previously the SVD Schur back-projection
    dropped the constraint recovery term, producing discontinuous eigenvalues.
    """
    xi = 1.0
    dm_crit = np.sqrt(xi) / 2.0  # = 0.5
    dm_near = dm_crit - 0.01  # 0.49, invertible M

    spec_near = _make_asymmetric_m_spec(dm_near, xi)
    spec_crit = _make_asymmetric_m_spec(dm_crit, xi)

    layout_n, grid_n, ce_n, kg_n, rs_n = _build_eval_context(spec_near)
    layout_c, grid_c, ce_c, kg_c, rs_c = _build_eval_context(spec_crit)

    A_near, _, _, _, _, _, _, _ = _build_evolution_matrices(
        spec_near, layout_n, grid_n, ce_n, kg_n, rs_n,
    )
    A_crit, _, _, _, _, _, _, _ = _build_evolution_matrices(
        spec_crit, layout_c, grid_c, ce_c, kg_c, rs_c,
    )

    # Compare eigenvalues at a representative non-DC mode
    m_idx = 2
    w_near = np.asarray(
        np.sort_complex(np.linalg.eigvals(A_near[m_idx])),  # type: ignore[reportUnknownArgumentType]
        dtype=np.complex128,
    )
    w_crit = np.asarray(
        np.sort_complex(np.linalg.eigvals(A_crit[m_idx])),  # type: ignore[reportUnknownArgumentType]
        dtype=np.complex128,
    )

    # At the critical point, one DOF is eliminated, so A_crit may have
    # some zero eigenvalues.  Compare the non-zero eigenvalues.
    mag_near = np.abs(w_near)
    mag_crit = np.abs(w_crit)

    # Physical modes: eigenvalues with |λ| > 0.01
    phys_near = sorted(w_near[mag_near > 0.01], key=lambda z: -abs(z.imag))
    phys_crit = sorted(w_crit[mag_crit > 0.01], key=lambda z: -abs(z.imag))

    # Both should have the same number of physical modes (the critical
    # point eliminates one field DOF = 2 eigenvalues from the first-order
    # system, but the remaining physical modes should persist).
    assert len(phys_crit) >= 2, (
        f"Critical point should retain at least one physical oscillation mode "
        f"(2 eigenvalues ±iω), got {len(phys_crit)} non-zero eigenvalues. "
        f"All eigenvalues: {w_crit}"
    )

    # The largest physical eigenvalue magnitudes should be close
    # (continuity check).  Allow 20% relative tolerance for the
    # discrete dm step.
    max_imag_near = max(abs(z.imag) for z in phys_near)
    max_imag_crit = max(abs(z.imag) for z in phys_crit)
    rel_diff = abs(max_imag_near - max_imag_crit) / max(max_imag_near, 1e-15)
    assert rel_diff < 0.2, (
        f"Eigenvalue discontinuity at critical point: "
        f"|ω_near| = {max_imag_near:.6f}, |ω_crit| = {max_imag_crit:.6f}, "
        f"relative diff = {rel_diff:.3f}. "
        f"SVD Schur back-projection may be incorrect (#260)."
    )


# ---------------------------------------------------------------------------
# Test 4: rank-deficient K_cc handled by pseudoinverse (#260, CDT case)
# ---------------------------------------------------------------------------


def _make_cdt_like_spec(deltam: float, xi: float = 1.0) -> EquationSystem:
    """Build a 4-field spec mimicking the CDT dark-photon rank structure.

    Four fields: a (photon), t0, t1, t2 (torsion trace components).
    Each torsion field has d2_t cross-couplings to a and to each other,
    producing a 4×4 mass matrix with rank that drops at deltam = √ξ/2.

    At generic deltam, M has rank 3 (1 null direction from the trace
    structure).  At deltam = √ξ/2, M drops to rank 2 (additional null
    direction from the kinetic-mixing critical point), making K_cc
    expand from 1×1 to 2×2 and potentially rank-deficient.

    This tests the np.linalg.pinv(K_cc) path that handles rank-deficient
    constraint stiffness blocks.
    """
    from tidal.symbolic.json_loader import normalize_kinetic_coefficients

    # Trace coupling: each torsion field has d2_t to both other torsion fields
    # and to the photon.  M structure:
    #   M[a, t_i] = -2*deltam   (kinetic mixing)
    #   M[t_i, a] = +2*deltam   (antisymmetric)
    #   M[t_i, t_j] = -xi       (trace structure, i≠j)
    #   M[a, a] = 1, M[t_i, t_i] = -xi
    fields = [
        {"name": "a", "index": 0, "is_dynamical": True},
        {"name": "t0", "index": 1, "is_dynamical": True},
        {"name": "t1", "index": 2, "is_dynamical": True},
        {"name": "t2", "index": 3, "is_dynamical": True},
    ]

    equations = []

    # Photon equation: d2_t(a) + sum_i 2*deltam * d2_t(t_i) = -k^2 * a
    equations.append(
        {
            "field": "a",
            "lhs": {
                "expression": "d2_t(a)",
                "order": {"time": 2, "space": 0},
                "kinetic_coefficient_symbolic": "-1",
            },
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    *[
                        {
                            "coefficient": 2 * deltam,
                            "operator": "d2_t",
                            "field": f"t{i}",
                        }
                        for i in range(3)
                    ],
                    {
                        "coefficient": -1.0,
                        "operator": "laplacian_x",
                        "field": "a",
                    },
                ],
            },
        },
    )

    # Torsion equations: (-xi)*d2_t(t_i) - 2*deltam*d2_t(a)
    #   + sum_{j≠i} (-xi)*d2_t(t_j) = -k^2 * t_i
    for i in range(3):
        terms = [
            {
                "coefficient": -2 * deltam,
                "operator": "d2_t",
                "field": "a",
            },
            *[
                {
                    "coefficient": -xi,
                    "operator": "d2_t",
                    "field": f"t{j}",
                }
                for j in range(3)
                if j != i
            ],
            {
                "coefficient": -1.0,
                "operator": "laplacian_x",
                "field": f"t{i}",
            },
        ]
        equations.append(
            {
                "field": f"t{i}",
                "lhs": {
                    "expression": f"d2_t(t{i})",
                    "order": {"time": 2, "space": 0},
                    "kinetic_coefficient_symbolic": str(-xi),
                },
                "rhs": {"type": "linear_combination", "terms": terms},
            },
        )

    spec_dict: dict[str, object] = {
        "metadata": {"source": "inline-test"},
        "spacetime": {
            "dimension": 2,
            "signature": [-1, 1],
            "coordinates": ["t", "x"],
        },
        "fields": fields,
        "equations": equations,
    }
    spec = cast(
        "EquationSystem",
        EquationSystem.from_dict(spec_dict),  # type: ignore[arg-type]
    )
    return normalize_kinetic_coefficients(spec, {})


def test_kcc_rank_deficient_pinv_continuity() -> None:
    """K_cc pseudoinverse maintains continuity when K_cc becomes rank-deficient.

    The CDT dark-photon model's mass matrix has a null space that expands
    at the kinetic-mixing critical point (deltam = sqrt(xi)/2).  At this
    point, the rotated stiffness constraint block K_cc becomes rank-deficient.
    The pseudoinverse (np.linalg.pinv) must correctly handle undetermined
    constraint modes to produce continuous eigenvalues.

    Regression test for #260 (CDT case): previously, 1e-14*I regularization
    of K_cc corrupted the constraint recovery for undetermined modes.
    """
    xi = 1.0
    dm_crit = np.sqrt(xi) / 2.0
    dm_near = dm_crit - 0.01

    spec_near = _make_cdt_like_spec(dm_near, xi)
    spec_crit = _make_cdt_like_spec(dm_crit, xi)

    layout_n, grid_n, ce_n, kg_n, rs_n = _build_eval_context(spec_near)
    layout_c, grid_c, ce_c, kg_c, rs_c = _build_eval_context(spec_crit)

    A_near, _, _, _, _, _, _, _ = _build_evolution_matrices(
        spec_near, layout_n, grid_n, ce_n, kg_n, rs_n,
    )
    A_crit, _, _, _, _, _, _, _ = _build_evolution_matrices(
        spec_crit, layout_c, grid_c, ce_c, kg_c, rs_c,
    )

    m_idx = 2
    w_near = np.asarray(np.linalg.eigvals(A_near[m_idx]), dtype=np.complex128)
    w_crit = np.asarray(np.linalg.eigvals(A_crit[m_idx]), dtype=np.complex128)

    # No tachyonic eigenvalues at either parameter point
    max_re_near = float(np.max(np.abs(w_near.real)))
    max_re_crit = float(np.max(np.abs(w_crit.real)))
    assert max_re_near < 1.0, (
        f"Near-critical point has tachyonic eigenvalue: max |Re(λ)| = {max_re_near:.3e}"
    )
    assert max_re_crit < 1.0, (
        f"Critical point has tachyonic eigenvalue: max |Re(λ)| = {max_re_crit:.3e}. "
        f"K_cc pseudoinverse may not be handling rank-deficient constraint stiffness."
    )

    # Physical eigenvalues (|λ| > 0.01) should be continuous
    phys_near: list[complex] = sorted(
        w_near[np.abs(w_near) > 0.01], key=lambda z: -abs(z.imag),
    )
    phys_crit: list[complex] = sorted(
        w_crit[np.abs(w_crit) > 0.01], key=lambda z: -abs(z.imag),
    )

    if phys_near and phys_crit:
        max_imag_near = max(abs(z.imag) for z in phys_near)
        max_imag_crit = max(abs(z.imag) for z in phys_crit)
        rel_diff = abs(max_imag_near - max_imag_crit) / max(max_imag_near, 1e-15)
        assert rel_diff < 0.5, (
            f"Eigenvalue discontinuity at K_cc rank-deficient critical point: "
            f"|ω_near| = {max_imag_near:.6f}, |ω_crit| = {max_imag_crit:.6f}, "
            f"relative diff = {rel_diff:.3f}. "
            f"K_cc pseudoinverse may be incorrect (#260, CDT case)."
        )
