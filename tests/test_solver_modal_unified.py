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
    StateLayout, GridInfo, CoefficientEvaluator, list[np.ndarray], tuple[int, ...]
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

    A_rhs, B_lhs, _, _, _, _ = _build_evolution_matrices(
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
        import scipy.linalg as sla

        w = cast(
            "np.ndarray[Any, Any]",
            sla.eig(A_rhs[m_idx], B_lhs[m_idx], right=False),
        )
    else:
        w = np.linalg.eigvals(A_rhs[m_idx])

    w_finite = w[np.isfinite(w)]
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
    import scipy.linalg as sla

    spec = _make_rank_deficient_spec()
    layout, grid, coeff_eval, k_grid, rfft_shape = _build_eval_context(spec)

    A_rhs, B_lhs, _, _, _, _ = _build_evolution_matrices(
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
