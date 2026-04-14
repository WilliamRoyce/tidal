"""Regression tests for the unified modal-solver evolution-matrix builder (#256).

Pins the failure mode and the correct behaviour before and after the
unification of `build_constraint_eliminated_matrices` and
`_build_generalized_evolution_matrices` into one path.

The core failure: a theory with algebraic constraint fields **and**
cross-`d2_t` couplings on dynamical equations (e.g. the TorsionCDT-trace
dark photon) routes through `build_constraint_eliminated_matrices` which
silently mis-classifies `d2_t(other_field)` RHS terms as spatial operators
(because it uses `_EXACT_MULTIPLIERS[op]` which returns only the spatial
part, `ones_like(k) = 1` for `d2_t`). The resulting evolution matrix has
an implicit diagonal mass matrix instead of the actual rank-deficient
pattern, producing spurious tachyons.

The correct path is `_build_generalized_evolution_matrices`, which
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
    _build_generalized_evolution_matrices,  # type: ignore[attr-defined]
    _build_k_axes,  # type: ignore[attr-defined]
    _build_k_grid,  # type: ignore[attr-defined]
    build_constraint_eliminated_matrices,
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
# Test 1: constraint path misclassifies d2_t as identity
# ---------------------------------------------------------------------------


def test_constraint_path_misclassifies_d2_t_cross_couplings() -> None:
    """Pins the #256 bug: constraint path treats d2_t as identity on the field slot.

    After normalization, each equation's RHS has two `d2_t(t_j)` terms
    with coefficient `-1` (from +1 moved to LHS, divided by -1 LHS kin).
    The constraint path's iteration at modal.py:521 does:

        dj = dyn_slot_map[term.field]
        A_dd[:, vel_slot, dj] += coeff * mult

    With `mult = _EXACT_MULTIPLIERS["d2_t"](k) = ones_like(k) = 1`, this
    adds `-1 * 1 = -1` to `A_dd[vel_i, field_j]` — i.e., treats `d2_t(t_j)`
    as if it were `identity(t_j)` and puts it in the stiffness block.

    The correct treatment would put `-1` into a mass matrix `M[i, j]`,
    not the stiffness, and then invert M to get the true evolution.
    """
    spec = _make_rank_deficient_spec()
    layout, grid, coeff_eval, k_grid, rfft_shape = _build_eval_context(spec)

    A_reduced, _, _, _, _ = build_constraint_eliminated_matrices(
        spec,
        layout,
        grid,
        coeff_eval,
        k_grid,
        rfft_shape,
    )

    # A_reduced is (n_modes, 6, 6) for 3 fields × 2 (field+velocity).
    # The constraint path's mis-classification shows up as:
    #   A_dd[vel_i, field_j] = -1  for j != i   (silently wrong)
    # Whereas a correct handling would have A_reduced dependent on
    # M_inv · K, which for M = [[1,1,1],[1,1,1],[1,1,1]] and K = diag(-k^2)
    # would NOT produce these cross-field entries.
    #
    # Specifically: the velocity rows of A_reduced should couple to
    # field_j only through a mass-inverted stiffness matrix. Let's check
    # mode 2 (k != 0):
    m_idx = 2  # first non-DC mode

    # Find the velocity slot indices
    vel_slots = sorted(layout.velocity_slot_map.values())
    field_slots = sorted(layout.field_slot_map.values())
    assert len(vel_slots) == 3, f"expected 3 velocity slots, got {len(vel_slots)}"
    assert len(field_slots) == 3, f"expected 3 field slots, got {len(field_slots)}"

    # The (vel_0, field_1) entry carries the cross-d2_t coupling.
    # In the silently-wrong constraint path, this should be exactly -1
    # (the normalized coefficient of the d2_t RHS term treated as identity).
    cross = A_reduced[m_idx, vel_slots[0], field_slots[1]]
    # The silent bug produces exactly -1 at this entry. A correct
    # mass-inverted handling would produce a k^2-dependent value.
    assert np.isclose(cross.real, -1.0, atol=1e-12), (
        f"Expected silently-misclassified cross-coupling = -1 "
        f"(the #256 bug), got {cross.real:.6f}. "
        f"If this test fails, either the bug has been fixed (good — "
        f"delete this xfail) or the constraint path's iteration "
        f"behaviour has changed in an unexpected way."
    )


# ---------------------------------------------------------------------------
# Test 2: generalized path correctly eliminates the mass-null subspace
# ---------------------------------------------------------------------------


def test_generalized_path_handles_rank_deficient_mass() -> None:
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

    A_rhs, B_lhs, _, _, _, _ = _build_generalized_evolution_matrices(
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
# Test 3: pre-solved evolution matrix is stable across all modes
# ---------------------------------------------------------------------------


def test_presolved_unified_matrix_is_tachyon_free() -> None:
    """The unified builder (D.2) will pre-solve B_lhs · A = A_final.

    Verify that the pre-solved form, which is what D.2 will return from
    the unified builder, produces a first-order evolution matrix whose
    eigenvalues are tachyon-free across every mode. This locks in the
    behaviour that D.2 will make the default across all specs.
    """
    spec = _make_rank_deficient_spec()
    layout, grid, coeff_eval, k_grid, rfft_shape = _build_eval_context(spec)

    A_rhs, B_lhs, _, _, _, _ = _build_generalized_evolution_matrices(
        spec,
        layout,
        grid,
        coeff_eval,
        k_grid,
        rfft_shape,
    )

    # Pre-solve B_lhs · d' = A_rhs · d, mirroring D.2.
    if B_lhs is not None:
        A_final = np.zeros_like(A_rhs)
        for m in range(A_rhs.shape[0]):
            try:
                A_final[m] = np.linalg.solve(B_lhs[m], A_rhs[m])
            except np.linalg.LinAlgError:
                A_final[m] = cast(
                    "np.ndarray[Any, Any]",
                    np.linalg.lstsq(B_lhs[m], A_rhs[m], rcond=None)[0],
                )
    else:
        A_final = A_rhs

    # Scan every mode for spurious tachyons (positive Re part).
    max_re_across_modes = 0.0
    for m in range(A_final.shape[0]):
        w = cast(
            "np.ndarray[Any, Any]",
            np.linalg.eigvals(A_final[m]),
        )
        w_finite = w[np.isfinite(w)]
        if len(w_finite) > 0:
            max_re_across_modes = max(max_re_across_modes, float(np.max(w_finite.real)))

    assert max_re_across_modes < 1e-6, (
        f"Pre-solved unified matrix has a tachyonic eigenvalue "
        f"(max Re(λ) = {max_re_across_modes:.3e}) across all modes. "
        f"The unified builder (D.2) must produce tachyon-free dynamics."
    )
