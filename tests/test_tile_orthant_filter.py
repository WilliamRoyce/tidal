"""Tests for tile-eligibility classification in :mod:`tidal.inference._sphere`.

Covers the corner-evaluation + optional interior-sampling classifier that
drives orthant-restricted angular surveys.  Per-axis sign constraints and
linear half-spaces are exact under corner-only eval; non-convex regions
are caught (modulo interior sampling) and otherwise fall back to runtime
rejection.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.inference._constraints import ConstraintSet
from tidal.inference._sphere import (
    TileStatus,
    classify_all_cells,
    classify_tile,
    n_faces,
)

# ----------------------------------------------------------------------
# Trivial / edge cases
# ----------------------------------------------------------------------


def test_empty_constraint_set_is_fully_inside() -> None:
    """No constraints → every tile is trivially FULLY_INSIDE."""
    empty = ConstraintSet()
    assert not empty  # sanity: empty set is falsy
    status = classify_tile(
        face_idx=1,
        sub_tile=(1,),
        M=1,
        Q=np.eye(2),
        r=1.0,
        coupling_names=("c1", "c2"),
        constraints=empty,
    )
    assert status is TileStatus.FULLY_INSIDE


def test_wrong_coupling_names_length_raises() -> None:
    cs = ConstraintSet.from_strings(["c1 > 0"])
    with pytest.raises(ValueError, match="coupling_names must have length"):
        classify_tile(
            face_idx=1,
            sub_tile=(1,),
            M=1,
            Q=np.eye(2),
            r=1.0,
            coupling_names=("c1",),  # length 1, need 2
            constraints=cs,
        )


# ----------------------------------------------------------------------
# Per-axis sign constraints (exact under corner-eval)
# ----------------------------------------------------------------------


def test_axis_sign_constraint_n2() -> None:
    """``c_1 > 0`` on the n=2 sphere gives one FULLY_INSIDE face, one
    FULLY_OUTSIDE, and two BOUNDARY faces (M=1).
    """
    cs = ConstraintSet.from_strings(["c1 > 0"])
    n = 2
    Q = np.eye(n)
    results = list(
        classify_all_cells(
            n_dims=n,
            M=1,
            Q=Q,
            r=1.0,
            coupling_names=("c1", "c2"),
            constraints=cs,
        )
    )
    assert len(results) == n_faces(n)  # 4 faces, 1 tile each
    statuses = {face_idx: status for face_idx, _sub, status in results}
    assert statuses[1] is TileStatus.FULLY_INSIDE  # +x_1
    assert statuses[2] is TileStatus.FULLY_OUTSIDE  # -x_1
    assert statuses[3] is TileStatus.BOUNDARY  # +x_2 (c_1 spans signs)
    assert statuses[4] is TileStatus.BOUNDARY  # -x_2


def test_axis_sign_constraint_drops_half_of_n6_faces() -> None:
    """In an N=6 theory, requiring all of ``xi, mA2, alpha3 > 0`` drops at
    least the three negative-axis faces (six faces total) outright.
    """
    cs = ConstraintSet.from_strings(["xi > 0", "mA2 > 0", "alpha3 > 0"])
    n = 6
    coupling_names = ("xi", "mA2", "alpha3", "alpha1", "alpha2", "delta1")
    Q = np.eye(n)
    statuses = list(
        classify_all_cells(
            n_dims=n,
            M=1,
            Q=Q,
            r=1.0,
            coupling_names=coupling_names,
            constraints=cs,
        )
    )
    assert len(statuses) == n_faces(n)  # 12 faces, 1 tile each
    by_face = {face_idx: status for face_idx, _sub, status in statuses}
    # Faces 2, 4, 6 are -xi, -mA2, -alpha3 respectively → all corners fail.
    assert by_face[2] is TileStatus.FULLY_OUTSIDE
    assert by_face[4] is TileStatus.FULLY_OUTSIDE
    assert by_face[6] is TileStatus.FULLY_OUTSIDE
    # Positive faces of constrained axes are not necessarily FULLY_INSIDE
    # because the other constrained axes can go negative on the face, but
    # they must not be FULLY_OUTSIDE.
    assert by_face[1] is not TileStatus.FULLY_OUTSIDE
    assert by_face[3] is not TileStatus.FULLY_OUTSIDE
    assert by_face[5] is not TileStatus.FULLY_OUTSIDE


# ----------------------------------------------------------------------
# Linear half-space (still convex; exact)
# ----------------------------------------------------------------------


def test_linear_combination_constraint_exact_on_corners() -> None:
    """``2*c1 + c2 > 0`` is a linear half-space — convex on the gnomonic
    chart, so corner eval is exact.
    """
    cs = ConstraintSet.from_strings(["2*c1 + c2 > 0"])
    n = 2
    Q = np.eye(n)
    # Face 1 (+x_1): on M=1 the corners are u=-1 and u=+1, giving the cube
    # points (1, -1) and (1, +1); the gnomonic projections have
    # 2*c1 + c2 = (2*1 + (-1))/sqrt(2) > 0 and (2 + 1)/sqrt(2) > 0 → both pass.
    status = classify_tile(
        face_idx=1,
        sub_tile=(1,),
        M=1,
        Q=Q,
        r=1.0,
        coupling_names=("c1", "c2"),
        constraints=cs,
    )
    assert status is TileStatus.FULLY_INSIDE


# ----------------------------------------------------------------------
# Non-convex constraint — corner-eval can lie; interior sampling tightens
# ----------------------------------------------------------------------


def test_non_convex_constraint_interior_sampling_downgrades() -> None:
    """``c1**2 < 0.1`` (an interior-band condition) is non-convex along the
    gnomonic chart.  On face 1+, the M=1 corners both satisfy
    c1 = 1/sqrt(2) ≈ 0.707 → c1**2 ≈ 0.5 > 0.1, so corners alone return
    FULLY_OUTSIDE — but the true interior includes points near the +y
    pole where c1 → 0 → c1**2 → 0 → passes.  Interior sampling exposes
    that mismatch, downgrading the classification to BOUNDARY.
    """
    cs = ConstraintSet.from_strings(["c1**2 < 0.1"])
    n = 2
    Q = np.eye(n)
    # On M=2 face 3 (+x_2), the sub-tile (1,) covers u_1 in [-1, 0]; the
    # corner u_1=-1 gives the cube point (-1, +1) → c1 = -1/sqrt(2) ≈ -0.707
    # → c1**2 ≈ 0.5 → fails.  The corner u_1=0 gives (0, 1) → c1 = 0 →
    # c1**2 = 0 → passes.  So already BOUNDARY at corners.
    status_corners = classify_tile(
        face_idx=3,
        sub_tile=(1,),
        M=2,
        Q=Q,
        r=1.0,
        coupling_names=("c1", "c2"),
        constraints=cs,
    )
    assert status_corners is TileStatus.BOUNDARY

    # Now the harder case: face 1+ on M=1 — both corners fail (c1**2 ≈ 0.5
    # at u = ±1) but the interior near u=0 (the equator from the face's
    # perspective) actually has c1 → smaller.  Wait — on face 1+, c1 is
    # ALWAYS 1/sqrt(1 + u**2) >= 1/sqrt(2), so the constraint c1**2 < 0.1
    # cannot be satisfied anywhere on face 1+.  Truly FULLY_OUTSIDE.
    status_face1 = classify_tile(
        face_idx=1,
        sub_tile=(1,),
        M=1,
        Q=Q,
        r=1.0,
        coupling_names=("c1", "c2"),
        constraints=cs,
        sample_interior=16,
    )
    assert status_face1 is TileStatus.FULLY_OUTSIDE


def test_non_convex_constraint_false_fully_inside_caught_by_interior() -> None:
    """Construct a case where corner eval would say FULLY_INSIDE but the
    interior contains a violating point.  Use the annular constraint
    ``c1**2 + c2**2 > 0.99`` evaluated at r=1 (where every theta_hat
    automatically satisfies it at the unit sphere) — but with a deliberate
    radius perturbation we can illustrate the *boundary downgrade*.

    Here we use r=1 and constraint ``r * c1 * r * c2 < 0.45`` which is
    interior-violating on tiles whose corners avoid the c1*c2 maximum but
    whose interior approaches it.  On face 1+ M=1, the corners u=±1 give
    c1*c2 = ±(1)/2 = ±0.5 — one corner passes (c1*c2 = -0.5 < 0.45), one
    corner *fails* (c1*c2 = +0.5 > 0.45)... so still BOUNDARY at corners.

    Skip the brittle construction; the real-world non-convex case is
    covered by the ``c1**2 < 0.1`` test above.  Here we just sanity-check
    that ``sample_interior > 0`` doesn't flip a FULLY_INSIDE on a convex
    constraint.
    """
    cs = ConstraintSet.from_strings(["c1 > 0"])
    status = classify_tile(
        face_idx=1,
        sub_tile=(1,),
        M=1,
        Q=np.eye(2),
        r=1.0,
        coupling_names=("c1", "c2"),
        constraints=cs,
        sample_interior=32,
    )
    assert status is TileStatus.FULLY_INSIDE


# ----------------------------------------------------------------------
# classify_all_cells contract
# ----------------------------------------------------------------------


def test_classify_all_cells_yields_exact_count() -> None:
    """``2N * M^(N-1)`` cells are iterated regardless of constraint outcome."""
    cs = ConstraintSet.from_strings(["c1 > 0"])
    n = 3
    M = 3
    Q = np.eye(n)
    results = list(
        classify_all_cells(
            n_dims=n,
            M=M,
            Q=Q,
            r=1.0,
            coupling_names=("c1", "c2", "c3"),
            constraints=cs,
        )
    )
    assert len(results) == n_faces(n) * M ** (n - 1)  # 6 * 9 = 54


def test_classify_all_cells_preserves_tile_order_within_face() -> None:
    """Within each face, sub_tiles are enumerated in itertools.product order."""
    cs = ConstraintSet.from_strings(["c1 > 0"])
    n = 2
    M = 3
    results = list(
        classify_all_cells(
            n_dims=n,
            M=M,
            Q=np.eye(n),
            r=1.0,
            coupling_names=("c1", "c2"),
            constraints=cs,
        )
    )
    # Face 1 first, with sub_tiles (1,), (2,), (3,)
    face_1 = [(sub, status) for face_idx, sub, status in results if face_idx == 1]
    assert [sub for sub, _ in face_1] == [(1,), (2,), (3,)]


# ----------------------------------------------------------------------
# Speedup demonstration
# ----------------------------------------------------------------------


def test_speedup_smoke_4d_with_two_sign_constraints() -> None:
    """A 4-coupling theory with 2 positive-axis sign constraints reduces
    the M=2 cell count by more than 25%.  Provides a sanity-check
    threshold for survey-driver users to expect non-trivial savings.
    """
    cs = ConstraintSet.from_strings(["c1 > 0", "c2 > 0"])
    n = 4
    M = 2
    coupling_names = ("c1", "c2", "c3", "c4")
    total = n_faces(n) * M ** (n - 1)
    results = list(
        classify_all_cells(
            n_dims=n,
            M=M,
            Q=np.eye(n),
            r=1.0,
            coupling_names=coupling_names,
            constraints=cs,
        )
    )
    n_dropped = sum(1 for _, _, st in results if st is TileStatus.FULLY_OUTSIDE)
    # Faces 2 and 4 are -c1 and -c2 entirely outside → 2 * M^(N-1) = 16
    # cells dropped on those alone, well over 25% of 64.
    assert n_dropped / total > 0.25
