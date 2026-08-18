"""Cubed-sphere geometry tests for :mod:`tidal.inference._sphere`."""

from __future__ import annotations

import numpy as np
import pytest

from tidal.inference._sphere import (
    cell_label,
    enumerate_cells,
    enumerate_tiles,
    face_label,
    face_label_math,
    face_to_axis_sign,
    face_to_direction,
    n_faces,
    random_rotation,
    tile_bounds,
    tile_label,
)

# ----------------------------------------------------------------------
# Face indexing
# ----------------------------------------------------------------------


@pytest.mark.parametrize(("n_dims", "expected"), [(2, 4), (3, 6), (4, 8), (6, 12)])
def test_n_faces(n_dims: int, expected: int) -> None:
    assert n_faces(n_dims) == expected


def test_n_faces_rejects_one_dim() -> None:
    with pytest.raises(ValueError, match="n_dims must be >= 2"):
        n_faces(1)


@pytest.mark.parametrize(
    ("face_idx", "expected"),
    [
        (1, (1, +1.0)),
        (2, (1, -1.0)),
        (3, (2, +1.0)),
        (4, (2, -1.0)),
        (11, (6, +1.0)),
        (12, (6, -1.0)),
    ],
)
def test_face_to_axis_sign(face_idx: int, expected: tuple[int, float]) -> None:
    assert face_to_axis_sign(face_idx) == expected


def test_face_to_axis_sign_rejects_zero() -> None:
    with pytest.raises(ValueError, match="face_idx must be >= 1"):
        face_to_axis_sign(0)


@pytest.mark.parametrize(
    ("face_idx", "ascii_lbl", "math_lbl"),
    [
        (1, "01p", r"1\uparrow"),
        (2, "01m", r"1\downarrow"),
        (5, "03p", r"3\uparrow"),
        (12, "06m", r"6\downarrow"),
    ],
)
def test_face_labels(face_idx: int, ascii_lbl: str, math_lbl: str) -> None:
    assert face_label(face_idx) == ascii_lbl
    assert face_label_math(face_idx) == math_lbl


# ----------------------------------------------------------------------
# Gnomonic projection
# ----------------------------------------------------------------------


def test_face_to_direction_unit_norm() -> None:
    n = 6
    Q = np.eye(n)
    rng = np.random.default_rng(42)
    for face_idx in range(1, n_faces(n) + 1):
        for _ in range(20):
            u = rng.uniform(-1.0, 1.0, size=n - 1)
            x = face_to_direction(face_idx, u, Q)
            assert x.shape == (n,)
            np.testing.assert_allclose(np.linalg.norm(x), 1.0, atol=1e-12)


def test_face_to_direction_face_center_is_axis() -> None:
    """``u = 0`` on face k+ projects to ``+e_k``; on face k- to ``-e_k``."""
    n = 4
    Q = np.eye(n)
    u_center = np.zeros(n - 1)
    for k in range(1, n + 1):
        face_plus = 2 * (k - 1) + 1
        face_minus = 2 * (k - 1) + 2
        x_plus = face_to_direction(face_plus, u_center, Q)
        x_minus = face_to_direction(face_minus, u_center, Q)
        expected_plus = np.zeros(n)
        expected_plus[k - 1] = 1.0
        np.testing.assert_allclose(x_plus, expected_plus, atol=1e-12)
        np.testing.assert_allclose(x_minus, -expected_plus, atol=1e-12)


def test_face_to_direction_rejects_wrong_shape() -> None:
    n = 4
    Q = np.eye(n)
    with pytest.raises(ValueError, match=r"u must have shape"):
        face_to_direction(1, np.zeros(n), Q)


def test_face_to_direction_rejects_face_out_of_range() -> None:
    n = 3
    Q = np.eye(n)
    with pytest.raises(ValueError, match=r"face_idx \d+ exceeds"):
        face_to_direction(7, np.zeros(n - 1), Q)


def test_face_to_direction_q_rotation_is_orthogonal_transform() -> None:
    """Pre-rotating Q only changes the physical-axis basis, not norms or angles."""
    n = 5
    Q = random_rotation(n, seed=7)
    Q_id = np.eye(n)
    rng = np.random.default_rng(11)
    for face_idx in range(1, n_faces(n) + 1):
        u = rng.uniform(-1.0, 1.0, size=n - 1)
        x_id = face_to_direction(face_idx, u, Q_id)
        x_q = face_to_direction(face_idx, u, Q)
        # x_q is Q^T @ x_id by construction
        np.testing.assert_allclose(x_q, Q.T @ x_id, atol=1e-12)
        # Both unit norm
        np.testing.assert_allclose(np.linalg.norm(x_q), 1.0, atol=1e-12)


# ----------------------------------------------------------------------
# Sub-tile bounds
# ----------------------------------------------------------------------


def test_tile_bounds_m_one_returns_full_face() -> None:
    u_lo, u_hi = tile_bounds((1,) * 5, M=1)
    np.testing.assert_array_equal(u_lo, -np.ones(5))
    np.testing.assert_array_equal(u_hi, +np.ones(5))


@pytest.mark.parametrize("M", [2, 3, 4, 8])
def test_tile_bounds_union_covers_face(M: int) -> None:
    """All M^(N-1) sub-tiles of a face partition ``[-1, 1]^(N-1)`` exactly."""
    n_dims = 4
    total = 0.0
    full_volume = 2.0 ** (n_dims - 1)
    for sub_tile in enumerate_tiles(n_dims, M):
        u_lo, u_hi = tile_bounds(sub_tile, M)
        widths = u_hi - u_lo
        total += float(np.prod(widths))
    np.testing.assert_allclose(total, full_volume, atol=1e-12)


def test_tile_bounds_rejects_out_of_range() -> None:
    with pytest.raises(ValueError, match=r"sub_tile entries must be in"):
        tile_bounds((0, 1), M=2)
    with pytest.raises(ValueError, match=r"sub_tile entries must be in"):
        tile_bounds((1, 3), M=2)


def test_tile_bounds_rejects_zero_m() -> None:
    with pytest.raises(ValueError, match="M must be >= 1"):
        tile_bounds((1, 1), M=0)


# ----------------------------------------------------------------------
# Enumeration
# ----------------------------------------------------------------------


@pytest.mark.parametrize(("n_dims", "M"), [(2, 1), (3, 2), (4, 3), (6, 2)])
def test_enumerate_tiles_count(n_dims: int, M: int) -> None:
    tiles = list(enumerate_tiles(n_dims, M))
    assert len(tiles) == M ** (n_dims - 1)
    # All tile indices should be 1-indexed in [1, M].
    for t in tiles:
        assert len(t) == n_dims - 1
        assert all(1 <= s <= M for s in t)


@pytest.mark.parametrize(("n_dims", "M"), [(2, 1), (3, 2), (4, 1), (5, 2)])
def test_enumerate_cells_count(n_dims: int, M: int) -> None:
    cells = list(enumerate_cells(n_dims, M))
    assert len(cells) == 2 * n_dims * M ** (n_dims - 1)
    face_indices = {f for f, _ in cells}
    assert face_indices == set(range(1, 2 * n_dims + 1))


# ----------------------------------------------------------------------
# Random rotation
# ----------------------------------------------------------------------


@pytest.mark.parametrize("n_dims", [2, 3, 5, 6])
def test_random_rotation_orthogonal(n_dims: int) -> None:
    Q = random_rotation(n_dims, seed=42)
    np.testing.assert_allclose(Q.T @ Q, np.eye(n_dims), atol=1e-12)
    np.testing.assert_allclose(Q @ Q.T, np.eye(n_dims), atol=1e-12)


@pytest.mark.parametrize("n_dims", [2, 3, 5, 6])
def test_random_rotation_proper(n_dims: int) -> None:
    """Det must be +1 (proper rotation, not improper)."""
    Q = random_rotation(n_dims, seed=42)
    np.testing.assert_allclose(np.linalg.det(Q), 1.0, atol=1e-12)


def test_random_rotation_seed_reproducible() -> None:
    Q1 = random_rotation(4, seed=99)
    Q2 = random_rotation(4, seed=99)
    np.testing.assert_array_equal(Q1, Q2)


def test_random_rotation_different_seeds_differ() -> None:
    Q1 = random_rotation(4, seed=1)
    Q2 = random_rotation(4, seed=2)
    assert not np.allclose(Q1, Q2)


# ----------------------------------------------------------------------
# Cell labeling
# ----------------------------------------------------------------------


def test_tile_label_format() -> None:
    assert tile_label((1,)) == "1"
    assert tile_label((2, 3)) == "2_3"
    assert tile_label((1, 1, 1, 1, 1)) == "1_1_1_1_1"


def test_cell_label_format() -> None:
    assert cell_label(1, (1,)) == "01p_tile1"
    assert cell_label(4, (2, 3)) == "02m_tile2_3"
    assert cell_label(11, (1, 1, 1, 1, 1)) == "06p_tile1_1_1_1_1"
