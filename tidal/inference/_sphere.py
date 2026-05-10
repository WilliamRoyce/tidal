r"""Cubed-sphere geometry for the radial-angular coupling-space prior.

Provides the angular chart used by :class:`RadialAngularPrior` in
:mod:`tidal.inference._prior`.  The chart wraps a cube ``[-1, 1]^N``
around the origin and projects each cube face radially out to the unit
sphere ``S^(N-1)``; the per-face cube coordinates ``u in [-1, 1]^(N-1)``
become a quasi-rectilinear chart on a sphere region with no polar
crowding.

The geometry follows the standard cubed-sphere construction (gnomonic
projection — Ronchi, Iacono & Paolucci 1996).  Face indexing, sub-tile
labelling, and orientation conventions match those of the supervisor's
reference implementation so per-tile outputs and labels are
interchangeable for cross-checking.

Public API
----------
``n_faces(n_dims)``                    — number of faces (= 2 * n_dims).
``face_to_axis_sign(face_idx)``        — 1-indexed face -> (axis k, sign).
``face_label(face_idx)``               — ``"01p"`` / ``"03m"`` ASCII tag.
``face_label_math(face_idx)``          — LaTeX ``r"1\\uparrow"``.
``face_to_direction(face_idx, u, Q)``  — gnomonic + Q^T to physical axes.
``tile_bounds(sub_tile, M)``           — sub-tile cube bounds.
``enumerate_tiles(n_dims, M)``         — tile-index tuples per face.
``enumerate_cells(n_dims, M)``         — (face_idx, tile_indices) pairs.
``random_rotation(n_dims, seed)``      — proper rotation via QR (det +1).
``tile_label(sub_tile)``               — ``"2_3"`` ASCII tag.
``cell_label(face_idx, sub_tile)``     — ``"01p_tile2_3"`` directory key.
"""

from __future__ import annotations

import itertools
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Iterator

    from numpy.typing import NDArray


def n_faces(n_dims: int) -> int:
    """Return the number of faces of an n_dims-dimensional cube (= 2 * n_dims)."""
    if n_dims < 2:
        msg = f"n_dims must be >= 2 (got {n_dims})"
        raise ValueError(msg)
    return 2 * n_dims


def face_to_axis_sign(face_idx: int) -> tuple[int, float]:
    """Decode a 1-indexed face id into ``(axis k, sign s)``.

    Convention: face 1 = +x_1, face 2 = -x_1, face 3 = +x_2, ...
    """
    if face_idx < 1:
        msg = f"face_idx must be >= 1 (got {face_idx})"
        raise ValueError(msg)
    k = (face_idx - 1) // 2 + 1  # 1-indexed axis
    s = +1.0 if (face_idx - 1) % 2 == 0 else -1.0
    return k, s


def face_label(face_idx: int) -> str:
    """ASCII face label, e.g. ``"01p"`` for face 1+, ``"03m"`` for face 3-."""
    k, s = face_to_axis_sign(face_idx)
    return f"{k:02d}{'p' if s > 0 else 'm'}"


def face_label_math(face_idx: int) -> str:
    r"""LaTeX face label, e.g. ``r"1\\uparrow"`` / ``r"6\\downarrow"``."""
    k, s = face_to_axis_sign(face_idx)
    arrow = r"\uparrow" if s > 0 else r"\downarrow"
    return f"{k}{arrow}"


def face_to_direction(
    face_idx: int,
    u: NDArray[np.float64],
    Q: NDArray[np.float64],  # noqa: N803 (Q matches mathematical convention)
) -> NDArray[np.float64]:
    r"""Map face index + face coordinates ``u in [-1, 1]^(N-1)`` to a unit N-vector on ``S^(N-1)``.

    The gnomonic projection is

    .. math::

        x = \\frac{v}{\\|v\\|}, \\qquad v_k = s, \\; v_j = u_{j'}

    where ``k = (face_idx - 1) // 2 + 1`` (1-indexed dominant axis),
    ``s = +1`` for faces 1, 3, ... and ``-1`` for 2, 4, ..., and ``j'``
    is the index of ``u`` corresponding to the j-th non-dominant axis.
    The optional rotation ``Q`` is applied by ``Q^T x`` to map to the
    physical (un-rotated) axes; pass ``np.eye(n_dims)`` for the
    identity (the default in higher-level callers).

    Parameters
    ----------
    face_idx : int
        1-indexed face id in ``[1, 2 * n_dims]``.
    u : array, shape (N - 1,)
        Cube-face coordinates in ``[-1, 1]``.
    Q : array, shape (N, N)
        Pre-rotation matrix.  ``Q^T x`` gives the physical-coords vector.

    Returns
    -------
    x : array, shape (N,)
        Unit-norm direction on ``S^(N-1)`` in physical coordinates.
    """
    q_arr = np.asarray(Q)
    n = q_arr.shape[0]
    if q_arr.shape != (n, n):
        msg = f"Q must be square; got shape {q_arr.shape}"
        raise ValueError(msg)
    u = np.asarray(u, dtype=np.float64)
    if u.shape != (n - 1,):
        msg = f"u must have shape ({n - 1},); got {u.shape}"
        raise ValueError(msg)
    k, s = face_to_axis_sign(face_idx)
    if k > n:
        msg = f"face_idx {face_idx} exceeds 2 * n_dims = {2 * n}"
        raise ValueError(msg)

    x = np.zeros(n, dtype=np.float64)
    x[k - 1] = s
    j = 0
    for i in range(n):
        if i == k - 1:
            continue
        x[i] = u[j]
        j += 1
    norm = float(np.linalg.norm(x))
    x /= norm
    return q_arr.T @ x


def tile_bounds(
    sub_tile: tuple[int, ...],
    M: int,  # noqa: N803 (M matches mathematical convention for subdivision count)
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    r"""Return ``(u_lo, u_hi)`` for one M-subdivided tile of a face.

    ``sub_tile`` is a length-(N-1) tuple of 1-indexed integers in
    ``[1, M]``; ``M`` is the per-axis subdivision count.  The bounds are

    .. math::

        u_{lo,i} = -1 + 2 (s_i - 1) / M, \\qquad
        u_{hi,i} = -1 + 2 s_i / M.

    ``M = 1`` returns ``(-1, +1)`` per axis (the whole face).
    """
    if M < 1:
        msg = f"M must be >= 1 (got {M})"
        raise ValueError(msg)
    t = np.asarray(sub_tile, dtype=int)
    if t.ndim != 1:
        msg = f"sub_tile must be 1-D; got shape {t.shape}"
        raise ValueError(msg)
    if (t < 1).any() or (t > M).any():
        msg = f"sub_tile entries must be in [1, {M}]; got {sub_tile}"
        raise ValueError(msg)
    u_lo = -1.0 + 2.0 * (t - 1) / M
    u_hi = -1.0 + 2.0 * t / M
    return u_lo.astype(np.float64), u_hi.astype(np.float64)


def tile_label(sub_tile: tuple[int, ...]) -> str:
    """Underscore-joined sub-tile label, e.g. ``(2, 3) -> "2_3"``."""
    return "_".join(str(i) for i in sub_tile)


def cell_label(face_idx: int, sub_tile: tuple[int, ...]) -> str:
    """Per-cell directory key, e.g. ``"01p_tile2_3"`` for face 1+, sub-tile (2, 3)."""
    return f"{face_label(face_idx)}_tile{tile_label(sub_tile)}"


def enumerate_tiles(
    n_dims: int,
    M: int,  # noqa: N803
) -> Iterator[tuple[int, ...]]:
    """Yield all length-(N-1) tile-index tuples for an M-subdivided face."""
    if n_dims < 2:
        msg = f"n_dims must be >= 2 (got {n_dims})"
        raise ValueError(msg)
    if M < 1:
        msg = f"M must be >= 1 (got {M})"
        raise ValueError(msg)
    return itertools.product(*([range(1, M + 1)] * (n_dims - 1)))


def enumerate_cells(
    n_dims: int,
    M: int,  # noqa: N803
) -> Iterator[tuple[int, tuple[int, ...]]]:
    """Yield ``(face_idx, sub_tile)`` for all ``2N * M^(N-1)`` cells."""
    for face_idx in range(1, n_faces(n_dims) + 1):
        for sub_tile in enumerate_tiles(n_dims, M):
            yield face_idx, sub_tile


def random_rotation(
    n_dims: int,
    seed: int = 0,
) -> NDArray[np.float64]:
    """Generate a proper rotation matrix ``Q`` (``n_dims x n_dims``, ``det = +1``) via QR decomposition of an N(0, I) random matrix.

    Pre-rotating coupling space by ``Q`` breaks alignment between cube
    faces and theory symmetry axes (a "trick for more readable plots"
    per supervisor 10 May 2026; not load-bearing physics).  The default
    in higher-level callers is the identity, which keeps face boundaries
    aligned with the Lagrangian's natural axes.
    """
    if n_dims < 2:
        msg = f"n_dims must be >= 2 (got {n_dims})"
        raise ValueError(msg)
    rng = np.random.default_rng(seed)
    A = rng.standard_normal((n_dims, n_dims))  # noqa: N806
    Q, _ = np.linalg.qr(A)  # noqa: N806
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q.astype(np.float64)
