"""Shared utilities for mode-matching between dispersion curves."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray


def find_shared_modes(
    k_source: NDArray[np.float64],
    k_target: NDArray[np.float64],
    active_source: NDArray[np.bool_] | None = None,
    active_target: NDArray[np.bool_] | None = None,
    *,
    rtol: float = 1e-6,
) -> tuple[list[int], list[int]]:
    """Find matching wavenumber indices between two arrays.

    Returns lists of (source_indices, target_indices) where
    ``k_source[i] ≈ k_target[j]`` within relative tolerance.

    Parameters
    ----------
    k_source, k_target : ndarray
        Wavenumber arrays.
    active_source, active_target : ndarray of bool or None
        Optional masks for active modes. If None, all modes are active.
    rtol : float
        Relative tolerance for matching.

    Returns
    -------
    tuple[list[int], list[int]]
        Paired index lists ``(idx_source, idx_target)``.
    """
    idx_src: list[int] = []
    idx_tgt: list[int] = []
    for i, k_s in enumerate(k_source):
        if active_source is not None and not active_source[i]:
            continue
        mask = np.isclose(k_target, k_s, rtol=rtol)
        if active_target is not None:
            mask &= active_target
        matches = np.where(mask)[0]
        if len(matches) > 0:
            idx_src.append(i)
            idx_tgt.append(matches[0])
    return idx_src, idx_tgt
