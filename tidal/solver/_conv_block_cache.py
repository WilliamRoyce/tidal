"""Module-level cache for convolution-matrix blocks across PolyChord likelihood calls.

GH #384 Phase A′: for position-dependent BSM-separable RHS terms, the
``(n_modes × n_modes)`` complex convolution block produced by
``_add_convolution_coupling`` / ``_term_conv_block`` in ``tidal.solver.modal``
is identical across PolyChord calls modulo a multiplicative BSM scalar.
Cache it at chain start, scale per call.

Cache architecture
------------------
- **Scope**: module-level dict, per-Python-process. PolyChord's MPI ranks
  are independent Python processes, so the cache is naturally per-rank
  (no locking needed).
- **Lifetime**: until process exit or explicit ``clear()``. Keyed on a
  string-based hash that's stable across function calls in the same process.
- **Memory**: each entry is an ``(n_modes, n_modes)`` complex128 array.
  For a 38-field PGT theory at N=128: ~30 MB per rank for ~400 cached terms.

Key contents
------------
The cache key carries everything that affects the block contents:

- ``coeff_geom_str``: the coefficient expression with BSM symbols already
  substituted out (or set to 1). Stable across PolyChord calls because
  PolyChord only varies BSM symbols.
- ``operator_name``: e.g. ``"laplacian_x"``. Determines the Fourier
  multiplier applied per-mode.
- ``grid_shape``: tuple of ints. Determines ``n_modes`` and the FFT size.
- ``geometry_hash``: hash of the geometric parameters (everything in
  ``parameters`` that's NOT in ``bsm_symbols``). Same chain → same hash.

Public API
----------
- ``get_or_compute(key, compute_fn)``: standard memoise.
- ``clear()``: drop all cached blocks (test isolation).
- ``stats()``: return ``{"hits": ..., "misses": ..., "bytes": ...}``.

Linked: GH #379 (the constraint-Schur path being optimised),
GH #384 (this optimisation), the [[gh384-conv-block-cache]] memory.
"""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Callable

    import numpy as np
    from numpy.typing import NDArray


CacheKey = tuple[str, str, tuple[int, ...], int]
"""(coeff_geom_str, operator_name, grid_shape, geometry_hash)"""


_CACHE: dict[CacheKey, NDArray[np.complex128]] = {}
_STATS = {"hits": 0, "misses": 0}
_LOCK = threading.Lock()  # defensive; single-rank Python is single-threaded anyway


def make_geometry_hash(
    parameters: dict[str, float],
    bsm_symbols: frozenset[str] | set[str] | tuple[str, ...],
) -> int:
    """Stable hash of the non-BSM (i.e. geometry) parameters.

    Two PolyChord calls with the same geometry and different BSM samplings
    produce the same hash. A rerun with different geometry produces a
    different hash and the cache invalidates naturally.
    """
    bsm_set = frozenset(bsm_symbols)
    geom_items = sorted((k, v) for k, v in parameters.items() if k not in bsm_set)
    return hash(tuple(geom_items))


def make_key(
    coeff_geom_str: str,
    operator_name: str,
    grid_shape: tuple[int, ...],
    geometry_hash: int,
) -> CacheKey:
    """Build a cache key. See module docstring for components."""
    return (coeff_geom_str, operator_name, tuple(grid_shape), int(geometry_hash))


def get_or_compute(
    key: CacheKey,
    compute_fn: Callable[[], NDArray[np.complex128]],
) -> NDArray[np.complex128]:
    """Return the cached block, computing and storing it on miss."""
    with _LOCK:
        cached = _CACHE.get(key)
        if cached is not None:
            _STATS["hits"] += 1
            return cached
    # Compute outside the lock — compute_fn may be expensive.
    block = compute_fn()
    with _LOCK:
        # Re-check in case another path raced us
        cached = _CACHE.get(key)
        if cached is not None:
            _STATS["hits"] += 1
            return cached
        _CACHE[key] = block
        _STATS["misses"] += 1
    return block


def clear() -> None:
    """Drop all cached blocks and reset stats.

    Call between independent chains in tests; not needed for normal use.
    """
    with _LOCK:
        _CACHE.clear()
        _STATS["hits"] = 0
        _STATS["misses"] = 0


def stats() -> dict[str, int]:
    """Return current cache statistics.

    Keys:
        hits: number of get_or_compute calls that hit the cache.
        misses: number that computed a new block.
        size: number of cached blocks.
        bytes: total memory footprint of cached blocks.
    """
    with _LOCK:
        total_bytes = sum(blk.nbytes for blk in _CACHE.values())
        return {
            "hits": _STATS["hits"],
            "misses": _STATS["misses"],
            "size": len(_CACHE),
            "bytes": total_bytes,
        }
