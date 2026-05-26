"""Tests for tidal.solver._conv_block_cache.

GH #384 Phase A′ cache infrastructure: get/compute, hit/miss accounting,
key stability under BSM-only variation, key change under geometry change.
"""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver._conv_block_cache import (
    clear,
    get_or_compute,
    make_geometry_hash,
    make_key,
    stats,
)


@pytest.fixture(autouse=True)
def _isolate_cache() -> None:
    """Reset the module-level cache between every test."""
    clear()


def _block(n: int = 4) -> np.ndarray:
    """Sample dense complex block for testing."""
    rng = np.random.default_rng(seed=12345)
    return rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))


def test_first_call_is_miss_subsequent_are_hits() -> None:
    key = make_key("x", "laplacian_x", (8,), 42)
    n_computes = [0]

    def compute() -> np.ndarray:
        n_computes[0] += 1
        return _block()

    out1 = get_or_compute(key, compute)
    out2 = get_or_compute(key, compute)
    out3 = get_or_compute(key, compute)

    assert n_computes[0] == 1, "compute_fn must run exactly once"
    assert out1 is out2 is out3, "cache must return the same object"
    s = stats()
    assert s["hits"] == 2
    assert s["misses"] == 1
    assert s["size"] == 1


def test_distinct_keys_distinct_blocks() -> None:
    block_a = _block()
    block_b = _block(n=8)  # different shape on purpose

    out_a = get_or_compute(make_key("x", "laplacian_x", (8,), 1), lambda: block_a)
    out_b = get_or_compute(make_key("y", "laplacian_x", (8,), 1), lambda: block_b)
    assert out_a is block_a
    assert out_b is block_b
    s = stats()
    assert s["size"] == 2
    assert s["misses"] == 2


def test_geometry_hash_stable_under_bsm_change() -> None:
    """Two parameter dicts that differ only in BSM symbols must hash identically."""
    bsm = {"alpha1", "delta1"}
    p1 = {"alpha1": 0.1, "delta1": -0.3, "Bpeak": 0.01, "sigB": 5.0}
    p2 = {"alpha1": 0.9, "delta1": +0.7, "Bpeak": 0.01, "sigB": 5.0}
    assert make_geometry_hash(p1, bsm) == make_geometry_hash(p2, bsm)


def test_geometry_hash_changes_when_geometry_changes() -> None:
    bsm = {"alpha1"}
    p_low = {"alpha1": 0.1, "Bpeak": 0.01, "sigB": 5.0}
    p_high = {"alpha1": 0.1, "Bpeak": 0.02, "sigB": 5.0}  # geometry varied
    assert make_geometry_hash(p_low, bsm) != make_geometry_hash(p_high, bsm)


def test_cross_call_hits_across_BSM_change() -> None:
    """The headline GH #384 use case: two solve_modal calls at different
    BSM couplings must share the cached block.
    """
    bsm = {"alpha1"}
    h1 = make_geometry_hash({"alpha1": 0.1, "Bpeak": 0.01}, bsm)
    h2 = make_geometry_hash({"alpha1": 0.9, "Bpeak": 0.01}, bsm)
    assert h1 == h2  # geometry unchanged

    n_computes = [0]

    def compute() -> np.ndarray:
        n_computes[0] += 1
        return _block()

    # Two calls with BSM-only variation
    key1 = make_key("Bpeak * x", "laplacian_x", (8,), h1)
    key2 = make_key("Bpeak * x", "laplacian_x", (8,), h2)
    out1 = get_or_compute(key1, compute)
    out2 = get_or_compute(key2, compute)
    assert out1 is out2
    assert n_computes[0] == 1


def test_clear_resets_state() -> None:
    key = make_key("x", "laplacian_x", (8,), 1)
    get_or_compute(key, _block)
    assert stats()["size"] == 1
    clear()
    s = stats()
    assert s["size"] == 0
    assert s["hits"] == 0
    assert s["misses"] == 0


def test_bytes_counter_tracks_memory() -> None:
    key = make_key("x", "laplacian_x", (8,), 1)
    block = _block(n=16)
    get_or_compute(key, lambda: block)
    s = stats()
    # complex128 → 16 bytes/entry × 16×16 = 4096 bytes
    assert s["bytes"] == block.nbytes
