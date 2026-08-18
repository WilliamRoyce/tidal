"""Determinism of anesthetic-backed rendering and bootstrap statistics.

anesthetic draws from NumPy's *legacy* global RNG in two places that
matter here: ``compress_weights`` / ``triangular_sample_compression_2d``
during KDE plotting, and ``logX`` inside ``NestedSamples.stats``.  Neither
accepts a ``Generator``, so the only lever is the global state.

``_deterministic_render`` seeds it and — the part that matters — restores
it afterwards, so a render cannot leave the process in a fixed state for
whatever runs next.  See issue #388.
"""

from __future__ import annotations

import numpy as np

from tidal.inference._visualize import _deterministic_render


class TestDeterministicRender:
    def test_restores_global_state(self) -> None:
        """The context manager must not leak a seeded state to its caller."""
        np.random.seed(12345)  # noqa: NPY002
        before = np.random.get_state()  # noqa: NPY002

        with _deterministic_render():
            np.random.rand(100)  # noqa: NPY002  # consume, as anesthetic does

        after = np.random.get_state()  # noqa: NPY002
        assert before[0] == after[0]
        assert np.array_equal(before[1], after[1])
        assert before[2] == after[2]

    def test_inner_draws_reproduce_across_entries(self) -> None:
        with _deterministic_render():
            first = np.random.rand(4)  # noqa: NPY002

        np.random.rand(999)  # noqa: NPY002  # perturb the outer stream

        with _deterministic_render():
            second = np.random.rand(4)  # noqa: NPY002

        assert np.array_equal(first, second)

    def test_outer_stream_unaffected_by_inner_consumption(self) -> None:
        """How much the block consumes must not shift the caller's stream."""
        np.random.seed(7)  # noqa: NPY002
        expected = np.random.rand()  # noqa: NPY002

        np.random.seed(7)  # noqa: NPY002
        with _deterministic_render():
            np.random.rand(50)  # noqa: NPY002
        actual = np.random.rand()  # noqa: NPY002

        assert expected == actual

    def test_distinct_seeds_give_distinct_draws(self) -> None:
        """Negative control: a degenerate implementation would pass the rest."""
        with _deterministic_render(seed=0):
            a = np.random.rand(4)  # noqa: NPY002
        with _deterministic_render(seed=1):
            b = np.random.rand(4)  # noqa: NPY002

        assert not np.array_equal(a, b)
