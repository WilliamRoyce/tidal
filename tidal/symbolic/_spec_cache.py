"""Process-level cache for parsed equation systems.

Several hot paths (the inference likelihood pre-flight, the sweep
pre-flight, the prior-stability sweep) parsed the equation JSON on
every call.  At ~1–3 ms per parse and ~10⁵ likelihood evaluations per
chain, that's a few minutes of pure JSON-decoding overhead per worker
per run.  This cache amortizes the cost: the second time a path is
loaded within the process, the parsed :class:`EquationSystem` is
returned directly.

The cache is keyed on ``(absolute_path, mtime)`` so editing the JSON
on disk invalidates the entry automatically — important when a user
re-derives a spec mid-session.

The cache is intentionally module-level (not per-call) and not
weakref'd: the dominant use case is one or two specs per worker held
for the lifetime of a chain, and ``EquationSystem`` is small enough
that retaining it across the whole run costs nothing meaningful.

Issue #327 — performance plan ("Phase 1 — likelihood spec reuse").
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    from tidal.symbolic.json_loader import EquationSystem


# Global cache.  Module-level: shared across all callers in the
# process.  Cleared automatically when the process exits.
_SPEC_CACHE: dict[tuple[str, float], EquationSystem] = {}


def load_spec_cached(path: Path) -> EquationSystem:
    """Return the parsed :class:`EquationSystem` for *path*, cached.

    Cache key: ``(absolute_path, mtime)`` — editing the JSON on disk
    invalidates the entry on the next call.

    Parameters
    ----------
    path : Path
        Path to the equation JSON.  Resolved to an absolute path
        before it becomes the cache key, so callers may pass either
        absolute or relative paths.
    """
    from tidal.symbolic.json_loader import load_equation_system

    abs_path = path.resolve()
    mtime = abs_path.stat().st_mtime
    key = (str(abs_path), mtime)
    cached = _SPEC_CACHE.get(key)
    if cached is None:
        cached = load_equation_system(abs_path)
        _SPEC_CACHE[key] = cached
    return cached


def clear_cache() -> None:
    """Drop all cached entries.

    Used by tests that mutate spec JSON files; production code never
    needs this because mtime-based invalidation handles edits.
    """
    _SPEC_CACHE.clear()
