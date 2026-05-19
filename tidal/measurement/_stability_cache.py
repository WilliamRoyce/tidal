"""Process-level cache of parameter-independent probe structural state.

The pre-flight stability probe in :mod:`tidal.measurement._stability`
recomputes a handful of structural objects on every call:

* :class:`tidal.solver.state.StateLayout` from ``spec`` and grid size.
* The Fourier wavenumber grid ``k_vals``, ``rfft_shape``.
* The source-field slot lookup against the layout.

None of these depend on the parameter dict.  In an MCMC chain with
~10⁵ probe calls per worker, recomputing them on every call is pure
overhead.  This module memoises them per ``(id(spec), grid-tuple,
source)`` so subsequent calls are O(dict-lookup).

Performance budget contributions (measured at the
``dark_photon_plasma`` fixture, N=32, ~27 ms total per probe call):

* ``StateLayout.from_spec`` ≈ 0.1–0.2 ms.
* k-grid construction ≈ <0.1 ms.
* Slot lookup decision tree ≈ negligible.

Total saving: ~0.3–0.5 ms per call.  Modest but free-of-risk.

Issue #327 — performance plan ("Phase 2 — structural cache").
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver.grid import GridInfo
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem


@dataclass(frozen=True)
class ProbeStructuralContext:
    """Parameter-independent context for ``check_conversion_stability``.

    Computed once per ``(spec, grid, source)`` and cached.  Members
    are intentionally narrow: anything that depends on the parameter
    dict (e.g. the per-mode evolution matrices) is *not* here.
    """

    layout: StateLayout
    """The layout for ``(spec, grid.num_points)``."""

    k_vals: NDArray[np.float64]
    """Real-fft wavenumber grid for the probe's 1-D Fourier scan."""

    rfft_shape: tuple[int, ...]
    """rfft output shape ``(N // 2 + 1,)`` for the 1-D probe."""

    ic_k_default: float
    """Default IC wavenumber ``2π / L`` if the caller passes no
    ``ic_wavevector``."""

    src_slot_full: int | None
    """Full-layout slot index for the source field; ``None`` if the
    source field is not registered (early-return path)."""

    early_return_reason: str | None
    """If non-``None``, the probe should immediately return a stable
    result with this diagnostic message — covers the "field not
    found" case.  Other early-return cases (constraint field, no
    coupled block) require ``mapping`` from
    ``_build_evolution_matrices`` and are handled by the caller."""


# Module-level cache.  Key: ``(id(spec), grid-tuple, source)``.  The
# ``id(spec)`` half is safe because :class:`EquationSystem` is a
# frozen dataclass and the spec cache in
# ``tidal.symbolic._spec_cache`` keeps a stable reference for the
# lifetime of the process.  Cleared automatically at process exit.
_STRUCTURAL_CACHE: dict[
    tuple[int, tuple[int, ...], tuple[tuple[float, float], ...], tuple[bool, ...], str],
    ProbeStructuralContext,
] = {}


def get_structural_context(
    spec: EquationSystem,
    grid: GridInfo,
    source: str,
) -> ProbeStructuralContext:
    """Return the cached :class:`ProbeStructuralContext` for the inputs.

    Builds and caches on first call; returns the cached entry on
    subsequent calls.  The cache key is
    ``(id(spec), grid.shape, grid.bounds, grid.periodic, source)`` —
    safe because both ``spec`` and ``grid`` are frozen dataclasses
    and the spec cache keeps spec identities stable across MCMC
    evaluations.
    """
    from tidal.solver.state import StateLayout

    key = (
        id(spec),
        tuple(grid.shape),
        tuple(grid.bounds),
        tuple(grid.periodic),
        source,
    )
    cached = _STRUCTURAL_CACHE.get(key)
    if cached is not None:
        return cached

    layout = StateLayout.from_spec(spec, grid.num_points)

    n = int(grid.shape[0])
    dx = float(grid.dx[0])
    k_vals = np.asarray(2 * math.pi * np.fft.rfftfreq(n, dx), dtype=np.float64)
    rfft_shape: tuple[int, ...] = (n // 2 + 1,)

    L = float(grid.bounds[0][1] - grid.bounds[0][0])
    ic_k_default = 2 * math.pi / L

    src_slot_full: int | None
    early_return_reason: str | None
    try:
        src_slot_full = int(layout.field_slot_map[source])
        early_return_reason = None
    except KeyError:
        src_slot_full = None
        early_return_reason = f"Source field {source!r} not found in reduced system."

    ctx = ProbeStructuralContext(
        layout=layout,
        k_vals=k_vals,
        rfft_shape=rfft_shape,
        ic_k_default=ic_k_default,
        src_slot_full=src_slot_full,
        early_return_reason=early_return_reason,
    )
    _STRUCTURAL_CACHE[key] = ctx
    return ctx


def clear_cache() -> None:
    """Drop all cached entries.

    Used by tests that mutate spec or grid identities.  Production
    code never needs this.
    """
    _STRUCTURAL_CACHE.clear()
