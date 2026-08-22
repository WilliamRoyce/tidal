"""Integration test for GH #384 Phase A′ convolution-block cache.

Verifies that two ``solve_modal`` calls at different BSM couplings produce
results identical to a no-cache baseline, and that the cache reports
high hit rate on the broken Phase E theory.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from tidal.solver._conv_block_cache import clear as clear_cache
from tidal.solver._conv_block_cache import stats as cache_stats
from tidal.solver.grid import GridInfo
from tidal.solver.modal import solve_modal
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

BROKEN_THEORY_PATH = Path(
    "examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"
)

# Skip whole module if the broken-theory JSON isn't present (e.g. in some
# reduced test environments).
pytestmark = [
    # GH #455: this module drives the nonminimal dual-Gaussian spec through
    # solve_modal; the corrected #444 velocity-row scale exposes the
    # ungauged-gravity near-singular composition (operator abscissa
    # 0.51 -> 1594), making these runs diverge or hang. Cache-mechanics
    # coverage remains in test_conv_block_cache.py. Re-enable once #455
    # adjudicates (WS2 constraint-residual oracle).
    pytest.mark.skip(
        reason="GH #455: near-singular ungauged-gravity "
        "composition under the corrected #444 scale — see issue"
    ),
    pytest.mark.skipif(
        not BROKEN_THEORY_PATH.exists(),
        reason=f"Broken-theory JSON missing: {BROKEN_THEORY_PATH}",
    ),
]


@pytest.fixture
def broken_theory_setup() -> dict:
    """Load the broken theory at canonical Phase E geometry."""
    with Path(BROKEN_THEORY_PATH).open(encoding="utf-8") as f:
        spec_raw = EquationSystem.from_dict(json.load(f))
    spec_with_bsm = dataclasses.replace(
        spec_raw,
        metadata={
            **spec_raw.metadata,
            "_inference_sampled_params": (
                "alpha1",
                "alpha2",
                "alpha3",
                "delta1",
            ),
        },
    )
    base_params = {
        "kappa": 1.0,
        "Bpeak": 0.01,
        "sigB": 5.0,
        "zc1": 50.0,
        "zc2": 150.0,
    }
    grid = GridInfo(
        shape=(64,),
        bounds=((0.0, 200.0),),
        periodic=(True,),
    )
    layout = StateLayout.from_spec(spec_raw, grid.num_points)
    y0 = np.zeros(layout.total_size)
    # IC on h_5 (Phase E sampled field)
    if "h_5" in layout.field_slot_map:
        h5_slot = layout.field_slot_map["h_5"]
        h5_start = h5_slot * grid.num_points
        x = np.linspace(0.0, 200.0, grid.num_points, endpoint=False)
        y0[h5_start : h5_start + grid.num_points] = 1e-4 * np.exp(
            -((x - 100.0) ** 2) / 100.0,
        )
    return {
        "spec_raw": spec_raw,
        "spec_with_bsm": spec_with_bsm,
        "params_base": base_params,
        "grid": grid,
        "y0": y0,
    }


def _run(spec, params, grid, y0):
    return solve_modal(
        spec=spec,
        grid=grid,
        y0=y0,
        t_span=(0.0, 0.2),
        parameters=params,
        num_snapshots=2,
        return_eigendata=False,
    )


def test_cached_matches_uncached_bit_for_bit(broken_theory_setup: dict) -> None:
    """Two solve_modal calls — one with cache engaged via BSM metadata,
    one without — must produce identical state vectors.
    """
    clear_cache()
    params = {
        **broken_theory_setup["params_base"],
        "alpha1": 0.1,
        "alpha2": -0.1,
        "alpha3": 0.2,
        "delta1": -0.05,
    }
    res_no_cache = _run(
        broken_theory_setup["spec_raw"],
        params,
        broken_theory_setup["grid"],
        broken_theory_setup["y0"],
    )
    clear_cache()  # isolate the second run from any prior cache state
    res_cached = _run(
        broken_theory_setup["spec_with_bsm"],
        params,
        broken_theory_setup["grid"],
        broken_theory_setup["y0"],
    )
    y_no = res_no_cache["y"][-1]
    y_yes = res_cached["y"][-1]
    assert y_no.shape == y_yes.shape
    diff = np.abs(y_no - y_yes).max()
    scale = max(np.abs(y_no).max(), 1e-30)
    assert diff / scale < 1e-10, (
        f"Cached run diverged from no-cache: max rel diff {diff / scale:.2e}"
    )


def test_cross_call_cache_hits_on_BSM_only_variation(
    broken_theory_setup: dict,
) -> None:
    """Two cached calls at the SAME geometry but different BSM samples
    should reuse the cached blocks — hit rate climbs above 50% on the
    second call.
    """
    clear_cache()
    p1 = {
        **broken_theory_setup["params_base"],
        "alpha1": 0.1,
        "alpha2": -0.1,
        "alpha3": 0.2,
        "delta1": -0.05,
    }
    p2 = {
        **broken_theory_setup["params_base"],
        "alpha1": 0.9,
        "alpha2": +0.4,
        "alpha3": -0.3,
        "delta1": +0.7,
    }
    _run(
        broken_theory_setup["spec_with_bsm"],
        p1,
        broken_theory_setup["grid"],
        broken_theory_setup["y0"],
    )
    s1 = cache_stats()
    _run(
        broken_theory_setup["spec_with_bsm"],
        p2,
        broken_theory_setup["grid"],
        broken_theory_setup["y0"],
    )
    s2 = cache_stats()
    new_hits = s2["hits"] - s1["hits"]
    new_misses = s2["misses"] - s1["misses"]
    hit_rate = new_hits / max(new_hits + new_misses, 1)
    assert hit_rate > 0.5, (
        f"Cross-call cache hit rate {hit_rate:.2%} too low; "
        f"new hits={new_hits}, new misses={new_misses}"
    )
