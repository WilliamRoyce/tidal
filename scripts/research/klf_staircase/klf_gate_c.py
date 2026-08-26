"""Gate (c): the staircase prototype on the E.cal pencil (GH #473 stage 1).

usage: klf_gate_c.py N tau
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))  # research script: sibling imports
from klf_port import klf_deflate

from tidal.solver import modal
from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_convolution_matrix_with_constraints,
    _build_k_axes_full,
    _build_k_grid,
    _pencil_deflate,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

REPO = Path("/workspaces/torsion-gertsenshtein")
PARAMS = {"kappa": 1.0, "Bpeak": 0.01, "sigB": 5.0, "zc1": 25.0, "zc2": 75.0}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 16
TAU = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-6
OBS = ("h_5", "a_1")

spec = EquationSystem.from_dict(
    json.loads(
        (REPO / "examples/data/gertsenshtein_ungauged_e_dual_gaussian.json").read_text()
    )
)
grid = GridInfo(bounds=((0.0, 100.0),), shape=(N,), periodic=(True,))
layout = StateLayout.from_spec(spec, N)
ce = CoefficientEvaluator(spec, grid, PARAMS)
k_grid = _build_k_grid(_build_k_axes_full(grid))
cap: dict[str, np.ndarray] = {}


def _capture(A, B, *, context="", diagnostics=None, tag=("pencil", 0)):
    cap["A"], cap["B"] = np.array(A, copy=True), np.array(B, copy=True)
    return np.asarray(A, dtype=np.complex128), np.eye(A.shape[0], dtype=np.complex128)


modal._pencil_deflate = _capture
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    _, _, _, o2r, _ = _build_convolution_matrix_with_constraints(
        spec, layout, grid, ce, k_grid, tuple(grid.shape)
    )
modal._pencil_deflate = _pencil_deflate
A0, B0 = cap["A"], cap["B"]
n_tot = A0.shape[0]

obs_idx = np.concatenate(
    [
        np.arange(N) + o2r[m[f]] * N
        for f in OBS
        for m in (layout.field_slot_map, layout.velocity_slot_map)
    ]
)
mask = np.zeros(n_tot, dtype=bool)
mask[obs_idx] = True

t0 = time.perf_counter()
A_eff, proj, pinned, info = klf_deflate(A0, B0, TAU)
dt = time.perf_counter() - t0

blk_rel = np.linalg.norm(
    A_eff[np.ix_(mask, mask)] - A0[np.ix_(mask, mask)]
) / np.linalg.norm(A0[np.ix_(mask, mask)])
cross = max(
    np.linalg.norm(A_eff[np.ix_(mask, ~mask)]),
    np.linalg.norm(A_eff[np.ix_(~mask, mask)]),
)
ev = np.linalg.eigvals(A_eff)
print(
    f"N={N} dim={n_tot} tau={TAU:.0e}: contract={info['contract']:.2e}  "
    f"pinned={info['n_pinned']}  finite={info['n_finite']}  "
    f"maxRe={info['max_re']:+.3e}  |Im|max={float(np.max(np.abs(ev.imag))):.3f}  "
    f"obs-block rel={blk_rel:.2e} cross={cross:.2e}  [{dt:.1f}s]"
)
