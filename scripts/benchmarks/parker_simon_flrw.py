"""Parker–Simon-style perturbative-reduction calibration.

Calibrates the perturbative-reduction Pass 0 + Pass 1 machinery on a
Parker--Simon-style higher-derivative correction (Parker & Simon 1993,
gr-qc/9211002). The base theory is a 1+1D massive Klein--Gordon wave;
the order-ε correction is a higher-derivative biharmonic operator,

    ∂_t² φ = ∂_x² φ - m² φ - ε · ∂_x⁴ φ

which is the canonical structure of a Parker--Simon-style reduction
of higher-derivative gravity to second-order field equations: an
order-zero base wave operator plus a small higher-derivative
correction. The full dispersion is ω² = k² + m² + ε·k⁴.

Pass 0 solves the base modal evolution (ε=0); Pass 1 evaluates the
closed-form Duhamel kernel for the higher-derivative correction
source on the Pass 0 trajectory. The combined Pass 0 + Pass 1
solution matches the full-ε modal solution to O(ε²) by construction.

This test is *distinct* from `forced_oscillator.py` (which probes
the order-zero mass-shift correction): the correction operator here
is fourth-order in space, exercising a different code path through
the perturbative-reduction source assembly.

The 1+0D scalar a^(1) recovery of Parker & Simon eqs. 3.39–3.44 is
documented in the planned `tests/test_perturbative_parker_simon.py`
(Task 6.2 in PERTURBATIVE_REDUCTION_IMPLEMENTATION.md); the present
benchmark provides the spatial-field-theoretic analogue.

Serves:   manuscript/sections/appendices/validation.tex (App D §4)
Consumes: scripts/figures/figD_parker_simon.py
Writes:   benchmark_results/canonical/parker_simon_flrw.json
"""

from __future__ import annotations

import argparse
import copy
import datetime
import json
import platform
import socket
import subprocess  # noqa: S404
import sys
from pathlib import Path
from typing import Any, cast

import numpy as np
import scipy

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.grid import GridInfo  # noqa: E402
from tidal.solver.modal import solve_modal, solve_modal_pass1  # noqa: E402
from tidal.solver.state import StateLayout  # noqa: E402
from tidal.symbolic.json_loader import EquationSystem  # noqa: E402

DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "parker_simon_flrw.json"

# 1+1D Klein-Gordon with biharmonic-x order-ε higher-derivative correction.
_PS_BASE: dict[str, Any] = {
    "metadata": {"source": "inline-benchmark", "parameters": {"m2": 1.0, "eps": 0.0}},
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
    "equations": [
        {
            "field": "phi_0",
            "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-m2",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                    {
                        "coefficient": -1.0,
                        "operator": "biharmonic",
                        "field": "phi_0",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        },
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2"]]},
}

# eps ladder restricted to the regime where err_combined is above the
# round-off floor (~1e-7 in this setup). Points below the floor are
# excluded to avoid misleading log-log slopes.
FULL_EPS_VALUES = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005]
SMOKE_EPS_VALUES = [0.1, 0.05, 0.01]
M2 = 1.0
T_END = 2.0
N_GRID = 64
K_MODE = 2.0  # higher mode amplifies biharmonic ε·k⁴ vs identity ε·1 contrast
N_SNAP = 21


def _git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _metadata(parameters: dict) -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "parameters": parameters,
    }


def _setup() -> tuple[EquationSystem, GridInfo, np.ndarray]:
    spec = EquationSystem.from_dict(copy.deepcopy(_PS_BASE))
    n = N_GRID
    length = 2 * np.pi
    grid = GridInfo(shape=(n,), bounds=((0.0, length),), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    x = np.linspace(0.0, length, n, endpoint=False)
    y0 = np.zeros(layout.num_slots * grid.num_points)
    y0[:n] = np.sin(K_MODE * x)
    return spec, grid, y0


def _run_eps(eps: float) -> dict[str, float]:
    spec, grid, y0 = _setup()
    base_spec = spec.filter_by_order(0)
    correction_spec = spec.filter_by_order(1)

    pass0 = cast(
        "dict[str, Any]",
        solve_modal(
            base_spec,
            grid,
            y0,
            t_span=(0.0, T_END),
            parameters={"m2": M2, "eps": eps},
            num_snapshots=N_SNAP,
            return_eigendata=True,
        ),
    )
    pass1 = solve_modal_pass1(
        pass0["eigendata"],
        correction_spec,
        grid,
        pass0["t"],
        parameters={"m2": M2, "eps": eps},
    )
    q_pass0 = pass0["y"]
    q_total = pass0["y"] + pass1["y"]

    full = cast(
        "dict[str, Any]",
        solve_modal(
            spec,
            grid,
            y0,
            t_span=(0.0, T_END),
            parameters={"m2": M2, "eps": eps},
            num_snapshots=N_SNAP,
        ),
    )
    q_full = full["y"]

    err_pass0 = float(np.max(np.abs(q_pass0[-1, :N_GRID] - q_full[-1, :N_GRID])))
    err_combined = float(np.max(np.abs(q_total[-1, :N_GRID] - q_full[-1, :N_GRID])))
    return {
        "eps": eps,
        "err_pass0_vs_full": err_pass0,
        "err_combined_vs_full": err_combined,
        "improvement_ratio": (err_pass0 / err_combined)
        if err_combined > 0
        else float("inf"),
        "scaled_err_combined_over_eps2": err_combined / max(eps * eps, 1e-30),
    }


def run(*, smoke: bool) -> dict:
    eps_values = SMOKE_EPS_VALUES if smoke else FULL_EPS_VALUES
    rows: list[dict] = []
    for eps in eps_values:
        print(f"[parker_simon_flrw] eps={eps}", flush=True)
        rows.append(_run_eps(eps))

    eps_arr = np.array([r["eps"] for r in rows])
    err_arr = np.array([r["err_combined_vs_full"] for r in rows])
    mask = (eps_arr > 0) & (err_arr > 1e-15) & np.isfinite(err_arr)
    slope = None
    if mask.sum() >= 2:
        slope, _ = np.polyfit(np.log10(eps_arr[mask]), np.log10(err_arr[mask]), 1)
        slope = float(slope)

    summary = {
        "n_eps_points": len(rows),
        "fitted_slope_err_vs_eps": slope,
        "expected_slope": 2.0,
        "max_err_combined": float(max(r["err_combined_vs_full"] for r in rows)),
        "min_err_combined": float(min(r["err_combined_vs_full"] for r in rows)),
        "median_improvement_ratio": float(
            np.median([r["improvement_ratio"] for r in rows])
        ),
        "correction_operator": "biharmonic",
        "reference": "Parker & Simon 1993, gr-qc/9211002 (canonical higher-derivative perturbative-reduction structure)",
    }
    return {
        "metadata": _metadata(
            {
                "m2": M2,
                "t_end": T_END,
                "n_grid": N_GRID,
                "k_mode": K_MODE,
                "n_snapshots": N_SNAP,
                "domain": [0.0, 2 * np.pi],
                "eps_values": eps_values,
                "smoke": smoke,
                "correction": "ε · ∂_x^4 (phi) [derivative_4_x]",
                "full_dispersion": "omega² = k² + m² + ε k⁴",
            }
        ),
        "summary": summary,
        "results": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    data = run(smoke=args.smoke)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
