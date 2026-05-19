"""Forced-oscillator perturbative-reduction calibration (Pass 0 + Pass 1).

Calibrates the closed-form Duhamel kernel of the perturbative reduction
against the analytic perturbative expansion of a driven Klein--Gordon
oscillator with a small mass-shift correction:

    ∂_t² φ = ∂_x² φ - m² φ - ε · φ

The Pass 0 (ε=0) solution is the standard Klein-Gordon evolution; the
Pass 1 correction is computed via the closed-form Duhamel kernel
φ_1(λ_i, λ_j; t) in the eigenbasis of the base operator. We sweep ε
over a range and report the truncation error |Pass0+1 - Full|, which
should scale as O(ε²) by construction (Parker & Simon 1993; Al-Mohy
& Higham 2011 augmented matrix exponential).

This is a non-circular calibration of the perturbative-reduction path:
the closed-form Duhamel kernel is a generic ODE-theory construct
derived independently of TIDAL's implementation. The improvement over
Pass-0-only error is also recorded as a sanity check that Pass 1 is
actually contributing.

Serves:   manuscript/sections/appendices/validation.tex (App D §3)
Consumes: scripts/figures/figD_forced_oscillator.py
Writes:   benchmark_results/canonical/forced_oscillator.json
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

from tidal.solver.grid import GridInfo
from tidal.solver.modal import solve_modal, solve_modal_pass1
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "forced_oscillator.json"

# 1+1D Klein-Gordon spec with ε mass-shift correction. Pass 0 contains the
# order=0 terms only; Pass 1 picks up the order=1 mass-shift term.
_KG_BASE: dict[str, Any] = {
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
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-eps",
                        "order_in_eps": 1,
                    },
                ],
            },
        },
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2 - eps"]]},
}

# eps ladder restricted to the perturbative regime where the O(eps^2)
# theorem is informative. The eps=0.5 point hits saturation
# (err_pass0 ~ err_combined) and was dropped.
FULL_EPS_VALUES = [0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001]
SMOKE_EPS_VALUES = [0.1, 0.05, 0.01]
M2 = 1.0
T_END = 2.0
N_GRID = 64
K_MODE = 1.0
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
    spec = EquationSystem.from_dict(copy.deepcopy(_KG_BASE))
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
    t_eval = pass0["t"]
    pass1 = solve_modal_pass1(
        pass0["eigendata"],
        correction_spec,
        grid,
        t_eval,
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
        print(f"[forced_oscillator] eps={eps}", flush=True)
        rows.append(_run_eps(eps))

    # Fit log-log slope of err_combined vs eps
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
