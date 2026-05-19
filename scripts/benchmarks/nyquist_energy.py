"""Nyquist-mode energy drift with and without the IC zeroing fix.

Sweeps N ∈ {128, 256, 512, 1024} for ``coupled_scalars`` (massless) and
``coupled_scalars_massive`` (omegaP2=1), measuring the relative energy drift
|dE/E| under
the Fourier modal solver with ``_zero_nyquist=True`` (production default)
and ``_zero_nyquist=False`` (reproduces the pre-fix behaviour).

The IC is a broad Gaussian plus a small Nyquist component (amplitude
``NYQ_AMP``), so that the Nyquist energy fraction is ~1e-5 of total.
With the fix the drift falls to the floating-point floor; without it the
drift is proportional to the Nyquist energy fraction and is
resolution-independent.

Serves:   manuscript/sections/appendices/numerical.tex (tab:NyquistEnergy)
Consumes: scripts/tables/tab_nyquist_energy.py
Writes:   benchmark_results/canonical/nyquist_energy.json by default
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import socket
import subprocess  # noqa: S404
import sys
from pathlib import Path

import numpy as np
import scipy

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.measurement._energy import compute_energy_timeseries
from tidal.measurement._io import SimulationData
from tidal.solver.grid import GridInfo
from tidal.solver.modal import solve_modal
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

THEORIES: list[tuple[str, Path, dict]] = [
    (
        "coupled_scalars",
        REPO_ROOT / "examples" / "data" / "coupled_scalars.json",
        {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
    ),
    (
        # Same JSON, with plasma mass omegaP2=1 — different dispersion (ω²=k²+m²).
        # The gertsenshtein theory has growing modes (Re λ=k_nyq) in the modal
        # basis (gauge-field ghost mode), so it is not modal-solver compatible.
        "coupled_scalars_massive",
        REPO_ROOT / "examples" / "data" / "coupled_scalars.json",
        {"kappa": 1.0, "B0": 0.1, "omegaP2": 1.0, "mg2": 0.0},
    ),
]
N_VALUES = [128, 256, 512, 1024, 2048, 4096]
DOMAIN = (0.0, 100.0)  # matches examples/*/run.sh --bounds 0:100
T_END = 10.0
NUM_SNAPSHOTS = 51

# Gaussian IC parameters
GAUSS_AMPLITUDE = 0.1
GAUSS_WIDTH = 10.0  # broad Gaussian on [0,100] domain — minimal natural Nyquist content

# Nyquist component amplitude relative to GAUSS_AMPLITUDE.
# For energy fraction ~1e-5 in Nyquist mode, eps = NYQ_AMP * GAUSS_AMPLITUDE:
#   Gaussian mean-sq ≈ A^2 * σ*sqrt(π)/L ≈ 0.01 * 10*1.77/100 ≈ 1.8e-3
#   Nyquist mean-sq = (eps)^2 / 2; fraction = eps^2 / (2 * 1.8e-3)
#   For fraction ~1e-5: eps ≈ sqrt(1e-5 * 2 * 1.8e-3) ≈ 1.9e-4 → NYQ_AMP ≈ 1.9e-3
NYQ_AMP = 2e-3  # relative to GAUSS_AMPLITUDE → eps ≈ 2e-4; Nyquist fraction ~1e-5


def _git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            stderr=subprocess.DEVNULL,
        )
        return out.decode().strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def _metadata() -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "parameters": {
            "theories": [name for name, _, _ in THEORIES],
            "n_values": N_VALUES,
            "domain": list(DOMAIN),
            "t_end": T_END,
            "gauss_amplitude": GAUSS_AMPLITUDE,
            "gauss_width": GAUSS_WIDTH,
            "nyq_relative_amplitude": NYQ_AMP,
        },
    }


def _build_ic(spec: EquationSystem, grid: GridInfo) -> np.ndarray:
    """Gaussian on the first field plus a small pure-Nyquist component."""
    layout = StateLayout.from_spec(spec, grid.num_points)
    y0 = np.zeros(layout.num_slots * grid.num_points)
    N = grid.num_points
    x = np.linspace(DOMAIN[0], DOMAIN[1], N, endpoint=False)
    center = (DOMAIN[0] + DOMAIN[1]) / 2.0

    # Broad Gaussian — minimal natural Nyquist content
    gaussian = GAUSS_AMPLITUDE * np.exp(-((x - center) ** 2) / (2.0 * GAUSS_WIDTH**2))

    # Pure Nyquist grid mode: (-1)^n alternating pattern.
    # In continuous form this is cos(π·N/L·x); on the grid x[n]=n·L/N it
    # reduces to cos(π·n) = (-1)^n.  The rfft of (-1)^n has energy only in
    # the Nyquist bin (k=N//2), which is what we want to test.
    nyquist = NYQ_AMP * GAUSS_AMPLITUDE * (-1) ** np.arange(N)

    y0[:N] = gaussian + nyquist
    return y0


def _energy_drift(
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    parameters: dict,
    zero_nyquist: bool,
) -> float:
    """Return |ΔE/E| = |E(t_end) - E(0)| / |E(0)| for the given IC and flag."""
    result = solve_modal(
        spec,
        grid,
        y0,
        (0.0, T_END),
        parameters=parameters,
        num_snapshots=NUM_SNAPSHOTS,
        _zero_nyquist=zero_nyquist,
    )
    data = SimulationData.from_result(result, spec, grid, parameters=parameters)
    _times, _per_field, _interaction, total = compute_energy_timeseries(data)
    e0 = float(total[0])
    if abs(e0) < 1e-30:
        return float("nan")
    return float(abs(total[-1] - e0) / abs(e0))


# ----------------------------------------------------------------------------
# Module-level callable for the HPC master driver (`run_all_hpc.py`).
#
# `run_one_nyquist` is a single (theory, n, zero_nyquist) measurement.
# The Nyquist drift is deterministic so we only need one call per
# (theory, n, projection) combination — no rep loop, no error bars.
# ----------------------------------------------------------------------------


def run_one_nyquist(
    theory_name: str,
    json_path: Path,
    n: int,
    parameters: dict,
    zero_nyquist: bool,
) -> dict:
    """Run one modal solve and return its energy drift."""
    with json_path.open(encoding="utf-8") as fh:
        spec = EquationSystem.from_dict(json.load(fh))
    grid = GridInfo(bounds=(DOMAIN,), shape=(n,), periodic=(True,))
    y0 = _build_ic(spec, grid)
    drift = _energy_drift(spec, grid, y0, parameters, zero_nyquist=zero_nyquist)
    return {
        "theory": theory_name,
        "n": n,
        "zero_nyquist": zero_nyquist,
        "abs_dE_over_E": drift,
    }


def run() -> dict:
    """Standalone-CLI path: dev-container smoke test only.

    Canonical-quality data must come from
    ``scripts/benchmarks/run_all_hpc.py`` on a CSD3 sapphire compute node.
    """
    rows: list[dict] = []
    for theory_name, json_path, params in THEORIES:
        for n in N_VALUES:
            print(f"  {theory_name}  N={n}", flush=True)
            r_fix = run_one_nyquist(
                theory_name, json_path, n, params, zero_nyquist=True
            )
            r_nofix = run_one_nyquist(
                theory_name, json_path, n, params, zero_nyquist=False
            )
            drift_fixed = r_fix["abs_dE_over_E"]
            drift_nofixed = r_nofix["abs_dE_over_E"]
            rows.append(
                {
                    "theory": theory_name,
                    "n": n,
                    "abs_dE_over_E_fixed": drift_fixed,
                    "abs_dE_over_E_nofixed": drift_nofixed,
                }
            )
            print(
                f"    with fix: {drift_fixed:.2e}   without fix: {drift_nofixed:.2e}",
                flush=True,
            )

    return {"metadata": _metadata(), "results": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "benchmark_results" / "canonical" / "nyquist_energy.json",
    )
    args = parser.parse_args()

    print("Running Nyquist energy drift benchmark...")
    data = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
