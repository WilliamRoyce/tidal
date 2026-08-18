"""Padé scaling-and-squaring vs eigendecomposition wall-time across theories.

Serves:   manuscript/sections/appendices/numerical.tex (tab:PadeBenchmark)
Consumes: scripts/tables/tab_pade_benchmark.py
Writes:   benchmark_results/canonical/pade_vs_eig.json by default

Times scipy.linalg.expm (Padé path, production default) against the
retired eigendecomposition path V @ diag(exp(D·t)) @ V_inv on the same
per-mode matrices built from three representative theories:

  coupled_scalars    N=256 1D  — 4-slot per-mode blocks,  128 modes
  gertsenshtein      N=512 1D  — 12-slot per-mode blocks, 256 modes
  massive_gravity_3d N=32  2D  — 12-slot per-mode blocks, 544 modes

For each theory the per-mode matrices are built via the modal solver's
internal _build_per_mode_matrices() and the per-mode wall time is
averaged over all active modes and N_REPS repetitions.
"""

from __future__ import annotations

import argparse
import datetime
import json
import platform
import socket
import subprocess  # noqa: S404
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.linalg

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_k_axes,
    _build_k_grid,
    _build_per_mode_matrices,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# (theory_name, json_path, grid_shapes, bounds, params, n_labels)
# grid_shapes is a list of shapes to sweep; n_labels are the matching display
# labels (LaTeX-safe; "32^2" renders as $32^2$ via the figure script).
THEORIES: list[
    tuple[
        str,
        Path,
        list[tuple[int, ...]],
        tuple[tuple[float, float], ...],
        dict,
        list[str],
    ]
] = [
    (
        "coupled_scalars",
        REPO_ROOT / "examples" / "data" / "coupled_scalars.json",
        [(128,), (256,), (512,), (1024,), (2048,), (4096,)],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
        ["128", "256", "512", "1024", "2048", "4096"],
    ),
    (
        "gertsenshtein",
        REPO_ROOT / "examples" / "data" / "gertsenshtein.json",
        [(128,), (256,), (512,), (1024,), (2048,), (4096,)],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.3},
        ["128", "256", "512", "1024", "2048", "4096"],
    ),
    (
        "massive_gravity_3d",
        REPO_ROOT / "examples" / "data" / "massive_gravity_3d.json",
        [(16, 16), (24, 24), (32, 32), (48, 48), (64, 64)],
        ((0.0, 50.0), (0.0, 50.0)),
        {"m2": 1.0},
        ["16^2", "24^2", "32^2", "48^2", "64^2"],
    ),
]

DT = 0.2  # timestep for matrix exponential (matches T_END=10, num_snapshots=51)
N_WARMUP = 3  # warm-up calls before timing
N_REPS = 20  # timing repetitions (averaged)


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
            "theories": [name for name, *_ in THEORIES],
            "dt": DT,
            "n_warmup": N_WARMUP,
            "n_reps": N_REPS,
        },
    }


def _eig_expm(M: np.ndarray, dt: float) -> np.ndarray:
    """Eigendecomposition path: V @ diag(exp(d·dt)) @ V_inv."""
    d, V = np.linalg.eig(M)
    V_inv = np.linalg.inv(V)
    return V @ (np.exp(d * dt)[:, None] * V_inv)


# ----------------------------------------------------------------------------
# Module-level callables for the HPC master driver (`run_all_hpc.py`).
#
# `build_pade_setup` builds the per-mode matrices for one (theory, n_label).
# This is the expensive O(n_modes × n_slots²) construction cost; the master
# driver caches the result per (theory, n_label) inside each worker so
# multiple timed reps amortize the build.
#
# `time_one_pade_rep` times exactly one (Padé, Eig) pair on a pre-built
# setup. The master driver invokes this once per (theory, n_label, rep_idx).
# ----------------------------------------------------------------------------


def build_pade_setup(
    theory_name: str,
    json_path: Path,
    grid_shape: tuple[int, ...],
    bounds: tuple[tuple[float, float], ...],
    params: dict,
) -> dict:
    """Build per-mode matrices for one (theory, grid_shape). Heavy."""
    with json_path.open(encoding="utf-8") as fh:
        spec = EquationSystem.from_dict(json.load(fh))

    grid = GridInfo(
        bounds=bounds, shape=grid_shape, periodic=tuple(True for _ in grid_shape)
    )
    layout = StateLayout.from_spec(spec, grid.num_points)
    coeff_eval = CoefficientEvaluator(spec, grid, params)
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape = (*grid_shape[:-1], grid_shape[-1] // 2 + 1)

    A_modes = _build_per_mode_matrices(
        spec, layout, grid, coeff_eval, k_grid, rfft_shape
    )
    raw = [A_modes[m] for m in range(A_modes.shape[0])]
    return {
        "theory": theory_name,
        "raw": raw,
        "n_modes": A_modes.shape[0],
        "n_slots": A_modes.shape[1],
    }


def time_one_pade_rep(
    setup: dict,
    dt: float = DT,
    batch_size: int = 1,
) -> dict:
    """One timed rep on a pre-built setup. Returns per-mode (Padé, Eig) times.

    With ``batch_size > 1`` the inner timing loop is repeated ``batch_size``
    times before averaging; this is the adaptive-batching path used to push
    per-batch wall time above the OS-jitter floor on small-N matrices.
    """
    raw = setup["raw"]
    n_modes = setup["n_modes"]
    denom = n_modes * batch_size
    t0 = time.perf_counter()
    for _ in range(batch_size):
        for M in raw:
            scipy.linalg.expm(M * dt)
    pade_s = (time.perf_counter() - t0) / denom
    t0 = time.perf_counter()
    for _ in range(batch_size):
        for M in raw:
            _eig_expm(M, dt)
    eig_s = (time.perf_counter() - t0) / denom
    return {"pade_s_per_mode": pade_s, "eig_s_per_mode": eig_s}


def probe_pade_batch_size(setup: dict, target_s: float = 0.1, dt: float = DT) -> int:
    """Pick a batch size so one (Padé+Eig) measurement runs for ``target_s``
    or longer. Brings sub-millisecond timings above the OS-jitter floor.
    """
    raw = setup["raw"]
    t0 = time.perf_counter()
    for M in raw:
        scipy.linalg.expm(M * dt)
    probe_s = time.perf_counter() - t0
    if probe_s >= target_s:
        return 1
    return max(1, int(np.ceil(target_s / probe_s)))


def warmup_pade(setup: dict, n_warmup: int = N_WARMUP, dt: float = DT) -> None:
    """Untimed warmup runs to settle BLAS / branch prediction / page tables."""
    raw = setup["raw"]
    for _ in range(n_warmup):
        for M in raw:
            scipy.linalg.expm(M * dt)
            _eig_expm(M, dt)


def _time_theory(
    theory_name: str,
    json_path: Path,
    grid_shape: tuple[int, ...],
    bounds: tuple[tuple[float, float], ...],
    params: dict,
    n_label: str,
) -> dict:
    """Standalone-CLI path: build matrices once + run N_REPS sequentially.

    This is the dev-container smoke-test path. Canonical-quality timing
    must be produced by ``scripts/benchmarks/run_all_hpc.py`` on a CSD3
    sapphire compute node — see plan file.
    """
    setup = build_pade_setup(theory_name, json_path, grid_shape, bounds, params)
    n_modes = setup["n_modes"]
    n_slots = setup["n_slots"]
    print(
        f"  {theory_name}  N={n_label}  modes={n_modes}  slots={n_slots}x{n_slots}",
        flush=True,
    )

    warmup_pade(setup)
    batch_size = probe_pade_batch_size(setup)
    pade_per_mode_reps: list[float] = []
    eig_per_mode_reps: list[float] = []
    for _ in range(N_REPS):
        rep = time_one_pade_rep(setup, batch_size=batch_size)
        pade_per_mode_reps.append(rep["pade_s_per_mode"])
        eig_per_mode_reps.append(rep["eig_s_per_mode"])

    pade_mean = float(np.mean(pade_per_mode_reps))
    pade_std = float(np.std(pade_per_mode_reps, ddof=1))
    eig_mean = float(np.mean(eig_per_mode_reps))
    eig_std = float(np.std(eig_per_mode_reps, ddof=1))
    ratio = pade_mean / eig_mean

    print(
        f"    Padé: {pade_mean:.4g}±{pade_std:.2g} s/mode  "
        f"Eig: {eig_mean:.4g}±{eig_std:.2g} s/mode  Ratio: {ratio:.2f}",
        flush=True,
    )

    return {
        "theory": theory_name,
        "n": n_label,
        "n_modes": n_modes,
        "n_slots": n_slots,
        "pade_s": pade_mean,
        "eig_s": eig_mean,
        "ratio": ratio,
        "pade_s_mean": pade_mean,
        "pade_s_std": pade_std,
        "pade_s_reps": pade_per_mode_reps,
        "eig_s_mean": eig_mean,
        "eig_s_std": eig_std,
        "eig_s_reps": eig_per_mode_reps,
    }


def run() -> dict:
    rows: list[dict] = []
    for theory_name, json_path, grid_shapes, bounds, params, n_labels in THEORIES:
        for grid_shape, n_label in zip(grid_shapes, n_labels, strict=True):
            row = _time_theory(
                theory_name, json_path, grid_shape, bounds, params, n_label
            )
            rows.append(row)
    return {"metadata": _metadata(), "results": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "benchmark_results" / "canonical" / "pade_vs_eig.json",
    )
    args = parser.parse_args()

    print("Running Padé vs eigendecomposition benchmark...")
    data = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
