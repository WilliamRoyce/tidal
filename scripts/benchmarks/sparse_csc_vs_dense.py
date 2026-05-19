"""Sparse-CSC vs dense ``expm_multiply`` wall-time on localised Gertsenshtein.

Serves:   manuscript/sections/appendices/numerical.tex (tab:SparseCSCBenchmark)
Consumes: scripts/tables/tab_sparse_csc.py
Writes:   benchmark_results/canonical/sparse_csc_vs_dense.json by default

Builds the full convolution evolution matrix for the localised
Gertsenshtein system (Gaussian B-field background), then times
``scipy.sparse.linalg.expm_multiply`` in dense and CSC formats at
N ∈ {128, 256, 512}.  Density is reported alongside the wall times so
the reader can verify when the heuristic threshold is met.
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
import scipy.sparse
import scipy.sparse.linalg

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.modal import (
    _build_convolution_matrix,
    _build_k_axes,
    _build_k_grid,
)
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

THEORY_JSON = REPO_ROOT / "examples" / "data" / "gertsenshtein_localized.json"
PARAMS = {"kappa": 1.0, "Bpeak": 0.1, "R": 5.0}
DOMAIN = (-50.0, 50.0)  # background Gaussian centred at 0, width R=5

# Dense `expm_multiply` becomes infeasible above N≈1024 (matrix exceeds
# 300 MB at N=2048). N=4096 is capped because building the convolution
# matrix needs the dense intermediate (~9.7 GB) even for the sparse-only
# path. Sparse-CSC reaches N=2048 here; production goes higher because
# the matrix is built mode-by-mode without ever fully materialising.
N_DENSE = [64, 128, 256, 512, 1024]
N_SPARSE = [64, 128, 256, 512, 1024, 2048]

DT = 0.2  # matches T_END=10, num_snapshots=51
THRESHOLD_REL = 1e-14  # match production threshold in _evolve_full_matrix
N_WARMUP = 2
# 10 reps is enough for stable mean/std on the cheap N values; dense
# N=1024 (~10 s/rep) makes the total acceptable (~100 s).
N_REPS = 10


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
            "theory": "gertsenshtein_localized",
            "n_dense": N_DENSE,
            "n_sparse": N_SPARSE,
            "domain": list(DOMAIN),
            "params": PARAMS,
            "dt": DT,
            "threshold_relative": THRESHOLD_REL,
            "n_reps": N_REPS,
        },
    }


# ----------------------------------------------------------------------------
# Module-level callables for the HPC master driver (`run_all_hpc.py`).
#
# `build_sparse_setup` builds the convolution matrix and both its dense
# (`A_dense_dt = A * DT`) and CSC (`A_csc`) representations.  Memory: the
# dense matrix is the dominant cost (~9.7 GB at N=4096 complex128).  The
# master driver caches one setup per (n, path) in each worker so multiple
# timed reps reuse the build.
#
# `time_one_sparse_rep(setup)` runs ONE `expm_multiply` call and returns
# its wall time.
# ----------------------------------------------------------------------------


def build_sparse_setup(n: int, path: str) -> dict:
    """Build matrix + IC for one (n, path).

    path ∈ {"dense", "sparse"} — only the requested representation is
    materialised, so memory cost matches the actual path being timed.
    """
    with THEORY_JSON.open(encoding="utf-8") as fh:
        spec = EquationSystem.from_dict(json.load(fh))

    grid = GridInfo(bounds=(DOMAIN,), shape=(n,), periodic=(True,))
    layout = StateLayout.from_spec(spec, n)
    coeff_eval = CoefficientEvaluator(spec, grid, PARAMS)
    k_axes = _build_k_axes(grid)
    k_grid = _build_k_grid(k_axes)
    rfft_shape = (n // 2 + 1,)

    A_full = _build_convolution_matrix(
        spec, layout, grid, coeff_eval, k_grid, rfft_shape
    )
    threshold = float(np.max(np.abs(A_full))) * THRESHOLD_REL
    A_full[np.abs(A_full) < threshold] = 0
    density = float(np.count_nonzero(A_full)) / A_full.size
    matrix_size = A_full.shape[0]

    rng = np.random.default_rng(seed=42)
    y0 = (
        rng.standard_normal(matrix_size) + 1j * rng.standard_normal(matrix_size)
    ).astype(np.complex128)

    if path == "dense":
        A_op = A_full * DT
    elif path == "sparse":
        A_op = scipy.sparse.csc_array(A_full * DT)
        # Release the dense original (it's the bulk of the memory).
        A_full = None
    else:
        msg = f"Unknown sparse_csc path: {path!r}"
        raise ValueError(msg)

    return {
        "n": n,
        "path": path,
        "matrix_size": matrix_size,
        "density_pct": density * 100.0,
        "A_op": A_op,
        "y0": y0,
        "note": "above_density_cutoff" if density >= 0.30 else "below_density_cutoff",
    }


def time_one_sparse_rep(setup: dict, batch_size: int = 1) -> dict:
    """One timed ``expm_multiply`` measurement on a pre-built setup.

    With ``batch_size > 1`` the call is repeated ``batch_size`` times before
    averaging; this is the adaptive-batching path used to push per-batch
    wall time above the OS-jitter floor on small-N matrices.
    """
    A_op = setup["A_op"]
    y0 = setup["y0"]
    t0 = time.perf_counter()
    for _ in range(batch_size):
        scipy.sparse.linalg.expm_multiply(A_op, y0)
    return {"wall_s": (time.perf_counter() - t0) / batch_size}


def probe_sparse_batch_size(setup: dict, target_s: float = 0.1) -> int:
    """Pick a batch size so one measurement runs for ``target_s`` or longer."""
    t0 = time.perf_counter()
    scipy.sparse.linalg.expm_multiply(setup["A_op"], setup["y0"])
    probe_s = time.perf_counter() - t0
    if probe_s >= target_s:
        return 1
    return max(1, int(np.ceil(target_s / probe_s)))


def warmup_sparse(setup: dict, n_warmup: int = N_WARMUP) -> None:
    for _ in range(n_warmup):
        scipy.sparse.linalg.expm_multiply(setup["A_op"], setup["y0"])


def _time_one(n: int, *, time_dense: bool, time_sparse: bool) -> list[dict]:
    """Standalone-CLI path: time requested paths for one N. Smoke-test only.

    Canonical-quality timing must come from
    ``scripts/benchmarks/run_all_hpc.py`` on a CSD3 sapphire compute node.
    """
    rows: list[dict] = []

    def _run_path(path: str) -> dict:
        setup = build_sparse_setup(n, path)
        print(
            f"  N={n}  matrix {setup['matrix_size']}x{setup['matrix_size']}  "
            f"density={setup['density_pct']:.2f}%  path={path}",
            flush=True,
        )
        warmup_sparse(setup)
        batch_size = probe_sparse_batch_size(setup)
        reps_s = [
            time_one_sparse_rep(setup, batch_size=batch_size)["wall_s"]
            for _ in range(N_REPS)
        ]
        mean = float(np.mean(reps_s))
        std = float(np.std(reps_s, ddof=1))
        print(
            f"    {path}: {mean:.4g}±{std:.2g} s  ({setup['note']})",
            flush=True,
        )
        return {
            "n": n,
            "path": path,
            "matrix_size": setup["matrix_size"],
            "density_pct": setup["density_pct"],
            "wall_s": mean,
            "wall_s_mean": mean,
            "wall_s_std": std,
            "wall_s_reps": reps_s,
            "n_reps": N_REPS,
            "note": setup["note"],
        }

    if time_dense:
        rows.append(_run_path("dense"))
    if time_sparse:
        rows.append(_run_path("sparse"))

    return rows


def run() -> dict:
    rows: list[dict] = []
    n_all = sorted(set(N_DENSE) | set(N_SPARSE))
    for n in n_all:
        rows.extend(
            _time_one(n, time_dense=(n in N_DENSE), time_sparse=(n in N_SPARSE))
        )
    return {"metadata": _metadata(), "results": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT
        / "benchmark_results"
        / "canonical"
        / "sparse_csc_vs_dense.json",
    )
    args = parser.parse_args()

    print("Running sparse-CSC vs dense expm_multiply benchmark...")
    data = run()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
