r"""Per-Newton-iteration wall time: analytical Jacobian tiers vs FD proxy.

Measures all three analytical-Jacobian delivery tiers (dense ``jacfn``,
sparse-CSC ``jacfn``, GMRES ``jactimes``) at the **same** n_total values across
three representative theories, bypassing the production auto-selection in
``try_analytical_jacobian()``.  The callbacks are called directly after
``build_jacobian_matrices()`` — no IDA/CVODE solver instantiation required.

The FD proxy cost is estimated as ``n_colors × one_matvec``, where
``n_colors = max_nnz_per_column`` of the union sparsity pattern (an upper
bound on the chromatic number of the column-intersection graph, exact for
banded systems).  ``one_matvec`` is the time to compute
``J_union @ y`` for a random vector ``y``, which equals the analytical cost of
one IDA residual evaluation for a linear time-independent system.

All callbacks are timed with adaptive batching (≥100 ms per batch) and
N_REPS=20 outer repetitions.  Trials are shuffled across all (theory, N_grid)
configurations per rep to average out monotonic drift.

Usage on the compute node, after ``hpc_shuttle.sh attach <jobid>``:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
    uv run python scripts/benchmarks/jac_speedup.py \\
        --reps 20 --seed 42 --out benchmark_results/canonical/jac_speedup.json

Output is written to benchmark_results/canonical/jac_speedup.json by default.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import platform
import random
import socket
import subprocess  # noqa: S404
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver._types import DENSE_THRESHOLD, SPARSE_THRESHOLD  # noqa: E402
from tidal.solver.analytical_jacobian import (  # noqa: E402
    _create_ida_jactimes,
    _create_jacfn,
    _create_sparse_jacfn,
    build_jacobian_matrices,
)
from tidal.solver.grid import GridInfo  # noqa: E402
from tidal.solver.state import StateLayout  # noqa: E402
from tidal.symbolic.json_loader import EquationSystem  # noqa: E402

# ---------------------------------------------------------------------------
# Benchmark configuration
# ---------------------------------------------------------------------------

# Dense callback is skipped when n_total exceeds this — allocating two n×n
# float64 matrices costs 2×n²×8 bytes (40 GB at n=50K, fine on sapphire).
DENSE_BENCH_LIMIT: int = 50_000

N_REPS: int = 20
N_WARMUP: int = 3
_PROBE_TARGET_S: float = 0.1  # adaptive batch target (seconds)

# (theory_name, json_path, grid_shapes, bounds, params)
# grid_shapes is a list of shapes to sweep; for 2-D theories each entry is
# (N, N).  n_total is computed at runtime from layout.total_size.
THEORIES: list[
    tuple[
        str,
        Path,
        list[tuple[int, ...]],
        tuple[tuple[float, float], ...],
        dict,
    ]
] = [
    (
        "coupled_scalars",
        REPO_ROOT / "examples" / "data" / "coupled_scalars.json",
        [
            (64,),
            (128,),
            (256,),
            (512,),
            (1024,),
            (2048,),
            (4096,),
            (8192,),
            (16384,),
            (32768,),
            (65536,),
            (131072,),
        ],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
    ),
    (
        "gertsenshtein",
        REPO_ROOT / "examples" / "data" / "gertsenshtein.json",
        [
            (16,),
            (32,),
            (64,),
            (128,),
            (256,),
            (512,),
            (1024,),
            (2048,),
            (4096,),
            (8192,),
            (16384,),
            (32768,),
        ],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.3},
    ),
    (
        # massive_gravity_3d has d2_t in constraint equations — incompatible
        # with the analytical Jacobian builder.  Use navier_cauchy_2d instead:
        # 2 displacement fields (+ 2 velocity slots) on an N×N grid give
        # n_total = 4N², the same O(N²) 2-D scaling the figure aims to show.
        "navier_cauchy_2d",
        REPO_ROOT / "examples" / "data" / "navier_cauchy_2d.json",
        [
            (8, 8),
            (12, 12),
            (16, 16),
            (24, 24),
            (32, 32),
            (48, 48),
            (64, 64),
            (96, 96),
            (128, 128),
            (160, 160),
            (192, 192),
            (256, 256),
        ],
        ((0.0, 1.0), (0.0, 1.0)),
        {"rho": 1.0, "lam": 1.0, "mu": 1.0},
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


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


def _metadata(reps: int, seed: int) -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "parameters": {
            "theories": [t[0] for t in THEORIES],
            "n_reps": reps,
            "seed": seed,
            "dense_bench_limit": DENSE_BENCH_LIMIT,
            "dense_threshold": DENSE_THRESHOLD,
            "sparse_threshold": SPARSE_THRESHOLD,
        },
    }


def _n_colors_estimate(J_csc: scipy.sparse.csc_matrix) -> int:
    """Upper bound on FD coloring count: max nonzeros per column.

    Equals the exact chromatic number for banded matrices (tight for FD
    stencils).  For general sparse patterns it is a conservative upper bound.
    """
    nnz_per_col = np.diff(J_csc.indptr)
    return int(nnz_per_col.max()) if J_csc.shape[1] > 0 else 1


# ---------------------------------------------------------------------------
# Setup (heavy — one-time per (theory, N_grid))
# ---------------------------------------------------------------------------


def build_jac_setup(
    theory_name: str,
    json_path: Path,
    grid_shape: tuple[int, ...],
    bounds: tuple[tuple[float, float], ...],
    params: dict,
) -> dict:
    """Build all callbacks for one (theory, N_grid) point.

    Returns a setup dict consumed by ``time_one_jac_rep`` and ``warmup_jac``.
    The dense callback is None when n_total > DENSE_BENCH_LIMIT.
    """
    with json_path.open(encoding="utf-8") as fh:
        spec = EquationSystem.from_dict(json.load(fh))

    grid = GridInfo(
        bounds=bounds,
        shape=grid_shape,
        periodic=tuple(True for _ in grid_shape),
    )
    layout = StateLayout.from_spec(spec, grid.num_points)
    bc = None  # periodic — GridInfo.periodic propagates periodicity

    dF_dy, dF_dyp = build_jacobian_matrices(spec, layout, grid, bc, params)
    n_total = layout.total_size

    rng = np.random.default_rng(42)
    y0 = np.zeros(n_total)
    yp0 = np.zeros(n_total)
    v_rand = rng.standard_normal(n_total)
    jv = np.empty(n_total)

    # --- Dense analytical callback ---
    if n_total <= DENSE_BENCH_LIMIT:
        jacfn_dense = _create_jacfn(dF_dy, dF_dyp)
        JJ_dense = np.empty((n_total, n_total))

        def call_dense() -> None:
            jacfn_dense(0.0, y0, yp0, y0, 1.0, JJ_dense)

    else:
        call_dense = None  # type: ignore[assignment]

    # --- Sparse analytical callback ---
    jacfn_sparse, sparsity = _create_sparse_jacfn(dF_dy, dF_dyp)
    JJ_sparse = np.empty(sparsity.nnz)

    def call_sparse() -> None:
        jacfn_sparse(0.0, y0, yp0, y0, 1.0, JJ_sparse)

    # --- GMRES jactimes callback ---
    solvefn_gmres = _create_ida_jactimes(dF_dy, dF_dyp)

    def call_gmres() -> None:
        solvefn_gmres(0.0, y0, yp0, y0, v_rand, jv, 1.0)

    # --- FD proxy: n_colors × one sparse mat-vec ---
    J_union = (abs(dF_dy) + abs(dF_dyp)).tocsc()
    J_union.eliminate_zeros()
    y_fd = rng.standard_normal(n_total)
    out_fd = np.empty(n_total)
    n_colors = _n_colors_estimate(J_union)

    def call_matvec() -> None:
        out_fd[:] = J_union @ y_fd

    return {
        "theory": theory_name,
        "grid_shape": grid_shape,
        "n_total": n_total,
        "n_colors": n_colors,
        "nnz": int(J_union.nnz),
        "call_dense": call_dense,
        "call_sparse": call_sparse,
        "call_gmres": call_gmres,
        "call_matvec": call_matvec,
    }


# ---------------------------------------------------------------------------
# Adaptive batch sizing and warmup
# ---------------------------------------------------------------------------


def probe_jac_batch_size(setup: dict, target_s: float = _PROBE_TARGET_S) -> int:
    """Pick a batch size so one sparse-callback measurement ≥ target_s."""
    t0 = time.perf_counter()
    setup["call_sparse"]()
    probe_s = time.perf_counter() - t0
    if probe_s >= target_s:
        return 1
    return max(1, math.ceil(target_s / probe_s))


def warmup_jac(setup: dict, n_warmup: int = N_WARMUP) -> None:
    """Untimed warmup to settle caches and branch prediction."""
    for _ in range(n_warmup):
        if setup["call_dense"] is not None:
            setup["call_dense"]()
        setup["call_sparse"]()
        setup["call_gmres"]()
        setup["call_matvec"]()


# ---------------------------------------------------------------------------
# Per-rep timing
# ---------------------------------------------------------------------------


def time_one_jac_rep(setup: dict, batch_size: int = 1) -> dict:
    """Time one rep of all four paths. Returns per-call wall times in seconds.

    The FD proxy is timed as ``n_colors`` sequential mat-vec calls (simulating
    the n_colors residual evaluations that SUNDIALS performs per Jacobian fill).
    """
    results: dict[str, float] = {}
    n_colors = setup["n_colors"]

    # Dense
    if setup["call_dense"] is not None:
        t0 = time.perf_counter()
        for _ in range(batch_size):
            setup["call_dense"]()
        results["dense_s"] = (time.perf_counter() - t0) / batch_size
    else:
        results["dense_s"] = float("nan")

    # Sparse
    t0 = time.perf_counter()
    for _ in range(batch_size):
        setup["call_sparse"]()
    results["sparse_s"] = (time.perf_counter() - t0) / batch_size

    # GMRES jactimes (one JVP call)
    t0 = time.perf_counter()
    for _ in range(batch_size):
        setup["call_gmres"]()
    results["gmres_s"] = (time.perf_counter() - t0) / batch_size

    # FD proxy: n_colors mat-vec calls in a loop (= n_colors residual calls)
    t0 = time.perf_counter()
    for _ in range(batch_size):
        for _ in range(n_colors):
            setup["call_matvec"]()
    results["fd_s"] = (time.perf_counter() - t0) / batch_size

    return results


# ---------------------------------------------------------------------------
# Standalone timing for one (theory, N_grid)
# ---------------------------------------------------------------------------


def _time_one_config(
    theory_name: str,
    json_path: Path,
    grid_shape: tuple[int, ...],
    bounds: tuple[tuple[float, float], ...],
    params: dict,
    reps: int,
) -> dict:
    """Build + warm up + time one (theory, N_grid) config. Returns result row."""
    shape_str = "x".join(str(s) for s in grid_shape)
    print(f"  {theory_name}  shape={shape_str}", flush=True, end="")

    setup = build_jac_setup(theory_name, json_path, grid_shape, bounds, params)
    n_total = setup["n_total"]
    n_colors = setup["n_colors"]
    tier = (
        "dense"
        if n_total <= DENSE_THRESHOLD
        else ("sparse" if n_total <= SPARSE_THRESHOLD else "gmres")
    )
    print(
        f"  n_total={n_total}  n_colors={n_colors}  nnz={setup['nnz']}"
        f"  auto_tier={tier}",
        flush=True,
    )

    warmup_jac(setup)
    batch_size = probe_jac_batch_size(setup)

    dense_reps: list[float] = []
    sparse_reps: list[float] = []
    gmres_reps: list[float] = []
    fd_reps: list[float] = []

    for _ in range(reps):
        rep = time_one_jac_rep(setup, batch_size=batch_size)
        dense_reps.append(rep["dense_s"])
        sparse_reps.append(rep["sparse_s"])
        gmres_reps.append(rep["gmres_s"])
        fd_reps.append(rep["fd_s"])

    def _stats(vals: list[float]) -> tuple[float, float]:
        arr = [v for v in vals if not math.isnan(v)]
        if not arr:
            return float("nan"), float("nan")
        return float(np.mean(arr)), float(np.std(arr, ddof=1))

    dense_mean, dense_std = _stats(dense_reps)
    sparse_mean, sparse_std = _stats(sparse_reps)
    gmres_mean, gmres_std = _stats(gmres_reps)
    fd_mean, fd_std = _stats(fd_reps)

    # Speedup: FD / best-analytical (whichever analytical tier is auto-selected)
    auto_mean = {"dense": dense_mean, "sparse": sparse_mean, "gmres": gmres_mean}[tier]
    speedup = (
        fd_mean / auto_mean
        if not math.isnan(auto_mean) and auto_mean > 0
        else float("nan")
    )

    print(
        f"    dense={dense_mean:.3g}s  sparse={sparse_mean:.3g}s"
        f"  gmres={gmres_mean:.3g}s  FD={fd_mean:.3g}s  speedup={speedup:.2f}x",
        flush=True,
    )

    return {
        "theory": theory_name,
        "n_grid": int(np.prod(grid_shape)),
        "grid_shape": list(grid_shape),
        "n_total": n_total,
        "n_colors": n_colors,
        "nnz": setup["nnz"],
        "auto_tier": tier,
        "dense_s_mean": dense_mean,
        "dense_s_std": dense_std,
        "dense_s_reps": dense_reps,
        "sparse_s_mean": sparse_mean,
        "sparse_s_std": sparse_std,
        "sparse_s_reps": sparse_reps,
        "gmres_s_mean": gmres_mean,
        "gmres_s_std": gmres_std,
        "gmres_s_reps": gmres_reps,
        "fd_s_mean": fd_mean,
        "fd_s_std": fd_std,
        "fd_s_reps": fd_reps,
        "speedup_auto_vs_fd": speedup,
    }


# ---------------------------------------------------------------------------
# Standalone orchestrator (shuffled across all configs)
# ---------------------------------------------------------------------------


def run(reps: int = N_REPS, seed: int = 42) -> dict:
    """Run the full benchmark.

    Trials are shuffled across all (theory, N_grid) configurations per rep
    to average out monotonic thermal drift.
    """
    rng = random.Random(seed)  # noqa: S311

    # Build all configs
    configs = []
    for theory_name, json_path, grid_shapes, bounds, params in THEORIES:
        configs.extend(
            (theory_name, json_path, grid_shape, bounds, params)
            for grid_shape in grid_shapes
        )

    print(
        f"jac_speedup: {len(configs)} configs x {reps} reps  seed={seed}",
        flush=True,
    )

    # Run each config sequentially (setup is heavy; shuffle is per-config not per-rep
    # because we re-use setup across reps within a config).
    shuffled = list(configs)
    rng.shuffle(shuffled)

    rows: list[dict] = []
    for theory_name, json_path, grid_shape, bounds, params in shuffled:
        row = _time_one_config(theory_name, json_path, grid_shape, bounds, params, reps)
        rows.append(row)

    return {"metadata": _metadata(reps, seed), "results": rows}


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=REPO_ROOT / "benchmark_results" / "canonical" / "jac_speedup.json",
    )
    parser.add_argument("--reps", type=int, default=N_REPS)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    data = run(reps=args.reps, seed=args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
