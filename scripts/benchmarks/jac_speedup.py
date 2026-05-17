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
        # 4 slots, 1-D → n_total = 4N.  Dense tier (n<2000): N≤500.
        # Densified at low N to resolve the dense-tier speedup region.
        "coupled_scalars",
        REPO_ROOT / "examples" / "data" / "coupled_scalars.json",
        [
            (32,),
            (48,),
            (64,),
            (96,),
            (128,),
            (192,),
            (256,),
            (384,),
            (500,),
            (768,),
            (1024,),
            (1536,),
            (2048,),
            (3072,),
            (4096,),
            (6144,),
            (8192,),
            (12288,),
            (16384,),
            (24576,),
            (32768,),
            (49152,),
            (65536,),
        ],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
    ),
    (
        # 12 slots, 1-D → n_total = 12N.  Dense tier (n<2000): N≤166.
        # Densified at low N; upper end trimmed (was N=32768 ≈ 393k).
        "gertsenshtein",
        REPO_ROOT / "examples" / "data" / "gertsenshtein.json",
        [
            (4,),
            (8,),
            (12,),
            (16,),
            (24,),
            (32,),
            (48,),
            (64,),
            (96,),
            (128,),
            (160,),
            (192,),
            (256,),
            (384,),
            (512,),
            (768,),
            (1024,),
            (1536,),
            (2048,),
            (3072,),
            (4096,),
            (6144,),
            (8192,),
            (12288,),
            (16384,),
        ],
        ((0.0, 100.0),),
        {"kappa": 1.0, "B0": 0.3},
    ),
    (
        # massive_gravity_3d has d2_t in constraint equations — incompatible
        # with the analytical Jacobian builder.  Use navier_cauchy_2d instead:
        # 2 displacement fields (+ 2 velocity slots) on an N×N grid give
        # n_total = 4N², the same O(N²) 2-D scaling the figure aims to show.
        # Densified at low N; upper end trimmed (was N=256 ≈ 262k).
        "navier_cauchy_2d",
        REPO_ROOT / "examples" / "data" / "navier_cauchy_2d.json",
        [
            (4, 4),
            (6, 6),
            (8, 8),
            (10, 10),
            (12, 12),
            (14, 14),
            (16, 16),
            (20, 20),
            (24, 24),
            (32, 32),
            (40, 40),
            (48, 48),
            (56, 56),
            (64, 64),
            (80, 80),
            (96, 96),
            (112, 112),
            (128, 128),
            (160, 160),
            (192, 192),
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


def _probe_batch_size(call, target_s: float = _PROBE_TARGET_S) -> int:
    """Pick a batch size for one callback so probe time ≥ target_s.

    Probing per-callback is essential: the dense path can be 1e5x slower than
    sparse, so using a sparse-derived batch_size for dense over-batches by 5
    orders of magnitude and stalls the job.
    """
    t0 = time.perf_counter()
    call()
    probe_s = time.perf_counter() - t0
    if probe_s >= target_s:
        return 1
    return max(1, math.ceil(target_s / probe_s))


def probe_jac_batch_size(
    setup: dict, target_s: float = _PROBE_TARGET_S
) -> dict[str, int]:
    """Per-callback adaptive batch sizes."""
    return {
        "dense": _probe_batch_size(setup["call_dense"], target_s)
        if setup["call_dense"]
        else 1,
        "sparse": _probe_batch_size(setup["call_sparse"], target_s),
        "gmres": _probe_batch_size(setup["call_gmres"], target_s),
        "matvec": _probe_batch_size(setup["call_matvec"], target_s),
    }


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


def time_one_jac_rep(setup: dict, batch_sizes: dict[str, int] | int = 1) -> dict:
    """Time one rep of all four paths. Returns per-call wall times in seconds.

    The FD proxy is timed as ``n_colors`` sequential mat-vec calls (simulating
    the n_colors residual evaluations that SUNDIALS performs per Jacobian fill).
    ``batch_sizes`` may be a single int (legacy) or a per-callback dict from
    ``probe_jac_batch_size`` — the latter is essential because dense and
    sparse callbacks differ by orders of magnitude in cost.
    """
    if isinstance(batch_sizes, int):
        bs = {
            "dense": batch_sizes,
            "sparse": batch_sizes,
            "gmres": batch_sizes,
            "matvec": batch_sizes,
        }
    else:
        bs = batch_sizes

    results: dict[str, float] = {}
    n_colors = setup["n_colors"]

    # Dense
    if setup["call_dense"] is not None:
        b = bs["dense"]
        t0 = time.perf_counter()
        for _ in range(b):
            setup["call_dense"]()
        results["dense_s"] = (time.perf_counter() - t0) / b
    else:
        results["dense_s"] = float("nan")

    # Sparse
    b = bs["sparse"]
    t0 = time.perf_counter()
    for _ in range(b):
        setup["call_sparse"]()
    results["sparse_s"] = (time.perf_counter() - t0) / b

    # GMRES jactimes (one JVP call)
    b = bs["gmres"]
    t0 = time.perf_counter()
    for _ in range(b):
        setup["call_gmres"]()
    results["gmres_s"] = (time.perf_counter() - t0) / b

    # FD proxy: time one mat-vec (batched), then scale by n_colors.  Scaling
    # is exact for linear FD (the colored evaluations are independent) and
    # avoids paying n_colors × batch_size matvecs in the inner loop.
    b = bs["matvec"]
    t0 = time.perf_counter()
    for _ in range(b):
        setup["call_matvec"]()
    matvec_s = (time.perf_counter() - t0) / b
    results["matvec_s"] = matvec_s
    results["fd_s"] = matvec_s * n_colors

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
    batch_sizes = probe_jac_batch_size(setup)

    dense_reps: list[float] = []
    sparse_reps: list[float] = []
    gmres_reps: list[float] = []
    fd_reps: list[float] = []

    for _ in range(reps):
        rep = time_one_jac_rep(setup, batch_sizes=batch_sizes)
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


def _config_key(
    theory_name: str, grid_shape: tuple[int, ...]
) -> tuple[str, tuple[int, ...]]:
    return (theory_name, tuple(grid_shape))


def _load_existing_rows(path: Path) -> dict[tuple[str, tuple[int, ...]], dict]:
    """Load (theory, grid_shape) -> row mapping from an existing JSON file."""
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    out: dict[tuple[str, tuple[int, ...]], dict] = {}
    for row in data.get("results", []):
        key = (row["theory"], tuple(row["grid_shape"]))
        out[key] = row
    return out


def run(
    reps: int = N_REPS,
    seed: int = 42,
    out_path: Path | None = None,
    *,
    skip_existing: bool = False,
    order: str = "shuffle",
) -> dict:
    """Run the full benchmark.

    ``order`` controls the traversal sequence: ``shuffle`` (default) randomises
    across all (theory, N_grid) configurations to average out monotonic
    thermal drift; ``ascending`` sorts by n_total so cheap small-N configs
    complete first — useful when you want rapid early checkpoints.  After
    each config completes the full JSON is rewritten to ``out_path`` so a
    job killed by walltime still leaves partial data on disk.  If
    ``skip_existing`` is set, configs already present in ``out_path`` are
    skipped — useful for splitting the sweep across multiple INTR rounds.
    """  # noqa: DOC501  (ValueError raised for invalid `order`)
    rng = random.Random(seed)  # noqa: S311

    # Build all configs
    configs = []
    for theory_name, json_path, grid_shapes, bounds, params in THEORIES:
        configs.extend(
            (theory_name, json_path, grid_shape, bounds, params)
            for grid_shape in grid_shapes
        )

    existing: dict[tuple[str, tuple[int, ...]], dict] = {}
    if out_path is not None and skip_existing:
        existing = _load_existing_rows(out_path)
        if existing:
            print(
                f"jac_speedup: found {len(existing)} existing rows in {out_path}",
                flush=True,
            )

    print(
        f"jac_speedup: {len(configs)} configs x {reps} reps  seed={seed}  order={order}",
        flush=True,
    )

    # Run each config sequentially (setup is heavy; order is per-config not
    # per-rep because we re-use setup across reps within a config).
    ordered = list(configs)
    if order == "shuffle":
        rng.shuffle(ordered)
    elif order == "ascending":
        # Sort by n_total (= slots × product(grid_shape)).  Slot count requires
        # loading the spec; cheap to do once.  Use product(grid_shape) as a
        # proxy (correct ordering within a theory; close-enough across).
        from math import prod

        ordered.sort(key=lambda c: prod(c[2]))
    else:
        msg = f"unknown order={order!r}; use 'shuffle' or 'ascending'"
        raise ValueError(msg)
    shuffled = ordered

    rows: list[dict] = list(existing.values())
    for theory_name, json_path, grid_shape, bounds, params in shuffled:
        if skip_existing and _config_key(theory_name, grid_shape) in existing:
            print(
                f"  SKIP {theory_name} shape={grid_shape} (already in JSON)", flush=True
            )
            continue
        row = _time_one_config(theory_name, json_path, grid_shape, bounds, params, reps)
        rows.append(row)

        # Per-config checkpoint: rewrite JSON immediately so a job killed by
        # the INTR walltime still leaves partial results on disk.
        if out_path is not None:
            data = {"metadata": _metadata(reps, seed), "results": rows}
            out_path.parent.mkdir(parents=True, exist_ok=True)
            tmp = out_path.with_suffix(out_path.suffix + ".tmp")
            with tmp.open("w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2)
            tmp.replace(out_path)

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
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip configs already present in --out (resume mode).",
    )
    parser.add_argument(
        "--order",
        choices=["shuffle", "ascending"],
        default="shuffle",
        help="Config traversal order: 'shuffle' (default, randomises) or "
        "'ascending' (sort by n_total — small/cheap configs first).",
    )
    args = parser.parse_args()

    data = run(
        reps=args.reps,
        seed=args.seed,
        out_path=args.out,
        skip_existing=args.skip_existing,
        order=args.order,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
