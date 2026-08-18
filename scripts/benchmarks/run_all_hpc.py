r"""HPC master driver — runs all three App-C benchmarks in one shuffled,
isolated, parallel sweep.

Designed for the CSD3 sapphire interactive INTR session (112 cores, 1 h).
Each (benchmark × config × rep) tuple is submitted as an independent task
to a `ProcessPoolExecutor` whose workers run with single-threaded BLAS and
are pinned to private CPU cores.  The task list is shuffled with a fixed
seed so time-correlated noise is decorrelated from any single
(benchmark, theory, N) axis.

CRITICAL: this module sets `OMP_NUM_THREADS=1` etc. **before** importing
numpy. Do not insert any numpy / scipy / tidal imports above the env-var
setup block at the top of `__main__`.

Usage on the compute node, after `hpc_shuttle.sh attach <jobid>`:

    cd $REMOTE_ROOT
    nohup uv run python scripts/benchmarks/run_all_hpc.py \\
        --tier baseline --reps 30 --parallel 112 \\
        > run.log 2>&1 &
    # tail -f run.log to monitor
    # (or use tmux session for resilience across SSH drops)

Once the canonical JSONs land in `benchmark_results/canonical/`,
they ride the `hpc_shuttle.sh pull <jobid>` whitelist back to local.
"""

from __future__ import annotations

import argparse
import datetime
import json
import operator
import os
import platform
import random
import socket
import subprocess  # noqa: S404
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))


# ---------------------------------------------------------------------------
# Tier definitions: which (theory, N) configs to time at each tier.
# ---------------------------------------------------------------------------

# Padé 1D theories share the same N ladder.
_PADE_1D_BASELINE = [
    128,
    192,
    256,
    384,
    512,
    768,
    1024,
    1536,
    2048,
    3072,
    4096,
    6144,
    8192,
]
_PADE_1D_EXTENDED = [*_PADE_1D_BASELINE, 12288, 16384]
_PADE_1D_AMBITIOUS = sorted(
    {*_PADE_1D_EXTENDED, 2304, 2816, 3584, 4608, 5120, 5632, 7168}
)

_PADE_2D_BASELINE = [16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96]
_PADE_2D_EXTENDED = [*_PADE_2D_BASELINE, 112, 128]
_PADE_2D_AMBITIOUS = sorted({*_PADE_2D_EXTENDED, 36, 44, 52, 60, 72, 88})

_SPARSE_BASELINE = [64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048]
_SPARSE_SPARSE_BASELINE = [*_SPARSE_BASELINE, 3072, 4096]
_SPARSE_EXTENDED_DENSE = [*_SPARSE_BASELINE, 3072]
_SPARSE_AMBITIOUS_FILL = sorted(
    {*_SPARSE_BASELINE, 192, 320, 448, 576, 704, 896, 1152, 1280, 1408, 1664, 1792}
)

_NYQUIST_BASELINE = _PADE_1D_BASELINE
_NYQUIST_EXTENDED = _PADE_1D_EXTENDED
_NYQUIST_AMBITIOUS = _PADE_1D_AMBITIOUS


def _tier_pade_1d(tier: str) -> list[int]:
    return {
        "baseline": _PADE_1D_BASELINE,
        "extended": _PADE_1D_EXTENDED,
        "ambitious": _PADE_1D_AMBITIOUS,
    }[tier]


def _tier_pade_2d(tier: str) -> list[int]:
    return {
        "baseline": _PADE_2D_BASELINE,
        "extended": _PADE_2D_EXTENDED,
        "ambitious": _PADE_2D_AMBITIOUS,
    }[tier]


def _tier_sparse_dense(tier: str) -> list[int]:
    return {
        "baseline": _SPARSE_BASELINE,
        "extended": _SPARSE_EXTENDED_DENSE,
        "ambitious": sorted(set(_SPARSE_EXTENDED_DENSE + _SPARSE_AMBITIOUS_FILL)),
    }[tier]


def _tier_sparse_sparse(tier: str) -> list[int]:
    return {
        "baseline": _SPARSE_SPARSE_BASELINE,
        "extended": _SPARSE_SPARSE_BASELINE,
        "ambitious": sorted(set(_SPARSE_SPARSE_BASELINE + _SPARSE_AMBITIOUS_FILL)),
    }[tier]


def _tier_nyquist(tier: str) -> list[int]:
    return {
        "baseline": _NYQUIST_BASELINE,
        "extended": _NYQUIST_EXTENDED,
        "ambitious": _NYQUIST_AMBITIOUS,
    }[tier]


# ---------------------------------------------------------------------------
# Per-config concurrency caps (memory safety).
# Key: (benchmark, n, path_or_None) → max workers running this config in parallel.
# ---------------------------------------------------------------------------

# Both dense and sparse paths build the dense convolution matrix as an
# intermediate (~1.2 GB at N=2048, ~4.8 GB at N=4096, ~10.8 GB at N=4096
# during the dense-times-DT * sparse conversion peak). Cap BOTH paths
# identically so total peak memory stays within ~150 GB across all
# concurrently-running workers on a sapphire node (≈256 GB total RAM,
# leaving headroom for Padé/Nyquist setups + OS).
_CONFIG_CAPS = {
    ("sparse_csc_vs_dense", 4096, "dense"): 12,
    ("sparse_csc_vs_dense", 4096, "sparse"): 12,
    ("sparse_csc_vs_dense", 3072, "dense"): 20,
    ("sparse_csc_vs_dense", 3072, "sparse"): 20,
    ("sparse_csc_vs_dense", 2048, "dense"): 32,
    ("sparse_csc_vs_dense", 2048, "sparse"): 32,
    ("sparse_csc_vs_dense", 1536, "dense"): 48,
    ("sparse_csc_vs_dense", 1536, "sparse"): 48,
    ("sparse_csc_vs_dense", 1024, "dense"): 80,
    ("sparse_csc_vs_dense", 1024, "sparse"): 80,
}


def _config_cap(benchmark: str, n: int, path: str | None) -> int:
    return _CONFIG_CAPS.get((benchmark, n, path), 999)


# ---------------------------------------------------------------------------
# Worker initializer: pin to a core; cache setups per (benchmark, key).
# ---------------------------------------------------------------------------

_WORKER_CORE: int | None = None
# OrderedDict-backed LRU; capped to limit per-worker memory growth from
# accumulated sparse-matrix builds.
from collections import OrderedDict

_SETUP_CACHE: OrderedDict[tuple, Any] = OrderedDict()
_SETUP_CACHE_MAX = 1  # one cached setup per worker — matrices can be huge,
# and with N_workers × LRU × peak_setup_size as the memory budget, the
# previous LRU=2 setting pushed sapphire over the OOM line on shuffled
# task lists that span the largest 2D massive_gravity_3d configs.


def _cache_get(key):
    """Return cached setup if present, moving it to the MRU end."""
    if key in _SETUP_CACHE:
        _SETUP_CACHE.move_to_end(key)
        return _SETUP_CACHE[key]
    return None


def _cache_put(key, value) -> None:
    """Insert; evict LRU when over capacity."""
    _SETUP_CACHE[key] = value
    _SETUP_CACHE.move_to_end(key)
    while len(_SETUP_CACHE) > _SETUP_CACHE_MAX:
        _SETUP_CACHE.popitem(last=False)


def _worker_init(core_id: int) -> None:
    global _WORKER_CORE
    _WORKER_CORE = core_id
    try:
        os.sched_setaffinity(0, {core_id})
    except (AttributeError, PermissionError, OSError):
        pass  # macOS dev, or restricted env — degrade gracefully


def _setup_pade(theory_name, json_path_str, grid_shape, bounds, params):
    from scripts.benchmarks import pade_vs_eig as bm

    key = ("pade", theory_name, tuple(grid_shape))
    cached = _cache_get(key)
    if cached is not None:
        return cached
    setup = bm.build_pade_setup(
        theory_name,
        Path(json_path_str),
        tuple(grid_shape),
        tuple(tuple(b) for b in bounds),
        params,
    )
    bm.warmup_pade(setup)
    setup["batch_size"] = bm.probe_pade_batch_size(setup)
    _cache_put(key, setup)
    return setup


def _setup_sparse(n, path):
    from scripts.benchmarks import sparse_csc_vs_dense as bm

    key = ("sparse", n, path)
    cached = _cache_get(key)
    if cached is not None:
        return cached
    setup = bm.build_sparse_setup(n, path)
    bm.warmup_sparse(setup)
    setup["batch_size"] = bm.probe_sparse_batch_size(setup)
    _cache_put(key, setup)
    return setup


def _do_pade_rep(task: dict) -> dict:
    from scripts.benchmarks import pade_vs_eig as bm

    try:
        setup = _setup_pade(
            task["theory"],
            task["json_path"],
            task["grid_shape"],
            task["bounds"],
            task["params"],
        )
        rep = bm.time_one_pade_rep(setup, batch_size=setup["batch_size"])
        return {
            "benchmark": "pade_vs_eig",
            "theory": task["theory"],
            "n_label": task["n_label"],
            "n_modes": setup["n_modes"],
            "n_slots": setup["n_slots"],
            "batch_size": setup["batch_size"],
            "rep_idx": task["rep_idx"],
            "pade_s_per_mode": rep["pade_s_per_mode"],
            "eig_s_per_mode": rep["eig_s_per_mode"],
            "worker_core": _WORKER_CORE,
        }
    except Exception as exc:
        return {
            "benchmark": "pade_vs_eig",
            "theory": task["theory"],
            "n_label": task["n_label"],
            "rep_idx": task["rep_idx"],
            "error": f"{type(exc).__name__}: {exc}",
        }


def _do_sparse_rep(task: dict) -> dict:
    from scripts.benchmarks import sparse_csc_vs_dense as bm

    try:
        setup = _setup_sparse(task["n"], task["path"])
        rep = bm.time_one_sparse_rep(setup, batch_size=setup["batch_size"])
        return {
            "benchmark": "sparse_csc_vs_dense",
            "n": task["n"],
            "path": task["path"],
            "matrix_size": setup["matrix_size"],
            "density_pct": setup["density_pct"],
            "note": setup["note"],
            "batch_size": setup["batch_size"],
            "rep_idx": task["rep_idx"],
            "wall_s": rep["wall_s"],
            "worker_core": _WORKER_CORE,
        }
    except Exception as exc:
        return {
            "benchmark": "sparse_csc_vs_dense",
            "n": task["n"],
            "path": task["path"],
            "rep_idx": task["rep_idx"],
            "error": f"{type(exc).__name__}: {exc}",
        }


def _do_nyquist(task: dict) -> dict:
    from scripts.benchmarks import nyquist_energy as bm

    try:
        return {
            "benchmark": "nyquist_energy",
            **bm.run_one_nyquist(
                task["theory"],
                Path(task["json_path"]),
                task["n"],
                task["params"],
                zero_nyquist=task["zero_nyquist"],
            ),
            "worker_core": _WORKER_CORE,
        }
    except Exception as exc:
        return {
            "benchmark": "nyquist_energy",
            "theory": task["theory"],
            "n": task["n"],
            "zero_nyquist": task["zero_nyquist"],
            "error": f"{type(exc).__name__}: {exc}",
        }


_DISPATCH = {
    "pade_vs_eig": _do_pade_rep,
    "sparse_csc_vs_dense": _do_sparse_rep,
    "nyquist_energy": _do_nyquist,
}


def _execute(task: dict) -> dict:
    return _DISPATCH[task["benchmark"]](task)


# ---------------------------------------------------------------------------
# Task list construction.
# ---------------------------------------------------------------------------


def _build_task_list(tier: str, n_reps: int) -> list[dict]:
    """Union of all (benchmark × config × rep) tuples for the given tier."""
    tasks: list[dict] = []

    # Padé 1D theories.
    pade_1d_theories = [
        (
            "coupled_scalars",
            str(REPO_ROOT / "examples" / "data" / "coupled_scalars.json"),
            ((0.0, 100.0),),
            {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
        ),
        (
            "gertsenshtein",
            str(REPO_ROOT / "examples" / "data" / "gertsenshtein.json"),
            ((0.0, 100.0),),
            {"kappa": 1.0, "B0": 0.3},
        ),
    ]
    for theory, json_path, bounds, params in pade_1d_theories:
        for n in _tier_pade_1d(tier):
            tasks.extend(
                {
                    "benchmark": "pade_vs_eig",
                    "theory": theory,
                    "json_path": json_path,
                    "grid_shape": (n,),
                    "bounds": bounds,
                    "params": params,
                    "n_label": str(n),
                    "rep_idx": rep,
                }
                for rep in range(n_reps)
            )

    # Padé 2D theory.
    for n in _tier_pade_2d(tier):
        tasks.extend(
            {
                "benchmark": "pade_vs_eig",
                "theory": "massive_gravity_3d",
                "json_path": str(
                    REPO_ROOT / "examples" / "data" / "massive_gravity_3d.json"
                ),
                "grid_shape": (n, n),
                "bounds": ((0.0, 50.0), (0.0, 50.0)),
                "params": {"m2": 1.0},
                "n_label": f"{n}^2",
                "rep_idx": rep,
            }
            for rep in range(n_reps)
        )

    # Sparse-CSC dense.
    for n in _tier_sparse_dense(tier):
        tasks.extend(
            {
                "benchmark": "sparse_csc_vs_dense",
                "n": n,
                "path": "dense",
                "rep_idx": rep,
            }
            for rep in range(n_reps)
        )

    # Sparse-CSC sparse.
    for n in _tier_sparse_sparse(tier):
        tasks.extend(
            {
                "benchmark": "sparse_csc_vs_dense",
                "n": n,
                "path": "sparse",
                "rep_idx": rep,
            }
            for rep in range(n_reps)
        )

    # Nyquist — one task per (theory, n, projection); no reps.
    nyquist_theories = [
        (
            "coupled_scalars",
            str(REPO_ROOT / "examples" / "data" / "coupled_scalars.json"),
            {"kappa": 1.0, "B0": 0.1, "omegaP2": 0.0, "mg2": 0.0},
        ),
        (
            "coupled_scalars_massive",
            str(REPO_ROOT / "examples" / "data" / "coupled_scalars.json"),
            {"kappa": 1.0, "B0": 0.1, "omegaP2": 1.0, "mg2": 0.0},
        ),
    ]
    for theory, json_path, params in nyquist_theories:
        for n in _tier_nyquist(tier):
            tasks.extend(
                {
                    "benchmark": "nyquist_energy",
                    "theory": theory,
                    "json_path": json_path,
                    "n": n,
                    "params": params,
                    "zero_nyquist": zero,
                    "rep_idx": 0,
                }
                for zero in (True, False)
            )

    return tasks


# ---------------------------------------------------------------------------
# Per-config concurrency throttling at submission time.
#
# We can't easily put a multiprocessing.Semaphore inside the worker (no shared
# state across fork in some kernels), so we throttle in the SUBMITTER: only
# advance to the next task when the active task budget for its (bench,n,path)
# allows it. Cheap configs default to cap=999 (no throttle).
# ---------------------------------------------------------------------------


def _config_key(task: dict) -> tuple:
    return (
        task["benchmark"],
        task.get("n", task.get("n_label")),
        task.get("path"),
    )


def _config_cap_for(task: dict) -> int:
    n = task.get("n")
    if n is None:
        try:
            n = int(task.get("n_label", "0").split("^")[0])
        except ValueError:
            n = 0
    return _config_cap(task["benchmark"], n, task.get("path"))


# ---------------------------------------------------------------------------
# Aggregation: collapse the flat rep list into the canonical JSON schemas.
# ---------------------------------------------------------------------------


def _load_existing_pade_reps(path: Path) -> dict[tuple[str, str], dict]:
    """Read an existing pade_vs_eig.json and return per-config seed reps."""
    if not path.exists():
        return {}
    try:
        d = json.load(path.open(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    seeds: dict[tuple[str, str], dict] = {}
    for r in d.get("results", []):
        key = (r["theory"], r["n"])
        seeds[key] = {
            "pade_reps": list(r.get("pade_s_reps", [])),
            "eig_reps": list(r.get("eig_s_reps", [])),
            "n_modes": r.get("n_modes"),
            "n_slots": r.get("n_slots"),
        }
    return seeds


def _load_existing_sparse_reps(path: Path) -> dict[tuple[int, str], dict]:
    if not path.exists():
        return {}
    try:
        d = json.load(path.open(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    seeds: dict[tuple[int, str], dict] = {}
    for r in d.get("results", []):
        key = (int(r["n"]), r["path"])
        seeds[key] = {
            "wall_reps": list(r.get("wall_s_reps", [])),
            "matrix_size": r.get("matrix_size"),
            "density_pct": r.get("density_pct"),
            "note": r.get("note", "below_density_cutoff"),
        }
    return seeds


def _aggregate_pade(rep_rows: list[dict], seeds: dict | None = None) -> list[dict]:
    """Group Padé reps by (theory, n_label) → mean/std rows.

    If ``seeds`` is provided (from an existing canonical JSON), its per-rep
    arrays are PREPENDED so new reps APPEND rather than overwrite previous
    runs' data — supports cross-host accumulation (local + HPC + …).
    """
    import numpy as np  # imported here to keep top-of-module BLAS-clean

    seeds = seeds or {}
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in rep_rows:
        if "error" in r:
            continue
        grouped[r["theory"], r["n_label"]].append(r)

    all_keys = set(grouped) | set(seeds)
    rows: list[dict] = []
    for theory, n_label in all_keys:
        reps = grouped.get((theory, n_label), [])
        seed = seeds.get((theory, n_label), {})
        pade_reps = list(seed.get("pade_reps", [])) + [
            r["pade_s_per_mode"] for r in reps
        ]
        eig_reps = list(seed.get("eig_reps", [])) + [r["eig_s_per_mode"] for r in reps]
        if not pade_reps:
            continue
        pade_mean = float(np.mean(pade_reps))
        pade_std = float(np.std(pade_reps, ddof=1)) if len(pade_reps) > 1 else 0.0
        eig_mean = float(np.mean(eig_reps))
        eig_std = float(np.std(eig_reps, ddof=1)) if len(eig_reps) > 1 else 0.0
        n_modes = reps[0]["n_modes"] if reps else seed.get("n_modes")
        n_slots = reps[0]["n_slots"] if reps else seed.get("n_slots")
        rows.append(
            {
                "theory": theory,
                "n": n_label,
                "n_modes": n_modes,
                "n_slots": n_slots,
                "pade_s": pade_mean,
                "eig_s": eig_mean,
                "ratio": pade_mean / eig_mean if eig_mean > 0 else float("nan"),
                "pade_s_mean": pade_mean,
                "pade_s_std": pade_std,
                "pade_s_reps": pade_reps,
                "eig_s_mean": eig_mean,
                "eig_s_std": eig_std,
                "eig_s_reps": eig_reps,
            }
        )
    return sorted(rows, key=lambda r: (r["theory"], _n_to_float(r["n"])))


def _aggregate_sparse(rep_rows: list[dict], seeds: dict | None = None) -> list[dict]:
    import numpy as np

    seeds = seeds or {}
    grouped: dict[tuple[int, str], list[dict]] = defaultdict(list)
    for r in rep_rows:
        if "error" in r:
            continue
        grouped[r["n"], r["path"]].append(r)
    all_keys = set(grouped) | set(seeds)
    rows: list[dict] = []
    for n, path in all_keys:
        reps = grouped.get((n, path), [])
        seed = seeds.get((n, path), {})
        wall_reps = list(seed.get("wall_reps", [])) + [r["wall_s"] for r in reps]
        if not wall_reps:
            continue
        mean = float(np.mean(wall_reps))
        std = float(np.std(wall_reps, ddof=1)) if len(wall_reps) > 1 else 0.0
        matrix_size = reps[0]["matrix_size"] if reps else seed.get("matrix_size")
        density = reps[0]["density_pct"] if reps else seed.get("density_pct")
        note = reps[0]["note"] if reps else seed.get("note", "below_density_cutoff")
        rows.append(
            {
                "n": n,
                "path": path,
                "matrix_size": matrix_size,
                "density_pct": density,
                "wall_s": mean,
                "wall_s_mean": mean,
                "wall_s_std": std,
                "wall_s_reps": wall_reps,
                "n_reps": len(wall_reps),
                "note": note,
            }
        )
    return sorted(rows, key=operator.itemgetter("path", "n"))


def _load_existing_nyquist_rows(path: Path) -> dict[tuple[str, int], dict]:
    """Read an existing nyquist_energy.json and return per-config rows by key."""
    if not path.exists():
        return {}
    try:
        d = json.load(path.open(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return {}
    return {(r["theory"], int(r["n"])): r for r in d.get("results", [])}


def _aggregate_nyquist(rep_rows: list[dict], seeds: dict | None = None) -> list[dict]:
    """Collapse pairs of (zero=T, zero=F) into one row per (theory, n).

    Merges with ``seeds`` (existing canonical-JSON rows) so a partial
    re-run only overwrites the (theory, n) keys it actually measured —
    older entries for other N values survive untouched. The Nyquist
    measurement is deterministic so there are no per-rep arrays to
    accumulate; merging is "newer-row wins" only for keys re-measured.
    """
    seeds = seeds or {}
    grouped: dict[tuple[str, int], dict] = defaultdict(dict)
    for r in rep_rows:
        if "error" in r:
            continue
        grouped[r["theory"], int(r["n"])][r["zero_nyquist"]] = r["abs_dE_over_E"]

    # Start with existing seeds; overwrite per-key only where the new
    # run measured both projections.
    merged: dict[tuple[str, int], dict] = dict(seeds)
    for (theory, n), d in grouped.items():
        row = {
            "theory": theory,
            "n": n,
            "abs_dE_over_E_fixed": d.get(True, float("nan")),
            "abs_dE_over_E_nofixed": d.get(False, float("nan")),
        }
        # Don't clobber a fully-populated seed with a half-populated new row.
        if (theory, n) in seeds and (
            not isinstance(row["abs_dE_over_E_fixed"], float)
            or row["abs_dE_over_E_fixed"] != row["abs_dE_over_E_fixed"]  # NaN
            or row["abs_dE_over_E_nofixed"] != row["abs_dE_over_E_nofixed"]
        ):
            continue
        merged[theory, n] = row

    return sorted(merged.values(), key=operator.itemgetter("theory", "n"))


def _n_to_float(label: Any) -> float:
    """For sort ordering: '256' → 256.0, '32^2' → 1024.0, 256 → 256.0."""
    if isinstance(label, (int, float)):
        return float(label)
    s = str(label)
    if "^" in s:
        base, exp = s.split("^")
        return float(base) ** float(exp)
    return float(s)


# ---------------------------------------------------------------------------
# Metadata.
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


def _metadata(tier: str, n_reps: int, seed: int, n_workers: int) -> dict:
    import numpy as np
    import scipy

    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "driver": "run_all_hpc.py",
        "parameters": {
            "tier": tier,
            "n_reps": n_reps,
            "shuffle_seed": seed,
            "n_workers": n_workers,
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
            "mkl_num_threads": os.environ.get("MKL_NUM_THREADS"),
            "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        },
    }


# ---------------------------------------------------------------------------
# Main driver.
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tier",
        choices=("baseline", "extended", "ambitious"),
        default="baseline",
    )
    parser.add_argument(
        "--reps", type=int, default=30, help="N_reps per Padé/Sparse config"
    )
    parser.add_argument("--parallel", type=int, default=112, help="Worker pool size")
    parser.add_argument("--seed", type=int, default=42, help="Shuffle seed")
    parser.add_argument(
        "--repeat-passes",
        type=int,
        default=1,
        help="Run the whole tier this many times, concatenating reps",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=REPO_ROOT / "benchmark_results" / "canonical",
    )
    parser.add_argument(
        "--no-merge",
        dest="merge",
        action="store_false",
        default=True,
        help="Overwrite canonical JSONs instead of appending new reps onto existing per-rep arrays",
    )
    parser.add_argument(
        "--max-n",
        type=int,
        default=None,
        help="Skip any (theory, N) config with N > max_n. For local dev-container "
        "runs where memory or wall-time forbids the heavy end of the sweep.",
    )
    parser.add_argument(
        "--min-n",
        type=int,
        default=None,
        help="Skip any (theory, N) config with N < min_n. Pair with --merge to "
        "resume a timed-out run by collecting only the missing high-N reps.",
    )
    args = parser.parse_args()

    print(
        f"run_all_hpc: tier={args.tier} reps={args.reps} "
        f"parallel={args.parallel} passes={args.repeat_passes} "
        f"seed={args.seed}",
        flush=True,
    )

    # Load merge seeds ONCE before any tasks run.  Passing them explicitly to
    # every _write_canonical call prevents the periodic-checkpoint path from
    # re-reading the file it just wrote as fresh seeds and accumulating reps
    # exponentially across checkpoint writes.
    initial_pade_seeds = (
        _load_existing_pade_reps(args.out_dir / "pade_vs_eig.json")
        if args.merge
        else {}
    )
    initial_sparse_seeds = (
        _load_existing_sparse_reps(args.out_dir / "sparse_csc_vs_dense.json")
        if args.merge
        else {}
    )
    initial_nyquist_seeds = (
        _load_existing_nyquist_rows(args.out_dir / "nyquist_energy.json")
        if args.merge
        else {}
    )

    # Accumulate per-rep rows across passes (different seeds per pass).
    all_rep_rows: list[dict] = []
    pass_seeds: list[int] = []

    for pass_idx in range(args.repeat_passes):
        pass_seed = args.seed + pass_idx
        pass_seeds.append(pass_seed)
        tasks = _build_task_list(args.tier, args.reps)
        if args.max_n is not None or args.min_n is not None:

            def _task_n(t):
                n = t.get("n")
                if n is None:
                    try:
                        n = int(t.get("n_label", "0").split("^")[0])
                    except ValueError:
                        n = 0
                return int(n)

            if args.max_n is not None:
                tasks = [t for t in tasks if _task_n(t) <= args.max_n]
            if args.min_n is not None:
                tasks = [t for t in tasks if _task_n(t) >= args.min_n]
        rng = random.Random(pass_seed)
        rng.shuffle(tasks)
        print(
            f"  pass {pass_idx + 1}/{args.repeat_passes}: "
            f"{len(tasks)} tasks, seed={pass_seed}"
            f"{', max_n=' + str(args.max_n) if args.max_n else ''}",
            flush=True,
        )

        # Throttled submission: keep at most _config_cap_for(t) of each
        # (bench, n, path) in flight at once.
        active_per_key: dict[tuple, int] = defaultdict(int)
        in_flight: set = set()
        pending = list(tasks)
        pending.reverse()  # pop() yields the shuffled head

        t_start = time.time()
        # Build worker→core mapping; only use the first `parallel` cores.
        avail_cores = (
            sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else list(range(os.cpu_count() or 1))
        )
        worker_cores = (
            avail_cores * ((args.parallel // max(1, len(avail_cores))) + 1)
        )[: args.parallel]

        # Use spawn to ensure a clean BLAS-thread state in each child
        # (the parent has already set the env vars in __main__).
        import multiprocessing as mp

        ctx = mp.get_context("fork")  # fork inherits parent env vars

        # Per-worker SETUP_CACHE growth is bounded inside the worker via an
        # LRU cap (see `_setup_pade`/`_setup_sparse` callers) — workers
        # therefore hold at most a few matrix-build artifacts at a time,
        # preventing the multi-GB cache bloat that killed the first pass.
        with ProcessPoolExecutor(
            max_workers=args.parallel,
            mp_context=ctx,
            initializer=_worker_init,
            initargs=(worker_cores[0],),
        ) as pool:
            futures = {}

            def _try_submit() -> bool:
                """Submit as many pending tasks as caps allow. Returns True if anything submitted."""
                submitted = False
                i = len(pending) - 1
                while i >= 0 and len(in_flight) < args.parallel:
                    t = pending[i]
                    key = _config_key(t)
                    if active_per_key[key] >= _config_cap_for(t):
                        i -= 1
                        continue
                    pending.pop(i)
                    fut = pool.submit(_execute, t)
                    futures[fut] = (t, key)
                    in_flight.add(fut)
                    active_per_key[key] += 1
                    submitted = True
                    i -= 1
                return submitted

            _try_submit()
            n_completed = 0
            last_checkpoint = time.time()
            while in_flight:
                done = next(as_completed(in_flight))
                in_flight.remove(done)
                t, key = futures.pop(done)
                active_per_key[key] -= 1
                try:
                    row = done.result()
                    all_rep_rows.append(row)
                except Exception as exc:
                    all_rep_rows.append(
                        {
                            "benchmark": t["benchmark"],
                            "error": f"future-exception: {type(exc).__name__}: {exc}",
                            **{k: v for k, v in t.items() if k != "params"},
                        }
                    )
                n_completed += 1
                _try_submit()

                # Periodic progress + checkpoint.
                now = time.time()
                if now - last_checkpoint > 30.0:
                    elapsed = now - t_start
                    print(
                        f"    [{elapsed:6.1f}s] completed {n_completed} / "
                        f"{len(tasks)} pass-{pass_idx + 1} tasks "
                        f"({len(pending)} pending, {len(in_flight)} in flight)",
                        flush=True,
                    )
                    _write_partial(
                        args.out_dir,
                        all_rep_rows,
                        args.tier,
                        args.reps,
                        args.seed,
                        args.parallel,
                        pass_seeds,
                    )
                    # Also write canonical JSONs so a Slurm timeout only loses
                    # the last 30 s of work rather than the entire run.
                    # Seeds are passed explicitly so this call does NOT re-read
                    # from disk (which would accumulate reps exponentially).
                    _write_canonical(
                        args.out_dir,
                        all_rep_rows,
                        args.tier,
                        args.reps,
                        args.seed,
                        args.parallel,
                        pass_seeds,
                        merge=args.merge,
                        pade_seeds=initial_pade_seeds,
                        sparse_seeds=initial_sparse_seeds,
                        nyquist_seeds=initial_nyquist_seeds,
                    )
                    last_checkpoint = now

        print(
            f"  pass {pass_idx + 1} complete: {time.time() - t_start:.1f}s wall",
            flush=True,
        )

    # Final write — same pre-loaded seeds as checkpoints.
    _write_canonical(
        args.out_dir,
        all_rep_rows,
        args.tier,
        args.reps,
        args.seed,
        args.parallel,
        pass_seeds,
        merge=args.merge,
        pade_seeds=initial_pade_seeds,
        sparse_seeds=initial_sparse_seeds,
        nyquist_seeds=initial_nyquist_seeds,
    )
    print(
        f"run_all_hpc: wrote canonical JSONs to {args.out_dir}",
        flush=True,
    )
    n_errors = sum(1 for r in all_rep_rows if "error" in r)
    if n_errors:
        print(
            f"  WARNING: {n_errors} / {len(all_rep_rows)} tasks errored "
            "(see `error` field in JSON metadata)",
            flush=True,
        )
    return 0


def _write_partial(
    out_dir: Path,
    rows: list[dict],
    tier: str,
    reps: int,
    seed: int,
    n_workers: int,
    pass_seeds: list[int],
) -> None:
    """Crash-resistant checkpoint of the in-progress run."""
    out_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "metadata": {
            **_metadata(tier, reps, seed, n_workers),
            "pass_seeds_completed_or_in_progress": pass_seeds,
            "status": "partial",
        },
        "rep_rows": rows,
    }
    path = out_dir / "run_all_hpc.partial.json"
    with path.open("w") as fh:
        json.dump(payload, fh, indent=1)


def _write_canonical(
    out_dir: Path,
    rows: list[dict],
    tier: str,
    reps: int,
    seed: int,
    n_workers: int,
    pass_seeds: list[int],
    *,
    merge: bool = True,
    pade_seeds: dict | None = None,
    sparse_seeds: dict | None = None,
    nyquist_seeds: dict | None = None,
) -> None:
    """Aggregate per-rep rows into the three canonical per-benchmark JSONs.

    When ``merge=True`` (default), the existing canonical JSONs' per-rep
    arrays are loaded as SEEDS that the new reps APPEND onto. This makes
    cross-host data collection additive: a local pass + an HPC pass + a
    later re-run all accumulate into one larger sample, rather than
    overwriting each other.

    Pass pre-loaded ``pade_seeds``/``sparse_seeds``/``nyquist_seeds`` to
    avoid re-reading from disk on every periodic checkpoint call. Without
    this, each checkpoint would re-read the file it just wrote as fresh
    seeds, accumulating reps exponentially across checkpoint writes.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    errors = [r for r in rows if "error" in r]

    pade_rows = [r for r in rows if r["benchmark"] == "pade_vs_eig"]
    sparse_rows = [r for r in rows if r["benchmark"] == "sparse_csc_vs_dense"]
    nyquist_rows = [r for r in rows if r["benchmark"] == "nyquist_energy"]

    if pade_seeds is None:
        pade_seeds = (
            _load_existing_pade_reps(out_dir / "pade_vs_eig.json") if merge else {}
        )
    if sparse_seeds is None:
        sparse_seeds = (
            _load_existing_sparse_reps(out_dir / "sparse_csc_vs_dense.json")
            if merge
            else {}
        )
    if nyquist_seeds is None:
        nyquist_seeds = (
            _load_existing_nyquist_rows(out_dir / "nyquist_energy.json")
            if merge
            else {}
        )

    meta_base = _metadata(tier, reps, seed, n_workers)
    meta_base["parameters"]["pass_seeds"] = pass_seeds
    meta_base["parameters"]["merge_mode"] = merge
    if errors:
        meta_base["parameters"]["n_errors"] = len(errors)
        meta_base["parameters"]["errors_sample"] = errors[:5]

    payload_pade = {
        "metadata": meta_base,
        "results": _aggregate_pade(pade_rows, seeds=pade_seeds),
    }
    with (out_dir / "pade_vs_eig.json").open("w") as fh:
        json.dump(payload_pade, fh, indent=2)

    payload_sparse = {
        "metadata": meta_base,
        "results": _aggregate_sparse(sparse_rows, seeds=sparse_seeds),
    }
    with (out_dir / "sparse_csc_vs_dense.json").open("w") as fh:
        json.dump(payload_sparse, fh, indent=2)

    payload_nyq = {
        "metadata": meta_base,
        "results": _aggregate_nyquist(nyquist_rows, seeds=nyquist_seeds),
    }
    with (out_dir / "nyquist_energy.json").open("w") as fh:
        json.dump(payload_nyq, fh, indent=2)

    # Drop the partial file once the canonical write succeeds.
    partial = out_dir / "run_all_hpc.partial.json"
    if partial.exists():
        partial.unlink()


if __name__ == "__main__":
    # Set thread limits FIRST, before any numpy/scipy/tidal import path is
    # exercised. The benchmark modules are only imported inside worker tasks
    # via the late `from scripts.benchmarks import ...` statements above.
    for var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        os.environ.setdefault(var, "1")
    sys.exit(main())
