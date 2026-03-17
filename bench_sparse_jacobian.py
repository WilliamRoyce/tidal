"""Benchmark: sparse analytical Jacobian vs FD colored Jacobian.

Compares wall time for IDA/CVODE simulations in the sparse tier
(2,000 < N ≤ 200,000) with and without analytical sparse jacfn.

Usage:
    uv run python bench_sparse_jacobian.py
"""

from __future__ import annotations

import json
import statistics
import time
from pathlib import Path
from typing import Any
from unittest.mock import patch

import numpy as np

# ---------------------------------------------------------------------------
# Benchmark configuration
# ---------------------------------------------------------------------------

EXAMPLES_DIR = Path("examples/data")

# (json_name, grid_shapes, t_end, extra_params)
# Parameters: coupled_scalars needs B0, kappa, mg2, omegaP2
#             coupled_scattering needs g0, R, mPhi2, mChi2
#             em_3d has no free parameters
BENCHMARK_CASES: list[tuple[str, list[tuple[int, ...]], float, dict[str, float]]] = [
    # CVODE examples (pure ODE, no constraints)
    (
        "coupled_scalars.json",
        [(32, 32), (48, 48)],
        1.0,
        {"B0": 1.0, "kappa": 0.1, "mg2": 1.0, "omegaP2": 1.0},
    ),
    (
        "coupled_scattering.json",
        [(32, 32), (48, 48)],
        1.0,
        {"g0": 0.5, "R": 2.0, "mPhi2": 1.0, "mChi2": 1.0},
    ),
    # IDA example (constraints present)
    ("em_3d.json", [(8, 8, 8)], 0.5, {}),
]

N_REPEATS = 3  # median of N runs


def _make_patched_try():
    """Create a patched version of try_analytical_jacobian that forces FD for sparse tier."""
    import tidal.solver.analytical_jacobian as aj_mod

    orig_fn = aj_mod.try_analytical_jacobian

    def _patched(options, spec, layout, grid, bc, parameters, *, solver="ida"):
        from tidal.solver._types import DENSE_THRESHOLD, SPARSE_THRESHOLD

        n = layout.total_size
        if DENSE_THRESHOLD < n <= SPARSE_THRESHOLD:
            return False  # Force FD path
        return orig_fn(options, spec, layout, grid, bc, parameters, solver=solver)

    return _patched


def _run_simulation(
    json_path: Path,
    grid_shape: tuple[int, ...],
    t_end: float,
    params: dict[str, float],
    *,
    force_fd: bool = False,
) -> tuple[float, np.ndarray]:
    """Run a simulation and return (wall_time, final_state)."""
    import tidal.solver.analytical_jacobian as aj_mod
    from tidal.solver.grid import GridInfo
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem

    spec_data: dict[str, Any] = json.loads(json_path.read_text(encoding="utf-8"))
    spec = EquationSystem.from_dict(spec_data)

    ndim = len(grid_shape)
    bounds = tuple((0.0, 10.0) for _ in range(ndim))
    grid = GridInfo(
        bounds=bounds, shape=grid_shape, periodic=tuple(True for _ in range(ndim))
    )
    bc = "periodic"
    layout = StateLayout.from_spec(spec, grid.num_points)
    n_state = layout.total_size

    # Determine solver from equation structure
    has_constraints = any(
        eq.get("lhs", {}).get("order", {}).get("time", 2) == 0
        for eq in spec_data.get("equations", [])
    )

    # Build IC: small Gaussian on first field
    y0 = np.zeros(n_state)
    n_pts = grid.num_points
    center = np.array([5.0] * ndim)
    coords = np.meshgrid(
        *[
            np.linspace(b[0], b[1], s, endpoint=False)
            for b, s in zip(bounds, grid_shape, strict=False)
        ],
        indexing="ij",
    )
    r2 = sum((c - c0) ** 2 for c, c0 in zip(coords, center, strict=False))
    ic = 0.1 * np.exp(-r2 / 2.0)
    y0[:n_pts] = ic.ravel()

    t_span = (0.0, t_end)

    ctx = (
        patch.object(aj_mod, "try_analytical_jacobian", _make_patched_try())
        if force_fd
        else _nullcontext()
    )

    if has_constraints:
        from tidal.solver.ida import solve_ida

        with ctx:
            t0 = time.perf_counter()
            result = solve_ida(
                spec, grid, y0, t_span, bc=bc, parameters=params, num_snapshots=10
            )
            wall = time.perf_counter() - t0
    else:
        from tidal.solver.cvode import solve_cvode

        with ctx:
            t0 = time.perf_counter()
            result = solve_cvode(
                spec, grid, y0, t_span, bc=bc, parameters=params, num_snapshots=10
            )
            wall = time.perf_counter() - t0

    return wall, result["y"][-1]


class _nullcontext:
    """Minimal no-op context manager."""

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def main() -> None:
    print("=" * 78)
    print("Sparse Analytical Jacobian Benchmark")
    print("=" * 78)
    print()

    results: list[dict[str, Any]] = []

    for json_name, grid_shapes, t_end, params in BENCHMARK_CASES:
        json_path = EXAMPLES_DIR / json_name
        if not json_path.exists():
            print(f"SKIP {json_name} (not found)")
            continue

        for grid_shape in grid_shapes:
            # Compute expected state count
            spec_data = json.loads(json_path.read_text(encoding="utf-8"))
            n_slots = 0
            for eq in spec_data["equations"]:
                to = eq["lhs"]["order"]["time"]
                n_slots += 1
                if to >= 2:
                    n_slots += 1  # velocity slot
            n_pts = 1
            for s in grid_shape:
                n_pts *= s
            n_state = n_slots * n_pts

            label = (
                f"{json_name} @ {'x'.join(str(s) for s in grid_shape)} (N={n_state})"
            )
            print(f"\n--- {label} ---")

            # Warmup (1 run each)
            try:
                _, y_analytical = _run_simulation(
                    json_path, grid_shape, t_end, params, force_fd=False
                )
                _, y_fd = _run_simulation(
                    json_path, grid_shape, t_end, params, force_fd=True
                )
            except Exception as e:
                print(f"  ERROR during warmup: {e}")
                import traceback

                traceback.print_exc()
                continue

            # Correctness check
            max_diff = float(np.max(np.abs(y_analytical - y_fd)))
            print(f"  Solution max-diff (analytical vs FD): {max_diff:.2e}")

            # Timed runs
            times_analytical = []
            times_fd = []
            for _i in range(N_REPEATS):
                t_a, _ = _run_simulation(
                    json_path, grid_shape, t_end, params, force_fd=False
                )
                t_f, _ = _run_simulation(
                    json_path, grid_shape, t_end, params, force_fd=True
                )
                times_analytical.append(t_a)
                times_fd.append(t_f)

            med_a = statistics.median(times_analytical)
            med_f = statistics.median(times_fd)
            speedup = med_f / med_a if med_a > 0 else float("inf")

            print(f"  Analytical sparse: {med_a:.3f}s (median of {N_REPEATS})")
            print(f"  FD colored sparse: {med_f:.3f}s (median of {N_REPEATS})")
            print(f"  Speedup: {speedup:.2f}x")

            results.append(
                {
                    "example": json_name,
                    "grid": list(grid_shape),
                    "n_state": n_state,
                    "analytical_s": round(med_a, 4),
                    "fd_s": round(med_f, 4),
                    "speedup": round(speedup, 2),
                    "max_diff": max_diff,
                }
            )

    # Summary table
    print("\n" + "=" * 78)
    print("SUMMARY")
    print("=" * 78)
    print(
        f"{'Example':<30} {'Grid':<12} {'N':>7} {'Analytical':>10} {'FD':>10} {'Speedup':>8} {'MaxDiff':>10}"
    )
    print("-" * 92)
    for r in results:
        grid_str = "x".join(str(g) for g in r["grid"])
        print(
            f"{r['example']:<30} {grid_str:<12} {r['n_state']:>7} "
            f"{r['analytical_s']:>9.3f}s {r['fd_s']:>9.3f}s "
            f"{r['speedup']:>7.2f}x {r['max_diff']:>10.2e}"
        )

    # Acceptance criteria
    print()
    if results:
        speedups = [r["speedup"] for r in results]
        avg_speedup = statistics.mean(speedups)
        print(f"Average speedup: {avg_speedup:.2f}x")

        max_diffs = [r["max_diff"] for r in results]
        if max(max_diffs) < 1e-5:
            print(f"PASS: All solutions match to < 1e-5 (max: {max(max_diffs):.2e})")
        else:
            print(f"WARNING: Max solution diff = {max(max_diffs):.2e}")


if __name__ == "__main__":
    main()
