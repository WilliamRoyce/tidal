#!/usr/bin/env python3
"""Local benchmark: scipy vs JAX modal solver.

Runs both backends on the coupled-scalars and KG test specs and reports:
  - scipy mean wall time per solve call
  - JAX JIT compile time (first call)
  - JAX mean wall time per solve call (warm, post-JIT)
  - Max relative error between backends on y[-1]

Usage::

    uv run python scripts/local_jax_benchmark.py
    uv run python scripts/local_jax_benchmark.py --n-grid 128 --n-runs 20
    uv run python scripts/local_jax_benchmark.py --theory examples/data/coupled_scalars.json

Exits non-zero if max relative error > 1e-10.

Profiling evidence (from code analysis and benchmark runs):
  - Hot path: scipy.linalg.expm called per active mode in _evolve_per_mode_pade
    (modal.py:1996). For N=128 grid, coupled 4-slot system: ~129 rfft modes,
    4×4 expm per mode → O(N) matrix-exponential calls dominate wall time.
  - JAX path: jax.vmap(jax.scipy.linalg.expm) over all modes in one fused
    kernel + jax.lax.scan for the time loop → typical 3-10× speedup on CPU
    for N>=64, with JIT amortised over multiple likelihood calls in PolyChord.
  - Sparse-IC optimisation (scipy path): skips inactive modes. JAX vmap is
    always dense (all modes computed). For dense ICs the two paths are
    comparable; for single-mode plane-wave IC the scipy path wins on
    precompute cost. JAX wins on the scan evolution loop regardless.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Inline test specs (no JSON file required)
# ---------------------------------------------------------------------------

_KG_1D_SPEC: dict[str, Any] = {
    "metadata": {"source": "benchmark", "parameters": {"m2": 1.0}},
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
                ],
            },
        }
    ],
    "coupling": {"mass_matrix_symbolic": [["-m2"]]},
}

_COUPLED_SPEC: dict[str, Any] = {
    "metadata": {
        "source": "benchmark",
        "parameters": {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5},
    },
    "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
    "fields": [
        {"name": "phi_0", "index": 0, "is_dynamical": True},
        {"name": "chi_0", "index": 1, "is_dynamical": True},
    ],
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
                        "coefficient_symbolic": "-mPhi2",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ],
            },
        },
        {
            "field": "chi_0",
            "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "chi_0",
                        "coefficient_symbolic": "-mChi2",
                    },
                    {
                        "coefficient": -1.0,
                        "operator": "identity",
                        "field": "phi_0",
                        "coefficient_symbolic": "-gCpl",
                    },
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                ],
            },
        },
    ],
    "coupling": {"mass_matrix_symbolic": [["-mPhi2", "-gCpl"], ["-gCpl", "-mChi2"]]},
}


def _gaussian_ic(spec: object, grid: object) -> np.ndarray:
    from tidal.solver.state import StateLayout

    layout = StateLayout.from_spec(spec, grid.num_points)  # type: ignore[union-attr]
    y0 = np.zeros(layout.num_slots * grid.num_points)  # type: ignore[union-attr]
    x = np.linspace(
        grid.bounds[0][0],
        grid.bounds[0][1],
        grid.shape[0],
        endpoint=False,  # type: ignore[union-attr]
    )
    center = (grid.bounds[0][0] + grid.bounds[0][1]) / 2  # type: ignore[union-attr]
    y0[: grid.num_points] = 0.1 * np.exp(-((x - center) ** 2) / (2 * 1.5**2))  # type: ignore[union-attr]
    return y0


def _run_benchmark(
    spec: object,
    grid: object,
    y0: np.ndarray,
    params: dict[str, float],
    *,
    n_runs: int,
    num_snapshots: int,
    t_end: float,
    label: str,
) -> dict[str, Any]:
    """Benchmark scipy vs JAX on a single spec+grid combination."""
    from tidal.solver.modal import solve_modal

    try:
        import jax

        jax.config.update("jax_enable_x64", True)  # noqa: FBT003
        from tidal.solver.modal_jax import solve_modal_jax

        has_jax = True
    except ImportError:
        has_jax = False
        print(f"  [WARNING] JAX not installed — skipping JAX benchmark for {label}")
        print("           Install: uv sync --extra jax")

    print(f"\n{'=' * 60}")
    print(f"  {label}  (N={grid.shape}, snapshots={num_snapshots})")  # type: ignore[union-attr]
    print(f"{'=' * 60}")

    # --- JAX JIT warmup (done first to avoid contaminating scipy timing) ---
    if not has_jax:
        result_scipy = None
        result_jax = None
        jit_time = 0.0
    else:
        print("\n  JAX JIT warmup:")
        t_jit_start = time.perf_counter()
        result_jax = solve_modal_jax(  # type: ignore[arg-type]
            spec, grid, y0, (0.0, t_end), parameters=params, num_snapshots=num_snapshots
        )
        jit_time = time.perf_counter() - t_jit_start
        print(
            f"    JIT compile time: {jit_time * 1000:.0f} ms  (first call, not timed)"
        )

    # --- scipy warmup (after JAX JIT so caches are at steady state) --------
    _ = solve_modal(
        spec,
        grid,
        y0,
        (0.0, t_end),
        parameters=params,  # type: ignore[arg-type]
        num_snapshots=num_snapshots,
    )

    # --- Interleaved timing: alternating scipy/JAX to cancel monotonic drift -
    # Both backends see the same cache state / thermal conditions at each pair.
    times_scipy: list[float] = []
    times_jax: list[float] = []

    for _ in range(n_runs):
        # scipy
        t0 = time.perf_counter()
        result_scipy = solve_modal(
            spec,
            grid,
            y0,
            (0.0, t_end),
            parameters=params,  # type: ignore[arg-type]
            num_snapshots=num_snapshots,
        )
        times_scipy.append(time.perf_counter() - t0)

        if has_jax:
            t0 = time.perf_counter()
            result_jax = solve_modal_jax(  # type: ignore[arg-type]
                spec,
                grid,
                y0,
                (0.0, t_end),
                parameters=params,
                num_snapshots=num_snapshots,
            )
            times_jax.append(time.perf_counter() - t0)

    def _fmt_ms(vals: list[float]) -> str:
        arr = np.array(vals) * 1000.0
        p25, p50, p75 = np.percentile(arr, [25, 50, 75])
        return (
            f"median={p50:.2f} ms  IQR=[{p25:.2f},{p75:.2f}]  mean={arr.mean():.2f} ms"
            f"  (n={len(vals)})"
        )

    print(f"\n  scipy backend:  {_fmt_ms(times_scipy)}")

    if not has_jax:
        return {
            "label": label,
            "scipy_mean_ms": float(np.mean(times_scipy)) * 1000,
            "jax_available": False,
        }

    print(f"  JAX backend:    {_fmt_ms(times_jax)}  [post-JIT]")
    print(f"    JIT compile time: {jit_time * 1000:.0f} ms")

    mean_scipy = float(np.mean(times_scipy))
    mean_jax = float(np.mean(times_jax))
    median_scipy = float(np.median(times_scipy))
    median_jax = float(np.median(times_jax))

    # Speedup from medians (more robust than means for skewed distributions)
    speedup = median_scipy / median_jax
    print(f"\n  Speedup (scipy/JAX, median): {speedup:.2f}x")

    # --- Numerical agreement ---------------------------------------------
    y_s = result_scipy["y"][-1]
    y_j = result_jax["y"][-1]
    max_ref = float(np.max(np.abs(y_s)))
    abs_floor = 1e-15  # below this, treat solution as zero
    rel_tol = 1e-10
    if max_ref < abs_floor:
        max_err = float(np.max(np.abs(y_j - y_s)))
        print(f"  Max abs error (zero-IC): {max_err:.2e}")
        rel_err = max_err
    else:
        rel_err = float(np.max(np.abs(y_j - y_s))) / max_ref
        print(f"  Max relative error: {rel_err:.2e}  (tolerance: {rel_tol:.0e})")

    if rel_err > rel_tol:
        print(f"  [FAIL] Relative error {rel_err:.2e} exceeds {rel_tol:.0e} threshold")
    else:
        print("  [PASS] Relative error within tolerance")

    return {
        "label": label,
        "scipy_mean_ms": mean_scipy * 1000,
        "scipy_median_ms": median_scipy * 1000,
        "jax_jit_ms": jit_time * 1000,
        "jax_mean_ms": mean_jax * 1000,
        "jax_median_ms": median_jax * 1000,
        "speedup": speedup,
        "max_rel_error": rel_err,
        "pass": rel_err <= rel_tol,
        "jax_available": True,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--n-grid",
        type=int,
        default=128,
        help="Grid points per dimension (default: 128)",
    )
    parser.add_argument(
        "--n-runs",
        type=int,
        default=50,
        help="Interleaved benchmark pairs per backend (default: 50)",
    )
    parser.add_argument(
        "--snapshots",
        type=int,
        default=101,
        help="Snapshots per solve call (default: 101)",
    )
    parser.add_argument(
        "--t-end", type=float, default=10.0, help="Simulation end time (default: 10.0)"
    )
    parser.add_argument(
        "--theory",
        type=str,
        default=None,
        help="Optional: path to a TIDAL JSON spec to benchmark",
    )
    parser.add_argument(
        "--params",
        type=str,
        default=None,
        help="Parameters for --theory as key=val,key=val",
    )
    args = parser.parse_args(argv)

    sys.path.insert(0, str(Path(__file__).parent.parent))
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem

    results: list[dict[str, Any]] = []
    any_fail = False

    if args.theory:
        # Benchmark a user-supplied JSON spec
        with Path(args.theory).open(encoding="utf-8") as f:
            data = json.load(f)
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(
            shape=(args.n_grid,),
            bounds=[(0.0, 10.0)],
            periodic=[True],
        )
        params: dict[str, float] = {}
        if args.params:
            for kv in args.params.split(","):
                k, v = kv.split("=")
                params[k.strip()] = float(v.strip())
        from tidal.solver.state import StateLayout

        layout = StateLayout.from_spec(spec, grid.num_points)
        y0 = np.zeros(layout.num_slots * grid.num_points)
        x = np.linspace(0, 10, args.n_grid, endpoint=False)
        y0[: args.n_grid] = 0.1 * np.exp(-((x - 5.0) ** 2) / (2 * 1.5**2))
        r = _run_benchmark(
            spec,
            grid,
            y0,
            params,
            n_runs=args.n_runs,
            num_snapshots=args.snapshots,
            t_end=args.t_end,
            label=Path(args.theory).stem,
        )
        results.append(r)
        if r.get("jax_available") and not r.get("pass", True):
            any_fail = True
    else:
        # Built-in benchmarks
        for spec_dict, params, label in [
            (_KG_1D_SPEC, {"m2": 1.0}, f"KG-1D (N={args.n_grid})"),
            (
                _COUPLED_SPEC,
                {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5},
                f"Coupled scalars (N={args.n_grid})",
            ),
        ]:
            spec = EquationSystem.from_dict(spec_dict)
            grid = GridInfo(
                shape=(args.n_grid,),
                bounds=[(0.0, 10.0)],
                periodic=[True],
            )
            y0 = _gaussian_ic(spec, grid)
            r = _run_benchmark(
                spec,
                grid,
                y0,
                params,
                n_runs=args.n_runs,
                num_snapshots=args.snapshots,
                t_end=args.t_end,
                label=label,
            )
            results.append(r)
            if r.get("jax_available") and not r.get("pass", True):
                any_fail = True

    # Summary
    print(f"\n{'=' * 60}")
    print("  SUMMARY")
    print(f"{'=' * 60}")
    for r in results:
        if r.get("jax_available"):
            status = "PASS" if r.get("pass") else "FAIL"
            print(
                f"  [{status}] {r['label']}: "
                f"scipy_med={r['scipy_median_ms']:.2f}ms  "
                f"jax_med={r['jax_median_ms']:.2f}ms  "
                f"speedup={r['speedup']:.2f}x  "
                f"rel_err={r['max_rel_error']:.1e}"
            )
        else:
            print(f"  [SKIP] {r['label']}: JAX not available")

    if any_fail:
        print(
            "\n[FAIL] One or more benchmarks exceeded the 1e-10 relative error threshold."
        )
        return 1
    print("\n[PASS] All benchmarks within relative error tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
