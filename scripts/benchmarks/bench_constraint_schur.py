r"""Per-likelihood-call profile of solve_modal on the broken Phase E theory.

Headline metric: per-call wall time of ``solve_modal()`` on
``torsion_gertsenshtein_nonminimal_e_dual_gaussian.json`` at the canonical
Phase E geometry, varying BSM couplings (alpha1, alpha2, alpha3, delta1)
across calls.

Used as the GH #384 gate decision: which phase of the per-call work
dominates? Specifically:

- (a) ``CoefficientEvaluator.__init__`` + ``resolve()`` (coefficient
  evaluation, including the L2 spatial precompute).
- (b) ``_build_convolution_matrix_with_constraints`` fill phase
  (constraint and dynamical matrix assembly).
- (c) The LU solve at modal.py:2261 (``np.linalg.solve(K_cc, K_cd)``).
- (d) ``expm_multiply`` time evolution.

If (c) is < 20% of per-call time → fill is dominant; Phase A (coefficient
separable cache) is the right optimization to pursue.

If (c) is > 70% of per-call time → solve dominates; Phase A only saves a
small fraction. Pause and replan: switch to sparse LU
(``scipy.sparse.linalg.splu``) or preconditioned iterative solver instead.

Usage
-----
    uv run python scripts/benchmarks/bench_constraint_schur.py \\
        --calls 30 --output benchmark_results/bench_constraint_schur.json

Output: JSON with per-call timing, cProfile breakdown, headline ratios.
"""

from __future__ import annotations

import argparse
import cProfile
import json
import pstats
import time
from pathlib import Path

import numpy as np


def _load_spec_and_params() -> tuple[object, dict[str, float], object, object]:
    """Load the broken Phase E theory + initial state setup.

    Returns (spec, base_params, grid, y0).
    """
    from tidal.solver.grid import GridInfo
    from tidal.solver.state import StateLayout
    from tidal.symbolic.json_loader import EquationSystem

    json_path = Path(
        "examples/data/torsion_gertsenshtein_nonminimal_e_dual_gaussian.json"
    )
    with open(json_path) as f:
        spec_dict = json.load(f)
    spec = EquationSystem.from_dict(spec_dict)

    # GH #384 Phase A′: stash sampled-BSM-symbols on spec.metadata so the
    # modal solver's convolution-block cache engages on PolyChord-style
    # calls (varying alpha1..delta1, geometry fixed). Without this, the
    # cache misses on every call and we measure baseline performance.
    import dataclasses as _dc

    spec = _dc.replace(
        spec,
        metadata={
            **spec.metadata,
            "_inference_sampled_params": ("alpha1", "alpha2", "alpha3", "delta1"),
        },
    )

    # Canonical Phase E geometry (FROZEN — see scripts/hpc_submit_drafts/v3e_localised/_geometry.env)
    base_params = {
        "kappa": 1.0,
        "Bpeak": 0.01,
        "sigB": 5.0,
        "zc1": 50.0,
        "zc2": 150.0,
        # BSM couplings — sampled per-call; arbitrary initial values
        "alpha1": 0.0,
        "alpha2": 0.0,
        "alpha3": 0.0,
        "delta1": 0.0,
    }

    grid = GridInfo(shape=(128,), bounds=((0.0, 200.0),), periodic=(True,))
    layout = StateLayout.from_spec(spec, grid.num_points)
    y0 = np.zeros(layout.total_size, dtype=np.float64)
    # Small gaussian IC on h_5 (canonical Phase E sampled field)
    if "h_5" in layout.field_slot_map:
        h5_slot = layout.field_slot_map["h_5"]
        h5_start = h5_slot * grid.num_points
        x = np.linspace(0.0, 200.0, grid.num_points, endpoint=False)
        y0[h5_start : h5_start + grid.num_points] = 1e-4 * np.exp(
            -((x - 100.0) ** 2) / 100.0
        )
    return spec, base_params, grid, y0


def _time_single_call(spec, params, grid, y0) -> dict[str, float]:
    """One solve_modal call, returning timing dict."""
    from tidal.solver.modal import solve_modal

    t_start = time.perf_counter()
    _ = solve_modal(
        spec=spec,
        grid=grid,
        y0=y0,
        t_span=(0.0, 1.0),
        parameters=params,
        num_snapshots=3,
        return_eigendata=False,
    )
    t_wall = time.perf_counter() - t_start
    return {"wall_s": t_wall}


def _profile_breakdown(spec, params, grid, y0) -> dict[str, float]:
    """Run a single solve_modal under cProfile, return phase breakdown.

    Buckets functions by:
      - coeff_eval: CoefficientEvaluator and dependencies
      - fill: matrix builders (_build_convolution_*, _emit_term)
      - solve: np.linalg.solve / lu_factor / lu_solve
      - evolve: expm_multiply, expm
      - other: everything else
    """
    from tidal.solver.modal import solve_modal

    prof = cProfile.Profile()
    prof.enable()
    _ = solve_modal(
        spec=spec,
        grid=grid,
        y0=y0,
        t_span=(0.0, 1.0),
        parameters=params,
        num_snapshots=3,
        return_eigendata=False,
    )
    prof.disable()

    stats = pstats.Stats(prof)
    buckets = {
        "coeff_eval": 0.0,
        "fill": 0.0,
        "solve": 0.0,
        "evolve": 0.0,
        "other": 0.0,
    }
    total_tt = 0.0
    for func_key, (_cc, _nc, tt, _ct, _callers) in stats.stats.items():  # type: ignore[attr-defined]
        path, _line, fn_name = func_key
        total_tt += tt
        path_str = str(path)
        # Bucket assignment by function/path heuristic
        if "coefficients.py" in path_str or "_eval_utils" in path_str:
            buckets["coeff_eval"] += tt
        elif (
            "_build_convolution" in fn_name
            or "_emit_term" in fn_name
            or "_build_per_mode" in fn_name
            or "_build_evolution_matrices" in fn_name
        ):
            buckets["fill"] += tt
        elif fn_name in {"solve", "lu_factor", "lu_solve"} or (
            "linalg" in path_str and "expm" not in fn_name
        ):
            buckets["solve"] += tt
        elif "expm_multiply" in fn_name or "expm" in fn_name:
            buckets["evolve"] += tt
        else:
            buckets["other"] += tt

    # Total wall time from cProfile = sum of own-time (tt)
    buckets["_total"] = total_tt
    return buckets


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--calls", type=int, default=20, help="Number of solve_modal calls to time"
    )
    p.add_argument("--output", type=Path, default=None, help="JSON output path")
    p.add_argument(
        "--seed", type=int, default=42, help="RNG seed for varying BSM couplings"
    )
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    print("Loading spec and setting up...")
    spec, base_params, grid, y0 = _load_spec_and_params()
    print(f"  spec: {len(spec.equations)} equations")
    print(f"  grid: shape={grid.shape}, bounds={grid.bounds}")
    print(f"  layout total size: {y0.shape[0]}")

    # Warm-up call (JIT, caches)
    print("\nWarm-up call...")
    t_warm = _time_single_call(spec, base_params, grid, y0)
    print(f"  wall: {t_warm['wall_s']:.3f} s")

    # Time many calls at varying BSM couplings
    print(f"\nTiming {args.calls} calls with varying BSM couplings...")
    timings = []
    for i in range(args.calls):
        params = dict(base_params)
        params["alpha1"] = float(rng.uniform(-0.5, 0.5))
        params["alpha2"] = float(rng.uniform(-0.5, 0.5))
        params["alpha3"] = float(rng.uniform(-0.5, 0.5))
        params["delta1"] = float(rng.uniform(-0.5, 0.5))
        t = _time_single_call(spec, params, grid, y0)
        timings.append(t["wall_s"])
        if (i + 1) % 5 == 0:
            print(f"  {i + 1}/{args.calls}: {t['wall_s']:.3f} s")

    times = np.array(timings)
    summary = {
        "n_calls": int(args.calls),
        "mean_s": float(times.mean()),
        "p50_s": float(np.percentile(times, 50)),
        "p95_s": float(np.percentile(times, 95)),
        "min_s": float(times.min()),
        "max_s": float(times.max()),
        "std_s": float(times.std()),
    }
    print("\nPer-call wall time:")
    print(f"  mean = {summary['mean_s']:.3f} s")
    print(f"  p50  = {summary['p50_s']:.3f} s")
    print(f"  p95  = {summary['p95_s']:.3f} s")
    print(f"  min  = {summary['min_s']:.3f} s, max = {summary['max_s']:.3f} s")

    # cProfile breakdown on a single call
    print("\nProfiling one call for phase breakdown...")
    params = dict(base_params)
    params["alpha1"] = float(rng.uniform(-0.5, 0.5))
    params["delta1"] = float(rng.uniform(-0.5, 0.5))
    buckets = _profile_breakdown(spec, params, grid, y0)
    total = buckets.pop("_total")
    print(f"  cProfile own-time total: {total:.3f} s")
    print("  Phase breakdown (own-time, % of total):")
    for name, t in sorted(buckets.items(), key=lambda kv: -kv[1]):
        pct = 100 * t / total if total > 0 else 0
        print(f"    {name:12s}: {t:.3f} s ({pct:5.1f}%)")

    # Gate decision per plan
    solve_pct = 100 * buckets["solve"] / total if total > 0 else 0
    fill_pct = 100 * buckets["fill"] / total if total > 0 else 0
    coeff_pct = 100 * buckets["coeff_eval"] / total if total > 0 else 0
    print("\nGate decision:")
    print(f"  solve = {solve_pct:.1f}%")
    if solve_pct < 20:
        print("  → Phase A (fill-side optimization) is the right target. PROCEED.")
    elif solve_pct > 70:
        print("  → Solve dominates. Phase A only saves a small fraction.")
        print("    Pause and replan: target sparse LU or preconditioned iterative.")
    else:
        print(
            "  → Mixed. Phase A worthwhile but expect modest gains. Continue with caveat."
        )

    result = {
        "config": {
            "theory": "torsion_gertsenshtein_nonminimal_e_dual_gaussian",
            "n_grid": int(grid.shape[0]),
            "L": float(grid.bounds[0][1] - grid.bounds[0][0]),
            "t_end": 1.0,
            "n_snapshots": 3,
            "calls": int(args.calls),
        },
        "summary": summary,
        "breakdown_pct": {
            "coeff_eval": coeff_pct,
            "fill": fill_pct,
            "solve": solve_pct,
            "evolve": 100 * buckets["evolve"] / total if total > 0 else 0,
            "other": 100 * buckets["other"] / total if total > 0 else 0,
        },
        "breakdown_s": buckets | {"_total": total},
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2))
        print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()
