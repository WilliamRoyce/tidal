r"""In-process modal vs IDA vs CVODE timing on the GH #367 reproducer.

Question this benchmark answers: how much faster would the modal solver be
than IDA/CVODE on the position-dependent + periodic Gertsenshtein case if
the discretization artifact were fixed? The current Phase 1 auto-route to
IDA is correct but potentially slow at PolyChord scale (10^4 likelihood
calls per chain) — this measurement gates the algorithm-research work
(Phase B-D of the plan).

Run on the safe modal regime (t_end small enough that k_max*t < 2) where
the modal solver is still accurate; that establishes the upper bound on
speedup if the algorithm fix lands. Also runs at the broken regime
(t_end=20) to document the wall-time-vs-correctness tradeoff.

Usage:
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
        uv run python scripts/benchmarks/modal_vs_ida_inprocess.py
"""

from __future__ import annotations

import json
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tidal.solver.cvode import solve_cvode
from tidal.solver.grid import GridInfo
from tidal.solver.ida import solve_ida
from tidal.solver.modal import solve_modal
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

THEORY_JSON = REPO_ROOT / "examples/data/gertsenshtein_e0_dual_gaussian.json"
PARAMS: dict[str, float] = {
    "kappa": 1.0,
    "Bpeak": 0.01,
    "sigB": 5.0,
    "zc1": 25.0,
    "zc2": 75.0,
}
IC_AMPLITUDE = 1e-2
IC_WIDTH = 5.0
IC_CENTER = 25.0
DOMAIN_LENGTH = 100.0
IC_SLOT = "h_5"

N_REPS = 5
N_SNAPSHOTS = 11


def build_y0(spec: EquationSystem, grid: GridInfo) -> np.ndarray:
    """Build the Gaussian-in-h_5 initial state matching the GH #367 reproducer."""
    layout = StateLayout.from_spec(spec, grid.num_points)
    y0 = np.zeros(layout.total_size, dtype=float)
    x_axis = np.linspace(0.0, DOMAIN_LENGTH, grid.num_points, endpoint=False)
    gauss = IC_AMPLITUDE * np.exp(-((x_axis - IC_CENTER) ** 2) / (2 * IC_WIDTH**2))
    slot_idx = layout.slot_name_to_idx[IC_SLOT]
    y0[layout.slot_slice(slot_idx)] = gauss
    return y0


def build_grid(n: int) -> GridInfo:
    return GridInfo(
        shape=(n,),
        bounds=[(0.0, DOMAIN_LENGTH)],
        periodic=(True,),
    )


def k_max_for_n(n: int) -> float:
    """Nyquist wavenumber for periodic discretisation of [0, L]."""
    dx = DOMAIN_LENGTH / n
    return math.pi / dx


def time_solver(
    fn: Any,
    spec: EquationSystem,
    grid: GridInfo,
    y0: np.ndarray,
    t_span: tuple[float, float],
    parameters: dict[str, float],
    n_reps: int,
) -> tuple[float, float, dict[str, Any] | None]:
    """Return (median_seconds, h5_peak, result) over n_reps timed calls.

    h5_peak is taken from the LAST call so we can sanity-check correctness.
    If the solver raises SimulationDivergedError (modal at broken regime),
    return (median_time_until_failure, nan, None).
    """
    from tidal.solver._exceptions import SimulationDivergedError

    layout = StateLayout.from_spec(spec, grid.num_points)
    times: list[float] = []
    result: dict[str, Any] | None = None
    diverged = False
    for _ in range(n_reps):
        t0 = time.perf_counter()
        try:
            result = fn(
                spec,
                grid,
                y0,
                t_span=t_span,
                parameters=parameters,
                num_snapshots=N_SNAPSHOTS,
            )
        except SimulationDivergedError:
            diverged = True
        times.append(time.perf_counter() - t0)
    if diverged or result is None:
        return statistics.median(times), float("nan"), None
    h5_slice = layout.slot_slice(layout.slot_name_to_idx[IC_SLOT])
    h5_peak = float(np.max(np.abs(result["y"][-1][h5_slice])))
    return statistics.median(times), h5_peak, result


def main() -> None:
    print(f"Loading {THEORY_JSON.relative_to(REPO_ROOT)}")
    with Path(THEORY_JSON).open(encoding="utf-8") as f:
        spec = EquationSystem.from_dict(json.load(f))

    rows: list[dict[str, Any]] = []
    configs = [
        # (N, t_end, label)
        (64, 0.4, "safe-modal-N64"),  # k_max*t = 2.01 * 0.4 = 0.8 (well within safe)
        (64, 1.0, "boundary-modal-N64"),  # k_max*t = 2.01 (4% error per handoff table)
        (64, 20.0, "broken-modal-N64"),  # k_max*t = 40 (1e9x error per handoff)
        (128, 0.2, "safe-modal-N128"),  # k_max*t = 4.02 * 0.2 = 0.8
        (128, 1.0, "boundary-modal-N128"),  # k_max*t = 4.02
        (256, 0.1, "safe-modal-N256"),  # k_max*t = 8.04 * 0.1 = 0.8
    ]

    for n, t_end, label in configs:
        kmax = k_max_for_n(n)
        kt = kmax * t_end
        print(f"\n=== {label}: N={n}, t_end={t_end}, k_max*t={kt:.2f} ===")
        grid = build_grid(n)
        y0 = build_y0(spec, grid)

        modal_t, modal_h5, _ = time_solver(
            solve_modal, spec, grid, y0, (0.0, t_end), PARAMS, N_REPS
        )
        print(f"  modal  : {modal_t * 1000:7.2f} ms   h5_peak = {modal_h5:.4e}")

        cvode_t, cvode_h5, _ = time_solver(
            solve_cvode, spec, grid, y0, (0.0, t_end), PARAMS, N_REPS
        )
        print(f"  cvode  : {cvode_t * 1000:7.2f} ms   h5_peak = {cvode_h5:.4e}")

        ida_t, ida_h5, _ = time_solver(
            solve_ida, spec, grid, y0, (0.0, t_end), PARAMS, N_REPS
        )
        print(f"  ida    : {ida_t * 1000:7.2f} ms   h5_peak = {ida_h5:.4e}")

        if math.isnan(modal_h5):
            print("  modal/cvode err: DIVERGED (safety net caught the run)")
            modal_err = float("nan")
        else:
            modal_err = abs(modal_h5 - cvode_h5) / max(cvode_h5, 1e-300)
            print(f"  modal/cvode err: {modal_err * 100:.2f}%")
        print(
            f"  ratios:  modal/cvode = {modal_t / cvode_t:.2f}x   "
            f"modal/ida = {modal_t / ida_t:.2f}x   "
            f"cvode/ida = {cvode_t / ida_t:.2f}x"
        )

        rows.append(
            {
                "label": label,
                "N": n,
                "t_end": t_end,
                "kmax_t": kt,
                "modal_ms": modal_t * 1000,
                "cvode_ms": cvode_t * 1000,
                "ida_ms": ida_t * 1000,
                "modal_h5": modal_h5,
                "cvode_h5": cvode_h5,
                "ida_h5": ida_h5,
                "modal_err_pct": modal_err * 100,
                "modal_vs_cvode": modal_t / cvode_t,
                "modal_vs_ida": modal_t / ida_t,
            }
        )

    out_path = REPO_ROOT / "benchmark_results" / "modal_vs_ida_inprocess.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with Path(out_path).open("w", encoding="utf-8") as f:
        json.dump({"reps": N_REPS, "rows": rows}, f, indent=2)
    print(f"\nWrote {out_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
