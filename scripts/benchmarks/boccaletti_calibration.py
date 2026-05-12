r"""Boccaletti-kernel calibration for the Einstein--Maxwell baseline.

Calibrates TIDAL's reproduction of the Gertsenshtein photon-graviton
conversion probability $P = \\sin^2(\\kappa B_0 D / 2)$ in two regimes:

  (a) B0 sweep at fixed propagation distance ($D = t_\\mathrm{end}$);
  (b) multi-resolution convergence at a single regime point.

For each sweep cell the simulated $P_\\mathrm{final}$ and $P_\\mathrm{max}$
are recorded with their residual against the analytic baselines:

  * $P_\\mathrm{final}^\\mathrm{theory} = \\sin^2(\\kappa B_0 t/2)$
  * $P_\\mathrm{max}^\\mathrm{theory}   = \\sin^2(\\kappa B_0 t/2) \\,
    k^2 / (k^2 + \\kappa^2 B_0^2)$

The $P_\\mathrm{max}$ correction comes from the graviton effective-mass
contribution $\\kappa^2 B_0^2 / 2$ that the massless photon cannot absorb;
see comments in `examples/gertsenshtein/sweep_B0.sh`.

A planned follow-up extends this script to (c) the localised-$B_0$ regime
via `theory_localized.toml`, once that variant's JSON spec is derived.

Serves:   manuscript/sections/appendices/validation.tex (BoccalettiCalibration)
Consumes: scripts/figures/figD_boccaletti_calibration.py,
          scripts/figures/fig1_boccaletti_validation.py
Writes:   benchmark_results/canonical/boccaletti_calibration.json
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
import math
import operator
import platform
import shutil
import socket
import subprocess  # noqa: S404
import tempfile
from pathlib import Path

import numpy as np
import scipy

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = REPO_ROOT / "examples"
DEFAULT_OUT = (
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_calibration.json"
)

# Physical / numerical parameters (held constant across sweep cells).
KAPPA = 1.0
KWAVE = 2.0106  # 32 wavelengths in the [0, 100] domain
BOUNDS = (0.0, 100.0)
T_END = 50.0
B0_FIXED = 0.05  # regime point for convergence study

# Full HPC-scale sweep parameters.
FULL_B0_LO, FULL_B0_HI, FULL_B0_N = 0.005, 0.25, 40
FULL_GRID_N = 512
FULL_CONVERGE_N = [64, 128, 256, 512, 1024, 2048]

# Smoke-mode parameters (local sanity check, < 2 min).
SMOKE_B0_LO, SMOKE_B0_HI, SMOKE_B0_N = 0.05, 0.20, 4
SMOKE_GRID_N = 128
SMOKE_CONVERGE_N = [64, 128]


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


def _metadata(parameters: dict) -> dict:
    return {
        "timestamp": datetime.datetime.now(datetime.UTC).isoformat(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "git_sha": _git_sha(),
        "parameters": parameters,
    }


def _read_results_csv(run_dir: Path) -> list[dict]:
    csv_path = run_dir / "results.csv"
    if not csv_path.exists():
        msg = f"sweep result not found: {csv_path}"
        raise FileNotFoundError(msg)
    with csv_path.open() as fh:
        return list(csv.DictReader(fh))


def _common_args(out_dir: Path, *, grid_n: int, parallel: int) -> list[str]:
    return [
        "tidal",
        "sweep",
        str(EXAMPLES / "data" / "gertsenshtein.json"),
        "--measure",
        "conversion",
        "--source",
        "h_7",
        "--target",
        "a_2",
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{BOUNDS[0]}:{BOUNDS[1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        f"{KWAVE}",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        "h_7",
        "--t-end",
        f"{T_END}",
        "--param",
        f"kappa={KAPPA}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--parallel",
        str(parallel),
        "--force",
    ]


def _run_b0_sweep(
    out_dir: Path, *, b0_lo: float, b0_hi: float, b0_n: int, grid_n: int, parallel: int
) -> list[dict]:
    cmd = _common_args(out_dir, grid_n=grid_n, parallel=parallel)
    cmd += ["--sweep", f"B0={b0_lo}:{b0_hi}:{b0_n}"]
    print(f"[boccaletti_calibration] B0 sweep: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return _read_results_csv(out_dir)


def _run_convergence(
    out_dir: Path, *, n_values: list[int], parallel: int
) -> list[dict]:
    cmd = _common_args(out_dir, grid_n=n_values[0], parallel=parallel)
    cmd += ["--converge", ",".join(str(n) for n in n_values)]
    cmd += ["--param", f"B0={B0_FIXED}"]
    print(f"[boccaletti_calibration] convergence: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return _read_results_csv(out_dir)


def _analytic_bare(b0: float, t_end: float = T_END) -> float:
    return float(math.sin(0.5 * KAPPA * b0 * t_end) ** 2)


def _analytic_pmax_corrected(b0: float, t_end: float = T_END) -> float:
    bare = _analytic_bare(b0, t_end)
    cap = KWAVE**2 / (KWAVE**2 + (KAPPA * b0) ** 2)
    return bare * cap


def _f(row: dict, *keys: str, default: float | None = None) -> float | None:
    for k in keys:
        v = row.get(k)
        if not v:
            continue
        try:
            return float(v)
        except ValueError:
            continue
    return default


def _b0_residuals(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    for row in rows:
        b0 = _f(row, "B0")
        p_final = _f(row, "P_final")
        p_max = _f(row, "P_max")
        if b0 is None or p_final is None or p_max is None:
            continue
        bare = _analytic_bare(b0)
        corrected = _analytic_pmax_corrected(b0)
        out.append(
            {
                "B0": b0,
                "t_end": T_END,
                "P_final_sim": p_final,
                "P_max_sim": p_max,
                "P_final_analytic": bare,
                "P_max_analytic_corrected": corrected,
                "residual_final": p_final - bare,
                "residual_max": p_max - corrected,
            }
        )
    return out


def _convergence_results(rows: list[dict]) -> list[dict]:
    out: list[dict] = []
    p_final_target = _analytic_bare(B0_FIXED)
    for row in rows:
        n = _f(row, "grid_shape", "grid_n", "N", "Nx", "Ngrid")
        p_max = _f(row, "P_max")
        p_final = _f(row, "P_final")
        if n is None or p_final is None:
            continue
        out.append(
            {
                "N": int(n),
                "B0": B0_FIXED,
                "t_end": T_END,
                "P_max_sim": p_max,
                "P_final_sim": p_final,
                "P_final_analytic": p_final_target,
                "abs_error_final": abs(p_final - p_final_target),
            }
        )
    out.sort(key=operator.itemgetter("N"))
    return out


def _summary(b0_rows: list[dict], conv_rows: list[dict]) -> dict:
    res_final = np.array([r["residual_final"] for r in b0_rows])
    res_max = np.array([r["residual_max"] for r in b0_rows])
    return {
        "max_abs_residual_final": float(np.max(np.abs(res_final)))
        if res_final.size
        else None,
        "rms_residual_final": float(np.sqrt(np.mean(res_final**2)))
        if res_final.size
        else None,
        "max_abs_residual_max": float(np.max(np.abs(res_max)))
        if res_max.size
        else None,
        "rms_residual_max": float(np.sqrt(np.mean(res_max**2)))
        if res_max.size
        else None,
        "n_b0_points": int(res_final.size),
        "min_abs_error_at_highest_N": (
            float(min(r["abs_error_final"] for r in conv_rows)) if conv_rows else None
        ),
        "highest_N_in_convergence": (
            int(max(r["N"] for r in conv_rows)) if conv_rows else None
        ),
    }


def run(*, smoke: bool, parallel: int, work_dir: Path) -> dict:
    b0_lo, b0_hi, b0_n = (
        (SMOKE_B0_LO, SMOKE_B0_HI, SMOKE_B0_N)
        if smoke
        else (FULL_B0_LO, FULL_B0_HI, FULL_B0_N)
    )
    grid_n = SMOKE_GRID_N if smoke else FULL_GRID_N
    n_values = SMOKE_CONVERGE_N if smoke else FULL_CONVERGE_N

    work_dir.mkdir(parents=True, exist_ok=True)

    sweep_dir = work_dir / "b0_sweep"
    if sweep_dir.exists():
        shutil.rmtree(sweep_dir)
    b0_raw = _run_b0_sweep(
        sweep_dir, b0_lo=b0_lo, b0_hi=b0_hi, b0_n=b0_n, grid_n=grid_n, parallel=parallel
    )
    b0_rows = _b0_residuals(b0_raw)

    conv_dir = work_dir / "convergence"
    if conv_dir.exists():
        shutil.rmtree(conv_dir)
    conv_raw = _run_convergence(
        conv_dir, n_values=n_values, parallel=max(1, parallel // 2)
    )
    conv_rows = _convergence_results(conv_raw)

    return {
        "metadata": _metadata(
            {
                "b0_range": [b0_lo, b0_hi],
                "b0_n_points": b0_n,
                "convergence_n_values": n_values,
                "grid_n_sweep": grid_n,
                "kappa": KAPPA,
                "kwave": KWAVE,
                "bounds": list(BOUNDS),
                "t_end": T_END,
                "b0_convergence_point": B0_FIXED,
                "fd_order": 4,
                "smoke": smoke,
            }
        ),
        "summary": _summary(b0_rows, conv_rows),
        "b0_sweep": b0_rows,
        "convergence": conv_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=DEFAULT_OUT, help="output JSON path"
    )
    parser.add_argument(
        "--smoke", action="store_true", help="run a small sweep for local verification"
    )
    parser.add_argument(
        "--parallel", type=int, default=4, help="parallel workers for tidal sweep"
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=None,
        help="scratch dir for sweep outputs (defaults to a tempdir)",
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="boccaletti_calibration_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, parallel=args.parallel, work_dir=work)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
