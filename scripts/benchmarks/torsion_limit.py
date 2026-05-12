r"""Einstein--Maxwell limit recovery from the torsion-extended theory.

Calibrates the symbolic-and-numerical pipeline jointly through a sharp
limit: the torsion-Gertsenshtein Lagrangian must reduce to pure
Einstein--Maxwell as the torsion-coupling parameter $\\xi$ approaches
zero. We sweep $\\xi$ over a logarithmic ladder and record

  $\\Delta(\\xi) = |P_\\mathrm{max}^\\mathrm{torsion}(\\xi)
                 - P_\\mathrm{max}^\\mathrm{EM}|$

against a single pure-Einstein--Maxwell reference run at the same
$(\\Bzero, t_\\mathrm{end})$. The expected behaviour is $\\Delta \\to 0$
linearly in $\\xi$ for the propagating-PGT theory. The test engages every
stage of the pipeline: the symbolic derivation must drop torsion couplings
analytically in the limit, the numerical solver must reproduce the
reduction.

Serves:   manuscript/sections/appendices/validation.tex (App D, calibration 2)
Consumes: scripts/figures/figD_torsion_limit.py
Writes:   benchmark_results/canonical/torsion_limit.json
"""

from __future__ import annotations

import argparse
import csv
import datetime
import json
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "torsion_limit.json"

TORSION_SPEC = EXAMPLES / "data" / "torsion_gertsenshtein_minimal_propagating.json"
EM_SPEC = EXAMPLES / "data" / "gertsenshtein.json"

# Physics parameters held constant across the sweep.
KAPPA = 1.0
B0_FIXED = 0.05
KWAVE = 2.0106
BOUNDS = (0.0, 100.0)
T_END = 50.0
DELTA1_FIXED = 0.0  # propagating-PGT defaults to delta1 = 0
BETA_DEFAULTS = {"beta1": 1.0, "beta2": 1.0, "beta3": 1.0}

# xi-sweep: log-spaced; smoke mode is coarse for local sanity check.
FULL_XI_VALUES = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
FULL_GRID_N = 512
SMOKE_XI_VALUES = [1e-1, 1e-3, 1e-5]
SMOKE_GRID_N = 128


def _git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, stderr=subprocess.DEVNULL
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


def _f(row: dict, *keys: str) -> float | None:
    for k in keys:
        v = row.get(k)
        if not v:
            continue
        try:
            return float(v)
        except ValueError:
            continue
    return None


def _common_args(
    out_dir: Path,
    *,
    json_spec: Path,
    grid_n: int,
    parallel: int,
    ic_component: str,
    source: str,
    target: str,
    extra_params: dict[str, float] | None = None,
) -> list[str]:
    cmd = [
        "tidal",
        "sweep",
        str(json_spec),
        "--measure",
        "conversion",
        "--source",
        source,
        "--target",
        target,
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
        ic_component,
        "--t-end",
        f"{T_END}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_FIXED}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--parallel",
        str(parallel),
        "--force",
    ]
    if extra_params:
        for k, v in extra_params.items():
            cmd.extend(["--param", f"{k}={v}"])
    return cmd


def _run_torsion_xi_sweep(
    out_dir: Path, *, xi_values: list[float], grid_n: int, parallel: int
) -> list[dict]:
    extra = dict(BETA_DEFAULTS)
    extra["delta1"] = DELTA1_FIXED
    cmd = _common_args(
        out_dir,
        json_spec=TORSION_SPEC,
        grid_n=grid_n,
        parallel=parallel,
        ic_component="h_5",  # h_x channel — common to EM and propagating PGT
        source="h_5",
        target="a_1",
        extra_params=extra,
    )
    xi_csv = ",".join(f"{x}" for x in xi_values)
    cmd += ["--sweep", f"xi={xi_csv}"]
    print(f"[torsion_limit] xi sweep: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return _read_results_csv(out_dir)


def _run_em_reference(out_dir: Path, *, grid_n: int) -> dict:
    """Single Einstein-Maxwell run + measure, no sweep."""
    sim_cmd = [
        "tidal",
        "simulate",
        str(EM_SPEC),
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
        "h_5",
        "--t-end",
        f"{T_END}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_FIXED}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    print(f"[torsion_limit] EM simulate: {' '.join(sim_cmd)}", flush=True)
    subprocess.run(sim_cmd, check=True, cwd=REPO_ROOT)

    meas_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "conversion",
        "--source",
        "h_5",
        "--target",
        "a_1",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_FIXED}",
        "--json",
        "--quiet",
    ]
    print(f"[torsion_limit] EM measure: {' '.join(meas_cmd)}", flush=True)
    out = subprocess.check_output(meas_cmd, cwd=REPO_ROOT)
    meas = json.loads(out)
    conv = meas.get("conversion", {})
    pk = float(conv.get("peak_probability", 0.0))
    return {"P_max": pk, "P_final": pk}


def _analyse(
    torsion_rows: list[dict], em_ref: dict, xi_values: list[float]
) -> list[dict]:
    out: list[dict] = []
    em_pmax = em_ref["P_max"] or 0.0
    em_pfinal = em_ref["P_final"] or 0.0
    # Map xi to the row that has it; tidal sweep may emit rows ordered.
    rows_by_xi = {}
    for row in torsion_rows:
        xi = _f(row, "xi")
        if xi is None:
            continue
        rows_by_xi[round(xi, 16)] = row
    for xi in sorted(xi_values, reverse=True):
        row = rows_by_xi.get(round(xi, 16))
        if row is None:
            continue
        p_max = _f(row, "P_max") or 0.0
        p_final = _f(row, "P_final") or 0.0
        out.append(
            {
                "xi": xi,
                "P_max_sim": p_max,
                "P_final_sim": p_final,
                "P_max_EM": em_pmax,
                "P_final_EM": em_pfinal,
                "abs_diff_P_max": abs(p_max - em_pmax),
                "abs_diff_P_final": abs(p_final - em_pfinal),
            }
        )
    return out


def _summary(rows: list[dict]) -> dict:
    if not rows:
        return {"n_xi_points": 0}
    diffs_max = np.array([r["abs_diff_P_max"] for r in rows])
    diffs_final = np.array([r["abs_diff_P_final"] for r in rows])
    xis = np.array([r["xi"] for r in rows])
    # Linear-in-xi fit: log|diff| ~ slope * log(xi) + intercept
    mask = (xis > 0) & (diffs_max > 0)
    slope_max = None
    if mask.sum() >= 2:
        slope_max, _ = np.polyfit(np.log10(xis[mask]), np.log10(diffs_max[mask]), 1)
        slope_max = float(slope_max)
    return {
        "n_xi_points": len(rows),
        "min_xi": float(xis.min()),
        "max_xi": float(xis.max()),
        "diff_P_max_at_smallest_xi": float(diffs_max[np.argmin(xis)]),
        "diff_P_final_at_smallest_xi": float(diffs_final[np.argmin(xis)]),
        "log_log_slope_P_max_vs_xi": slope_max,
    }


def run(*, smoke: bool, parallel: int, work_dir: Path) -> dict:
    xi_values = SMOKE_XI_VALUES if smoke else FULL_XI_VALUES
    grid_n = SMOKE_GRID_N if smoke else FULL_GRID_N
    work_dir.mkdir(parents=True, exist_ok=True)

    torsion_dir = work_dir / "torsion_xi_sweep"
    if torsion_dir.exists():
        shutil.rmtree(torsion_dir)
    torsion_rows = _run_torsion_xi_sweep(
        torsion_dir, xi_values=xi_values, grid_n=grid_n, parallel=parallel
    )

    em_dir = work_dir / "em_reference"
    if em_dir.exists():
        shutil.rmtree(em_dir)
    em_ref = _run_em_reference(em_dir, grid_n=grid_n)

    rows = _analyse(torsion_rows, em_ref, xi_values)
    return {
        "metadata": _metadata(
            {
                "xi_values": xi_values,
                "grid_n": grid_n,
                "kappa": KAPPA,
                "B0": B0_FIXED,
                "kwave": KWAVE,
                "bounds": list(BOUNDS),
                "t_end": T_END,
                "delta1": DELTA1_FIXED,
                "beta_defaults": BETA_DEFAULTS,
                "fd_order": 4,
                "smoke": smoke,
            }
        ),
        "summary": _summary(rows),
        "em_reference": em_ref,
        "xi_sweep": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="torsion_limit_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, parallel=args.parallel, work_dir=work)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
