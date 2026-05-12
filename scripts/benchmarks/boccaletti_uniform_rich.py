"""Boccaletti-kernel reproduction in the uniform-B0 regime, three-overlay.

Calibrates the Einstein--Maxwell Gertsenshtein simulation against three
distinct analytic predictions of the conversion probability on a
uniform background:

  (bare)         P_final = sin^2(kappa B0 t/2)
  (energy-cap)   P_max   = sin^2(kappa B0 t/2) * k^2/(k^2 + kappa^2 B0^2 / 2)
                            -- accounts for graviton effective-mass energy
                               the massless photon cannot absorb;
                               derivation in examples/gertsenshtein/sweep_B0.sh
                               and docs/tex/gertsenshtein_formula.tex eq. P-corrected
  (Raffelt-Stodolsky)
                 P_RS    = g^2/(g^2 + Delta^2) * sin^2(sqrt(g^2+Delta^2) D)
                            with g = kappa B0 / 2 and Delta = -kappa^2 B0^2 / (4 omega)
                            -- the two-mode detuning formula from
                               manuscript perturbative_regime.tex eq. 344--347.

Sweeps (B0, t_end) on a 2D grid (HPC-parallel) and records P_final, P_max,
and the residual against each of the three analytic predictions per cell.
Adds a multi-resolution N-sweep at a regime point to certify the
discretisation does not dominate the residual.

Serves:   manuscript/sections/appendices/validation.tex (App D §1)
Consumes: scripts/figures/figD_boccaletti_uniform.py
Writes:   benchmark_results/canonical/boccaletti_uniform_rich.json
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
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_uniform_rich.json"
)

KAPPA = 1.0
KWAVE = 2.0106  # 32 wavelengths in [0, 100]
BOUNDS = (0.0, 100.0)
SOURCE = "h_5"
TARGET = "a_1"
IC_COMPONENT = "h_5"

# Full HPC scale.
FULL_B0_LO, FULL_B0_HI, FULL_B0_N = 0.005, 0.20, 16
FULL_T_END_VALUES = [25.0, 50.0, 75.0, 100.0]
FULL_GRID_N = 1024
FULL_N_CONVERGE = [128, 256, 512, 1024, 2048]

# Smoke mode.
SMOKE_B0_LO, SMOKE_B0_HI, SMOKE_B0_N = 0.05, 0.20, 4
SMOKE_T_END_VALUES = [25.0, 50.0]
SMOKE_GRID_N = 128
SMOKE_N_CONVERGE = [64, 128]


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


def _sweep_b0(
    *,
    t_end: float,
    b0_lo: float,
    b0_hi: float,
    b0_n: int,
    grid_n: int,
    parallel: int,
    out_dir: Path,
) -> list[dict]:
    cmd = [
        "tidal",
        "sweep",
        str(EXAMPLES / "data" / "gertsenshtein.json"),
        "--measure",
        "conversion",
        "--source",
        SOURCE,
        "--target",
        TARGET,
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
        IC_COMPONENT,
        "--t-end",
        f"{t_end}",
        "--param",
        f"kappa={KAPPA}",
        "--fd-order",
        "4",
        "--sweep",
        f"B0={b0_lo}:{b0_hi}:{b0_n}",
        "--output",
        str(out_dir),
        "--parallel",
        str(parallel),
        "--force",
    ]
    print(
        f"[boccaletti_uniform_rich] t_end={t_end}: {' '.join(cmd[:6])} ...", flush=True
    )
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return _read_results_csv(out_dir)


def _converge(
    *, n_values: list[int], b0: float, t_end: float, parallel: int, out_dir: Path
) -> list[dict]:
    cmd = [
        "tidal",
        "sweep",
        str(EXAMPLES / "data" / "gertsenshtein.json"),
        "--measure",
        "conversion",
        "--source",
        SOURCE,
        "--target",
        TARGET,
        "--grid-shape",
        str(n_values[0]),
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
        IC_COMPONENT,
        "--t-end",
        f"{t_end}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={b0}",
        "--fd-order",
        "4",
        "--converge",
        ",".join(str(n) for n in n_values),
        "--output",
        str(out_dir),
        "--parallel",
        str(parallel),
        "--force",
    ]
    print(f"[boccaletti_uniform_rich] N-converge: {' '.join(cmd[:6])} ...", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)
    return _read_results_csv(out_dir)


def _analytic_bare(b0: float, t: float) -> float:
    return float(math.sin(0.5 * KAPPA * b0 * t) ** 2)


def _analytic_energy_cap(b0: float, t: float) -> float:
    bare = _analytic_bare(b0, t)
    cap = KWAVE**2 / (KWAVE**2 + 0.5 * (KAPPA * b0) ** 2)
    return bare * cap


def _analytic_raffelt_stodolsky(b0: float, t: float) -> float:
    """Two-mode detuned formula from perturbative_regime.tex.

    g  = kappa B0 / 2          (mixing strength, ∝ B0)
    Delta = -kappa^2 B0^2 / (4 omega)  with omega ≈ k (massless limit)
                               (detuning from graviton effective mass, ∝ B0²)
    P  = g²/(g² + Delta²) sin²(sqrt(g² + Delta²) t)
    """
    g = 0.5 * KAPPA * b0
    omega = KWAVE
    delta = -(KAPPA**2) * (b0**2) / (4.0 * omega)
    g2 = g * g
    d2 = delta * delta
    omega_m = math.sqrt(g2 + d2)
    if omega_m == 0.0:
        return 0.0
    return float(g2 / (g2 + d2) * math.sin(omega_m * t) ** 2)


def _b0_residuals(rows: list[dict], t_end: float) -> list[dict]:
    out: list[dict] = []
    for row in rows:
        b0 = _f(row, "B0")
        p_final = _f(row, "P_final")
        p_max = _f(row, "P_max")
        if b0 is None or p_final is None or p_max is None:
            continue
        bare = _analytic_bare(b0, t_end)
        capped = _analytic_energy_cap(b0, t_end)
        rs = _analytic_raffelt_stodolsky(b0, t_end)
        out.append(
            {
                "B0": b0,
                "t_end": t_end,
                "P_final_sim": p_final,
                "P_max_sim": p_max,
                "P_bare": bare,
                "P_energy_cap": capped,
                "P_raffelt_stodolsky": rs,
                "residual_final_vs_bare": p_final - bare,
                "residual_max_vs_capped": p_max - capped,
                "residual_max_vs_RS": p_max - rs,
            }
        )
    return out


def _converge_rows(rows: list[dict], b0: float, t_end: float) -> list[dict]:
    out: list[dict] = []
    bare = _analytic_bare(b0, t_end)
    for row in rows:
        n = _f(row, "grid_shape", "grid_n", "N")
        p_final = _f(row, "P_final")
        p_max = _f(row, "P_max")
        if n is None or p_final is None:
            continue
        out.append(
            {
                "N": int(n),
                "B0": b0,
                "t_end": t_end,
                "P_final_sim": p_final,
                "P_max_sim": p_max,
                "P_bare": bare,
                "abs_error_final": abs(p_final - bare),
            }
        )
    out.sort(key=operator.itemgetter("N"))
    return out


def run(*, smoke: bool, parallel: int, work_dir: Path) -> dict:
    b0_lo, b0_hi, b0_n = (
        (SMOKE_B0_LO, SMOKE_B0_HI, SMOKE_B0_N)
        if smoke
        else (FULL_B0_LO, FULL_B0_HI, FULL_B0_N)
    )
    t_values = SMOKE_T_END_VALUES if smoke else FULL_T_END_VALUES
    grid_n = SMOKE_GRID_N if smoke else FULL_GRID_N
    n_converge = SMOKE_N_CONVERGE if smoke else FULL_N_CONVERGE

    work_dir.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    for t_end in t_values:
        sub = work_dir / f"b0_t{str(t_end).replace('.', 'p')}"
        if sub.exists():
            shutil.rmtree(sub)
        raw = _sweep_b0(
            t_end=t_end,
            b0_lo=b0_lo,
            b0_hi=b0_hi,
            b0_n=b0_n,
            grid_n=grid_n,
            parallel=parallel,
            out_dir=sub,
        )
        all_rows.extend(_b0_residuals(raw, t_end))

    # Convergence at a regime point.
    converge_dir = work_dir / "converge"
    if converge_dir.exists():
        shutil.rmtree(converge_dir)
    b0_converge = b0_lo + 0.5 * (b0_hi - b0_lo)
    t_converge = t_values[0]
    converge_raw = _converge(
        n_values=n_converge,
        b0=b0_converge,
        t_end=t_converge,
        parallel=max(1, parallel // 2),
        out_dir=converge_dir,
    )
    converge_rows = _converge_rows(converge_raw, b0=b0_converge, t_end=t_converge)

    # Summary
    rfin = np.array([r["residual_final_vs_bare"] for r in all_rows])
    rcap = np.array([r["residual_max_vs_capped"] for r in all_rows])
    rrs = np.array([r["residual_max_vs_RS"] for r in all_rows])
    summary = {
        "max_abs_residual_final_bare": float(np.max(np.abs(rfin)))
        if rfin.size
        else None,
        "rms_residual_final_bare": float(np.sqrt(np.mean(rfin**2)))
        if rfin.size
        else None,
        "max_abs_residual_max_capped": float(np.max(np.abs(rcap)))
        if rcap.size
        else None,
        "rms_residual_max_capped": float(np.sqrt(np.mean(rcap**2)))
        if rcap.size
        else None,
        "max_abs_residual_max_RS": float(np.max(np.abs(rrs))) if rrs.size else None,
        "rms_residual_max_RS": float(np.sqrt(np.mean(rrs**2))) if rrs.size else None,
        "n_grid_cells": int(rfin.size),
        "n_t_values": len(t_values),
        "min_abs_error_at_highest_N": (
            float(min(r["abs_error_final"] for r in converge_rows))
            if converge_rows
            else None
        ),
    }

    return {
        "metadata": _metadata(
            {
                "b0_range": [b0_lo, b0_hi],
                "b0_n_points": b0_n,
                "t_end_values": t_values,
                "grid_n": grid_n,
                "n_converge_values": n_converge,
                "b0_converge_point": b0_converge,
                "t_converge_point": t_converge,
                "kappa": KAPPA,
                "kwave": KWAVE,
                "bounds": list(BOUNDS),
                "fd_order": 4,
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "results": all_rows,
        "grid": all_rows,
        "convergence": converge_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--parallel", type=int, default=4)
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="boccaletti_uniform_rich_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, parallel=args.parallel, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
