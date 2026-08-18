"""Boccaletti path-integrated reproduction (localized B-field).

For a finite magnetized region the Boccaletti formula generalizes to
the path-integrated form

  P = sin^2( (kappa / 2) * int_{-inf}^{+inf} B_x(z) dz )

and for the Gaussian profile B_x(z) = Bpeak * exp(-z^2 / (2 R^2)),

  P = sin^2( kappa * Bpeak * R * sqrt(pi/2) ).

This benchmark sweeps (Bpeak, R) on a 2D grid using the
gertsenshtein_localized.json spec (derived from theory_localized.toml).
Each cell runs a Gaussian graviton wave-packet on h_7 that traverses
the magnetized region at z=0, then measures the h_7 -> a_2 conversion
probability. Residuals against the path-integrated form are recorded
per cell.

IC convention follows examples/gertsenshtein/run_localized.sh:
  domain [-100, 100], Gaussian IC at z=-50 width 5, wavevector 2.0,
  t_end=120 (one passage through the B-field region).

Note: energy measurement is not reliable on position-dependent
coefficients (known limitation); conversion is the primary metric.

Serves:   manuscript/sections/appendices/validation.tex (App D §2)
Consumes: scripts/figures/figD_boccaletti_localised.py
Writes:   benchmark_results/canonical/boccaletti_localised.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
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
    REPO_ROOT / "benchmark_results" / "canonical" / "boccaletti_localised.json"
)

JSON_SPEC = EXAMPLES / "data" / "gertsenshtein_localized.json"

KAPPA = 1.0
BOUNDS = (-100.0, 100.0)
T_END = 120.0
IC_CENTER = -50.0
IC_WIDTH = 5.0
IC_KWAVE = 2.0
IC_COMPONENT = "h_7"
SOURCE = "h_7"
TARGET = "a_2"

# Full HPC scale: 8 x 5 = 40 cells spanning the perturbative regime
# (kappa*Bpeak*R*sqrt(pi/2) < pi/2 i.e. Bpeak*R < 1.25 at kappa=1).
FULL_BPEAK_VALUES = [0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18]
FULL_R_VALUES = [3.0, 4.5, 6.0, 7.5, 9.0]
FULL_GRID_N = 1024

# Smoke.
SMOKE_BPEAK_VALUES = [0.06, 0.10, 0.14]
SMOKE_R_VALUES = [4.0, 6.0]
SMOKE_GRID_N = 256


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


def _analytic_pathintegrated(bpeak: float, r: float) -> float:
    """P = sin^2(kappa * Bpeak * R * sqrt(pi/2))."""
    return float(math.sin(KAPPA * bpeak * r * math.sqrt(math.pi / 2.0)) ** 2)


def _run_cell(bpeak: float, r: float, *, grid_n: int, out_dir: Path) -> dict:
    sim_cmd = [
        "tidal",
        "simulate",
        str(JSON_SPEC),
        "--grid-shape",
        str(grid_n),
        f"--bounds={BOUNDS[0]}:{BOUNDS[1]}",
        "--no-periodic",
        "--bc",
        "neumann",
        "--ic",
        "gaussian",
        "--ic-wavevector",
        f"{IC_KWAVE}",
        "--ic-amplitude",
        "0.1",
        "--ic-width",
        f"{IC_WIDTH}",
        f"--ic-center={IC_CENTER}",
        "--ic-component",
        IC_COMPONENT,
        "--t-end",
        f"{T_END}",
        "--fd-order",
        "4",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"Bpeak={bpeak}",
        "--param",
        f"R={r}",
        "--output",
        str(out_dir),
        "--force",
    ]
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "Bpeak": bpeak,
            "R": r,
            "ok": False,
            "error": err[-1] if err else "sim failed",
        }

    meas_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "conversion",
        "--source",
        SOURCE,
        "--target",
        TARGET,
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"Bpeak={bpeak}",
        "--param",
        f"R={r}",
        "--json",
        "--quiet",
    ]
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "Bpeak": bpeak,
            "R": r,
            "ok": False,
            "error": err[-1] if err else "meas failed",
        }
    meas = json.loads(res.stdout)
    p_sim = float(meas.get("conversion", {}).get("peak_probability", 0.0))
    p_analytic = _analytic_pathintegrated(bpeak, r)
    return {
        "Bpeak": bpeak,
        "R": r,
        "ok": True,
        "P_sim": p_sim,
        "P_analytic": p_analytic,
        "residual": p_sim - p_analytic,
        "abs_rel_diff": abs(p_sim - p_analytic) / max(abs(p_analytic), 1e-12),
    }


def run(*, smoke: bool, work_dir: Path) -> dict:
    bpeaks = SMOKE_BPEAK_VALUES if smoke else FULL_BPEAK_VALUES
    rs = SMOKE_R_VALUES if smoke else FULL_R_VALUES
    grid_n = SMOKE_GRID_N if smoke else FULL_GRID_N
    work_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for bpeak in bpeaks:
        for r in rs:
            sub = (
                work_dir
                / f"bp{str(bpeak).replace('.', 'p')}_R{str(r).replace('.', 'p')}"
            )
            if sub.exists():
                shutil.rmtree(sub)
            print(f"[boccaletti_localised] Bpeak={bpeak}, R={r}", flush=True)
            rows.append(_run_cell(bpeak, r, grid_n=grid_n, out_dir=sub))

    ok = [r for r in rows if r.get("ok")]
    res = np.array([r["residual"] for r in ok]) if ok else np.array([])
    summary = {
        "n_cells": len(rows),
        "n_ok": len(ok),
        "max_abs_residual": float(np.max(np.abs(res))) if res.size else None,
        "rms_residual": float(np.sqrt(np.mean(res**2))) if res.size else None,
        "max_abs_rel_diff": (max(r["abs_rel_diff"] for r in ok) if ok else None),
    }
    return {
        "metadata": _metadata(
            {
                "Bpeak_values": bpeaks,
                "R_values": rs,
                "grid_n": grid_n,
                "kappa": KAPPA,
                "bounds": list(BOUNDS),
                "t_end": T_END,
                "ic_center": IC_CENTER,
                "ic_width": IC_WIDTH,
                "ic_wavevector": IC_KWAVE,
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "results": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="boccaletti_localised_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
