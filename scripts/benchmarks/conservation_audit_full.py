"""Energy-conservation audit across the example library.

Runs `tidal simulate` + `tidal measure --what conservation` on a curated
list of canonical examples and records the maximum relative energy
drift |dE/E|_max per example. The audit catches pipeline regressions
that point-validations would miss: a constraint solver that returns
inconsistent initial conditions, a symbolic stage that loses
cross-coupling terms in field decomposition, or a time integrator
that drifts past its predicted order.

The selection here mirrors the documented energy-conservation table
in memory/energy-conservation.md (2026-02-27 baseline). Two outliers
are tracked separately:

  - de_sitter_kg: 30.9% drift is physical (Hubble friction)
  - massive_3form: Hamiltonian skipped (rank-3 Lagrangian decomposition)

Serves:   manuscript/sections/appendices/validation.tex (App D §8)
Consumes: scripts/figures/figD_conservation.py
Writes:   benchmark_results/canonical/conservation_audit_full.json
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
    REPO_ROOT / "benchmark_results" / "canonical" / "conservation_audit_full.json"
)

# Each entry: (label, json filename, ic_component, bounds, t_end,
#  grid_n_full, grid_n_smoke, params, notes).
# Grid sizes deliberately modest to keep the audit fast — the metric
# is structural integrity, not per-example precision.
EXAMPLES_LIST: list[dict] = [
    {
        "label": "gertsenshtein",
        "json": "gertsenshtein.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 50.0,
        "grid_n_full": 256,
        "grid_n_smoke": 64,
        "params": {"kappa": 1.0, "B0": 0.05},
    },
    {
        "label": "coupled_scalars",
        "json": "coupled_scalars.json",
        "ic_component": "h_0",
        "bounds": (0.0, math.tau),
        "t_end": 10.0,
        "grid_n_full": 128,
        "grid_n_smoke": 32,
        "params": {"kappa": 1.0, "B0": 0.05, "omegaP2": 0.0, "mg2": 0.0},
    },
    {
        "label": "gertsenshtein_proca",
        "json": "gertsenshtein_proca.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 50.0,
        "grid_n_full": 256,
        "grid_n_smoke": 64,
        "params": {"kappa": 1.0, "B0": 0.05, "mA2": 1.0},
    },
    {
        "label": "coupled_scattering",
        "json": "coupled_scattering.json",
        "ic_component": "chi_0",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 256,
        "grid_n_smoke": 64,
        "params": {},
    },
    {
        "label": "chern_simons_3d",
        "json": "chern_simons_3d.json",
        "ic_component": "a_0",
        "bounds": (0.0, 10.0),
        "t_end": 5.0,
        "grid_n_full": 32,
        "grid_n_smoke": 16,
        "params": {},
    },
    {
        "label": "conformal_kg_static",
        "json": "conformal_kg_static.json",
        "ic_component": "phi",
        "bounds": (0.0, 10.0),
        "t_end": 5.0,
        "grid_n_full": 128,
        "grid_n_smoke": 32,
        "params": {},
    },
    {
        "label": "cylindrical_kg_1d",
        "json": "cylindrical_kg_1d.json",
        "ic_component": "phi",
        "bounds": (0.5, 10.0),
        "t_end": 5.0,
        "grid_n_full": 128,
        "grid_n_smoke": 32,
        "params": {},
    },
    {
        "label": "navier_cauchy_2d",
        "json": "navier_cauchy_2d.json",
        "ic_component": "u_0",
        "bounds": (0.0, 10.0),
        "t_end": 3.0,
        "grid_n_full": 32,
        "grid_n_smoke": 16,
        "params": {},
    },
    {
        "label": "gw_plane_wave_1d",
        "json": "gw_plane_wave_1d.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 256,
        "grid_n_smoke": 64,
        "params": {"kappa": 1.0},
    },
    {
        "label": "torsion_gertsenshtein",
        "json": "torsion_gertsenshtein.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 128,
        "grid_n_smoke": 64,
        "params": {
            "kappa": 1.0,
            "B0": 0.05,
            "alpha1": 0.0,
            "alpha2": 0.0,
            "alpha3": 0.0,
            "b5": 0.0,
        },
    },
    {
        "label": "torsion_gertsenshtein_nonminimal",
        "json": "torsion_gertsenshtein_nonminimal.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 128,
        "grid_n_smoke": 64,
        "params": {
            "kappa": 1.0,
            "B0": 0.05,
            "alpha1": 0.0,
            "alpha2": 0.0,
            "alpha3": 0.0,
            "delta1": 0.0,
        },
    },
    {
        "label": "dark_photon_plasma",
        "json": "dark_photon_plasma.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 64,
        "grid_n_smoke": 32,
        "params": {
            "kappa": 1.0,
            "B0": 0.05,
            "alpha3": 0.123,
            "xi": 0.274,
            "deltam": 0.01,
            "mA2": 0.955,
        },
    },
    {
        "label": "torsion_dark_photon_fv",
        "json": "torsion_dark_photon_fv.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 30.0,
        "grid_n_full": 128,
        "grid_n_smoke": 64,
        "params": {
            "kappa": 1.0,
            "B0": 0.05,
            "mT2": 0.246,
            "xi": 0.274,
            "deltam": 0.01,
            "mA2": 0.955,
        },
    },
]


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


def _run_example(ex: dict, *, grid_n: int, out_dir: Path) -> dict:
    bounds = ex["bounds"]
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / "data" / ex["json"]),
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{bounds[0]}:{bounds[1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        ex["ic_component"],
        "--t-end",
        f"{ex['t_end']}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in ex.get("params", {}).items():
        sim_cmd.extend(["--param", f"{k}={v}"])

    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "label": ex["label"],
            "ok": False,
            "stage": "simulate",
            "error": err[-1] if err else "non-zero exit",
        }

    meas_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "conservation",
        "--json",
        "--quiet",
    ]
    for k, v in ex.get("params", {}).items():
        meas_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "label": ex["label"],
            "ok": False,
            "stage": "measure",
            "error": err[-1] if err else "non-zero exit",
        }
    meas = json.loads(res.stdout)
    conv = meas.get("conservation", {})
    drift = (
        conv.get("max_relative_error")
        or conv.get("max_relative_energy_drift")
        or conv.get("dE_over_E_max")
    )
    return {
        "label": ex["label"],
        "ok": True,
        "max_abs_dE_over_E": float(drift) if drift is not None else None,
        "is_conserved": bool(conv.get("is_conserved", False)),
        "threshold": float(conv.get("threshold", 1e-3)),
    }


def run(*, smoke: bool, work_dir: Path) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    examples = EXAMPLES_LIST[:4] if smoke else EXAMPLES_LIST
    for ex in examples:
        grid_n = ex["grid_n_smoke"] if smoke else ex["grid_n_full"]
        sub = work_dir / ex["label"]
        if sub.exists():
            shutil.rmtree(sub)
        rows.append(_run_example(ex, grid_n=grid_n, out_dir=sub))

    ok = [r for r in rows if r.get("ok")]
    drifts = [
        r["max_abs_dE_over_E"] for r in ok if r.get("max_abs_dE_over_E") is not None
    ]
    summary = {
        "n_examples": len(rows),
        "n_ok": len(ok),
        "n_failed": len(rows) - len(ok),
        "max_drift": max(drifts) if drifts else None,
        "median_drift": float(np.median(drifts)) if drifts else None,
        "min_drift": min(drifts) if drifts else None,
    }
    return {
        "metadata": _metadata(
            {
                "smoke": smoke,
                "examples": [ex["label"] for ex in examples],
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
    with tempfile.TemporaryDirectory(prefix="conservation_audit_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
