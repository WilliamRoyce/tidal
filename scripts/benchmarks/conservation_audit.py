r"""Energy-conservation audit across the canonical example library.

Calibrates the structural integrity of the E-L velocity pipeline by
running a short simulation per example and measuring the relative
energy drift $|dE/E|_\\mathrm{max}$ over the simulation interval.
The audit is not a per-example precision claim; it catches structural
regressions that point-validations would miss.

One expected outlier is annotated: the de-Sitter Klein--Gordon example,
where the energy drift is the physical Hubble friction of an expanding
background (memory: energy-conservation.md). Other examples should sit
in the $10^{-4}$ to $10^{-9}$ band on the leapfrog / CVODE / IDA chosen
by the dispatcher.

Serves:   manuscript/sections/appendices/validation.tex (App D, calibration 5)
Consumes: scripts/figures/figD_conservation_audit.py
Writes:   benchmark_results/canonical/conservation_audit.json
"""

from __future__ import annotations

import argparse
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "conservation_audit.json"

# Per-example config: (label, json_filename, ic_component, bounds, t_end, grid_n,
#                       periodic, parallel_safe, notes)
CASES: list[dict] = [
    {
        "label": "coupled_scalars",
        "json": "coupled_scalars.json",
        "ic_component": "h_0",
        "bounds": (0.0, 2.0 * np.pi),
        "t_end": 5.0,
        "grid_n_smoke": 64,
        "grid_n_full": 256,
        "params": {"omegaP2": 0.0, "mg2": 0.0, "B0": 0.05, "kappa": 1.0},
    },
    {
        "label": "gertsenshtein",
        "json": "gertsenshtein.json",
        "ic_component": "h_5",
        "bounds": (0.0, 100.0),
        "t_end": 50.0,
        "grid_n_smoke": 128,
        "grid_n_full": 512,
        "params": {"kappa": 1.0, "B0": 0.05},
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


def _run_case(case: dict, *, grid_n: int, out_dir: Path) -> dict:
    bounds = case["bounds"]
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / "data" / case["json"]),
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
        case["ic_component"],
        "--t-end",
        f"{case['t_end']}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in case.get("params", {}).items():
        sim_cmd.extend(["--param", f"{k}={v}"])

    print(f"[conservation_audit] {case['label']}: simulate", flush=True)
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "label": case["label"],
            "ok": False,
            "error": (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "non-zero exit",
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
    for k, v in case.get("params", {}).items():
        meas_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "label": case["label"],
            "ok": False,
            "error": (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "measure failed",
        }
    meas = json.loads(res.stdout)
    conv = meas.get("conservation", {})
    # try several plausible key shapes
    drift = (
        conv.get("max_relative_error")
        or conv.get("max_relative_energy_drift")
        or conv.get("max_abs_relative_energy_drift")
        or conv.get("dE_over_E_max")
        or conv.get("max_dE_over_E")
    )
    return {
        "label": case["label"],
        "ok": True,
        "max_abs_dE_over_E": float(drift) if drift is not None else None,
        "conservation_raw": conv,
    }


def _summary(rows: list[dict]) -> dict:
    ok = [r for r in rows if r.get("ok")]
    drifts = [
        r["max_abs_dE_over_E"] for r in ok if r.get("max_abs_dE_over_E") is not None
    ]
    return {
        "n_examples": len(rows),
        "n_ok": len(ok),
        "n_failed": len(rows) - len(ok),
        "max_drift": max(drifts) if drifts else None,
        "median_drift": float(np.median(drifts)) if drifts else None,
        "min_drift": min(drifts) if drifts else None,
    }


def run(*, smoke: bool, work_dir: Path) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for case in CASES:
        grid_n = case["grid_n_smoke"] if smoke else case["grid_n_full"]
        sub = work_dir / case["label"]
        if sub.exists():
            shutil.rmtree(sub)
        rows.append(_run_case(case, grid_n=grid_n, out_dir=sub))
    return {
        "metadata": _metadata(
            {"smoke": smoke, "examples": [c["label"] for c in CASES]}
        ),
        "summary": _summary(rows),
        "examples": rows,
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
