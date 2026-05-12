r"""Cross-backend self-consistency calibration.

Calibrates the TIDAL solver dispatcher by running the same physical
problem under multiple solver backends (modal, CVODE, IDA, leapfrog) and
recording pairwise relative differences in the conversion probability
$P_\\mathrm{final}$. A consistent dispatcher must produce results that
agree to the regression tolerance of the strictest backend.

Two representative theories are exercised:

  * `coupled_scalars` -- flat, no constraints; covers leapfrog / CVODE / modal
  * `gertsenshtein` -- Einstein--Maxwell baseline; covers modal-eligible flow

For each (theory, backend) pair we record the simulation's `P_final` (or
the peak field amplitude at $t_\\mathrm{end}$ when conversion is not the
natural metric). Pairwise relative differences are computed in the
analysis step.

Serves:   manuscript/sections/appendices/validation.tex (App D, calibration 3)
Consumes: scripts/figures/figD_cross_backend.py
Writes:   benchmark_results/canonical/cross_backend.json
"""

from __future__ import annotations

import argparse
import datetime
import itertools
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "cross_backend.json"

# Backends covered by the dispatcher. `modal` and `leapfrog` need
# periodic BCs; `cvode` and `ida` are general-purpose.
BACKENDS = ["modal", "cvode", "ida", "leapfrog"]

# Theory cases: (label, json_path, ic_component, source, target,
#                 B0_value_or_None, t_end, kwave, grid_n_full, grid_n_smoke)
CASES_FULL = [
    {
        "label": "gertsenshtein",
        "json_spec": "data/gertsenshtein.json",
        "ic_component": "h_5",
        "source": "h_5",
        "target": "a_1",
        "params": {"kappa": 1.0, "B0": 0.05},
        "t_end": 50.0,
        "kwave": 2.0106,
        "bounds": (0.0, 100.0),
        "grid_n_full": 512,
        "grid_n_smoke": 128,
    },
    {
        "label": "coupled_scalars",
        "json_spec": "data/coupled_scalars.json",
        "ic_component": "h_0",
        "source": "h_0",
        "target": "a_0",
        "params": {"omegaP2": 0.0, "mg2": 0.0, "B0": 0.05, "kappa": 1.0},
        "t_end": 10.0,
        "kwave": 1.0,
        "bounds": (0.0, 2.0 * np.pi),
        "grid_n_full": 256,
        "grid_n_smoke": 64,
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


def _run_case(case: dict, backend: str, *, grid_n: int, out_dir: Path) -> dict:
    """Run one (theory, backend) configuration; return P_final via `tidal measure`."""
    bounds = case["bounds"]
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / case["json_spec"]),
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{bounds[0]}:{bounds[1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        f"{case['kwave']}",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        case["ic_component"],
        "--t-end",
        f"{case['t_end']}",
        "--scheme",
        backend,
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in case["params"].items():
        sim_cmd.extend(["--param", f"{k}={v}"])
    print(f"[cross_backend] {case['label']}/{backend}: simulate", flush=True)
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "label": case["label"],
            "backend": backend,
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
        "conversion",
        "--source",
        case["source"],
        "--target",
        case["target"],
        "--json",
        "--quiet",
    ]
    for k, v in case["params"].items():
        meas_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "label": case["label"],
            "backend": backend,
            "ok": False,
            "error": (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "measure failed",
        }
    meas = json.loads(res.stdout)
    conv = meas.get("conversion", {})
    return {
        "label": case["label"],
        "backend": backend,
        "ok": True,
        "P_max": float(conv.get("peak_probability", 0.0)),
        "P_peak_time": float(conv.get("peak_time", 0.0)),
    }


def _pairwise(rows: list[dict]) -> list[dict]:
    """Compute pairwise |Pa - Pb| / max(|Pa|, |Pb|, 1e-12) within each theory."""
    by_label: dict[str, dict[str, float]] = {}
    for r in rows:
        if not r.get("ok"):
            continue
        by_label.setdefault(r["label"], {})[r["backend"]] = r["P_max"]
    out: list[dict] = []
    for label, vals in by_label.items():
        for a, b in itertools.combinations(sorted(vals), 2):
            pa, pb = vals[a], vals[b]
            denom = max(abs(pa), abs(pb), 1e-12)
            out.append(
                {
                    "label": label,
                    "backend_a": a,
                    "backend_b": b,
                    "P_a": pa,
                    "P_b": pb,
                    "rel_diff": abs(pa - pb) / denom,
                }
            )
    return out


def _summary(rows: list[dict], pairs: list[dict]) -> dict:
    ok = [r for r in rows if r.get("ok")]
    failures = [
        (r["label"], r["backend"], r.get("error", "?")) for r in rows if not r.get("ok")
    ]
    rels = [p["rel_diff"] for p in pairs]
    return {
        "n_runs": len(rows),
        "n_ok": len(ok),
        "n_failed": len(rows) - len(ok),
        "failures": failures,
        "max_pairwise_rel_diff": max(rels) if rels else None,
        "median_pairwise_rel_diff": float(np.median(rels)) if rels else None,
    }


def run(*, smoke: bool, work_dir: Path) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for case in CASES_FULL:
        grid_n = case["grid_n_smoke"] if smoke else case["grid_n_full"]
        for backend in BACKENDS:
            sub = work_dir / f"{case['label']}_{backend}"
            if sub.exists():
                shutil.rmtree(sub)
            rows.append(_run_case(case, backend, grid_n=grid_n, out_dir=sub))
    pairs = _pairwise(rows)
    return {
        "metadata": _metadata(
            {
                "backends": BACKENDS,
                "theories": [c["label"] for c in CASES_FULL],
                "smoke": smoke,
            }
        ),
        "summary": _summary(rows, pairs),
        "runs": rows,
        "pairwise": pairs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="cross_backend_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
