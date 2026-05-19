"""Multi-method agreement: 8 backends × 2 theories + Hamiltonian-vs-EOM.

Broader than the old cross_backend benchmark. Records pairwise relative
agreement in P_max for every applicable solver backend on each of two
representative theories, plus a Hamiltonian-derived vs EOM-derived
energy cross-check on one canonical theory.

Backends covered: modal, cvode, ida, leapfrog (Yoshida 2 and 4 via
--leapfrog-order), scipy DOP853 / Radau / BDF (via --scheme scipy
--method NAME).

Theories:
  - gertsenshtein           — Einstein-Maxwell baseline (modal-eligible)
  - coupled_scalars         — Raffelt-Stodolsky-style two-channel mixing
                              (cite raffelt1988mixing in App D §6)

The Hamiltonian-vs-EOM measurement-path cross-check runs `tidal measure
--what energy` (Wolfram Hamiltonian, Parseval gradient energy on periodic
domains; memory: spectral_energy_parseval.md) and `--what conservation`
(EOM-based energy drift), recording the agreement between the two
independent code paths.

Serves:   manuscript/sections/appendices/validation.tex (App D §6)
Consumes: scripts/figures/figD_multi_method.py
Writes:   benchmark_results/canonical/multi_method.json
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "multi_method.json"

# Each backend entry: (label, scheme_args). scheme_args is a list of
# extra argv tokens for `tidal simulate`. Modal is auto-only -- pass it
# explicitly via --scheme modal.
BACKENDS: list[tuple[str, list[str]]] = [
    ("modal", ["--scheme", "modal"]),
    ("cvode", ["--scheme", "cvode"]),
    ("ida", ["--scheme", "ida"]),
    ("leapfrog_Y2", ["--scheme", "leapfrog", "--leapfrog-order", "2"]),
    ("leapfrog_Y4", ["--scheme", "leapfrog", "--leapfrog-order", "4"]),
    ("scipy_DOP853", ["--scheme", "scipy", "--method", "DOP853"]),
    # scipy Radau and BDF are implicit and hang convergence-iterating
    # on these wave problems at our grid sizes (~40 min per cell with
    # no return). Excluded; the implicit family is already represented
    # by IDA and CVODE.
]
# Per-cell wall-clock timeout — safety net against future backend
# regressions silently locking up the benchmark.
CELL_TIMEOUT = 240

THEORIES: list[dict] = [
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
        "grid_n_full": 256,
        "grid_n_smoke": 64,
    },
    {
        "label": "coupled_scalars",
        "json_spec": "data/coupled_scalars.json",
        "ic_component": "h_0",
        "source": "h_0",
        "target": "a_0",
        "params": {"kappa": 1.0, "B0": 0.05, "omegaP2": 0.0, "mg2": 0.0},
        "t_end": 10.0,
        "kwave": 1.0,
        "bounds": (0.0, 2.0 * np.pi),
        "grid_n_full": 128,
        "grid_n_smoke": 32,
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


def _run_one(
    theory: dict, backend: tuple[str, list[str]], *, grid_n: int, out_dir: Path
) -> dict:
    label, scheme_args = backend
    bounds = theory["bounds"]
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / theory["json_spec"]),
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{bounds[0]}:{bounds[1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        f"{theory['kwave']}",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        theory["ic_component"],
        "--t-end",
        f"{theory['t_end']}",
        "--fd-order",
        "4",
        *scheme_args,
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in theory["params"].items():
        sim_cmd.extend(["--param", f"{k}={v}"])

    try:
        res = subprocess.run(
            sim_cmd,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
            check=False,
            timeout=CELL_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return {
            "theory": theory["label"],
            "backend": label,
            "ok": False,
            "error": f"sim timed out after {CELL_TIMEOUT}s",
        }
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "theory": theory["label"],
            "backend": label,
            "ok": False,
            "error": err[-1] if err else "non-zero exit",
        }

    meas_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "conversion",
        "--source",
        theory["source"],
        "--target",
        theory["target"],
        "--json",
        "--quiet",
    ]
    for k, v in theory["params"].items():
        meas_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        err = (res.stderr or res.stdout or "").splitlines()
        return {
            "theory": theory["label"],
            "backend": label,
            "ok": False,
            "error": err[-1] if err else "meas failed",
        }
    meas = json.loads(res.stdout)
    p = float(meas.get("conversion", {}).get("peak_probability", 0.0))
    return {"theory": theory["label"], "backend": label, "ok": True, "P_max": p}


def _measurement_path_cross_check(theory: dict, *, grid_n: int, out_dir: Path) -> dict:
    """Hamiltonian-vs-EOM energy agreement on one theory under modal."""
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / theory["json_spec"]),
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{theory['bounds'][0]}:{theory['bounds'][1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        f"{theory['kwave']}",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        theory["ic_component"],
        "--t-end",
        f"{theory['t_end']}",
        "--fd-order",
        "4",
        "--scheme",
        "modal",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in theory["params"].items():
        sim_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {"theory": theory["label"], "ok": False, "error": "sim failed"}

    # Energy measurement (Hamiltonian-derived).
    e_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "energy",
        "--json",
        "--quiet",
    ]
    for k, v in theory["params"].items():
        e_cmd.extend(["--param", f"{k}={v}"])
    e_res = subprocess.run(
        e_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    energy_meas = json.loads(e_res.stdout) if e_res.returncode == 0 else {}

    # Conservation diagnostic (EOM-based).
    c_cmd = [
        "tidal",
        "measure",
        str(out_dir),
        "--what",
        "conservation",
        "--json",
        "--quiet",
    ]
    for k, v in theory["params"].items():
        c_cmd.extend(["--param", f"{k}={v}"])
    c_res = subprocess.run(
        c_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    cons_meas = json.loads(c_res.stdout) if c_res.returncode == 0 else {}

    return {
        "theory": theory["label"],
        "ok": e_res.returncode == 0 and c_res.returncode == 0,
        "energy_payload": energy_meas.get("energy", {}),
        "conservation_payload": cons_meas.get("conservation", {}),
    }


def _pairwise(rows: list[dict]) -> list[dict]:
    by_theory: dict[str, dict[str, float]] = {}
    for r in rows:
        if not r.get("ok"):
            continue
        by_theory.setdefault(r["theory"], {})[r["backend"]] = r["P_max"]
    out: list[dict] = []
    for theory, vals in by_theory.items():
        for a, b in itertools.combinations(sorted(vals), 2):
            pa, pb = vals[a], vals[b]
            denom = max(abs(pa), abs(pb), 1e-30)
            out.append(
                {
                    "theory": theory,
                    "backend_a": a,
                    "backend_b": b,
                    "P_a": pa,
                    "P_b": pb,
                    "rel_diff": abs(pa - pb) / denom,
                }
            )
    return out


def run(*, smoke: bool, work_dir: Path) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    backends = BACKENDS[:4] if smoke else BACKENDS
    rows: list[dict] = []
    for theory in THEORIES:
        grid_n = theory["grid_n_smoke"] if smoke else theory["grid_n_full"]
        for backend in backends:
            sub = work_dir / f"{theory['label']}_{backend[0]}"
            if sub.exists():
                shutil.rmtree(sub)
            rows.append(_run_one(theory, backend, grid_n=grid_n, out_dir=sub))

    pairs = _pairwise(rows)

    # Hamiltonian-vs-EOM cross-check on gertsenshtein only.
    ham_dir = work_dir / "ham_vs_eom"
    if ham_dir.exists():
        shutil.rmtree(ham_dir)
    grid_n = THEORIES[0]["grid_n_smoke"] if smoke else THEORIES[0]["grid_n_full"]
    ham_eom = _measurement_path_cross_check(THEORIES[0], grid_n=grid_n, out_dir=ham_dir)

    ok = [r for r in rows if r.get("ok")]
    failures = [
        (r["theory"], r["backend"], r.get("error", "?"))
        for r in rows
        if not r.get("ok")
    ]
    rels = [p["rel_diff"] for p in pairs]
    summary = {
        "n_runs": len(rows),
        "n_ok": len(ok),
        "n_failed": len(rows) - len(ok),
        "failures": failures,
        "max_pairwise_rel_diff": max(rels) if rels else None,
        "median_pairwise_rel_diff": float(np.median(rels)) if rels else None,
        "backends_per_theory": len(backends),
        "theories": [t["label"] for t in THEORIES],
    }
    return {
        "metadata": _metadata(
            {
                "backends": [b[0] for b in backends],
                "theories": [t["label"] for t in THEORIES],
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "results": rows,
        "pairwise": pairs,
        "hamiltonian_vs_eom": ham_eom,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="multi_method_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
