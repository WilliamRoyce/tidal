r"""Time-integration order calibration on the Gertsenshtein baseline.

Calibrates the time-integration schemes used by TIDAL by exercising
leapfrog (Yoshida-2 and Yoshida-4) and CVODE on a problem with a known
closed-form $P_\\mathrm{final}$ (the Boccaletti kernel at low $B_0$ in
the perturbative regime). We sweep the time step $\\Delta t$ and record
the L1 error $|P_\\mathrm{final}^\\mathrm{sim}(\\Delta t)
                  - P_\\mathrm{final}^\\mathrm{analytic}|$.
For a $p$-th order scheme the error should scale as $\\Delta t^p$ until
the spatial discretisation error or round-off floor dominates.

Schemes calibrated:
  * leapfrog --leapfrog-order 2 (Yoshida-2, expect O(dt^2))
  * leapfrog --leapfrog-order 4 (Yoshida-4, expect O(dt^4))
  * cvode (adaptive; expect monotonic improvement under tightening rtol)

Serves:   manuscript/sections/appendices/validation.tex (App D, calibration 4)
Consumes: scripts/figures/figD_dispersion_time_order.py
Writes:   benchmark_results/canonical/time_integration_order.json
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
    REPO_ROOT / "benchmark_results" / "canonical" / "time_integration_order.json"
)

# Physics: stay in the linear regime (low B0), short integration.
KAPPA = 1.0
B0_FIXED = 0.05
KWAVE = 2.0106
BOUNDS = (0.0, 100.0)
T_END = 50.0
SOURCE = "h_5"
TARGET = "a_1"
IC_COMPONENT = "h_5"

# Reference / analytic: bare sin² in the perturbative regime.
# kappa*B0*t/2 = 0.5*0.05*50 = 1.25 < pi/2; first quarter-period not reached.
ANALYTIC_PFINAL = math.sin(0.5 * KAPPA * B0_FIXED * T_END) ** 2

# Time-step ladder. Smaller dt → higher CFL safety; the asymptote shape
# is the diagnostic, not any single error value.
FULL_DT_LADDER = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]
SMOKE_DT_LADDER = [0.1, 0.05, 0.025]
GRID_N_FULL = 512
GRID_N_SMOKE = 128


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
    *, scheme: str, dt: float, grid_n: int, leapfrog_order: int | None, out_dir: Path
) -> dict:
    sim_cmd = [
        "tidal",
        "simulate",
        str(EXAMPLES / "data" / "gertsenshtein.json"),
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
        f"{T_END}",
        "--dt",
        f"{dt}",
        "--scheme",
        scheme,
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
    if leapfrog_order is not None:
        sim_cmd.extend(["--leapfrog-order", str(leapfrog_order)])
    label = f"{scheme}{leapfrog_order or ''}@dt={dt}"
    print(f"[time_integration_order] {label}: simulate", flush=True)
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "scheme": scheme,
            "leapfrog_order": leapfrog_order,
            "dt": dt,
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
        SOURCE,
        "--target",
        TARGET,
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_FIXED}",
        "--json",
        "--quiet",
    ]
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {
            "scheme": scheme,
            "leapfrog_order": leapfrog_order,
            "dt": dt,
            "ok": False,
            "error": (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "measure failed",
        }
    meas = json.loads(res.stdout)
    conv = meas.get("conversion", {})
    p_final = float(conv.get("peak_probability", 0.0))
    return {
        "scheme": scheme,
        "leapfrog_order": leapfrog_order,
        "dt": dt,
        "ok": True,
        "P_final_sim": p_final,
        "P_final_analytic": ANALYTIC_PFINAL,
        "abs_error": abs(p_final - ANALYTIC_PFINAL),
    }


def _fit_order(rows: list[dict]) -> float | None:
    """Log-log slope of (abs_error vs dt) on the cleanly-converging segment."""
    ok = [r for r in rows if r.get("ok") and r["abs_error"] > 1e-14]
    if len(ok) < 2:
        return None
    dts = np.log10([r["dt"] for r in ok])
    errs = np.log10([r["abs_error"] for r in ok])
    slope, _ = np.polyfit(dts, errs, 1)
    return float(slope)


def run(*, smoke: bool, work_dir: Path) -> dict:
    dt_ladder = SMOKE_DT_LADDER if smoke else FULL_DT_LADDER
    grid_n = GRID_N_SMOKE if smoke else GRID_N_FULL
    work_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    # leapfrog Yoshida-2 and Yoshida-4
    for lf_order in (2, 4):
        for dt in dt_ladder:
            sub = work_dir / f"lf{lf_order}_dt{str(dt).replace('.', 'p')}"
            if sub.exists():
                shutil.rmtree(sub)
            rows.append(
                _run_one(
                    scheme="leapfrog",
                    dt=dt,
                    grid_n=grid_n,
                    leapfrog_order=lf_order,
                    out_dir=sub,
                )
            )
    # cvode adaptive (dt here is the max step)
    for dt in dt_ladder:
        sub = work_dir / f"cvode_dt{str(dt).replace('.', 'p')}"
        if sub.exists():
            shutil.rmtree(sub)
        rows.append(
            _run_one(
                scheme="cvode", dt=dt, grid_n=grid_n, leapfrog_order=None, out_dir=sub
            )
        )

    # Slope-fit per scheme/order
    slopes = {
        "leapfrog_2": _fit_order(
            [
                r
                for r in rows
                if r.get("scheme") == "leapfrog" and r.get("leapfrog_order") == 2
            ]
        ),
        "leapfrog_4": _fit_order(
            [
                r
                for r in rows
                if r.get("scheme") == "leapfrog" and r.get("leapfrog_order") == 4
            ]
        ),
        "cvode": _fit_order([r for r in rows if r.get("scheme") == "cvode"]),
    }
    summary = {
        "n_runs": len(rows),
        "n_ok": sum(1 for r in rows if r.get("ok")),
        "fitted_slopes": slopes,
        "expected_slopes": {"leapfrog_2": 2.0, "leapfrog_4": 4.0, "cvode": "adaptive"},
    }
    return {
        "metadata": _metadata(
            {
                "dt_ladder": dt_ladder,
                "grid_n": grid_n,
                "kappa": KAPPA,
                "B0": B0_FIXED,
                "kwave": KWAVE,
                "bounds": list(BOUNDS),
                "t_end": T_END,
                "analytic_P_final": ANALYTIC_PFINAL,
                "fd_order": 4,
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "runs": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="time_integration_order_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
