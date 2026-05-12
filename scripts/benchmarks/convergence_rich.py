"""Rich convergence diagnostics: dispersion + Rabi-N sweep + multi-scheme dt.

Combines three convergence calibrations into one HPC-parallel campaign:

  (a) Proca-mass dispersion ω(k) = sqrt(k² + m²) at three values of m²
      (0.1, 1.0, 10.0), recovered by spatial-FFT mode identification +
      time-FFT with Hann window and quadratic peak refinement.

  (b) Rabi-frequency grid convergence: Ω_eff / Ω_theory at six grid
      resolutions N ∈ {128, 256, 512, 1024, 2048, 4096} on the
      Einstein-Maxwell Gertsenshtein baseline. Expected (kΔx)² scaling
      of the FD-4 stencil error.

  (c) Time-integration order: leapfrog Y2/Y4, CVODE, IDA, scipy
      DOP853/Radau/BDF on a low-B0 Einstein-Maxwell run; recovers the
      predicted algebraic order from a log-log fit of |P_final - sin²|
      vs Δt.

Serves:   manuscript/sections/appendices/validation.tex (App D §7)
Consumes: scripts/figures/figD_convergence.py
Writes:   benchmark_results/canonical/convergence_rich.json
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "convergence_rich.json"

KAPPA = 1.0
B0_RABI = 0.05
KWAVE = 2.0106
BOUNDS = (0.0, 100.0)
T_END = 50.0

# (a) Proca dispersion params
PROCA_SPEC = EXAMPLES / "data" / "gertsenshtein_proca.json"
PROCA_MASS_SQ = [0.1, 1.0, 10.0]
PROCA_KMULT_FULL = [4, 8, 16, 24, 32, 48]
PROCA_KMULT_SMOKE = [4, 8, 16]
PROCA_T_END_FULL = 100.0
PROCA_T_END_SMOKE = 30.0
PROCA_GRID_N_FULL = 1024
PROCA_GRID_N_SMOKE = 256
PROCA_SNAP_DT_FULL = 0.05
PROCA_SNAP_DT_SMOKE = 0.05

# (b) Rabi N-sweep params
RABI_N_FULL = [128, 256, 512, 1024, 2048, 4096]
RABI_N_SMOKE = [128, 256]

# (c) Time-integration order params
TIO_SCHEMES_FULL: list[tuple[str, list[str]]] = [
    ("leapfrog_Y2", ["--scheme", "leapfrog", "--leapfrog-order", "2"]),
    ("leapfrog_Y4", ["--scheme", "leapfrog", "--leapfrog-order", "4"]),
    ("cvode", ["--scheme", "cvode"]),
    ("ida", ["--scheme", "ida"]),
    ("scipy_DOP853", ["--scheme", "scipy", "--method", "DOP853"]),
    # scipy Radau and BDF: implicit solvers hang convergence-iterating
    # on these wave problems at our grid sizes. Excluded; implicit
    # behaviour is already represented by IDA and CVODE.
]
TIO_SCHEMES_SMOKE = TIO_SCHEMES_FULL[:3]
TIO_CELL_TIMEOUT = 240  # seconds; safety net per (scheme, dt) cell
TIO_DT_FULL = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625]
TIO_DT_SMOKE = [0.1, 0.05, 0.025]
TIO_GRID_N_FULL = 512
TIO_GRID_N_SMOKE = 128
TIO_ANALYTIC = math.sin(0.5 * KAPPA * B0_RABI * T_END) ** 2


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


# -------- (a) Proca dispersion --------


def _proca_simulate(
    *,
    k_request: float,
    mass2: float,
    grid_n: int,
    snapshot_dt: float,
    t_end: float,
    out_dir: Path,
) -> None:
    cmd = [
        "tidal",
        "simulate",
        str(PROCA_SPEC),
        "--grid-shape",
        str(grid_n),
        "--bounds",
        f"{BOUNDS[0]}:{BOUNDS[1]}",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        f"{k_request}",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        "a_1",
        "--t-end",
        f"{t_end}",
        "--snapshots",
        f"{snapshot_dt}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        "B0=0.0",
        "--param",
        f"mA2={mass2}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    print(f"[convergence_rich] proca: m2={mass2}, k_req={k_request:.4f}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def _proca_omega(out_dir: Path, component: str = "a_1") -> tuple[float, float]:
    times = np.load(out_dir / "times.npy")
    field = np.load(out_dir / f"{component}.npy")
    nx = field.shape[1]
    x_axis_len = BOUNDS[1] - BOUNDS[0]
    spec0 = np.fft.rfft(field[0])
    kx_axis = np.fft.rfftfreq(nx, d=x_axis_len / nx) * 2.0 * math.pi
    mode_idx = int(np.argmax(np.abs(spec0)))
    k_real = float(kx_axis[mode_idx])
    mode_amp = np.array(
        [np.fft.rfft(field[i])[mode_idx] for i in range(field.shape[0])]
    )
    series = mode_amp.real - mode_amp.real.mean()
    window = np.hanning(series.size)
    windowed = series * window
    dt = float(times[1] - times[0]) if times.size >= 2 else 1.0
    spec = np.fft.rfft(windowed)
    spec[0] = 0
    mag = np.abs(spec)
    freqs = np.fft.rfftfreq(series.size, d=dt)
    peak_idx = int(np.argmax(mag))
    f_peak = float(freqs[peak_idx])
    if 1 <= peak_idx < mag.size - 1:
        y0, y1, y2 = mag[peak_idx - 1], mag[peak_idx], mag[peak_idx + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            delta = 0.5 * (y0 - y2) / denom
            df = float(freqs[1] - freqs[0])
            f_peak = float(freqs[peak_idx]) + delta * df
    omega = 2.0 * math.pi * f_peak
    return omega, k_real


def _run_dispersion(*, smoke: bool, work_dir: Path) -> list[dict]:
    mass_values = PROCA_MASS_SQ
    k_mults = PROCA_KMULT_SMOKE if smoke else PROCA_KMULT_FULL
    grid_n = PROCA_GRID_N_SMOKE if smoke else PROCA_GRID_N_FULL
    snapshot_dt = PROCA_SNAP_DT_SMOKE if smoke else PROCA_SNAP_DT_FULL
    t_end = PROCA_T_END_SMOKE if smoke else PROCA_T_END_FULL
    base_k = 2.0 * math.pi / (BOUNDS[1] - BOUNDS[0])
    rows: list[dict] = []
    for m2 in mass_values:
        for kmult in k_mults:
            k_req = base_k * kmult
            sub = work_dir / f"proca_m{str(m2).replace('.', 'p')}_k{kmult}"
            if sub.exists():
                shutil.rmtree(sub)
            try:
                _proca_simulate(
                    k_request=k_req,
                    mass2=m2,
                    grid_n=grid_n,
                    snapshot_dt=snapshot_dt,
                    t_end=t_end,
                    out_dir=sub,
                )
                omega_sim, k_real = _proca_omega(sub)
            except (
                subprocess.CalledProcessError,
                FileNotFoundError,
                ValueError,
            ) as exc:
                rows.append(
                    {"mass2": m2, "k_requested": k_req, "ok": False, "error": str(exc)}
                )
                continue
            omega_ana = math.sqrt(k_real**2 + m2)
            rows.append(
                {
                    "mass2": m2,
                    "k_requested": k_req,
                    "k_realised": k_real,
                    "ok": True,
                    "omega_sim": omega_sim,
                    "omega_analytic": omega_ana,
                    "rel_error": abs(omega_sim - omega_ana)
                    / max(abs(omega_ana), 1e-30),
                }
            )
    return rows


# -------- (b) Rabi-frequency grid convergence --------


def _rabi_simulate(*, grid_n: int, out_dir: Path, snapshot_dt: float) -> None:
    cmd = [
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
        "h_5",
        "--t-end",
        f"{T_END}",
        "--snapshots",
        f"{snapshot_dt}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_RABI}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    print(f"[convergence_rich] rabi: N={grid_n}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def _rabi_omega_eff(out_dir: Path) -> float:
    times = np.load(out_dir / "times.npy")
    h5 = np.load(out_dir / "h_5.npy")
    # global amplitude time series via max-abs across space
    series = np.array([np.max(np.abs(h5[i])) for i in range(h5.shape[0])])
    series -= series.mean()
    window = np.hanning(series.size)
    dt = float(times[1] - times[0]) if times.size >= 2 else 1.0
    spec = np.abs(np.fft.rfft(series * window))
    spec[0] = 0
    freqs = np.fft.rfftfreq(series.size, d=dt)
    peak_idx = int(np.argmax(spec))
    f_peak = float(freqs[peak_idx])
    if 1 <= peak_idx < spec.size - 1:
        y0, y1, y2 = spec[peak_idx - 1], spec[peak_idx], spec[peak_idx + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            delta = 0.5 * (y0 - y2) / denom
            df = float(freqs[1] - freqs[0])
            f_peak = float(freqs[peak_idx]) + delta * df
    return 2.0 * math.pi * f_peak


def _run_rabi(*, smoke: bool, work_dir: Path) -> list[dict]:
    n_values = RABI_N_SMOKE if smoke else RABI_N_FULL
    omega_theory = KAPPA * B0_RABI  # Rabi angular frequency for sin^2(omega t/2)
    snapshot_dt = 0.1
    rows: list[dict] = []
    for n in n_values:
        sub = work_dir / f"rabi_N{n}"
        if sub.exists():
            shutil.rmtree(sub)
        try:
            _rabi_simulate(grid_n=n, out_dir=sub, snapshot_dt=snapshot_dt)
            omega_eff = _rabi_omega_eff(sub)
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as exc:
            rows.append({"N": n, "ok": False, "error": str(exc)})
            continue
        rows.append(
            {
                "N": n,
                "ok": True,
                "omega_eff": omega_eff,
                "omega_theory": omega_theory,
                "ratio": omega_eff / max(omega_theory, 1e-30),
                "k_dx": KWAVE * (BOUNDS[1] - BOUNDS[0]) / n,
            }
        )
    return rows


# -------- (c) Time-integration order --------


def _tio_run_one(
    scheme_label: str, scheme_args: list[str], dt: float, *, grid_n: int, out_dir: Path
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
        "h_5",
        "--t-end",
        f"{T_END}",
        "--dt",
        f"{dt}",
        *scheme_args,
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0_RABI}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    try:
        res = subprocess.run(
            sim_cmd,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
            check=False,
            timeout=TIO_CELL_TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return {
            "scheme": scheme_label,
            "dt": dt,
            "ok": False,
            "error": f"sim timed out after {TIO_CELL_TIMEOUT}s",
        }
    if res.returncode != 0:
        return {"scheme": scheme_label, "dt": dt, "ok": False, "error": "sim failed"}

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
        f"B0={B0_RABI}",
        "--json",
        "--quiet",
    ]
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return {"scheme": scheme_label, "dt": dt, "ok": False, "error": "meas failed"}
    meas = json.loads(res.stdout)
    p = float(meas.get("conversion", {}).get("peak_probability", 0.0))
    return {
        "scheme": scheme_label,
        "dt": dt,
        "ok": True,
        "P_final_sim": p,
        "P_final_analytic": TIO_ANALYTIC,
        "abs_error": abs(p - TIO_ANALYTIC),
    }


def _run_tio(*, smoke: bool, work_dir: Path) -> list[dict]:
    schemes = TIO_SCHEMES_SMOKE if smoke else TIO_SCHEMES_FULL
    dts = TIO_DT_SMOKE if smoke else TIO_DT_FULL
    grid_n = TIO_GRID_N_SMOKE if smoke else TIO_GRID_N_FULL
    rows: list[dict] = []
    for label, args in schemes:
        for dt in dts:
            sub = work_dir / f"tio_{label}_dt{str(dt).replace('.', 'p')}"
            if sub.exists():
                shutil.rmtree(sub)
            rows.append(_tio_run_one(label, args, dt, grid_n=grid_n, out_dir=sub))
    return rows


def _fit_slope(rows: list[dict], scheme: str) -> float | None:
    ok = [
        r
        for r in rows
        if r.get("ok")
        and r.get("scheme") == scheme
        and r.get("abs_error") is not None
        and math.isfinite(r["abs_error"])
        and r["abs_error"] > 1e-14
    ]
    if len(ok) < 2:
        return None
    dts = np.log10([r["dt"] for r in ok])
    errs = np.log10([r["abs_error"] for r in ok])
    slope, _ = np.polyfit(dts, errs, 1)
    return float(slope)


def run(*, smoke: bool, work_dir: Path) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    dispersion_rows = _run_dispersion(smoke=smoke, work_dir=work_dir)
    rabi_rows = _run_rabi(smoke=smoke, work_dir=work_dir)
    tio_rows = _run_tio(smoke=smoke, work_dir=work_dir)

    # Slope fits
    schemes_done = sorted({r["scheme"] for r in tio_rows if r.get("ok")})
    slopes = {s: _fit_slope(tio_rows, s) for s in schemes_done}

    dispersion_ok = [r for r in dispersion_rows if r.get("ok")]
    rels = [r["rel_error"] for r in dispersion_ok]
    summary = {
        "dispersion_n_runs": len(dispersion_rows),
        "dispersion_n_ok": len(dispersion_ok),
        "dispersion_max_rel_error": max(rels) if rels else None,
        "dispersion_median_rel_error": float(np.median(rels)) if rels else None,
        "rabi_n_runs": len(rabi_rows),
        "tio_n_runs": len(tio_rows),
        "tio_fitted_slopes": slopes,
    }

    return {
        "metadata": _metadata(
            {
                "kappa": KAPPA,
                "kwave": KWAVE,
                "B0_rabi": B0_RABI,
                "bounds": list(BOUNDS),
                "t_end_default": T_END,
                "proca_mass_sq": PROCA_MASS_SQ,
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "results": tio_rows,
        "dispersion": dispersion_rows,
        "rabi_convergence": rabi_rows,
        "time_integration_order": tio_rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="convergence_rich_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
