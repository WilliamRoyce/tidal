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
import multiprocessing as mp
import os
import platform
import shutil
import socket
import subprocess  # noqa: S404
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
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

# (b) Rabi-frequency grid convergence. Use a stronger background
# B0_RABI = 0.2 so the Rabi period 2*pi/(kappa*B0) = 31.4 fits comfortably
# inside t_end = 200 (six full periods); the FFT then resolves the slow
# Rabi modulation cleanly, and Omega_eff/Omega_theory exposes the
# (k*dx)^2 FD-stencil error as N varies.
RABI_B0 = 0.2
RABI_T_END = 200.0
RABI_N_FULL = [128, 256, 512, 1024, 2048, 4096]
RABI_N_SMOKE = [128, 256]
RABI_SNAPSHOT_DT = 0.1  # ~2000 snapshots over t_end=200

# (c) Time-integration order params. Restricted to the symplectic family
# (leapfrog Y2, Y4) which is the only family where the slope is a
# meaningful convergence-order diagnostic on a wave problem at fixed
# spatial discretisation. Implicit schemes (CVODE, IDA, scipy DOP853)
# saturate at their rtol setting and exhibit a flat horizontal in this
# panel; their behaviour is already covered by the multi-method
# agreement figure.
TIO_SCHEMES_FULL: list[tuple[str, list[str]]] = [
    ("leapfrog_Y2", ["--scheme", "leapfrog", "--leapfrog-order", "2"]),
    ("leapfrog_Y4", ["--scheme", "leapfrog", "--leapfrog-order", "4"]),
]
TIO_SCHEMES_SMOKE = TIO_SCHEMES_FULL
TIO_CELL_TIMEOUT = 240  # seconds; safety net per (scheme, dt) cell
# 60 log-spaced dt values within the CFL-stable region (dt <= dx = 0.049
# at N=2048, L=100; safety cap at 0.024).  60 cells/scheme x 2 schemes =
# 120 parallel cells on sapphire (112 cores → 2 waves, ~60 s total).
TIO_DT_FULL: list[float] = [
    round(float(x), 7) for x in np.logspace(np.log10(0.001), np.log10(0.024), 60)
]
TIO_DT_SMOKE = [0.2, 0.1, 0.05]
TIO_GRID_N_FULL = 2048
TIO_GRID_N_SMOKE = 256
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


def _worker_init(avail_cores: list[int]) -> None:
    """Pin each worker to a private core to minimise cache thrashing."""
    worker_id = mp.current_process()._identity  # type: ignore[attr-defined]
    idx = (worker_id[0] - 1) % len(avail_cores) if worker_id else 0
    try:
        os.sched_setaffinity(0, {avail_cores[idx]})
    except (AttributeError, PermissionError, OSError):
        pass  # macOS dev or restricted env — degrade gracefully


def _dispersion_cell(cell: dict) -> dict:
    """Top-level function for ProcessPoolExecutor: one Proca dispersion cell."""
    m2 = cell["mass2"]
    k_req = cell["k_req"]
    sub = Path(cell["sub"])
    if sub.exists():
        shutil.rmtree(sub)
    try:
        _proca_simulate(
            k_request=k_req,
            mass2=m2,
            grid_n=cell["grid_n"],
            snapshot_dt=cell["snapshot_dt"],
            t_end=cell["t_end"],
            out_dir=sub,
        )
        omega_sim, k_real = _proca_omega(sub)
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as exc:
        return {"mass2": m2, "k_requested": k_req, "ok": False, "error": str(exc)}
    omega_ana = math.sqrt(k_real**2 + m2)
    return {
        "mass2": m2,
        "k_requested": k_req,
        "k_realised": k_real,
        "ok": True,
        "omega_sim": omega_sim,
        "omega_analytic": omega_ana,
        "rel_error": abs(omega_sim - omega_ana) / max(abs(omega_ana), 1e-30),
    }


def _run_dispersion(*, smoke: bool, work_dir: Path, max_workers: int = 1) -> list[dict]:
    mass_values = PROCA_MASS_SQ
    k_mults = PROCA_KMULT_SMOKE if smoke else PROCA_KMULT_FULL
    grid_n = PROCA_GRID_N_SMOKE if smoke else PROCA_GRID_N_FULL
    snapshot_dt = PROCA_SNAP_DT_SMOKE if smoke else PROCA_SNAP_DT_FULL
    t_end = PROCA_T_END_SMOKE if smoke else PROCA_T_END_FULL
    base_k = 2.0 * math.pi / (BOUNDS[1] - BOUNDS[0])
    cells = [
        {
            "mass2": m2,
            "k_req": base_k * kmult,
            "grid_n": grid_n,
            "snapshot_dt": snapshot_dt,
            "t_end": t_end,
            "sub": str(work_dir / f"proca_m{str(m2).replace('.', 'p')}_k{kmult}"),
        }
        for m2 in mass_values
        for kmult in k_mults
    ]
    if max_workers > 1:
        avail = (
            sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else list(range(os.cpu_count() or 1))
        )
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=min(max_workers, len(cells)),
            mp_context=ctx,
            initializer=_worker_init,
            initargs=(avail,),
        ) as ex:
            return list(ex.map(_dispersion_cell, cells))
    return [_dispersion_cell(c) for c in cells]


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
        f"{RABI_T_END}",
        "--snapshots",
        f"{snapshot_dt}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={RABI_B0}",
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


def _rabi_cell(cell: dict) -> dict:
    """Top-level function for ProcessPoolExecutor: one Rabi-N cell."""
    n = cell["N"]
    sub = Path(cell["sub"])
    omega_theory = cell["omega_theory"]
    if sub.exists():
        shutil.rmtree(sub)
    try:
        _rabi_simulate(grid_n=n, out_dir=sub, snapshot_dt=cell["snapshot_dt"])
        omega_eff = _rabi_omega_eff(sub)
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as exc:
        return {"N": n, "ok": False, "error": str(exc)}
    return {
        "N": n,
        "ok": True,
        "omega_eff": omega_eff,
        "omega_theory": omega_theory,
        "ratio": omega_eff / max(omega_theory, 1e-30),
        "k_dx": KWAVE * (BOUNDS[1] - BOUNDS[0]) / n,
    }


def _run_rabi(*, smoke: bool, work_dir: Path, max_workers: int = 1) -> list[dict]:
    n_values = RABI_N_SMOKE if smoke else RABI_N_FULL
    omega_theory = KAPPA * RABI_B0
    cells = [
        {
            "N": n,
            "omega_theory": omega_theory,
            "snapshot_dt": RABI_SNAPSHOT_DT,
            "sub": str(work_dir / f"rabi_N{n}"),
        }
        for n in n_values
    ]
    if max_workers > 1:
        avail = (
            sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else list(range(os.cpu_count() or 1))
        )
        ctx = mp.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=min(max_workers, len(cells)),
            mp_context=ctx,
            initializer=_worker_init,
            initargs=(avail,),
        ) as ex:
            return list(ex.map(_rabi_cell, cells))
    return [_rabi_cell(c) for c in cells]


# -------- (c) Time-integration order --------


def _tio_modal_reference(
    *, grid_n: int, out_dir: Path, t_end: float = T_END
) -> float | None:
    """Run modal at t_end to obtain a per-cell aligned P reference.

    Pass t_end = math.ceil(T_END / dt) * dt so modal reports P at the same
    simulation time as the leapfrog cell, eliminating the snapshot-misalignment
    term (dt · dP/dt) that otherwise dominates the small-dt residual and masks
    the Y4 slope-4 signal.
    """
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
        f"{t_end}",
        "--scheme",
        "modal",
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
        return None
    if res.returncode != 0:
        return None
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
        return None
    meas = json.loads(res.stdout)
    return float(meas.get("conversion", {}).get("peak_probability", 0.0))


def _tio_run_one(
    scheme_label: str,
    scheme_args: list[str],
    dt: float,
    *,
    grid_n: int,
    out_dir: Path,
    modal_reference: float | None = None,
    t_align: float | None = None,
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
        "peak_conversion",
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
    # P_final = probability at last snapshot = P(t_align); robust regardless of
    # snapshot count (unlike peak_probability which is argmax and can pick an
    # intermediate snapshot when Y4's oscillatory error temporarily overshoots).
    p = float(meas.get("peak_conversion", {}).get("P_final", 0.0))
    err_modal = abs(p - modal_reference) if modal_reference is not None else None
    t_ref = t_align if t_align is not None else T_END
    p_analytic_aligned = math.sin(0.5 * KAPPA * B0_RABI * t_ref) ** 2
    return {
        "scheme": scheme_label,
        "dt": dt,
        "t_align": t_ref,
        "ok": True,
        "P_final_sim": p,
        "P_final_analytic": TIO_ANALYTIC,
        "P_final_analytic_aligned": p_analytic_aligned,
        "P_final_modal": modal_reference,
        "abs_error": abs(p - TIO_ANALYTIC),
        "abs_error_analytic_aligned": abs(p - p_analytic_aligned),
        "abs_error_vs_modal": err_modal,
    }


def _tio_cell(cell: dict) -> dict:
    """Top-level function for ProcessPoolExecutor: one (scheme, dt) TIO cell.

    Runs modal at t_align = ceil(T_END/dt)*dt, then leapfrog at T_END with
    the same dt.  Using a scheme-specific modal dir avoids write conflicts when
    Y2 and Y4 cells for the same dt run concurrently.
    """
    label = cell["label"]
    scheme_args = cell["scheme_args"]
    dt = cell["dt"]
    grid_n = cell["grid_n"]
    work_dir = Path(cell["work_dir"])

    t_align = math.ceil(T_END / dt) * dt
    dt_tag = str(dt).replace(".", "p")
    modal_dir = work_dir / f"tio_modal_{label}_dt{dt_tag}"
    if modal_dir.exists():
        shutil.rmtree(modal_dir)
    p_modal = _tio_modal_reference(grid_n=grid_n, out_dir=modal_dir, t_end=t_align)
    print(
        f"[tio] {label} dt={dt:.6f} t_align={t_align:.6f} P_modal={p_modal!r}",
        flush=True,
    )
    sub = work_dir / f"tio_{label}_dt{dt_tag}"
    if sub.exists():
        shutil.rmtree(sub)
    return _tio_run_one(
        label,
        scheme_args,
        dt,
        grid_n=grid_n,
        out_dir=sub,
        modal_reference=p_modal,
        t_align=t_align,
    )


def _run_tio(*, smoke: bool, work_dir: Path, max_workers: int = 1) -> list[dict]:
    schemes = TIO_SCHEMES_SMOKE if smoke else TIO_SCHEMES_FULL
    dts = TIO_DT_SMOKE if smoke else TIO_DT_FULL
    grid_n = TIO_GRID_N_SMOKE if smoke else TIO_GRID_N_FULL
    cells = [
        {
            "label": label,
            "scheme_args": scheme_args,
            "dt": dt,
            "grid_n": grid_n,
            "work_dir": str(work_dir),
        }
        for label, scheme_args in schemes
        for dt in dts
    ]
    if max_workers <= 1:
        return [_tio_cell(c) for c in cells]

    ctx = mp.get_context("fork")
    rows: list[dict] = []
    avail = (
        sorted(os.sched_getaffinity(0))
        if hasattr(os, "sched_getaffinity")
        else list(range(os.cpu_count() or 1))
    )
    n_workers = min(max_workers, len(cells))
    with ProcessPoolExecutor(
        max_workers=n_workers,
        mp_context=ctx,
        initializer=_worker_init,
        initargs=(avail,),
    ) as pool:
        futures = {pool.submit(_tio_cell, c): c for c in cells}
        for n_done, fut in enumerate(as_completed(futures), start=1):
            rows.append(fut.result())
            if n_done % 10 == 0 or n_done == len(cells):
                print(f"[tio] {n_done}/{len(cells)} cells done", flush=True)
    return rows


def _fit_slope(
    rows: list[dict], scheme: str, *, key: str = "abs_error_vs_modal"
) -> float | None:
    """Fit log-log slope of the given error metric vs dt.

    Default key is the modal-reference error so spatial/IC floors cancel;
    falls back to "abs_error" (vs sin² analytic) if modal reference is
    unavailable.
    """
    ok = [
        r
        for r in rows
        if r.get("ok")
        and r.get("scheme") == scheme
        and r.get(key) is not None
        and math.isfinite(r[key])
        and r[key] > 1e-14
    ]
    if len(ok) < 2:
        # Fall back to vs-analytic if vs-modal is empty.
        if key != "abs_error":
            return _fit_slope(rows, scheme, key="abs_error")
        return None
    dts = np.log10([r["dt"] for r in ok])
    errs = np.log10([r[key] for r in ok])
    slope, _ = np.polyfit(dts, errs, 1)
    return float(slope)


def run(*, smoke: bool, work_dir: Path, max_workers: int = 1) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    dispersion_rows = _run_dispersion(
        smoke=smoke, work_dir=work_dir, max_workers=max_workers
    )
    rabi_rows = _run_rabi(smoke=smoke, work_dir=work_dir, max_workers=max_workers)
    tio_rows = _run_tio(smoke=smoke, work_dir=work_dir, max_workers=max_workers)

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


def append_tio(*, out: Path, work_dir: Path, max_workers: int = 1) -> None:
    """Append missing (scheme, dt) time_integration_order rows.

    Loads the existing JSON, identifies (scheme, dt) pairs in
    TIO_SCHEMES_FULL x TIO_DT_FULL that are not already present, runs only
    those (via _tio_cell so per-cell aligned modal reference is applied), and
    merges into the row list. Existing rows survive untouched.

    Seeds are loaded ONCE before any cells run; the output file is only
    written after all new rows are collected, so re-reading the file mid-run
    cannot cause double-counting.

    Raises
    ------
    FileNotFoundError
        If ``out`` does not exist (full run required first).
    """
    if not out.exists():
        msg = f"--append-tio requires existing {out}; run the full benchmark first."
        raise FileNotFoundError(msg)
    # Load seeds once; never re-read until after all new rows are ready.
    with out.open(encoding="utf-8") as fh:
        data = json.load(fh)
    existing = data.get("time_integration_order", [])
    have = {(r.get("scheme"), float(r["dt"])) for r in existing}

    work_dir.mkdir(parents=True, exist_ok=True)
    grid_n = TIO_GRID_N_FULL
    cells = [
        {
            "label": label,
            "scheme_args": scheme_args,
            "dt": dt,
            "grid_n": grid_n,
            "work_dir": str(work_dir),
        }
        for label, scheme_args in TIO_SCHEMES_FULL
        for dt in TIO_DT_FULL
        if (label, dt) not in have
    ]
    if not cells:
        print("append_tio: nothing to do (all cells already present)")
        return

    new_rows: list[dict] = []
    if max_workers <= 1:
        new_rows = [_tio_cell(c) for c in cells]
    else:
        ctx = mp.get_context("fork")
        avail = (
            sorted(os.sched_getaffinity(0))
            if hasattr(os, "sched_getaffinity")
            else list(range(os.cpu_count() or 1))
        )
        n_workers = min(max_workers, len(cells))
        with ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=ctx,
            initializer=_worker_init,
            initargs=(avail,),
        ) as pool:
            futures = {pool.submit(_tio_cell, c): c for c in cells}
            for n_done, fut in enumerate(as_completed(futures), start=1):
                new_rows.append(fut.result())
                if n_done % 10 == 0 or n_done == len(cells):
                    print(f"[append_tio] {n_done}/{len(cells)} cells done", flush=True)

    merged = list(existing) + new_rows
    merged.sort(key=lambda r: (r.get("scheme", ""), float(r["dt"])))
    data["time_integration_order"] = merged
    with out.open("w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)
    print(f"appended {len(new_rows)} TIO rows to {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Run benchmark cells concurrently using N worker processes "
            "(fork context; set to 112 on a sapphire INTR node)."
        ),
    )
    parser.add_argument(
        "--append-tio",
        action="store_true",
        help=(
            "Run only the time_integration_order sweep for (scheme, dt) "
            "pairs missing from the existing JSON; preserves all other rows."
        ),
    )
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="convergence_rich_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        if args.append_tio:
            append_tio(out=args.out, work_dir=work, max_workers=args.parallel)
            return
        data = run(smoke=args.smoke, work_dir=work, max_workers=args.parallel)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    # Force single-threaded BLAS in every worker subprocess so that
    # N workers x M BLAS threads don't fight over the same cores.
    # Must be set before ProcessPoolExecutor forks workers (fork inherits env).
    for _var in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "BLIS_NUM_THREADS",
    ):
        os.environ.setdefault(_var, "1")
    main()
