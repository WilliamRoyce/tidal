r"""Proca dispersion calibration: $\\omega^2 = k^2 + m^2$.

Calibrates TIDAL's time evolution of a massive vector (Proca) plane-wave
against the closed-form dispersion relation. We use the
`gertsenshtein_proca` example at $\\Bzero = 0$ (decoupled from gravity)
to obtain the free Proca equation $\\partial_t^2 a = \\partial_x^2 a -
m^2 a$. For each $k$ in the sweep we run a plane-wave IC at that
wavevector, sample the time series of one field component finely, FFT
to recover the dominant angular frequency $\\omega(k)$, and compare to
the analytic prediction $\\omega(k) = \\sqrt{k^2 + m^2}$.

Serves:   manuscript/sections/appendices/validation.tex (App D, calibration 4)
Consumes: scripts/figures/figD_dispersion_time_order.py
Writes:   benchmark_results/canonical/proca_dispersion.json
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "proca_dispersion.json"

PROCA_SPEC = EXAMPLES / "data" / "gertsenshtein_proca.json"

# Physics parameters.
KAPPA = 1.0
B0 = 0.0  # decouple a from h
MASS2_FULL = 1.0
MASS2_SMOKE = 1.0
T_END_FULL = 100.0  # ~16 periods of slowest mode at m=1
T_END_SMOKE = 8.0
BOUNDS = (0.0, 100.0)  # large enough that several wavelengths fit

# Snapshots: dense enough to resolve the highest expected omega.
SNAPSHOT_DT_FULL = 0.05  # 2000 snapshots over t_end=100
SNAPSHOT_DT_SMOKE = 0.05  # 160 snapshots over t_end=8

# k-sweep, in units of 2*pi/L base modes — converted to physical k below.
FULL_K_MULTIPLIERS = [4, 8, 16, 24, 32, 48]
SMOKE_K_MULTIPLIERS = [4, 8, 16]

GRID_N_FULL = 1024
GRID_N_SMOKE = 256


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


def _simulate(
    *,
    k_request: float,
    grid_n: int,
    snapshot_dt: float,
    t_end: float,
    mass2: float,
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
        "a_1",  # photon component (a_x)
        "--t-end",
        f"{t_end}",
        "--snapshots",
        f"{snapshot_dt}",
        "--param",
        f"kappa={KAPPA}",
        "--param",
        f"B0={B0}",
        "--param",
        f"mA2={mass2}",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    print(f"[proca_dispersion] k_req={k_request}: simulate", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def _omega_from_snapshots(out_dir: Path, component: str = "a_1") -> tuple[float, float]:
    """Return (omega, k_realised) extracted from snapshots.

    omega: angular frequency from FFT of the spatial-mean(|field|) time series.
    k_realised: wavenumber of the dominant spatial Fourier mode at t=0.
    """
    times = np.load(out_dir / "times.npy")
    field = np.load(out_dir / f"{component}.npy")  # shape (Nt, Nx)
    bounds = BOUNDS
    nx = field.shape[1]
    x_axis_len = bounds[1] - bounds[0]

    # Realised k from t=0 snapshot
    spec0 = np.fft.rfft(field[0])
    kx_axis = np.fft.rfftfreq(nx, d=x_axis_len / nx) * 2.0 * math.pi
    k_realised = float(kx_axis[int(np.argmax(np.abs(spec0)))])

    # Time series: amplitude of the dominant spatial mode over t.
    mode_idx = int(np.argmax(np.abs(spec0)))
    mode_amp = np.array(
        [np.fft.rfft(field[i])[mode_idx] for i in range(field.shape[0])]
    )
    # Use the real part (cosine response) as the time series.
    series = mode_amp.real

    # Detrend (remove DC) and apply Hann window before FFT.
    series -= series.mean()
    window = np.hanning(series.size)
    windowed = series * window

    # FFT of the time series — note times may not be uniformly spaced in
    # general, but tidal simulate uses uniform snapshots so dt is stable.
    dt = float(times[1] - times[0]) if times.size >= 2 else 1.0
    spec = np.fft.rfft(windowed)
    freqs = np.fft.rfftfreq(series.size, d=dt)  # Hz (1/T)
    # Mask DC bin
    spec[0] = 0
    mag = np.abs(spec)
    peak_idx = int(np.argmax(mag))
    f_peak = float(freqs[peak_idx])
    # Quadratic interpolation around the peak (sub-bin refinement).
    if 1 <= peak_idx < mag.size - 1:
        y0, y1, y2 = mag[peak_idx - 1], mag[peak_idx], mag[peak_idx + 1]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            delta = 0.5 * (y0 - y2) / denom
            df = float(freqs[1] - freqs[0]) if freqs.size >= 2 else 0.0
            f_peak = float(freqs[peak_idx]) + delta * df
    omega = 2.0 * math.pi * f_peak
    return omega, k_realised


def run(*, smoke: bool, work_dir: Path) -> dict:
    k_multipliers = SMOKE_K_MULTIPLIERS if smoke else FULL_K_MULTIPLIERS
    grid_n = GRID_N_SMOKE if smoke else GRID_N_FULL
    snapshot_dt = SNAPSHOT_DT_SMOKE if smoke else SNAPSHOT_DT_FULL
    t_end = T_END_SMOKE if smoke else T_END_FULL
    mass2 = MASS2_SMOKE if smoke else MASS2_FULL

    base_k = 2.0 * math.pi / (BOUNDS[1] - BOUNDS[0])  # 2π/L
    k_requests = [base_k * m for m in k_multipliers]

    work_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    for k_req in k_requests:
        sub = work_dir / f"k_{k_req:.4f}".replace(".", "p")
        if sub.exists():
            shutil.rmtree(sub)
        try:
            _simulate(
                k_request=k_req,
                grid_n=grid_n,
                snapshot_dt=snapshot_dt,
                t_end=t_end,
                mass2=mass2,
                out_dir=sub,
            )
            omega, k_real = _omega_from_snapshots(sub)
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as exc:
            rows.append(
                {
                    "k_requested": k_req,
                    "ok": False,
                    "error": str(exc),
                }
            )
            continue
        omega_analytic = math.sqrt(k_real**2 + mass2)
        rows.append(
            {
                "k_requested": k_req,
                "k_realised": k_real,
                "ok": True,
                "omega_sim": omega,
                "omega_analytic": omega_analytic,
                "rel_error": (
                    abs(omega - omega_analytic) / max(abs(omega_analytic), 1e-12)
                ),
            }
        )

    ok = [r for r in rows if r.get("ok")]
    rels = [r["rel_error"] for r in ok]
    summary = {
        "n_k_points": len(rows),
        "n_ok": len(ok),
        "max_rel_error": max(rels) if rels else None,
        "median_rel_error": float(np.median(rels)) if rels else None,
    }

    return {
        "metadata": _metadata(
            {
                "k_multipliers": k_multipliers,
                "grid_n": grid_n,
                "snapshot_dt": snapshot_dt,
                "t_end": t_end,
                "mass2": mass2,
                "bounds": list(BOUNDS),
                "kappa": KAPPA,
                "B0": B0,
                "smoke": smoke,
            }
        ),
        "summary": summary,
        "k_sweep": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--work-dir", type=Path, default=None)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="proca_dispersion_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
