"""Fundamental-vector (FV) ↔ TorsionCDT bit-exact equivalence.

For complex multi-field PGT theories where no analytic closed form is
available, the only available consistency check is to derive the same
continuum action through two mathematically distinct field
formulations — a 10-field fundamental-vector Proca representation and
an 18-field constrained-dyad (TorsionCDT) representation — and verify
that the two produce identical observables (P_max) at the IEEE round-off
limit.

Equivalence map (post-2026-04-24 NEW convention):
    m_T²(FV) = 2 · α_3(CDT)
    α_3 > 0  ↔  m_T² > 0   (stable Proca)

This benchmark runs both formulations at five canonical parameter
points spanning the dark-photon-plasma sector and reports the
pairwise P_max relative difference. Bit-exact agreement
(Δ/P_max ≲ 10⁻¹⁴) is the criterion.

The methodology has surfaced four distinct bugs (GitHub #305 spectral
leakage, #320 eigendecomposition ill-conditioning, #264 QZ
rank-deficiency, modal template-cache linearity violation), each
described in the App D §5 prose.

Serves:   manuscript/sections/appendices/validation.tex (App D §5)
Consumes: scripts/figures/figD_fv_cdt.py
Writes:   benchmark_results/canonical/fv_cdt_equivalence.json
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
DEFAULT_OUT = REPO_ROOT / "benchmark_results" / "canonical" / "fv_cdt_equivalence.json"

FV_SPEC = EXAMPLES / "data" / "torsion_dark_photon_fv.json"
CDT_SPEC = EXAMPLES / "data" / "dark_photon_plasma.json"

KAPPA = 1.0
B0 = 0.05
KWAVE = 2.0106
BOUNDS = (0.0, 100.0)
T_END = 50.0
SOURCE = "h_5"
TARGET = "a_1"
IC_COMPONENT = "h_5"

# Five canonical test points. The first is the original equivalence-
# verified point from memory: fv_cdt_equivalence_verified.md.
# Map convention: m_T² = 2·α_3 (post 2026-04-24 flip, GH #318).
FULL_POINTS = [
    {"label": "canonical", "mA2": 0.955, "deltam": 0.01, "xi": 0.274, "alpha3": 0.123},
    {"label": "low_mass", "mA2": 0.500, "deltam": 0.01, "xi": 0.500, "alpha3": 0.050},
    {"label": "high_mass", "mA2": 1.000, "deltam": 0.005, "xi": 0.300, "alpha3": 0.100},
    {"label": "wide_mix", "mA2": 0.700, "deltam": 0.020, "xi": 0.200, "alpha3": 0.150},
    {"label": "tight_mix", "mA2": 0.300, "deltam": 0.005, "xi": 0.400, "alpha3": 0.075},
]
SMOKE_POINTS = FULL_POINTS[:2]
FULL_GRID_N = 64
SMOKE_GRID_N = 32


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


def _simulate_and_measure(
    *, json_spec: Path, params: dict[str, float], grid_n: int, out_dir: Path
) -> tuple[bool, float, str]:
    sim_cmd = [
        "tidal",
        "simulate",
        str(json_spec),
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
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in params.items():
        sim_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        sim_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return (
            False,
            0.0,
            (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "sim failed",
        )

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
        "--json",
        "--quiet",
    ]
    for k, v in params.items():
        meas_cmd.extend(["--param", f"{k}={v}"])
    res = subprocess.run(
        meas_cmd, capture_output=True, text=True, cwd=REPO_ROOT, check=False
    )
    if res.returncode != 0:
        return (
            False,
            0.0,
            (res.stderr or res.stdout or "").splitlines()[-1]
            if (res.stderr or res.stdout)
            else "meas failed",
        )
    meas = json.loads(res.stdout)
    p = float(meas.get("conversion", {}).get("peak_probability", 0.0))
    return True, p, ""


def _run_point(point: dict, *, grid_n: int, work_dir: Path) -> dict:
    mA2 = point["mA2"]
    deltam = point["deltam"]
    xi = point["xi"]
    alpha3 = point["alpha3"]
    mT2 = 2.0 * alpha3  # NEW convention equivalence map

    fv_dir = work_dir / f"{point['label']}_fv"
    cdt_dir = work_dir / f"{point['label']}_cdt"
    for d in (fv_dir, cdt_dir):
        if d.exists():
            shutil.rmtree(d)

    common = {"kappa": KAPPA, "B0": B0, "mA2": mA2, "deltam": deltam, "xi": xi}
    fv_params = {**common, "mT2": mT2}
    cdt_params = {**common, "alpha3": alpha3}

    print(
        f"[fv_cdt_equivalence] {point['label']}: FV mT2={mT2}, CDT alpha3={alpha3}",
        flush=True,
    )
    ok_fv, p_fv, err_fv = _simulate_and_measure(
        json_spec=FV_SPEC, params=fv_params, grid_n=grid_n, out_dir=fv_dir
    )
    ok_cdt, p_cdt, err_cdt = _simulate_and_measure(
        json_spec=CDT_SPEC, params=cdt_params, grid_n=grid_n, out_dir=cdt_dir
    )
    if not (ok_fv and ok_cdt):
        return {
            "label": point["label"],
            "params": point,
            "ok": False,
            "fv_error": err_fv,
            "cdt_error": err_cdt,
        }
    denom = max(abs(p_fv), abs(p_cdt), 1e-30)
    return {
        "label": point["label"],
        "params": point,
        "mT2_equivalence": mT2,
        "ok": True,
        "P_max_FV": p_fv,
        "P_max_CDT": p_cdt,
        "abs_diff": abs(p_fv - p_cdt),
        "rel_diff": abs(p_fv - p_cdt) / denom,
    }


def run(*, smoke: bool, work_dir: Path) -> dict:
    points = SMOKE_POINTS if smoke else FULL_POINTS
    grid_n = SMOKE_GRID_N if smoke else FULL_GRID_N
    work_dir.mkdir(parents=True, exist_ok=True)

    rows = [_run_point(p, grid_n=grid_n, work_dir=work_dir) for p in points]
    ok = [r for r in rows if r.get("ok")]
    rels = [r["rel_diff"] for r in ok]

    summary = {
        "n_points": len(rows),
        "n_ok": len(ok),
        "max_rel_diff": max(rels) if rels else None,
        "median_rel_diff": float(np.median(rels)) if rels else None,
        "min_rel_diff": min(rels) if rels else None,
        "ieee_floor": 2.22e-16,
        "convention_map": "m_T² = 2·α_3 (post 2026-04-24, GH #318)",
    }

    return {
        "metadata": _metadata(
            {
                "n_points": len(points),
                "grid_n": grid_n,
                "kappa": KAPPA,
                "B0": B0,
                "kwave": KWAVE,
                "t_end": T_END,
                "bounds": list(BOUNDS),
                "fv_spec": str(FV_SPEC.relative_to(REPO_ROOT)),
                "cdt_spec": str(CDT_SPEC.relative_to(REPO_ROOT)),
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
    with tempfile.TemporaryDirectory(prefix="fv_cdt_") as tmp:
        work = Path(args.work_dir) if args.work_dir else Path(tmp)
        data = run(smoke=args.smoke, work_dir=work)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        json.dump(data, fh, indent=2)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
