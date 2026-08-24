"""In-suite FV ↔ TorsionCDT equivalence pin (coverage-audit closure, GH #457).

The fundamental-vector (10-field Proca) and TorsionCDT (18-field
constrained-dyad) formulations derive the SAME continuum physics through
mathematically distinct field representations, related by the
post-2026-04-24 convention map ``mT2 = 2·alpha3``. Their P_max agreement
at the round-off limit is the adjudicating instrument that historically
surfaced #305, #320, #264 and the template-cache violation — yet until
now it lived only in ``scripts/benchmarks/fv_cdt_equivalence.py`` with
no in-suite pin (2026-08-26 historical audit on GH #457).

This test runs the benchmark's canonical point at reduced size and pins
the equivalence. MEASURED under the pencil-engine operator (commit
372e4a73, N=32, t_end=20): relative difference 4.6e-14 — the map is
representation-level and survives the corrected operator exactly (both
formulations changed together; the CDT side additionally exercises the
gauge quotient for its redundant constraint triplets, GH #465). The
1e-12 assertion keeps ~40x headroom over the measured value while still
catching any representation-splitting regression at the first digit
that could matter.

A FAILURE here means the two formulations' operators diverged — a
finding for GH #449/#468 adjudication, not a tolerance to loosen.
"""

from __future__ import annotations

import contextlib
import io
import json
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FV_SPEC = REPO_ROOT / "examples/data/torsion_dark_photon_fv.json"
CDT_SPEC = REPO_ROOT / "examples/data/dark_photon_plasma.json"

#: Canonical point of scripts/benchmarks/fv_cdt_equivalence.py.
COMMON = {"kappa": 1.0, "B0": 0.05, "mA2": 0.955, "deltam": 0.01, "xi": 0.274}
ALPHA3 = 0.123
MT2 = 2.0 * ALPHA3  # the NEW-convention equivalence map

pytestmark = pytest.mark.slow


def _simulate_and_measure(
    spec: Path, params: dict[str, float], out_dir: Path
) -> float:
    from tidal.cli import main

    sim_args = [
        "simulate",
        str(spec),
        "--grid-shape",
        "32",
        "--bounds",
        "0:100",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        "2.0106",
        "--ic-amplitude",
        "0.1",
        "--ic-component",
        "h_5",
        "--t-end",
        "20",
        "--fd-order",
        "4",
        "--output",
        str(out_dir),
        "--force",
    ]
    for k, v in params.items():
        sim_args += ["--param", f"{k}={v}"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(sim_args) == 0

    meas_args = [
        "measure",
        str(out_dir),
        "--what",
        "conversion",
        "--source",
        "h_5",
        "--target",
        "a_1",
        "--json",
        "--quiet",
    ]
    for k, v in params.items():
        meas_args += ["--param", f"{k}={v}"]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(meas_args) == 0
    payload = json.loads(buf.getvalue())
    return float(payload["conversion"]["peak_probability"])


def test_fv_cdt_pmax_equivalence(tmp_path: Path) -> None:
    if not (FV_SPEC.exists() and CDT_SPEC.exists()):
        pytest.skip("FV/CDT spec pair not present")
    p_fv = _simulate_and_measure(FV_SPEC, {**COMMON, "mT2": MT2}, tmp_path / "fv")
    p_cdt = _simulate_and_measure(
        CDT_SPEC, {**COMMON, "alpha3": ALPHA3}, tmp_path / "cdt"
    )
    assert p_cdt > 1e-4  # the point genuinely converts — no trivial-zero pass
    rel = abs(p_fv - p_cdt) / p_cdt
    assert rel < 1e-12, (
        f"FV/CDT representations diverged: P_fv={p_fv:.16e} "
        f"P_cdt={p_cdt:.16e} rel={rel:.3e} — adjudicate before touching "
        f"this tolerance (GH #449/#468)"
    )
