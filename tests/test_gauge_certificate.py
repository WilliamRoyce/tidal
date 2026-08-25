"""GH #468 "pin + certify": the measurement-time gauge certificate.

The uniform dark-photon CDT spec ships byte-identical constraint triplets
(GH #465), so the modal engine's gauge quotient pins undetermined
directions on every mode. The run's ``metadata.json`` carries, per state
slot, how much of the pinned subspace lives there
(``solver_diagnostics``); ``tidal measure`` turns that into a
per-observable certificate:

* the h_5 → a_1 conversion reads no pinned content → ``certified`` (the
  FV↔CDT agreement at 4.6e-14 is the manual proof this automates);
* the whole-state energy reads the pinned torsion slots → ``flagged``
  with the magnitude — reported, never silently presented as invariant.
"""

from __future__ import annotations

import contextlib
import io
import json
import warnings
from pathlib import Path

import pytest

from tidal.cli import main

REPO_ROOT = Path(__file__).resolve().parent.parent
CDT = REPO_ROOT / "examples/data/dark_photon_plasma.json"
PARAMS = {
    "kappa": 1.0,
    "B0": 0.05,
    "mA2": 0.955,
    "deltam": 0.01,
    "xi": 0.274,
    "alpha3": 0.123,
}

pytestmark = pytest.mark.slow


def _params() -> list[str]:
    out: list[str] = []
    for k, v in PARAMS.items():
        out += ["--param", f"{k}={v}"]
    return out


def _measure_json(run: Path, what: str, extra: list[str]) -> dict:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rc = main(
            [
                "measure",
                str(run),
                "--what",
                what,
                "--json",
                "--quiet",
                *extra,
                *_params(),
            ]
        )
    assert rc == 0
    return json.loads(buf.getvalue())


@pytest.fixture(scope="module")
def cdt_run(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not CDT.exists():
        pytest.skip("dark_photon_plasma spec not present")
    out = tmp_path_factory.mktemp("cdt_cert") / "run"
    args = [
        "simulate",
        str(CDT),
        "--grid-shape",
        "16",
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
        "5",
        "--output",
        str(out),
        "--force",
        *_params(),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(args) == 0
    return out


class TestGaugeCertificate:
    def test_run_records_pinned_directions(self, cdt_run: Path) -> None:
        meta = json.loads((cdt_run / "metadata.json").read_text())
        diag = meta["solver_diagnostics"]
        assert diag["pinned_dims"] > 0  # the #465 triplets are quotiented
        overlap = dict(zip(diag["slot_names"], diag["pin_overlap_max"], strict=True))
        # the pins live on the redundant torsion sector, not on the channel
        assert overlap["h_5"] < 1e-10
        assert overlap["a_1"] < 1e-10
        assert overlap["v_h_5"] < 1e-10
        assert overlap["v_a_1"] < 1e-10
        assert max(overlap.values()) > 0.1

    def test_conversion_is_certified(self, cdt_run: Path) -> None:
        payload = _measure_json(
            cdt_run, "conversion", ["--source", "h_5", "--target", "a_1"]
        )
        cert = payload["gauge_certificate"]
        assert cert["verdict"] == "certified"
        assert cert["max_pin_overlap"] < 1e-10
        assert cert["pinned_dims"] > 0

    def test_whole_state_energy_is_flagged(self, cdt_run: Path) -> None:
        payload = _measure_json(cdt_run, "energy", [])
        cert = payload["gauge_certificate"]
        assert cert["verdict"] == "flagged"
        assert cert["max_pin_overlap"] > 0.1
        assert cert["pinned_slots_in_support"]

    def test_text_output_carries_the_certificate(
        self, cdt_run: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rc = main(
                [
                    "measure",
                    str(cdt_run),
                    "--what",
                    "conversion",
                    "--source",
                    "h_5",
                    "--target",
                    "a_1",
                    *_params(),
                ]
            )
        assert rc == 0
        out = capsys.readouterr().out
        assert "Gauge certificate" in out
        assert "certified" in out
