"""GH #468 "pin + certify": the measurement-time gauge certificate.

The uniform dark-photon CDT spec ships byte-identical constraint triplets
(GH #465), so the modal engine's gauge quotient pins undetermined
directions on every mode. The run's ``metadata.json`` carries, per state
slot, how much of the pinned subspace lives there
(``solver_diagnostics``); ``tidal measure`` turns that into a
per-observable certificate:

* the h_5 → a_1 conversion is ``certified`` by the direct
  perturb-and-remeasure test (the FV↔CDT agreement at 4.6e-14 is the
  manual proof this automates);
* the whole-state energy READS the pinned torsion slots (slot overlap
  would flag it) yet is invariant under the direct test — pure-gauge
  directions carry no energy — so it is ``certified``;
* an observable that reads a pinned field's raw value is ``flagged``
  with its sensitivity — reported, never presented as invariant.
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

    def test_run_stores_pinned_subspace_probes(self, cdt_run: Path) -> None:
        import numpy as np

        probes = np.load(cdt_run / "pin_probes.npy")
        meta = json.loads((cdt_run / "metadata.json").read_text())
        n_probe_axes = 3  # (K, n_slots, N) in 1-D
        assert probes.ndim == n_probe_axes
        assert probes.shape[1] == len(meta["solver_diagnostics"]["slot_names"])
        assert np.all(np.isfinite(probes))
        assert np.max(np.abs(probes)) == 1.0  # normalized

    def test_conversion_is_certified(self, cdt_run: Path) -> None:
        """Direct test (arc C2): perturbing the state along the pinned
        subspace does not move the h_5 → a_1 peak conversion.
        """
        payload = _measure_json(
            cdt_run, "conversion", ["--source", "h_5", "--target", "a_1"]
        )
        cert = payload["gauge_certificate"]
        assert cert["method"] == "perturb-and-remeasure"
        assert cert["verdict"] == "certified"
        assert cert["max_sensitivity"] < 1e-6
        assert cert["pinned_dims"] > 0

    def test_whole_state_energy_is_certified_pins_are_pure_gauge(
        self, cdt_run: Path
    ) -> None:
        """The energy reads the pinned torsion slots, yet the direct test
        finds it INVARIANT: the pinned #465 redundancy directions carry no
        energy, as pure-gauge directions must. The slot-overlap method
        flagged this case — the over-flagging the direct test exists to
        remove — and a passing energy certificate doubles as the check
        that the pinned directions really were gauge.
        """
        payload = _measure_json(cdt_run, "energy", [])
        cert = payload["gauge_certificate"]
        assert cert["method"] == "perturb-and-remeasure"
        assert cert["verdict"] == "certified"
        assert cert["max_sensitivity"] < 1e-6
        assert cert["max_pin_overlap"] > 0.1  # slot overlap alone would flag

    def test_reading_a_pinned_field_directly_is_flagged(self, cdt_run: Path) -> None:
        """Positive control: an observable that reads a pinned field's raw
        value (the peak of the redundant torsion component) DOES move under
        the perturbation — the certificate must say so.
        """
        import numpy as np

        from tidal.cli._measure import _gauge_sensitivity
        from tidal.measurement import SimulationData
        from tidal.symbolic import load_equation_system

        meta = json.loads((cdt_run / "metadata.json").read_text())
        data = SimulationData.load(
            cdt_run, load_equation_system(Path(meta["spec_path"]))
        )
        slot_names = meta["solver_diagnostics"]["slot_names"]
        overlap = dict(
            zip(slot_names, meta["solver_diagnostics"]["pin_overlap_max"], strict=True)
        )
        pinned_field = max(
            (s for s in slot_names if not s.startswith("v_")), key=lambda s: overlap[s]
        )

        def _peak_of_pinned_field(d: SimulationData) -> float:
            return float(np.max(np.abs(d.fields[pinned_field])))

        sens = _gauge_sensitivity(
            data,
            slot_names,
            cdt_run / "pin_probes.npy",
            _peak_of_pinned_field,
            f"peak |{pinned_field}|",
        )
        assert sens is not None
        assert sens["max_sensitivity"] > 1e-3

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
