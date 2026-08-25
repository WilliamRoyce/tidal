"""End-to-end contract of the observable-sector closure (GH #468 route 3).

On the localized implicit-dynamical class (E.cal), whose full pencil the
modal builder refuses at double precision (GH #468/#474), ``tidal
simulate`` restricts the run to the exactly closed sector excited by the
IC and writes it as a first-class ``restricted_spec.json`` that becomes
the run's ``spec_path``. The honest-data rules are the contract:

* fields outside the closure are ABSENT from the outputs (never zero,
  never present-but-wrong), and the provenance is recorded;
* ``tidal measure`` works unchanged inside the evolved sector;
* measurements that would read omitted fields — total energy with
  dropped Hamiltonian terms, conversion on an omitted source/target —
  error with the sector alternative named.
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
ECAL = REPO_ROOT / "examples/data/gertsenshtein_ungauged_e_dual_gaussian.json"
PARAMS = {"kappa": 1.0, "Bpeak": 0.01, "sigB": 5.0, "zc1": 25.0, "zc2": 75.0}

pytestmark = pytest.mark.slow


def _params() -> list[str]:
    out: list[str] = []
    for k, v in PARAMS.items():
        out += ["--param", f"{k}={v}"]
    return out


@pytest.fixture(scope="module")
def restricted_run(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if not ECAL.exists():
        pytest.skip("E.cal spec not present")
    out = tmp_path_factory.mktemp("ecal_closure") / "run"
    args = [
        "simulate",
        str(ECAL),
        "--grid-shape",
        "16",
        "--bounds",
        "0:100",
        "--periodic",
        "--ic",
        "plane-wave",
        "--ic-wavevector",
        "0.19",
        "--ic-amplitude",
        "0.01",
        "--ic-component",
        "h_5",
        "--t-end",
        "5",
        "--output",
        str(out),
        *_params(),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert main(args) == 0
    return out


class TestClosureRestrictedRun:
    def test_restricted_spec_artifact_and_absence(self, restricted_run: Path) -> None:
        spec_file = restricted_run / "restricted_spec.json"
        assert spec_file.exists()
        meta = json.loads((restricted_run / "metadata.json").read_text())
        assert Path(meta["spec_path"]) == spec_file
        # omitted fields are ABSENT from the run's outputs, not zero
        assert set(meta["fields"]) == {"h_5", "a_1"}
        record = json.loads(spec_file.read_text())["metadata"]["restriction"]
        assert record["evolved"] == sorted(record["evolved"]) or set(
            record["evolved"]
        ) == {"h_5", "a_1"}
        assert "h_3" in record["omitted"]
        assert record["seeds"] == ["h_5"]
        assert Path(record["parent_spec"]) == ECAL
        assert record["dropped_hamiltonian_terms"] > 0

    def test_measure_inside_the_sector_works(self, restricted_run: Path) -> None:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rc = main(
                [
                    "measure",
                    str(restricted_run),
                    "--what",
                    "conversion",
                    "--source",
                    "h_5",
                    "--target",
                    "a_1",
                    "--json",
                    "--quiet",
                    *_params(),
                ]
            )
        assert rc == 0
        payload = json.loads(buf.getvalue())
        assert payload["conversion"]["peak_probability"] >= 0.0

    def test_total_energy_refuses_with_sector_offer(
        self, restricted_run: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main(["measure", str(restricted_run), "--what", "energy", *_params()])
        assert rc == 1
        err = capsys.readouterr().err + capsys.readouterr().out
        assert "closure-restricted" in err
        assert "conversion" in err

    def test_omitted_field_refuses(
        self, restricted_run: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        rc = main(
            [
                "measure",
                str(restricted_run),
                "--what",
                "conversion",
                "--source",
                "h_3",
                "--target",
                "a_1",
                *_params(),
            ]
        )
        assert rc == 1
        err = capsys.readouterr().err + capsys.readouterr().out
        assert "not evolved" in err
        assert "h_3" in err
