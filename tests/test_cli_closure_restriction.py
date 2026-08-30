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


# ---------------------------------------------------------------------------
# The physics pin: isolate the graviton–photon sector, evolve only it, and
# get the RIGHT conversion (GH #468/#449, arc C of the 2026-08-27 plan).
# ---------------------------------------------------------------------------

#: Frozen Phase E geometry (scripts/hpc_submit_drafts/v3e_localised/_geometry.env).
PHASE_E_ARGS = [
    "--grid-shape",
    "128",
    "--bounds",
    "0:100",
    "--periodic",
    "--ic",
    "gaussian",
    "--ic-component",
    "h_5",
    "--ic-center=8",
    "--ic-width",
    "3",
    "--ic-wavevector",
    "2",
    "--ic-amplitude",
    "0.01",
]


def _boccaletti_prediction() -> float:
    """Path-integrated Boccaletti P = sin²(κ·Bpeak·σ_B·√(2π)/2) per Gaussian."""
    import math

    arg = (
        PARAMS["kappa"] * PARAMS["Bpeak"] * PARAMS["sigB"] * math.sqrt(2 * math.pi) / 2
    )
    return math.sin(arg) ** 2


def _phase_e_run(out: Path, t_end: float) -> Path:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert (
            main(
                [
                    "simulate",
                    str(ECAL),
                    *PHASE_E_ARGS,
                    "--t-end",
                    str(t_end),
                    "--output",
                    str(out),
                    "--force",
                    *_params(),
                ]
            )
            == 0
        )
    return out


def _conversion_series(run: Path):
    from tidal.measurement import SimulationData
    from tidal.measurement._conversion import compute_conversion_probability
    from tidal.symbolic import load_equation_system

    spec_path = Path(json.loads((run / "metadata.json").read_text())["spec_path"])
    data = SimulationData.load(run, load_equation_system(spec_path))
    return compute_conversion_probability(data, "h_5", "a_1")


@pytest.fixture(scope="module")
def phase_e_runs(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    if not ECAL.exists():
        pytest.skip("E.cal spec not present")
    base = tmp_path_factory.mktemp("ecal_phase_e")
    return _phase_e_run(base / "t40", 40.0), _phase_e_run(base / "t80", 80.0)


class TestClosureReproducesBoccaletti:
    """The load-bearing physics claim of the closure route, pinned.

    Measured 2026-08-27 on the corrected operator: P_peak = 0.003912 vs
    sin²(0.0627) = 0.003922 (ratio 0.9975), P_peak(t_end=80) = 0.003914
    at t = 44, final P after the sign-flipped second Gaussian 7.8e-6.
    The archived Phase E value 0.0036 (defective operator) is ~8% low —
    a 1% tolerance separates the two cleanly.
    """

    def test_only_the_graviton_photon_sector_is_evolved(
        self, phase_e_runs: tuple[Path, Path]
    ) -> None:
        run40, _ = phase_e_runs
        meta = json.loads((run40 / "metadata.json").read_text())
        assert set(meta["fields"]) == {"h_5", "a_1"}
        assert Path(meta["spec_path"]).name == "restricted_spec.json"

    def test_peak_conversion_matches_boccaletti(
        self, phase_e_runs: tuple[Path, Path]
    ) -> None:
        run40, _ = phase_e_runs
        conv = _conversion_series(run40)
        p_peak = float(conv.probability.max())
        p_th = _boccaletti_prediction()
        assert abs(p_peak - p_th) / p_th < 1e-2, (
            f"P_peak={p_peak:.6e} vs Boccaletti {p_th:.6e} "
            f"(ratio {p_peak / p_th:.4f}) — the archived defective-operator "
            f"value was 0.0036"
        )

    def test_peak_is_t_end_independent(self, phase_e_runs: tuple[Path, Path]) -> None:
        run40, run80 = phase_e_runs
        p40 = float(_conversion_series(run40).probability.max())
        p80 = float(_conversion_series(run80).probability.max())
        assert abs(p80 - p40) / p40 < 5e-3, (
            f"A(80)/A(40) = {p80 / p40:.5f} — growth would be the #455 "
            f"spurious-operator signature"
        )

    def test_second_gaussian_cancels_the_conversion(
        self, phase_e_runs: tuple[Path, Path]
    ) -> None:
        # B0(x) = Bpeak·(g1 − g2): ∫B dz = 0 over the box, so after the
        # packet crosses the sign-flipped second Gaussian the accumulated
        # phase cancels and P returns to ~0 — the dual-Gaussian periodic
        # construction's intended numerical property.
        _, run80 = phase_e_runs
        conv = _conversion_series(run80)
        assert float(conv.probability[-1]) < 1e-4 * _boccaletti_prediction() * 1e2
        assert float(conv.probability[-1]) < 1e-5

    def test_measurement_is_gauge_certified(
        self, phase_e_runs: tuple[Path, Path]
    ) -> None:
        run40, _ = phase_e_runs
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rc = main(
                [
                    "measure",
                    str(run40),
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
        assert payload["gauge_certificate"]["verdict"] == "certified"
        assert payload["gauge_certificate"]["pinned_dims"] == 0
