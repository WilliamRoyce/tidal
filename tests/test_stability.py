"""Regression tests for the canonical conversion-channel stability probe.

These tests pin the empirical Stage C truth-table results
(:file:`stage_c_truth_table.csv` at repo root): the unit-IC all-k Padé
probe at ``threshold=0.15`` and ``t_test=20`` (post-#341) must catch
the contamination samples found by the post-hoc t-independence audit
on the D1 amplify chain (HPC job 28520217), and must not over-reject
samples the audit confirmed perturbative.

Threshold history:
* Pre-#341: ``threshold=0.3`` caught 56/57 contamination — sample 391
  (γ_eff(t=20)=0.272) leaked through; caught only by the inline
  Hwang–Noh ``P_max > 0.5`` gate after the simulation.
* Post-#341: ``threshold=0.15`` catches 57/57 (sample 391 included);
  motivated by Phase 6 hi-res D1 amp MAP at γ_eff=0.281, where the
  γ ∈ [0.15, 0.30] regime produced exponential A(t) ∝ exp(0.27·t)
  and failed the t_end-independence cross-check (#340).

See ``docs/tex/stability_probe.tex`` for the architecture and #323
Stage C investigation for the empirical foundations.
"""

from __future__ import annotations

import csv
import math
import time
import warnings
from argparse import Namespace
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pytest

from tidal.measurement._stability import (
    PROBE_PROFILE_NAME,
    check_conversion_stability,
)

if TYPE_CHECKING:
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


REPO_ROOT = Path(__file__).resolve().parents[1]
TRUTH_TABLE = REPO_ROOT / "stage_c_truth_table.csv"
D1_SPEC = REPO_ROOT / "examples/data/torsion_gertsenshtein_nonminimal.json"
T1_SPEC = REPO_ROOT / "examples/data/dark_photon_plasma.json"

# Stage C IC parameters (matching the audited D1 amplify chain).
IC_K = 0.06283185307179587
KAPPA = 1.0
B0 = 0.01

# Known-residual samples not caught by the probe.  Empty post-#341:
# sample 391 (γ_eff=0.272) is now caught by the threshold=0.15 probe.
# Kept as a sentinel set so future regressions can be flagged here
# without changing test logic.
KNOWN_PROBE_RESIDUAL: set[str] = set()

# A complete parameter set for T1 (dark_photon_plasma): the spec has four
# free couplings plus kappa/B0, and an incomplete set makes the probe
# raise rather than answer.
T1_PARAMS: dict[str, float] = {
    "kappa": 1.0,
    "B0": 0.01,
    "mA2": 0.05,
    "deltam": 0.1,
    "xi": 0.3,
    "alpha3": 0.1,
}


def _load_d1() -> tuple[EquationSystem, GridInfo]:
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import load_equation_system

    spec = load_equation_system(D1_SPEC)
    grid = GridInfo(shape=(256,), bounds=((0.0, 100.0),), periodic=(True,))
    return spec, grid


def _load_t1() -> tuple[EquationSystem, GridInfo]:
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import load_equation_system

    spec = load_equation_system(T1_SPEC)
    grid = GridInfo(shape=(32,), bounds=((0.0, 50.0),), periodic=(True,))
    return spec, grid


# ---------------------------------------------------------------------
# Stage C truth-table regression
# ---------------------------------------------------------------------


@pytest.mark.slow
class TestStageCTruthTable:
    """Probe verdict must match Stage C sim verdict on every audited row."""

    @pytest.fixture(scope="class")
    def truth_rows(self) -> list[dict[str, str]]:
        if not TRUTH_TABLE.exists():
            pytest.skip(f"truth table missing: {TRUTH_TABLE}")
        with TRUTH_TABLE.open() as f:
            return [
                row
                for row in csv.DictReader(f)
                if row["sim_verdict"] in {"sim_perturbative", "sim_non_perturbative"}
            ]

    @pytest.fixture(scope="class")
    def d1(self) -> tuple[EquationSystem, GridInfo]:
        return _load_d1()

    def test_contamination_caught(
        self,
        truth_rows: list[dict[str, str]],
        d1: tuple[EquationSystem, GridInfo],
    ) -> None:
        """Every audited contamination sample must be rejected (≤1 residual)."""
        spec, grid = d1
        contaminated = [
            r for r in truth_rows if r["sim_verdict"] == "sim_non_perturbative"
        ]
        assert contaminated, "truth table missing sim_non_perturbative rows"

        misses: list[str] = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for row in contaminated:
                params = {
                    "kappa": KAPPA,
                    "B0": B0,
                    "alpha1": float(row["alpha1"]),
                    "alpha2": float(row["alpha2"]),
                    "alpha3": float(row["alpha3"]),
                    "delta1": float(row["delta1"]),
                }
                try:
                    res = check_conversion_stability(
                        spec,
                        grid,
                        params,
                        source="h_5",
                        target="a_1",
                        ic_wavevector=IC_K,
                    )
                    rejected = not res.stable
                except Exception:  # noqa: BLE001
                    rejected = True
                if not rejected:
                    misses.append(row["sample_index"])

        unexpected = set(misses) - KNOWN_PROBE_RESIDUAL
        assert not unexpected, (
            f"Probe missed contamination samples beyond known residual "
            f"{sorted(KNOWN_PROBE_RESIDUAL)}: {sorted(unexpected)}"
        )

    @pytest.mark.skip(
        reason="GH #467: the Stage-C sample set runs the nonminimal-class "
        "operator at points that now refuse honestly (near constraint-"
        "index breakdown, deflation-contract violation) — re-derive the "
        "truth table under the corrected operator + refusal semantics"
    )
    def test_perturbative_not_overrejected(
        self,
        truth_rows: list[dict[str, str]],
        d1: tuple[EquationSystem, GridInfo],
    ) -> None:
        """Audited perturbative samples must not be falsely rejected (≤2)."""
        spec, grid = d1
        perturbative = [r for r in truth_rows if r["sim_verdict"] == "sim_perturbative"]
        assert perturbative, "truth table missing sim_perturbative rows"

        false_rejects: list[str] = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for row in perturbative:
                params = {
                    "kappa": KAPPA,
                    "B0": B0,
                    "alpha1": float(row["alpha1"]),
                    "alpha2": float(row["alpha2"]),
                    "alpha3": float(row["alpha3"]),
                    "delta1": float(row["delta1"]),
                }
                try:
                    res = check_conversion_stability(
                        spec,
                        grid,
                        params,
                        source="h_5",
                        target="a_1",
                        ic_wavevector=IC_K,
                    )
                except Exception:  # noqa: BLE001
                    false_rejects.append(row["sample_index"])
                    continue
                if not res.stable:
                    false_rejects.append(row["sample_index"])

        # Allow up to 2 false rejections — the probe is conservative on
        # samples sitting near the perturbative boundary.
        assert len(false_rejects) <= 2, (
            f"Probe falsely rejected too many perturbative samples "
            f"({len(false_rejects)}/{len(perturbative)}): {false_rejects[:10]}"
        )


# ---------------------------------------------------------------------
# Decoupling sanity & known-tachyon
# ---------------------------------------------------------------------


@pytest.mark.slow
class TestProbeBehavior:
    @pytest.fixture(scope="class")
    def d1(self) -> tuple[EquationSystem, GridInfo]:
        return _load_d1()

    def test_d1_working_point_accepts(
        self, d1: tuple[EquationSystem, GridInfo]
    ) -> None:
        """D1 documented working point (α₁=0, α₂=−0.6, α₃=0.5) is stable."""
        spec, grid = d1
        params = {
            "kappa": KAPPA,
            "B0": B0,
            "alpha1": 0.0,
            "alpha2": -0.6,
            "alpha3": 0.5,
            "delta1": 0.0,
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = check_conversion_stability(
                spec,
                grid,
                params,
                source="h_5",
                target="a_1",
                ic_wavevector=IC_K,
            )
        assert res.stable, f"working point should be stable: {res.message}"
        assert res.profile_name == PROBE_PROFILE_NAME

    def test_known_tachyon_rejected(self, d1: tuple[EquationSystem, GridInfo]) -> None:
        """Stage C sample 2808 (v1 γ_eff≈18.5, sim non-perturbative) → reject."""
        spec, grid = d1
        # Stage C truth-table row 2808: a strong tachyonic point that
        # the v1 probe catches at γ_eff ≈ 18.5; sim verdict
        # ``sim_non_perturbative``.  Concrete known-tachyon regression
        # for the canonical probe.
        params = {
            "kappa": KAPPA,
            "B0": B0,
            "alpha1": -0.292296681808213,
            "alpha2": -1.88963964322586,
            "alpha3": 0.578692478437655,
            "delta1": -0.105453324639565,
        }
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = check_conversion_stability(
                spec,
                grid,
                params,
                source="h_5",
                target="a_1",
                ic_wavevector=IC_K,
            )
        assert not res.stable, (
            f"known-tachyon Stage C sample 2808 should be rejected: {res.message}"
        )
        assert res.max_excess > 0.3


# ---------------------------------------------------------------------
# Performance
# ---------------------------------------------------------------------


@pytest.mark.slow
class TestProbePerformance:
    @pytest.fixture(scope="class")
    def t1(self) -> tuple[EquationSystem, GridInfo]:
        return _load_t1()

    def test_probe_perf(self, t1: tuple[EquationSystem, GridInfo]) -> None:
        """Probe wall stays within the recorded budget (T1, N=32).

        Guard only -- the asserted budget below is the operative number.
        Wall-clock here is load-sensitive, so this catches gross regressions
        rather than pinning a precise figure.
        """
        spec, grid = t1

        rng = np.random.default_rng(0)
        # T1 dark_photon_plasma has 4 free coupling parameters.
        timings: list[float] = []
        n_iter = 1000
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for _ in range(n_iter):
                params = {
                    "kappa": 1.0,
                    "B0": 0.01,
                    "mA2": float(rng.uniform(0.001, 1.0)),
                    "deltam": float(rng.uniform(-0.5, 0.5)),
                    "xi": float(rng.uniform(0.05, 1.0)),
                    "alpha3": float(rng.uniform(-0.5, 0.5)),
                }
                t0 = time.perf_counter()
                try:
                    check_conversion_stability(
                        spec,
                        grid,
                        params,
                        source="h_5",
                        target="a_1",
                        ic_wavevector=2.0 * np.pi / 50.0,
                    )
                except Exception:  # noqa: BLE001, S112
                    continue
                timings.append(time.perf_counter() - t0)

        assert timings, "no valid probe runs"
        median_ms = float(np.median(timings)) * 1e3
        # Empirical baseline on T1 dark_photon_plasma (N=32, 4 free
        # couplings) under the canonical unit-IC all-k probe at
        # t_test=20.  Two coupled wins (issue #327):
        #
        #  * Phases 1-2 (commit 350b53b et al.): per-process spec cache
        #    + parameter-independent structural cache for the probe
        #    saved ~2 ms/call.
        #  * Phase 3 (commit 49c9d447): ``np.einsum`` → ``np.matmul``
        #    refactor inside ``_build_evolution_matrices``.  Batched
        #    3-D matmul dispatches to BLAS gemm; the 3-way einsums
        #    that replicated ``(U @ K) @ V`` were ~70× slower without
        #    ``optimize=True``.  Saved ~12 ms/call.
        #
        # Combined, as measured at the time: ≈ 27.6 ms → ≈ 13 ms (2.13×).
        # (Historical -- superseded by the post-#341 figure below.)
        #
        # Post-#341 (threshold 0.30 → 0.15): the spectral-radius
        # prefilter cutoff is ``threshold * t_test`` so dropping the
        # threshold by 2× admits roughly twice as many modes through
        # the ``expm`` call — observed median rose from ~13 ms to
        # ~62 ms (T1 N=32). The probe is still fast in absolute terms
        # (∼60 ms vs typical sample-time of seconds) and the wall-time
        # cost is the unavoidable price of tightening the threshold
        # without weakening probe correctness. Budget set ∼25 % above
        # the observed median; a 50 % rise here would indicate a real
        # regression beyond the expected #341 cost.
        assert median_ms <= 80.0, f"median probe wall {median_ms:.2f} ms > 80 ms"


# ---------------------------------------------------------------------
# Shared probe stage and the sweep/sample policy (GH #454)
# ---------------------------------------------------------------------


class TestProbeForRunSharedStage:
    """``probe_for_run`` is the ONE probe stage both entry points use.

    ``tidal sweep`` and ``tidal sample`` each carried a copy of this
    preamble, and the copies drifted into two different tachyonic
    policies for four months (GH #454). These tests pin the properties
    that make one shared stage worth having.
    """

    @staticmethod
    def _args(**overrides: object) -> Namespace:
        base: dict[str, object] = {
            "param": [],
            "grid_shape": 32,
            "bounds": "0:50",
            "ic_wavevector": None,
        }
        base.update(overrides)
        return Namespace(**base)

    def test_not_applicable_without_conversion_measurement(self) -> None:
        from tidal.measurement._stability import probe_for_run

        result, meta = probe_for_run(
            T1_SPEC,
            self._args(),
            {},
            source=("a_1",),
            target=("a_2",),
            measurements={"energy"},
        )
        assert result is None
        assert meta == {}

    def test_not_applicable_without_source_and_target(self) -> None:
        from tidal.measurement._stability import probe_for_run

        result, meta = probe_for_run(
            T1_SPEC,
            self._args(),
            {},
            source=None,
            target=None,
            measurements={"conversion"},
        )
        assert result is None
        assert meta == {}

    def test_metadata_key_set_is_fixed(self) -> None:
        """Schema parity: a sweep row and a chain sample must be comparable.

        The absence of this invariant is what let sweep record three
        columns (one under a different name) while inference recorded
        five, so nothing could be joined across the two.
        """
        from tidal.measurement._stability import (
            PROBE_METADATA_KEYS,
            probe_for_run,
        )

        _result, meta = probe_for_run(
            T1_SPEC,
            self._args(),
            T1_PARAMS,
            source=("a_1",),
            target=("a_2",),
            measurements={"conversion"},
        )
        assert set(meta) == set(PROBE_METADATA_KEYS)
        # k_tachyonic is NaN, never absent, when nothing is tachyonic —
        # so the column exists in every row.
        assert "k_tachyonic" in meta

    def test_missing_bounds_no_longer_silently_disables_the_probe(self) -> None:
        """GH #454: sweeps without ``--bounds`` ran no probe at all.

        The old sweep code passed ``bounds=None`` into ``GridInfo``,
        which raises ``TypeError``, and a bare ``except Exception``
        swallowed it at DEBUG. Since ``--bounds`` defaults to None in
        every subparser, those sweeps got neither a verdict nor the
        diagnostic columns.
        """
        from tidal.measurement._stability import probe_for_run

        result, meta = probe_for_run(
            T1_SPEC,
            self._args(bounds=None),
            T1_PARAMS,
            source=("a_1",),
            target=("a_2",),
            measurements={"conversion"},
        )
        assert result is not None, "probe silently skipped without --bounds"
        assert meta

    def test_probe_grid_matches_the_grid_the_simulation_will_use(self) -> None:
        """GH #479: the probe must describe the system that gets evolved.

        Three private fallbacks (256, 64, 256) stood in for this and all
        three disagreed with the simulation's own default of 64 points on
        (0, 10), giving the probe a Nyquist 2.5-5x BELOW the solver's —
        so modes the solver evolves were never examined. That is a false
        negative in a probe whose design principle is never to risk one.
        """
        from tidal.measurement._run_stages import parse_bounds, parse_grid_shape
        from tidal.symbolic._spec_cache import load_spec_cached

        captured: dict[str, GridInfo] = {}
        import tidal.measurement._stability as stability_mod

        original = stability_mod.check_conversion_stability

        def _capture(
            _spec: object, grid: GridInfo, *_a: object, **_k: object
        ) -> object:
            captured["grid"] = grid
            return original(_spec, grid, *_a, **_k)  # type: ignore[arg-type]

        stability_mod.check_conversion_stability = _capture  # type: ignore[assignment]
        try:
            stability_mod.probe_for_run(
                T1_SPEC,
                self._args(grid_shape=None, bounds=None),
                T1_PARAMS,
                source=("a_1",),
                target=("a_2",),
                measurements={"conversion"},
            )
        finally:
            stability_mod.check_conversion_stability = original  # type: ignore[assignment]

        spec = load_spec_cached(T1_SPEC)
        sim_shape = parse_grid_shape(None, spec.spatial_dimension)
        sim_bounds = parse_bounds(None, spec.spatial_dimension)

        probe_grid = captured["grid"]
        assert probe_grid.shape[0] == sim_shape[0]
        assert probe_grid.bounds[0] == sim_bounds[0]

        # The property that actually matters: the probe must not examine a
        # narrower band of k than the solver will evolve.
        probe_k_max = (
            math.pi
            * probe_grid.shape[0]
            / (probe_grid.bounds[0][1] - probe_grid.bounds[0][0])
        )
        sim_k_max = math.pi * sim_shape[0] / (sim_bounds[0][1] - sim_bounds[0][0])
        assert probe_k_max >= sim_k_max, (
            f"probe Nyquist {probe_k_max:.2f} < solver Nyquist {sim_k_max:.2f}: "
            "tachyonic modes in the gap would be invisible to the probe"
        )

    def test_unavailable_meta_hook_is_used(self) -> None:
        """The GH #421 posdep-kinetic refusal keeps its warn-once marker."""
        import tidal.measurement._stability as stability_mod

        calls: list[BaseException] = []

        def _marker(exc: BaseException) -> dict[str, object]:
            calls.append(exc)
            return {"stability_profile": "unavailable-posdep-kinetic"}

        def _boom(*_a: object, **_k: object) -> object:
            msg = "position-dependent kinetic"
            raise NotImplementedError(msg)

        original = stability_mod.check_conversion_stability
        stability_mod.check_conversion_stability = _boom  # type: ignore[assignment]
        try:
            result, meta = stability_mod.probe_for_run(
                T1_SPEC,
                self._args(),
                {},
                source=("a_1",),
                target=("a_2",),
                measurements={"conversion"},
                unavailable_meta=_marker,
            )
        finally:
            stability_mod.check_conversion_stability = original  # type: ignore[assignment]

        assert result is None
        assert meta == {"stability_profile": "unavailable-posdep-kinetic"}
        assert len(calls) == 1


class TestSweepTachyonicPolicy:
    """``tidal sweep`` records the probe verdict; it does not block (GH #454).

    Until v0.49.5 the sweep path returned early with
    ``run_status="tachyonic"`` and never simulated, while ``tidal
    sample`` recorded the same verdict as metadata and continued. That
    was the April-2026 policy left behind by the 2026-05-10 v3 change,
    not a decision.
    """

    @staticmethod
    def _unstable() -> object:
        from tidal.measurement._stability import ConversionStabilityResult

        return ConversionStabilityResult(
            stable=False,
            max_excess=0.42,
            k_tachyonic=0.31,
            n_tachyonic_modes=7,
            message="tachyonic at k=0.31",
        )

    def _patch(
        self,
        monkeypatch: pytest.MonkeyPatch,
        simulated: list[bool],
    ) -> None:
        """Patch the shared stages, not the caller.

        ``_run_single`` delegates the sequence to ``run_point`` (GH #480),
        so the stages are patched where they now live.  The behavior under
        test is unchanged — that the probe verdict is recorded and the
        point simulated anyway — only the seam moved.
        """
        import tidal.measurement._run_stages as stages
        from tidal.measurement._stability import probe_metadata

        unstable = self._unstable()

        def _probe(*_a: object, **_k: object) -> tuple[object, dict[str, object]]:
            return unstable, probe_metadata(unstable)  # type: ignore[arg-type]

        def _simulate(*_a: object, **_k: object) -> tuple[int, float, None]:
            simulated.append(True)
            return 0, 1.25, None

        def _measure(*_a: object, **_k: object) -> dict[str, object]:
            return {"P_max": 0.017}

        monkeypatch.setattr(stages, "probe_for_run", _probe)
        monkeypatch.setattr(stages, "simulate_run", _simulate)
        monkeypatch.setattr(stages, "measure_run", _measure)

    @staticmethod
    def _args() -> Namespace:
        return Namespace(
            param=[],
            grid_shape=32,
            bounds="0:50",
            ic_wavevector=None,
        )

    def test_default_simulates_and_records_the_verdict(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
    ) -> None:
        import tidal.cli._sweep as sweep_mod

        simulated: list[bool] = []
        self._patch(monkeypatch, simulated)

        row = sweep_mod._run_single(
            self._args(),
            T1_SPEC,
            {},
            tmp_path,
            {"conversion"},
            ("a_1",),
            ("a_2",),
            0.99,
        )

        assert simulated == [True], "unstable point was not simulated"
        assert row["run_status"] == "success"
        assert row["P_max"] == 0.017
        assert row["wall_time_s"] > 0.0
        # The verdict is recorded, not acted on.
        assert row["tachyonic_excess"] == 0.42
        assert row["n_tachyonic_modes"] == 7
        assert row["k_tachyonic"] == 0.31

    def test_there_is_no_configuration_that_blocks(self) -> None:
        """The probe verdict is unconditionally a diagnostic (v0.50.0).

        ``--gated`` reproduced the pre-v0.49.5 blocking row for archived
        sweeps; it was removed once rejection on tachyonic growth was
        abandoned as policy, so the abandoned behavior is now structurally
        unreachable rather than merely off by default.
        """
        from tidal.cli import _build_parser

        with pytest.raises(SystemExit):
            _build_parser().parse_args(
                [
                    "sweep",
                    "x.json",
                    "--sweep",
                    "p=0:1:2",
                    "--measure",
                    "conversion",
                    "--gated",
                ]
            )
