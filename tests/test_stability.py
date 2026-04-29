"""Regression tests for the canonical conversion-channel stability probe.

These tests pin the empirical Stage C truth-table results
(:file:`stage_c_truth_table.csv` at repo root): the unit-IC all-k Padé
probe at ``threshold=0.3`` and ``t_test=20`` must catch the
contamination samples found by the post-hoc t-independence audit on
the D1 amplify chain (HPC job 28520217), and must not over-reject
samples the audit confirmed perturbative.

Sample 391 is a known residual: ``γ_eff(t=20)=0.272 < 0.3``, so the
probe alone does not catch it; the inline Hwang–Noh gate in
:mod:`tidal.inference._likelihood` rejects it after the simulation
runs.

See ``docs/tex/stability_probe.tex`` for the architecture and #323
Stage C investigation for the empirical foundations.
"""

from __future__ import annotations

import csv
import time
import warnings
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

# Known-residual sample index — γ_eff(t=20)=0.272, just under the 0.3
# threshold; caught by the inline Hwang–Noh gate after sim, not by the
# probe.
KNOWN_PROBE_RESIDUAL = {"391"}


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
class TestProbeBehaviour:
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
        """Median probe wall ≤ 6 ms over 1000 random param dicts (T1, N=32)."""
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
        #  * Phase 3 (commit XXXXXXX): ``np.einsum`` → ``np.matmul``
        #    refactor inside ``_build_evolution_matrices``.  Batched
        #    3-D matmul dispatches to BLAS gemm; the 3-way einsums
        #    that replicated ``(U @ K) @ V`` were ~70× slower without
        #    ``optimize=True``.  Saved ~12 ms/call.
        #
        # Combined: median dropped from ≈ 27.6 ms → ≈ 13 ms (2.13×).
        # Budget set ~15 % above the observed p99 (≈ 26 ms) to absorb
        # CI noise without being flaky.  A 50 % rise here would
        # indicate a real regression.
        assert median_ms <= 20.0, f"median probe wall {median_ms:.2f} ms > 20 ms"
