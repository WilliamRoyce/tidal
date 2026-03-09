"""Tests for ensemble/replicate runs in ``tidal sweep``.

Tests cover:
- Replicate expansion logic
- Parameter noise injection
- IC perturbation
- SweepResults aggregation methods
- CLI argument parsing
- TOML config parsing
- Seed determinism
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from tidal.cli._sweep import (
    _apply_param_noise_to_plan,
    _expand_replicates,
    _parse_param_noise,
)
from tidal.measurement._sweep_results import SweepResults

# ------------------------------------------------------------------
# _parse_param_noise
# ------------------------------------------------------------------


class TestParseParamNoise:
    """Test ``--param-noise`` CLI argument parsing."""

    def test_none_returns_empty(self) -> None:
        assert _parse_param_noise(None) == {}

    def test_empty_list_returns_empty(self) -> None:
        assert _parse_param_noise([]) == {}

    def test_single_param(self) -> None:
        result = _parse_param_noise(["B0=0.01"])
        assert result == {"B0": 0.01}

    def test_multiple_params(self) -> None:
        result = _parse_param_noise(["B0=0.01", "m2=0.05"])
        assert result == {"B0": 0.01, "m2": 0.05}

    def test_invalid_no_equals(self) -> None:
        with pytest.raises(ValueError, match="PARAM=SCALE"):
            _parse_param_noise(["B0_0.01"])


# ------------------------------------------------------------------
# _expand_replicates
# ------------------------------------------------------------------


def _make_plan(
    swept: dict[str, float] | None = None,
    overrides: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Create a minimal run plan for testing."""
    swept = swept or {"g0": 1.0}
    overrides = overrides or {"g0": 1.0}
    return {
        "index": 0,
        "param_overrides": dict(overrides),
        "swept_vals": dict(swept),
        "subdir": "g0_1",
        "run_dir": Path("/tmp/test/g0_1"),
        "grid_override": None,
    }


class TestExpandReplicates:
    """Test replicate expansion of run plans."""

    def test_basic_expansion(self) -> None:
        plans = [_make_plan()]
        expanded = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={}
        )
        assert len(expanded) == 3
        for i, rp in enumerate(expanded):
            assert rp["replicate"] == i
            assert rp["seed"] == 42 + i
            assert rp["subdir"] == f"g0_1/rep_{i}"
            assert rp["run_dir"] == Path(f"/tmp/test/g0_1/rep_{i}")
            assert rp["index"] == i

    def test_multi_point_expansion(self) -> None:
        plans = [
            _make_plan({"g0": 1.0}, {"g0": 1.0}),
            _make_plan({"g0": 2.0}, {"g0": 2.0}),
        ]
        plans[1]["subdir"] = "g0_2"
        plans[1]["run_dir"] = Path("/tmp/test/g0_2")
        plans[1]["index"] = 1
        expanded = _expand_replicates(
            plans, n_replicates=2, base_seed=0, param_noise={}
        )
        assert len(expanded) == 4
        # First two should be replicates of point 1
        assert expanded[0]["swept_vals"]["g0"] == 1.0
        assert expanded[0]["replicate"] == 0
        assert expanded[1]["swept_vals"]["g0"] == 1.0
        assert expanded[1]["replicate"] == 1
        # Next two should be replicates of point 2
        assert expanded[2]["swept_vals"]["g0"] == 2.0
        assert expanded[2]["replicate"] == 0
        assert expanded[3]["swept_vals"]["g0"] == 2.0
        assert expanded[3]["replicate"] == 1

    def test_no_mutation_of_originals(self) -> None:
        plans = [_make_plan()]
        original_swept = dict(plans[0]["swept_vals"])
        _expand_replicates(plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1})
        assert plans[0]["swept_vals"] == original_swept

    def test_indices_sequential(self) -> None:
        plans = [_make_plan({"g0": v}) for v in [1.0, 2.0, 3.0]]
        for i, p in enumerate(plans):
            p["index"] = i
            p["subdir"] = f"g0_{int(p['swept_vals']['g0'])}"
            p["run_dir"] = Path(f"/tmp/test/{p['subdir']}")
        expanded = _expand_replicates(
            plans, n_replicates=2, base_seed=10, param_noise={}
        )
        indices = [rp["index"] for rp in expanded]
        assert indices == list(range(6))


# ------------------------------------------------------------------
# Parameter noise
# ------------------------------------------------------------------


class TestParamNoise:
    """Test parameter noise injection in replicate plans."""

    def test_noise_applied(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        expanded = _expand_replicates(
            plans, n_replicates=5, base_seed=42, param_noise={"g0": 0.1}
        )
        # Each replicate should have a different g0 value
        g0_vals = [rp["swept_vals"]["g0"] for rp in expanded]
        assert len(set(g0_vals)) == 5  # All distinct
        # All should be close to 1.0 but not exactly 1.0
        for v in g0_vals:
            assert v != 1.0  # Noise was applied
            assert abs(v - 1.0) < 1.0  # Reasonable noise scale

    def test_nominal_vals_stored(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        expanded = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1}
        )
        for rp in expanded:
            assert rp["nominal_vals"] == {"g0": 1.0}

    def test_noise_deterministic(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        exp1 = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1}
        )
        exp2 = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1}
        )
        for a, b in zip(exp1, exp2, strict=True):
            assert a["swept_vals"]["g0"] == b["swept_vals"]["g0"]

    def test_noise_not_applied_when_empty(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        expanded = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={}
        )
        for rp in expanded:
            assert rp["swept_vals"]["g0"] == 1.0
            assert "nominal_vals" not in rp

    def test_apply_param_noise_to_plan(self) -> None:
        plan: dict[str, Any] = {
            "swept_vals": {"g0": 1.0, "m2": 2.0},
            "param_overrides": {"g0": 1.0, "m2": 2.0},
        }
        _apply_param_noise_to_plan(plan, {"g0": 0.1}, seed=42)
        assert plan["swept_vals"]["g0"] != 1.0
        assert plan["swept_vals"]["m2"] == 2.0  # Not noised
        assert plan["nominal_vals"]["g0"] == 1.0


# ------------------------------------------------------------------
# SweepResults replicate properties
# ------------------------------------------------------------------


def _make_results(
    n_points: int = 3,
    n_reps: int = 3,
    *,
    with_noise: bool = False,
) -> SweepResults:
    """Create synthetic SweepResults with replicate rows."""
    rng = np.random.default_rng(42)
    rows: list[dict[str, Any]] = []
    param_values = [1.0, 2.0, 3.0][:n_points]

    for p_val in param_values:
        for rep in range(n_reps):
            noise = rng.standard_normal() * 0.1 if with_noise else 0.0
            rows.append(
                {
                    "g0": p_val,
                    "P_max": p_val * 0.5 + noise,
                    "max_energy_error": 1e-5 + abs(noise) * 1e-5,
                    "wall_time_s": 1.0 + noise * 0.1,
                    "run_status": "success",
                    "error_message": None,
                    "solver_exit_code": 0,
                    "replicate": rep,
                    "seed": 42 + rep,
                }
            )

    return SweepResults(
        swept_params={"g0": param_values},
        fixed_params={"m2": 1.0},
        sim_settings={
            "grid_shape": 64,
            "t_end": 10.0,
            "dt": "auto",
            "scheme": "auto",
            "bc": "periodic",
        },
        rows=rows,
        run_dirs=[],
        spec_path="test.json",
        measurements=["peak_conversion"],
        metadata={"n_replicates": n_reps, "base_seed": 42},
    )


def _make_results_no_replicates() -> SweepResults:
    """Create SweepResults without replicates (backward compat)."""
    rows = [
        {"g0": 1.0, "P_max": 0.5, "run_status": "success"},
        {"g0": 2.0, "P_max": 1.0, "run_status": "success"},
        {"g0": 3.0, "P_max": 1.5, "run_status": "success"},
    ]
    return SweepResults(
        swept_params={"g0": [1.0, 2.0, 3.0]},
        fixed_params={},
        sim_settings={"grid_shape": 64},
        rows=rows,
        run_dirs=[],
        spec_path="test.json",
        measurements=["peak_conversion"],
    )


class TestSweepResultsReplicates:
    """Test replicate-related properties and methods on SweepResults."""

    def test_has_replicates_true(self) -> None:
        r = _make_results()
        assert r.has_replicates is True

    def test_has_replicates_false(self) -> None:
        r = _make_results_no_replicates()
        assert r.has_replicates is False

    def test_n_replicates(self) -> None:
        r = _make_results(n_reps=5)
        assert r.n_replicates == 5

    def test_n_replicates_no_replicates(self) -> None:
        r = _make_results_no_replicates()
        assert r.n_replicates == 1

    def test_group_by_point(self) -> None:
        r = _make_results(n_points=2, n_reps=3)
        groups = r.group_by_point()
        assert len(groups) == 2
        for rep_rows in groups.values():
            assert len(rep_rows) == 3
            # All rows in group should have same g0
            g0_vals = {row["g0"] for row in rep_rows}
            assert len(g0_vals) == 1


class TestSweepResultsAggregate:
    """Test the aggregate() method for computing ensemble statistics."""

    def test_aggregate_basic(self) -> None:
        r = _make_results(n_points=3, n_reps=3)
        agg = r.aggregate()
        assert len(agg.rows) == 3  # One per parameter point
        for row in agg.rows:
            assert "P_max_mean" in row
            assert "P_max_std" in row
            assert "P_max_sem" in row
            assert "P_max_ci95_lo" in row
            assert "P_max_ci95_hi" in row
            assert "P_max_median" in row
            assert "P_max_q25" in row
            assert "P_max_q75" in row
            assert "P_max_n" in row
            # Plain metric = mean (backward compat)
            assert row["P_max"] == row["P_max_mean"]

    def test_aggregate_correct_mean(self) -> None:
        # Use with_noise=False so mean = exact value
        r = _make_results(n_points=2, n_reps=4, with_noise=False)
        agg = r.aggregate()
        # g0=1.0 -> P_max = 0.5 for all reps
        row1 = next(row for row in agg.rows if row["g0"] == 1.0)
        assert row1["P_max_mean"] == pytest.approx(0.5)
        assert row1["P_max_std"] == pytest.approx(0.0)
        assert row1["P_max_sem"] == pytest.approx(0.0)

    def test_aggregate_noisy(self) -> None:
        r = _make_results(n_points=2, n_reps=10, with_noise=True)
        agg = r.aggregate()
        for row in agg.rows:
            assert row["P_max_std"] > 0
            assert row["P_max_sem"] < row["P_max_std"]
            assert row["P_max_ci95_lo"] < row["P_max_mean"]
            assert row["P_max_ci95_hi"] > row["P_max_mean"]
            assert row["P_max_n"] == 10

    def test_aggregate_no_replicates_returns_self(self) -> None:
        r = _make_results_no_replicates()
        agg = r.aggregate()
        assert agg is r  # Same object

    def test_aggregate_preserves_params(self) -> None:
        r = _make_results()
        agg = r.aggregate()
        g0_values = sorted(row["g0"] for row in agg.rows)
        assert g0_values == [1.0, 2.0, 3.0]

    def test_aggregate_with_failed_rows(self) -> None:
        """Failed rows should be excluded from metric aggregation."""
        r = _make_results(n_points=1, n_reps=3, with_noise=False)
        # Mark one row as failed
        r.rows[1]["run_status"] = "diverged"
        r.rows[1]["P_max"] = None
        agg = r.aggregate()
        assert agg.rows[0]["P_max_n"] == 2  # Only 2 successful

    def test_aggregate_single_replicate(self) -> None:
        r = _make_results(n_points=2, n_reps=1)
        agg = r.aggregate()
        for row in agg.rows:
            assert row["P_max_std"] == 0.0
            assert row["P_max_ci95_lo"] == row["P_max_mean"]

    def test_aggregate_specific_metrics(self) -> None:
        r = _make_results()
        agg = r.aggregate(metrics=["P_max"])
        for row in agg.rows:
            assert "P_max_mean" in row
            # max_energy_error should NOT be aggregated
            assert "max_energy_error_mean" not in row

    def test_aggregate_with_nominal_vals(self) -> None:
        """When param noise is used, grouping should use nominal values."""
        rows: list[dict[str, Any]] = [
            {
                "g0": 1.0 + rep * 0.01,  # Noised values
                "g0_nominal": 1.0,  # Nominal
                "P_max": 0.5 + rep * 0.1,
                "run_status": "success",
                "replicate": rep,
                "seed": rep,
            }
            for rep in range(3)
        ]
        r = SweepResults(
            swept_params={"g0": [1.0]},
            fixed_params={},
            sim_settings={},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )
        groups = r.group_by_point()
        assert len(groups) == 1  # All grouped under nominal=1.0


# ------------------------------------------------------------------
# IC perturbation
# ------------------------------------------------------------------


class TestICPerturbation:
    """Test _apply_ic_perturbation function."""

    def test_no_perturbation_when_none(self) -> None:
        from argparse import Namespace

        from tidal.cli._simulate import _apply_ic_perturbation

        args = Namespace(ic_perturbation=None, ic_amplitude=1.0)
        y0 = np.ones(64)
        # ic_perturbation=None → no-op, doesn't need spec/grid
        result = _apply_ic_perturbation(y0, args, None, None)  # type: ignore[arg-type]
        np.testing.assert_array_equal(result, y0)

    def test_no_perturbation_when_zero(self) -> None:
        from argparse import Namespace

        from tidal.cli._simulate import _apply_ic_perturbation

        args = Namespace(ic_perturbation=0.0, ic_amplitude=1.0)
        y0 = np.ones(64)
        result = _apply_ic_perturbation(y0, args, None, None)  # type: ignore[arg-type]
        np.testing.assert_array_equal(result, y0)


# ------------------------------------------------------------------
# TOML config parsing
# ------------------------------------------------------------------


class TestSweepConfigReplicates:
    """Test TOML config parsing for replicate settings."""

    def test_parse_execution_replicates(self, tmp_path: Path) -> None:
        toml_content = """
spec = "test.json"

[sweep.g0]
start = 0.1
stop = 1.0
count = 5

[measurement]
types = ["peak_conversion"]

[execution]
n_replicates = 10
base_seed = 123
"""
        toml_file = tmp_path / "sweep.toml"
        toml_file.write_text(toml_content)
        # Create the spec file
        (tmp_path / "test.json").write_text("{}")

        from tidal.cli._sweep_config import load_sweep_config

        config = load_sweep_config(toml_file)
        assert config.n_replicates == 10
        assert config.base_seed == 123

    def test_parse_ensemble_section(self, tmp_path: Path) -> None:
        toml_content = """
spec = "test.json"

[sweep.g0]
start = 0.1
stop = 1.0
count = 5

[measurement]
types = ["peak_conversion"]

[ensemble]
ic_perturbation = 0.05
param_noise = { g0 = 0.01, m2 = 0.05 }
"""
        toml_file = tmp_path / "sweep.toml"
        toml_file.write_text(toml_content)
        (tmp_path / "test.json").write_text("{}")

        from tidal.cli._sweep_config import load_sweep_config

        config = load_sweep_config(toml_file)
        assert config.ic_perturbation == 0.05
        assert config.param_noise == {"g0": 0.01, "m2": 0.05}

    def test_defaults_without_ensemble(self, tmp_path: Path) -> None:
        toml_content = """
spec = "test.json"

[sweep.g0]
start = 0.1
stop = 1.0
count = 5

[measurement]
types = ["peak_conversion"]
"""
        toml_file = tmp_path / "sweep.toml"
        toml_file.write_text(toml_content)
        (tmp_path / "test.json").write_text("{}")

        from tidal.cli._sweep_config import load_sweep_config

        config = load_sweep_config(toml_file)
        assert config.n_replicates == 1
        assert config.base_seed == 42
        assert config.ic_perturbation is None
        assert config.param_noise is None


# ------------------------------------------------------------------
# CSV serialization with replicates
# ------------------------------------------------------------------


class TestReplicateSerialization:
    """Test that replicate data round-trips through CSV/JSON."""

    def test_csv_includes_replicate_column(self) -> None:
        r = _make_results(n_points=2, n_reps=2)
        csv_str = r.to_csv()
        assert "replicate" in csv_str
        assert "seed" in csv_str

    def test_json_includes_replicate(self) -> None:
        import json

        r = _make_results(n_points=2, n_reps=2)
        json_str = r.to_json()
        data = json.loads(json_str)
        assert "replicate" in data["columns"]
        assert all("replicate" in row for row in data["rows"])

    def test_aggregated_csv_no_replicate(self) -> None:
        r = _make_results(n_points=2, n_reps=3)
        agg = r.aggregate()
        csv_str = agg.to_csv()
        assert "P_max_mean" in csv_str
        assert "P_max_std" in csv_str
        # Aggregated rows should not have replicate column
        assert "replicate" not in csv_str


# ------------------------------------------------------------------
# Seed determinism
# ------------------------------------------------------------------


class TestSeedDeterminism:
    """Test that the same base_seed produces identical ensembles."""

    def test_expand_deterministic(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        noise = {"g0": 0.5}

        exp_a = _expand_replicates(plans, 5, base_seed=99, param_noise=noise)
        exp_b = _expand_replicates(plans, 5, base_seed=99, param_noise=noise)

        for a, b in zip(exp_a, exp_b, strict=True):
            assert a["swept_vals"]["g0"] == b["swept_vals"]["g0"]
            assert a["seed"] == b["seed"]

    def test_different_seeds_different_noise(self) -> None:
        plans = [_make_plan({"g0": 1.0}, {"g0": 1.0})]
        noise = {"g0": 0.5}

        exp_a = _expand_replicates(plans, 3, base_seed=1, param_noise=noise)
        exp_b = _expand_replicates(plans, 3, base_seed=100, param_noise=noise)

        vals_a = [rp["swept_vals"]["g0"] for rp in exp_a]
        vals_b = [rp["swept_vals"]["g0"] for rp in exp_b]
        assert vals_a != vals_b


# ------------------------------------------------------------------
# Bug fix A1: Param noise RNG independence per parameter
# ------------------------------------------------------------------


class TestParamNoiseRNGIndependence:
    """Adding/removing a noised parameter should not change other draws."""

    def test_independent_per_param(self) -> None:
        """g0 noise should be identical whether m2 noise is present or not."""
        plans = [_make_plan({"g0": 1.0, "m2": 2.0}, {"g0": 1.0, "m2": 2.0})]
        plans[0]["subdir"] = "g0_1_m2_2"
        plans[0]["run_dir"] = Path("/tmp/test/g0_1_m2_2")

        # With only g0 noise
        exp_g0_only = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1}
        )
        g0_vals_only = [rp["swept_vals"]["g0"] for rp in exp_g0_only]

        # With both g0 and m2 noise
        exp_both = _expand_replicates(
            plans, n_replicates=3, base_seed=42, param_noise={"g0": 0.1, "m2": 0.05}
        )
        g0_vals_both = [rp["swept_vals"]["g0"] for rp in exp_both]

        # g0 draws should be identical regardless of m2 presence
        for a, b in zip(g0_vals_only, g0_vals_both, strict=True):
            assert a == b


# ------------------------------------------------------------------
# Bug fix A3: CSV column ordering
# ------------------------------------------------------------------


class TestCSVColumnOrdering:
    """Replicate/seed columns should appear before metrics."""

    def test_replicate_before_metrics(self) -> None:
        from tidal.cli._sweep import _build_row

        row = _build_row(
            swept_vals={"g0": 1.0},
            fixed_params={"m2": 2.0},
            sim_settings={"grid_shape": 64},
            metrics={"P_max": 0.5, "wall_time_s": 1.0, "run_status": "success"},
            grid_override=None,
            replicate=0,
            seed=42,
        )
        keys = list(row.keys())
        # replicate and seed should come before P_max
        assert keys.index("replicate") < keys.index("P_max")
        assert keys.index("seed") < keys.index("P_max")
        # swept params should come first
        assert keys.index("g0") < keys.index("replicate")

    def test_nominal_vals_after_swept(self) -> None:
        from tidal.cli._sweep import _build_row

        row = _build_row(
            swept_vals={"g0": 1.01},
            fixed_params={},
            sim_settings={},
            metrics={"P_max": 0.5},
            grid_override=None,
            replicate=0,
            seed=42,
            nominal_vals={"g0": 1.0},
        )
        keys = list(row.keys())
        assert keys.index("g0") < keys.index("g0_nominal")
        assert keys.index("g0_nominal") < keys.index("replicate")


# ------------------------------------------------------------------
# Bootstrap CI
# ------------------------------------------------------------------


class TestBootstrapCI:
    """Test non-parametric bootstrap confidence intervals."""

    def test_bootstrap_basic(self) -> None:
        r = _make_results(n_points=2, n_reps=10, with_noise=True)
        boot = r.bootstrap_ci("P_max", n_bootstrap=1000, seed=42)
        assert len(boot.rows) == 2
        for row in boot.rows:
            assert "P_max_boot_lo" in row
            assert "P_max_boot_hi" in row
            assert row["P_max_boot_lo"] < row["P_max_boot_hi"]

    def test_bootstrap_no_replicates_returns_self(self) -> None:
        r = _make_results_no_replicates()
        boot = r.bootstrap_ci("P_max")
        assert boot is r

    def test_bootstrap_deterministic(self) -> None:
        r = _make_results(n_points=1, n_reps=10, with_noise=True)
        b1 = r.bootstrap_ci("P_max", n_bootstrap=500, seed=123)
        b2 = r.bootstrap_ci("P_max", n_bootstrap=500, seed=123)
        assert b1.rows[0]["P_max_boot_lo"] == b2.rows[0]["P_max_boot_lo"]
        assert b1.rows[0]["P_max_boot_hi"] == b2.rows[0]["P_max_boot_hi"]

    def test_bootstrap_single_replicate(self) -> None:
        r = _make_results(n_points=1, n_reps=1)
        boot = r.bootstrap_ci("P_max")
        assert boot.rows[0]["P_max_boot_lo"] == boot.rows[0]["P_max_boot_hi"]

    def test_bootstrap_ci_contains_mean(self) -> None:
        r = _make_results(n_points=2, n_reps=20, with_noise=True)
        boot = r.bootstrap_ci("P_max", n_bootstrap=5000, seed=42)
        for row in boot.rows:
            assert row["P_max_boot_lo"] <= row["P_max"]
            assert row["P_max_boot_hi"] >= row["P_max"]


# ------------------------------------------------------------------
# Replicate convergence plot (smoke test)
# ------------------------------------------------------------------


class TestReplicateConvergencePlot:
    """Smoke test for render_replicate_convergence."""

    def test_render_runs_without_error(self) -> None:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt

        from tidal.cli._sweep_panels import render_replicate_convergence

        r = _make_results(n_points=2, n_reps=5, with_noise=True)
        fig, ax = plt.subplots()
        render_replicate_convergence(ax, r, "P_max")
        # Should not raise; check axis labels
        assert ax.get_xlabel() == "Number of replicates (k)"
        assert "P_max" in ax.get_ylabel()
        plt.close(fig)
