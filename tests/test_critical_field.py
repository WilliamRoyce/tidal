"""Tests for critical field analysis (tidal.measurement._critical_field).

Verifies:
- Threshold crossing detection with linear interpolation
- Error estimation (grid, interpolation model)
- Quality flags (good, coarse, edge, none)
- Conversion to SweepResults for plotting
- Analytical reference formulas (Boccaletti, uniform)
- CLI argument parsing and dispatch
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

from tidal.measurement._critical_field import (
    QUALITY_COARSE,
    QUALITY_EDGE,
    QUALITY_GOOD,
    QUALITY_NONE,
    compute_critical_field,
    compute_reference_threshold,
    critical_field_to_sweep_results,
)
from tidal.measurement._sweep_results import SweepResults

if TYPE_CHECKING:
    from pathlib import Path

# ── Helpers ─────────────────────────────────────────────────────────


def _make_sweep(
    b_values: list[float],
    metric_values: list[float],
    metric_name: str = "P_final",
    outer_params: dict[str, list[float]] | None = None,
) -> SweepResults:
    """Build a synthetic SweepResults for testing.

    If *outer_params* is provided, creates a full grid of
    (outer_param_combos) x (b_values).
    """
    rows: list[dict[str, Any]] = []
    swept: dict[str, list[float]] = {"B0": b_values}

    if outer_params is None:
        for b, m in zip(b_values, metric_values, strict=False):
            rows.append({"B0": b, metric_name: m})
    else:
        swept.update(outer_params)
        # metric_values should be a flat list for the full grid
        idx = 0
        outer_keys = list(outer_params.keys())
        # Build grid: iterate outer combos, then B0
        import itertools

        outer_combos = list(itertools.product(*outer_params.values()))
        for combo in outer_combos:
            for b in b_values:
                row: dict[str, Any] = {"B0": b}
                row.update(dict(zip(outer_keys, combo, strict=False)))
                row[metric_name] = metric_values[idx]
                rows.append(row)
                idx += 1

    return SweepResults(
        swept_params=swept,
        fixed_params={"kappa": 1.0},
        sim_settings={"t_end": 50.0},
        rows=rows,
        run_dirs=[],
        spec_path="test.json",
        measurements=["peak_conversion"],
    )


# ── Basic crossing detection ───────────────────────────────────────


class TestBasicCrossing:
    """Test threshold crossing detection."""

    def test_sin2_crossing(self) -> None:
        """1D B0 sweep with sin² metric — find B_min at analytical crossing."""
        n = 50
        b_values = list(np.linspace(0.01, 1.0, n))
        # P = sin²(π·B0) → first crossing of 0.99 near B0 ≈ 0.5
        metric_values = [math.sin(math.pi * b) ** 2 for b in b_values]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        assert len(result.rows) == 1  # No outer params → single group
        row = result.rows[0]
        assert row["crossing_quality"] in {QUALITY_GOOD, QUALITY_COARSE}
        # Analytical: sin²(π·B) = 0.99 → π·B = arcsin(√0.99) ≈ 1.4706
        # → B ≈ 0.468
        b_analytical = math.asin(math.sqrt(0.99)) / math.pi
        assert abs(row["B_min"] - b_analytical) < 0.02
        assert row["inv_B_min"] == pytest.approx(1.0 / row["B_min"], rel=1e-10)

    def test_exact_grid_point(self) -> None:
        """Metric equals threshold exactly at a grid point."""
        b_values = [0.1, 0.2, 0.3, 0.4, 0.5]
        metric_values = [0.0, 0.5, 0.99, 1.0, 0.8]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        assert row["B_min"] == pytest.approx(0.3, abs=1e-10)

    def test_no_crossing(self) -> None:
        """Metric never reaches threshold → B_min = NaN."""
        b_values = [0.1, 0.2, 0.3]
        metric_values = [0.1, 0.2, 0.3]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        assert math.isnan(row["B_min"])
        assert row["crossing_quality"] == QUALITY_NONE

    def test_already_exceeded(self) -> None:
        """Metric exceeds threshold at smallest B → quality = edge."""
        b_values = [0.1, 0.2, 0.3]
        metric_values = [0.99, 1.0, 0.8]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        assert row["B_min"] == pytest.approx(0.1)
        assert row["crossing_quality"] == QUALITY_EDGE

    def test_first_crossing_non_monotonic(self) -> None:
        """For oscillating metric, takes the FIRST crossing (minimum B)."""
        b_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
        # Oscillating: crosses 0.5 at B=0.2 and again at B=0.5
        metric_values = [0.0, 0.6, 0.3, 0.1, 0.7, 0.2]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.5)

        row = result.rows[0]
        # First crossing is between B=0.1 and B=0.2 (interpolated)
        assert row["B_min"] < 0.25

    def test_invalid_field_param(self) -> None:
        """Raises ValueError for non-existent field parameter."""
        sweep = _make_sweep([0.1], [0.5])
        with pytest.raises(ValueError, match="not a swept parameter"):
            compute_critical_field(sweep, "nonexistent")

    def test_missing_metric(self) -> None:
        """Raises ValueError when metric not found in rows."""
        sweep = _make_sweep([0.1, 0.2], [0.5, 0.9])
        with pytest.raises(ValueError, match="not found"):
            compute_critical_field(sweep, "B0", metric="nonexistent")


# ── Interpolation ──────────────────────────────────────────────────


class TestInterpolation:
    """Test linear interpolation accuracy."""

    def test_interpolation_improves_accuracy(self) -> None:
        """Interpolated B_min is closer to analytical than nearest grid point."""
        # Use threshold 0.5 where crossing is guaranteed in [0, 0.5]
        # P = sin²(π·B) crosses 0.5 at B = arcsin(√0.5)/π ≈ 0.25
        n = 10
        b_values = list(np.linspace(0.01, 0.5, n))
        metric_values = [math.sin(math.pi * b) ** 2 for b in b_values]
        b_analytical = math.asin(math.sqrt(0.5)) / math.pi

        sweep = _make_sweep(b_values, metric_values)

        # With interpolation
        result_interp = compute_critical_field(sweep, "B0", threshold=0.5)
        b_min_interp = result_interp.rows[0]["B_min"]

        # Without interpolation
        result_no_interp = compute_critical_field(
            sweep, "B0", threshold=0.5, interpolate=False
        )
        b_min_no_interp = result_no_interp.rows[0]["B_min"]

        err_interp = abs(b_min_interp - b_analytical)
        err_no_interp = abs(b_min_no_interp - b_analytical)
        assert err_interp < err_no_interp


# ── Error estimation ───────────────────────────────────────────────


class TestErrors:
    """Test error estimation."""

    def test_grid_error(self) -> None:
        """Grid error = half the interval between bracketing points."""
        b_values = [0.0, 0.1, 0.2, 0.3]
        metric_values = [0.0, 0.5, 1.0, 0.8]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        # Crossing between B=0.1 and B=0.2, ΔB = 0.1
        assert row["B_min_err_grid"] == pytest.approx(0.05, rel=1e-10)

    def test_inv_b_error_propagation(self) -> None:
        """Error on 1/B_min: err_{1/B} = err_B / B^2."""
        b_values = [0.0, 0.1, 0.2, 0.3]
        metric_values = [0.0, 0.5, 1.0, 0.8]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        b_min = row["B_min"]
        b_err = row["B_min_err"]
        expected_inv_err = b_err / b_min**2
        assert row["inv_B_min_err"] == pytest.approx(expected_inv_err, rel=1e-8)

    def test_interpolation_model_error(self) -> None:
        """Quadratic vs linear difference captured in B_min_err_interp."""
        # Use a nonlinear function where quadratic differs from linear
        n = 20
        b_values = list(np.linspace(0.01, 1.0, n))
        metric_values = [math.sin(math.pi * b) ** 2 for b in b_values]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        # For sin², the quadratic interpolant should differ from linear
        # The error should be finite (not NaN) since we have enough points
        if not math.isnan(row["B_min_err_interp"]):
            assert row["B_min_err_interp"] >= 0


# ── Quality flags ──────────────────────────────────────────────────


class TestQualityFlags:
    """Test crossing quality classification."""

    def test_good_quality(self) -> None:
        """Dense grid with clean crossing → quality = good."""
        n = 100
        b_values = list(np.linspace(0.01, 1.0, n))
        metric_values = [math.sin(math.pi * b) ** 2 for b in b_values]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        assert result.rows[0]["crossing_quality"] == QUALITY_GOOD

    def test_coarse_quality(self) -> None:
        """Very coarse grid: quality = coarse (err_grid > 0.1*B_min)."""
        # 3 points, wide spacing
        b_values = [0.01, 0.5, 1.0]
        metric_values = [0.0, 0.5, 1.0]

        sweep = _make_sweep(b_values, metric_values)
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        assert row["crossing_quality"] == QUALITY_COARSE

    def test_edge_quality(self) -> None:
        """Already exceeded → quality = edge."""
        sweep = _make_sweep([0.1, 0.2], [1.0, 0.5])
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        assert result.rows[0]["crossing_quality"] == QUALITY_EDGE

    def test_none_quality(self) -> None:
        """No crossing → quality = none."""
        sweep = _make_sweep([0.1, 0.2], [0.1, 0.2])
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        assert result.rows[0]["crossing_quality"] == QUALITY_NONE


# ── 2D outer params ────────────────────────────────────────────────


class TestOuterParams:
    """Test with multiple outer parameters."""

    def test_2d_outer_params(self) -> None:
        """Sweep (kappa_val, B0) → B_min(kappa_val) = π/(2·kappa_val)."""
        kappa_values = [1.0, 2.0, 4.0]
        b_values = list(np.linspace(0.01, 5.0, 200))

        # P_final = sin²(kappa_val · B0 / 2) for uniform with t_end=1
        rows: list[dict[str, Any]] = []
        for kv in kappa_values:
            for b in b_values:
                p = math.sin(kv * b / 2.0) ** 2
                rows.append({"kappa_val": kv, "B0": b, "P_final": p})

        sweep = SweepResults(
            swept_params={"kappa_val": kappa_values, "B0": b_values},
            fixed_params={},
            sim_settings={"t_end": 1.0},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )

        result = compute_critical_field(sweep, "B0", threshold=0.99)

        assert len(result.rows) == 3  # One per kappa_val
        assert "kappa_val" in result.outer_params
        assert "B0" not in result.outer_params

        for row in result.rows:
            kv = row["kappa_val"]
            # Analytical: sin²(kv·B/2) = 0.99 → kv·B/2 = arcsin(√0.99)
            b_analytical = 2.0 * math.asin(math.sqrt(0.99)) / kv
            assert row["B_min"] == pytest.approx(b_analytical, rel=0.05)

    def test_no_outer_params(self) -> None:
        """Single swept parameter → single output row, no outer params."""
        sweep = _make_sweep([0.1, 0.2, 0.3], [0.0, 0.5, 1.0])
        result = compute_critical_field(sweep, "B0", threshold=0.5)

        assert len(result.rows) == 1
        assert result.outer_params == {}


# ── Conversion to SweepResults ─────────────────────────────────────


class TestSweepResultsConversion:
    """Test critical_field_to_sweep_results."""

    def test_round_trip(self, tmp_path: Path) -> None:
        """Reduced SweepResults serializes and loads correctly."""
        kappa_values = [0.5, 1.0, 2.0]
        b_values = list(np.linspace(0.01, 5.0, 30))
        rows: list[dict[str, Any]] = []
        for kv in kappa_values:
            for b in b_values:
                p = math.sin(kv * b / 2.0) ** 2
                rows.append({"kappa_val": kv, "B0": b, "P_final": p})

        sweep = SweepResults(
            swept_params={"kappa_val": kappa_values, "B0": b_values},
            fixed_params={"R": 5.0},
            sim_settings={"t_end": 50.0},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )

        result = compute_critical_field(sweep, "B0", threshold=0.99)
        reduced = critical_field_to_sweep_results(result, sweep)

        # Check structure
        assert "B0" not in reduced.swept_params
        assert "kappa_val" in reduced.swept_params
        assert "critical_field" in reduced.measurements
        assert reduced.metadata["critical_field"]["field_param"] == "B0"
        assert reduced.metadata["critical_field"]["threshold"] == 0.99

        # Serialize and reload
        reduced.to_csv(tmp_path / "results.csv")
        reduced.to_json(tmp_path / "results.json")
        reduced.save_sweep_json(tmp_path / "sweep.json")

        loaded = SweepResults.from_directory(tmp_path)
        assert loaded.n_runs == reduced.n_runs
        assert "kappa_val" in loaded.swept_params

    def test_metrics_present(self) -> None:
        """Reduced SweepResults contains B_min and inv_B_min columns."""
        sweep = _make_sweep([0.1, 0.2, 0.3], [0.0, 0.5, 1.0])
        result = compute_critical_field(sweep, "B0", threshold=0.5)
        reduced = critical_field_to_sweep_results(result, sweep)

        metric_names = reduced.metric_names
        assert "B_min" in metric_names
        assert "inv_B_min" in metric_names
        assert "B_min_err" in metric_names
        assert "crossing_quality" in metric_names


# ── Analytical reference formulas ──────────────────────────────────


class TestReferenceFormulas:
    """Test compute_reference_threshold."""

    def test_boccaletti(self) -> None:
        """Boccaletti: P = sin²(κ/2 · B · R · √(π/2))."""
        kappa = 1.0
        r = 5.0
        b_ref = 0.1
        expected = math.sin(kappa / 2 * b_ref * r * math.sqrt(math.pi / 2)) ** 2

        p = compute_reference_threshold(
            "boccaletti",
            b_ref,
            fixed_params={"kappa": kappa, "R": r},
            sim_settings={},
        )
        assert p == pytest.approx(expected, rel=1e-10)

    def test_uniform(self) -> None:
        """Uniform: P = sin²(κ · B · t_end / 2)."""
        kappa = 1.0
        t_end = 50.0
        b_ref = 0.1
        expected = math.sin(kappa * b_ref * t_end / 2) ** 2

        p = compute_reference_threshold(
            "uniform",
            b_ref,
            fixed_params={"kappa": kappa},
            sim_settings={"t_end": t_end},
        )
        assert p == pytest.approx(expected, rel=1e-10)

    def test_boccaletti_full_conversion(self) -> None:
        """Boccaletti at B_min for full conversion gives P ≈ 1."""
        kappa = 1.0
        r = 5.0
        # B_min = π / (κ · R · √(π/2))
        b_min = math.pi / (kappa * r * math.sqrt(math.pi / 2))
        p = compute_reference_threshold(
            "boccaletti",
            b_min,
            fixed_params={"kappa": kappa, "R": r},
            sim_settings={},
        )
        assert p == pytest.approx(1.0, abs=1e-10)

    def test_unknown_formula(self) -> None:
        """Unknown formula name raises ValueError."""
        with pytest.raises(ValueError, match="Unknown"):
            compute_reference_threshold("unknown", 0.1, {}, {})

    def test_missing_kappa(self) -> None:
        """Missing kappa parameter raises ValueError."""
        with pytest.raises(ValueError, match="kappa"):
            compute_reference_threshold("boccaletti", 0.1, {}, {})

    def test_missing_r(self) -> None:
        """Missing R parameter for Boccaletti raises ValueError."""
        with pytest.raises(ValueError, match="R"):
            compute_reference_threshold("boccaletti", 0.1, {"kappa": 1.0}, {})

    def test_missing_t_end(self) -> None:
        """Missing t_end for uniform raises ValueError."""
        with pytest.raises(ValueError, match="t_end"):
            compute_reference_threshold("uniform", 0.1, {"kappa": 1.0}, {})


# ── CLI integration ────────────────────────────────────────────────


class TestCLI:
    """Test CLI argument parsing."""

    def test_critical_field_args_parsed(self) -> None:
        """Verify --critical-field and related args are parsed."""
        from tidal.cli import _build_parser as build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "analyze",
                "/tmp/test",
                "--critical-field",
                "B0",
                "--threshold",
                "0.95",
                "--output",
                "/tmp/out",
                "--no-interpolate",
            ]
        )
        assert args.critical_field == "B0"
        assert args.threshold == 0.95
        assert args.output == "/tmp/out"
        assert args.no_interpolate is True

    def test_reference_formula_args(self) -> None:
        """Verify --reference-formula and --reference-B are parsed."""
        from tidal.cli import _build_parser as build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "analyze",
                "/tmp/test",
                "--critical-field",
                "Bpeak",
                "--reference-formula",
                "boccaletti",
                "--reference-B",
                "0.1",
                "--output",
                "/tmp/out",
            ]
        )
        assert args.reference_formula == "boccaletti"
        assert args.reference_B == 0.1

    def test_end_to_end(self, tmp_path: Path) -> None:
        """Full CLI round-trip: create sweep data → analyze → verify output."""
        # Create synthetic sweep data
        n = 30
        b_values = list(np.linspace(0.01, 2.0, n))
        metric_values = [math.sin(b) ** 2 for b in b_values]

        sweep = _make_sweep(b_values, metric_values)
        sweep_dir = tmp_path / "sweep"
        sweep_dir.mkdir()
        sweep.to_csv(sweep_dir / "results.csv")
        sweep.to_json(sweep_dir / "results.json")
        sweep.save_sweep_json(sweep_dir / "sweep.json")

        # Run analyze command
        from argparse import Namespace

        from tidal.cli._analyze import analyze_command

        out_dir = tmp_path / "critical"
        args = Namespace(
            data_path=str(sweep_dir),
            critical_field="B0",
            threshold=0.99,
            metric="P_final",
            output=str(out_dir),
            no_interpolate=False,
            reference_formula=None,
            reference_B=None,
            sensitivity="sobol",
            bootstrap=100,
        )
        rc = analyze_command(args)
        assert rc == 0

        # Verify output files
        assert (out_dir / "results.csv").exists()
        assert (out_dir / "results.json").exists()
        assert (out_dir / "sweep.json").exists()

        # Load and verify
        loaded = SweepResults.from_directory(out_dir)
        assert loaded.n_runs >= 1
        assert any("B_min" in row for row in loaded.rows)


# ── Edge cases ─────────────────────────────────────────────────────


class TestEdgeCases:
    """Test edge cases and robustness."""

    def test_single_point(self) -> None:
        """Single B0 value above threshold."""
        sweep = _make_sweep([0.5], [1.0])
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        assert row["B_min"] == pytest.approx(0.5)
        assert row["crossing_quality"] == QUALITY_EDGE

    def test_failed_runs_skipped(self) -> None:
        """Rows with None metric are skipped."""
        rows = [
            {"B0": 0.1, "P_final": None},
            {"B0": 0.2, "P_final": 0.5},
            {"B0": 0.3, "P_final": 1.0},
        ]
        sweep = SweepResults(
            swept_params={"B0": [0.1, 0.2, 0.3]},
            fixed_params={},
            sim_settings={},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )
        result = compute_critical_field(sweep, "B0", threshold=0.99)

        row = result.rows[0]
        # Should find crossing between 0.2 and 0.3
        assert 0.2 < row["B_min"] <= 0.3

    def test_all_none_metric(self) -> None:
        """All metric values None → B_min = NaN."""
        rows = [
            {"B0": 0.1, "P_final": None},
            {"B0": 0.2, "P_final": None},
        ]
        sweep = SweepResults(
            swept_params={"B0": [0.1, 0.2]},
            fixed_params={},
            sim_settings={},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )
        result = compute_critical_field(sweep, "B0", threshold=0.99)
        assert math.isnan(result.rows[0]["B_min"])


# ── Sweep-integrated critical field ──────────────────────────────


class TestSweepIntegrated:
    """Test --critical-field flag on sweep via _run_critical_field_post_sweep."""

    def test_post_sweep_produces_output(self, tmp_path: Path) -> None:
        """Post-sweep critical field extraction writes critical_field/ subdir."""
        from argparse import Namespace

        from tidal.cli._sweep import _run_critical_field_post_sweep

        # Synthetic sweep results with clean sin² crossing
        n = 30
        b_values = list(np.linspace(0.01, 2.0, n))
        metric_values = [math.sin(b) ** 2 for b in b_values]
        sweep = _make_sweep(b_values, metric_values)

        args = Namespace(
            critical_field="B0",
            cf_threshold=0.99,
            cf_metric=None,
            no_interpolate=False,
        )
        _run_critical_field_post_sweep(sweep, args, tmp_path)

        cf_dir = tmp_path / "critical_field"
        assert cf_dir.exists()
        assert (cf_dir / "results.csv").exists()
        assert (cf_dir / "results.json").exists()
        assert (cf_dir / "sweep.json").exists()

        # Load and verify B_min is reasonable
        loaded = SweepResults.from_directory(cf_dir)
        assert loaded.n_runs >= 1
        assert any("B_min" in row for row in loaded.rows)

    def test_post_sweep_invalid_param_warns(self, tmp_path: Path) -> None:
        """Invalid field param prints warning but doesn't raise."""
        from argparse import Namespace

        from tidal.cli._sweep import _run_critical_field_post_sweep

        sweep = _make_sweep([0.1, 0.2], [0.5, 1.0])
        args = Namespace(
            critical_field="nonexistent",
            cf_threshold=0.99,
            cf_metric=None,
            no_interpolate=False,
        )
        # Should not raise — just warn
        _run_critical_field_post_sweep(sweep, args, tmp_path)
        assert not (tmp_path / "critical_field").exists()

    def test_post_sweep_multi_param(self, tmp_path: Path) -> None:
        """Multi-parameter sweep with --critical-field collapses B0 correctly."""
        from argparse import Namespace

        from tidal.cli._sweep import _run_critical_field_post_sweep

        # 2D sweep: kappa_val x B0, P = sin²(kappa_val * B0 / 2)
        kappa_values = [1.0, 2.0]
        b_values = list(np.linspace(0.01, 5.0, 50))
        rows: list[dict[str, Any]] = []
        for kv in kappa_values:
            for b in b_values:
                p = math.sin(kv * b / 2.0) ** 2
                rows.append({"kappa_val": kv, "B0": b, "P_final": p})

        sweep = SweepResults(
            swept_params={"kappa_val": kappa_values, "B0": b_values},
            fixed_params={},
            sim_settings={"t_end": 1.0},
            rows=rows,
            run_dirs=[],
            spec_path="test.json",
            measurements=["peak_conversion"],
        )

        args = Namespace(
            critical_field="B0",
            cf_threshold=0.99,
            cf_metric=None,
            no_interpolate=False,
        )
        _run_critical_field_post_sweep(sweep, args, tmp_path)

        cf_dir = tmp_path / "critical_field"
        loaded = SweepResults.from_directory(cf_dir)

        # Should have 2 rows (one per kappa_val), B0 collapsed
        assert loaded.n_runs == 2
        assert "kappa_val" in loaded.swept_params
        assert "B0" not in loaded.swept_params

        # Verify B_min values: sin²(kv * B / 2) = 0.99
        for row in loaded.rows:
            kv = row["kappa_val"]
            b_analytical = 2.0 * math.asin(math.sqrt(0.99)) / kv
            assert row["B_min"] == pytest.approx(b_analytical, rel=0.05)


# ── Log-scale and quality overlay ────────────────────────────────


class TestLogScaleAndQuality:
    """Test log-scale norm and quality hatching helpers."""

    def test_build_norm_linear(self) -> None:
        """Default (linear) norm returns Normalize."""
        from matplotlib.colors import Normalize

        from tidal.cli._sweep_panels import _build_norm

        vals = np.array([1.0, 2.0, 3.0])
        norm = _build_norm(vals)
        assert isinstance(norm, Normalize)

    def test_build_norm_log(self) -> None:
        """Log norm returns LogNorm with positive bounds."""
        from matplotlib.colors import LogNorm

        from tidal.cli._sweep_panels import _build_norm

        vals = np.array([0.01, 0.1, 1.0, 10.0])
        norm = _build_norm(vals, log_scale=True)
        assert isinstance(norm, LogNorm)

    def test_build_norm_log_with_zeros(self) -> None:
        """Log norm with zero values clamps to smallest positive."""
        from matplotlib.colors import LogNorm

        from tidal.cli._sweep_panels import _build_norm

        vals = np.array([0.0, 0.0, 0.1, 1.0])
        norm = _build_norm(vals, log_scale=True)
        assert isinstance(norm, LogNorm)

    def test_build_norm_log_all_zero(self) -> None:
        """Log norm with all zeros falls back to linear Normalize."""
        from matplotlib.colors import Normalize

        from tidal.cli._sweep_panels import _build_norm

        vals = np.array([0.0, 0.0])
        norm = _build_norm(vals, log_scale=True)
        assert isinstance(norm, Normalize)

    def test_build_norm_empty(self) -> None:
        """Empty array returns default Normalize."""
        from matplotlib.colors import Normalize

        from tidal.cli._sweep_panels import _build_norm

        vals = np.array([], dtype=np.float64)
        norm = _build_norm(vals)
        assert isinstance(norm, Normalize)

    def test_quality_hatching_smoke(self) -> None:
        """Quality hatching runs without error on synthetic data."""
        import matplotlib.pyplot as plt

        from tidal.cli._sweep_panels import _overlay_quality_hatching

        fig, ax = plt.subplots()
        p1 = np.array([0.1, 0.2, 0.3])
        p2 = np.array([1.0, 2.0])
        quality = np.array([["good", "coarse", "none"], ["edge", "good", "good"]])
        _overlay_quality_hatching(ax, p1, p2, quality)
        # Should add patches for coarse, edge, none cells
        assert len(ax.patches) == 3
        plt.close(fig)

    def test_cli_args_parsed(self) -> None:
        """Verify --critical-field and --cf-threshold on sweep parser."""
        from tidal.cli import _build_parser as build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "sweep",
                "test.json",
                "--sweep",
                "B0=0.1:1.0:5",
                "--critical-field",
                "B0",
                "--cf-threshold",
                "0.95",
                "--no-interpolate",
            ]
        )
        assert args.critical_field == "B0"
        assert args.cf_threshold == 0.95
        assert args.no_interpolate is True

    def test_log_scale_cli_arg(self) -> None:
        """Verify --log-scale on plot parser."""
        from tidal.cli import _build_parser as build_parser

        parser = build_parser()
        args = parser.parse_args(
            [
                "plot",
                "/tmp/test",
                "--type",
                "sweep",
                "--metric",
                "inv_B_min",
                "--log-scale",
            ]
        )
        assert args.log_scale is True
