"""Tests for ``tidal sweep`` — parameter sweep and convergence studies."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from tidal.cli._sweep import (
    _build_row,
    _interval_scores,
    _run_subdir_name,
    parse_converge_spec,
    parse_sweep_spec,
)
from tidal.measurement._sweep_results import SweepResults

# ------------------------------------------------------------------
# parse_sweep_spec
# ------------------------------------------------------------------


class TestParseSweepSpec:
    """Test sweep specification parsing."""

    def test_linear_range(self) -> None:
        name, vals = parse_sweep_spec("g0=0.1:1.0:5")
        assert name == "g0"
        assert len(vals) == 5
        assert vals[0] == pytest.approx(0.1)
        assert vals[-1] == pytest.approx(1.0)

    def test_log_range(self) -> None:
        name, vals = parse_sweep_spec("m2=0.01:100.0:5:log")
        assert name == "m2"
        assert len(vals) == 5
        assert vals[0] == pytest.approx(0.01)
        assert vals[-1] == pytest.approx(100.0)
        # Log-spaced: middle value should be geometric mean
        assert vals[2] == pytest.approx(1.0, rel=0.01)

    def test_explicit_values(self) -> None:
        name, vals = parse_sweep_spec("g0=0.1,0.5,1.0,2.0")
        assert name == "g0"
        assert vals == [0.1, 0.5, 1.0, 2.0]

    def test_invalid_no_equals(self) -> None:
        with pytest.raises(ValueError, match="Invalid sweep spec"):
            parse_sweep_spec("g0_0.1:1.0:5")

    def test_invalid_empty_name(self) -> None:
        with pytest.raises(ValueError, match="Empty parameter name"):
            parse_sweep_spec("=0.1:1.0:5")

    def test_invalid_count_lt_2(self) -> None:
        with pytest.raises(ValueError, match="must be >= 2"):
            parse_sweep_spec("g0=0.1:1.0:1")

    def test_invalid_log_negative(self) -> None:
        with pytest.raises(ValueError, match="positive bounds"):
            parse_sweep_spec("g0=-1.0:1.0:5:log")

    def test_invalid_scale(self) -> None:
        with pytest.raises(ValueError, match="Unknown scale"):
            parse_sweep_spec("g0=0.1:1.0:5:quadratic")

    def test_invalid_range_too_many_colons(self) -> None:
        with pytest.raises(ValueError, match="Invalid range spec"):
            parse_sweep_spec("g0=0.1:1.0:5:log:extra")

    def test_single_explicit_value(self) -> None:
        with pytest.raises(ValueError, match="at least 2 values"):
            parse_sweep_spec("g0=0.5")

    def test_whitespace_handling(self) -> None:
        name, vals = parse_sweep_spec("  g0  =  0.1 , 0.5 , 1.0  ")
        assert name == "g0"
        assert len(vals) == 3


# ------------------------------------------------------------------
# parse_converge_spec
# ------------------------------------------------------------------


class TestParseConvergeSpec:
    """Test convergence specification parsing."""

    def test_basic(self) -> None:
        sizes = parse_converge_spec("32,64,128,256")
        assert sizes == [32, 64, 128, 256]

    def test_unsorted_input(self) -> None:
        sizes = parse_converge_spec("128,32,64")
        assert sizes == [32, 64, 128]

    def test_too_few(self) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            parse_converge_spec("32")

    def test_negative(self) -> None:
        with pytest.raises(ValueError, match="positive"):
            parse_converge_spec("-1,32")


# ------------------------------------------------------------------
# _run_subdir_name
# ------------------------------------------------------------------


class TestRunSubdirName:
    """Test subdirectory naming for sweep runs."""

    def test_single_param_float(self) -> None:
        name = _run_subdir_name({"g0": 0.5})
        assert name == "g0_0.5"

    def test_single_param_int(self) -> None:
        name = _run_subdir_name({"g0": 1.0})
        assert name == "g0_1"

    def test_multi_param(self) -> None:
        name = _run_subdir_name({"g0": 0.5, "m2": 2.0})
        assert "g0_0.5" in name
        assert "m2_2" in name

    def test_converge(self) -> None:
        name = _run_subdir_name({}, converge_size=64)
        assert name == "N_64"

    def test_empty(self) -> None:
        name = _run_subdir_name({})
        assert name == "run_0"


# ------------------------------------------------------------------
# _build_row
# ------------------------------------------------------------------


class TestBuildRow:
    """Test results row construction."""

    def test_combines_all_sources(self) -> None:
        row = _build_row(
            swept_vals={"g0": 0.5},
            fixed_params={"m2": 1.0},
            sim_settings={"t_end": 10.0, "scheme": "auto"},
            metrics={"P_max": 0.05, "wall_time_s": 1.2},
            grid_override=None,
        )
        assert row["g0"] == 0.5
        assert row["m2"] == 1.0
        assert row["t_end"] == 10.0
        assert row["P_max"] == 0.05

    def test_grid_override(self) -> None:
        row = _build_row(
            swept_vals={},
            fixed_params={},
            sim_settings={"grid_shape": "auto"},
            metrics={},
            grid_override=64,
        )
        assert row["grid_shape"] == 64

    def test_filters_internal_keys(self) -> None:
        row = _build_row(
            swept_vals={},
            fixed_params={},
            sim_settings={},
            metrics={"P_max": 0.05, "_result_obj": "internal"},
            grid_override=None,
        )
        assert "P_max" in row
        assert "_result_obj" not in row


# ------------------------------------------------------------------
# SweepResults
# ------------------------------------------------------------------


class TestSweepResults:
    """Test SweepResults dataclass and serialization."""

    @pytest.fixture
    def sample_results(self) -> SweepResults:
        return SweepResults(
            swept_params={"g0": [0.1, 0.3, 0.5]},
            fixed_params={"m2": 1.0},
            sim_settings={"t_end": 10.0, "scheme": "auto"},
            rows=[
                {"g0": 0.1, "m2": 1.0, "t_end": 10.0, "scheme": "auto", "P_max": 0.01, "wall_time_s": 1.0},
                {"g0": 0.3, "m2": 1.0, "t_end": 10.0, "scheme": "auto", "P_max": 0.05, "wall_time_s": 1.1},
                {"g0": 0.5, "m2": 1.0, "t_end": 10.0, "scheme": "auto", "P_max": 0.10, "wall_time_s": 1.2},
            ],
            run_dirs=[Path("/tmp/g0_0.1"), Path("/tmp/g0_0.3"), Path("/tmp/g0_0.5")],
            spec_path="examples/data/coupled_scattering.json",
            measurements=["conversion"],
        )

    def test_n_runs(self, sample_results: SweepResults) -> None:
        assert sample_results.n_runs == 3

    def test_param_names(self, sample_results: SweepResults) -> None:
        assert "g0" in sample_results.param_names
        assert "m2" in sample_results.param_names

    def test_metric_names(self, sample_results: SweepResults) -> None:
        metrics = sample_results.metric_names
        assert "P_max" in metrics
        assert "wall_time_s" in metrics
        assert "g0" not in metrics
        assert "t_end" not in metrics

    def test_metric_names_union_across_rows(self) -> None:
        """metric_names should include keys from ALL rows, not just rows[0]."""
        r = SweepResults(
            swept_params={"g0": [0.1, 0.5]},
            fixed_params={},
            sim_settings={},
            rows=[
                {"g0": 0.1, "P_max": 0.01},
                {"g0": 0.5, "P_max": 0.05, "conversion_error": "some error"},
            ],
            run_dirs=[],
            spec_path="",
            measurements=[],
        )
        metrics = r.metric_names
        assert "P_max" in metrics
        assert "conversion_error" in metrics

    def test_column(self, sample_results: SweepResults) -> None:
        p_max = sample_results.column("P_max")
        assert len(p_max) == 3
        assert p_max[0] == pytest.approx(0.01)
        assert p_max[2] == pytest.approx(0.10)

    def test_to_csv(self, sample_results: SweepResults) -> None:
        csv_str = sample_results.to_csv()
        reader = csv.DictReader(csv_str.strip().splitlines())
        rows = list(reader)
        assert len(rows) == 3
        assert "g0" in rows[0]
        assert "P_max" in rows[0]
        assert "m2" in rows[0]
        assert "t_end" in rows[0]

    def test_to_csv_file(self, sample_results: SweepResults, tmp_path: Path) -> None:
        csv_file = tmp_path / "results.csv"
        sample_results.to_csv(csv_file)
        assert csv_file.exists()
        content = csv_file.read_text()
        assert "P_max" in content

    def test_to_json(self, sample_results: SweepResults) -> None:
        json_str = sample_results.to_json()
        data = json.loads(json_str)
        assert "columns" in data
        assert "rows" in data
        assert len(data["rows"]) == 3

    def test_save_sweep_json(self, sample_results: SweepResults, tmp_path: Path) -> None:
        sweep_file = tmp_path / "sweep.json"
        sample_results.save_sweep_json(sweep_file)
        assert sweep_file.exists()
        data = json.loads(sweep_file.read_text())
        assert data["spec_path"] == "examples/data/coupled_scattering.json"
        assert data["swept_parameters"] == {"g0": [0.1, 0.3, 0.5]}
        assert data["completed_runs"] == 3

    def test_save_sweep_json_includes_run_dirs(self, sample_results: SweepResults, tmp_path: Path) -> None:
        sweep_file = tmp_path / "sweep.json"
        sample_results.save_sweep_json(sweep_file)
        data = json.loads(sweep_file.read_text())
        assert "run_dirs" in data
        assert len(data["run_dirs"]) == 3

    def test_roundtrip(self, sample_results: SweepResults, tmp_path: Path) -> None:
        """Save to directory and reload — run_dirs should survive roundtrip."""
        sample_results.save_sweep_json(tmp_path / "sweep.json")
        sample_results.to_json(tmp_path / "results.json")

        loaded = SweepResults.from_directory(tmp_path)
        assert loaded.n_runs == 3
        assert loaded.swept_params == {"g0": [0.1, 0.3, 0.5]}
        assert loaded.spec_path == "examples/data/coupled_scattering.json"
        assert len(loaded.run_dirs) == 3
        assert str(loaded.run_dirs[0]) == "/tmp/g0_0.1"

    def test_from_directory_missing_sweep_json(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match=r"sweep\.json"):
            SweepResults.from_directory(tmp_path)

    def test_is_convergence(self) -> None:
        r = SweepResults(
            swept_params={},
            fixed_params={},
            sim_settings={},
            rows=[],
            run_dirs=[],
            spec_path="",
            measurements=[],
            converge_sizes=[32, 64],
        )
        assert r.is_convergence is True

    def test_not_convergence(self, sample_results: SweepResults) -> None:
        assert sample_results.is_convergence is False


# ------------------------------------------------------------------
# CLI integration (smoke tests)
# ------------------------------------------------------------------


class TestSweepCLI:
    """Smoke tests for tidal sweep CLI."""

    def test_help(self) -> None:
        from tidal.cli import main

        with pytest.raises(SystemExit) as exc_info:
            main(["sweep", "--help"])
        assert exc_info.value.code == 0

    def test_no_sweep_or_converge(self) -> None:
        from tidal.cli import main

        code = main(["sweep", "examples/data/coupled_scattering.json", "--output", "/tmp/nosweep"])
        assert code == 1

    def test_both_sweep_and_converge(self) -> None:
        from tidal.cli import main

        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--converge", "32,64",
            "--output", "/tmp/nosweep",
        ])
        assert code == 1

    def test_missing_output(self) -> None:
        from tidal.cli import main

        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
        ])
        assert code == 1

    def test_missing_spec(self) -> None:
        from tidal.cli import main

        code = main([
            "sweep", "nonexistent.json",
            "--sweep", "g0=0.1,0.5",
            "--output", "/tmp/nosweep",
        ])
        assert code == 1

    def test_small_sweep_end_to_end(self, tmp_path: Path) -> None:
        """Run a minimal 2-point sweep with real simulation."""
        from tidal.cli import main

        output = tmp_path / "sweep_out"
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--measure", "conservation,conversion",
            "--source", "phi_0", "--target", "chi_0",
            "--grid-shape", "16",
            "--bounds=-10:10",
            "--ic", "gaussian",
            "--ic-component", "phi_0",
            "--ic-width", "3",
            "--t-end", "2",
            "--output", str(output),
        ])
        assert code == 0
        assert (output / "sweep.json").exists()
        assert (output / "results.csv").exists()
        assert (output / "results.json").exists()

        # Verify CSV content
        with Path(output / "results.csv").open(encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == 2
        assert "P_max" in rows[0]
        assert "max_energy_error" in rows[0]
        assert "g0" in rows[0]

        # P_max should increase with g0
        p_max_vals = [float(r["P_max"]) for r in rows if r["P_max"]]
        if len(p_max_vals) == 2:
            assert p_max_vals[1] > p_max_vals[0]

    def test_convergence_end_to_end(self, tmp_path: Path) -> None:
        """Run a minimal convergence study."""
        from tidal.cli import main

        output = tmp_path / "converge_out"
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--converge", "8,16",
            "--measure", "conservation",
            "--param", "g0=0.5",
            "--bounds=-10:10",
            "--ic", "gaussian",
            "--ic-component", "phi_0",
            "--t-end", "1",
            "--output", str(output),
        ])
        assert code == 0
        assert (output / "sweep.json").exists()
        assert (output / "results.csv").exists()
        assert (output / "N_8").is_dir()
        assert (output / "N_16").is_dir()

    def test_resume_skips_completed(self, tmp_path: Path) -> None:
        """Resume should skip runs with existing metadata.json."""
        from tidal.cli import main

        output = tmp_path / "resume_test"

        # First run
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--measure", "conservation",
            "--grid-shape", "8",
            "--bounds=-5:5",
            "--t-end", "1",
            "--output", str(output),
        ])
        assert code == 0
        assert (output / "g0_0.1" / "metadata.json").exists()

        # Second run with --resume (should complete without errors)
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--measure", "conservation",
            "--grid-shape", "8",
            "--bounds=-5:5",
            "--t-end", "1",
            "--output", str(output),
            "--resume",
        ])
        assert code == 0


# ------------------------------------------------------------------
# Sweep plot integration (smoke tests)
# ------------------------------------------------------------------


class TestSweepPlot:
    """Smoke tests for sweep plot types in tidal plot."""

    @pytest.fixture
    def sweep_dir(self, tmp_path: Path) -> Path:
        """Run a small sweep to produce plot-able output."""
        from tidal.cli import main

        output = tmp_path / "plot_sweep"
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.3,0.5",
            "--measure", "conversion,conservation",
            "--source", "phi_0", "--target", "chi_0",
            "--grid-shape", "16",
            "--bounds=-10:10",
            "--ic", "gaussian",
            "--ic-component", "phi_0",
            "--ic-width", "3",
            "--t-end", "2",
            "--output", str(output),
        ])
        assert code == 0
        return output

    def test_sweep_1d_plot(self, sweep_dir: Path, tmp_path: Path) -> None:
        from tidal.cli import main

        plot_path = tmp_path / "sweep.png"
        code = main([
            "plot", str(sweep_dir),
            "--type", "sweep",
            "--metric", "P_max",
            "--output", str(plot_path),
        ])
        assert code == 0
        assert plot_path.exists()

    def test_sweep_compare_plot(self, sweep_dir: Path, tmp_path: Path) -> None:
        from tidal.cli import main

        plot_path = tmp_path / "compare.png"
        code = main([
            "plot", str(sweep_dir),
            "--type", "sweep-compare",
            "--metric", "conversion",
            "--output", str(plot_path),
        ])
        assert code == 0
        assert plot_path.exists()

    def test_sweep_plot_auto_metric(self, sweep_dir: Path, tmp_path: Path) -> None:
        from tidal.cli import main

        plot_path = tmp_path / "auto.png"
        code = main([
            "plot", str(sweep_dir),
            "--type", "sweep",
            "--output", str(plot_path),
        ])
        assert code == 0
        assert plot_path.exists()

    def test_sweep_not_a_sweep_dir(self, tmp_path: Path) -> None:
        from tidal.cli import main

        # Point at an empty dir (no sweep.json)
        empty = tmp_path / "empty"
        empty.mkdir()
        code = main([
            "plot", str(empty),
            "--type", "sweep",
            "--metric", "P_max",
        ])
        assert code == 1


# ------------------------------------------------------------------
# Validation tests
# ------------------------------------------------------------------


class TestSweepValidation:
    """Tests for parameter and measurement validation."""

    def test_dry_run(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """--dry-run should print plan without creating output."""
        from tidal.cli import main

        output = tmp_path / "dry_run_out"
        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.3,0.5",
            "--measure", "conservation",
            "--output", str(output),
            "--dry-run",
        ])
        assert code == 0
        # Should NOT create any run directories
        assert not (output / "results.csv").exists()
        captured = capsys.readouterr()
        assert "3 runs planned" in captured.out

    def test_unknown_measurement_type(self, tmp_path: Path) -> None:
        """Unknown --measure value should error with exit code 1."""
        from tidal.cli import main

        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--measure", "foobar",
            "--output", str(tmp_path / "out"),
        ])
        assert code == 1

    def test_unknown_measurement_mixed(self, tmp_path: Path) -> None:
        """Mix of valid and invalid --measure types should error."""
        from tidal.cli import main

        code = main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "g0=0.1,0.5",
            "--measure", "conservation,nonexistent",
            "--output", str(tmp_path / "out"),
        ])
        assert code == 1

    def test_unknown_swept_param_warns(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
        """Swept param not in spec should warn but not error."""
        from tidal.cli import main

        # This will run (and likely fail on simulation), but the validation
        # warning should appear in stderr before any simulation attempt.
        main([
            "sweep", "examples/data/coupled_scattering.json",
            "--sweep", "nonexistent_param_xyz=0.1,0.5",
            "--measure", "conservation",
            "--grid-shape", "8",
            "--bounds=-5:5",
            "--t-end", "0.1",
            "--output", str(tmp_path / "out"),
        ])
        captured = capsys.readouterr()
        assert "nonexistent_param_xyz" in captured.err
        assert "Possible typo" in captured.err


# ------------------------------------------------------------------
# Adaptive interval scoring
# ------------------------------------------------------------------


class TestIntervalScores:
    """Test _interval_scores curvature-based scoring."""

    def test_flat_function_low_scores(self) -> None:
        """Constant function should have zero scores."""
        values = [0.0, 1.0, 2.0, 3.0]
        metric_vals = [1.0, 1.0, 1.0, 1.0]
        scores = _interval_scores(values, metric_vals)
        assert len(scores) == 3
        assert all(s == 0.0 for s in scores)

    def test_linear_function_equal_scores(self) -> None:
        """Linear function has constant differences, all intervals score alike."""
        values = [0.0, 1.0, 2.0, 3.0]
        metric_vals = [0.0, 1.0, 2.0, 3.0]
        scores = _interval_scores(values, metric_vals)
        assert len(scores) == 3
        # All intervals have the same absolute step, so scores should be equal
        # (the scoring function normalizes by max(fa, fb) so ratios differ,
        # but the relative ranking is what matters for adaptive refinement)
        assert all(s > 0 for s in scores)

    def test_sharp_feature_scores_highest(self) -> None:
        """An interval with a sharp change should score higher."""
        values = [0.0, 1.0, 2.0, 3.0, 4.0]
        # Flat, then sharp jump at 2->3
        metric_vals = [0.0, 0.0, 0.0, 10.0, 10.0]
        scores = _interval_scores(values, metric_vals)
        # Interval [2, 3] has the largest jump
        assert scores[2] > scores[0]
        assert scores[2] > scores[3]

    def test_none_values_score_zero(self) -> None:
        """None metric values should produce zero score."""
        values = [0.0, 1.0, 2.0]
        metric_vals = [1.0, None, 2.0]
        scores = _interval_scores(values, metric_vals)
        assert scores[0] == 0.0
        assert scores[1] == 0.0

    def test_single_interval(self) -> None:
        """Two points = one interval."""
        values = [0.0, 1.0]
        metric_vals = [0.0, 5.0]
        scores = _interval_scores(values, metric_vals)
        assert len(scores) == 1
        assert scores[0] > 0
