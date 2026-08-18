"""Tests for ``tidal sweep --config`` TOML configuration."""

from __future__ import annotations

import argparse
from typing import TYPE_CHECKING

import pytest

from tidal.cli._sweep_config import (
    SweepConfig,
    _parse_sweep_section,
    apply_config_to_args,
    load_sweep_config,
)

if TYPE_CHECKING:
    from pathlib import Path

# ------------------------------------------------------------------
# _parse_sweep_section
# ------------------------------------------------------------------


class TestParseSweepSection:
    """Test individual [sweep.PARAM] section parsing."""

    def test_linear_range(self) -> None:
        values, adaptive = _parse_sweep_section(
            "g0",
            {"start": 0.1, "stop": 1.0, "count": 5},
        )
        assert len(values) == 5
        assert values[0] == pytest.approx(0.1)
        assert values[-1] == pytest.approx(1.0)
        assert adaptive == {}

    def test_log_range(self) -> None:
        values, adaptive = _parse_sweep_section(
            "m2",
            {"start": 0.01, "stop": 100.0, "count": 5, "scale": "log"},
        )
        assert len(values) == 5
        assert values[0] == pytest.approx(0.01)
        assert values[-1] == pytest.approx(100.0)
        assert adaptive == {}

    def test_explicit_values(self) -> None:
        values, adaptive = _parse_sweep_section("g0", {"values": [0.1, 0.5, 1.0]})
        assert values == [0.1, 0.5, 1.0]
        assert adaptive == {}

    def test_adaptive_mode(self) -> None:
        values, adaptive = _parse_sweep_section(
            "g0",
            {
                "start": 0.0,
                "stop": 1.0,
                "scale": "adaptive",
                "initial_count": 3,
                "max_count": 20,
                "metric": "P_max",
                "threshold": 0.05,
            },
        )
        assert len(values) == 3
        assert values[0] == pytest.approx(0.0)
        assert values[-1] == pytest.approx(1.0)
        assert adaptive["bounds"] == (0.0, 1.0)
        assert adaptive["max_count"] == 20
        assert adaptive["metric"] == "P_max"

    def test_missing_count_errors(self) -> None:
        with pytest.raises(ValueError, match="missing required key 'count'"):
            _parse_sweep_section("g0", {"start": 0.0, "stop": 1.0})

    def test_count_too_small(self) -> None:
        with pytest.raises(ValueError, match="must be >= 2"):
            _parse_sweep_section("g0", {"start": 0.0, "stop": 1.0, "count": 1})

    def test_log_negative_bounds(self) -> None:
        with pytest.raises(ValueError, match="positive bounds"):
            _parse_sweep_section(
                "g0",
                {"start": -1.0, "stop": 1.0, "count": 5, "scale": "log"},
            )

    def test_unknown_scale(self) -> None:
        with pytest.raises(ValueError, match="unknown scale"):
            _parse_sweep_section(
                "g0",
                {"start": 0.0, "stop": 1.0, "count": 5, "scale": "cubic"},
            )

    def test_explicit_values_too_few(self) -> None:
        with pytest.raises(ValueError, match="must be a list with >= 2"):
            _parse_sweep_section("g0", {"values": [0.5]})


# ------------------------------------------------------------------
# load_sweep_config
# ------------------------------------------------------------------


class TestLoadSweepConfig:
    """Test full TOML file loading."""

    def test_minimal_config(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        # Create a dummy spec file
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n[sweep.g0]\nstart = 0.1\nstop = 1.0\ncount = 5\n',
        )
        config = load_sweep_config(toml)
        assert config.spec_path == spec
        assert "g0" in config.swept_params
        assert len(config.swept_params["g0"]) == 5

    def test_full_config(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n'
            "\n"
            "[parameters]\n"
            "mPhi2 = 1.0\n"
            "mChi2 = 4.0\n"
            "\n"
            "[sweep.g0]\n"
            "start = 0.01\n"
            "stop = 1.0\n"
            "count = 8\n"
            'scale = "log"\n'
            "\n"
            "[simulation]\n"
            "grid_shape = 128\n"
            "bounds = [0, 100]\n"
            "periodic = true\n"
            "t_end = 20.0\n"
            "\n"
            "[measurement]\n"
            'types = ["conversion", "mixing"]\n'
            'source = "phi_0"\n'
            'target = "chi_0"\n'
            "\n"
            "[output]\n"
            'path = "sweep_out"\n'
            "\n"
            "[execution]\n"
            "parallel = 4\n"
            "resume = true\n",
        )
        config = load_sweep_config(toml)
        assert config.fixed_params == {"mPhi2": 1.0, "mChi2": 4.0}
        assert len(config.swept_params["g0"]) == 8
        assert config.sim_settings["grid_shape"] == 128
        assert config.sim_settings["periodic"] is True
        assert config.measurements == ["conversion", "mixing"]
        assert config.source == "phi_0"
        assert config.target == "chi_0"
        assert config.output == tmp_path / "sweep_out"
        assert config.parallel == 4
        assert config.resume is True

    def test_convergence_mode(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n[convergence]\ngrid_sizes = [32, 64, 128, 256]\n',
        )
        config = load_sweep_config(toml)
        assert config.converge_sizes == [32, 64, 128, 256]
        assert config.swept_params == {}

    def test_missing_spec_errors(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        toml.write_text("[sweep.g0]\nstart = 0.1\nstop = 1.0\ncount = 5\n")
        with pytest.raises(ValueError, match="Missing required key 'spec'"):
            load_sweep_config(toml)

    def test_no_sweep_or_converge_errors(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text('spec = "spec.json"\n')
        with pytest.raises(ValueError, match=r"No .sweep"):
            load_sweep_config(toml)

    def test_sweep_and_converge_mutually_exclusive(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n'
            "[sweep.g0]\n"
            "start = 0.1\n"
            "stop = 1.0\n"
            "count = 5\n"
            "[convergence]\n"
            "grid_sizes = [32, 64]\n",
        )
        with pytest.raises(ValueError, match="mutually exclusive"):
            load_sweep_config(toml)

    def test_path_resolution_relative_to_toml(self, tmp_path: Path) -> None:
        subdir = tmp_path / "configs"
        subdir.mkdir()
        toml = subdir / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "../spec.json"\n'
            "[sweep.g0]\n"
            "start = 0.1\n"
            "stop = 1.0\n"
            "count = 3\n"
            "[output]\n"
            'path = "../results"\n',
        )
        config = load_sweep_config(toml)
        assert config.spec_path == subdir / ".." / "spec.json"
        assert config.output == subdir / ".." / "results"

    def test_unknown_sections_warn(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n'
            "[sweep.g0]\n"
            "start = 0.1\n"
            "stop = 1.0\n"
            "count = 3\n"
            "[foobar]\n"
            "x = 1\n",
        )
        config = load_sweep_config(toml)
        captured = capsys.readouterr()
        assert "foobar" in captured.err
        # Still loads successfully
        assert "g0" in config.swept_params

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_sweep_config(tmp_path / "nonexistent.toml")

    def test_explicit_values_in_toml(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text('spec = "spec.json"\n[sweep.g0]\nvalues = [0.1, 0.5, 0.9]\n')
        config = load_sweep_config(toml)
        assert config.swept_params["g0"] == [0.1, 0.5, 0.9]

    def test_multi_sweep_params(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n'
            "[sweep.g0]\n"
            "start = 0.1\n"
            "stop = 1.0\n"
            "count = 3\n"
            "[sweep.mChi2]\n"
            "start = 0.5\n"
            "stop = 4.0\n"
            "count = 4\n",
        )
        config = load_sweep_config(toml)
        assert len(config.swept_params) == 2
        assert len(config.swept_params["g0"]) == 3
        assert len(config.swept_params["mChi2"]) == 4

    def test_sweep_strategy_parsed(self, tmp_path: Path) -> None:
        toml = tmp_path / "sweep.toml"
        spec = tmp_path / "spec.json"
        spec.write_text("{}")
        toml.write_text(
            'spec = "spec.json"\n'
            "[sweep]\n"
            'strategy = "latin_hypercube"\n'
            "n_samples = 50\n"
            "[sweep.g0]\n"
            "start = 0.1\n"
            "stop = 1.0\n"
            "count = 3\n"
            "[sweep.mChi2]\n"
            "start = 0.5\n"
            "stop = 4.0\n"
            "count = 4\n",
        )
        config = load_sweep_config(toml)
        assert config.sweep_strategy == "latin_hypercube"
        assert config.n_samples == 50
        assert len(config.swept_params) == 2


# ------------------------------------------------------------------
# apply_config_to_args
# ------------------------------------------------------------------


class TestApplyConfigToArgs:
    """Test merging TOML config into argparse Namespace."""

    def _make_config(self, tmp_path: Path) -> SweepConfig:
        return SweepConfig(
            spec_path=tmp_path / "spec.json",
            swept_params={"g0": [0.1, 0.5, 1.0]},
            fixed_params={"mPhi2": 1.0, "mChi2": 4.0},
            sim_settings={"grid_shape": 128, "t_end": 20.0, "bounds": [0, 100]},
            measurements=["conversion", "mixing"],
            source="phi_0",
            target="chi_0",
            output=tmp_path / "sweep_out",
            parallel=4,
            resume=True,
            energy_threshold=1e-3,
            converge_sizes=None,
            force_large_sweep=False,
            dry_run=False,
        )

    def test_fills_empty_args(self, tmp_path: Path) -> None:
        config = self._make_config(tmp_path)
        args = argparse.Namespace(
            json_path=None,
            output=None,
            param=[],
            measure=None,
            source=None,
            target=None,
            energy_threshold=None,
            parallel=None,
            resume=False,
            force_large_sweep=False,
            dry_run=False,
            grid_shape=None,
            t_end=None,
            bounds=None,
        )
        swept, converge = apply_config_to_args(config, args)
        assert args.json_path == str(tmp_path / "spec.json")
        assert args.output == str(tmp_path / "sweep_out")
        assert "mPhi2=1.0" in args.param
        assert "mChi2=4.0" in args.param
        assert args.measure == "conversion,mixing"
        assert args.source == "phi_0"
        assert args.target == "chi_0"
        assert args.parallel == 4
        assert args.resume is True
        assert swept == {"g0": [0.1, 0.5, 1.0]}
        assert converge is None

    def test_cli_overrides_toml(self, tmp_path: Path) -> None:
        config = self._make_config(tmp_path)
        args = argparse.Namespace(
            json_path="other.json",
            output="/tmp/other",
            param=["mPhi2=2.0"],
            measure="conservation",
            source="chi_0",
            target="phi_0",
            energy_threshold=None,
            parallel=8,
            resume=False,
            force_large_sweep=False,
            dry_run=False,
            grid_shape=None,
            t_end=None,
            bounds=None,
        )
        apply_config_to_args(config, args)
        # CLI values preserved
        assert args.json_path == "other.json"
        assert args.output == "/tmp/other"
        assert args.measure == "conservation"
        assert args.source == "chi_0"
        assert args.parallel == 8
        # mPhi2 from CLI, mChi2 from TOML
        assert "mPhi2=2.0" in args.param
        assert "mChi2=4.0" in args.param
        # mPhi2 should NOT be duplicated from TOML
        assert sum(1 for p in args.param if p.startswith("mPhi2=")) == 1

    def test_sim_settings_convert_bounds(self, tmp_path: Path) -> None:
        config = self._make_config(tmp_path)
        args = argparse.Namespace(
            json_path=None,
            output=None,
            param=[],
            measure=None,
            source=None,
            target=None,
            energy_threshold=None,
            parallel=None,
            resume=False,
            force_large_sweep=False,
            dry_run=False,
            grid_shape=None,
            t_end=None,
            bounds=None,
        )
        apply_config_to_args(config, args)
        # bounds = [0, 100] should become "0:100"
        assert args.bounds == "0:100"
        # grid_shape = 128 should become "128"
        assert args.grid_shape == "128"
        assert args.t_end == 20.0
