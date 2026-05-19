"""Tests for ``tidal plot`` subcommand.

Tests cover CLI argument parsing, all 6 plot types, error handling,
and interaction with the SimulationData / disk-backed output pipeline.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import pytest

from tidal.cli import main

if TYPE_CHECKING:
    from pathlib import Path


# ============================================================
# Fixtures
# ============================================================


@pytest.fixture(scope="session")
def inline_2d_spec(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Write an inline 2+1D Klein-Gordon JSON spec."""
    import json

    spec: dict[str, Any] = {
        "metadata": {
            "source": "inline-test-2d",
            "lagrangian_expr": "-1/2 (d phi)^2",
            "derived_from": "Euler-Lagrange",
            "gauge": "none",
            "linearized": False,
            "parameters": {},
        },
        "spacetime": {
            "dimension": 3,
            "signature": [-1, 1, 1],
            "coordinates": ["t", "x", "y"],
        },
        "fields": [
            {"name": "phi_0", "index": 0, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            },
        ],
        "coupling": {},
    }
    d = tmp_path_factory.getbasetemp() / "inline_specs"
    d.mkdir(exist_ok=True)
    p = d / "kg_2d.json"
    if not p.exists():
        p.write_text(json.dumps(spec, indent=2))
    return p


@pytest.fixture(scope="session")
def two_d_sim_dir(
    inline_2d_spec: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Run a short 2D simulation for 2D plot tests."""
    output = tmp_path_factory.mktemp("plot2d") / "sim_out"
    ret = main(
        [
            "simulate",
            str(inline_2d_spec),
            "--t-end",
            "2.0",
            "--grid-shape",
            "16,16",
            "--bounds",
            "0:10,0:10",
            "--output",
            str(output),
            "--quiet",
        ],
    )
    assert ret == 0, "2D simulation failed"
    assert output.is_dir()
    return output


# ============================================================
# Basic CLI tests
# ============================================================


class TestPlotHelp:
    def test_plot_help(self) -> None:
        with pytest.raises(SystemExit, match="0"):
            main(["plot", "--help"])

    def test_plot_requires_type(self) -> None:
        """Missing --type should error."""
        with pytest.raises(SystemExit):
            main(["plot", "/nonexistent"])

    def test_plot_nonexistent_dir(self, tmp_path: Path) -> None:
        ret = main(["plot", str(tmp_path / "nope"), "--type", "heatmap"])
        assert ret == 1

    def test_plot_invalid_type(self, coupled_scalars_dir: Path) -> None:
        # argparse catches invalid choice and exits
        with pytest.raises(SystemExit):
            main(["plot", str(coupled_scalars_dir), "--type", "bogus"])


# ============================================================
# 1D plot types (using coupled_scalars_dir)
# ============================================================


class TestHeatmap:
    def test_heatmap_default_field(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "heatmap.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "heatmap",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_heatmap_specific_field(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "heatmap_chi.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "heatmap",
                "--field",
                "chi_0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_heatmap_custom_cmap(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "heatmap_viridis.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "heatmap",
                "--cmap",
                "viridis",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_heatmap_unknown_field(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "heatmap",
                "--field",
                "nonexistent",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1


class TestSnapshot1D:
    def test_snapshot_last(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "snapshot",
                "--time-index",
                "-1",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_snapshot_first(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap0.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "snapshot",
                "--time-index",
                "0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestAmplitude:
    def test_amplitude_all_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "amp.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_amplitude_single_field(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "amp_phi.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--fields",
                "phi_0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_amplitude_with_overlay(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "amp_overlay.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--overlay",
                "exp(-0.1*t)",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_amplitude_bad_overlay(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--overlay",
                "__import__('os').system('echo bad')",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1


class TestEnergy:
    def test_energy_all_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "energy.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "energy",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_energy_selected_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "energy_phi.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "energy",
                "--fields",
                "phi_0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestProfile:
    def test_profile_default_times(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "profile.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "profile",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_profile_specific_times(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "profile_times.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "profile",
                "--time-indices",
                "0,5,-1",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestCompare:
    def test_compare_all_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "compare.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "compare",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


# ============================================================
# 2D plot types
# ============================================================


class TestSnapshot2D:
    def test_snapshot_2d(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap2d.png"
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "snapshot",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_heatmap_rejects_2d(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        """Heatmap requires 1D spatial data — should error for 2D."""
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "heatmap",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_profile_2d_requires_cross_section(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        """Profile on 2D data without --cross-section should error."""
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "profile",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_profile_2d_with_cross_section(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "profile_cs.png"
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "profile",
                "--cross-section",
                "y=5.0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_compare_2d_requires_cross_section(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        """Compare on 2D data without --cross-section should error."""
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "compare",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_compare_2d_with_cross_section(
        self,
        two_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "compare_cs.png"
        ret = main(
            [
                "plot",
                str(two_d_sim_dir),
                "--type",
                "compare",
                "--cross-section",
                "x=5.0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


# ============================================================
# 3D plot types
# ============================================================


@pytest.fixture(scope="session")
def inline_3d_spec(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Write an inline 3+1D Klein-Gordon JSON spec."""
    import json

    spec: dict[str, Any] = {
        "metadata": {
            "source": "inline-test-3d",
            "lagrangian_expr": "-1/2 (d phi)^2",
            "derived_from": "Euler-Lagrange",
            "gauge": "none",
            "linearized": False,
            "parameters": {},
        },
        "spacetime": {
            "dimension": 4,
            "signature": [-1, 1, 1, 1],
            "coordinates": ["t", "x", "y", "z"],
        },
        "fields": [
            {"name": "phi_0", "index": 0, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            },
        ],
        "coupling": {},
    }
    d = tmp_path_factory.getbasetemp() / "inline_specs"
    d.mkdir(exist_ok=True)
    p = d / "kg_3d.json"
    if not p.exists():
        p.write_text(json.dumps(spec, indent=2))
    return p


@pytest.fixture(scope="session")
def three_d_sim_dir(
    inline_3d_spec: Path,
    tmp_path_factory: pytest.TempPathFactory,
) -> Path:
    """Run a short 3D simulation for 3D plot tests."""
    output = tmp_path_factory.mktemp("plot3d") / "sim_out"
    ret = main(
        [
            "simulate",
            str(inline_3d_spec),
            "--t-end",
            "1.0",
            "--grid-shape",
            "4,4,8",
            "--bounds",
            "0:4,0:4,0:8",
            "--output",
            str(output),
            "--quiet",
        ],
    )
    assert ret == 0, "3D simulation failed"
    assert output.is_dir()
    return output


class TestSnapshot3D:
    def test_snapshot_3d_initial(
        self,
        three_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap3d_0.png"
        ret = main(
            [
                "plot",
                str(three_d_sim_dir),
                "--type",
                "snapshot",
                "--time-index",
                "0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_snapshot_3d_final(
        self,
        three_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap3d_last.png"
        ret = main(
            [
                "plot",
                str(three_d_sim_dir),
                "--type",
                "snapshot",
                "--time-index",
                "-1",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_snapshot_3d_specific_field(
        self,
        three_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "snap3d_phi.png"
        ret = main(
            [
                "plot",
                str(three_d_sim_dir),
                "--type",
                "snapshot",
                "--field",
                "phi_0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestProfile3D:
    def test_profile_3d(
        self,
        three_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "profile3d.png"
        ret = main(
            [
                "plot",
                str(three_d_sim_dir),
                "--type",
                "profile",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestCompare3D:
    def test_compare_3d(
        self,
        three_d_sim_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "compare3d.png"
        ret = main(
            [
                "plot",
                str(three_d_sim_dir),
                "--type",
                "compare",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


class TestSingleFieldSelection:
    """Test that single_field() prefers dynamical fields over constraints."""

    def test_prefers_dynamical_field(self) -> None:
        from tidal.cli._panels import single_field

        class _MockData:
            fields = {"h_0": None, "h_1": None, "h_4": None}
            dynamical_fields = ("h_4",)

        result = single_field(_MockData(), None)  # type: ignore[arg-type]
        assert result == "h_4"

    def test_prefers_dynamical_even_with_constraint_velocities(self) -> None:
        """IDA writes constraint velocities — single_field should still pick dynamical."""
        from tidal.cli._panels import single_field

        class _MockData:
            fields = {"h_0": None, "h_1": None, "h_4": None}
            dynamical_fields = ("h_4",)
            velocities = {"h_0": None, "h_1": None, "h_4": None}

        result = single_field(_MockData(), None)  # type: ignore[arg-type]
        assert result == "h_4"

    def test_falls_back_to_first_if_all_constraints(self) -> None:
        from tidal.cli._panels import single_field

        class _MockData:
            fields = {"h_0": None, "h_1": None}
            dynamical_fields = ()

        result = single_field(_MockData(), None)  # type: ignore[arg-type]
        assert result == "h_0"

    def test_picks_largest_amplitude_among_dynamical(self) -> None:
        """When multiple dynamical fields exist, pick the one with largest amplitude."""
        import numpy as np

        from tidal.cli._panels import single_field

        class _MockData:
            fields = {
                "A_0": np.zeros((10, 32, 32)),  # constraint
                "A_1": np.ones((10, 32, 32)) * 1e-13,  # near-zero (longitudinal)
                "A_2": np.ones((10, 32, 32)) * 0.5,  # actual wave
            }
            dynamical_fields = ("A_1", "A_2")

        result = single_field(_MockData(), None)  # type: ignore[arg-type]
        assert result == "A_2"


# ============================================================
# Output options
# ============================================================


class TestOutputOptions:
    def test_custom_title(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "titled.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--title",
                "My Custom Title",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_custom_dpi(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "hidpi.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--dpi",
                "72",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_custom_figsize(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "sized.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--figsize",
                "10,4",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_default_output_path(
        self,
        coupled_scalars_dir: Path,
    ) -> None:
        """Without --output, file goes to DATA_DIR/{type}.png."""
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--quiet",
            ],
        )
        assert ret == 0
        expected = coupled_scalars_dir / "amplitude.png"
        assert expected.exists()
        # Clean up
        expected.unlink()

    def test_verbose_output(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Without --quiet, progress messages should print."""
        output = tmp_path / "verbose.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--output",
                str(output),
            ],
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "Loaded" in captured.out
        assert "Saved" in captured.out


# ============================================================
# Argument parsing edge cases
# ============================================================


class TestArgParsing:
    def test_bad_time_indices(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "profile",
                "--time-indices",
                "abc",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_bad_cross_section_format(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "profile",
                "--cross-section",
                "badformat",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_bad_cross_section_axis(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "profile",
                "--cross-section",
                "w=5.0",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_bad_figsize(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--figsize",
                "10,abc",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_bad_figsize_not_two(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--figsize",
                "10",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_unknown_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "amplitude",
                "--fields",
                "phi_0,nonexistent",
                "--output",
                str(tmp_path / "bad.png"),
                "--quiet",
            ],
        )
        assert ret == 1

    def test_param_override(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        """--param should be accepted without error."""
        output = tmp_path / "param.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "energy",
                "--param",
                "mPhi2=2.0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()


# ============================================================
# Unit tests for argument parsers and panel helpers
# ============================================================


class TestPanelHelpers:
    def test_parse_fields_none(self) -> None:
        from tidal.cli._plot_command import _parse_fields

        assert _parse_fields(None) is None

    def test_parse_fields_csv(self) -> None:
        from tidal.cli._plot_command import _parse_fields

        assert _parse_fields("a, b , c") == ["a", "b", "c"]

    def test_parse_time_indices_none(self) -> None:
        from tidal.cli._plot_command import _parse_time_indices

        assert _parse_time_indices(None) is None

    def test_parse_time_indices_valid(self) -> None:
        from tidal.cli._plot_command import _parse_time_indices

        assert _parse_time_indices("0, 5, -1") == [0, 5, -1]

    def test_parse_time_indices_invalid(self) -> None:
        from tidal.cli._plot_command import _parse_time_indices

        with pytest.raises(ValueError, match="Invalid --time-indices"):
            _parse_time_indices("abc")

    def test_parse_cross_section_none(self) -> None:
        from tidal.cli._plot_command import _parse_cross_section

        assert _parse_cross_section(None) is None

    def test_parse_cross_section_valid(self) -> None:
        from tidal.cli._plot_command import _parse_cross_section

        assert _parse_cross_section("y=25.0") == ("y", 25.0)

    def test_parse_cross_section_no_equals(self) -> None:
        from tidal.cli._plot_command import _parse_cross_section

        with pytest.raises(ValueError, match="expected AXIS=VAL"):
            _parse_cross_section("y25")

    def test_parse_cross_section_bad_axis(self) -> None:
        from tidal.cli._plot_command import _parse_cross_section

        with pytest.raises(ValueError, match="expected x, y, or z"):
            _parse_cross_section("w=5.0")

    def test_parse_figsize_valid(self) -> None:
        from tidal.cli._plot_command import _parse_figsize

        assert _parse_figsize("8,6") == (8.0, 6.0)

    def test_parse_figsize_invalid_count(self) -> None:
        from tidal.cli._plot_command import _parse_figsize

        with pytest.raises(ValueError, match="expected W,H"):
            _parse_figsize("8")

    def test_resolve_time_indices_default(self) -> None:
        from tidal.cli._panels import resolve_time_indices

        # Create a minimal mock
        class _MockData:
            n_snapshots = 10

        result = resolve_time_indices(_MockData(), None)  # type: ignore[arg-type]
        assert len(result) == 5
        assert result[0] == 0
        assert result[-1] == 9

    def test_resolve_time_indices_explicit(self) -> None:
        from tidal.cli._panels import resolve_time_indices

        class _MockData:
            n_snapshots = 10

        result = resolve_time_indices(_MockData(), [0, -1])  # type: ignore[arg-type]
        assert result == [0, 9]

    def test_resolve_time_indices_out_of_range(self) -> None:
        from tidal.cli._panels import resolve_time_indices

        class _MockData:
            n_snapshots = 5

        with pytest.raises(ValueError, match="out of range"):
            resolve_time_indices(_MockData(), [10])  # type: ignore[arg-type]


# ============================================================
# Hamiltonian and Conservation plot types
# ============================================================


class TestHamiltonian:
    def test_hamiltonian_all_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "hamiltonian.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "hamiltonian",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_hamiltonian_selected_fields(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "hamiltonian_phi.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "hamiltonian",
                "--fields",
                "phi_0",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_hamiltonian_default_filename(
        self,
        coupled_scalars_dir: Path,
    ) -> None:
        """Without --output, saves to DATA_DIR/hamiltonian.png."""
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "hamiltonian",
                "--quiet",
            ],
        )
        assert ret == 0
        expected = coupled_scalars_dir / "hamiltonian.png"
        assert expected.exists()
        expected.unlink()


class TestConservation:
    def test_conservation_default_threshold(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "conservation.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "conservation",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_conservation_custom_threshold(
        self,
        coupled_scalars_dir: Path,
        tmp_path: Path,
    ) -> None:
        output = tmp_path / "conservation_tight.png"
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "conservation",
                "--threshold",
                "1e-6",
                "--output",
                str(output),
                "--quiet",
            ],
        )
        assert ret == 0
        assert output.exists()

    def test_conservation_default_filename(
        self,
        coupled_scalars_dir: Path,
    ) -> None:
        """Without --output, saves to DATA_DIR/conservation.png."""
        ret = main(
            [
                "plot",
                str(coupled_scalars_dir),
                "--type",
                "conservation",
                "--quiet",
            ],
        )
        assert ret == 0
        expected = coupled_scalars_dir / "conservation.png"
        assert expected.exists()
        expected.unlink()


# ============================================================
# Panel visibility tests (position-dependent error sentinel)
# ============================================================


class TestMeasurePlotEmptyPanels:
    """Verify that panels with error sentinels are hidden rather than showing
    empty 'No X data' text.
    """

    def _make_ax(self) -> Any:  # noqa: ANN401
        import matplotlib.pyplot as plt

        _, ax = plt.subplots()
        return ax

    def test_spectrum_panel_hidden_on_error(self) -> None:
        from tidal.cli._measure_plot import _plot_spectrum

        ax = self._make_ax()
        results: dict[str, Any] = {
            "spectrum": {"error": "position-dependent term detected"},
        }
        _plot_spectrum(ax, results)
        assert not ax.get_visible()

    def test_spectrum_panel_hidden_when_missing(self) -> None:
        from tidal.cli._measure_plot import _plot_spectrum

        ax = self._make_ax()
        _plot_spectrum(ax, {})
        assert not ax.get_visible()

    def test_dispersion_panel_hidden_on_error(self) -> None:
        from tidal.cli._measure_plot import _plot_dispersion

        ax = self._make_ax()
        results: dict[str, Any] = {
            "dispersion": {"error": "requires spatially uniform system"},
        }
        _plot_dispersion(ax, results)
        assert not ax.get_visible()

    def test_spectrum_panel_visible_with_valid_data(self) -> None:
        """When data is present and correct, the panel should be visible."""
        import numpy as np

        from tidal.cli._measure_plot import _plot_spectrum

        ax = self._make_ax()
        wn = np.linspace(0.0, 1.0, 5)
        power = np.ones(5)
        results: dict[str, Any] = {
            "spectrum": {
                "phi_0": {
                    "initial": {"wavenumbers": wn, "power": power},
                    "final": {"wavenumbers": wn, "power": power},
                },
            },
        }
        _plot_spectrum(ax, results)
        assert ax.get_visible()
