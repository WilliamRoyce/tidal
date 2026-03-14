"""Tests for the ``tidal`` command-line interface."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import ClassVar

import pytest

from tidal.cli import main


class TestMainEntryPoint:
    def test_no_args_returns_zero(self) -> None:
        assert main([]) == 0

    def test_help_exits_zero(self) -> None:
        with pytest.raises(SystemExit, match="0"):
            main(["--help"])

    def test_version(self) -> None:
        with pytest.raises(SystemExit, match="0"):
            main(["--version"])

    def test_module_invocation(self) -> None:
        """``python -m tidal.cli`` should be importable."""
        import tidal.cli.__main__  # noqa: F401  # pyright: ignore[reportUnusedImport]

    def test_get_version_fallback(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """_get_version() should return 'unknown' when package is not installed."""
        from importlib.metadata import PackageNotFoundError

        from tidal.cli import _get_version

        def _raise(_name: str) -> str:
            raise PackageNotFoundError

        monkeypatch.setattr("importlib.metadata.version", _raise)
        assert _get_version() == "unknown"

    def test_entry_point_subprocess(self) -> None:
        """The ``tidal`` entry point should be callable as a subprocess."""
        result = subprocess.run(
            [sys.executable, "-m", "tidal.cli", "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0
        assert "tidal" in result.stdout


class TestInspectCommand:
    @pytest.mark.parametrize(
        ("fixture_name", "dim_label", "fields", "expected_strs"),
        [
            ("inline_kg_1d_json", "1+1D", 1, ["phi_0"]),
            ("massive_3form_json", "3+1D", 4, ["C_0", "C_3", "m2"]),
            ("inline_coupled_scalars_json", "1+1D", 2, ["phi_0", "chi_0"]),
        ],
    )
    def test_inspect_specs(
        self,
        fixture_name: str,
        dim_label: str,
        fields: int,
        expected_strs: list[str],
        capsys: pytest.CaptureFixture[str],
        request: pytest.FixtureRequest,
    ) -> None:
        json_path: Path = request.getfixturevalue(fixture_name)
        ret = main(["inspect", str(json_path)])
        assert ret == 0

        out = capsys.readouterr().out
        assert dim_label in out
        if fields > 1:
            assert f"{fields} components" in out
        for s in expected_strs:
            assert s in out

    def test_inspect_nonexistent_file(self) -> None:
        ret = main(["inspect", "/nonexistent/file.json"])
        assert ret == 1

    def test_inspect_with_params_flag(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(["inspect", str(inline_kg_1d_json), "--params"])
        assert ret == 0

    def test_inspect_json_output(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--json should output valid, parseable JSON."""
        import json

        ret = main(["inspect", str(inline_kg_1d_json), "--json"])
        assert ret == 0

        out = capsys.readouterr().out
        data = json.loads(out)
        assert data["spacetime"]["dimension"] == 2
        assert data["spacetime"]["spatial_dimension"] == 1
        assert "phi_0" in data["fields"]
        assert len(data["equations"]) == 1
        assert "mass_matrix" in data
        assert "coupling_matrix" in data

    def test_inspect_json_with_params(
        self, inline_coupled_scalars_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--json --params should include default parameter values."""
        import json

        ret = main(["inspect", str(inline_coupled_scalars_json), "--json", "--params"])
        assert ret == 0

        data = json.loads(capsys.readouterr().out)
        assert "required" in data["parameters"]
        assert len(data["fields"]) == 2


class TestListCommand:
    def test_list_default_dir(self, capsys: pytest.CaptureFixture[str]) -> None:
        ret = main(["list"])
        # Returns 0 even with empty examples/data/
        assert ret == 0

    def test_list_custom_dir(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(["list", "--dir", str(inline_kg_1d_json.parent)])
        assert ret == 0

        out = capsys.readouterr().out
        assert "klein_gordon_1d.json" in out
        assert "specifications found" in out

    def test_list_nonexistent_dir(self) -> None:
        ret = main(["list", "--dir", "/nonexistent/dir"])
        assert ret == 1


class TestSimulateCommand:
    def test_simulate_1d_summary(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "1.0",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "phi_0" in out

    def test_simulate_with_params(
        self, klein_gordon_3d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(klein_gordon_3d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--grid-shape",
                "8",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_png_output(self, inline_kg_1d_json: Path, tmp_path: Path) -> None:
        output = tmp_path / "test_output.png"
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--output",
                str(output),
            ]
        )
        assert ret == 0
        assert output.exists()

    def test_simulate_2d_plot_output(
        self,
        chern_simons_json: Path,
        tmp_path: Path,
    ) -> None:
        """2D spec should produce a non-empty PNG file (exercises plot_2d path)."""
        output = tmp_path / "cs_2d.png"
        ret = main(
            [
                "simulate",
                str(chern_simons_json),
                "--grid-shape",
                "8",
                "--t-end",
                "0.2",
                "--output",
                str(output),
            ]
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    # --- IC types (parametrized) ---

    @pytest.mark.parametrize("ic_type", ["gaussian", "plane-wave", "zero"])
    def test_simulate_ic_types(
        self,
        inline_kg_1d_json: Path,
        ic_type: str,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                ic_type,
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_gaussian_custom(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "gaussian",
                "--ic-width",
                "1.5",
                "--ic-amplitude",
                "2.0",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_off_center_gaussian(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "gaussian",
                "--ic-center",
                "30.0",
                "--ic-width",
                "3.0",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_chern_simons_constraint(
        self, chern_simons_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Test simulation of a system with constraint (time_order=0) + dynamical fields."""
        ret = main(
            [
                "simulate",
                str(chern_simons_json),
                "--t-end",
                "0.5",
                "--grid-shape",
                "8",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "A_0" in out

    def test_simulate_invalid_param_format(self, inline_kg_1d_json: Path) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "bad_no_equals",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_invalid_ic_component(self, inline_kg_1d_json: Path) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--ic-component",
                "nonexistent_field",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_nonexistent_file(self) -> None:
        ret = main(["simulate", "/nonexistent/file.json", "--no-plot"])
        assert ret == 1

    def test_simulate_custom_grid(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--grid-shape",
                "32",
                "--bounds",
                "0:20",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    # --- Feature: --bc (mixed boundary conditions) ---

    def test_simulate_bc_single(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--bc",
                "neumann",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_bc_mixed_2d(
        self, polar_kg_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(polar_kg_json),
                "--bc",
                "neumann,periodic",
                "--grid-shape",
                "16",
                "--bounds",
                "0.5:8,0:6.28",
                "--ic",
                "gaussian",
                "--ic-center",
                "3.0,3.14",
                "--ic-width",
                "0.5",
                "--t-end",
                "0.5",
                "--dt",
                "0.01",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_bc_invalid_type(self, inline_kg_1d_json: Path) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--bc",
                "invalid_bc",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_bc_wrong_count(self, inline_kg_1d_json: Path) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--bc",
                "neumann,periodic",  # 2 values for 1D
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_bc_dirichlet_accepted(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Dirichlet BC is supported by the native solver."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--bc",
                "dirichlet",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    # --- Feature: --ic formula ---

    def test_simulate_formula_ic(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "formula",
                "--ic-formula",
                "np.exp(-((x - 5)**2) / 2)",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_formula_ic_2d(
        self, polar_kg_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(polar_kg_json),
                "--ic",
                "formula",
                "--ic-formula",
                "np.exp(-((x - 3)**2 + (y - pi)**2) / 0.5**2)",
                "--grid-shape",
                "16",
                "--bounds",
                "0.5:8,0:6.28",
                "--bc",
                "neumann,periodic",
                "--t-end",
                "0.5",
                "--dt",
                "0.01",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_formula_ic_missing_expr(self, inline_kg_1d_json: Path) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--ic",
                "formula",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_formula_ic_constant(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "formula",
                "--ic-formula",
                "0.5",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_formula_ic_undefined_var(self, inline_kg_1d_json: Path) -> None:
        """Formula with undefined variable should fail."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--ic",
                "formula",
                "--ic-formula",
                "badvar * 2",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_formula_ic_attribute_access_rejected(
        self,
        inline_kg_1d_json: Path,
    ) -> None:
        """Formula with attribute access should be rejected by AST validator."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--ic",
                "formula",
                "--ic-formula",
                "x.__class__.__name__",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    # --- Feature: --mode constraint ---

    def test_simulate_constraint_mode(
        self, inline_constraint_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_constraint_json),
                "--mode",
                "constraint",
                "--grid-shape",
                "16",
                "--bc",
                "neumann",
                "--ic",
                "gaussian",
                "--ic-component",
                "rho",
                "--no-plot",
            ]
        )
        assert ret == 0

        out = capsys.readouterr().out
        assert "Constraint solve complete" in out

    # --- Solver options ---

    def test_simulate_explicit_dt(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--dt",
                "0.01",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_default_scheme(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify default scheme (auto) auto-selects and logs solver choice."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        # Inline KG spec lacks canonical section → auto selects IDA
        assert "Auto-selected solver:" in captured.out
        assert "Scheme:" in captured.out

    def test_auto_selects_modal_for_periodic_constraints(
        self, inline_em_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Periodic constraint equations should auto-select modal (Schur)."""
        ret = main(
            [
                "simulate",
                str(inline_em_1d_json),
                "--ic",
                "plane-wave",
                "--ic-component",
                "A_2",
                "--t-end",
                "0.1",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "Auto-selected solver: modal" in captured.out

    def test_simulate_ida_scheme(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify IDA solver runs end-to-end on Klein-Gordon."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--scheme",
                "ida",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "IDA solver" in captured.out
        assert "snapshots stored" in captured.out

    def test_simulate_ida_plane_wave(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """IDA with plane-wave IC should preserve amplitude."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--scheme",
                "ida",
                "--ic",
                "plane-wave",
                "--t-end",
                "1.0",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_leapfrog_scheme(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify leapfrog solver runs end-to-end on Klein-Gordon."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--scheme",
                "leapfrog",
                "--t-end",
                "0.5",
                "--dt",
                "0.01",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "leapfrog" in captured.out.lower()
        assert "snapshots stored" in captured.out

    def test_simulate_custom_snapshots(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--snapshots",
                "0.25",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    # --- 3D simulation ---

    def test_simulate_3d(
        self, klein_gordon_3d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(
            [
                "simulate",
                str(klein_gordon_3d_json),
                "--param",
                "m2=1.0",
                "--grid-shape",
                "8",
                "--t-end",
                "0.2",
                "--no-plot",
            ]
        )
        assert ret == 0

    # --- Edge-case validation ---

    def test_simulate_explicit_format_flag(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Explicit --format flag should work for all formats."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--format",
                "summary",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_periodic_flag(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--no-periodic should produce non-periodic (Neumann) BCs."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--no-periodic",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_wavevector_custom(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--ic-wavevector should override default wavevector for plane-wave."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "plane-wave",
                "--ic-wavevector",
                "0.5",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0

    def test_simulate_ic_center_wrong_dim(self, inline_kg_1d_json: Path) -> None:
        """--ic-center with wrong dimension count should fail."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic-center",
                "5.0,5.0",  # 2 values for 1D
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_snapshots_nonpositive(self, inline_kg_1d_json: Path) -> None:
        """--snapshots with zero or negative value should fail."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--snapshots",
                "0",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_dt_nonpositive(self, inline_kg_1d_json: Path) -> None:
        """--dt with zero or negative value should fail."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--dt",
                "-0.1",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 1

    def test_simulate_t_end_nonpositive(self, inline_kg_1d_json: Path) -> None:
        """--t-end with zero or negative value should fail."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0",
                "--no-plot",
            ]
        )
        assert ret == 1

    # --- Feature: --quiet ---

    def test_simulate_quiet_suppresses_progress(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--quiet should suppress progress messages but keep results."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--no-plot",
                "--quiet",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        # Results should still appear
        assert "Results:" in out
        # Progress messages should NOT appear
        assert "Loading equation specification" not in out
        assert "Building PDE" not in out
        assert "Running simulation" not in out

    def test_simulate_quiet_short_flag(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """-q should work as shorthand for --quiet."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--no-plot",
                "-q",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Results:" in out
        assert "Loading equation specification" not in out


class TestZeroEvolutionWarning:
    """Tests for the zero-evolution diagnostic warning."""

    def test_em_plane_wave_no_warning(
        self, inline_em_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Plane-wave IC on EM A_2 (transverse) provides non-zero π → no warning."""
        ret = main(
            [
                "simulate",
                str(inline_em_1d_json),
                "--ic",
                "plane-wave",
                "--ic-component",
                "A_2",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "all evolution rates are zero" not in captured.err

    def test_kg_gaussian_no_warning(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """KG Gaussian IC has non-zero dπ/dt from laplacian → no warning."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--ic",
                "gaussian",
                "--t-end",
                "0.5",
                "--no-plot",
            ]
        )
        assert ret == 0
        captured = capsys.readouterr()
        assert "all evolution rates are zero" not in captured.err

    def test_em_plane_wave_amplitude_stable(
        self,
        inline_em_1d_json: Path,
        tmp_path: Path,
    ) -> None:
        """EM plane-wave in A_2 (transverse): amplitude stays bounded (no growth)."""
        import numpy as np

        out_dir = tmp_path / "em_stable"
        ret = main(
            [
                "simulate",
                str(inline_em_1d_json),
                "--ic",
                "plane-wave",
                "--ic-component",
                "A_2",
                "--t-end",
                "5.0",
                "--scheme",
                "ida",
                "--periodic",
                "--output",
                str(out_dir),
            ]
        )
        assert ret == 0

        # Read final snapshot from disk-backed storage
        a2_data = np.load(str(out_dir / "A_2.npy"), mmap_mode="r")
        final_peak = float(np.max(np.abs(a2_data[-1])))
        initial_peak = 1.0  # plane-wave amplitude
        # Amplitude should stay within 20% of initial (generous tolerance for numerics)
        assert final_peak < 1.2 * initial_peak, (
            f"A_2 amplitude grew from {initial_peak:.3f} to {final_peak:.3f} — "
            f"constraint solver may not be active"
        )


class TestRequireStable:
    """Tests for the --require-stable flag end-to-end."""

    def test_require_stable_aborts_on_unstable_spec(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Unstable coupled spec + --require-stable should sys.exit(1)."""
        import json

        # g0^2 = 25 > mPhi2 * mChi2 = 0.01  =>  strongly unstable
        spec_data = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [
                {"name": "phi_0", "index": 0, "is_dynamical": True},
                {"name": "chi_0", "index": 1, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -0.1,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": -5.0,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -0.1,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                            {
                                "coefficient": -5.0,
                                "operator": "identity",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "chi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec_path = tmp_path / "unstable.json"
        spec_path.write_text(json.dumps(spec_data, indent=2))

        with pytest.raises(SystemExit) as exc_info:
            main(
                [
                    "simulate",
                    str(spec_path),
                    "--t-end",
                    "0.5",
                    "--no-plot",
                    "--require-stable",
                ]
            )
        assert exc_info.value.code == 1

        err = capsys.readouterr().err
        assert "eigenvalue" in err.lower()

    def test_stable_spec_passes_with_require_stable(
        self, inline_kg_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Stable spec + --require-stable should complete successfully."""
        ret = main(
            [
                "simulate",
                str(inline_kg_1d_json),
                "--param",
                "m2=1.0",
                "--t-end",
                "0.5",
                "--no-plot",
                "--require-stable",
            ]
        )
        assert ret == 0
        err = capsys.readouterr().err
        assert "eigenvalue" not in err.lower()


class TestDeriveCommand:
    def test_derive_toml_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Test Scalar"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - m2/2 phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "wolframscript" in out
        assert "DefManifold" in out
        assert "DefConstantSymbol[m2]" in out
        assert "EulerLagrangeEquation" in out

    def test_derive_toml_dry_run_with_parameters(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Dry-run with [parameters] should inject metadata into generated WLS."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Test With Params"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - m2/2 phi[]^2"

[parameters]
m2 = 1.0

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0
        out = capsys.readouterr().out
        assert '"m2" -> 1.0' in out
        assert "parameters" in out.lower()

    def test_derive_toml_save_script(self, tmp_path: Path) -> None:
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Test Vector"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/4 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]] - m2/2 A[-a] eta[a,b] A[-b]"

[output]
path = "output.json"
""")
        script_path = tmp_path / "generated.wls"
        # Just save the script, don't run (wolframscript may not be available)
        main(["derive", str(config), "--save-script", str(script_path)])
        # May return 1 if wolframscript not found, but script should be saved
        assert script_path.exists()

        content = script_path.read_text()
        assert "DefTensor" in content
        assert "tsA" in content or "tvA" in content  # Prefixed field name

    def test_derive_nonexistent_file(self) -> None:
        ret = main(["derive", "/nonexistent/theory.toml"])
        assert ret == 1

    def test_derive_bad_extension(self) -> None:
        ret = main(["derive", "/tmp/test.xyz"])
        assert ret == 1

    def test_derive_multi_field_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        config = tmp_path / "coupled.toml"
        config.write_text("""
[theory]
name = "Coupled Scalars"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[fields]]
name = "chi"
type = "scalar"

[constants]
names = ["g"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - 1/2 CD[-a][chi[]] eta[a,b] CD[-b][chi[]] - g phi[] chi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Multi-field should use VarD per field
        assert "VarD" in out
        assert "BuildMultiFieldJSONStructure" in out

    def test_derive_tensor_field_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        config = tmp_path / "tensor.toml"
        config.write_text("""
[theory]
name = "Massive Three Form"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "C"
type = "tensor"
rank = 3
symmetry = "antisymmetric"

[constants]
names = ["m2"]

[lagrangian]
expression = "-1/12 CD[-d][C[-a, -b, -c]] eta[d, h] eta[a, e] eta[b, f] eta[c, g] CD[-h][C[-e, -f, -g]] - m2/12 C[-a, -b, -c] eta[a, e] eta[b, f] eta[c, g] C[-e, -f, -g]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Antisymmetric" in out
        assert "DefConstantSymbol[m2]" in out

    def test_derive_single_field_passes_metric_matrix(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Single-field DecomposeToComponents must include MetricMatrix option."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Test Scalar"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert '"MetricMatrix"' in out
        assert "tsMetricMatrix" in out

    def test_derive_multi_field_passes_metric_matrix(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Multi-field DecomposeToComponents must include MetricMatrix option."""
        config = tmp_path / "coupled.toml"
        config.write_text("""
[theory]
name = "Coupled Scalars"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[fields]]
name = "chi"
type = "scalar"

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - 1/2 CD[-a][chi[]] eta[a,b] CD[-b][chi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert '"MetricMatrix"' in out
        assert "csMetricMatrix" in out

    def test_derive_curvilinear_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Curvilinear metric with coordinate-dependent entries should generate valid WLS."""
        config = tmp_path / "polar.toml"
        config.write_text("""
[theory]
name = "Polar Klein-Gordon"

[spacetime]
dimension = 3
metric = "diagonal"
diagonal = [-1, 1, "x[]^2"]

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["polm2"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - polm2/2 phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Should have coordinate-dependent metric entry
        assert "x[]^2" in out
        # Should use DiagonalMatrix
        assert "DiagonalMatrix" in out
        # Should pass MetricMatrix to DecomposeToComponents
        assert '"MetricMatrix"' in out
        assert "DefConstantSymbol[polm2]" in out

    def test_derive_chart_placeholder_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Component-derivative notation with -chart placeholder should be substituted."""
        config = tmp_path / "elasticity.toml"
        config.write_text("""
[theory]
name = "Navier Cauchy"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "ux"
type = "scalar"

[[fields]]
name = "uy"
type = "scalar"

[constants]
names = ["rho", "lam", "mu"]

[lagrangian]
expression = "rho/2 * (CD[{0, -chart}][ux[]]^2 + CD[{0, -chart}][uy[]]^2) - lam/2 * (CD[{1, -chart}][ux[]] + CD[{2, -chart}][uy[]])^2 - mu * (CD[{1, -chart}][ux[]]^2 + CD[{2, -chart}][uy[]]^2 + 1/2 * (CD[{2, -chart}][ux[]] + CD[{1, -chart}][uy[]])^2)"

[parameters]
rho = 1.0
lam = 1.0
mu = 1.0

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Lagrangian assignment should have chart placeholder substituted
        assert "ncCD[{0, -ncCart}][ncUx[]]" in out
        assert "ncCD[{1, -ncCart}][ncUx[]]" in out
        # Metadata lagrangian_expr stores original expression (unsubstituted) — that's OK
        # Field names should be prefixed
        assert "ncUx[]" in out
        assert "ncUy[]" in out
        # Constants should be defined
        assert "DefConstantSymbol[rho]" in out
        assert "DefConstantSymbol[lam]" in out
        assert "DefConstantSymbol[mu]" in out
        # Multi-field: should have VarD for both fields
        assert "VarD[ncUx[]" in out
        assert "VarD[ncUy[]" in out

    def test_derive_linearization_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Linearization (xPert) TOML should generate valid WLS."""
        config = tmp_path / "gw.toml"
        config.write_text("""
[theory]
name = "Linearized Gravity"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[linearization]
expression = "Einstein[CD][-a, -b]"
perturbation_field = "h"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Should load xPert
        assert "xAct`xPert`" in out
        assert "Linearize.wl" in out
        # Should have SetupMetricPerturbation
        assert "SetupMetricPerturbation" in out
        # Should have LinearizeTensorExpression
        assert "LinearizeTensorExpression" in out
        # Should have notation conversion rule
        assert "LI[1]" in out
        # Should have DecomposeToComponents
        assert "DecomposeToComponents" in out
        # CD in Einstein should be prefixed
        assert "Einstein[lgCD]" in out
        # Should NOT have VarD or Lagrangian
        assert "VarD" not in out
        # Metadata should indicate linearized
        assert '"linearized" -> True' in out

    def test_derive_massive_gravity_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Fierz-Pauli linearization with bg tensor generates valid WLS."""
        config = tmp_path / "mg.toml"
        config.write_text("""
[theory]
name = "Massive Gravity"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[constants]
names = ["m2"]

[linearization]
expression = "Einstein[CD][-a, -b] - m2 (eta[-a, -b] - bg[-a, -b]) + m2 eta[-a, -b] eta[c, d] (eta[-c, -d] - bg[-c, -d])"
perturbation_field = "h"

[parameters]
m2 = 1.0

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Should load xPert + Linearize.wl
        assert "xAct`xPert`" in out
        assert "Linearize.wl" in out
        # Constant should be defined
        assert "DefConstantSymbol[m2]" in out
        # Background metric tensor defined (not perturbed by xPert)
        assert "DefTensor[mgBg[-a, -b]" in out
        # Explicit zero perturbation rule for bg (Unprotect/Protect required)
        assert "Unprotect[Perturbation]" in out
        assert "Perturbation[mgBg[__], ___] := 0" in out
        assert "Protect[Perturbation]" in out
        # xPert setup
        assert "SetupMetricPerturbation" in out
        assert "LinearizeTensorExpression" in out
        # Metric reference should be substituted (eta -> mgEta, bg -> mgBg)
        assert "mgEta[-a, -b]" in out
        assert "mgBg[-a, -b]" in out
        # bg -> metric replacement before decomposition
        assert "linExprPlain = linExprPlain /. mgBg -> mgEta" in out
        # Notation conversion
        assert "LI[1]" in out
        # Decomposition
        assert "DecomposeToComponents" in out
        # Metadata
        assert '"linearized" -> True' in out
        # No Euler-Lagrange path
        assert "VarD" not in out
        # Parameter default value in metadata
        assert '"m2" -> 1.0' in out


class TestGaugeFixing:
    """Tests for ``[[gauge]]`` TOML → validation and WLS generation."""

    _BASE_TOML = """\
[theory]
name = "EM"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"

[output]
path = "output.json"
"""

    def test_gauge_lorenz_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Lorenz gauge on vector field → GaugeFix.wl loaded + builder in WLS."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "lorenz"
xi = 1.0
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "GaugeFix.wl" in out
        assert "BuildLorenzGaugeTerm" in out
        assert "AddGaugeFixingTerm" in out

    def test_gauge_custom_lagrangian_term_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Custom gauge expression → appears in WLS output."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
expression = "-(1/2) * eta[a,b] CD[-a][A[-c]] eta[c,d] CD[-b][A[-d]]"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "GaugeTerm" in out
        assert "AddGaugeFixingTerm" in out

    def test_gauge_metadata_lorenz(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Lorenz gauge metadata string in WLS."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "lorenz"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0
        out = capsys.readouterr().out
        assert '"gauge" -> "lorenz(A)"' in out

    def test_gauge_metadata_none_default(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """No [[gauge]] → 'gauge -> none' in metadata."""
        config = tmp_path / "theory.toml"
        config.write_text(self._BASE_TOML)
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0
        out = capsys.readouterr().out
        assert '"gauge" -> "none"' in out

    def test_gauge_bad_field_ref(self, tmp_path: Path) -> None:
        """Reject field='Z' not declared in [[fields]]."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "Z"
type = "lorenz"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_gauge_bad_type(self, tmp_path: Path) -> None:
        """Reject unknown gauge type."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "unknown"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_gauge_scalar_not_allowed(self, tmp_path: Path) -> None:
        """Lorenz gauge on scalar field → exit code 1."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "phi[]^2"

[output]
path = "output.json"

[[gauge]]
field = "phi"
type = "lorenz"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_gauge_duplicate_field(self, tmp_path: Path) -> None:
        """Two gauge entries for same field → exit code 1."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "lorenz"

[[gauge]]
field = "A"
type = "lorenz"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_gauge_custom_missing_expression(self, tmp_path: Path) -> None:
        """Custom gauge without expression key is rejected."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "custom"
mechanism = "lagrangian_term"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_gauge_custom_missing_mechanism(self, tmp_path: Path) -> None:
        """Custom gauge without mechanism key is rejected."""
        config = tmp_path / "theory.toml"
        config.write_text(
            self._BASE_TOML
            + """
[[gauge]]
field = "A"
type = "custom"
expression = "eta[a,b] CD[-a][A[-b]]"
"""
        )
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1


class TestDeriveAbsolutePaths:
    """Verify generated WLS scripts use absolute paths (not $InputFileName-relative)."""

    _MINIMAL_CONFIG: ClassVar[dict[str, object]] = {
        "theory": {"name": "Test Scalar"},
        "spacetime": {"dimension": 2, "metric": "minkowski"},
        "fields": [{"name": "phi", "type": "scalar"}],
        "lagrangian": {"expression": "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"},
        "output": {"path": "output.json"},
    }

    def test_pipeline_path_is_absolute(self) -> None:
        from tidal.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG)
        assert "$InputFileName" not in script
        assert "/tidal/wolfram" in script

    def test_output_path_is_absolute(self) -> None:
        from tidal.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG, config_dir=Path("/project"))
        # The output path "output.json" should be resolved to /project/output.json
        assert '"/project/output.json"' in script

    def test_output_override_absolute(self) -> None:
        from tidal.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG, output_override="/tmp/test.json")
        assert '"/tmp/test.json"' in script

    def test_output_path_resolves_relative_to_config_dir(self, tmp_path: Path) -> None:
        from tidal.cli._derive import generate_wls

        # Create config_dir = tmp_path/sub/
        config_dir = tmp_path / "sub"
        config_dir.mkdir()

        config = {
            **self._MINIMAL_CONFIG,
            "output": {"path": "../data/out.json"},
        }
        script = generate_wls(config, config_dir=config_dir)
        # ../data/out.json relative to tmp_path/sub → tmp_path/data/out.json
        expected = str(tmp_path / "data" / "out.json")
        assert f'"{expected}"' in script


class TestDeriveValidation:
    def test_missing_spacetime(self, tmp_path: Path) -> None:
        config = tmp_path / "bad.toml"
        config.write_text("""
[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "phi[]^2"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_missing_fields(self, tmp_path: Path) -> None:
        config = tmp_path / "bad.toml"
        config.write_text("""
[spacetime]
dimension = 2

[lagrangian]
expression = "phi[]^2"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_bad_field_name(self, tmp_path: Path) -> None:
        config = tmp_path / "bad.toml"
        config.write_text("""
[spacetime]
dimension = 2

[[fields]]
name = "phi-1"
type = "scalar"

[lagrangian]
expression = "phi-1[]^2"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_tensor_rank_exceeds_dimension(self, tmp_path: Path) -> None:
        config = tmp_path / "bad.toml"
        config.write_text("""
[spacetime]
dimension = 2

[[fields]]
name = "C"
type = "tensor"
rank = 5
symmetry = "antisymmetric"

[lagrangian]
expression = "C[-a,-b,-c,-d,-e] C[-a,-b,-c,-d,-e]"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_bad_dimension(self, tmp_path: Path) -> None:
        config = tmp_path / "bad.toml"
        config.write_text("""
[spacetime]
dimension = 99

[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "phi[]^2"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1


class TestValidateCommand:
    """Tests for ``tidal validate``."""

    def test_validate_valid_spec(
        self,
        inline_kg_1d_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        ret = main(["validate", str(inline_kg_1d_json)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "OK" in out

    def test_validate_nonexistent_file(self) -> None:
        ret = main(["validate", "nonexistent_file.json"])
        assert ret == 1

    def test_validate_invalid_json(self, tmp_path: Path) -> None:
        bad_json = tmp_path / "bad.json"
        bad_json.write_text("{ not valid json }", encoding="utf-8")
        ret = main(["validate", str(bad_json)])
        assert ret == 1

    def test_validate_warns_on_missing_param_defaults(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """A spec with symbolic coefficient but no parameter defaults should warn."""
        import json

        spec = {
            "metadata": {"source": "test"},
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-m2",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "phi_0",
                            },
                        ],
                    },
                }
            ],
            "coupling": {"mass_matrix": [[1.0]], "coupling_matrix": [[0.0]]},
        }
        spec_path = tmp_path / "no_defaults.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")
        main(["validate", str(spec_path)])
        err = capsys.readouterr().err
        assert "WARNING" in err
        assert "m2" in err

    def test_validate_coupled_scalars(
        self,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Coupled scalars spec should validate successfully."""
        ret = main(["validate", str(inline_coupled_scalars_json)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "OK" in out

    def test_validate_unknown_operator(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """JSON with unknown operator should fail validation."""
        import json

        spec = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "operator": "nonexistent_op",
                                "field": "phi_0",
                                "coefficient": 1.0,
                            }
                        ],
                    },
                }
            ],
        }
        spec_path = tmp_path / "bad_op.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")
        ret = main(["validate", str(spec_path)])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Unknown operator" in err
        assert "nonexistent_op" in err

    def test_validate_bad_field_reference(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """JSON with unknown field reference should fail validation."""
        import json

        spec = {
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2, "space": 0},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "operator": "identity",
                                "field": "chi_99",
                                "coefficient": 1.0,
                            }
                        ],
                    },
                }
            ],
        }
        spec_path = tmp_path / "bad_ref.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")
        ret = main(["validate", str(spec_path)])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Unknown field reference" in err
        assert "chi_99" in err

    def test_validate_directory_instead_of_file(
        self,
        tmp_path: Path,
    ) -> None:
        """Passing a directory path should fail validation."""
        ret = main(["validate", str(tmp_path)])
        assert ret == 1

    def test_validate_stability_stable_spec(
        self,
        inline_kg_1d_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--stability on a stable spec should return 0 and say [stable]."""
        ret = main(["validate", str(inline_kg_1d_json), "--stability", "--param", "m2=1.0"])
        assert ret == 0
        out = capsys.readouterr().out
        assert "stable" in out.lower()

    def test_validate_stability_tachyon_detection(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--stability should detect a tachyonic mode (negative mass²)."""
        import json

        spec = {
            "metadata": {"source": "test", "parameters": {"m2": -1.0}},
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,   # positive identity → negative eigenvalue
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-m2",
                            },
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                        ],
                    },
                }
            ],
            "coupling": {"mass_matrix_symbolic": [["-m2"]]},
        }
        spec_path = tmp_path / "tachyon.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")
        ret = main(["validate", str(spec_path), "--stability"])
        assert ret == 1
        err = capsys.readouterr().err
        assert "ERROR" in err

    def test_validate_stability_param_override_triggers_instability(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--param overrides should affect stability check (large coupling → tachyon)."""
        import json

        # coupled_scalars: stable when gCpl² < mPhi2 * mChi2 = 1*4 = 4
        # With gCpl=3.0 it becomes unstable (g² = 9 > 4)
        spec = {
            "metadata": {"source": "test", "parameters": {"mPhi2": 1.0, "mChi2": 4.0, "gCpl": 0.5}},
            "spacetime": {
                "dimension": 2,
                "signature": [-1, 1],
                "coordinates": ["t", "x"],
            },
            "fields": [
                {"name": "phi_0", "index": 0, "is_dynamical": True},
                {"name": "chi_0", "index": 1, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "phi_0", "coefficient_symbolic": "-mPhi2"},
                            {"coefficient": -0.5, "operator": "identity", "field": "chi_0", "coefficient_symbolic": "-gCpl"},
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -4.0, "operator": "identity", "field": "chi_0", "coefficient_symbolic": "-mChi2"},
                            {"coefficient": -0.5, "operator": "identity", "field": "phi_0", "coefficient_symbolic": "-gCpl"},
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "chi_0"},
                        ],
                    },
                },
            ],
            "coupling": {"mass_matrix_symbolic": [["-mPhi2", None], [None, "-mChi2"]]},
        }
        spec_path = tmp_path / "coupled.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")

        # Default params: gCpl=0.5 → stable
        ret_stable = main(["validate", str(spec_path), "--stability"])
        assert ret_stable == 0

        # Override: gCpl=3.0 → unstable (g²=9 > mPhi2*mChi2=4)
        ret_unstable = main(["validate", str(spec_path), "--stability", "--param", "gCpl=3.0"])
        assert ret_unstable == 1


class TestMeasureCommand:
    def test_measure_help(self) -> None:
        with pytest.raises(SystemExit, match="0"):
            main(["measure", "--help"])

    def test_measure_nonexistent_path(self) -> None:
        ret = main(
            ["measure", "/nonexistent/data_dir", "--spec", "/nonexistent/spec.json"]
        )
        assert ret == 1

    def test_measure_summary_text(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Default mode should print full summary with energy + conservation."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Energy Conservation:" in out
        assert "Per-Field Energy" in out

    def test_measure_json_output(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--json should produce valid JSON with simulation metadata."""
        import json as json_mod

        capsys.readouterr()  # flush fixture output
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        raw = capsys.readouterr().out
        data = json_mod.loads(raw)
        assert "simulation" in data
        assert "energy" in data or "conservation" in data

    def test_measure_energy_only(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=energy should compute only energy (no conversion/mixing)."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "energy",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Per-Field Energy" in out
        assert "Conversion" not in out

    def test_measure_conversion_requires_source(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
    ) -> None:
        """--what=conversion without --source should error."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "conversion",
            ]
        )
        assert ret == 1

    def test_measure_conversion_with_fields(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Explicit --source and --target should compute conversion."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "conversion",
                "--source",
                "phi_0",
                "--target",
                "chi_0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Conversion" in out
        assert "Peak P(t)" in out

    def test_measure_png_output(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        tmp_path: Path,
    ) -> None:
        """--output with .png extension should create a plot file."""
        output = tmp_path / "measurement.png"
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--output",
                str(output),
            ]
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_measure_spec_autodiscovery(
        self,
        coupled_scalars_dir: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Snapshot directory should contain metadata.json with spec_path for auto-discovery."""
        ret = main(["measure", str(coupled_scalars_dir)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "Energy Conservation:" in out

    def test_measure_quiet(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--quiet should suppress progress messages but keep results."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--quiet",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Loading:" not in out
        assert "Energy Conservation:" in out

    def test_measure_invalid_what(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
    ) -> None:
        """Unknown --what value should return 1."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "nonexistent_measure",
            ]
        )
        assert ret == 1

    def test_measure_param_override(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--param should override parameter values from snapshot data."""
        capsys.readouterr()  # flush fixture output
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "energy",
                "--param",
                "mPhi2=2.0",
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        import json as json_mod

        data = json_mod.loads(capsys.readouterr().out)
        assert data["simulation"]["parameters"]["mPhi2"] == 2.0

    def test_measure_mixing_standalone(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=mixing (without explicit conversion) should auto-detect source/target."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "mixing",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Mixing Length:" in out or "Mixing Length: not extracted" in out
        # Conversion should also have been computed (as a dependency)
        assert "Conversion" in out

    def test_measure_spectrum_individual(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=spectrum alone should compute power spectrum."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectrum",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Spectrum" in out
        # Should NOT include conversion or mixing sections
        assert "Conversion" not in out
        assert "Mixing" not in out

    def test_measure_combined_energy_conversion(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=energy,conversion should compute both but not mixing."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "energy,conversion",
                "--source",
                "phi_0",
                "--target",
                "chi_0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Per-Field Energy" in out
        assert "Conversion" in out
        assert "Mixing" not in out

    def test_measure_spectral_conversion_auto_target(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--source without --target should auto-select remaining dynamical fields."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectral_conversion",
                "--source",
                "phi_0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Spectral Conversion" in out

    def test_measure_spectral_conversion(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=spectral_conversion should compute P(k,t)."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectral_conversion",
                "--source",
                "phi_0",
                "--target",
                "chi_0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Spectral Conversion" in out
        assert "Active k-modes:" in out

    def test_measure_spectral_conversion_json(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=spectral_conversion --json should include spectral_conversion key."""
        import json as json_mod

        capsys.readouterr()  # flush
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectral_conversion",
                "--source",
                "phi_0",
                "--target",
                "chi_0",
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        data = json_mod.loads(capsys.readouterr().out)
        assert "spectral_conversion" in data
        assert data["spectral_conversion"]["n_active_modes"] > 0

    def test_measure_spectral_conversion_requires_source(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
    ) -> None:
        """--what=spectral_conversion without --source should error."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectral_conversion",
            ]
        )
        assert ret == 1

    def test_measure_spectral_conversion_plot(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        tmp_path: Path,
    ) -> None:
        """--what=spectral_conversion --output .png should create plot."""
        output = tmp_path / "sc_plot.png"
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "spectral_conversion",
                "--source",
                "phi_0",
                "--target",
                "chi_0",
                "--output",
                str(output),
            ]
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_measure_dispersion_text(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=dispersion should compute dispersion relation."""
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "dispersion",
                "--source",
                "phi_0",
            ]
        )
        assert ret == 0
        out = capsys.readouterr().out
        assert "Dispersion" in out

    def test_measure_dispersion_json(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=dispersion --json should include 'dispersion' key."""
        import json as json_mod

        capsys.readouterr()
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "dispersion",
                "--source",
                "phi_0",
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        data = json_mod.loads(capsys.readouterr().out)
        assert "dispersion" in data

    def test_measure_dispersion_plot(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        tmp_path: Path,
    ) -> None:
        """--what=dispersion --output .png should create plot."""
        output = tmp_path / "disp_plot.png"
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "dispersion",
                "--source",
                "phi_0",
                "--output",
                str(output),
            ]
        )
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_measure_dispersion_group_field_name(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--source phi_0,chi_0 should sum group power and join field names."""
        import json as json_mod

        capsys.readouterr()
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "dispersion",
                "--source",
                "phi_0,chi_0",
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        data = json_mod.loads(capsys.readouterr().out)
        assert "dispersion" in data
        assert data["dispersion"]["field"] == "phi_0, chi_0"

    def test_measure_dispersion_no_source_uses_all_dynamical(
        self,
        coupled_scalars_dir: Path,
        inline_coupled_scalars_json: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """--what=dispersion without --source uses all dynamical fields."""
        import json as json_mod

        capsys.readouterr()
        ret = main(
            [
                "measure",
                str(coupled_scalars_dir),
                "--spec",
                str(inline_coupled_scalars_json),
                "--what",
                "dispersion",
                "--json",
                "--quiet",
            ]
        )
        assert ret == 0
        data = json_mod.loads(capsys.readouterr().out)
        assert "dispersion" in data
        # Both phi_0 and chi_0 are dynamical — all should be included
        assert "phi_0" in data["dispersion"]["field"]
        assert "chi_0" in data["dispersion"]["field"]


class TestBackgroundFields:
    """Tests for [[background_fields]] TOML feature — WLS generation dry-runs."""

    def test_background_scalar_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Scalar background field should generate DefTensor + ComponentValue."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Scalar in External Potential"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - V[] * phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Background fields" in out
        assert "ComponentValue" in out
        assert "V0" in out
        assert "EulerLagrangeEquation" in out

    def test_background_vector_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Vector background field should generate ComponentValue for each component."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "EM in External B Field"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["B0val"]

[[background_fields]]
name = "B"
type = "vector"
components = [0, 0, "B0val"]

[lagrangian]
expression = "-1/2 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert out.count("ComponentValue") == 3
        assert "B0val" in out
        assert "{0," in out
        assert "{1," in out
        assert "{2," in out

    def test_background_tensor_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Rank-2 tensor background field should generate dim^2 ComponentValues."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Tensor Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "G"
type = "tensor"
rank = 2
symmetry = "symmetric"
components = [-1, 0, 0, 1]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert out.count("ComponentValue") == 4

    def test_background_field_in_lagrangian_substitution(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Background field name should be substituted in Lagrangian expression."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Substitution Test"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "V"
type = "scalar"
components = [1.0]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - V[] * phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Theory name "Substitution Test" -> prefix "st"
        assert "stV[]" in out

    def test_background_missing_components(self, tmp_path: Path) -> None:
        """Background field without components should fail validation."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "V"
type = "scalar"

[lagrangian]
expression = "phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_background_wrong_component_count(self, tmp_path: Path) -> None:
        """Background field with wrong number of components should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad Count"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "B"
type = "vector"
components = [0, 1]

[lagrangian]
expression = "phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_background_name_collision_dynamical(self, tmp_path: Path) -> None:
        """Background field with same name as dynamical field should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Collision"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "phi"
type = "scalar"
components = [1.0]

[lagrangian]
expression = "phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_background_name_collision_derived(self, tmp_path: Path) -> None:
        """Background field with same name as derived field should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Collision"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[[background_fields]]
name = "F"
type = "scalar"
components = [1.0]

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_no_background_fields_backward_compat(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Existing TOMLs without [[background_fields]] should still work."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "No Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Background fields" not in out
        assert "ComponentValue" not in out
        assert "EulerLagrangeEquation" in out

    def test_background_numeric_components(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Background field with numeric (int/float) components should work."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Numeric Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "G"
type = "vector"
components = [0, 1.5]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "1.5" in out

    def test_background_multi_field_coupling_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Two dynamical fields coupled through a background field."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Background Coupling"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "psi"
type = "scalar"

[[fields]]
name = "chi"
type = "scalar"

[constants]
names = ["gCpl", "Bval"]

[[background_fields]]
name = "B"
type = "scalar"
components = ["Bval"]

[lagrangian]
expression = "-1/2 CD[-a][psi[]] eta[a, b] CD[-b][psi[]] - 1/2 CD[-a][chi[]] eta[a, b] CD[-b][chi[]] + gCpl * B[] * psi[] * chi[]"

[parameters]
gCpl = 0.3
Bval = 1.0

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Background fields" in out
        assert "ComponentValue" in out
        assert "Bval" in out
        assert "DecomposeToComponents" in out
        assert "bcPsi[]" in out
        assert "bcChi[]" in out
        assert '"Bval" -> 1.0' in out

    def test_background_localized_unitstep_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Spatially-varying background with UnitStep should parse."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Localized Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0 * UnitStep[x[] - 30] * UnitStep[70 - x[]]"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - V[]/2 * phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "UnitStep" in out
        assert "ComponentValue" in out

    def test_scalar_background_explicit_substitution(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Scalar background should generate ReplaceAll before DecomposeToComponents."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Scalar Potential Well"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - V[]/2 * phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "spwV[] -> V0" in out
        decompose_pos = out.index("DecomposeToComponents")
        replace_pos = out.index("spwV[] -> V0")
        assert replace_pos < decompose_pos

    def test_scalar_background_multi_field_substitution(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Multi-field with scalar background should substitute in all EOMs."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Background Coupling"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "psi"
type = "scalar"

[[fields]]
name = "chi"
type = "scalar"

[constants]
names = ["Bval"]

[[background_fields]]
name = "B"
type = "scalar"
components = ["Bval"]

[lagrangian]
expression = "-1/2 CD[-a][psi[]] eta[a, b] CD[-b][psi[]] - 1/2 CD[-a][chi[]] eta[a, b] CD[-b][chi[]] + B[] * psi[] * chi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Both EOM variables should have the substitution
        assert out.count("bcB[] -> Bval") >= 2

    def test_scalar_background_localized_substitution(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Localized scalar background should substitute the full expression."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Localized Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0 * UnitStep[x[] - 30] * UnitStep[70 - x[]]"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - V[]/2 * phi[]^2"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "lbV[] -> V0 * UnitStep[x[] - 30] * UnitStep[70 - x[]]" in out

    def test_vector_background_explicit_substitution(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Vector background generates ReplaceAll after decomposition."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Vector No Sub"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[background_fields]]
name = "B"
type = "vector"
components = [0, 0, 1.0]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Vector backgrounds should NOT use scalar-style substitution
        assert "Substitute scalar background" not in out
        # But SHOULD use vector/tensor-style substitution after decomposition
        assert "Substitute vector/tensor background" in out

    def test_coupled_scattering_multifield_background(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Multi-field + scalar background with Exp should produce correct WLS."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Coupled Scalar Scattering"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[fields]]
name = "chi"
type = "scalar"

[constants]
names = ["mPhi2", "mChi2", "g0", "R"]

[[background_fields]]
name = "G"
type = "scalar"
components = ["g0 * Exp[-(x[]^2 + y[]^2) / (2 * R^2)]"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - 1/2 CD[-a][chi[]] eta[a, b] CD[-b][chi[]] - mPhi2/2 * phi[]^2 - mChi2/2 * chi[]^2 - G[] * phi[] * chi[]"

[parameters]
mPhi2 = 0.0
mChi2 = 1.0
g0 = 2.0
R = 8.0

[output]
path = "coupled_scattering.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Prefix from "Coupled Scalar Scattering" → "css"
        # Background G should be defined
        assert "cssG" in out
        # Scalar background should get ReplaceAll substitution
        assert "Substitute scalar background" in out
        # Both fields should appear in Euler-Lagrange
        assert "cssPhi" in out
        assert "cssChi" in out
        # Constants should be declared
        assert "DefConstantSymbol" in out
        assert "g0" in out
        assert "mChi2" in out


    def test_background_vector_curved_metric_contravariant(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Vector background in curved diagonal metric: contravariant A^μ = A_μ / g_{μμ}.

        For spherical metric diag[-1, 1, r², r²sin²θ], the θ-component (index 2)
        has g_{22} = x[]^2 (= r²).  Given covariant A_θ = B/(2*x[]^2), the correct
        contravariant value is A^θ = B/(2*x[]^4), NOT B/(2*x[]^2).

        Bug 2: prior to this fix, the same covariant value was used for the
        contravariant substitution in _wls_vector_background_substitution, giving
        wrong coupling coefficients in curved-spacetime theories.
        """
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Curved Metric Vector Background"

[spacetime]
dimension = 4
metric = "diagonal"
diagonal = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["Bpeak", "z0"]

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "Bpeak*z0^3/(2*x[]^2)", "0"]

[lagrangian]
expression = "-1/4 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Covariant substitution (index 2, -Cart): covariant value unchanged
        assert "Bpeak*z0^3/(2*x[]^2)" in out
        # Contravariant substitution (index 2, +Cart): must be A_θ / g_{θθ}
        # = Bpeak*z0^3/(2*x[]^2) / x[]^2 = Bpeak*z0^3/(2*x[]^4), computed via Simplify
        assert "Simplify[(Bpeak*z0^3/(2*x[]^2)) / (x[]^2)]" in out

    def test_background_vector_minkowski_covariant_eq_contravariant(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Vector background in Minkowski metric: A^μ = A_μ (spatial components equal).

        For Minkowski metric, g_{ii} = 1 for spatial components so A^i = A_i.
        The covariant == contravariant shortcut must be preserved (no Simplify wrapping).
        """
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Minkowski Vector Background"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["B0"]

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "B0", "0"]

[lagrangian]
expression = "-1/4 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # For Minkowski (empty metric_diagonal), contravariant = covariant
        assert "[{2, -" in out and "-> B0" in out
        assert "[{2, " in out and "-> B0" in out
        # Must NOT generate Simplify[] wrapping for Minkowski (old shortcut preserved)
        assert "Simplify[(B0)" not in out


class TestExceptionHandling:
    def test_value_error_shows_clean_message(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """ValueError should produce a clean 'Error:' message, not a traceback."""
        ret = main(["simulate", "/nonexistent/spec.json"])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Error:" in err

    def test_derive_derived_field_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Derived field F should generate DefTensor, MakeRule, and substitution."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Electromagnetism"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Derived field tensor definition
        assert "DefTensor" in out
        assert "Antisymmetric" in out
        # MakeRule for substitution
        assert "MakeRule" in out
        assert "MetricOn" in out
        # Substitution applied to Lagrangian
        assert "/." in out
        assert "ToCanonical" in out
        assert "ContractMetric" in out
        # VarD should vary w.r.t. A, not F
        assert "EulerLagrangeEquation" in out

    def test_derive_derived_scalar_dry_run(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Derived scalar field should generate DefTensor and MakeRule."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Test Derived Scalar"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[[derived_fields]]
name = "V"
type = "scalar"
definition = "phi[]^2"

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]] - V[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "MakeRule" in out
        assert "Rules" in out

    def test_derive_derived_field_name_collision(self, tmp_path: Path) -> None:
        """Derived field with same name as fundamental field should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad Collision"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[derived_fields]]
name = "A"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]]"

[lagrangian]
expression = "A[-a] eta[a, b] A[-b]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_derive_derived_field_missing_definition(self, tmp_path: Path) -> None:
        """Derived field without definition should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Missing Defn"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"

[lagrangian]
expression = "F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_derive_no_derived_fields_backward_compat(
        self,
        tmp_path: Path,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        """Existing TOMLs without [[derived_fields]] should still work."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Plain Scalar"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a, b] CD[-b][phi[]]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

        out = capsys.readouterr().out
        # Should NOT contain derived field artifacts
        assert "MakeRule" not in out
        assert "Derived field" not in out
        # Normal flow should work
        assert "EulerLagrangeEquation" in out

    def test_derive_rank_exceeds_max(self, tmp_path: Path) -> None:
        """Tensor rank exceeding max should be caught at config validation."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad Rank"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "T"
type = "tensor"
rank = 99

[lagrangian]
expression = "T[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    # -- [linearization] validation tests --

    def test_linearization_missing_expression(self, tmp_path: Path) -> None:
        """[linearization] without expression should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[linearization]
perturbation_field = "h"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_linearization_missing_perturbation_field(self, tmp_path: Path) -> None:
        """[linearization] without perturbation_field should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[linearization]
expression = "Einstein[CD][-a, -b]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_linearization_unknown_perturbation_field(self, tmp_path: Path) -> None:
        """perturbation_field not matching any [[fields]] entry should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[linearization]
expression = "Einstein[CD][-a, -b]"
perturbation_field = "g"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_linearization_and_lagrangian_coexist(self, tmp_path: Path) -> None:
        """[lagrangian] + [linearization] is valid (Lagrangian-first)."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Lagrangian-first linearization"

[spacetime]
dimension = 3
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[lagrangian]
expression = "RicciScalarCD[]"

[linearization]
perturbation_field = "h"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

    def test_linearization_or_lagrangian_required(self, tmp_path: Path) -> None:
        """Config with neither [lagrangian] nor [linearization] should fail."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1


class TestMatterPerturbations:
    """Tests for [[linearization.matter_perturbations]] multi-field support."""

    def test_matter_perturbation_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Einstein-Maxwell config with matter perturbation generates correct WLS."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Einstein-Maxwell"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[fields]]
name = "A"
type = "vector"

[[fields]]
name = "a"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "-B0 * z[]", "0"]

[constants]
names = ["kappa", "B0"]

[lagrangian]
expression = "(1/kappa^2) RicciScalarCD[] - 1/4 F[-a,-b] eta[a,c] eta[b,d] F[-c,-d]"

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[[gauge]]
field = "h"
type = "tt"

[[gauge]]
field = "a"
type = "lorenz"

[reduction]
type = "plane_wave"
propagation_axis = "z"

[parameters]
kappa = 1.0
B0 = 0.3

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0
        wls = capsys.readouterr().out
        # Should contain DefTensorPerturbation for matter field
        assert "DefTensorPerturbation" in wls
        # Should contain metric perturbation setup
        assert "SetupMetricPerturbation" in wls
        # Should have VarD for at least two fields
        assert wls.count("VarD[") >= 2
        # Should have DecomposeToComponents
        assert "DecomposeToComponents" in wls
        # Should have MakeRule for derived field F
        assert "MakeRule" in wls
        # Background ComponentValue for A
        assert "ComponentValue" in wls
        assert "-B0" in wls
        # Regression: field name "a" must not corrupt "eta[" substitution
        assert "Eta[" in wls, "eta[ should be substituted to {prefix}Eta["
        # The prefixed Eta should not be further corrupted by field "a" replacement
        assert "EttidalA[" not in wls, "field name 'a' corrupted Eta substitution"
        assert "ettidalA" not in wls, "field name 'a' corrupted eta substitution"
        # Background field DownValues + derivative evaluation in pipeline
        assert "SetBackgroundFieldDownValues" in wls
        assert '"BackgroundFieldRules"' in wls

    def test_matter_perturbation_vector_2d_dry_run(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """2D vector matter perturbation with position-dependent background."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Maxwell Perturbation 2D"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[[fields]]
name = "a"
type = "vector"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "-B0 * x[]"]

[constants]
names = ["m2", "B0"]

[lagrangian]
expression = "-1/4 (CD[-a][A[-b]] - CD[-b][A[-a]]) eta[a,c] eta[b,d] (CD[-c][A[-d]] - CD[-d][A[-c]]) - m2/2 A[-a] eta[a,b] A[-b]"

[linearization]

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[[gauge]]
field = "a"
type = "lorenz"

[parameters]
m2 = 1.0
B0 = 0.5

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0
        wls = capsys.readouterr().out
        # Matter-only perturbation (metric pert registered but h terms dropped)
        assert "DefTensorPerturbation" in wls
        assert "SetupMetricPerturbation" in wls
        assert "Matter-only: drop all metric perturbation orders" in wls
        # Background field evaluation pipeline
        assert "SetBackgroundFieldDownValues" in wls
        assert '"BackgroundFieldRules"' in wls
        assert "-B0 * x[]" in wls
        # Lorenz gauge on perturbation
        assert "BuildLorenzGaugeTerm" in wls

    def test_matter_perturbation_missing_field(self, tmp_path: Path) -> None:
        """Matter perturbation referencing non-existent field fails validation."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "0", "0"]

[lagrangian]
expression = "RicciScalarCD[]"

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "NoSuchField"
perturbation_name = "a"
background = "Abar"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_matter_perturbation_missing_background(self, tmp_path: Path) -> None:
        """Matter perturbation referencing non-existent background fails."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[fields]]
name = "A"
type = "vector"

[lagrangian]
expression = "RicciScalarCD[]"

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "NoSuchBg"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_matter_perturbation_without_lagrangian(self, tmp_path: Path) -> None:
        """Matter perturbations require a Lagrangian section."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[fields]]
name = "A"
type = "vector"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "0", "0"]

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_matter_perturbation_duplicate_field(self, tmp_path: Path) -> None:
        """Same field listed twice in matter_perturbations fails."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Bad"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[fields]]
name = "A"
type = "vector"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "0", "0"]

[lagrangian]
expression = "RicciScalarCD[]"

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "b"
background = "Abar"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 1

    def test_matter_perturbation_gauge_on_perturbation_name(
        self, tmp_path: Path
    ) -> None:
        """Gauge can target a perturbation_name (not just a [[fields]] name)."""
        config = tmp_path / "theory.toml"
        config.write_text("""
[theory]
name = "Gauge on pert"

[spacetime]
dimension = 4
metric = "minkowski"

[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"

[[fields]]
name = "A"
type = "vector"

[[fields]]
name = "a"
type = "vector"

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[[background_fields]]
name = "Abar"
type = "vector"
components = ["0", "0", "0", "0"]

[lagrangian]
expression = "(1/kappa^2) RicciScalarCD[] - 1/4 F[-a,-b] eta[a,c] eta[b,d] F[-c,-d]"

[constants]
names = ["kappa"]

[linearization]
perturbation_field = "h"

[[linearization.matter_perturbations]]
field = "A"
perturbation_name = "a"
background = "Abar"

[[gauge]]
field = "a"
type = "lorenz"

[reduction]
type = "plane_wave"
propagation_axis = "z"

[parameters]
kappa = 1.0

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0


class TestCoordValuesPreEvaluation:
    """Tests for _apply_coord_values and early coordinate-value substitution.

    This covers the optimization that eliminates killed-coordinate factors
    (e.g. Sin[y[]]^2) from Python-generated Wolfram strings *before* emission,
    so they never enter xPert / Christoffel / TT-traceless computations.
    """

    def test_apply_coord_values_empty(self) -> None:
        """Empty coord_values returns input unchanged."""
        from tidal.cli._derive import _apply_coord_values

        exprs = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]
        assert _apply_coord_values(exprs, {}) == exprs

    def test_apply_coord_values_substitutes(self) -> None:
        """y[] -> Pi/2 is substituted in all strings."""
        from tidal.cli._derive import _apply_coord_values

        exprs = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]
        result = _apply_coord_values(exprs, {"y": "Pi/2"})
        assert result == ["-1", "1", "x[]^2", "x[]^2*Sin[Pi/2]^2"]

    def test_apply_coord_values_multiple_coords(self) -> None:
        """Multiple coord substitutions are all applied."""
        from tidal.cli._derive import _apply_coord_values

        exprs = ["y[]*z[]", "Sin[y[]]*Cos[z[]]"]
        result = _apply_coord_values(exprs, {"y": "Pi/2", "z": "0"})
        assert result == ["Pi/2*0", "Sin[Pi/2]*Cos[0]"]

    def test_apply_coord_values_no_match(self) -> None:
        """Strings without the coord function are returned unchanged."""
        from tidal.cli._derive import _apply_coord_values

        exprs = ["x[]^2", "Bpeak*z0^3/(2*x[]^2)"]
        result = _apply_coord_values(exprs, {"y": "Pi/2"})
        assert result == exprs

    def test_metric_diagonal_not_pre_evaluated_in_wls(self, tmp_path: Path) -> None:
        """Metric matrix must NOT be pre-evaluated with coordinate_values.

        Substituting y→Pi/2 into the metric before xPert's L^(2) expansion
        creates a non-flat background: diag(-1,1,r²,r²) has non-zero Riemann,
        while the true spherical metric diag(-1,1,r²,r²sin²θ) is flat.
        The coordinate evaluation is deferred to _wls_plane_wave_coordinate_evaluation
        which applies "fieldEquations /. {y[] -> Pi/2}" after derivation.
        """
        from tidal.cli._derive import generate_wls
        import tomllib

        config_text = """
[theory]
name = "Spherical test"
[spacetime]
dimension = 4
metric = "diagonal"
diagonal = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]
[[fields]]
name = "phi"
type = "scalar"
[lagrangian]
expression = "-1/2 CD[-a][phi[]] CD[a][phi[]]"
[reduction]
type = "plane_wave"
propagation_axis = "x"
coordinate_values = {y = "Pi/2"}
[output]
path = "out.json"
"""
        config_path = tmp_path / "theory.toml"
        config_path.write_text(config_text)
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
        wls_text = generate_wls(config, str(config_path))
        # Metric matrix must NOT be pre-evaluated (creates non-flat background)
        assert "MetricMatrix /. {y[] -> Pi/2}" not in wls_text
        # But coordinate_values should still be applied to field equations after derivation
        assert "fieldEquations = fieldEquations /. {y[] -> Pi/2}" in wls_text

    def test_tt_traceless_weights_no_sin_y(self, tmp_path: Path) -> None:
        """TT-traceless substitution weights use Sin[Pi/2], not Sin[y[]]."""
        from tidal.cli._derive import generate_wls
        import tomllib

        config_text = """
[theory]
name = "GW spherical"
[spacetime]
dimension = 4
metric = "diagonal"
diagonal = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]
[[fields]]
name = "h"
type = "tensor"
rank = 2
symmetry = "symmetric"
[lagrangian]
expression = "(1/kappa^2) RicciScalarCD[]"
[constants]
names = ["kappa"]
[linearization]
perturbation_field = "h"
[[gauge]]
field = "h"
type = "tt"
[reduction]
type = "plane_wave"
propagation_axis = "x"
coordinate_values = {y = "Pi/2"}
[output]
path = "out.json"
"""
        config_path = tmp_path / "theory.toml"
        config_path.write_text(config_text)
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
        wls_text = generate_wls(config, str(config_path))
        wls_lines = wls_text.splitlines()
        # TT traceless substitution lines contain "args___" (RuleDelayed pattern)
        # The last-diagonal field name is prefix-dependent (e.g. gsH9, gedH9)
        tt_lines = [l for l in wls_lines if "args___" in l and "Sin[" in l]
        assert tt_lines, "Expected TT traceless substitution lines with Sin[] weight"
        for line in tt_lines:
            assert "Sin[y[]]" not in line, f"Sin[y[]] in TT traceless: {line}"
            assert "Sin[Pi/2]" in line, f"Expected Sin[Pi/2] in TT traceless: {line}"

    def test_coordinate_values_plane_wave_uses_expand(self, tmp_path: Path) -> None:
        """Plane-wave coordinate eval step uses Expand not Simplify."""
        from tidal.cli._derive import generate_wls
        import tomllib

        config_text = """
[theory]
name = "Coord eval speed test"
[spacetime]
dimension = 4
metric = "diagonal"
diagonal = ["-1", "1", "x[]^2", "x[]^2*Sin[y[]]^2"]
[[fields]]
name = "phi"
type = "scalar"
[lagrangian]
expression = "-1/2 CD[-a][phi[]] CD[a][phi[]]"
[reduction]
type = "plane_wave"
propagation_axis = "x"
coordinate_values = {y = "Pi/2"}
[output]
path = "out.json"
"""
        config_path = tmp_path / "theory.toml"
        config_path.write_text(config_text)
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
        wls_text = generate_wls(config, str(config_path))
        assert "fieldEquations = fieldEquations /. {y[] -> Pi/2}" in wls_text
        assert "Expand[fieldEquations[[k, 2]]]" in wls_text

    def test_no_coordinate_values_noop(self, tmp_path: Path) -> None:
        """Absent coordinate_values emits no early metric substitution."""
        from tidal.cli._derive import generate_wls
        import tomllib

        config_text = """
[theory]
name = "No coord values"
[spacetime]
dimension = 4
metric = "minkowski"
[[fields]]
name = "phi"
type = "scalar"
[lagrangian]
expression = "-1/2 CD[-a][phi[]] CD[a][phi[]]"
[reduction]
type = "plane_wave"
propagation_axis = "z"
[output]
path = "out.json"
"""
        config_path = tmp_path / "theory.toml"
        config_path.write_text(config_text)
        with open(config_path, "rb") as f:
            config = tomllib.load(f)
        wls_text = generate_wls(config, str(config_path))
        # No coordinate_values → no early metric ReplaceAll
        assert "MetricMatrix /. {" not in wls_text
        # No coordinate_values → no killed-coordinate field-equation substitution
        # (the plane-wave reduction does emit "fieldEquations /. {" for derivative-zeroing,
        # so we check specifically for a coordinate-value substitution like "y[] ->")
        assert "-> Pi/2" not in wls_text
