"""Tests for the ``tg`` command-line interface."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import ClassVar

import pytest

from torsion_gertsenshtein.cli import main

EXAMPLES_DIR = Path(__file__).parent.parent / "examples" / "data"


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
        """``python -m torsion_gertsenshtein.cli`` should be importable."""
        import torsion_gertsenshtein.cli.__main__  # noqa: F401

    def test_get_version_fallback(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """_get_version() should return 'unknown' when package is not installed."""
        from importlib.metadata import PackageNotFoundError

        from torsion_gertsenshtein.cli import _get_version

        def _raise(_name: str) -> str:
            raise PackageNotFoundError

        monkeypatch.setattr("importlib.metadata.version", _raise)
        assert _get_version() == "unknown"

    def test_entry_point_subprocess(self) -> None:
        """The ``tg`` entry point should be callable as a subprocess."""
        result = subprocess.run(
            [sys.executable, "-m", "torsion_gertsenshtein.cli", "--version"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0
        assert "tg" in result.stdout


class TestInspectCommand:
    @pytest.mark.parametrize(
        ("fixture_name", "dim_label", "fields", "expected_strs"),
        [
            ("klein_gordon_1d_json", "1+1D", 1, ["phi_0"]),
            ("massive_3form_json", "3+1D", 4, ["C_0", "C_3", "m2"]),
            ("coupled_scalars_json", "1+1D", 2, ["phi_0", "chi_0"]),
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
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(["inspect", str(klein_gordon_1d_json), "--params"])
        assert ret == 0

    def test_inspect_json_output(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--json should output valid, parseable JSON."""
        import json

        ret = main(["inspect", str(klein_gordon_1d_json), "--json"])
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
        self, coupled_scalars_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--json --params should include default parameter values."""
        import json

        ret = main(["inspect", str(coupled_scalars_json), "--json", "--params"])
        assert ret == 0

        data = json.loads(capsys.readouterr().out)
        assert "required" in data["parameters"]
        assert len(data["fields"]) == 2


class TestListCommand:
    def test_list_default_dir(self, capsys: pytest.CaptureFixture[str]) -> None:
        ret = main(["list"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "klein_gordon_1d.json" in out
        assert "specifications found" in out

    def test_list_custom_dir(self, capsys: pytest.CaptureFixture[str]) -> None:
        ret = main(["list", "--dir", str(EXAMPLES_DIR)])
        assert ret == 0

        out = capsys.readouterr().out
        assert "specifications found" in out

    def test_list_nonexistent_dir(self) -> None:
        ret = main(["list", "--dir", "/nonexistent/dir"])
        assert ret == 1


class TestSimulateCommand:
    def test_simulate_1d_summary(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main(["simulate", str(klein_gordon_1d_json), "--param", "m2=1.0", "--t-end", "1.0", "--no-plot"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "phi_0" in out

    def test_simulate_with_params(
        self, klein_gordon_3d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_3d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--grid-shape", "8",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_png_output(
        self, klein_gordon_1d_json: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "test_output.png"
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()

    def test_simulate_2d_plot_output(
        self, chern_simons_json: Path, tmp_path: Path,
    ) -> None:
        """2D spec should produce a non-empty PNG file (exercises plot_2d path)."""
        output = tmp_path / "cs_2d.png"
        ret = main([
            "simulate", str(chern_simons_json),
            "--grid-shape", "8",
            "--t-end", "0.2",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()
        assert output.stat().st_size > 0

    def test_simulate_npz_output(
        self, klein_gordon_1d_json: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "test_output.npz"
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()

    # --- IC types (parametrized) ---

    @pytest.mark.parametrize("ic_type", ["gaussian", "plane-wave", "zero"])
    def test_simulate_ic_types(
        self,
        klein_gordon_1d_json: Path,
        ic_type: str,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", ic_type,
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_gaussian_custom(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", "gaussian",
            "--ic-width", "1.5",
            "--ic-amplitude", "2.0",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_off_center_gaussian(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", "gaussian",
            "--ic-center", "30.0",
            "--ic-width", "3.0",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_chern_simons_constraint(
        self, chern_simons_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Test simulation of a system with constraint (time_order=0) + dynamical fields."""
        ret = main([
            "simulate", str(chern_simons_json),
            "--t-end", "0.5",
            "--grid-shape", "8",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "A_0" in out

    def test_simulate_invalid_param_format(self, klein_gordon_1d_json: Path) -> None:
        ret = main(["simulate", str(klein_gordon_1d_json), "--param", "bad_no_equals", "--no-plot"])
        assert ret == 1

    def test_simulate_invalid_ic_component(self, klein_gordon_1d_json: Path) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--ic-component", "nonexistent_field",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_nonexistent_file(self) -> None:
        ret = main(["simulate", "/nonexistent/file.json", "--no-plot"])
        assert ret == 1

    def test_simulate_custom_grid(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--grid-shape", "32",
            "--bounds", "0:20",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    # --- Feature: --bc (mixed boundary conditions) ---

    def test_simulate_bc_single(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--bc", "neumann",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_bc_mixed_2d(
        self, polar_kg_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(polar_kg_json),
            "--bc", "neumann,periodic",
            "--grid-shape", "16",
            "--bounds", "0.5:8,0:6.28",
            "--ic", "gaussian",
            "--ic-center", "3.0,3.14",
            "--ic-width", "0.5",
            "--t-end", "0.5",
            "--dt", "0.01",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_bc_invalid_type(self, klein_gordon_1d_json: Path) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--bc", "invalid_bc",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_bc_wrong_count(self, klein_gordon_1d_json: Path) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--bc", "neumann,periodic",  # 2 values for 1D
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_bc_dirichlet_rejected(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Dirichlet BC is not supported by py-pde; should fail with clear error."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--bc", "dirichlet",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Dirichlet" in err
        assert "not supported" in err

    # --- Feature: --ic formula ---

    def test_simulate_formula_ic(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", "formula",
            "--ic-formula", "np.exp(-((x - 5)**2) / 2)",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_formula_ic_2d(
        self, polar_kg_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(polar_kg_json),
            "--ic", "formula",
            "--ic-formula", "np.exp(-((x - 3)**2 + (y - pi)**2) / 0.5**2)",
            "--grid-shape", "16",
            "--bounds", "0.5:8,0:6.28",
            "--bc", "neumann,periodic",
            "--t-end", "0.5",
            "--dt", "0.01",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_formula_ic_missing_expr(self, klein_gordon_1d_json: Path) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--ic", "formula",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_formula_ic_constant(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", "formula",
            "--ic-formula", "0.5",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_formula_ic_undefined_var(self, klein_gordon_1d_json: Path) -> None:
        """Formula with undefined variable should fail."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--ic", "formula",
            "--ic-formula", "badvar * 2",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_formula_ic_attribute_access_rejected(
        self, klein_gordon_1d_json: Path,
    ) -> None:
        """Formula with attribute access should be rejected by AST validator."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--ic", "formula",
            "--ic-formula", "x.__class__.__name__",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    # --- Feature: --mode constraint ---

    def test_simulate_constraint_mode(
        self, electrostatics_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(electrostatics_json),
            "--mode", "constraint",
            "--grid-shape", "16",
            "--bc", "neumann",
            "--ic", "gaussian",
            "--ic-component", "rho",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Constraint solve complete" in out

    # --- Solver options ---

    def test_simulate_explicit_dt(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--dt", "0.01",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_scipy_scheme(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Verify scipy solver uses py-pde ScipySolver."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--scheme", "scipy",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_custom_snapshots(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--snapshots", "0.25",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    # --- 3D simulation ---

    def test_simulate_3d(
        self, klein_gordon_3d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_3d_json),
            "--param", "m2=1.0",
            "--grid-shape", "8",
            "--t-end", "0.2",
            "--no-plot",
        ])
        assert ret == 0

    # --- Edge-case validation ---

    def test_simulate_explicit_format_flag(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """Explicit --format flag should work for all formats."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--format", "summary",
        ])
        assert ret == 0
        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_periodic_flag(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--no-periodic should produce non-periodic (Neumann) BCs."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--no-periodic",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_wavevector_custom(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--ic-wavevector should override default wavevector for plane-wave."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic", "plane-wave",
            "--ic-wavevector", "0.5",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_ic_center_wrong_dim(self, klein_gordon_1d_json: Path) -> None:
        """--ic-center with wrong dimension count should fail."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--ic-center", "5.0,5.0",  # 2 values for 1D
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_snapshots_nonpositive(self, klein_gordon_1d_json: Path) -> None:
        """--snapshots with zero or negative value should fail."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--snapshots", "0",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_dt_nonpositive(self, klein_gordon_1d_json: Path) -> None:
        """--dt with zero or negative value should fail."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--dt", "-0.1",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_t_end_nonpositive(self, klein_gordon_1d_json: Path) -> None:
        """--t-end with zero or negative value should fail."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0",
            "--no-plot",
        ])
        assert ret == 1

    # --- Feature: --quiet ---

    def test_simulate_quiet_suppresses_progress(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """--quiet should suppress progress messages but keep results."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--no-plot",
            "--quiet",
        ])
        assert ret == 0
        out = capsys.readouterr().out
        # Results should still appear
        assert "Results:" in out
        # Progress messages should NOT appear
        assert "Loading equation specification" not in out
        assert "Building PDE" not in out
        assert "Running simulation" not in out

    def test_simulate_quiet_short_flag(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """-q should work as shorthand for --quiet."""
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--no-plot",
            "-q",
        ])
        assert ret == 0
        out = capsys.readouterr().out
        assert "Results:" in out
        assert "Loading equation specification" not in out


class TestDeriveCommand:
    def test_derive_toml_dry_run(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
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

    def test_derive_multi_field_dry_run(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
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

    def test_derive_tensor_field_dry_run(self, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
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
        from torsion_gertsenshtein.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG)
        assert "$InputFileName" not in script
        assert "/torsion_gertsenshtein/wolfram" in script

    def test_output_path_is_absolute(self) -> None:
        from torsion_gertsenshtein.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG, config_dir=Path("/project"))
        # The output path "output.json" should be resolved to /project/output.json
        assert '"/project/output.json"' in script

    def test_output_override_absolute(self) -> None:
        from torsion_gertsenshtein.cli._derive import generate_wls

        script = generate_wls(self._MINIMAL_CONFIG, output_override="/tmp/test.json")
        assert '"/tmp/test.json"' in script

    def test_output_path_resolves_relative_to_config_dir(self, tmp_path: Path) -> None:
        from torsion_gertsenshtein.cli._derive import generate_wls

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
    """Tests for ``tg validate``."""

    def test_validate_valid_spec(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        ret = main(["validate", str(klein_gordon_1d_json)])
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        """A spec with symbolic coefficient but no parameter defaults should warn."""
        import json

        spec = {
            "metadata": {"source": "test"},
            "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [{
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {"type": "linear_combination", "terms": [
                    {"coefficient": -1.0, "operator": "identity", "field": "phi_0",
                     "coefficient_symbolic": "-m2"},
                    {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                ]},
            }],
            "coupling": {"mass_matrix": [[1.0]], "coupling_matrix": [[0.0]]},
        }
        spec_path = tmp_path / "no_defaults.json"
        spec_path.write_text(json.dumps(spec), encoding="utf-8")
        main(["validate", str(spec_path)])
        err = capsys.readouterr().err
        assert "WARNING" in err
        assert "m2" in err

    def test_validate_coupled_scalars(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Coupled scalars spec should validate successfully."""
        spec_path = EXAMPLES_DIR / "coupled_scalars.json"
        if not spec_path.exists():
            pytest.skip("coupled_scalars.json not found")
        ret = main(["validate", str(spec_path)])
        assert ret == 0
        out = capsys.readouterr().out
        assert "OK" in out

    def test_validate_unknown_operator(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        """JSON with unknown operator should fail validation."""
        import json

        spec = {
            "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"operator": "nonexistent_op", "field": "phi_0", "coefficient": 1.0}
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
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str],
    ) -> None:
        """JSON with unknown field reference should fail validation."""
        import json

        spec = {
            "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
            "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"operator": "identity", "field": "chi_99", "coefficient": 1.0}
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
        self, tmp_path: Path,
    ) -> None:
        """Passing a directory path should fail validation."""
        ret = main(["validate", str(tmp_path)])
        assert ret == 1


class TestExceptionHandling:
    def test_value_error_shows_clean_message(self, capsys: pytest.CaptureFixture[str]) -> None:
        """ValueError should produce a clean 'Error:' message, not a traceback."""
        ret = main(["simulate", "/nonexistent/spec.json"])
        assert ret == 1
        err = capsys.readouterr().err
        assert "Error:" in err

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
