"""Tests for the ``tg`` command-line interface."""

from __future__ import annotations

from pathlib import Path

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
        ret = main(["simulate", str(klein_gordon_1d_json), "--t-end", "1.0", "--no-plot"])
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
            "--t-end", "0.5",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()

    def test_simulate_npz_output(
        self, klein_gordon_1d_json: Path, tmp_path: Path
    ) -> None:
        output = tmp_path / "test_output.npz"
        ret = main([
            "simulate", str(klein_gordon_1d_json),
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

    def test_simulate_bc_dirichlet(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
            "--bc", "dirichlet",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    # --- Feature: --ic formula ---

    def test_simulate_formula_ic(
        self, klein_gordon_1d_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(klein_gordon_1d_json),
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

    # --- Feature: --mode constraint ---

    def test_simulate_constraint_mode(
        self, electrostatics_json: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        ret = main([
            "simulate", str(electrostatics_json),
            "--mode", "constraint",
            "--grid-shape", "16",
            "--bc", "dirichlet",
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


class TestDeriveAbsolutePaths:
    """Verify generated WLS scripts use absolute paths (not $InputFileName-relative)."""

    _MINIMAL_CONFIG: dict[str, object] = {
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
