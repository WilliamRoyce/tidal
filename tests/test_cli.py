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


class TestInspectCommand:
    def test_inspect_klein_gordon_1d(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main(["inspect", str(json_path)])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Spacetime:" in out
        assert "1+1D" in out
        assert "phi_0" in out

    def test_inspect_massive_3form(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "massive_3form.json"
        if not json_path.exists():
            pytest.skip("massive_3form.json not found")

        ret = main(["inspect", str(json_path)])
        assert ret == 0

        out = capsys.readouterr().out
        assert "4 components" in out
        assert "C_0" in out
        assert "C_3" in out
        assert "m2" in out
        assert "3+1D" in out

    def test_inspect_coupled_scalars(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "coupled_scalars.json"
        if not json_path.exists():
            pytest.skip("coupled_scalars.json not found")

        ret = main(["inspect", str(json_path)])
        assert ret == 0

        out = capsys.readouterr().out
        assert "2 components" in out
        assert "phi_0" in out
        assert "chi_0" in out

    def test_inspect_nonexistent_file(self) -> None:
        ret = main(["inspect", "/nonexistent/file.json"])
        assert ret == 1

    def test_inspect_with_params_flag(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main(["inspect", str(json_path), "--params"])
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
    def test_simulate_1d_summary(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main(["simulate", str(json_path), "--t-end", "1.0", "--no-plot"])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "phi_0" in out

    def test_simulate_with_params(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_3d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_3d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--param", "m2=1.0",
            "--t-end", "0.5",
            "--grid-shape", "8",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_png_output(self, tmp_path: Path) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        output = tmp_path / "test_output.png"
        ret = main([
            "simulate", str(json_path),
            "--t-end", "0.5",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()

    def test_simulate_npz_output(self, tmp_path: Path) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        output = tmp_path / "test_output.npz"
        ret = main([
            "simulate", str(json_path),
            "--t-end", "0.5",
            "--output", str(output),
        ])
        assert ret == 0
        assert output.exists()

    def test_simulate_zero_ic(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--ic", "zero",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_plane_wave_ic(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--ic", "plane-wave",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

    def test_simulate_gaussian_ic(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--ic", "gaussian",
            "--ic-width", "1.5",
            "--ic-amplitude", "2.0",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out

    def test_simulate_chern_simons_constraint(self, capsys: pytest.CaptureFixture[str]) -> None:
        """Test simulation of a system with constraint (time_order=0) + dynamical fields."""
        json_path = EXAMPLES_DIR / "chern_simons_3d.json"
        if not json_path.exists():
            pytest.skip("chern_simons_3d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--t-end", "0.5",
            "--grid-shape", "8",
            "--no-plot",
        ])
        assert ret == 0

        out = capsys.readouterr().out
        assert "Results:" in out
        assert "A_0" in out

    def test_simulate_invalid_param_format(self) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main(["simulate", str(json_path), "--param", "bad_no_equals", "--no-plot"])
        assert ret == 1

    def test_simulate_invalid_ic_component(self) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--ic-component", "nonexistent_field",
            "--t-end", "0.5",
            "--no-plot",
        ])
        assert ret == 1

    def test_simulate_nonexistent_file(self) -> None:
        ret = main(["simulate", "/nonexistent/file.json", "--no-plot"])
        assert ret == 1

    def test_simulate_custom_grid(self, capsys: pytest.CaptureFixture[str]) -> None:
        json_path = EXAMPLES_DIR / "klein_gordon_1d.json"
        if not json_path.exists():
            pytest.skip("klein_gordon_1d.json not found")

        ret = main([
            "simulate", str(json_path),
            "--grid-shape", "32",
            "--bounds", "0:20",
            "--t-end", "0.5",
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
        ret = main(["derive", str(config), "--save-script", str(script_path)])
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
