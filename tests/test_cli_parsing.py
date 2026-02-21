"""Unit tests for CLI parsing helpers.

Tests the pure-function parsers in ``tidal.cli._simulate``
independently from the full simulation pipeline.
"""

from __future__ import annotations

from argparse import Namespace
from types import SimpleNamespace

import pytest

from tidal.cli._simulate import (
    _infer_output_format,
    _parse_bc,
    _parse_bounds,
    _parse_grid_shape,
    _parse_params,
    _parse_single_bound,
)

# ==================== _parse_grid_shape ====================


class TestParseGridShape:
    @pytest.mark.parametrize(
        ("raw", "dim", "expected"),
        [
            (None, 1, [64]),
            (None, 2, [32, 32]),
            (None, 3, [16, 16, 16]),
            ("32", 1, [32]),
            ("32", 2, [32, 32]),
            ("32", 3, [32, 32, 32]),
            ("32,64", 2, [32, 64]),
            ("8,16,32", 3, [8, 16, 32]),
        ],
    )
    def test_valid(self, raw: str | None, dim: int, expected: list[int]) -> None:
        assert _parse_grid_shape(raw, dim) == expected

    @pytest.mark.parametrize(
        ("raw", "dim"),
        [
            ("32,64", 1),     # 2 values for 1D
            ("32,64", 3),     # 2 values for 3D
            ("8,16,32", 2),   # 3 values for 2D
        ],
    )
    def test_dimension_mismatch(self, raw: str, dim: int) -> None:
        with pytest.raises(ValueError, match="--grid-shape"):
            _parse_grid_shape(raw, dim)


# ==================== _parse_bounds ====================


class TestParseBounds:
    @pytest.mark.parametrize(
        ("raw", "dim", "expected"),
        [
            (None, 1, [(0.0, 10.0)]),
            (None, 2, [(0.0, 10.0), (0.0, 10.0)]),
            ("0:20", 1, [(0.0, 20.0)]),
            ("0:20", 2, [(0.0, 20.0), (0.0, 20.0)]),
            ("0:20,0:10", 2, [(0.0, 20.0), (0.0, 10.0)]),
            ("-5:5,0:10,0:20", 3, [(-5.0, 5.0), (0.0, 10.0), (0.0, 20.0)]),
        ],
    )
    def test_valid(
        self,
        raw: str | None,
        dim: int,
        expected: list[tuple[float, float]],
    ) -> None:
        assert _parse_bounds(raw, dim) == expected

    @pytest.mark.parametrize(
        ("raw", "dim"),
        [
            ("0:20,0:10", 1),   # 2 values for 1D
            ("0:20,0:10", 3),   # 2 values for 3D
        ],
    )
    def test_dimension_mismatch(self, raw: str, dim: int) -> None:
        with pytest.raises(ValueError, match="--bounds"):
            _parse_bounds(raw, dim)


class TestParseSingleBound:
    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ("0:20", (0.0, 20.0)),
            ("-5:5", (-5.0, 5.0)),
            ("0.5:8", (0.5, 8.0)),
        ],
    )
    def test_valid(self, raw: str, expected: tuple[float, float]) -> None:
        assert _parse_single_bound(raw) == expected

    def test_missing_colon(self) -> None:
        with pytest.raises(ValueError, match="LO:HI"):
            _parse_single_bound("20")

    def test_non_numeric(self) -> None:
        with pytest.raises(ValueError, match="could not convert"):
            _parse_single_bound("abc:def")

    def test_inverted_bounds_rejected(self) -> None:
        """Inverted bounds (lo >= hi) should raise ValueError."""
        with pytest.raises(ValueError, match=r"lower.*must be less than upper"):
            _parse_single_bound("20:0")

    def test_equal_bounds_rejected(self) -> None:
        """Equal bounds (lo == hi) should raise ValueError."""
        with pytest.raises(ValueError, match=r"lower.*must be less than upper"):
            _parse_single_bound("5:5")


# ==================== _parse_bc ====================


class TestParseBC:
    @pytest.mark.parametrize(
        ("raw", "periodic", "dim", "expected"),
        [
            (None, True, 1, True),
            (None, False, 1, False),
            ("periodic", True, 1, [True]),
            ("periodic", True, 2, [True, True]),
            ("neumann", True, 1, [False]),
            ("neumann,periodic", False, 2, [False, True]),
            ("periodic,neumann,neumann", True, 3, [True, False, False]),
        ],
    )
    def test_valid(
        self,
        raw: str | None,
        periodic: bool,
        dim: int,
        expected: bool | list[bool],
    ) -> None:
        assert _parse_bc(raw, periodic=periodic, spatial_dim=dim) == expected

    def test_dirichlet_rejected(self) -> None:
        """Dirichlet BC is not supported by py-pde and should raise a clear error."""
        with pytest.raises(ValueError, match=r"Dirichlet.*not supported"):
            _parse_bc("dirichlet", periodic=True, spatial_dim=1)

    def test_invalid_bc_type(self) -> None:
        with pytest.raises(ValueError, match="Invalid boundary condition"):
            _parse_bc("robin", periodic=True, spatial_dim=1)

    def test_wrong_count(self) -> None:
        with pytest.raises(ValueError, match="--bc"):
            _parse_bc("neumann,periodic", periodic=True, spatial_dim=1)

    def test_bc_overrides_periodic_flag(self) -> None:
        """--bc should override --periodic flag."""
        result = _parse_bc("neumann", periodic=True, spatial_dim=1)
        assert result == [False]


# ==================== _parse_params ====================


def _make_spec_stub(metadata: dict[str, object] | None = None) -> object:
    """Create a minimal spec-like object with a metadata dict."""
    return SimpleNamespace(metadata=metadata or {})


class TestParseParams:
    def test_empty_list_no_metadata(self) -> None:
        spec = _make_spec_stub()
        assert _parse_params([], spec) == {}  # type: ignore[arg-type]

    def test_metadata_defaults_loaded(self) -> None:
        spec = _make_spec_stub({"parameters": {"m2": 1.0, "g": 0.5}})
        result = _parse_params([], spec)  # type: ignore[arg-type]
        assert result == {"m2": 1.0, "g": 0.5}

    def test_cli_overrides_metadata(self) -> None:
        spec = _make_spec_stub({"parameters": {"m2": 1.0}})
        result = _parse_params(["m2=2.0"], spec)  # type: ignore[arg-type]
        assert result == {"m2": 2.0}

    def test_missing_equals_raises(self) -> None:
        spec = _make_spec_stub()
        with pytest.raises(ValueError, match="KEY=VALUE"):
            _parse_params(["bad_no_equals"], spec)  # type: ignore[arg-type]

    def test_non_numeric_value_raises(self) -> None:
        spec = _make_spec_stub()
        with pytest.raises(ValueError, match="Must be a number"):
            _parse_params(["m2=abc"], spec)  # type: ignore[arg-type]

    def test_non_numeric_metadata_skipped(self, capsys: pytest.CaptureFixture[str]) -> None:
        spec = _make_spec_stub({"parameters": {"note": "not a number", "m2": 1.0}})
        result = _parse_params([], spec)  # type: ignore[arg-type]
        assert result == {"m2": 1.0}
        assert "Warning" in capsys.readouterr().err

    def test_no_parameters_key_in_metadata(self) -> None:
        spec = _make_spec_stub({"source": "xAct"})
        assert _parse_params([], spec) == {}  # type: ignore[arg-type]

    def test_warns_unknown_param(self, capsys: pytest.CaptureFixture[str]) -> None:
        spec = _make_spec_stub({"parameters": {"m2": 1.0}})
        result = _parse_params(["m2=2.0", "bogus=3.0"], spec)  # type: ignore[arg-type]
        assert result == {"m2": 2.0, "bogus": 3.0}
        err = capsys.readouterr().err
        assert "bogus" in err
        assert "not found" in err

    def test_no_warning_for_known_param(self, capsys: pytest.CaptureFixture[str]) -> None:
        spec = _make_spec_stub({"parameters": {"m2": 1.0}})
        _parse_params(["m2=2.0"], spec)  # type: ignore[arg-type]
        err = capsys.readouterr().err
        assert "not found" not in err


# ==================== _infer_output_format ====================


def _make_args(**kwargs: object) -> Namespace:
    """Create a Namespace with defaults for _infer_output_format."""
    defaults: dict[str, object] = {"no_plot": False, "output_format": None, "output": None}
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestInferOutputFormat:
    def test_no_plot_returns_summary(self) -> None:
        assert _infer_output_format(_make_args(no_plot=True)) == "summary"

    def test_explicit_format_wins(self) -> None:
        assert _infer_output_format(_make_args(output_format="png")) == "png"

    def test_npz_extension_raises(self) -> None:
        """NPZ format is no longer supported — must raise ValueError."""
        with pytest.raises(ValueError, match="no longer supported"):
            _infer_output_format(_make_args(output="foo.npz"))

    def test_png_extension(self) -> None:
        assert _infer_output_format(_make_args(output="foo.png")) == "png"

    def test_svg_extension(self) -> None:
        assert _infer_output_format(_make_args(output="foo.svg")) == "svg"

    def test_pdf_extension(self) -> None:
        assert _infer_output_format(_make_args(output="foo.pdf")) == "pdf"

    def test_jpg_extension(self) -> None:
        assert _infer_output_format(_make_args(output="foo.jpg")) == "jpg"

    def test_no_output_defaults_png(self) -> None:
        assert _infer_output_format(_make_args()) == "png"

    def test_no_plot_takes_priority(self) -> None:
        """--no-plot should win even if --format or --output is given."""
        assert _infer_output_format(_make_args(no_plot=True, output_format="png")) == "summary"


class TestValidateFormulaAst:
    """Tests for _validate_formula_ast AST-based formula validation."""

    def test_valid_simple_expression(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast("x + 1", {"x"})

    def test_valid_function_calls(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast("sin(x) + cos(y)", {"sin", "cos", "x", "y"})

    def test_valid_numpy_attribute(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast("np.exp(-x**2)", {"np", "x"})

    def test_rejects_unknown_name(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(ValueError, match="Disallowed name 'badvar'"):
            _validate_formula_ast("badvar * 2", {"x"})

    def test_rejects_attribute_access(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(ValueError, match="Attribute access not allowed"):
            _validate_formula_ast("x.__class__", {"x"})

    def test_rejects_import_name(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(ValueError, match="Disallowed name '__import__'"):
            _validate_formula_ast("__import__('os')", {"x"})

    def test_rejects_lambda(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(TypeError, match=r"Disallowed construct.*Lambda"):
            _validate_formula_ast("(lambda: 1)()", {"x"})

    def test_valid_ternary(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast("x if x > 0 else -x", {"x"})

    def test_valid_complex_math(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast(
            "exp(-((x - 5)**2 + (y - pi)**2) / 0.5**2)",
            {"exp", "x", "y", "pi"},
        )

    def test_rejects_list_comprehension(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(TypeError, match=r"Disallowed construct.*ListComp"):
            _validate_formula_ast("[i for i in x]", {"x", "i"})

    def test_rejects_walrus_operator(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(TypeError, match=r"Disallowed construct.*NamedExpr"):
            _validate_formula_ast("(y := x + 1)", {"x", "y"})

    def test_rejects_nested_attribute(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        with pytest.raises(ValueError, match="Attribute access not allowed"):
            _validate_formula_ast("x.a.b", {"x"})

    def test_allows_subscript_slicing(self) -> None:
        from tidal.cli._simulate import _validate_formula_ast

        _validate_formula_ast("x[0:5]", {"x"})


# ==================== generate_wls constraint_solver ====================


class TestConstraintSolverToml:
    """Tests for [constraint_solver] TOML → WLS metadata generation."""

    def _generate(self, config: dict[str, object]) -> str:
        """Generate WLS from minimal config with constraint_solver."""
        from tidal.cli._derive import generate_wls

        base: dict[str, object] = {
            "theory": {"name": "Test"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "lagrangian": {
                "expression": "CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"
            },
            "output": {"path": "out.json"},
        }
        base.update(config)
        return generate_wls(base)

    def test_constraint_solver_toml_in_wls(self) -> None:
        """[constraint_solver] in TOML → solve_constraints -> True in WLS."""
        wls = self._generate({"constraint_solver": {"enabled": True}})
        assert '"solve_constraints" -> True' in wls

    def test_constraint_solver_absent_no_flag(self) -> None:
        """No [constraint_solver] → no solve_constraints in WLS."""
        wls = self._generate({})
        assert "solve_constraints" not in wls

    def test_constraint_solver_dirichlet_bcs(self) -> None:
        """Dirichlet BC config → correct Wolfram Association syntax."""
        wls = self._generate({
            "constraint_solver": {
                "enabled": True,
                "boundary_conditions": {
                    "x": {"type": "dirichlet", "value": 0.0},
                    "y": {"type": "dirichlet", "value": 0.0},
                },
            }
        })
        assert '"type" -> "dirichlet"' in wls
        assert '"value" -> 0.0' in wls

    def test_constraint_solver_periodic_bcs(self) -> None:
        """Periodic BC config → Wolfram Association without value key."""
        wls = self._generate({
            "constraint_solver": {
                "enabled": True,
                "boundary_conditions": {
                    "x": {"type": "periodic"},
                    "y": {"type": "periodic"},
                },
            }
        })
        assert '"type" -> "periodic"' in wls
        # Periodic BCs should NOT have a "value" key
        assert '"value"' not in wls


class TestGaugeToml:
    """Tests for ``[[gauge]]`` TOML → WLS generation."""

    def _generate(self, config: dict[str, object]) -> str:
        from tidal.cli._derive import generate_wls

        base: dict[str, object] = {
            "theory": {"name": "Test"},
            "spacetime": {"dimension": 2, "metric": "minkowski"},
            "fields": [{"name": "A", "type": "vector"}],
            "derived_fields": [
                {
                    "name": "F",
                    "type": "tensor",
                    "rank": 2,
                    "symmetry": "antisymmetric",
                    "definition": "CD[-a][A[-b]] - CD[-b][A[-a]]",
                }
            ],
            "lagrangian": {
                "expression": "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d]"
            },
            "output": {"path": "out.json"},
        }
        base.update(config)
        return generate_wls(base)

    def test_lorenz_preset_in_wls(self) -> None:
        """Lorenz preset → BuildLorenzGaugeTerm + AddGaugeFixingTerm + GaugeFix.wl."""
        wls = self._generate({"gauge": [{"field": "A", "type": "lorenz"}]})
        assert "BuildLorenzGaugeTerm" in wls
        assert "AddGaugeFixingTerm" in wls
        assert "GaugeFix.wl" in wls

    def test_lorenz_xi_parameter(self) -> None:
        """Parameter xi=2.0 appears as fourth argument to BuildLorenzGaugeTerm."""
        import re

        wls = self._generate(
            {"gauge": [{"field": "A", "type": "lorenz", "xi": 2.0}]}
        )
        # xi=2.0 should appear as the fourth argument to the builder call
        assert re.search(r"BuildLorenzGaugeTerm\[.*,\s*2(\.0)?\s*\]", wls)

    def test_custom_expression_in_wls(self) -> None:
        """Custom gauge expression → GaugeTerm + AddGaugeFixingTerm in WLS."""
        wls = self._generate({
            "gauge": [{
                "field": "A",
                "type": "custom",
                "mechanism": "lagrangian_term",
                "expression": "-(1/2) * eta[a,b] CD[-a][A[-b]] eta[c,d] CD[-c][A[-d]]",
            }]
        })
        assert "GaugeTerm" in wls
        assert "AddGaugeFixingTerm" in wls

    def test_no_gauge_no_gaugefix(self) -> None:
        """No [[gauge]] → GaugeFix.wl not loaded."""
        wls = self._generate({})
        assert "GaugeFix.wl" not in wls
        assert '"gauge" -> "none"' in wls

    def test_gauge_metadata_string(self) -> None:
        """Lorenz gauge → metadata string 'lorenz(A)'."""
        wls = self._generate({"gauge": [{"field": "A", "type": "lorenz"}]})
        assert '"gauge" -> "lorenz(A)"' in wls

    def test_gauge_canonicalization_in_wls(self) -> None:
        """Lorenz gauge → ToCanonical + ContractMetric after AddGaugeFixingTerm."""
        wls = self._generate({"gauge": [{"field": "A", "type": "lorenz"}]})
        add_idx = wls.index("AddGaugeFixingTerm")
        canon_idx = wls.index("ToCanonical", add_idx)
        contract_idx = wls.index("ContractMetric", canon_idx)
        assert add_idx < canon_idx < contract_idx

    def test_temporal_substitution_in_wls(self) -> None:
        """Temporal gauge → fieldEquations /. {comp :> 0} in WLS output."""
        wls = self._generate({"gauge": [{"field": "A", "type": "temporal"}]})
        assert "fieldEquations /." in wls or "fieldEquations/." in wls
        assert ":> 0" in wls
        # GaugeFix.wl should NOT be loaded for Type B constraints
        assert "GaugeFix.wl" not in wls

    def test_coulomb_constraint_in_wls(self) -> None:
        """Coulomb gauge → AppendTo[fieldEquations, ...] with Derivative terms."""
        wls = self._generate({"gauge": [{"field": "A", "type": "coulomb"}]})
        assert "AppendTo[fieldEquations" in wls
        assert "Derivative" in wls
        assert "coulomb" in wls.lower()
        # GaugeFix.wl should NOT be loaded for Type B constraints
        assert "GaugeFix.wl" not in wls

    def test_axial_substitution_in_wls(self) -> None:
        """Axial gauge → fieldEquations /. {last-spatial-comp :> 0} in WLS output."""
        wls = self._generate({"gauge": [{"field": "A", "type": "axial"}]})
        assert "fieldEquations /." in wls or "fieldEquations/." in wls
        assert ":> 0" in wls

    def test_de_donder_in_wls(self) -> None:
        """De Donder gauge on tensor → BuildDeDonderGaugeTerm in WLS output."""
        wls = self._generate({
            "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "symmetric"}],
            "derived_fields": [],
            "lagrangian": {"expression": "h[-a, -b] eta[a, c] eta[b, d] h[-c, -d]"},
            "gauge": [{"field": "h", "type": "de_donder"}],
        })
        assert "BuildDeDonderGaugeTerm" in wls
        assert "GaugeFix.wl" in wls
        assert "ToCanonical" in wls

    def test_mixed_type_a_and_b_in_wls(self) -> None:
        """Mixed Type A + Type B: loads GaugeFix.wl for Type A, has substitution for Type B."""
        wls = self._generate({
            "fields": [
                {"name": "A", "type": "vector"},
                {"name": "B", "type": "vector"},
            ],
            "derived_fields": [],
            "lagrangian": {
                "expression": (
                    "-1/2 CD[-a][A[-b]] eta[a,c] eta[b,d] CD[-c][A[-d]] "
                    "- 1/2 CD[-a][B[-b]] eta[a,c] eta[b,d] CD[-c][B[-d]]"
                )
            },
            "gauge": [
                {"field": "A", "type": "lorenz"},
                {"field": "B", "type": "temporal"},
            ],
        })
        assert "BuildLorenzGaugeTerm" in wls
        assert "GaugeFix.wl" in wls
        assert ":> 0" in wls

    # --- TT gauge tests ---

    def test_tt_gauge_in_wls(self) -> None:
        """TT gauge on tensor → temporal zeroing, transverse constraints, traceless substitution."""
        wls = self._generate({
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "symmetric"}],
            "derived_fields": [],
            "lagrangian": {"expression": "h[-a, -b] eta[a, c] eta[b, d] h[-c, -d]"},
            "gauge": [{"field": "h", "type": "tt"}],
        })
        # Temporal zeroing: h_0, h_1, h_2 = 0 (dim=3, first 3 components)
        assert "TT-temporal" in wls
        assert ":> 0" in wls
        # Transverse constraints (must come BEFORE traceless so substitution catches them)
        assert "TT transverse" in wls
        assert "AppendTo[fieldEquations" in wls
        # Traceless substitution + Expand + constraint replacement
        assert "TT traceless" in wls
        assert "Expand[fieldEquations" in wls
        # Ordering: transverse before traceless (so h_{d-1,d-1} refs get substituted)
        assert wls.index("TT transverse") < wls.index("TT traceless")
        # GaugeFix.wl should NOT be loaded for Type B constraints
        assert "GaugeFix.wl" not in wls

    def test_tt_gauge_with_linearization(self) -> None:
        """TT gauge + linearization → gauge fixing applied in linearisation path."""
        wls = self._generate({
            "spacetime": {"dimension": 4, "metric": "minkowski"},
            "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "symmetric"}],
            "derived_fields": [],
            "lagrangian": {"expression": "RicciScalarCD[]"},
            "linearization": {"perturbation_field": "h"},
            "gauge": [{"field": "h", "type": "tt"}],
        })
        # Should contain both linearisation AND gauge-fixing code
        assert "Perturbation" in wls
        assert "TT-temporal" in wls
        assert "TT traceless" in wls
        assert "TT transverse" in wls
        # GaugeFix.wl should NOT be loaded (TT is Type B only)
        assert "GaugeFix.wl" not in wls

    def test_gauge_linearization_without_lagrangian_raises(self) -> None:
        """[[gauge]] + [linearization] without [lagrangian] → ValueError."""
        import pytest

        from tidal.cli._derive import generate_wls

        config: dict[str, object] = {
            "theory": {"name": "Test"},
            "spacetime": {"dimension": 2, "metric": "minkowski"},
            "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "symmetric"}],
            "linearization": {
                "perturbation_field": "h",
                "expression": "SomeExpression[]",
            },
            "gauge": [{"field": "h", "type": "tt"}],
            "output": {"path": "out.json"},
        }
        with pytest.raises(ValueError, match=r"requires \[lagrangian\]"):
            generate_wls(config)

    def test_gauge_metadata_tt(self) -> None:
        """TT gauge → metadata string 'tt(h)'."""
        wls = self._generate({
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "symmetric"}],
            "derived_fields": [],
            "lagrangian": {"expression": "h[-a, -b] eta[a, c] eta[b, d] h[-c, -d]"},
            "gauge": [{"field": "h", "type": "tt"}],
        })
        assert '"gauge" -> "tt(h)"' in wls

    def test_tt_gauge_rejects_vector(self) -> None:
        """TT gauge on vector field → ValueError."""
        import pytest

        with pytest.raises(ValueError, match="requires a tensor"):
            self._generate({"gauge": [{"field": "A", "type": "tt"}]})

    def test_tt_gauge_rejects_non_symmetric(self) -> None:
        """TT gauge on non-symmetric tensor → ValueError."""
        import pytest

        with pytest.raises(ValueError, match="requires a symmetric tensor"):
            self._generate({
                "fields": [{"name": "h", "type": "tensor", "rank": 2, "symmetry": "antisymmetric"}],
                "derived_fields": [],
                "lagrangian": {"expression": "h[-a, -b] eta[a, c] eta[b, d] h[-c, -d]"},
                "gauge": [{"field": "h", "type": "tt"}],
            })
