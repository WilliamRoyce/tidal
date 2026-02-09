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
        assert _infer_output_format(_make_args(output_format="npz")) == "npz"

    def test_npz_extension(self) -> None:
        assert _infer_output_format(_make_args(output="foo.npz")) == "npz"

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
        assert _infer_output_format(_make_args(no_plot=True, output_format="npz")) == "summary"


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
