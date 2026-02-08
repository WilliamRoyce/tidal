"""Unit tests for CLI parsing helpers.

Tests the pure-function parsers in ``torsion_gertsenshtein.cli._simulate``
independently from the full simulation pipeline.
"""

from __future__ import annotations

import pytest

from torsion_gertsenshtein.cli._simulate import (
    _parse_bc,
    _parse_bounds,
    _parse_grid_shape,
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

    def test_inverted_bounds_parsed(self) -> None:
        """Inverted bounds parse successfully (validation is elsewhere)."""
        lo, hi = _parse_single_bound("20:0")
        assert lo == 20.0
        assert hi == 0.0


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
            ("dirichlet", True, 1, [False]),
            ("neumann,periodic", False, 2, [False, True]),
            ("periodic,neumann,dirichlet", True, 3, [True, False, False]),
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
