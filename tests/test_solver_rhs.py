"""Tests for tidal.solver.rhs — RHSEvaluator unified operator application."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.fields import FieldSet
from tidal.solver.grid import GridInfo
from tidal.solver.rhs import RHSEvaluator
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(
    terms: list[dict[str, Any]],
    *,
    dim: int = 2,
    field_name: str = "phi_0",
) -> EquationSystem:
    """Build a minimal 2nd-order spec with given RHS terms."""
    signature = [-1] + [1] * (dim - 1)
    data: dict[str, Any] = {
        "spacetime": {"dimension": dim, "signature": signature},
        "fields": [{"name": field_name, "index": 0}],
        "equations": [
            {
                "field": field_name,
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {"type": "linear_combination", "terms": terms},
            }
        ],
    }
    return EquationSystem.from_dict(data)


def _make_grid(n: int = 32) -> GridInfo:
    return GridInfo(bounds=((0, 2 * np.pi),), shape=(n,), periodic=(True,))


def _make_evaluator(
    spec: EquationSystem,
    grid: GridInfo,
    params: dict[str, float] | None = None,
    bc: str | tuple[str, ...] | None = None,
) -> tuple[RHSEvaluator, CoefficientEvaluator]:
    """Build RHSEvaluator + CoefficientEvaluator pair."""
    coeff = CoefficientEvaluator(spec, grid, parameters=params)
    rhs = RHSEvaluator(spec, grid, coeff, bc=bc)
    return rhs, coeff


# ---------------------------------------------------------------------------
# Tests — basic evaluation
# ---------------------------------------------------------------------------


class TestRHSBasic:
    def test_single_laplacian(self) -> None:
        """Laplacian of sin(x) should give -sin(x)."""
        spec = _make_spec(
            [
                {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
            ]
        )
        n = 64
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        x = grid.axes_coords(0)
        fs["phi_0"] = np.sin(x)

        result = rhs_eval.evaluate(0, fs)
        expected = -np.sin(x)  # d²sin(x)/dx² = -sin(x)
        np.testing.assert_allclose(result, expected, atol=0.01)

    def test_multiple_terms(self) -> None:
        """Two terms: laplacian + identity (wave equation with mass)."""
        spec = _make_spec(
            [
                {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                {"coefficient": -2.0, "operator": "identity", "field": "phi_0"},
            ]
        )
        n = 64
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        x = grid.axes_coords(0)
        fs["phi_0"] = np.sin(x)

        result = rhs_eval.evaluate(0, fs)
        expected = -np.sin(x) + (-2.0) * np.sin(x)  # lap + mass
        np.testing.assert_allclose(result, expected, atol=0.01)

    def test_cross_field_term(self) -> None:
        """Term referencing a different field."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [
                {"name": "phi_0", "index": 0},
                {"name": "chi_0", "index": 1},
            ],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {
                        "expression": "d2_t(phi_0)",
                        "order": {"time": 2},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "phi_0",
                            },
                            {
                                "coefficient": 0.5,
                                "operator": "identity",
                                "field": "chi_0",
                            },
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {
                        "expression": "d2_t(chi_0)",
                        "order": {"time": 2},
                    },
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian",
                                "field": "chi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        n = 32
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        fs["chi_0"] = np.ones(n) * 3.0

        # eq 0: lap(phi) + 0.5 * chi; phi=0 → lap=0, chi=3 → 1.5
        result = rhs_eval.evaluate(0, fs)
        np.testing.assert_allclose(result, 1.5, atol=1e-10)

    def test_zero_field_gives_zero(self) -> None:
        """Zero-valued field produces zero RHS contribution."""
        spec = _make_spec(
            [
                {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
            ]
        )
        n = 16
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))  # phi_0 = 0

        result = rhs_eval.evaluate(0, fs)
        np.testing.assert_allclose(result, 0.0, atol=1e-14)


# ---------------------------------------------------------------------------
# Tests — coefficient integration
# ---------------------------------------------------------------------------


class TestRHSCoefficients:
    def test_parameter_override(self) -> None:
        """Parameter override flows through to RHS result."""
        spec = _make_spec(
            [
                {
                    "coefficient": 0.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "-m2",
                },
            ]
        )
        n = 16
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid, params={"m2": 4.0})

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        fs["phi_0"] = np.ones(n)

        result = rhs_eval.evaluate(0, fs)
        np.testing.assert_allclose(result, -4.0, atol=1e-14)

    def test_position_dependent(self) -> None:
        """Position-dependent coefficient multiplied correctly."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ]
        )
        n = 32
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        fs["phi_0"] = np.ones(n)

        result = rhs_eval.evaluate(0, fs)
        x = grid.axes_coords(0)
        np.testing.assert_allclose(result, np.cos(x), atol=1e-14)


# ---------------------------------------------------------------------------
# Tests — evaluate_by_field
# ---------------------------------------------------------------------------


class TestEvaluateByField:
    def test_lookup_by_name(self) -> None:
        """evaluate_by_field looks up the equation for a field name."""
        spec = _make_spec(
            [
                {"coefficient": 2.0, "operator": "identity", "field": "phi_0"},
            ]
        )
        n = 8
        grid = _make_grid(n)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (n,))
        fs["phi_0"] = np.ones(n)

        result = rhs_eval.evaluate_by_field("phi_0", fs)
        np.testing.assert_allclose(result, 2.0, atol=1e-14)

    def test_unknown_field_raises(self) -> None:
        """evaluate_by_field raises KeyError for unknown field."""
        spec = _make_spec(
            [
                {"coefficient": 1.0, "operator": "identity", "field": "phi_0"},
            ]
        )
        grid = _make_grid(8)
        rhs_eval, _ = _make_evaluator(spec, grid)

        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        with pytest.raises(KeyError, match="nonexistent"):
            rhs_eval.evaluate_by_field("nonexistent", fs)
