"""Tests for tidal.solver.validation — extracted spec validation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.validation import (
    check_cfl_stability,
    check_mass_sign,
    validate_field_references,
    validate_operator_dimensions,
)
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(
    terms: list[dict[str, Any]],
    *,
    dim: int = 2,
    field_name: str = "phi_0",
    order: int = 2,
) -> EquationSystem:
    """Build a minimal spec."""
    signature = [-1] + [1] * (dim - 1)
    lhs_expr = "d2_t(phi_0)" if order == 2 else "0"
    data: dict[str, Any] = {
        "spacetime": {"dimension": dim, "signature": signature},
        "fields": [{"name": field_name, "index": 0}],
        "equations": [
            {
                "field": field_name,
                "lhs": {"expression": lhs_expr, "order": {"time": order}},
                "rhs": {"type": "linear_combination", "terms": terms},
            }
        ],
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# Tests — operator dimensions
# ---------------------------------------------------------------------------


class TestOperatorDimensions:
    def test_laplacian_in_1d_passes(self) -> None:
        """Laplacian works in any dimension."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
            dim=2,
        )
        validate_operator_dimensions(spec)  # Should not raise

    def test_gradient_z_in_2d_raises(self) -> None:
        """gradient_z requires 3D but spec is 2D."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_z", "field": "phi_0"}],
            dim=3,  # 2+1D → spatial_dim = 2
        )
        with pytest.raises(ValueError, match="gradient_z"):
            validate_operator_dimensions(spec)

    def test_gradient_y_in_2d_passes(self) -> None:
        """gradient_y requires 2D, spec has spatial_dim=2."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_y", "field": "phi_0"}],
            dim=3,
        )
        validate_operator_dimensions(spec)  # Should not raise

    def test_gradient_x_in_1d_passes(self) -> None:
        """gradient_x requires 1D, spec has spatial_dim=1."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "gradient_x", "field": "phi_0"}],
            dim=2,
        )
        validate_operator_dimensions(spec)  # Should not raise


# ---------------------------------------------------------------------------
# Tests — field references
# ---------------------------------------------------------------------------


class TestFieldReferences:
    def test_valid_self_reference(self) -> None:
        """Term referencing its own field passes."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "identity", "field": "phi_0"}],
        )
        validate_field_references(spec)  # Should not raise

    def test_valid_cross_reference(self) -> None:
        """Cross-field reference between known fields passes."""
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
                                "operator": "identity",
                                "field": "phi_0",
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        validate_field_references(spec)  # Should not raise


# ---------------------------------------------------------------------------
# Tests — CFL stability
# ---------------------------------------------------------------------------


class TestCFLStability:
    def test_cfl_satisfied(self) -> None:
        """Small enough dt → no warnings."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
        )
        grid = GridInfo(
            bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,)
        )
        dx = grid.dx[0]
        dt = dx * 0.5  # Well within CFL
        warnings = check_cfl_stability(spec, grid, dt)
        assert len(warnings) == 0

    def test_cfl_violated(self) -> None:
        """Dt too large → warning returned."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
        )
        grid = GridInfo(
            bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,)
        )
        dt = 10.0  # Way too large
        warnings = check_cfl_stability(spec, grid, dt)
        assert len(warnings) == 1
        assert "CFL violation" in warnings[0]

    def test_constraint_ignored(self) -> None:
        """Constraint equations (order=0) are not CFL-checked."""
        spec = _make_spec(
            [{"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"}],
            order=0,
        )
        grid = GridInfo(
            bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,)
        )
        warnings = check_cfl_stability(spec, grid, dt=100.0)
        assert len(warnings) == 0


# ---------------------------------------------------------------------------
# Tests — mass sign diagnostic
# ---------------------------------------------------------------------------


class TestMassSign:
    def test_sign_change_detected(self) -> None:
        """Sin[x] mass term changes sign → warning returned."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sin[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(
            bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,)
        )
        import warnings as w

        with w.catch_warnings():
            w.simplefilter("ignore")
            coeff_eval = CoefficientEvaluator(spec, grid)
        result = check_mass_sign(coeff_eval, spec)
        assert len(result) == 1
        assert "changes sign" in result[0]

    def test_positive_definite_no_warning(self) -> None:
        """Always-positive mass term → no warning."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]^2 + 1",
                    "coordinate_dependent": ["x"],
                },
            ],
        )
        grid = GridInfo(
            bounds=((0, 2 * np.pi),), shape=(32,), periodic=(True,)
        )
        coeff_eval = CoefficientEvaluator(spec, grid)
        result = check_mass_sign(coeff_eval, spec)
        assert len(result) == 0
