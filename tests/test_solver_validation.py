"""Tests for tidal.solver.validation — extracted spec validation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.solver.operators import AxisBCSpec, SideBCSpec
from tidal.solver.validation import (
    check_cfl_stability,
    check_mass_sign,
    check_pointwise_mass_stability,
    check_robin_stability,
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


# ---------------------------------------------------------------------------
# check_robin_stability
# ---------------------------------------------------------------------------


class TestCheckRobinStability:
    """Tests for Robin BC ghost-cell stability check."""

    def test_no_axis_bcs_returns_empty(self) -> None:
        grid = GridInfo(bounds=((0, 10),), shape=(64,), periodic=(False,))
        assert check_robin_stability(grid) == []

    def test_stable_robin(self) -> None:
        """gamma*dx < 2 should produce no warnings."""
        side = SideBCSpec(kind="robin", gamma=1.0, value=0.0)
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),), shape=(64,), periodic=(False,),
            axis_bcs=(abc,),
        )
        result = check_robin_stability(grid)
        assert len(result) == 0

    def test_unstable_robin_warns(self) -> None:
        """gamma*dx >= 2 should produce a warning."""
        # dx = 10/4 = 2.5, gamma = 1.0 -> gamma*dx = 2.5 >= 2
        side = SideBCSpec(kind="robin", gamma=1.0, value=0.0)
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),), shape=(4,), periodic=(False,),
            axis_bcs=(abc,),
        )
        result = check_robin_stability(grid)
        assert len(result) > 0
        assert "unstable" in result[0].lower()

    def test_periodic_axis_skipped(self) -> None:
        """Periodic axes are not checked for Robin stability."""
        abc = AxisBCSpec(periodic=True)
        grid = GridInfo(
            bounds=((0, 10),), shape=(4,), periodic=(True,),
            axis_bcs=(abc,),
        )
        assert check_robin_stability(grid) == []

    def test_neumann_not_warned(self) -> None:
        """Neumann BC should not trigger Robin stability warning."""
        side = SideBCSpec(kind="neumann")
        abc = AxisBCSpec(periodic=False, low=side, high=side)
        grid = GridInfo(
            bounds=((0, 10),), shape=(4,), periodic=(False,),
            axis_bcs=(abc,),
        )
        assert check_robin_stability(grid) == []


# ---------------------------------------------------------------------------
# Tests — pointwise mass stability
# ---------------------------------------------------------------------------


def _make_coupled_spec(m_phi2: float, m_chi2: float, g0: float) -> EquationSystem:
    """Two coupled scalar fields with constant cross-coupling."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}, {"name": "chi_0", "index": 1}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": -m_phi2, "operator": "identity", "field": "phi_0"},
                        {"coefficient": -g0, "operator": "identity", "field": "chi_0"},
                    ],
                },
            },
            {
                "field": "chi_0",
                "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": -m_chi2, "operator": "identity", "field": "chi_0"},
                        {"coefficient": -g0, "operator": "identity", "field": "phi_0"},
                    ],
                },
            },
        ],
    }
    return EquationSystem.from_dict(data)


class TestPointwiseMassStability:
    def test_stable_diagonal(self) -> None:
        """Decoupled system with positive masses is stable."""
        spec = _make_coupled_spec(m_phi2=1.0, m_chi2=1.0, g0=0.0)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert result.errors == []
        assert result.notes == []

    def test_stable_coupled_constant(self) -> None:
        """Coupled system satisfying det > 0 (mPhi2*mChi2 > g0^2) is stable."""
        spec = _make_coupled_spec(m_phi2=1.0, m_chi2=4.0, g0=1.0)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert result.errors == []
        assert result.notes == []

    def test_unstable_coupled_constant(self) -> None:
        """Coupled system violating det > 0 (g0^2 > mPhi2*mChi2) is unstable."""
        # det = 1.0*0.5 - 1.0^2 = -0.5 < 0
        spec = _make_coupled_spec(m_phi2=1.0, m_chi2=0.5, g0=1.0)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert len(result.errors) == 1
        assert "eigenvalue" in result.errors[0]

    def test_unstable_boundary_condition(self) -> None:
        """g0 = sqrt(mPhi2*mChi2) exactly is on the stability boundary (det=0)."""
        import math
        g0 = math.sqrt(1.0 * 1.0)  # exactly on boundary → det = 0, min eigenvalue = 0
        spec = _make_coupled_spec(m_phi2=1.0, m_chi2=1.0, g0=g0)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        # det=0: zero eigenvalue — within tolerance, so no warning expected
        assert result.errors == []

    def test_strongly_unstable(self) -> None:
        """Large coupling completely dominates — clearly unstable."""
        spec = _make_coupled_spec(m_phi2=0.1, m_chi2=0.1, g0=5.0)
        grid = GridInfo(bounds=((0, 10),), shape=(4,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert len(result.errors) == 1
        assert "unstable" in result.errors[0].lower()

    def test_single_field_stable(self) -> None:
        """Single field with positive mass is stable (1x1 matrix)."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(8,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert result.errors == []

    def test_2d_grid(self) -> None:
        """Stability check works on a 2D grid (broadcast correctness)."""
        spec = _make_coupled_spec(m_phi2=1.0, m_chi2=1.0, g0=0.5)
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True),
        )
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        assert result.errors == []

    def test_constraint_equation_no_false_positive(self) -> None:
        """Constraint field (time_order=0) should not cause spurious warnings."""
        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [
                {"name": "A_0", "index": 0},
                {"name": "A_1", "index": 1},
            ],
            "equations": [
                {
                    "field": "A_0",
                    "lhs": {"expression": "0", "order": {"time": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"},
                        ],
                    },
                },
                {
                    "field": "A_1",
                    "lhs": {"expression": "d2_t(A_1)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "A_1"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(
            bounds=((0, 10), (0, 10)), shape=(4, 4), periodic=(True, True),
        )
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        # A_0 has no identity term; A_1 has positive mass. Should pass.
        assert result.errors == []

    def test_asymmetric_stable_no_instability_error(self) -> None:
        """Asymmetric but stable coupling: note emitted, no error."""
        # Matrix: [[1.0, 0.5], [0.3, 1.0]] — asymmetric, but eigenvalues
        # are 0.613 and 1.387, both positive → stable.
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}, {"name": "chi_0", "index": 1}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                            {"coefficient": -0.5, "operator": "identity", "field": "chi_0"},
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "chi_0"},
                            {"coefficient": -0.3, "operator": "identity", "field": "phi_0"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(4,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        # Asymmetry detected as informational note
        assert any("asymmetric" in n.lower() for n in result.notes)
        # But stable — no errors
        assert result.errors == []

    def test_asymmetric_unstable_detected(self) -> None:
        """Asymmetric and unstable coupling: both note and error emitted."""
        # Matrix: [[0.1, 0.5], [0.3, 0.1]] — asymmetric AND unstable
        # (coupling dominates diagonal).
        data: dict[str, Any] = {
            "spacetime": {"dimension": 2, "signature": [-1, 1]},
            "fields": [{"name": "phi_0", "index": 0}, {"name": "chi_0", "index": 1}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -0.1, "operator": "identity", "field": "phi_0"},
                            {"coefficient": -5.0, "operator": "identity", "field": "chi_0"},
                        ],
                    },
                },
                {
                    "field": "chi_0",
                    "lhs": {"expression": "d2_t(chi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -0.1, "operator": "identity", "field": "chi_0"},
                            {"coefficient": -3.0, "operator": "identity", "field": "phi_0"},
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(bounds=((0, 10),), shape=(4,), periodic=(True,))
        coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        # Asymmetry detected as informational note
        assert any("asymmetric" in n.lower() for n in result.notes)
        # Instability detected as error
        assert len(result.errors) == 1
        assert "unstable" in result.errors[0].lower()

    def test_2d_position_dependent_coupling(self) -> None:
        """Position-dependent coupling on a 2D grid (core use case)."""
        import warnings as w

        data: dict[str, Any] = {
            "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
            "fields": [{"name": "phi_0", "index": 0}, {"name": "chi_0", "index": 1}],
            "equations": [
                {
                    "field": "phi_0",
                    "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": -1.0, "operator": "identity", "field": "phi_0"},
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "chi_0",
                                "coefficient_symbolic": "-0.1*Exp[-(x[]^2 + y[]^2)/2]",
                                "coordinate_dependent": ["x", "y"],
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
                            {"coefficient": -1.0, "operator": "identity", "field": "chi_0"},
                            {
                                "coefficient": -1.0,
                                "operator": "identity",
                                "field": "phi_0",
                                "coefficient_symbolic": "-0.1*Exp[-(x[]^2 + y[]^2)/2]",
                                "coordinate_dependent": ["x", "y"],
                            },
                        ],
                    },
                },
            ],
        }
        spec = EquationSystem.from_dict(data)
        grid = GridInfo(
            bounds=((-5, 5), (-5, 5)), shape=(8, 8), periodic=(False, False),
        )
        with w.catch_warnings():
            w.simplefilter("ignore")
            coeff_eval = CoefficientEvaluator(spec, grid, {})
        result = check_pointwise_mass_stability(coeff_eval, spec, grid)
        # Gaussian coupling peak is 0.1 << 1.0 diagonal → stable
        assert result.errors == []
