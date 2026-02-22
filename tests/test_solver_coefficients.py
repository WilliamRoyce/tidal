"""Tests for tidal.solver.coefficients — CoefficientEvaluator."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest

from tidal.solver.coefficients import CoefficientEvaluator
from tidal.solver.grid import GridInfo
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_spec(
    terms: list[dict[str, Any]],
    *,
    dim: int = 2,
    signature: list[int] | None = None,
    field_name: str = "phi_0",
) -> EquationSystem:
    """Build a minimal EquationSystem with the given RHS terms."""
    if signature is None:
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


def _make_grid_1d(n: int = 32) -> GridInfo:
    return GridInfo(bounds=((0, 2 * np.pi),), shape=(n,), periodic=(True,))


def _make_grid_2d(n: int = 8) -> GridInfo:
    return GridInfo(
        bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
        shape=(n, n),
        periodic=(True, True),
    )


# ---------------------------------------------------------------------------
# Tests — constant coefficients
# ---------------------------------------------------------------------------


class TestConstantCoefficients:
    def test_numeric_only(self) -> None:
        """No symbolic → return term.coefficient directly."""
        spec = _make_spec([
            {"coefficient": 1.5, "operator": "laplacian", "field": "phi_0"},
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == 1.5

    def test_parameter_override(self) -> None:
        """Symbolic resolved from parameters dict."""
        spec = _make_spec([
            {
                "coefficient": 0.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "m2",
            },
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid, parameters={"m2": 2.5})

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == 2.5

    def test_negated_parameter(self) -> None:
        """Symbolic ``"-m2"`` resolves to negated parameter value."""
        spec = _make_spec([
            {
                "coefficient": 0.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "-m2",
            },
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid, parameters={"m2": 3.0})

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == -3.0

    def test_compound_expression(self) -> None:
        """Symbolic ``"-2*m2"`` resolves as compound expression."""
        spec = _make_spec([
            {
                "coefficient": 0.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "-2*m2",
            },
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid, parameters={"m2": 1.5})

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == -3.0


# ---------------------------------------------------------------------------
# Tests — position-dependent coefficients
# ---------------------------------------------------------------------------


class TestPositionDependent:
    def test_spatial_coefficient(self) -> None:
        """Position-dependent coeff returns grid-shaped array."""
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
            dim=2,
        )
        grid = _make_grid_1d(16)
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert isinstance(result, np.ndarray)
        assert result.shape == (16,)

        x = grid.axes_coords(0)
        np.testing.assert_allclose(result, np.sin(x), atol=1e-14)

    def test_spatial_cached(self) -> None:
        """Spatial-only coefficients are cached (same object on repeat)."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(8)
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        r1 = ev.resolve(term, eq_idx=0, term_idx=0)
        r2 = ev.resolve(term, eq_idx=0, term_idx=0)
        assert r1 is r2  # Same array object — cached

    def test_2d_spatial(self) -> None:
        """Position-dependent in 2D returns correctly shaped array."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[] + y[]",
                    "coordinate_dependent": ["x", "y"],
                },
            ],
            dim=3,
        )
        grid = _make_grid_2d(4)
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert isinstance(result, np.ndarray)
        assert result.shape == (4, 4)


# ---------------------------------------------------------------------------
# Tests — time-dependent coefficients
# ---------------------------------------------------------------------------


class TestTimeDependent:
    def test_varies_with_t(self) -> None:
        """Time-dependent coefficient changes with t."""
        spec = _make_spec([
            {
                "coefficient": 1.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "Exp[-t]",
                "time_dependent": True,
                "coordinate_dependent": ["t"],
            },
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        r0 = ev.resolve(term, t=0.0, eq_idx=0, term_idx=0)
        np.testing.assert_allclose(r0, 1.0, atol=1e-14)

        ev.begin_timestep(1.0)
        r1 = ev.resolve(term, t=1.0, eq_idx=0, term_idx=0)
        np.testing.assert_allclose(r1, np.exp(-1.0), atol=1e-14)


# ---------------------------------------------------------------------------
# Tests — error handling
# ---------------------------------------------------------------------------


class TestErrorHandling:
    def test_missing_parameter(self) -> None:
        """Missing parameter in symbolic → ValueError."""
        spec = _make_spec([
            {
                "coefficient": 0.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "unknown_param",
            },
        ])
        grid = _make_grid_1d()
        # Should raise during construction (pre-resolve attempt)
        with pytest.raises(ValueError, match="unknown_param"):
            CoefficientEvaluator(spec, grid)

    def test_nan_expression(self) -> None:
        """Expression producing NaN → ValueError."""
        spec = _make_spec([
            {
                "coefficient": 0.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "Sqrt[-1]",
            },
        ])
        grid = _make_grid_1d()
        # sqrt(-1) gives NaN for real-valued numpy
        with pytest.raises((ValueError, TypeError)):
            CoefficientEvaluator(spec, grid)

    def test_begin_timestep_clears_cache(self) -> None:
        """L3 timestep cache is cleared by begin_timestep()."""
        spec = _make_spec([
            {
                "coefficient": 1.0,
                "operator": "identity",
                "field": "phi_0",
                "coefficient_symbolic": "Cos[t]",
                "time_dependent": True,
                "coordinate_dependent": ["t"],
            },
        ])
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        ev.resolve(term, t=0.0, eq_idx=0, term_idx=0)
        assert len(ev._timestep_cache) == 1

        ev.begin_timestep(1.0)
        assert len(ev._timestep_cache) == 0


# ---------------------------------------------------------------------------
# Tests — mass sign diagnostic
# ---------------------------------------------------------------------------


class TestMassSignDiagnostic:
    def test_sign_change_warns(self) -> None:
        """Mass term that changes sign across grid emits UserWarning."""
        # Sin[x] changes sign: positive on (0,π), negative on (π,2π)
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
            dim=2,
        )
        grid = _make_grid_1d(32)
        with pytest.warns(UserWarning, match="tachyonic"):
            CoefficientEvaluator(spec, grid)

    def test_no_warning_for_positive_mass(self) -> None:
        """Positive-definite mass term does not warn."""
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
            dim=2,
        )
        grid = _make_grid_1d(32)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            CoefficientEvaluator(spec, grid)  # Should not warn
