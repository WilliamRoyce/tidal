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
        spec = _make_spec(
            [
                {"coefficient": 1.5, "operator": "laplacian", "field": "phi_0"},
            ]
        )
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid)

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == 1.5

    def test_parameter_override(self) -> None:
        """Symbolic resolved from parameters dict."""
        spec = _make_spec(
            [
                {
                    "coefficient": 0.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "m2",
                },
            ]
        )
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid, parameters={"m2": 2.5})

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == 2.5

    def test_negated_parameter(self) -> None:
        """Symbolic ``"-m2"`` resolves to negated parameter value."""
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
        grid = _make_grid_1d()
        ev = CoefficientEvaluator(spec, grid, parameters={"m2": 3.0})

        term = spec.equations[0].rhs_terms[0]
        result = ev.resolve(term, eq_idx=0, term_idx=0)
        assert result == -3.0

    def test_compound_expression(self) -> None:
        """Symbolic ``"-2*m2"`` resolves as compound expression."""
        spec = _make_spec(
            [
                {
                    "coefficient": 0.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "-2*m2",
                },
            ]
        )
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
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Exp[-t]",
                    "time_dependent": True,
                    "coordinate_dependent": ["t"],
                },
            ]
        )
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
        spec = _make_spec(
            [
                {
                    "coefficient": 0.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "unknown_param",
                },
            ]
        )
        grid = _make_grid_1d()
        # Should raise during construction (pre-resolve attempt)
        with pytest.raises(ValueError, match="unknown_param"):
            CoefficientEvaluator(spec, grid)

    def test_nan_expression(self) -> None:
        """Expression producing NaN → ValueError."""
        spec = _make_spec(
            [
                {
                    "coefficient": 0.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sqrt[-1]",
                },
            ]
        )
        grid = _make_grid_1d()
        # sqrt(-1) gives NaN for real-valued numpy
        with pytest.raises((ValueError, TypeError)):
            CoefficientEvaluator(spec, grid)

    def test_begin_timestep_clears_cache(self) -> None:
        """L3 timestep cache is cleared by begin_timestep()."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "identity",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[t]",
                    "time_dependent": True,
                    "coordinate_dependent": ["t"],
                },
            ]
        )
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


# ---------------------------------------------------------------------------
# Tests — periodic coefficient continuity check
# ---------------------------------------------------------------------------


class TestPeriodicCoefficientContinuity:
    """Tests for check_periodic_coefficient_continuity().

    This check prevents the subtle O(1) energy non-conservation that occurs
    when position-dependent coefficients are discontinuous at periodic
    boundaries.  The integration-by-parts identity requires [β(L)-β(0)]=0;
    when it doesn't hold, the boundary term leaks energy regardless of
    grid resolution or solver tolerance.
    """

    def test_periodic_coefficient_no_warning(self) -> None:
        """Periodic coefficient (Cos) on periodic domain: no warning.

        Cos[x] on [0, 2π] is exactly periodic: Cos(0)=Cos(2π)=1.
        On a cell-centered grid [dx/2, ..., 2π-dx/2], the first and last
        values are Cos(dx/2) ≈ Cos(2π-dx/2) by symmetry.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "Cos[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(32)
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((True,))

    def test_nonperiodic_coefficient_warns(self) -> None:
        """Non-periodic coefficient with moderate jump (1-50%): warns.

        Use a coefficient that varies smoothly but has a small mismatch
        at the periodic boundary.  Sin[x] on cell-centered [0,2π] grid
        has sin(dx/2) ≈ 0.098 at left and sin(2π-dx/2) ≈ -0.098 at right,
        giving a ~20% relative jump (warn, not error).
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "Sin[x[]]",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(32)
        ev = CoefficientEvaluator(spec, grid)
        with pytest.warns(UserWarning, match="jump at the periodic boundary"):
            ev.check_periodic_coefficient_continuity((True,))

    def test_large_jump_raises_error(self) -> None:
        """Coefficient with >50% relative jump raises ValueError.

        x² on [0, 2π] has x(0)²≈0 but x(2π)²≈39.5 → 100% relative jump.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[]^2",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(32)
        ev = CoefficientEvaluator(spec, grid)
        with pytest.raises(ValueError, match="jump at the periodic boundary"):
            ev.check_periodic_coefficient_continuity((True,))

    def test_divergent_coefficient_raises_error(self) -> None:
        """1/x on [0.1, 10]: left=10 vs right=0.1, 99% jump → error.

        This models the original dipolar B∝1/r³ failure case.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "1/x[]",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = GridInfo(bounds=((0.1, 10.0),), shape=(32,), periodic=(True,))
        ev = CoefficientEvaluator(spec, grid)
        with pytest.raises(ValueError, match="jump at the periodic boundary"):
            ev.check_periodic_coefficient_continuity((True,))

    def test_no_check_when_not_periodic(self) -> None:
        """No warning even with non-periodic coefficient if BCs aren't periodic."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[]^2",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(32)
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((False,))

    def test_small_absolute_jump_benign(self) -> None:
        """Coefficient periodic in B(x) but with x² factor: small absolute jump.

        This is the centered-dipolar scenario: B(x) is perfectly periodic but
        the Wolfram expansion produces B²·x² terms where x² breaks periodicity.
        The key insight is that |coeff(boundary)| ≪ |coeff(center)|, so the
        leak metric (rel_jump * boundary_fraction) is negligible.
        """
        # Simulate B²·x² where B is tiny at boundaries
        # B = 0.001 at x=1, B = 1.0 at x=5 (center), B = 0.001 at x=9
        # coeff = B² * x² → at x=1: 1e-6, at x=9: 8.1e-5
        # Relative jump = |8.1e-5 - 1e-6| / max ≈ 0.3% — below 1% threshold
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    # Coefficient that's small at both ends, large at center
                    "coefficient_symbolic": ("x[]^2 / (1 + (x[] - 3.14159)^2)^3"),
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(64)
        ev = CoefficientEvaluator(spec, grid)
        # This coefficient has tiny values at both boundaries despite
        # containing x² — the localized denominator suppresses them.
        # Verify it doesn't raise (absolute values at boundary are tiny).
        ev.check_periodic_coefficient_continuity((True,))

    def test_antisymmetric_localized_no_warning(self) -> None:
        """Antisymmetric localized coefficient: large relative jump but negligible boundaries.

        x·exp(-x²/R²) is antisymmetric — opposite signs at boundaries.
        On [0, 2π] with R=1 (center at π): boundary values are ~10⁻⁴ of peak,
        so the relative jump appears large (~200%) but the absolute boundary
        magnitude is negligible.  The boundary-magnitude filter should suppress
        this false positive.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": (
                        "(x[] - 3.14159) * exp(-(x[] - 3.14159)^2)"
                    ),
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(64)
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((True,))

    def test_localized_gaussian_periodic_no_warning(self) -> None:
        """Gaussian exp(-x²/R²) on large periodic domain: zero at boundaries.

        This models the Gertsenshtein localized B-field on a periodic domain.
        The Gaussian decays to float64 zero at the boundaries, so no warning
        should trigger regardless of domain size.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "exp(-x[]^2 / 50)",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        # Domain [-50, 50] with R=5: exp(-50²/50) ≈ 10⁻²²
        grid = GridInfo(bounds=((-50.0, 50.0),), shape=(64,), periodic=(True,))
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((True,))

    def test_centered_dipolar_no_warning(self) -> None:
        """Centered dipolar B²·x²: localized with small boundary values.

        Models the Gertsenshtein false-positive scenario.  The coefficient
        peaks at ~10 in the centre but boundary values are ~0.02.  The old
        relative-jump check warned (rel_jump ~2-6%), but the leak metric
        (rel_jump * boundary_fraction) is ~8e-6, consistent with the
        observed |dE/E| = 4.53e-7 in actual simulations.
        """
        # Simulate B(x)² * x² where B is a Gaussian centred at x=π
        # on [0, 2π].  The x² factor breaks strict periodicity but both
        # boundary values are negligible (B(0)² * 0² ≈ 0, B(2π)² * (2π)² ≈ tiny).
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": ("x[]^2 * exp(-(x[] - 3.14159)^2 / 0.5)"),
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(128)
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((True,))

    def test_genuine_discontinuity_still_errors(self) -> None:
        """x² on [0, 2π] with large boundary values: still detected.

        Verifies the boundary-magnitude filter doesn't weaken detection of
        real discontinuities where boundary values are significant.
        x(0)²≈0 but x(2π)²≈39.5 → 100% relative jump with large |last|.
        """
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[]^2",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(32)
        ev = CoefficientEvaluator(spec, grid)
        with pytest.raises(ValueError, match="jump at the periodic boundary"):
            ev.check_periodic_coefficient_continuity((True,))

    def test_2d_checks_periodic_axis(self) -> None:
        """In 2D, check detects discontinuity along the periodic x-axis."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[]^2",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=3,
        )
        grid = _make_grid_2d(8)
        ev = CoefficientEvaluator(spec, grid)
        # x² on [0, 2π] has 100% relative jump → hard error
        with pytest.raises(ValueError, match="jump at the periodic boundary"):
            ev.check_periodic_coefficient_continuity((True, True))

    def test_2d_nonperiodic_axis_skipped(self) -> None:
        """In 2D, non-periodic x-axis is skipped even if coefficient varies."""
        spec = _make_spec(
            [
                {
                    "coefficient": 1.0,
                    "operator": "laplacian",
                    "field": "phi_0",
                    "coefficient_symbolic": "x[]^2",
                    "coordinate_dependent": ["x"],
                },
            ],
            dim=3,
        )
        grid = _make_grid_2d(8)
        ev = CoefficientEvaluator(spec, grid)
        # Only y-axis is periodic; x-axis not checked
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((False, True))

    def test_constant_coefficient_no_check(self) -> None:
        """Constant coefficients aren't in _spatial_cache, so no check needed."""
        spec = _make_spec(
            [
                {
                    "coefficient": 2.5,
                    "operator": "laplacian",
                    "field": "phi_0",
                },
            ],
            dim=2,
        )
        grid = _make_grid_1d(16)
        ev = CoefficientEvaluator(spec, grid)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            ev.check_periodic_coefficient_continuity((True,))
