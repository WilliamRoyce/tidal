"""Stress tests for background fields feature.

Tests the full pipeline from _mathematica_to_python conversion through
coefficient resolution to energy computation.  These tests exercise the
code paths that background field coefficients trigger — including
position-dependent expressions like UnitStep, Sign, Max, Min.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid

from tidal.measurement._io import SimulationData
from tidal.symbolic.json_loader import (
    ComponentEquation,
    ConstraintSolverConfig,
    EquationSystem,
    OperatorTerm,
)
from tidal.symbolic.pde_builder import PDEFromSpec

if TYPE_CHECKING:
    import pathlib


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_spec(
    *,
    coefficient_symbolic: str = "1.0",
    coordinate_dependent: tuple[str, ...] = (),
    operator: str = "laplacian",
    field: str = "phi",
    coefficient: float = 1.0,
) -> EquationSystem:
    """Build a minimal 1-field EquationSystem for coefficient testing."""
    return EquationSystem(
        n_components=1,
        dimension=3,
        spatial_dimension=2,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(
                        coefficient,
                        operator,
                        field,
                        coefficient_symbolic=coefficient_symbolic,
                        coordinate_dependent=coordinate_dependent,
                    ),
                ),
            ),
        ),
        mass_matrix=((0.0,),),
        coupling_matrix=((0.0,),),
        metadata={},
        coordinates=("t", "x", "y"),
    )


def _make_pde(
    spec: EquationSystem,
    parameters: dict[str, float] | None = None,
) -> PDEFromSpec:
    return PDEFromSpec(spec, parameters=parameters or {})


# ===========================================================================
# Group 1: _mathematica_to_python conversion — Bug A + Bug B fixes
# ===========================================================================


class TestMathematicaConversionBugFixes:
    """Verify UnitStep, HeavisideTheta, Sign, Max, Min conversions."""

    def make_pde(self) -> PDEFromSpec:
        spec = _make_spec()
        return _make_pde(spec)

    def test_unitstep_conversion(self) -> None:
        """UnitStep[x] → heaviside(x)."""
        pde = self.make_pde()
        result = pde._mathematica_to_python("UnitStep[x - 20]")
        assert result == "heaviside(x - 20)"

    def test_unitstep_nested(self) -> None:
        """Product of UnitStep calls for localized background."""
        pde = self.make_pde()
        result = pde._mathematica_to_python(
            "V0 * UnitStep[x() - 30] * UnitStep[70 - x()]"
        )
        assert "heaviside" in result
        assert "V0" in result

    def test_heavisidetheta_conversion(self) -> None:
        """HeavisideTheta[x] → heaviside(x)."""
        pde = self.make_pde()
        result = pde._mathematica_to_python("HeavisideTheta[x - 5]")
        assert result == "heaviside(x - 5)"

    def test_sign_conversion(self) -> None:
        """Sign[x] → sign(x), not np.sign(x)."""
        pde = self.make_pde()
        result = pde._mathematica_to_python("Sign[x]")
        assert result == "sign(x)"
        assert "np." not in result

    def test_max_conversion(self) -> None:
        """Max[x, 0] → maximum(x, 0), not np.maximum(x, 0)."""
        pde = self.make_pde()
        result = pde._mathematica_to_python("Max[x, 0]")
        assert result == "maximum(x, 0)"
        assert "np." not in result

    def test_min_conversion(self) -> None:
        """Min[x, 1] → minimum(x, 1), not np.minimum(x, 1)."""
        pde = self.make_pde()
        result = pde._mathematica_to_python("Min[x, 1]")
        assert result == "minimum(x, 1)"
        assert "np." not in result


# ===========================================================================
# Group 2: Coefficient resolution with background-like expressions
# ===========================================================================


class TestBackgroundCoefficientResolution:
    """Verify coefficient resolution for expressions that background fields produce."""

    def test_constant_symbolic_coefficient(self) -> None:
        """A pure symbolic constant (e.g., 'Bval') resolves to the parameter value."""
        spec = _make_spec(
            coefficient_symbolic="Bval",
            operator="identity",
            coefficient=0.0,
        )
        pde = _make_pde(spec, parameters={"Bval": 2.5})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=None)
        assert result == pytest.approx(2.5)

    def test_product_symbolic_coefficient(self) -> None:
        """Product of two symbolic constants resolves correctly."""
        spec = _make_spec(
            coefficient_symbolic="gCpl*Bval",
            operator="identity",
            coefficient=0.0,
        )
        pde = _make_pde(spec, parameters={"gCpl": 0.3, "Bval": 2.0})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=None)
        assert result == pytest.approx(0.6)

    def test_negated_symbolic_coefficient(self) -> None:
        """Negated symbolic coefficient (e.g., '-V0') resolves correctly."""
        spec = _make_spec(
            coefficient_symbolic="-V0",
            operator="identity",
            coefficient=0.0,
        )
        pde = _make_pde(spec, parameters={"V0": 4.0})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=None)
        assert result == pytest.approx(-4.0)

    def test_position_dependent_unitstep(self) -> None:
        """UnitStep coefficient evaluates to step-function array."""
        spec = _make_spec(
            coefficient_symbolic="V0*UnitStep[x() - 3]",
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=0.0,
        )
        grid = CartesianGrid(bounds=[(0, 10), (0, 10)], shape=[20, 20], periodic=True)
        pde = _make_pde(spec, parameters={"V0": 5.0})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

        assert isinstance(result, np.ndarray)
        assert result.shape == (20, 20)

        # x < 3: UnitStep = 0, so coefficient = 0
        x_coords = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]
        mask_below = x_coords < 3.0
        mask_above = x_coords > 3.0
        assert_allclose(result[mask_below], 0.0, atol=1e-10)
        assert_allclose(result[mask_above], 5.0, atol=1e-10)

    def test_position_dependent_double_unitstep(self) -> None:
        """Two UnitStep calls create a localized step-function region."""
        spec = _make_spec(
            coefficient_symbolic="V0*UnitStep[x() - 3]*UnitStep[7 - x()]",
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=0.0,
        )
        grid = CartesianGrid(bounds=[(0, 10), (0, 10)], shape=[40, 40], periodic=True)
        pde = _make_pde(spec, parameters={"V0": 2.0})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

        assert isinstance(result, np.ndarray)
        x_coords = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]

        # Inside [3, 7]: coefficient = V0 = 2.0
        mask_inside = (x_coords > 3.0) & (x_coords < 7.0)
        # Outside the well: coefficient = 0
        mask_outside = (x_coords < 3.0) | (x_coords > 7.0)
        assert_allclose(result[mask_inside], 2.0, atol=1e-10)
        assert_allclose(result[mask_outside], 0.0, atol=1e-10)

    def test_position_dependent_sign(self) -> None:
        """Sign coefficient evaluates to ±1 array."""
        spec = _make_spec(
            coefficient_symbolic="Sign[x() - 5]",
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=0.0,
        )
        grid = CartesianGrid(bounds=[(0, 10), (0, 10)], shape=[20, 20], periodic=True)
        pde = _make_pde(spec, parameters={})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

        assert isinstance(result, np.ndarray)
        x_coords = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]
        mask_below = x_coords < 5.0
        mask_above = x_coords > 5.0
        assert_allclose(result[mask_below], -1.0, atol=1e-10)
        assert_allclose(result[mask_above], 1.0, atol=1e-10)

    def test_position_dependent_maximum(self) -> None:
        """Max[x - 5, 0] evaluates to ReLU-like array."""
        spec = _make_spec(
            coefficient_symbolic="Max[x() - 5, 0]",
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=0.0,
        )
        grid = CartesianGrid(bounds=[(0, 10), (0, 10)], shape=[20, 20], periodic=True)
        pde = _make_pde(spec, parameters={})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)

        assert isinstance(result, np.ndarray)
        x_coords = np.asarray(grid.cell_coords[..., 0])  # pyright: ignore[reportUnknownArgumentType,reportUnknownVariableType]
        expected = np.maximum(x_coords - 5.0, 0.0)
        assert_allclose(result, expected, atol=1e-10)


# ===========================================================================
# Group 3: Energy computation with background-derived terms
# ===========================================================================


class TestBackgroundEnergy:
    """Verify energy computation with background-field-like coefficient patterns."""

    def test_constant_background_mass_resolves(self) -> None:
        """Constant symbolic mass from background field resolves for energy."""
        from tidal.measurement._energy import _resolve_mass_squared
        from tidal.measurement._io import SimulationData

        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(
                            -4.0,
                            "identity",
                            "phi",
                            coefficient_symbolic="-V0",
                        ),
                    ),
                ),
            ),
            mass_matrix=((4.0,),),
            coupling_matrix=((0.0,),),
            mass_matrix_symbolic=(("-V0",),),
            metadata={},
            coordinates=("t", "x", "y"),
        )

        # Minimal SimulationData — only need spec + parameters for mass resolution
        data = SimulationData(
            spec=spec,
            fields={"phi": np.zeros((2, 8, 8))},
            momenta={"pi_phi": np.zeros((2, 8, 8))},
            times=np.array([0.0, 1.0]),
            grid_spacing=(10.0 / 8, 10.0 / 8),
            grid_bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
            parameters={"V0": 4.0},
        )

        m2 = _resolve_mass_squared(data, 0)
        # Convention: matrix[i][j] = -(coefficient of identity(field_j))
        # mass_matrix_symbolic stores raw "-V0", resolution gives -V0 = -4.0
        # then _resolve_mass_squared negates: -(-4.0) = 4.0
        assert m2 == pytest.approx(4.0)

    def test_position_dependent_mass_in_virial_works(self) -> None:
        """Position-dependent identity term evaluates on grid in virial energy."""
        from tidal.measurement._energy import _compute_virial_potential
        from tidal.measurement._io import SimulationData

        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(
                            -1.0,
                            "identity",
                            "phi",
                            coefficient_symbolic="-V0*UnitStep[x() - 3]",
                            coordinate_dependent=("x",),
                        ),
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )

        data = SimulationData(
            spec=spec,
            fields={"phi": np.ones((2, 8, 8))},
            momenta={"pi_phi": np.zeros((2, 8, 8))},
            times=np.array([0.0, 1.0]),
            grid_spacing=(10.0 / 8, 10.0 / 8),
            grid_bounds=((0.0, 10.0), (0.0, 10.0)),
            periodic=(True, True),
            parameters={"V0": 4.0},
        )

        # Should now work (no longer raises) — position-dependent coefficients
        # are evaluated on the grid for the virial integral.
        result = _compute_virial_potential(data, t_idx=0)
        assert np.isfinite(result)
        # With constant phi=1 field, virial is non-trivial
        assert result != 0.0

    def test_position_dependent_mass_resolved_on_grid(self) -> None:
        """_resolve_mass_squared returns ndarray for position-dependent mass."""
        from tidal.measurement._energy import _resolve_mass_squared
        from tidal.measurement._io import SimulationData

        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(1.0, "laplacian", "phi"),
                        OperatorTerm(
                            -1.0,
                            "identity",
                            "phi",
                            coefficient_symbolic="-V0*exp(-(x()^2 + y()^2))",
                            coordinate_dependent=("x", "y"),
                        ),
                    ),
                ),
            ),
            mass_matrix=((1.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )

        data = SimulationData(
            spec=spec,
            fields={"phi": np.ones((2, 16, 16))},
            momenta={"pi_phi": np.zeros((2, 16, 16))},
            times=np.array([0.0, 1.0]),
            grid_spacing=(10.0 / 16, 10.0 / 16),
            grid_bounds=((-5.0, 5.0), (-5.0, 5.0)),
            periodic=(True, True),
            parameters={"V0": 3.0},
        )

        m2 = _resolve_mass_squared(data, 0)
        # Position-dependent → returns ndarray
        assert isinstance(m2, np.ndarray)
        assert m2.shape == (16, 16)
        # Center cell at (0.3125, 0.3125) — not exactly (0,0) so exp != 1
        # coefficient_symbolic = -V0*exp(-(x^2+y^2)), _resolve_mass_squared negates
        # Verify center has highest value and edges are smaller
        assert m2[8, 8] > m2[0, 0], "Center should be > corner"
        assert m2[8, 8] == pytest.approx(3.0, abs=0.6)  # close to V0=3.0
        assert m2[0, 0] < 0.01  # far from center, exp decays


# ===========================================================================
# Group 4: TOML validation edge cases
# ===========================================================================


class TestBackgroundFieldValidation:
    """TOML validation edge cases for [[background_fields]]."""

    def test_background_gradient_in_lagrangian_raises(
        self, tmp_path: pathlib.Path,
    ) -> None:
        """CD[-a][G[]] in Lagrangian should raise with clear error message."""
        from tidal.cli._derive import _validate_config

        config: dict[str, Any] = {
            "theory": {"name": "Bad Gradient"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "constants": {"names": ["g0"]},
            "background_fields": [
                {"name": "G", "type": "scalar", "components": ["g0"]},
            ],
            "lagrangian": {
                "expression": (
                    "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] "
                    "- CD[-a][G[]] eta[a,b] CD[-b][phi[]]"
                ),
            },
            "output": {"path": "out.json"},
        }
        with pytest.raises(ValueError, match=r"covariant derivative.*background.*G"):
            _validate_config(config)

    def test_background_no_gradient_passes(
        self, tmp_path: pathlib.Path,
    ) -> None:
        """G[] * phi[] (no gradient) should pass validation."""
        from tidal.cli._derive import _validate_config

        config: dict[str, Any] = {
            "theory": {"name": "OK Coupling"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "constants": {"names": ["g0"]},
            "background_fields": [
                {"name": "G", "type": "scalar", "components": ["g0"]},
            ],
            "lagrangian": {
                "expression": "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - G[] * phi[]^2",
            },
            "output": {"path": "out.json"},
        }
        # Should not raise
        _validate_config(config)

    def test_multiple_background_fields(
        self, tmp_path: pathlib.Path,
    ) -> None:
        """Two background fields in same config should work."""
        from tidal.cli import main

        config = tmp_path / "multi_bg.toml"
        config.write_text("""
[theory]
name = "Multi Background"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "phi"
type = "scalar"

[constants]
names = ["V0", "W0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[[background_fields]]
name = "W"
type = "scalar"
components = ["W0"]

[lagrangian]
expression = "-1/2 CD[-a][phi[]] eta[a,b] CD[-b][phi[]] - V[]*phi[]^2 - W[]*phi[]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0

    def test_background_with_derived_fields(
        self, tmp_path: pathlib.Path,
    ) -> None:
        """Background and derived fields should coexist."""
        from tidal.cli import main

        config = tmp_path / "bg_derived.toml"
        config.write_text("""
[theory]
name = "Background With Derived"

[spacetime]
dimension = 2
metric = "minkowski"

[[fields]]
name = "A"
type = "vector"

[constants]
names = ["V0"]

[[background_fields]]
name = "V"
type = "scalar"
components = ["V0"]

[[derived_fields]]
name = "F"
type = "tensor"
rank = 2
symmetry = "antisymmetric"
definition = "CD[-a][A[-b]] - CD[-b][A[-a]]"

[lagrangian]
expression = "-1/4 F[-a, -b] eta[a, c] eta[b, d] F[-c, -d] - V[] * A[-a] eta[a, b] A[-b]"

[output]
path = "output.json"
""")
        ret = main(["derive", str(config), "--dry-run"])
        assert ret == 0


# ---------------------------------------------------------------------------
# Piecewise coefficient resolution (Mathematica Simplify output)
# ---------------------------------------------------------------------------


class TestPiecewiseCoefficientResolution:
    """Test that Mathematica Piecewise expressions resolve correctly."""

    def test_piecewise_position_dependent(self) -> None:
        """Piecewise from scalar_potential_well JSON evaluates on a grid."""
        spec = _make_spec(
            coefficient_symbolic=(
                "Piecewise[{{-V0, Inequality[30, LessEqual, x[], LessEqual, 70]}}, 0]"
            ),
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=1.0,
        )
        grid = CartesianGrid([(0, 100)], 256, periodic=True)
        pde = PDEFromSpec(spec, parameters={"V0": 4.0})
        term = spec.equations[0].rhs_terms[0]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)
        arr = np.asarray(result)
        mid = arr[128]  # x=50.2 — inside [30, 70], coefficient is -V0 = -4.0
        assert mid == pytest.approx(-4.0)
        outside = arr[10]  # x≈3.9 — outside well, coefficient is 0
        assert outside == pytest.approx(0.0)

    def test_piecewise_equivalent_to_unitstep(self) -> None:
        """Piecewise form gives same result as UnitStep product form."""
        pw_expr = (
            "Piecewise[{{-V0, Inequality[30, LessEqual, x[], LessEqual, 70]}}, 0]"
        )
        us_expr = "-V0*UnitStep[x[] - 30]*UnitStep[70 - x[]]"
        grid = CartesianGrid([(0, 100)], 128, periodic=True)

        spec_pw = _make_spec(
            coefficient_symbolic=pw_expr,
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=1.0,
        )
        spec_us = _make_spec(
            coefficient_symbolic=us_expr,
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=1.0,
        )
        pde_pw = PDEFromSpec(spec_pw, parameters={"V0": 4.0})
        pde_us = PDEFromSpec(spec_us, parameters={"V0": 4.0})

        arr_pw = np.asarray(
            pde_pw._resolve_coefficient_at_point(
                spec_pw.equations[0].rhs_terms[0], t=0.0, grid=grid,
            ),
        )
        arr_us = np.asarray(
            pde_us._resolve_coefficient_at_point(
                spec_us.equations[0].rhs_terms[0], t=0.0, grid=grid,
            ),
        )
        # Interior and exterior should match (boundary may differ by 0.5
        # due to heaviside(0)=0.5 vs LessEqual giving full value)
        coords = np.asarray(grid.cell_coords[..., 0], dtype=float)
        interior = (coords > 31) & (coords < 69)
        assert_allclose(arr_pw[interior], arr_us[interior])
        exterior = (coords < 29) | (coords > 71)
        assert_allclose(arr_pw[exterior], arr_us[exterior])


# ---------------------------------------------------------------------------
# Gaussian Exp-based coupling coefficient on 2D grid
# ---------------------------------------------------------------------------


class TestGaussianCoupling2DGrid:
    """Test Gaussian Exp-based coupling coefficient evaluates correctly on a 2D grid."""

    def test_gaussian_coupling_peak_and_decay(self) -> None:
        """Gaussian coupling g0*E^(-r^2/(2R^2)) has peak at center, decays at edges."""
        # Build a 2-field spec with a Gaussian cross-field coupling term
        spec = EquationSystem(
            n_components=2,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi_0", "chi_0"),
            equations=(
                ComponentEquation(
                    field_name="phi_0",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0, "laplacian_x", "phi_0",
                        ),
                        OperatorTerm(
                            1.0, "laplacian_y", "phi_0",
                        ),
                        OperatorTerm(
                            -1.0, "identity", "chi_0",
                            coefficient_symbolic="g0 * E^(-(x[]^2 + y[]^2) / (2 * R^2))",
                            coordinate_dependent=("x", "y"),
                        ),
                    ),
                ),
                ComponentEquation(
                    field_name="chi_0",
                    field_index=1,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0, "laplacian_x", "chi_0",
                        ),
                        OperatorTerm(
                            1.0, "laplacian_y", "chi_0",
                        ),
                        OperatorTerm(
                            -1.0, "identity", "phi_0",
                            coefficient_symbolic="g0 * E^(-(x[]^2 + y[]^2) / (2 * R^2))",
                            coordinate_dependent=("x", "y"),
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 0.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        params = {"g0": 2.0, "R": 8.0}
        pde = PDEFromSpec(spec, parameters=params)
        grid = CartesianGrid([(-30, 30), (-30, 30)], 64, periodic=True)

        # Resolve the cross-field coupling term (index 2 in phi_0 equation)
        term = spec.equations[0].rhs_terms[2]
        result = pde._resolve_coefficient_at_point(term, t=0.0, grid=grid)
        arr = np.asarray(result)

        # No NaN or Inf
        assert np.all(np.isfinite(arr)), "Gaussian coupling has NaN/Inf"

        # Peak at center: _resolve_coefficient_at_point evaluates the symbolic
        # expression directly (g0 * exp(0) = 2.0); numeric coefficient is NOT applied.
        center = arr[32, 32]
        assert center == pytest.approx(2.0, rel=0.1)

        # Edge should be near zero (large r^2)
        edge = arr[0, 0]
        assert abs(edge) < 0.1, f"Edge value should decay to ~0, got {edge}"

        # Symmetric: arr[i,j] ≈ arr[j,i] (rotational symmetry)
        assert_allclose(arr[20, 32], arr[32, 20], atol=1e-10)


# ===========================================================================
# Group 6: Spatial coefficient caching (B6)
# ===========================================================================


class TestSpatialCoefficientCache:
    """Verify that spatial-only coefficients are cached after first use."""

    def test_spatial_cache_populated_on_first_rhs(self) -> None:
        """Spatial-only coefficient should be cached after _ensure_spatial_cache."""
        spec = _make_spec(
            coefficient_symbolic="g0 * exp(-(x()^2 + y()^2) / (2 * R^2))",
            coordinate_dependent=("x", "y"),
            operator="identity",
            coefficient=1.0,
        )
        pde = PDEFromSpec(spec, parameters={"g0": 2.0, "R": 8.0})
        grid = CartesianGrid([(-30, 30), (-30, 30)], 32, periodic=True)

        # Before first use, spatial cache should be empty
        assert not pde._spatial_cache_ready
        assert len(pde._spatial_cache) == 0

        # Trigger cache population
        pde._ensure_spatial_cache(grid)

        assert pde._spatial_cache_ready
        assert len(pde._spatial_cache) == 1
        cached_key = "g0 * exp(-(x()^2 + y()^2) / (2 * R^2))"
        assert cached_key in pde._spatial_cache
        arr = pde._spatial_cache[cached_key]
        assert arr.shape == (32, 32)
        # Peak at center
        assert arr[16, 16] == pytest.approx(2.0, rel=0.1)

    def test_spatial_cache_skips_time_dependent(self) -> None:
        """Time-dependent terms should NOT be cached in spatial cache."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0, "identity", "phi",
                            coefficient_symbolic="g0 * exp(-t()) * exp(-x()^2)",
                            coordinate_dependent=("x",),
                            time_dependent=True,
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={"g0": 1.0})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        pde._ensure_spatial_cache(grid)

        # Time-dependent term should NOT be in spatial cache; no spatial-only
        # terms exist so cache stays in initial state
        assert not pde._spatial_cache_ready
        assert len(pde._spatial_cache) == 0

    def test_time_dependent_background_evaluates_correctly(self) -> None:
        """Time-dependent background coefficient should vary with t."""
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0, "identity", "phi",
                            coefficient_symbolic="g0 * cos(t())",
                            time_dependent=True,
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={"g0": 2.0})
        term = spec.equations[0].rhs_terms[0]

        # t=0: cos(0) = 1, so coefficient = 2.0
        assert pde._resolve_coefficient_at_point(term, t=0.0, grid=None) == pytest.approx(2.0)
        # t=pi/2: cos(pi/2) = 0, so coefficient = 0.0
        assert pde._resolve_coefficient_at_point(term, t=np.pi / 2, grid=None) == pytest.approx(0.0, abs=1e-14)
        # t=pi: cos(pi) = -1, so coefficient = -2.0
        assert pde._resolve_coefficient_at_point(term, t=np.pi, grid=None) == pytest.approx(-2.0)

    def test_time_and_position_dependent_background(self) -> None:
        """Background with both time and position dependence evaluates correctly."""
        spec = _make_spec(
            coefficient_symbolic="g0 * exp(-x()^2) * cos(t())",
            coordinate_dependent=("x",),
        )
        # Need to mark as time-dependent too — rebuild the term
        spec = EquationSystem(
            n_components=1,
            dimension=3,
            spatial_dimension=2,
            component_names=("phi",),
            equations=(
                ComponentEquation(
                    field_name="phi",
                    field_index=0,
                    time_derivative_order=2,
                    rhs_terms=(
                        OperatorTerm(
                            1.0, "laplacian", "phi",
                            coefficient_symbolic="g0 * exp(-x()^2) * cos(t())",
                            coordinate_dependent=("x",),
                            time_dependent=True,
                        ),
                    ),
                ),
            ),
            mass_matrix=((0.0,),),
            coupling_matrix=((0.0,),),
            metadata={},
            coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={"g0": 3.0})
        grid = CartesianGrid([(-5, 5), (-5, 5)], 16, periodic=True)
        term = spec.equations[0].rhs_terms[0]

        # t=0: cos(0)=1, so result = 3*exp(-x^2)
        result_t0 = np.asarray(
            pde._resolve_coefficient_at_point(term, t=0.0, grid=grid),
        )
        assert result_t0[8, 8] == pytest.approx(3.0, rel=0.1)  # center: x≈0

        # t=pi: cos(pi)=-1, so result = -3*exp(-x^2)
        result_tpi = np.asarray(
            pde._resolve_coefficient_at_point(term, t=np.pi, grid=grid),
        )
        assert_allclose(result_tpi, -result_t0, atol=1e-12)

    def test_sign_change_warning_for_position_dependent_mass(self) -> None:
        """Position-dependent mass that changes sign emits UserWarning."""
        import warnings

        # Build spec with a mass coefficient that changes sign: sin(x)
        spec = _make_spec(
            coefficient_symbolic="sin(x())",
            coordinate_dependent=("x",),
            operator="identity",
            field="phi",
            coefficient=1.0,
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 32, periodic=True)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pde._ensure_spatial_cache(grid)

        tachyonic_warnings = [
            x for x in w
            if "tachyonic" in str(x.message).lower()
        ]
        assert len(tachyonic_warnings) == 1
        assert "changes sign" in str(tachyonic_warnings[0].message)

    def test_no_warning_for_positive_definite_mass(self) -> None:
        """Position-dependent mass that stays positive should not warn."""
        import warnings

        # exp(-x^2) is always positive
        spec = _make_spec(
            coefficient_symbolic="exp(-x()^2)",
            coordinate_dependent=("x",),
            operator="identity",
            field="phi",
            coefficient=1.0,
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 32, periodic=True)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            pde._ensure_spatial_cache(grid)

        tachyonic_warnings = [
            x for x in w
            if "tachyonic" in str(x.message).lower()
        ]
        assert len(tachyonic_warnings) == 0

    def test_spatial_cache_reused_across_calls(self) -> None:
        """Ensure cached value is the same object across multiple cache checks."""
        spec = _make_spec(
            coefficient_symbolic="g0 * exp(-x()^2)",
            coordinate_dependent=("x",),
            operator="identity",
            coefficient=1.0,
        )
        pde = PDEFromSpec(spec, parameters={"g0": 3.0})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        pde._ensure_spatial_cache(grid)
        cached_arr = pde._spatial_cache["g0 * exp(-x()^2)"]

        # Second call should be a no-op
        pde._ensure_spatial_cache(grid)
        assert pde._spatial_cache["g0 * exp(-x()^2)"] is cached_arr


# ===========================================================================
# Group 8: Measurement module hardening — FFT guards for position-dependent
# ===========================================================================


def _make_position_dependent_sim_data(
    *,
    position_dependent_coupling: bool = False,
    position_dependent_mass: bool = False,
) -> SimulationData:
    """Build synthetic SimulationData with position-dependent terms.

    Two 2+1D scalar fields phi_0, chi_0 with optional position-dependent
    mass or coupling terms.
    """
    phi_terms_list: list[OperatorTerm] = [
        OperatorTerm(coefficient=1.0, operator="laplacian", field="phi_0"),
        OperatorTerm(coefficient=-1.0, operator="identity", field="phi_0"),
    ]
    chi_terms_list: list[OperatorTerm] = [
        OperatorTerm(coefficient=1.0, operator="laplacian", field="chi_0"),
        OperatorTerm(coefficient=-1.0, operator="identity", field="chi_0"),
    ]

    if position_dependent_mass:
        # Override phi's self-mass to be position-dependent
        phi_terms_list[1] = OperatorTerm(
            coefficient=-1.0,
            operator="identity",
            field="phi_0",
            coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
            coordinate_dependent=("x", "y"),
        )

    if position_dependent_coupling:
        # Add position-dependent cross-field coupling
        phi_terms_list.append(OperatorTerm(
            coefficient=-0.1,
            operator="identity",
            field="chi_0",
            coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
            coordinate_dependent=("x", "y"),
        ))
        chi_terms_list.append(OperatorTerm(
            coefficient=-0.1,
            operator="identity",
            field="phi_0",
            coefficient_symbolic="-g0 * exp(-(x()^2 + y()^2))",
            coordinate_dependent=("x", "y"),
        ))

    eq_phi = ComponentEquation(
        field_name="phi_0", field_index=0, time_derivative_order=2,
        rhs_terms=tuple(phi_terms_list),
    )
    eq_chi = ComponentEquation(
        field_name="chi_0", field_index=1, time_derivative_order=2,
        rhs_terms=tuple(chi_terms_list),
    )

    spec = EquationSystem(
        n_components=2,
        dimension=3,
        spatial_dimension=2,
        component_names=("phi_0", "chi_0"),
        equations=(eq_phi, eq_chi),
        mass_matrix=((1.0, 0.0), (0.0, 1.0)),
        coupling_matrix=((0.0, 0.1), (0.1, 0.0)) if position_dependent_coupling else ((0.0, 0.0), (0.0, 0.0)),
        metadata={},
        coordinates=("t", "x", "y"),
    )

    n_grid = 16
    n_snap = 5
    dx = 10.0 / n_grid
    shape = (n_snap, n_grid, n_grid)
    rng = np.random.default_rng(42)

    return SimulationData(
        times=np.linspace(0.0, 1.0, n_snap),
        fields={
            "phi_0": rng.standard_normal(shape),
            "chi_0": rng.standard_normal(shape) * 0.01,
        },
        momenta={
            "phi_0": rng.standard_normal(shape),
            "chi_0": rng.standard_normal(shape) * 0.01,
        },
        grid_spacing=(dx, dx),
        grid_bounds=((-5.0, 5.0), (-5.0, 5.0)),
        periodic=(True, True),
        spec=spec,
        parameters={"g0": 1.0},
    )


class TestMeasurementFFTGuards:
    """Verify FFT-based methods reject position-dependent systems."""

    def test_spectral_energy_rejects_position_dependent_mass(self) -> None:
        """compute_spectral_energy raises TypeError for ndarray mass."""
        from tidal.measurement._spectral import compute_spectral_energy

        field = np.ones((16, 16))
        mom = np.zeros((16, 16))
        mass_array = np.ones((16, 16)) * 2.0  # ndarray, not scalar

        with pytest.raises(TypeError, match="spatially uniform mass"):
            compute_spectral_energy(field, mom, mass_array, (0.5, 0.5), (True, True))

    def test_spectral_conversion_rejects_position_dependent_mass(self) -> None:
        """Spectral conversion P(k,t) rejects position-dependent mass."""
        from tidal.measurement._spectral_conversion import compute_spectral_conversion

        data = _make_position_dependent_sim_data(position_dependent_mass=True)
        with pytest.raises(ValueError, match="position-dependent"):
            compute_spectral_conversion(data, "phi_0", "chi_0")

    def test_spectral_conversion_rejects_position_dependent_coupling(self) -> None:
        """Spectral conversion P(k,t) rejects position-dependent coupling."""
        from tidal.measurement._spectral_conversion import compute_spectral_conversion

        data = _make_position_dependent_sim_data(position_dependent_coupling=True)
        with pytest.raises(ValueError, match="position-dependent"):
            compute_spectral_conversion(data, "phi_0", "chi_0")

    def test_group_spectral_conversion_rejects_position_dependent(self) -> None:
        """Group spectral conversion also rejects position-dependent terms."""
        from tidal.measurement._spectral_conversion import (
            compute_group_spectral_conversion,
        )

        data = _make_position_dependent_sim_data(position_dependent_coupling=True)
        with pytest.raises(ValueError, match="position-dependent"):
            compute_group_spectral_conversion(data, "phi_0", "chi_0")

    def test_dispersion_rejects_position_dependent_coupling(self) -> None:
        """Dispersion relation omega(k) rejects position-dependent coupling."""
        from tidal.measurement._dispersion import compute_dispersion

        data = _make_position_dependent_sim_data(position_dependent_coupling=True)
        with pytest.raises(ValueError, match="position-dependent"):
            compute_dispersion(data, "phi_0")

    def test_conversion_probability_accepts_position_dependent_mass(self) -> None:
        """P(t) works with position-dependent mass (not FFT-based)."""
        from tidal.measurement._conversion import compute_conversion_probability

        data = _make_position_dependent_sim_data(position_dependent_mass=True)
        # Should NOT raise — P(t) uses real-space energy, not Fourier
        result = compute_conversion_probability(data, "phi_0", "chi_0")
        assert result.probability.shape == (5,)

    def test_energy_conservation_with_position_dependent_mass(self) -> None:
        """check_energy_conservation works with position-dependent mass."""
        from tidal.measurement._diagnostics import check_energy_conservation

        data = _make_position_dependent_sim_data(position_dependent_mass=True)
        result = check_energy_conservation(data)
        assert result.relative_error.shape == (5,)


# ===========================================================================
# Group 9: Proca + Scalar Background (Phase 2B)
# ===========================================================================


class TestProcaScalarBackground:
    """Tests for coupled Proca with Lorentzian scalar background."""

    TOML_PATH = (
        Path(__file__).parent.parent
        / "examples"
        / "proca_background"
        / "theory.toml"
    )

    def test_proca_background_wls_dry_run(self) -> None:
        """WLS generation succeeds and contains ReplaceAll for G."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)
        # Scalar BG uses ReplaceAll (not ComponentValue)
        assert "/." in wls or "ReplaceAll" in wls
        # Must contain the Lorentzian profile
        assert "Power" in wls
        # Must define background tensor G
        assert "pwlG[]" in wls or "DefTensor[pwlG[]" in wls

    def test_proca_background_lorentzian_coefficient_evaluates(self) -> None:
        """Lorentzian profile evaluates correctly: g0/(1 + r^2/R^2)."""
        from tidal.symbolic._eval_utils import evaluate_coefficient

        expr = "g0 * (1 + (x()^2 + y()^2) / R^2)**(-1)"
        params = {"g0": 1.0, "R": 8.0}
        x = np.linspace(-10, 10, 32)
        y = np.linspace(-10, 10, 32)
        xg, yg = np.meshgrid(x, y, indexing="ij")

        result = evaluate_coefficient(
            expr, params, ("t", "x", "y"),
            coord_arrays={"x": xg, "y": yg},
        )
        assert isinstance(result, np.ndarray)
        # Origin should be g0 = 1.0
        center = result[16, 16]
        assert_allclose(center, 1.0, atol=0.1)
        # Corners should be smaller (Lorentzian decay)
        corner = result[0, 0]
        assert corner < center
        # Algebraic tails — at r^2 ~ 200, G ~ 1/(1+200/64) ~ 0.24
        assert corner > 0.1  # heavier tails than Gaussian

    def test_proca_background_constraint_source_position_dependent(self) -> None:
        """Constraint source handles position-dependent cross-field coupling.

        Manually construct a spec where A_0 is a constraint with a
        position-dependent coupling to B_1, simulating what the Wolfram
        pipeline would produce for the Lorentzian background.
        """
        # A_0: constraint (time_derivative_order=0), coupled to B_1 via G(x,y)
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
            OperatorTerm(
                coefficient=-1.0, operator="identity", field="B_1",
                coefficient_symbolic="-gcoup * g0 * (1 + (x()^2 + y()^2) / R^2)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
        )
        # B_1: dynamical with constant mass
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        eq_a0 = ComponentEquation(
            field_name="A_0", field_index=0, time_derivative_order=0,
            rhs_terms=a0_terms,
        )
        eq_b1 = ComponentEquation(
            field_name="B_1", field_index=1, time_derivative_order=2,
            rhs_terms=b1_terms,
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(eq_a0, eq_b1),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={"gcoup": 0.5, "g0": 1.0, "R": 8.0})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        # The constraint source term should be resolved as an array
        from pde import FieldCollection, ScalarField
        b1_field = ScalarField(grid, data=np.ones((16, 16)))
        a0_field = ScalarField(grid, data=0.0)
        state = FieldCollection([a0_field, b1_field])

        source = pde._compute_constraint_source(
            0, state, "periodic", t=0.0,
            source_terms=[t for t in spec.equations[0].rhs_terms
                          if t.field != "A_0"],
        )
        # Source should be non-uniform because of position-dependent coupling
        source_data = source.data
        assert source_data.shape == (16, 16)
        # Center should have stronger coupling than corners
        assert abs(source_data[8, 8]) > abs(source_data[0, 0])

    def test_proca_background_spatial_cache_with_constraints(self) -> None:
        """Spatial cache populates for position-dependent constraint terms."""
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="B_1",
                coefficient_symbolic="-gcoup * g0 * (1 + (x()^2 + y()^2) / R^2)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
        )
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(
                ComponentEquation(field_name="A_0", field_index=0, time_derivative_order=0, rhs_terms=a0_terms),
                ComponentEquation(field_name="B_1", field_index=1, time_derivative_order=2, rhs_terms=b1_terms),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={"gcoup": 0.5, "g0": 1.0, "R": 8.0})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        pde._ensure_spatial_cache(grid)
        # The Lorentzian expression should be cached
        assert len(pde._spatial_cache) > 0

    def test_proca_background_energy_with_position_dependent_coupling(self) -> None:
        """compute_system_energy works for constraint + dynamic with pos-dep coupling."""
        from tidal.measurement._energy import compute_system_energy

        # Build 2-field spec: A_0 (constraint) + B_1 (dynamic) with G(x,y) coupling
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="B_1",
                coefficient_symbolic="-0.5 * (1 + (x()^2 + y()^2) / 64)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
        )
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(
                ComponentEquation(field_name="A_0", field_index=0, time_derivative_order=0, rhs_terms=a0_terms),
                ComponentEquation(field_name="B_1", field_index=1, time_derivative_order=2, rhs_terms=b1_terms),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        data = SimulationData(
            times=np.array([0.0, 1.0]),
            fields={
                "A_0": np.zeros((2, 16, 16)),
                "B_1": np.ones((2, 16, 16)) * 0.1,
            },
            momenta={"B_1": np.zeros((2, 16, 16))},
            grid_spacing=(1.25, 1.25),
            grid_bounds=((-10.0, 10.0), (-10.0, 10.0)),
            periodic=(True, True),
            spec=spec,
            parameters={},
        )
        result = compute_system_energy(data, 0)
        assert result.total >= 0.0

    def test_proca_background_measurement_cross_validation(self) -> None:
        """Manual Lorentzian energy vs virial for simple constraint+dynamic system."""
        from tidal.measurement._energy import compute_system_energy

        # B_1 dynamic field with constant mass; A_0 constraint (no energy contribution)
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
        )
        b1_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="B_1"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_1"),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("A_0", "B_1"),
            equations=(
                ComponentEquation(field_name="A_0", field_index=0, time_derivative_order=0, rhs_terms=a0_terms),
                ComponentEquation(field_name="B_1", field_index=1, time_derivative_order=2, rhs_terms=b1_terms),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        # Uniform B_1 = 1.0: KE=0, PE = m²/2 * integral of B_1^2
        n_grid = 16
        dx = 20.0 / n_grid
        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "A_0": np.zeros((1, n_grid, n_grid)),
                "B_1": np.ones((1, n_grid, n_grid)),
            },
            momenta={"B_1": np.zeros((1, n_grid, n_grid))},
            grid_spacing=(dx, dx),
            grid_bounds=((-10.0, 10.0), (-10.0, 10.0)),
            periodic=(True, True),
            spec=spec,
            parameters={},
        )
        result = compute_system_energy(data, 0)
        # Expected: m²=2, B_1=1 everywhere, V = 10.0 per grid cell
        # E_potential = 0.5 * m² * integral(B_1²) dV = 0.5 * 2 * 1² * 20² = 400
        expected_potential = 0.5 * 2.0 * (20.0 ** 2)
        assert_allclose(result.total, expected_potential, rtol=0.1)


# ===========================================================================
# Group 10: Vector Background (Phase 2C)
# ===========================================================================


class TestVectorBackground:
    """Tests for scalar-vector coupling with tanh domain wall vector background."""

    TOML_PATH = (
        Path(__file__).parent.parent
        / "examples"
        / "vector_background"
        / "theory.toml"
    )

    def test_vector_background_wls_dry_run(self) -> None:
        """WLS generation succeeds and contains ComponentValue lines."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)
        # Must contain Tanh profile
        assert "Tanh" in wls
        # Must define background tensor B
        assert "vbdB" in wls

    def test_vector_background_component_values_emitted(self) -> None:
        """WLS contains 3 ComponentValue lines for the vector background."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)

        # Count ComponentValue lines for the background B
        cv_lines = [
            line for line in wls.splitlines()
            if "ComponentValue[vbdB[" in line
        ]
        assert len(cv_lines) == 3, f"Expected 3 ComponentValue lines, got {len(cv_lines)}: {cv_lines}"
        # B_2 has tanh profile
        assert any("Tanh" in line for line in cv_lines)

    def test_vector_background_validation_check_emitted(self) -> None:
        """WLS contains ValidateNoUnresolvedBackgrounds call."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)

        assert "ValidateNoUnresolvedBackgrounds" in wls
        # Specifically validates the vector background B (not scalar backgrounds)
        assert "vbdB" in wls

    def test_vector_background_no_replace_all_for_vector(self) -> None:
        """Vector BG does NOT get scalar-style ReplaceAll substitution."""
        import tomllib

        from tidal.cli._derive import generate_wls

        with self.TOML_PATH.open("rb") as f:
            config = tomllib.load(f)
        wls = generate_wls(config, config_dir=self.TOML_PATH.parent)

        # Scalar BGs use /. (ReplaceAll) — vector BGs should NOT
        replace_all_lines = [
            line for line in wls.splitlines()
            if "/." in line and "vbdB" in line
        ]
        assert len(replace_all_lines) == 0, (
            f"Vector background B should NOT use ReplaceAll, but found: {replace_all_lines}"
        )

    def test_vector_background_tanh_coefficient_evaluates(self) -> None:
        """Tanh * Exp profile evaluates correctly on a grid."""
        from tidal.symbolic._eval_utils import evaluate_coefficient

        expr = "B0 * tanh(x / W) * exp(-(y**2) / (2 * R**2))"
        params: dict[str, float] = {"B0": 1.0, "W": 3.0, "R": 8.0}
        # Use 33 points so x=0.0 falls exactly on index 16
        n = 33
        x = np.linspace(-10, 10, n)
        y = np.linspace(-10, 10, n)
        xg, yg = np.meshgrid(x, y, indexing="ij")

        result = evaluate_coefficient(
            expr, params, ("t", "x", "y"),
            coord_arrays={"x": xg, "y": yg},
        )
        assert isinstance(result, np.ndarray)
        assert result.shape == (n, n)

        # Antisymmetric in x: B(x, y) = -B(-x, y) at y=0
        mid = n // 2  # index 16 = x=0
        assert_allclose(result[:, mid] + result[::-1, mid], 0.0, atol=1e-10)

        # At origin: tanh(0) = 0, exp(0) = 1 → B = 0
        assert_allclose(result[mid, mid], 0.0, atol=1e-10)

        # At x=+10, y=0: tanh(10/3) ~ 0.9997, exp(0) = 1
        assert result[-1, mid] > 0.9

        # At x=-10, y=0: tanh(-10/3) ~ -0.9997
        assert result[0, mid] < -0.9

    def test_vector_background_zero_component_constant(self) -> None:
        """B_0=0, B_1=0 components give constant zero, not position-dependent."""
        from tidal.symbolic._eval_utils import evaluate_coefficient

        # The zero components should evaluate to 0.0 (scalar), not an array
        result = evaluate_coefficient("0", {}, ("t", "x", "y"))
        assert result == 0.0

    def test_vector_background_wrong_component_count_raises(self) -> None:
        """2 components for a 3D vector raises ValueError."""
        from tidal.cli._derive import generate_wls

        config: dict[str, Any] = {
            "theory": {"name": "Test"},
            "spacetime": {"dimension": 3, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "constants": {"names": ["B0"]},
            "background_fields": [{"name": "B", "type": "vector", "components": ["0", "B0"]}],
            "lagrangian": {"expression": "CD[-a][phi[]] eta[a,b] CD[-b][phi[]]"},
            "output": {"path": "out.json"},
        }
        with pytest.raises(ValueError, match="components"):
            generate_wls(config)

    def test_vector_background_energy_computation(self) -> None:
        """compute_system_energy works with mixed-rank position-dependent spec."""
        from tidal.measurement._energy import compute_system_energy

        # phi: scalar dynamic, A_2: dynamic with position-dependent coupling to phi
        phi_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="phi_0"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="phi_0"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="A_2",
                coefficient_symbolic="-gBV * B0 * tanh(x / W) * exp(-(y**2) / (2 * R**2))",
                coordinate_dependent=("x", "y"),
            ),
        )
        a2_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="A_2"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="A_2"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="phi_0",
                coefficient_symbolic="-gBV * B0 * tanh(x / W) * exp(-(y**2) / (2 * R**2))",
                coordinate_dependent=("x", "y"),
            ),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("phi_0", "A_2"),
            equations=(
                ComponentEquation(field_name="phi_0", field_index=0, time_derivative_order=2, rhs_terms=phi_terms),
                ComponentEquation(field_name="A_2", field_index=1, time_derivative_order=2, rhs_terms=a2_terms),
            ),
            mass_matrix=((1.0, 0.0), (0.0, 2.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        n_grid = 16
        dx = 20.0 / n_grid
        data = SimulationData(
            times=np.array([0.0]),
            fields={
                "phi_0": np.ones((1, n_grid, n_grid)),
                "A_2": np.zeros((1, n_grid, n_grid)),
            },
            momenta={
                "phi_0": np.zeros((1, n_grid, n_grid)),
                "A_2": np.zeros((1, n_grid, n_grid)),
            },
            grid_spacing=(dx, dx),
            grid_bounds=((-10.0, 10.0), (-10.0, 10.0)),
            periodic=(True, True),
            spec=spec,
            parameters={"gBV": 0.5, "B0": 1.0, "W": 3.0, "R": 8.0},
        )
        result = compute_system_energy(data, 0)
        # Should succeed — position-dependent coupling evaluated on grid
        assert result.total >= 0.0
        assert "phi_0" in result.per_field


class TestConstraintAutoDetectionPositionDependent:
    """Tests that method='auto' correctly routes to matrix/Gauss-Seidel
    when constraint terms have position-dependent coefficients.
    """

    def test_coupled_constraints_auto_routes_to_gauss_seidel(self) -> None:
        """Coupled constraints with position-dependent cross-terms use Gauss-Seidel.

        Simulates the proca_background scenario: A_0 and B_0 are coupled
        constraints with position-dependent identity terms from a Lorentzian
        background. method='auto' on a periodic grid should route to
        Gauss-Seidel instead of FFT.
        """
        from pde import FieldCollection, ScalarField

        lorentzian_expr = "(1 + (x()**2 + y()**2) / 64)**(-1)"
        a0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="A_0"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="A_0"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="B_0",
                coefficient_symbolic=f"-0.5 * {lorentzian_expr}",
                coordinate_dependent=("x", "y"),
            ),
            # Source from dynamical field rho
            OperatorTerm(coefficient=1.0, operator="identity", field="rho"),
        )
        b0_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="B_0"),
            OperatorTerm(coefficient=-2.0, operator="identity", field="B_0"),
            OperatorTerm(
                coefficient=-0.5, operator="identity", field="A_0",
                coefficient_symbolic=f"-0.5 * {lorentzian_expr}",
                coordinate_dependent=("x", "y"),
            ),
        )
        # Source field to drive constraints
        rho_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="rho"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="rho"),
        )
        spec = EquationSystem(
            n_components=3, dimension=3, spatial_dimension=2,
            component_names=("A_0", "B_0", "rho"),
            equations=(
                ComponentEquation(
                    field_name="A_0", field_index=0, time_derivative_order=0,
                    rhs_terms=a0_terms,
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True, method="auto", max_iterations=30,
                    ),
                ),
                ComponentEquation(
                    field_name="B_0", field_index=1, time_derivative_order=0,
                    rhs_terms=b0_terms,
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True, method="auto", max_iterations=30,
                    ),
                ),
                ComponentEquation(
                    field_name="rho", field_index=2, time_derivative_order=2,
                    rhs_terms=rho_terms,
                ),
            ),
            mass_matrix=((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 1.0)),
            coupling_matrix=((0.0, 0.0, -1.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        # Gaussian source in rho
        x, y = grid.cell_coords[..., 0], grid.cell_coords[..., 1]
        gaussian = np.exp(-(x**2 + y**2) / 8)
        state = FieldCollection([
            ScalarField(grid, data=0.0),  # A_0
            ScalarField(grid, data=0.0),  # B_0
            ScalarField(grid, data=gaussian),  # rho
            ScalarField(grid, data=0.0),  # pi_rho
        ])

        # Should NOT raise ValueError about FFT incompatibility
        pde.evolution_rate(state, t=0.0)

        # A_0 should be non-zero (driven by rho source)
        a0_data = np.asarray(state[0].data)
        assert np.any(np.abs(a0_data) > 1e-12)

    def test_single_constraint_auto_routes_to_matrix(self) -> None:
        """Single constraint with position-dependent self-term uses matrix solver.

        A single constraint with a spatially-varying identity coefficient
        should auto-detect and use the matrix path, not FFT.
        """
        from pde import FieldCollection, ScalarField

        phi_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="phi"),
            OperatorTerm(
                coefficient=-1.0, operator="identity", field="phi",
                coefficient_symbolic="-(1 + (x()**2 + y()**2) / 64)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
            # Source from dynamical field rho
            OperatorTerm(coefficient=1.0, operator="identity", field="rho"),
        )
        rho_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="rho"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="rho"),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("phi", "rho"),
            equations=(
                ComponentEquation(
                    field_name="phi", field_index=0, time_derivative_order=0,
                    rhs_terms=phi_terms,
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True, method="auto",
                    ),
                ),
                ComponentEquation(
                    field_name="rho", field_index=1, time_derivative_order=2,
                    rhs_terms=rho_terms,
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 1.0)),
            coupling_matrix=((0.0, -1.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        x, y = grid.cell_coords[..., 0], grid.cell_coords[..., 1]
        gaussian = np.exp(-(x**2 + y**2) / 8)
        state = FieldCollection([
            ScalarField(grid, data=0.0),  # phi (constraint)
            ScalarField(grid, data=gaussian),  # rho
            ScalarField(grid, data=0.0),  # pi_rho
        ])

        # Should NOT raise ValueError about FFT incompatibility
        pde.evolution_rate(state, t=0.0)

        # Constraint should be solved
        phi_data = np.asarray(state[0].data)
        assert np.any(np.abs(phi_data) > 1e-12)

    def test_explicit_fft_with_position_dependent_still_raises(self) -> None:
        """Explicit method='fft' with position-dependent self-term still fails.

        The auto-detection gracefully routes to matrix, but if the user
        explicitly requests FFT, they should get a clear error message.
        """
        from pde import FieldCollection, ScalarField

        phi_terms = (
            OperatorTerm(coefficient=-1.0, operator="laplacian", field="phi"),
            OperatorTerm(
                coefficient=-1.0, operator="identity", field="phi",
                coefficient_symbolic="-(1 + (x()**2 + y()**2) / 64)**(-1)",
                coordinate_dependent=("x", "y"),
            ),
        )
        rho_terms = (
            OperatorTerm(coefficient=1.0, operator="laplacian", field="rho"),
            OperatorTerm(coefficient=-1.0, operator="identity", field="rho"),
        )
        spec = EquationSystem(
            n_components=2, dimension=3, spatial_dimension=2,
            component_names=("phi", "rho"),
            equations=(
                ComponentEquation(
                    field_name="phi", field_index=0, time_derivative_order=0,
                    rhs_terms=phi_terms,
                    constraint_solver=ConstraintSolverConfig(
                        enabled=True, method="fft",
                    ),
                ),
                ComponentEquation(
                    field_name="rho", field_index=1, time_derivative_order=2,
                    rhs_terms=rho_terms,
                ),
            ),
            mass_matrix=((0.0, 0.0), (0.0, 1.0)),
            coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
            metadata={}, coordinates=("t", "x", "y"),
        )
        pde = PDEFromSpec(spec, parameters={})
        grid = CartesianGrid([(-10, 10), (-10, 10)], 16, periodic=True)

        x, y = grid.cell_coords[..., 0], grid.cell_coords[..., 1]
        gaussian = np.exp(-(x**2 + y**2) / 8)
        state = FieldCollection([
            ScalarField(grid, data=0.0),
            ScalarField(grid, data=gaussian),
            ScalarField(grid, data=0.0),
        ])

        with pytest.raises(ValueError, match="FFT"):
            pde.evolution_rate(state, t=0.0)
