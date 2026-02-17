"""Stress tests for background fields feature.

Tests the full pipeline from _mathematica_to_python conversion through
coefficient resolution to energy computation.  These tests exercise the
code paths that background field coefficients trigger — including
position-dependent expressions like UnitStep, Sign, Max, Min.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid

from tidal.symbolic.json_loader import (
    ComponentEquation,
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

    def test_position_dependent_mass_in_virial_raises(self) -> None:
        """Position-dependent identity term raises ValueError in virial energy."""
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

        with pytest.raises(ValueError, match=r"position-dependent|Position-dependent"):
            _compute_virial_potential(data, t_idx=0)


# ===========================================================================
# Group 4: TOML validation edge cases
# ===========================================================================


class TestBackgroundFieldValidation:
    """TOML validation edge cases for [[background_fields]]."""

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
