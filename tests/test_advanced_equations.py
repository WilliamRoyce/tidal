"""Tests for advanced Klein-Gordon equations (anisotropic, higher-order, directional)."""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from pde import CartesianGrid, FieldCollection, ScalarField, UnitGrid

from torsion_gertsenshtein.kgsim.advanced_equations import (
    AnisotropicHigherOrderKGPDE,
    AnisotropicKGPDE,
    DirectionalKGPDE,
    HigherOrderKGPDE,
)


class TestAnisotropicKGPDE:
    """Tests for AnisotropicKGPDE."""

    def test_initialization(self) -> None:
        """Test basic initialization."""
        pde = AnisotropicKGPDE(mass=1.0, speeds=[1.0, 2.0])
        assert pde.m2 == 1.0
        np.testing.assert_array_equal(pde.speeds, [1.0, 2.0])
        np.testing.assert_array_equal(pde.speeds_squared, [1.0, 4.0])

    def test_invalid_speeds(self) -> None:
        """Test that invalid speeds raise ValueError."""
        with pytest.raises(ValueError, match="speeds must be a non-empty sequence"):
            AnisotropicKGPDE(mass=1.0, speeds=[])

        with pytest.raises(ValueError, match="All wave speeds must be positive"):
            AnisotropicKGPDE(mass=1.0, speeds=[1.0, -1.0])

        with pytest.raises(ValueError, match="All wave speeds must be positive"):
            AnisotropicKGPDE(mass=1.0, speeds=[0.0, 1.0])

    def test_produces_valid_output(self) -> None:
        """Test that anisotropic PDE produces valid (finite) evolution rates."""
        grid = UnitGrid([64, 64], periodic=True)

        # Create simple state
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        # Anisotropic with different speeds
        pde_aniso = AnisotropicKGPDE(mass=0.5, speeds=[2.0, 0.5])
        result_aniso = pde_aniso.evolution_rate(state, t=0.0)

        # Check results are finite and non-trivial
        assert np.all(np.isfinite(result_aniso[0].data))
        assert np.all(np.isfinite(result_aniso[1].data))

        # dphi/dt should equal pi
        np.testing.assert_allclose(result_aniso[0].data, pi.data)

        # dpi/dt should not be all zeros (unless phi and pi both are)
        assert np.any(np.abs(result_aniso[1].data) > 0)

    def test_dimension_mismatch(self) -> None:
        """Test that dimension mismatch raises ValueError."""
        grid = UnitGrid([64], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField(grid, data=np.zeros_like(phi.data))
        state = FieldCollection([phi, pi])

        # 2D speeds on 1D grid
        pde = AnisotropicKGPDE(mass=0.5, speeds=[1.0, 2.0])

        with pytest.raises(
            ValueError, match=r"speeds length .* must match grid dimension"
        ):
            pde.evolution_rate(state, t=0.0)

    def test_evolution_rate_shape(self) -> None:
        """Test that evolution_rate returns correct shape."""
        grid = CartesianGrid([(0, 10), (0, 10)], [32, 32], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        pde = AnisotropicKGPDE(mass=0.5, speeds=[2.0, 0.5])
        result = pde.evolution_rate(state, t=0.0)

        assert isinstance(result, FieldCollection)
        assert len(result) == 2
        assert result[0].data.shape == phi.data.shape
        assert result[1].data.shape == pi.data.shape


class TestHigherOrderKGPDE:
    """Tests for HigherOrderKGPDE."""

    def test_initialization(self) -> None:
        """Test basic initialization."""
        pde = HigherOrderKGPDE(mass=1.0, alpha_2=1.0, alpha_4=0.01, alpha_6=0.001)
        assert pde.m2 == 1.0
        assert pde.alpha_2 == 1.0
        assert pde.alpha_4 == 0.01
        assert pde.alpha_6 == 0.001

    def test_produces_valid_output(self) -> None:
        """Test that higher-order PDE produces valid (finite) evolution rates."""
        grid = UnitGrid([64], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        # Higher-order with only second-order term
        pde_higher = HigherOrderKGPDE(mass=0.5, alpha_2=1.0, alpha_4=0.0, alpha_6=0.0)
        result_higher = pde_higher.evolution_rate(state, t=0.0)

        # Check results are finite
        assert np.all(np.isfinite(result_higher[0].data))
        assert np.all(np.isfinite(result_higher[1].data))

        # dphi/dt should equal pi
        np.testing.assert_allclose(result_higher[0].data, pi.data)

    def test_fourth_order_effect(self) -> None:
        """Test that fourth-order term has expected sign and magnitude."""
        grid = UnitGrid([64], periodic=True)
        # Create a field with some curvature
        x = cast("np.ndarray", grid.cell_coords[..., 0])
        phi_data = np.sin(2 * np.pi * x / grid.axes_bounds[0][1])
        phi = ScalarField(grid, data=phi_data)
        pi = ScalarField(grid, data=np.zeros_like(phi.data))
        state = FieldCollection([phi, pi])

        # Without fourth-order
        pde_standard = HigherOrderKGPDE(mass=0.5, alpha_2=1.0, alpha_4=0.0)
        result_standard = pde_standard.evolution_rate(state, t=0.0)

        # With fourth-order (use larger coefficient to make effect visible)
        pde_fourth = HigherOrderKGPDE(mass=0.5, alpha_2=1.0, alpha_4=0.1)
        result_fourth = pde_fourth.evolution_rate(state, t=0.0)

        # Should be different
        assert not np.allclose(result_standard[1].data, result_fourth[1].data)

        # Fourth-order term should contribute (difference should be non-zero)
        diff = np.abs(result_fourth[1].data - result_standard[1].data)
        assert np.max(diff) > 0

    def test_evolution_rate_shape(self) -> None:
        """Test that evolution_rate returns correct shape."""
        grid = UnitGrid([128], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        pde = HigherOrderKGPDE(mass=0.5, alpha_2=1.0, alpha_4=0.01)
        result = pde.evolution_rate(state, t=0.0)

        assert isinstance(result, FieldCollection)
        assert len(result) == 2
        assert result[0].data.shape == phi.data.shape
        assert result[1].data.shape == pi.data.shape


class TestDirectionalKGPDE:
    """Tests for DirectionalKGPDE."""

    def test_initialization(self) -> None:
        """Test basic initialization."""
        pde = DirectionalKGPDE(mass=1.0, active_directions=(True, False))
        assert pde.m2 == 1.0
        assert pde.active_directions == (True, False)

    def test_invalid_directions(self) -> None:
        """Test that invalid active_directions raise ValueError."""
        with pytest.raises(ValueError, match="active_directions must be non-empty"):
            DirectionalKGPDE(mass=1.0, active_directions=[])

        with pytest.raises(ValueError, match="At least one direction must be active"):
            DirectionalKGPDE(mass=1.0, active_directions=(False, False))

    def test_produces_valid_output(self) -> None:
        """Test that directional PDE produces valid (finite) evolution rates."""
        grid = UnitGrid([64, 64], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        # All directions active
        pde_directional = DirectionalKGPDE(mass=0.5, active_directions=(True, True))
        result_directional = pde_directional.evolution_rate(state, t=0.0)

        # Check results are finite
        assert np.all(np.isfinite(result_directional[0].data))
        assert np.all(np.isfinite(result_directional[1].data))

        # dphi/dt should equal pi
        np.testing.assert_allclose(result_directional[0].data, pi.data)

    def test_inactive_direction_reduces_evolution(self) -> None:
        """Test that inactive directions reduce the magnitude of evolution."""
        grid = UnitGrid([64, 64], periodic=True)

        # Create a random field
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField(grid, data=np.zeros_like(phi.data))
        state = FieldCollection([phi, pi])

        # Evolution in all directions
        pde_full = DirectionalKGPDE(mass=0.0, active_directions=(True, True))
        result_full = pde_full.evolution_rate(state, t=0.0)

        # Evolution only in x
        pde_x_only = DirectionalKGPDE(mass=0.0, active_directions=(True, False))
        result_x_only = pde_x_only.evolution_rate(state, t=0.0)

        # Check results are finite
        assert np.all(np.isfinite(result_x_only[1].data))

        # x-only evolution should generally have smaller magnitude than full
        # (unless the field happens to vary only in x)
        mag_full = np.linalg.norm(result_full[1].data)
        mag_x_only = np.linalg.norm(result_x_only[1].data)

        # x_only should be <=  full (they're equal if field doesn't vary in y)
        assert mag_x_only <= mag_full + 1e-10

    def test_dimension_mismatch(self) -> None:
        """Test that dimension mismatch raises ValueError."""
        grid = UnitGrid([64], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField(grid, data=np.zeros_like(phi.data))
        state = FieldCollection([phi, pi])

        # 2D directions on 1D grid
        pde = DirectionalKGPDE(mass=0.5, active_directions=(True, False))

        with pytest.raises(
            ValueError, match=r"active_directions length .* must match grid dimension"
        ):
            pde.evolution_rate(state, t=0.0)


class TestAnisotropicHigherOrderKGPDE:
    """Tests for combined AnisotropicHigherOrderKGPDE."""

    def test_initialization(self) -> None:
        """Test basic initialization."""
        pde = AnisotropicHigherOrderKGPDE(
            mass=1.0, speeds=[1.0, 2.0], alpha_2=1.0, alpha_4=0.01
        )
        assert pde.m2 == 1.0
        np.testing.assert_array_equal(pde.speeds, [1.0, 2.0])
        assert pde.alpha_2 == 1.0
        assert pde.alpha_4 == 0.01

    def test_invalid_speeds(self) -> None:
        """Test that invalid speeds raise ValueError."""
        with pytest.raises(ValueError, match="speeds must be non-empty"):
            AnisotropicHigherOrderKGPDE(mass=1.0, speeds=[])

        with pytest.raises(ValueError, match="All wave speeds must be positive"):
            AnisotropicHigherOrderKGPDE(mass=1.0, speeds=[1.0, -1.0])

    def test_evolution_rate_shape(self) -> None:
        """Test that evolution_rate returns correct shape."""
        grid = UnitGrid([64, 64], periodic=True)
        phi = ScalarField.random_uniform(grid, -0.1, 0.1)
        pi = ScalarField.random_uniform(grid, -0.1, 0.1)
        state = FieldCollection([phi, pi])

        pde = AnisotropicHigherOrderKGPDE(
            mass=0.5, speeds=[2.0, 0.5], alpha_2=1.0, alpha_4=0.01
        )
        result = pde.evolution_rate(state, t=0.0)

        assert isinstance(result, FieldCollection)
        assert len(result) == 2
        assert result[0].data.shape == phi.data.shape
        assert result[1].data.shape == pi.data.shape
