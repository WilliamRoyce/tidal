"""Tests for the vectorfield package."""

from __future__ import annotations

from typing import cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import CartesianGrid, FieldCollection

from tidal.symbolic.json_loader import (
    ComponentEquation,
    EquationSystem,
    OperatorTerm,
)
from tidal.vectorfield.config import ComponentFieldParams
from tidal.vectorfield.initial_conditions import (
    ComponentGaussianPulse,
    ComponentPlaneWave,
    create_gaussian_pulse_state,
)

# === Fixtures ===


@pytest.fixture
def grid_1d() -> CartesianGrid:
    """1D periodic grid."""
    return CartesianGrid([(0, 100)], 128, periodic=True)


@pytest.fixture
def grid_2d() -> CartesianGrid:
    """2D periodic grid."""
    return CartesianGrid([(0, 50), (0, 50)], [64, 64], periodic=True)


@pytest.fixture
def em_spec() -> EquationSystem:
    """EM field specification (two components)."""
    return EquationSystem(
        n_components=2,
        dimension=2,
        spatial_dimension=1,
        component_names=("A_0", "A_1"),
        equations=(
            ComponentEquation(
                field_name="A_0",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_0"),
                ),
            ),
            ComponentEquation(
                field_name="A_1",
                field_index=1,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="A_1"),
                ),
            ),
        ),
        mass_matrix=((0.0, 0.0), (0.0, 0.0)),
        coupling_matrix=((0.0, 0.0), (0.0, 0.0)),
        metadata={"gauge": "lorenz"},
    )


@pytest.fixture
def kg_spec() -> EquationSystem:
    """Klein-Gordon specification (single component with mass)."""
    return EquationSystem(
        n_components=1,
        dimension=2,
        spatial_dimension=1,
        component_names=("phi",),
        equations=(
            ComponentEquation(
                field_name="phi",
                field_index=0,
                time_derivative_order=2,
                rhs_terms=(
                    OperatorTerm(coefficient=1.0, operator="laplacian", field="phi"),
                    OperatorTerm(coefficient=-1.0, operator="identity", field="phi"),
                ),
            ),
        ),
        mass_matrix=((1.0,),),
        coupling_matrix=((0.0,),),
        metadata={},
    )


# === ComponentFieldParams Tests ===


class TestComponentFieldParams:
    """Tests for ComponentFieldParams."""

    def test_from_equation_system_em(self, em_spec: EquationSystem) -> None:
        """Test creating params from EM specification."""
        params = ComponentFieldParams.from_equation_system(em_spec)

        num_em_components = 2
        assert params.n_components == num_em_components
        assert params.component_names == ("A_0", "A_1")
        assert params.spatial_dimension == 1
        assert params.is_massless is True

    def test_from_equation_system_kg(self, kg_spec: EquationSystem) -> None:
        """Test creating params from KG specification."""
        params = ComponentFieldParams.from_equation_system(kg_spec)

        assert params.n_components == 1
        assert params.component_names == ("phi",)
        assert params.is_massless is False

    def test_get_component_index(self, em_spec: EquationSystem) -> None:
        """Test getting component index by name."""
        params = ComponentFieldParams.from_equation_system(em_spec)

        assert params.get_component_index("A_0") == 0
        assert params.get_component_index("A_1") == 1

    def test_get_component_index_invalid(self, em_spec: EquationSystem) -> None:
        """Test that invalid component name raises."""
        params = ComponentFieldParams.from_equation_system(em_spec)

        with pytest.raises(ValueError, match="Unknown component"):
            params.get_component_index("A_2")

    def test_get_state_field_indices(self, em_spec: EquationSystem) -> None:
        """Test getting state array indices."""
        params = ComponentFieldParams.from_equation_system(em_spec)

        indices = params.get_state_field_indices()

        assert indices["A_0"] == (0, 1)  # field, momentum
        assert indices["A_1"] == (2, 3)

    def test_validation_n_components(self) -> None:
        """Test that invalid n_components raises."""
        with pytest.raises(ValueError, match="n_components must be at least 1"):
            ComponentFieldParams(
                n_components=0,
                component_names=(),
                spatial_dimension=1,
            )

    def test_validation_component_names_length(self) -> None:
        """Test that mismatched component_names length raises."""
        with pytest.raises(ValueError, match="component_names length"):
            ComponentFieldParams(
                n_components=2,
                component_names=("A_0",),
                spatial_dimension=1,
            )


# === ComponentGaussianPulse Tests ===


class TestComponentGaussianPulse:
    """Tests for ComponentGaussianPulse."""

    def test_create_single_component(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test creating Gaussian pulse in single component."""
        pulse = ComponentGaussianPulse(
            center=(50.0,),
            width=5.0,
            amplitude=2.0,
            active_components={"A_1": 1.0},
        )

        state = pulse.create(grid_1d, em_spec)

        assert isinstance(state, FieldCollection)
        num_em_states = 4  # A_0, Pi_0, A_1, Pi_1
        assert len(state) == num_em_states

        # A_0 should be zero
        assert_allclose(state[0].data, 0.0, atol=1e-10)
        assert_allclose(state[1].data, 0.0, atol=1e-10)

        # A_1 should have Gaussian profile
        assert np.max(state[2].data) > 0
        # Momentum should be zero (stationary pulse)
        assert_allclose(state[3].data, 0.0, atol=1e-10)

    def test_create_all_components(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test creating Gaussian pulse in all components."""
        pulse = ComponentGaussianPulse(
            center=(50.0,),
            width=5.0,
            amplitude=1.0,
            active_components=None,  # All components
        )

        state = pulse.create(grid_1d, em_spec)

        # Both A_0 and A_1 should have the pulse
        assert np.max(state[0].data) > 0
        assert np.max(state[2].data) > 0

        # Same profile in both
        assert_allclose(state[0].data, state[2].data)

    def test_gaussian_peak_at_center(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test that Gaussian peak is at the specified center."""
        center = 50.0
        pulse = ComponentGaussianPulse(
            center=(center,),
            width=5.0,
            amplitude=1.0,
            active_components={"A_0": 1.0},
        )

        state = pulse.create(grid_1d, em_spec)

        # Find peak location
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        peak_idx = np.argmax(state[0].data)
        peak_x = x[peak_idx]

        # Should be close to center (within one grid cell)
        dx = x[1] - x[0]
        assert abs(peak_x - center) < dx

    def test_gaussian_amplitude(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test that Gaussian has correct peak amplitude."""
        amplitude = 3.5
        pulse = ComponentGaussianPulse(
            center=(50.0,),
            width=5.0,
            amplitude=amplitude,
            active_components={"A_0": 1.0},
        )

        state = pulse.create(grid_1d, em_spec)

        # Peak should be approximately amplitude
        assert_allclose(np.max(state[0].data), amplitude, rtol=0.01)

    def test_2d_gaussian(self, grid_2d: CartesianGrid, em_spec: EquationSystem) -> None:
        """Test Gaussian in 2D grid."""
        pulse = ComponentGaussianPulse(
            center=(25.0, 25.0),
            width=5.0,
            amplitude=1.0,
            active_components={"A_0": 1.0},
        )

        state = pulse.create(grid_2d, em_spec)

        # Should have 2D shape
        assert state[0].data.shape == (64, 64)

        # Peak should be at center
        peak_idx = np.unravel_index(np.argmax(state[0].data), state[0].data.shape)
        # Center should be around index 32, 32
        assert abs(peak_idx[0] - 32) <= 1
        assert abs(peak_idx[1] - 32) <= 1

    def test_invalid_width_raises(self) -> None:
        """Test that non-positive width raises."""
        with pytest.raises(ValueError, match="width must be positive"):
            ComponentGaussianPulse(center=(0.0,), width=0.0)

        with pytest.raises(ValueError, match="width must be positive"):
            ComponentGaussianPulse(center=(0.0,), width=-1.0)


# === ComponentPlaneWave Tests ===


class TestComponentPlaneWave:
    """Tests for ComponentPlaneWave."""

    def test_create_plane_wave(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test creating plane wave initial condition."""
        k = 2 * np.pi / 100  # One wavelength across domain
        wave = ComponentPlaneWave(
            wavevector=(k,),
            amplitude=1.0,
            phase=0.0,
            active_components={"A_1": 1.0},
        )

        state = wave.create(grid_1d, em_spec)

        # A_1 field should be cos(kx)
        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        expected_field = np.cos(k * x)
        assert_allclose(state[2].data, expected_field, atol=1e-10)

        # A_1 momentum should be -k * sin(kx) (for right-moving wave)
        expected_momentum = -k * np.sin(k * x)
        assert_allclose(state[3].data, expected_momentum, atol=1e-10)

    def test_plane_wave_with_phase(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test plane wave with non-zero phase."""
        k = 0.1
        phase = np.pi / 4
        wave = ComponentPlaneWave(
            wavevector=(k,),
            amplitude=1.0,
            phase=phase,
            active_components={"A_0": 1.0},
        )

        state = wave.create(grid_1d, em_spec)

        x = cast("np.ndarray", grid_1d.cell_coords[..., 0])
        expected_field = np.cos(k * x + phase)
        assert_allclose(state[0].data, expected_field, atol=1e-10)


# === Convenience Function Tests ===


class TestCreateGaussianPulseState:
    """Tests for create_gaussian_pulse_state function."""

    def test_single_component(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test creating state with single component pulse."""
        state = create_gaussian_pulse_state(
            grid_1d,
            em_spec,
            center=(50.0,),
            width=5.0,
            amplitude=1.5,
            component="A_1",
        )

        # A_0 should be zero
        assert_allclose(state[0].data, 0.0)

        # A_1 should have pulse
        assert np.max(state[2].data) > 0

    def test_all_components(
        self, grid_1d: CartesianGrid, em_spec: EquationSystem
    ) -> None:
        """Test creating state with pulse in all components."""
        state = create_gaussian_pulse_state(
            grid_1d,
            em_spec,
            center=(50.0,),
            width=5.0,
            amplitude=1.0,
            component=None,  # All components
        )

        # Both components should have pulse
        assert np.max(state[0].data) > 0
        assert np.max(state[2].data) > 0
