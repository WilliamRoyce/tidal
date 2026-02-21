"""Tests for tidal.solver.operators — spatial FD operators on plain numpy."""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.operators import (
    apply_operator,
    biharmonic,
    cross_derivative,
    directional_laplacian,
    gradient,
    identity,
    laplacian,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _periodic_grid_1d(n: int = 64) -> GridInfo:
    return GridInfo(bounds=((0, 2 * np.pi),), shape=(n,), periodic=(True,))


def _periodic_grid_2d(n: int = 64) -> GridInfo:
    return GridInfo(
        bounds=((0, 2 * np.pi), (0, 2 * np.pi)),
        shape=(n, n),
        periodic=(True, True),
    )


def _neumann_grid_1d(n: int = 64) -> GridInfo:
    return GridInfo(bounds=((0, np.pi),), shape=(n,), periodic=(False,))


def _neumann_grid_2d(n: int = 64) -> GridInfo:
    return GridInfo(
        bounds=((0, np.pi), (0, np.pi)),
        shape=(n, n),
        periodic=(False, False),
    )


# ---------------------------------------------------------------------------
# Periodic 1D
# ---------------------------------------------------------------------------

class TestGradientPeriodic1D:
    def test_sine(self) -> None:
        g = _periodic_grid_1d()
        x = g.axes_coords(0)
        data = np.sin(x)
        result = gradient(data, 0, g)
        # d/dx sin(x) = cos(x), 2nd-order FD error ~ O(dx²)
        np.testing.assert_allclose(result, np.cos(x), atol=5e-3)

    def test_constant(self) -> None:
        g = _periodic_grid_1d()
        data = np.ones(g.shape)
        result = gradient(data, 0, g)
        np.testing.assert_allclose(result, 0, atol=1e-14)


class TestLaplacianPeriodic1D:
    def test_sine(self) -> None:
        g = _periodic_grid_1d()
        x = g.axes_coords(0)
        data = np.sin(x)
        result = laplacian(data, g)
        # ∇² sin(x) = -sin(x)
        np.testing.assert_allclose(result, -np.sin(x), atol=5e-3)

    def test_directional_equals_full(self) -> None:
        g = _periodic_grid_1d()
        x = g.axes_coords(0)
        data = np.sin(x)
        assert np.allclose(laplacian(data, g), directional_laplacian(data, 0, g))


# ---------------------------------------------------------------------------
# Periodic 2D
# ---------------------------------------------------------------------------

class TestOperatorsPeriodic2D:
    @pytest.fixture
    def setup(self) -> tuple[GridInfo, np.ndarray, np.ndarray, np.ndarray]:
        g = _periodic_grid_2d()
        xs, ys = g.coord_arrays()
        data = np.sin(xs) * np.cos(ys)  # sin(x)cos(y)
        return g, data, xs, ys

    def test_gradient_x(self, setup: tuple) -> None:
        g, data, xs, ys = setup
        result = gradient(data, 0, g)
        np.testing.assert_allclose(result, np.cos(xs) * np.cos(ys), atol=5e-3)

    def test_gradient_y(self, setup: tuple) -> None:
        g, data, xs, ys = setup
        result = gradient(data, 1, g)
        np.testing.assert_allclose(result, -np.sin(xs) * np.sin(ys), atol=5e-3)

    def test_laplacian(self, setup: tuple) -> None:
        g, data, _xs, _ys = setup
        result = laplacian(data, g)
        # ∇²[sin(x)cos(y)] = -2 sin(x)cos(y)
        np.testing.assert_allclose(result, -2 * data, atol=5e-2)

    def test_directional_laplacian_x(self, setup: tuple) -> None:
        g, data, xs, ys = setup
        result = directional_laplacian(data, 0, g)
        np.testing.assert_allclose(result, -np.sin(xs) * np.cos(ys), atol=5e-3)

    def test_cross_derivative(self, setup: tuple) -> None:
        g, data, xs, ys = setup
        result = cross_derivative(data, 0, 1, g)
        # ∂²/(∂x∂y) sin(x)cos(y) = -cos(x)sin(y)
        np.testing.assert_allclose(result, -np.cos(xs) * np.sin(ys), atol=5e-2)

    def test_biharmonic(self, setup: tuple) -> None:
        g, data, _xs, _ys = setup
        result = biharmonic(data, g)
        # ∇⁴[sin(x)cos(y)] = 4 sin(x)cos(y)
        np.testing.assert_allclose(result, 4 * data, atol=0.5)

    def test_identity(self, setup: tuple) -> None:
        g, data, _, _ = setup
        result = identity(data, g)
        assert result is data  # no copy


class TestNeumannBC:
    def test_gradient_linear(self) -> None:
        """Gradient of linear function should be constant."""
        g = GridInfo(bounds=((0, 1),), shape=(32,), periodic=(False,))
        x = g.axes_coords(0)
        data = 3.0 * x + 1.0
        result = gradient(data, 0, g, bc="neumann")
        # Interior cells should be exactly 3.0; boundary cells differ
        np.testing.assert_allclose(result[1:-1], 3.0, atol=1e-12)

    def test_laplacian_neumann_compatible(self) -> None:
        """cos(pi*x) satisfies Neumann on [0, pi]."""
        g = _neumann_grid_1d()
        x = g.axes_coords(0)
        data = np.cos(x)  # cos(x) on [0, pi], d/dx = -sin(x) = 0 at boundaries
        result = directional_laplacian(data, 0, g, bc="neumann")
        analytic = -np.cos(x)
        np.testing.assert_allclose(result, analytic, atol=5e-3)

    def test_2d_neumann(self) -> None:
        g = _neumann_grid_2d()
        xs, ys = g.coord_arrays()
        data = np.cos(xs) * np.cos(ys)
        result = laplacian(data, g, bc="neumann")
        analytic = -2 * np.cos(xs) * np.cos(ys)
        np.testing.assert_allclose(result, analytic, atol=5e-3)


# ---------------------------------------------------------------------------
# Mixed BCs (per-axis)
# ---------------------------------------------------------------------------

class TestMixedBC:
    def test_per_axis_bc(self) -> None:
        g = GridInfo(
            bounds=((0, 2 * np.pi), (0, np.pi)),
            shape=(32, 32),
            periodic=(True, False),
        )
        xs, ys = g.coord_arrays()
        data = np.sin(xs) * np.cos(ys)
        result = laplacian(data, g, bc=("periodic", "neumann"))
        np.testing.assert_allclose(result, -2 * data, atol=5e-2)


# ---------------------------------------------------------------------------
# BC validation
# ---------------------------------------------------------------------------

class TestBCValidation:
    def test_unknown_bc(self) -> None:
        g = _periodic_grid_1d()
        data = np.zeros(g.shape)
        with pytest.raises(ValueError, match="Unknown BC type"):
            gradient(data, 0, g, bc="absorbing")

    def test_bc_axis_count_mismatch(self) -> None:
        g = _periodic_grid_2d()
        data = np.zeros(g.shape)
        with pytest.raises(ValueError, match="Expected 2 BC entries"):
            gradient(data, 0, g, bc=("periodic",))


# ---------------------------------------------------------------------------
# Operator registry
# ---------------------------------------------------------------------------

class TestApplyOperator:
    def test_known_operator(self) -> None:
        g = _periodic_grid_1d()
        data = np.ones(g.shape)
        result = apply_operator("identity", data, g)
        assert result is data

    def test_unknown_operator(self) -> None:
        g = _periodic_grid_1d()
        data = np.ones(g.shape)
        with pytest.raises(ValueError, match="Unknown operator"):
            apply_operator("nonexistent", data, g)

    def test_all_registry_entries_callable(self) -> None:
        """Smoke test: every registered operator runs without error on 3D data."""
        from tidal.solver.operators import OPERATOR_REGISTRY

        g = GridInfo(
            bounds=((0, 1), (0, 1), (0, 1)),
            shape=(4, 4, 4),
            periodic=(True, True, True),
        )
        data = np.random.default_rng(42).standard_normal(g.shape)
        for name in OPERATOR_REGISTRY:
            result = apply_operator(name, data, g)
            assert result.shape == g.shape, f"{name} changed shape"


# ---------------------------------------------------------------------------
# Default BC inference from grid
# ---------------------------------------------------------------------------

class TestDefaultBC:
    def test_periodic_grid_infers_periodic(self) -> None:
        g = GridInfo(bounds=((0, 1),), shape=(10,), periodic=(True,))
        data = np.ones(g.shape)
        # Should not raise — periodic BC inferred from grid
        gradient(data, 0, g)

    def test_non_periodic_grid_infers_neumann(self) -> None:
        g = GridInfo(bounds=((0, 1),), shape=(10,), periodic=(False,))
        data = np.ones(g.shape)
        result = gradient(data, 0, g)
        np.testing.assert_allclose(result, 0, atol=1e-14)
