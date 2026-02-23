"""Tests for tidal.solver.operators — spatial FD operators on plain numpy."""

from __future__ import annotations

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.operators import (
    AxisBCSpec,
    SideBCSpec,
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


# ---------------------------------------------------------------------------
# SideBCSpec and AxisBCSpec dataclass tests
# ---------------------------------------------------------------------------


class TestSideBCSpec:
    """Unit tests for SideBCSpec ghost-cell formula."""

    def test_neumann_zero_ghost_params(self) -> None:
        """Homogeneous Neumann: ghost = interior (mirror)."""
        bc = SideBCSpec(kind="neumann")
        const, factor = bc.ghost_params(dx=0.1)
        assert const == 0.0
        assert factor == 1.0

    def test_neumann_nonzero_ghost_params(self) -> None:
        """Neumann with derivative D: ghost = dx*D + interior."""
        bc = SideBCSpec(kind="neumann", derivative=2.0)
        const, factor = bc.ghost_params(dx=0.5)
        assert const == pytest.approx(1.0)  # 0.5 * 2.0
        assert factor == 1.0

    def test_dirichlet_zero_ghost_params(self) -> None:
        """Homogeneous Dirichlet: ghost = -interior (antisymmetric)."""
        bc = SideBCSpec(kind="dirichlet")
        const, factor = bc.ghost_params(dx=0.1)
        assert const == 0.0
        assert factor == -1.0

    def test_dirichlet_nonzero_ghost_params(self) -> None:
        """Dirichlet with value V: ghost = 2*V - interior."""
        bc = SideBCSpec(kind="dirichlet", value=3.0)
        const, factor = bc.ghost_params(dx=0.1)
        assert const == pytest.approx(6.0)  # 2 * 3.0
        assert factor == -1.0

    def test_robin_ghost_params(self) -> None:
        """Robin: d_n f + gamma*f = beta."""
        bc = SideBCSpec(kind="robin", value=1.0, gamma=2.0)  # beta=1, gamma=2
        dx = 0.5
        const, factor = bc.ghost_params(dx)
        # denom = gamma*dx + 2 = 2*0.5 + 2 = 3
        # const = 2*dx*beta / denom = 2*0.5*1.0 / 3 = 1/3
        # factor = (2 - gamma*dx) / denom = (2 - 1.0) / 3 = 1/3
        assert const == pytest.approx(1.0 / 3.0)
        assert factor == pytest.approx(1.0 / 3.0)

    def test_robin_reduces_to_neumann(self) -> None:
        """Robin with gamma=0 is equivalent to Neumann."""
        bc_robin = SideBCSpec(kind="robin", value=0.0, gamma=0.0)
        bc_neumann = SideBCSpec(kind="neumann", derivative=0.0)
        dx = 0.25
        c_robin, f_robin = bc_robin.ghost_params(dx)
        c_neumann, f_neumann = bc_neumann.ghost_params(dx)
        assert c_robin == pytest.approx(c_neumann)
        assert f_robin == pytest.approx(f_neumann)

    def test_invalid_kind_raises(self) -> None:
        with pytest.raises(ValueError, match="Invalid BC kind"):
            SideBCSpec(kind="absorbing")

    def test_negative_robin_gamma_raises(self) -> None:
        with pytest.raises(ValueError, match="Robin gamma must be"):
            SideBCSpec(kind="robin", gamma=-1.0)


class TestAxisBCSpec:
    """Unit tests for AxisBCSpec construction and validation."""

    def test_periodic(self) -> None:
        bc = AxisBCSpec(periodic=True)
        assert bc.periodic
        assert bc.low is None
        assert bc.high is None

    def test_non_periodic_with_sides(self) -> None:
        lo = SideBCSpec(kind="dirichlet", value=1.0)
        hi = SideBCSpec(kind="neumann")
        bc = AxisBCSpec(periodic=False, low=lo, high=hi)
        assert not bc.periodic
        assert bc.low is lo
        assert bc.high is hi

    def test_periodic_with_low_raises(self) -> None:
        with pytest.raises(ValueError, match="Periodic axis cannot"):
            AxisBCSpec(periodic=True, low=SideBCSpec(kind="neumann"))

    def test_non_periodic_missing_low_raises(self) -> None:
        with pytest.raises(ValueError, match="requires both low and high"):
            AxisBCSpec(periodic=False, high=SideBCSpec(kind="neumann"))

    def test_non_periodic_missing_high_raises(self) -> None:
        with pytest.raises(ValueError, match="requires both low and high"):
            AxisBCSpec(periodic=False, low=SideBCSpec(kind="neumann"))


class TestStrToAxisBC:
    """Test legacy string -> AxisBCSpec conversion."""

    def test_periodic_string(self) -> None:
        from tidal.solver.operators import _str_to_axis_bc

        bc = _str_to_axis_bc("periodic")
        assert bc.periodic

    def test_neumann_string(self) -> None:
        from tidal.solver.operators import _str_to_axis_bc

        bc = _str_to_axis_bc("neumann")
        assert not bc.periodic
        assert bc.low is not None
        assert bc.low.kind == "neumann"
        assert bc.low.derivative == 0.0

    def test_dirichlet_string(self) -> None:
        from tidal.solver.operators import _str_to_axis_bc

        bc = _str_to_axis_bc("dirichlet")
        assert not bc.periodic
        assert bc.low is not None
        assert bc.low.kind == "dirichlet"
        assert bc.low.value == 0.0

    def test_invalid_string_raises(self) -> None:
        from tidal.solver.operators import _str_to_axis_bc

        with pytest.raises(ValueError, match="Unknown BC type"):
            _str_to_axis_bc("absorbing")


# ---------------------------------------------------------------------------
# Ghost-cell formula operator tests (non-zero values, per-side, Robin)
# ---------------------------------------------------------------------------


class TestNonZeroDirichletBC:
    """Test operators with non-zero Dirichlet boundary values."""

    def test_gradient_nonzero_dirichlet(self) -> None:
        """Gradient near boundary with Dirichlet(v=1.0).

        Field = 1.0 everywhere. With Dirichlet(1.0) on both sides,
        ghost = 2*1.0 - 1.0 = 1.0, so gradient should be zero everywhere.
        """
        g = GridInfo(bounds=((0, 1),), shape=(16,), periodic=(False,))
        data = np.ones(g.shape)
        side = SideBCSpec(kind="dirichlet", value=1.0)
        bc = (AxisBCSpec(periodic=False, low=side, high=side),)
        result = gradient(data, 0, g, bc=bc)
        np.testing.assert_allclose(result, 0.0, atol=1e-14)

    def test_laplacian_nonzero_dirichlet(self) -> None:
        """Laplacian of constant field with matching Dirichlet value = 0."""
        g = GridInfo(bounds=((0, 1),), shape=(16,), periodic=(False,))
        data = 2.0 * np.ones(g.shape)
        side = SideBCSpec(kind="dirichlet", value=2.0)
        bc = (AxisBCSpec(periodic=False, low=side, high=side),)
        result = laplacian(data, g, bc=bc)
        np.testing.assert_allclose(result, 0.0, atol=1e-14)


class TestNonZeroNeumannBC:
    """Test operators with non-zero Neumann derivative values."""

    def test_gradient_nonzero_neumann(self) -> None:
        """Field = constant. With Neumann(deriv=D), ghost = dx*D + interior.

        For f=0 everywhere, Neumann(deriv=1.0):
        ghost_left = dx*1.0 + 0 = dx
        gradient at boundary = (interior - ghost) / (2*dx) = -1/2
        """
        g = GridInfo(bounds=((0, 1),), shape=(16,), periodic=(False,))
        data = np.zeros(g.shape)
        side = SideBCSpec(kind="neumann", derivative=1.0)
        bc = (AxisBCSpec(periodic=False, low=side, high=side),)
        result = gradient(data, 0, g, bc=bc)
        # At left boundary: (f[1] - ghost_left) / (2*dx) = (0 - dx) / (2*dx) = -0.5
        assert result[0] == pytest.approx(-0.5, abs=1e-12)
        # At right boundary: (ghost_right - f[-2]) / (2*dx) = (dx - 0) / (2*dx) = 0.5
        assert result[-1] == pytest.approx(0.5, abs=1e-12)
        # Interior: all zero
        np.testing.assert_allclose(result[1:-1], 0.0, atol=1e-14)


class TestPerSideBC:
    """Test different BC types on left vs right side of same axis."""

    def test_dirichlet_left_neumann_right(self) -> None:
        """Dirichlet(0) on left, Neumann(0) on right."""
        g = GridInfo(bounds=((0, np.pi),), shape=(64,), periodic=(False,))
        x = g.axes_coords(0)
        # sin(x) satisfies: f(0) = 0 (Dirichlet), df/dx(pi) = -1 != 0 (not Neumann-compatible)
        # But sin(x) does satisfy Dirichlet(0) at left
        data = np.sin(x)

        left = SideBCSpec(kind="dirichlet", value=0.0)
        right = SideBCSpec(kind="neumann", derivative=0.0)
        bc = (AxisBCSpec(periodic=False, low=left, high=right),)

        # Check Laplacian: should be close to -sin(x) in interior
        result = directional_laplacian(data, 0, g, bc=bc)
        analytic = -np.sin(x)
        # Interior accuracy (away from boundaries)
        np.testing.assert_allclose(result[2:-2], analytic[2:-2], atol=5e-3)


class TestRobinBC:
    """Test Robin boundary condition: d_n f + gamma*f = beta."""

    def test_robin_ghost_cell_values(self) -> None:
        """Verify ghost cell values match the Robin formula directly."""
        from tidal.solver.operators import _pad_axis

        g = GridInfo(bounds=((0, 1),), shape=(8,), periodic=(False,))
        dx = g.dx[0]
        data = np.arange(8, dtype=float)

        side_lo = SideBCSpec(kind="robin", value=1.0, gamma=2.0)
        side_hi = SideBCSpec(kind="robin", value=0.5, gamma=1.0)
        bc = AxisBCSpec(periodic=False, low=side_lo, high=side_hi)

        padded = _pad_axis(data, 0, bc, dx)
        assert padded.shape == (10,)

        # Check left ghost
        c_lo, f_lo = side_lo.ghost_params(dx)
        expected_left = c_lo + f_lo * data[0]
        assert padded[0] == pytest.approx(expected_left)

        # Check right ghost
        c_hi, f_hi = side_hi.ghost_params(dx)
        expected_right = c_hi + f_hi * data[-1]
        assert padded[-1] == pytest.approx(expected_right)

    def test_robin_laplacian_1d(self) -> None:
        """Robin BC on a quadratic: Laplacian should be constant = 2."""
        g = GridInfo(bounds=((0, 1),), shape=(32,), periodic=(False,))
        x = g.axes_coords(0)
        data = x**2  # f = x^2, f'' = 2

        # f'(0) = 0, f(0) = 0 → Robin: f' + gamma*f = 0 → satisfied for any gamma at x=0
        # f'(1) = 2, f(1) = 1 → Robin: f' + gamma*f = beta → 2 + gamma = beta
        gamma = 1.0
        beta_right = 2.0 + gamma * 1.0  # = 3.0
        left = SideBCSpec(kind="robin", value=0.0, gamma=gamma)
        right = SideBCSpec(kind="robin", value=beta_right, gamma=gamma)
        bc = (AxisBCSpec(periodic=False, low=left, high=right),)

        result = directional_laplacian(data, 0, g, bc=bc)
        # Interior should be close to 2.0
        np.testing.assert_allclose(result[2:-2], 2.0, atol=5e-2)


class TestAxisBCSpecWithOperators:
    """Test passing AxisBCSpec tuples through the full operator pipeline."""

    def test_periodic_axis_bc_spec(self) -> None:
        """AxisBCSpec(periodic=True) should work like 'periodic' string."""
        g = GridInfo(bounds=((0, 2 * np.pi),), shape=(64,), periodic=(True,))
        x = g.axes_coords(0)
        data = np.sin(x)

        result_str = gradient(data, 0, g, bc="periodic")
        result_spec = gradient(data, 0, g, bc=(AxisBCSpec(periodic=True),))
        np.testing.assert_allclose(result_spec, result_str, atol=1e-14)

    def test_neumann_axis_bc_spec(self) -> None:
        """AxisBCSpec with Neumann sides should match 'neumann' string."""
        g = GridInfo(bounds=((0, np.pi),), shape=(32,), periodic=(False,))
        x = g.axes_coords(0)
        data = np.cos(x)

        result_str = directional_laplacian(data, 0, g, bc="neumann")
        side = SideBCSpec(kind="neumann")
        result_spec = directional_laplacian(
            data, 0, g, bc=(AxisBCSpec(periodic=False, low=side, high=side),)
        )
        np.testing.assert_allclose(result_spec, result_str, atol=1e-14)

    def test_dirichlet_axis_bc_spec(self) -> None:
        """AxisBCSpec with Dirichlet sides should match 'dirichlet' string."""
        g = GridInfo(bounds=((0, np.pi),), shape=(32,), periodic=(False,))
        x = g.axes_coords(0)
        data = np.sin(x)

        result_str = directional_laplacian(data, 0, g, bc="dirichlet")
        side = SideBCSpec(kind="dirichlet")
        result_spec = directional_laplacian(
            data, 0, g, bc=(AxisBCSpec(periodic=False, low=side, high=side),)
        )
        np.testing.assert_allclose(result_spec, result_str, atol=1e-14)

    def test_2d_mixed_axis_bc_spec(self) -> None:
        """2D: periodic on axis 0, Dirichlet on axis 1 via AxisBCSpec."""
        g = GridInfo(
            bounds=((0, 2 * np.pi), (0, np.pi)),
            shape=(32, 32),
            periodic=(True, False),
        )
        xs, ys = g.coord_arrays()
        data = np.sin(xs) * np.sin(ys)  # sin vanishes at y=0,pi

        side_d = SideBCSpec(kind="dirichlet")
        bc = (AxisBCSpec(periodic=True), AxisBCSpec(periodic=False, low=side_d, high=side_d))
        result = laplacian(data, g, bc=bc)
        analytic = -2 * data
        np.testing.assert_allclose(result, analytic, atol=5e-2)
