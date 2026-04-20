"""Tests for higher-order finite-difference stencils.

Verifies that the FD accuracy order matches the theoretical convergence
rate for gradient (∂f/∂x) and Laplacian (∂²f/∂x²) operators.

Uses the method of manufactured solutions with sin(kx) on a periodic
domain, where the exact derivatives are known analytically.

Reference: Fornberg, "Generation of Finite Difference Formulas on
Arbitrarily Spaced Grids", Mathematics of Computation 51(184), 1988.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

import numpy as np
import pytest

from tidal.solver.grid import GridInfo
from tidal.solver.operators import (
    cross_derivative,
    directional_laplacian,
    get_fd_order,
    gradient,
    laplacian,
    set_fd_order,
)
from tidal.solver.sparsity import operator_stencil_offsets

if TYPE_CHECKING:
    from collections.abc import Generator

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_fd_order() -> Generator[None, None, None]:  # pyright: ignore[reportUnusedFunction]
    """Reset FD order to default (2) after each test."""
    yield
    set_fd_order(2)


# ---------------------------------------------------------------------------
# Convergence order tests
# ---------------------------------------------------------------------------


def _convergence_slope(
    errors: dict[int, float],
    grid_sizes: list[int],
) -> float:
    """Compute the log-log convergence slope from error vs grid size.

    Uses least-squares fit on log₂(error) vs log₂(dx).
    Reference: LeVeque, "Finite Difference Methods for Ordinary and Partial
    Differential Equations", SIAM, 2007, §2.1 (convergence rate estimation).
    """
    ns = np.array(grid_sizes, dtype=float)
    log_dx = np.log2(1.0 / ns)  # dx proportional to 1/N
    log_err = np.log2(np.array([errors[n] for n in grid_sizes]))
    # Least-squares linear fit: log_err = slope * log_dx + intercept
    coeffs = np.polyfit(log_dx, log_err, 1)
    return float(coeffs[0])


class TestGradientConvergenceOrder:
    """Verify gradient convergence rate matches FD order.

    For sin(2πx) on [0, 2π] periodic, ∂f/∂x = 2π·cos(2πx).
    The L∞ error should scale as O(dx^p) where p is the FD order.
    """

    grid_sizes: ClassVar[list[int]] = [32, 64, 128, 256]

    def _compute_gradient_errors(self, order: int) -> dict[int, float]:
        set_fd_order(order)
        errors: dict[int, float] = {}
        for n in self.grid_sizes:
            grid = GridInfo(
                bounds=((0.0, 2 * np.pi),),
                shape=(n,),
                periodic=(True,),
            )
            x = grid.cell_coords[..., 0]
            f = np.sin(x)
            g = gradient(f, 0, grid)
            exact = np.cos(x)
            errors[n] = float(np.max(np.abs(g - exact)))
        return errors

    def test_order_2(self) -> None:
        """2nd-order gradient converges as O(dx²)."""
        errors = self._compute_gradient_errors(2)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(2.0, abs=0.2), (
            f"Expected slope ~2, got {slope:.2f}"
        )

    def test_order_4(self) -> None:
        """4th-order gradient converges as O(dx⁴)."""
        errors = self._compute_gradient_errors(4)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(4.0, abs=0.2), (
            f"Expected slope ~4, got {slope:.2f}"
        )

    def test_order_6(self) -> None:
        """6th-order gradient converges as O(dx⁶)."""
        errors = self._compute_gradient_errors(6)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(6.0, abs=0.2), (
            f"Expected slope ~6, got {slope:.2f}"
        )


class TestLaplacianConvergenceOrder:
    """Verify Laplacian convergence rate matches FD order.

    Uses sin(3x) on [0, 2π] periodic: ∂²f/∂x² = -9·sin(3x).
    Wavenumber k=3 keeps errors in the measurable range even at high
    order, avoiding premature machine-precision saturation.
    """

    # Order 6 with k=1 hits machine precision at N=256 (~2e-12),
    # corrupting the slope.  Use k=3 + smaller grids instead.
    grid_sizes: ClassVar[list[int]] = [16, 32, 64, 128]

    def _compute_laplacian_errors(self, order: int) -> dict[int, float]:
        set_fd_order(order)
        k = 3.0
        errors: dict[int, float] = {}
        for n in self.grid_sizes:
            grid = GridInfo(
                bounds=((0.0, 2 * np.pi),),
                shape=(n,),
                periodic=(True,),
            )
            x = grid.cell_coords[..., 0]
            f = np.sin(k * x)
            lap = directional_laplacian(f, 0, grid)
            exact = -(k**2) * np.sin(k * x)
            errors[n] = float(np.max(np.abs(lap - exact)))
        return errors

    def test_order_2(self) -> None:
        """2nd-order Laplacian converges as O(dx²)."""
        errors = self._compute_laplacian_errors(2)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(2.0, abs=0.2), (
            f"Expected slope ~2, got {slope:.2f}"
        )

    def test_order_4(self) -> None:
        """4th-order Laplacian converges as O(dx⁴)."""
        errors = self._compute_laplacian_errors(4)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(4.0, abs=0.2), (
            f"Expected slope ~4, got {slope:.2f}"
        )

    def test_order_6(self) -> None:
        """6th-order Laplacian converges as O(dx⁶)."""
        errors = self._compute_laplacian_errors(6)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(6.0, abs=0.3), (
            f"Expected slope ~6, got {slope:.2f}"
        )


# ---------------------------------------------------------------------------
# Accuracy comparison tests
# ---------------------------------------------------------------------------


class TestHigherOrderAccuracy:
    """Verify that higher-order stencils are more accurate at fixed N."""

    def test_gradient_accuracy_improves(self) -> None:
        """At N=64, order 4 is more accurate than order 2, and 6 > 4."""
        n = 64
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi),),
            shape=(n,),
            periodic=(True,),
        )
        x = grid.cell_coords[..., 0]
        f = np.sin(x)
        exact = np.cos(x)

        errors = {}
        for order in (2, 4, 6):
            set_fd_order(order)
            g = gradient(f, 0, grid)
            errors[order] = float(np.max(np.abs(g - exact)))

        # Each doubling of order should improve by ~2 orders of magnitude
        assert errors[4] < errors[2] * 0.01
        assert errors[6] < errors[4] * 0.01

    def test_laplacian_accuracy_improves(self) -> None:
        """At N=64, order 4 is more accurate than order 2, and 6 > 4."""
        n = 64
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi),),
            shape=(n,),
            periodic=(True,),
        )
        x = grid.cell_coords[..., 0]
        f = np.sin(x)
        exact = -np.sin(x)

        errors = {}
        for order in (2, 4, 6):
            set_fd_order(order)
            lap = directional_laplacian(f, 0, grid)
            errors[order] = float(np.max(np.abs(lap - exact)))

        assert errors[4] < errors[2] * 0.01
        assert errors[6] < errors[4] * 0.01


# ---------------------------------------------------------------------------
# 2D tests
# ---------------------------------------------------------------------------


class TestHigherOrder2D:
    """Test higher-order operators on 2D grids."""

    def test_laplacian_2d_periodic(self) -> None:
        """Full Laplacian of sin(x)·sin(y) = -2·sin(x)·sin(y)."""
        n = 32
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi), (0.0, 2 * np.pi)),
            shape=(n, n),
            periodic=(True, True),
        )
        x = grid.cell_coords[..., 0]
        y = grid.cell_coords[..., 1]
        f = np.sin(x) * np.sin(y)
        exact = -2.0 * f

        errors: dict[int, float] = {}
        for order in (2, 4, 6):
            set_fd_order(order)
            lap = laplacian(f, grid)
            err = float(np.max(np.abs(lap - exact)))
            errors[order] = err

        # Higher order should be more accurate (each step ≥ 10× better)
        assert errors[4] < errors[2] * 0.1, (
            f"Order 4 ({errors[4]:.2e}) not better than order 2 ({errors[2]:.2e})"
        )
        assert errors[6] < errors[4] * 0.1, (
            f"Order 6 ({errors[6]:.2e}) not better than order 4 ({errors[4]:.2e})"
        )

    def test_cross_derivative_2d(self) -> None:
        """Cross derivative ∂²f/(∂x∂y) of sin(x)·sin(y) = cos(x)·cos(y)."""
        n = 32
        grid = GridInfo(
            bounds=((0.0, 2 * np.pi), (0.0, 2 * np.pi)),
            shape=(n, n),
            periodic=(True, True),
        )
        x = grid.cell_coords[..., 0]
        y = grid.cell_coords[..., 1]
        f = np.sin(x) * np.sin(y)
        exact = np.cos(x) * np.cos(y)

        errors: dict[int, float] = {}
        for order in (2, 4, 6):
            set_fd_order(order)
            cd = cross_derivative(f, 0, 1, grid)
            err = float(np.max(np.abs(cd - exact)))
            errors[order] = err

        assert errors[4] < errors[2] * 0.1, (
            f"Order 4 ({errors[4]:.2e}) not better than order 2 ({errors[2]:.2e})"
        )
        assert errors[6] < errors[4] * 0.1, (
            f"Order 6 ({errors[6]:.2e}) not better than order 4 ({errors[4]:.2e})"
        )


# ---------------------------------------------------------------------------
# Convergence stress tests (issue #82)
# ---------------------------------------------------------------------------


class TestCrossDerivativeConvergenceOrder:
    """Verify cross-derivative convergence rate matches FD order.

    Uses sin(x)·sin(y) on 2D periodic domain: ∂²f/∂x∂y = cos(x)·cos(y).
    """

    grid_sizes: ClassVar[list[int]] = [16, 32, 64, 128]

    def _compute_cross_errors(self, order: int) -> dict[int, float]:
        set_fd_order(order)
        errors: dict[int, float] = {}
        for n in self.grid_sizes:
            grid = GridInfo(
                bounds=((0.0, 2 * np.pi), (0.0, 2 * np.pi)),
                shape=(n, n),
                periodic=(True, True),
            )
            x = grid.cell_coords[..., 0]
            y = grid.cell_coords[..., 1]
            f = np.sin(x) * np.sin(y)
            exact = np.cos(x) * np.cos(y)
            cd = cross_derivative(f, 0, 1, grid)
            errors[n] = float(np.max(np.abs(cd - exact)))
        return errors

    def test_order_2(self) -> None:
        """2nd-order cross derivative converges as O(dx²)."""
        errors = self._compute_cross_errors(2)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(2.0, abs=0.3), (
            f"Expected slope ~2, got {slope:.2f}"
        )

    def test_order_4(self) -> None:
        """4th-order cross derivative converges as O(dx⁴)."""
        errors = self._compute_cross_errors(4)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(4.0, abs=0.3), (
            f"Expected slope ~4, got {slope:.2f}"
        )

    def test_order_6(self) -> None:
        """6th-order cross derivative converges as O(dx⁶)."""
        errors = self._compute_cross_errors(6)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(6.0, abs=0.4), (
            f"Expected slope ~6, got {slope:.2f}"
        )


class TestDirectionalLaplacianConvergenceOrder2D:
    """Verify directional Laplacian convergence in 2D.

    Uses sin(2x)·cos(y): ∂²f/∂x² = -4·sin(2x)·cos(y).
    Tests that per-axis Laplacian converges correctly even
    when the other axis has nontrivial structure.
    """

    grid_sizes: ClassVar[list[int]] = [16, 32, 64, 128]

    def _compute_dir_lap_errors(self, order: int) -> dict[int, float]:
        set_fd_order(order)
        errors: dict[int, float] = {}
        for n in self.grid_sizes:
            grid = GridInfo(
                bounds=((0.0, 2 * np.pi), (0.0, 2 * np.pi)),
                shape=(n, n),
                periodic=(True, True),
            )
            x = grid.cell_coords[..., 0]
            y = grid.cell_coords[..., 1]
            f = np.sin(2 * x) * np.cos(y)
            exact = -4.0 * np.sin(2 * x) * np.cos(y)
            lap_x = directional_laplacian(f, 0, grid)
            errors[n] = float(np.max(np.abs(lap_x - exact)))
        return errors

    def test_order_2(self) -> None:
        """2nd-order directional Laplacian in 2D converges as O(dx²)."""
        errors = self._compute_dir_lap_errors(2)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(2.0, abs=0.2), (
            f"Expected slope ~2, got {slope:.2f}"
        )

    def test_order_4(self) -> None:
        """4th-order directional Laplacian in 2D converges as O(dx⁴)."""
        errors = self._compute_dir_lap_errors(4)
        slope = _convergence_slope(errors, self.grid_sizes)
        assert slope == pytest.approx(4.0, abs=0.3), (
            f"Expected slope ~4, got {slope:.2f}"
        )


class TestHighWavenumberStress:
    """Stress test: high-wavenumber functions near the Nyquist limit.

    Uses sin(k·x) with k = N/4, which is half the Nyquist frequency.
    This tests the regime where FD errors are largest and numerical
    dispersion is most pronounced.
    """

    def test_gradient_high_k(self) -> None:
        """Gradient of high-k sine at N=256 with order 6 has reasonable error."""
        n = 256
        k = n // 4
        grid = GridInfo(bounds=((0.0, 2 * np.pi),), shape=(n,), periodic=(True,))
        x = grid.cell_coords[..., 0]
        f = np.sin(k * x)
        set_fd_order(6)
        g = gradient(f, 0, grid)
        exact = k * np.cos(k * x)
        rel_err = float(np.max(np.abs(g - exact))) / k
        # At k = N/4, order-6 FD should still achieve < 5% relative error
        assert rel_err < 0.05, f"Relative error {rel_err:.4f} too large"

    def test_laplacian_high_k(self) -> None:
        """Laplacian of high-k sine at N=256 with order 6."""
        n = 256
        k = n // 4
        grid = GridInfo(bounds=((0.0, 2 * np.pi),), shape=(n,), periodic=(True,))
        x = grid.cell_coords[..., 0]
        f = np.sin(k * x)
        set_fd_order(6)
        lap = directional_laplacian(f, 0, grid)
        exact = -(k**2) * np.sin(k * x)
        rel_err = float(np.max(np.abs(lap - exact))) / (k**2)
        # Laplacian is harder at high k; allow 15% relative error
        assert rel_err < 0.15, f"Relative error {rel_err:.4f} too large"


# ---------------------------------------------------------------------------
# State management tests
# ---------------------------------------------------------------------------


class TestFDOrderState:
    """Test the module-level FD order state management."""

    def test_default_order(self) -> None:
        """Default FD order is 2."""
        set_fd_order(2)
        assert get_fd_order() == 2

    def test_set_order_4(self) -> None:
        """Setting order 4 returns 4."""
        set_fd_order(4)
        assert get_fd_order() == 4

    def test_set_order_6(self) -> None:
        """Setting order 6 returns 6."""
        set_fd_order(6)
        assert get_fd_order() == 6

    def test_invalid_order_raises(self) -> None:
        """Invalid FD order raises ValueError."""
        with pytest.raises(ValueError, match="FD order must be 2, 4, or 6"):
            set_fd_order(3)
        with pytest.raises(ValueError, match="FD order must be 2, 4, or 6"):
            set_fd_order(8)


# ---------------------------------------------------------------------------
# Sparsity pattern tests
# ---------------------------------------------------------------------------


class TestSparsityWidensWithOrder:
    """Verify sparsity pattern bandwidth increases with FD order.

    Higher-order stencils have wider footprints, producing more nonzero
    entries in the Jacobian sparsity pattern.
    """

    def test_gradient_offsets_widen(self) -> None:
        """gradient_x offsets grow: [-1,0,1] → [-2..2] → [-3..3]."""
        for order, expected_len in [(2, 3), (4, 5), (6, 7)]:
            set_fd_order(order)
            offsets = operator_stencil_offsets("gradient_x", 1)
            assert len(offsets) == expected_len, (
                f"Order {order}: expected {expected_len} offsets, got {len(offsets)}"
            )

    def test_laplacian_offsets_widen(self) -> None:
        """Full Laplacian in 2D widens with FD order."""
        for order, min_offsets in [(2, 5), (4, 9), (6, 13)]:
            set_fd_order(order)
            offsets = operator_stencil_offsets("laplacian", 2)
            assert len(offsets) >= min_offsets, (
                f"Order {order}: expected >= {min_offsets} offsets, got {len(offsets)}"
            )

    def test_biharmonic_offsets_widen(self) -> None:
        """Biharmonic offsets widen: radius 2 → 4 → 6.

        In 1D, biharmonic has 2*radius+1 offsets: 5, 9, 13.
        """
        for order, expected_1d in [(2, 5), (4, 9), (6, 13)]:
            set_fd_order(order)
            offsets = operator_stencil_offsets("biharmonic", 1)
            assert len(offsets) >= expected_1d, (
                f"Order {order}: expected >= {expected_1d} offsets, got {len(offsets)}"
            )
