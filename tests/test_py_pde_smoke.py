# tests/test_py_pde_smoke.py
"""Smoke tests for py-pde integration and Klein-Gordon simulations.

These tests verify basic functionality of the py-pde library and
the Klein-Gordon PDE implementation, including initial condition
generation and parameter validation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
import pytest
from numpy.testing import assert_allclose
from pde import PDE, CartesianGrid, FieldCollection, ScalarField

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
)
from torsion_gertsenshtein.kgsim.config import KGParameters
from torsion_gertsenshtein.kgsim.equations import KleinGordonPDE
from torsion_gertsenshtein.kgsim.initial_conditions import multi_gaussian

if TYPE_CHECKING:
    from collections.abc import Sequence


# ==================== Smoke Tests ====================


class TestPyPDESmoke:
    """Basic smoke tests for py-pde library integration."""

    def test_diffusion_smooths_field(self) -> None:
        """Verify diffusion equation reduces field variance over time."""
        grid = CartesianGrid([(0.0, 1.0)], 64)
        rng = np.random.default_rng(42)
        state = ScalarField(grid, rng.random(64))

        eq = PDE({"u": "laplace(u)"})

        out = eq.solve(state, t_range=0.02, dt=1e-4, tracker=None)
        if isinstance(out, tuple):
            result, _ = out
        else:
            result = out
        assert result is not None
        result = cast("ScalarField", result)

        assert result.data.shape == state.data.shape
        assert np.isfinite(result.data).all()

        # Diffusion should smooth things → variance should drop
        var0 = np.var(state.data)
        var1 = np.var(result.data)
        assert var1 < var0


# ==================== Klein-Gordon PDE Tests ====================


class TestKleinGordonPDE:
    """Tests for Klein-Gordon PDE expression and solution properties."""

    def test_expression_equivalence(self, grid_1d: CartesianGrid) -> None:
        """Verify KleinGordonPDE produces valid FieldCollection with periodic grid."""
        pde = KleinGordonPDE(KGParameters(mass=0.5))
        state = gaussian_pulse(grid_1d, amplitude=1.0, width=1.0)
        cfg = SimulationConfig(
            t_end=1.0, dt=0.1, solver="explicit", backend="numba", progress=False
        )
        out = pde.solve(
            state=state,
            t_range=cfg.t_end,
            dt=cfg.dt,
            solver=cfg.solver,
            backend=cfg.backend,
        )
        if isinstance(out, tuple):
            sol, _info = out
        else:
            sol = out
        assert sol is not None
        assert isinstance(sol, FieldCollection), "Expected FieldCollection from solve"

        first = sol[0]
        assert first is not None
        first = cast("ScalarField", first)
        periodic = getattr(first.grid, "periodic", None)
        if isinstance(periodic, (list, tuple)):
            periodic_seq = cast("Sequence[bool]", periodic)
            assert any(periodic_seq)
        else:
            assert bool(periodic)


# ==================== Multi-Gaussian Validation Tests ====================


class TestMultiGaussianValidation:
    """Tests for multi_gaussian parameter validation."""

    def test_non_1d_grid_raises_value_error(self, grid_2d: CartesianGrid) -> None:
        """Verify multi_gaussian rejects non-1D grids."""
        with pytest.raises(ValueError, match=r"(1D|1-dimension|one-dimensional)"):
            multi_gaussian(grid_2d, amplitudes=[1.0], widths=[1.0])

    @pytest.mark.parametrize("bad_width", [0.0, -1.0])
    def test_non_positive_widths_raise(
        self, grid_1d: CartesianGrid, bad_width: float
    ) -> None:
        """Verify multi_gaussian rejects non-positive width values."""
        with pytest.raises(ValueError, match=r"(width|positive|> ?0)"):
            multi_gaussian(grid_1d, amplitudes=[1.0], widths=[bad_width])

    def test_mismatched_parameter_lengths_raise(self, grid_1d: CartesianGrid) -> None:
        """Verify multi_gaussian rejects mismatched amplitudes/widths lengths."""
        with pytest.raises(ValueError, match=r"(mismatch|length|same length)"):
            multi_gaussian(grid_1d, amplitudes=[1.0, 2.0], widths=[1.0])

    def test_empty_amplitudes_raise(self, grid_1d: CartesianGrid) -> None:
        """Verify multi_gaussian rejects empty parameter lists."""
        with pytest.raises(ValueError, match=r"(empty|non-?empty|at least|>= ?1)"):
            multi_gaussian(grid_1d, amplitudes=[], widths=[])


# ==================== Gaussian Equivalence Tests ====================


class TestGaussianEquivalence:
    """Tests verifying equivalence between gaussian_pulse and multi_gaussian."""

    def test_gaussian_equals_multi_gaussian_single_field(self) -> None:
        """Verify gaussian_pulse matches multi_gaussian when N == 1."""
        grid_cfg = GridConfig(dim=1, shape=(64,), bounds=((0.0, 10.0),), periodic=True)
        grid = make_grid(grid_cfg)

        amplitude = 1.234
        width = 0.8
        center = [5.0]
        initial_velocity = 0.42

        fc_single = gaussian_pulse(
            grid,
            amplitude=amplitude,
            center=center,
            width=width,
            initial_velocity=initial_velocity,
        )
        fc_multi = multi_gaussian(
            grid,
            amplitudes=[amplitude],
            widths=[width],
            centers=center,
            velocities=[initial_velocity],
        )

        phi_single = np.asarray(fc_single[0].data).reshape(grid.shape)
        pi_single = np.asarray(fc_single[1].data).reshape(grid.shape)

        phi_multi = np.asarray(fc_multi[0].data).reshape(grid.shape)
        pi_multi = np.asarray(fc_multi[1].data).reshape(grid.shape)

        assert_allclose(
            phi_single, phi_multi, atol=0.0, rtol=0.0, err_msg="phi arrays do not match"
        )
        assert_allclose(
            pi_single, pi_multi, atol=0.0, rtol=0.0, err_msg="pi arrays do not match"
        )
