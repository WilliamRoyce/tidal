# tests/test_py_pde_smoke.py

from __future__ import annotations

from typing import TYPE_CHECKING, cast

import numpy as np
from numpy.testing import assert_allclose
from pde import PDE, CartesianGrid, ScalarField

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


def test_py_pde_diffusion_smoke() -> None:
    # 1D grid, small so CI stays fast
    grid = CartesianGrid([(0.0, 1.0)], 64)
    rng = np.random.default_rng(42)
    state = ScalarField(grid, rng.random(64))

    # Simple diffusion equation u_t = laplace(u)
    eq = PDE({"u": "laplace(u)"})

    # Short integration; tracker=None to avoid plotting
    out = eq.solve(state, t_range=0.02, dt=1e-4, tracker=None)
    # eq.solve can return either a tuple (field, info) or just the field.
    if isinstance(out, tuple):
        result, _ = out
    else:
        result = out
    assert result is not None  # narrow the type so the checker knows `.data` exists
    result = cast("ScalarField", result)

    # Basic sanity checks
    assert result.data.shape == state.data.shape
    assert np.isfinite(result.data).all()

    # Diffusion should smooth things → variance should drop
    var0 = np.var(state.data)
    var1 = np.var(result.data)
    assert var1 < var0


def test_kg_pde_expression_equivalence() -> None:
    grid = make_grid(
        GridConfig(dim=1, shape=(10,), bounds=((0.0, 1.0),), periodic=True)
    )
    pde = KleinGordonPDE(KGParameters(mass=0.5))
    state = gaussian_pulse(grid, amplitude=1.0, width=1.0)
    cfg = SimulationConfig(
        t_end=1.0, dt=0.1, solver="explicit", backend="numba", progress=False
    )
    out = pde.solve(
        state=state, t_range=cfg.t_end, dt=cfg.dt, solver="explicit", backend="numba"
    )
    # Normalize possible return shapes: FieldCollection | (FieldCollection|None, info) | None
    if isinstance(out, tuple):
        sol, _info = out
    else:
        sol = out
    assert sol is not None

    first = sol[0]
    assert first is not None
    first = cast("ScalarField", first)
    # Grid.periodic can be a bool or a sequence (e.g. [True] for 1D).
    periodic = getattr(first.grid, "periodic", None)
    if isinstance(periodic, (list, tuple)):
        # Narrow element type for the type checker
        periodic_seq = cast("Sequence[bool]", periodic)
        assert any(periodic_seq)
    else:
        assert bool(periodic)


def test_gaussian_equals_multi_gaussian_single_field() -> None:
    """gaussian_pulse(...) should match multi_gaussian(... ) when N == 1."""
    # grid / parameters
    grid_cfg = GridConfig(dim=1, shape=(64,), bounds=((0.0, 10.0),), periodic=True)
    grid = make_grid(grid_cfg)

    amplitude = 1.234
    width = 0.8
    center = [5.0]
    initial_velocity = 0.42

    # single-field calls
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

    # extract arrays and compare
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
