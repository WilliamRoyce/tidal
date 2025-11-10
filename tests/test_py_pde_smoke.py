# tests/test_py_pde_smoke.py

from __future__ import annotations

from typing import cast

import numpy as np
from pde import PDE, CartesianGrid, ScalarField


def test_py_pde_diffusion_smoke() -> None:
    # 1D grid, small so CI stays fast
    grid = CartesianGrid([(0.0, 1.0)], 64)
    rng = np.random.default_rng(42)
    state = ScalarField(grid, rng.random(64))

    # Simple diffusion equation u_t = laplace(u)
    eq = PDE({"u": "laplace(u)"})

    # Short integration; tracker=None to avoid plotting
    out = eq.solve(state, t_range=0.02, dt=1e-3, tracker=None)
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
