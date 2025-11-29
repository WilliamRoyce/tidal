# tests/test_step_covers_grid_matches_homogeneous.py
from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
from pde import FieldCollection, ScalarField

from torsion_gertsenshtein.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
)
from torsion_gertsenshtein.kgsim.equations import InhomogeneousKGPDE
from torsion_gertsenshtein.kgsim.profiles import step_region_1d

if TYPE_CHECKING:
    from collections.abc import Callable

# ---- helpers ----


def make_recorder(
    store: list[tuple[float, np.ndarray]], keep_every: int = 10
) -> Callable[[FieldCollection, float], dict[str, Any]]:
    """Return a callback(state, t) that keeps snapshots of phi every `keep_every` calls. We still get called every step (interval=1 in your runner), but only keep some."""
    counter = {"n": 0}

    def _callback(state: FieldCollection, t: float) -> dict[str, Any]:
        n = counter["n"]
        if n % keep_every == 0:
            arr = np.asarray(state[0].data)
            if not np.isfinite(arr).all():
                msg = f"Non-finite field at t={t}"
                raise RuntimeError(msg)
            store.append((t, arr.copy()))
        counter["n"] = n + 1
        return {}

    return _callback


# ---- the actual test ----


def test_setup(steps: int, dt: float) -> SimulationConfig:
    # --- configs (deterministic) ---
    t_end = steps * dt
    return SimulationConfig(
        t_end=t_end,
        dt=dt,  # fixed step
        solver="explicit",  # deterministic stepping
        backend="numpy",  # custom PDE runs on numpy; keep both on same backend
        progress=False,
    )


@pytest.mark.parametrize(
    ("nx", "L", "mass", "dt", "steps", "stride"), [(512, 200.0, 0.5, 0.05, 2000, 20)]
)
def test_step_covering_grid_matches_homogeneous(  # noqa: PLR0913, PLR0917
    grid_size: int, length: float, mass: float, dt: float, steps: int, stride: int
) -> None:
    """
    Run two simulations with periodic BCs.

      (A) Homogeneous mass m
      (B) Inhomogeneous PDE where the 'step' covers the whole grid
          with inside_value == outside_value == m^2
    They should match at each recorded snapshot to within small tolerance.
    """
    # --- grid & IC ---
    grid = make_grid(
        GridConfig(dim=1, shape=(grid_size,), bounds=((0.0, length),), periodic=True)
    )
    state0 = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

    # --- (A) homogeneous PDE ---
    pde_a = KleinGordonPDE(KGParameters(mass=mass))
    snapshots_a: list[tuple[float, np.ndarray]] = []
    recorder_a = make_recorder(snapshots_a, keep_every=stride)

    # clone initial state for B so both start identical
    state_a = FieldCollection(
        [
            ScalarField(grid, state0[0].data.copy()),
            ScalarField(grid, state0[1].data.copy()),
        ],
        labels=["phi", "pi"],
    )

    run(
        pde=pde_a,
        state=state_a,
        config=test_setup(steps=steps, dt=dt),
        extra_observer=recorder_a,
    )

    # --- (B) step covers entire grid -> constant m^2(x) ---
    m2_const = mass**2
    m2_field = step_region_1d(
        grid,
        x0=0.0,
        x1=length,  # “covers the whole grid”
        inside_value=m2_const,
        outside_value=m2_const,  # equal -> truly constant
    )
    pde_b = InhomogeneousKGPDE(m2_field=m2_field)

    snapshots_b: list[tuple[float, np.ndarray]] = []
    recorder_b = make_recorder(snapshots_b, keep_every=stride)

    state_b = FieldCollection(
        [
            ScalarField(grid, state0[0].data.copy()),
            ScalarField(grid, state0[1].data.copy()),
        ],
        labels=["phi", "pi"],
    )

    run(
        pde=pde_b,
        state=state_b,
        config=test_setup(steps=steps, dt=dt),
        extra_observer=recorder_b,
    )

    # --- compare time grids and fields ---
    assert len(snapshots_a) == len(snapshots_b) > 0
    # Allow tiny float jitter in recorded times
    times_a = np.array([t for t, _ in snapshots_a])
    times_b = np.array([t for t, _ in snapshots_b])
    assert np.allclose(times_a, times_b, rtol=0, atol=1e-12)

    # Compare fields at each snapshot
    # Tolerances: Explicit solver + different PDE implementations may differ at roundoff;
    # these values are tight but realistic for 1D KG with dt=0.05.
    rtol, atol = 1e-6, 1e-8
    for (t_a, phi_a), (_t_b, phi_b) in zip(snapshots_a, snapshots_b, strict=True):
        assert np.allclose(phi_a, phi_b, rtol=rtol, atol=atol), (
            f"mismatch at t≈{t_a:.6g}"
        )
