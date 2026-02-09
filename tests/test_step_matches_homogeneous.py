"""Tests verifying inhomogeneous KG with constant mass matches homogeneous KG.

This module runs two parallel simulations:
(A) Homogeneous KleinGordonPDE with mass m
(B) InhomogeneousKGPDE where the step covers the entire grid with constant m²

These should produce identical results, validating that the inhomogeneous
implementation correctly reduces to the homogeneous case.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import pytest
from pde import FieldCollection, ScalarField

from tidal.kgsim import (
    GridConfig,
    KGParameters,
    KleinGordonPDE,
    SimulationConfig,
    gaussian_pulse,
    make_grid,
    run,
)
from tidal.kgsim.equations import InhomogeneousKGPDE
from tidal.kgsim.profiles import step_region_1d

if TYPE_CHECKING:
    from collections.abc import Callable


# ==================== Helper Functions ====================


def make_recorder(
    store: list[tuple[float, np.ndarray]], keep_every: int = 10
) -> Callable[[FieldCollection, float], dict[str, Any]]:
    """Return a callback that keeps snapshots of phi every `keep_every` calls.

    Parameters
    ----------
    store : list[tuple[float, np.ndarray]]
        List to store (time, field_data) tuples.
    keep_every : int, optional
        Record every Nth callback, by default 10.

    Returns
    -------
    Callable[[FieldCollection, float], dict[str, Any]]
        Callback function for use with simulation runner.
    """
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


def _make_simulation_config(steps: int, dt: float) -> SimulationConfig:
    """Create a deterministic simulation configuration.

    Parameters
    ----------
    steps : int
        Number of time steps.
    dt : float
        Time step size.

    Returns
    -------
    SimulationConfig
        Configuration for explicit numpy-based solver.
    """
    t_end = steps * dt
    return SimulationConfig(
        t_end=t_end,
        dt=dt,
        solver="explicit",
        backend="numpy",
        progress=False,
    )


# ==================== Step-Homogeneous Equivalence Tests ====================


class TestStepMatchesHomogeneous:
    """Tests verifying inhomogeneous KG reduces to homogeneous when mass is constant."""

    # Test parameters: grid_size, length, mass, dt, steps, stride
    GRID_SIZE = 512
    LENGTH = 200.0
    MASS = 0.5
    DT = 0.05
    STEPS = 2000
    STRIDE = 20

    @pytest.mark.parametrize(
        ("grid_size", "length", "mass", "dt", "steps", "stride"),
        [(GRID_SIZE, LENGTH, MASS, DT, STEPS, STRIDE)],
        ids=["standard"],
    )
    def test_step_covering_grid_matches_homogeneous(
        self,
        grid_size: int,
        length: float,
        mass: float,
        dt: float,
        steps: int,
        stride: int,
    ) -> None:
        """Verify step covering entire grid with constant m² matches homogeneous PDE.

        When the step region covers the entire grid with inside_value == outside_value,
        the InhomogeneousKGPDE should produce identical results to KleinGordonPDE.
        """
        # --- grid & IC ---
        grid = make_grid(
            GridConfig(
                dim=1, shape=(grid_size,), bounds=((0.0, length),), periodic=True
            )
        )
        state0 = gaussian_pulse(grid, amplitude=1.0, width=5.0, initial_velocity=0.0)

        # --- (A) homogeneous PDE ---
        pde_a = KleinGordonPDE(KGParameters(mass=mass))
        snapshots_a: list[tuple[float, np.ndarray]] = []
        recorder_a = make_recorder(snapshots_a, keep_every=stride)

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
            config=_make_simulation_config(steps=steps, dt=dt),
            extra_observer=recorder_a,
        )

        # --- (B) step covers entire grid -> constant m²(x) ---
        m2_const = mass**2
        m2_field = step_region_1d(
            grid,
            x0=0.0,
            x1=length,
            inside_value=m2_const,
            outside_value=m2_const,
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
            config=_make_simulation_config(steps=steps, dt=dt),
            extra_observer=recorder_b,
        )

        # --- compare time grids and fields ---
        assert len(snapshots_a) == len(snapshots_b) > 0

        # Allow tiny float jitter in recorded times
        times_a = np.array([t for t, _ in snapshots_a])
        times_b = np.array([t for t, _ in snapshots_b])
        assert np.allclose(times_a, times_b, rtol=0, atol=1e-12)

        # Compare fields at each snapshot
        rtol, atol = 1e-6, 1e-8
        for (t_a, phi_a), (_t_b, phi_b) in zip(snapshots_a, snapshots_b, strict=True):
            assert np.allclose(phi_a, phi_b, rtol=rtol, atol=atol), (
                f"mismatch at t≈{t_a:.6g}"
            )
