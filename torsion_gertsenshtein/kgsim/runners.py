from __future__ import annotations

from typing import TYPE_CHECKING, Any

from pde import (
    ProgressTracker,
)
from pde.trackers import CallbackTracker

from torsion_gertsenshtein.kgsim.observers import total_energy_observer

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import (
        FieldCollection,
    )
    from pde.trackers.base import TrackerBase

    from torsion_gertsenshtein.kgsim.equations import KleinGordonPDE

    from .config import SimulationConfig


def run(
    *,
    pde: KleinGordonPDE,
    state: FieldCollection,
    config: SimulationConfig,
    extra_observer: Callable[[FieldCollection, float], dict[str, Any]] | None = None,
) -> FieldCollection:
    """
    Run the Klein-Gordon simulation with the configured solver and observers.

    This function constructs a set of runtime trackers (observers), selects and
    instantiates the requested time integrator, and then delegates the actual time
    integration to the provided PDE object's `solve` method.

    Behavior summary
    - Builds tracker items in order:
        - Optional progress tracker when cfg.progress is truthy.
        - A callback tracker that records the system's total energy (derived from
            pde.m2 when present).
        - An optional user-supplied callback observer (extra_observer).
        The trackers are wrapped into a tuple and passed to the PDE solver as the
        `tracker` argument.
    - Selects the solver implementation based on cfg.solver:
        - "scipy" -> ScipySolver(pde, method=cfg.method, backend=cfg.backend)
        - "explicit" -> ExplicitSolver(pde, backend=cfg.backend)
    - Calls pde.solve(...) with the assembled solver, tracker, time range and dt.

    Parameters
    ----------
    pde : KleinGordonPDE
            The PDE instance that implements the problem and exposes a `solve` method.
            If it defines an attribute `m2`, that value is used to compute a mass for
            the built-in total-energy observer.
    state : FieldCollection
            The initial field/state collection to pass to the solver.
    config : SimulationConfig
            Simulation configuration. Expected attributes:
                - solver: str, one of "scipy" or "explicit" (controls solver selection).
                - method: optional str, passed to ScipySolver when solver == "scipy".
                - backend: optional value, passed to solver constructors.
                - t_end: time range / endpoint to pass as `t_range` to pde.solve.
                - dt: timestep to pass to pde.solve.
                - progress: bool-like, if true a progress tracker is added.
    extra_observer : Callable[[FieldCollection, float], dict[str, Any]] | None, optional
            Optional user callback invoked periodically by a CallbackTracker. The
            callback should accept the current FieldCollection and a float time and
            return a mapping of observations. If provided it is appended to the
            built-in energy observer.

    Returns
    -------
    Any
            The return value from `pde.solve(...)`. The concrete type depends on the
            PDE/solver implementation (often a solution object or time-history).

    Raises
    ------
    ValueError
            If cfg.solver is not one of the supported solver identifiers.
    RuntimeError
            If the solver returns None instead of a FieldCollection result.

    Notes
    -----
    - The function always installs a total-energy observer (based on pde.m2 when
        available) in addition to any provided extra_observer.
    - The function constructs a tuple of tracker items and passes it as the
        `tracker` argument to `pde.solve`. The exact semantics and call frequency of
        CallbackTracker/ProgressTracker are governed by the tracking infrastructure.
    - Any solver-specific options should be supplied through the SimulationConfig
        instance (e.g. method/backend).
    """
    # --- trackers ---
    trackers: list[TrackerBase | str] = []
    if config.progress:
        trackers.append(ProgressTracker())  # or simply "progress"

    mass = float(getattr(pde, "m2", 0.0)) ** 0.5 if hasattr(pde, "m2") else 0.0
    trackers.append(CallbackTracker(total_energy_observer(mass=mass), interrupts=1))

    if extra_observer is not None:
        trackers.append(CallbackTracker(extra_observer, interrupts=1))

    tracker_param: tuple[TrackerBase | str, ...] = tuple(trackers)

    # --- solver selection: pass CLASS + kwargs (not an instance) ---
    if config.solver == "scipy":
        solver_name = "scipy"
        solver_kwargs: dict[str, Any] = {
            "method": config.method,
            "backend": config.backend,
        }
    elif config.solver == "explicit":
        solver_name = "explicit"
        solver_kwargs = {"backend": config.backend}
    else:
        # defensive runtime guard so the type checker knows `solver` is always bound
        msg = f"Unknown solver: {config.solver!r}"
        raise ValueError(msg)

    # --- run via PDEBase.solve ---
    raw = pde.solve(
        state=state,
        t_range=config.t_end,
        dt=config.dt,
        tracker=tracker_param,
        solver=solver_name,
        **solver_kwargs,
    )

    if isinstance(raw, tuple):
        result, _info = raw
    else:
        result = raw

    if result is None:
        msg = "solver returned None"
        raise RuntimeError(msg)

    return result
