from __future__ import annotations

import logging
from collections.abc import Callable  # add Callable import if not present
from typing import TYPE_CHECKING, Any

from pde import PDEBase, ProgressTracker
from pde.trackers import CallbackTracker

from torsion_gertsenshtein.kgsim.profiling import (
    Timer,
    first_tick_tracker,
    perf_counter,
    print_summary,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from pde import (
        FieldCollection,
    )
    from pde.trackers.base import TrackerBase

    from .config import SimulationConfig


def _build_observer_trackers(
    *,
    config: SimulationConfig,
    extra_observer: Callable[[Any, float], dict[str, Any]] | None,
    snapshot_interval: float | None,
    prof: dict[str, float],
    profile: bool = False,
) -> tuple[TrackerBase | str, ...]:
    """Assemble tracker objects (wrapped as TrackerBase instances)."""
    trackers: list[TrackerBase | str] = []

    if profile:
        ft = first_tick_tracker(prof, "t_first_tick")
        trackers.append(CallbackTracker(ft))

    if config.progress:
        trackers.append(ProgressTracker())

    if extra_observer is not None:
        trackers.append(
            CallbackTracker(
                extra_observer,
                interrupts=1 if snapshot_interval is None else snapshot_interval,
            )
        )

    return tuple(trackers)


def _select_solver_kwargs(
    *, pde: PDEBase, config: SimulationConfig
) -> tuple[str, dict[str, Any]]:
    """Select solver name and kwargs based on config and PDE capabilities.

    Parameters
    ----------
    pde : PDEBase
        The PDE instance to query for backend capabilities.
    config : SimulationConfig
        Simulation configuration object with attributes `solver`, `method`, and `backend`.

    Returns
    -------
    tuple[str, dict[str, Any]]
        A pair of solver name and keyword arguments to pass to the solver.

    Raises
    ------
    ValueError
        If `config.solver` is not one of the supported solver identifiers.
    """
    backend = config.backend
    if backend == "numba" and not hasattr(pde, "_make_pde_rhs_numba"):
        backend = "numpy"

    if config.solver == "scipy":
        solver_name = "scipy"
        solver_kwargs: dict[str, Any] = {
            "method": config.method,
            "backend": backend,
        }
    elif config.solver == "explicit":
        solver_name = "explicit"
        solver_kwargs = {"backend": backend}
    else:
        # defensive runtime guard so the type checker knows `solver` is always bound
        msg = f"Unknown solver: {config.solver!r}"
        raise ValueError(msg)

    return solver_name, solver_kwargs


def _call_solve_with_fallback(  # noqa: PLR0913
    *,
    pde: PDEBase,
    state: FieldCollection,
    t_end: float,
    dt: float | None,
    tracker: tuple[TrackerBase | str, ...],
    solver_name: str,
    solver_kwargs: dict[str, Any],
    prof: dict[str, float],
) -> FieldCollection | tuple[FieldCollection | None, dict[str, Any]] | None:
    """Call pde.solve and retry with numpy backend if numba backend is unsupported."""
    prof["t_call_solve"] = perf_counter()
    try:
        return pde.solve(
            state=state,
            t_range=t_end,
            dt=dt,
            tracker=tracker,
            solver=solver_name,
            **solver_kwargs,
        )
    except NotImplementedError:
        if solver_kwargs.get("backend") == "numba":
            solver_kwargs = {**solver_kwargs, "backend": "numpy"}
            return pde.solve(
                state=state,
                t_range=t_end,
                dt=dt,
                tracker=tracker,
                solver=solver_name,
                **solver_kwargs,
            )
        raise


def _log_profile(*, timer: Timer, prof: dict[str, float]) -> None:
    logger = logging.getLogger(__name__)
    init_delay = prof["t_first_tick"] - prof["t_call_solve"]
    logger.info("init delay until first step: %.3fs", init_delay)
    print_summary(timer.summary())


def run(  # noqa: PLR0913
    *,
    pde: PDEBase,
    state: FieldCollection,
    config: SimulationConfig,
    extra_observer: Callable[[FieldCollection, float], dict[str, Any]] | None = None,
    snapshot_interval: float | None = None,
    profile: bool = False,
) -> FieldCollection:
    """
    Run the Klein-Gordon simulation with the configured solver and observers.

    This function constructs a set of runtime trackers (observers), selects and
    instantiates the requested time integrator, and then delegates the actual time
    integration to the provided PDE object's `solve` method.

    Behavior summary
    - Builds tracker items in order:
        - Optional progress tracker when config.progress is truthy.
        - An optional user-supplied callback observer (extra_observer).
        The trackers are wrapped into a tuple and passed to the PDE solver as the
        `tracker` argument.
    - Selects the solver implementation based on config.solver:
        - "scipy" -> ScipySolver(pde, method=config.method, backend=config.backend)
        - "explicit" -> ExplicitSolver(pde, backend=config.backend)
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
    FieldCollection
            The return value from `pde.solve(...)`. This is always a FieldCollection
            object representing the final state of the simulation.

    Raises
    ------
    RuntimeError
            If the solver returns None instead of a FieldCollection result.

    Notes
    -----
    - The function constructs a tuple of tracker items and passes it as the
        `tracker` argument to `pde.solve`. The exact semantics and call frequency of
        CallbackTracker/ProgressTracker are governed by the tracking infrastructure.
    - Any solver-specific options should be supplied through the SimulationConfig
        instance (e.g. method/backend).
    """
    timer = Timer()

    # profiling bookkeeping
    prof = {"t_call_solve": 0.0, "t_first_tick": 0.0}

    # 1) build trackers
    tracker = _build_observer_trackers(
        config=config,
        extra_observer=extra_observer,
        snapshot_interval=snapshot_interval,
        profile=profile,
        prof=prof,
    )

    timer.mark("trackers_built")

    # 2) solver selection
    solver_name, solver_kwargs = _select_solver_kwargs(pde=pde, config=config)
    timer.mark("solver_selected")

    # 3) call solve (with backend fallback handled inside helper)
    solution_output = _call_solve_with_fallback(
        pde=pde,
        state=state,
        t_end=config.t_end,
        dt=config.dt,
        tracker=tracker,
        solver_name=solver_name,
        solver_kwargs=solver_kwargs,
        prof=prof,
    )
    timer.mark("solve_returned")

    if isinstance(solution_output, tuple):
        result, _info = solution_output
    else:
        result = solution_output

    if result is None:
        msg = "solver returned None"
        raise RuntimeError(msg)

    # 4) optional profiling log
    if profile:
        _log_profile(timer=timer, prof=prof)

    return result
