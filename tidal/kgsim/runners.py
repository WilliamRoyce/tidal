"""Simulation runner functions for executing Klein-Gordon time evolution.

This module provides high-level run() and run_with_snapshots() functions that
handle solver setup, tracker configuration, and execution.
"""

from __future__ import annotations

import logging
from collections.abc import Callable  # add Callable import if not present
from typing import TYPE_CHECKING, Any

from pde import (
    FieldCollection,
    MemoryStorage,
    PDEBase,
    ProgressTracker,
)
from pde.trackers import CallbackTracker

from tidal.kgsim.profiling import (
    Timer,
    first_tick_tracker,
    perf_counter,
    print_summary,
)

if TYPE_CHECKING:
    from collections.abc import Callable

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
    """Log profiling summary including initialization delay.

    Parameters
    ----------
    timer : Timer
        Timer instance with marked checkpoints.
    prof : dict[str, float]
        Profiling dictionary containing 't_call_solve' and 't_first_tick' keys.
    """
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

    # Type narrowing: assert result is FieldCollection after None check
    assert isinstance(result, FieldCollection), "Expected FieldCollection from solver"

    # 4) optional profiling log
    if profile:
        _log_profile(timer=timer, prof=prof)

    return result


def run_with_snapshots(
    *,
    pde: PDEBase,
    state: FieldCollection,
    config: SimulationConfig,
    snapshot_interval: float | None = None,
    profile: bool = False,
) -> tuple[FieldCollection, MemoryStorage]:
    """Run simulation with automatic snapshot storage using MemoryStorage.

    This is a convenience wrapper around `run()` that automatically sets up
    py-pde's MemoryStorage tracker to collect snapshots at regular intervals,
    eliminating the need for manual snapshot list management.

    Parameters
    ----------
    pde : PDEBase
        The PDE instance defining the evolution equations.
    state : FieldCollection
        Initial state (typically a FieldCollection of [phi, pi]).
    config : SimulationConfig
        Configuration object specifying solver, timestep, backend, etc.
    snapshot_interval : float | None, optional
        Time interval between snapshots. If None, uses config.dt or solver
        interrupts. Default is None.
    profile : bool, optional
        Whether to enable performance profiling. Default is False.

    Returns
    -------
    result : FieldCollection
        Final state after time evolution.
    storage : MemoryStorage
        Storage object containing all snapshots. Access snapshots via:
        - storage.times: array of snapshot times
        - storage[i]: state at snapshot i
        - storage.data: full trajectory data

    Raises
    ------
    RuntimeError
        If the solver returns None instead of a FieldCollection result.

    Examples
    --------
    >>> from tidal.kgsim import (
    ...     GridConfig, make_grid, gaussian_pulse,
    ...     InhomogeneousKGPDE, SimulationConfig, run_with_snapshots
    ... )
    >>> grid = make_grid(GridConfig(dim=1, shape=(128,), bounds=((0, 100),)))
    >>> state = gaussian_pulse(grid, amplitude=1.0, width=2.0)
    >>> pde = InhomogeneousKGPDE(m2_field=ScalarField(grid, data=0.5**2))
    >>> config = SimulationConfig(t_end=50.0, solver="scipy")
    >>> result, storage = run_with_snapshots(
    ...     pde=pde, state=state, config=config, snapshot_interval=1.0
    ... )
    >>> print(f"Collected {len(storage)} snapshots")
    >>> # Extract phi field from each snapshot
    >>> phi_snapshots = [storage[i][0].data for i in range(len(storage))]

    Notes
    -----
    This function replaces the common pattern of:

        snapshots = []
        def record_phi(state_coll, t):
            snapshots.append((t, state_coll[0].data.copy()))
            return {}
        run(..., extra_observer=record_phi)

    With the cleaner:

        result, storage = run_with_snapshots(...)
        snapshots = [(storage.times[i], storage[i][0].data)
                     for i in range(len(storage))]

    The MemoryStorage approach is more efficient and integrates better with
    py-pde's visualization and analysis tools.
    """
    # Create MemoryStorage tracker
    storage = MemoryStorage()
    storage_tracker = storage.tracker(
        interrupts=snapshot_interval if snapshot_interval is not None else 1
    )

    # Build other trackers
    trackers: list[TrackerBase | str] = []

    if profile:
        prof: dict[str, float] = {"t_first_tick": 0.0}
        ft = first_tick_tracker(prof, "t_first_tick")
        trackers.append(CallbackTracker(ft))

    if config.progress:
        trackers.append(ProgressTracker())

    # Add storage tracker
    trackers.append(storage_tracker)

    # Select solver
    solver_name, solver_kwargs = _select_solver_kwargs(pde=pde, config=config)

    # Call solve with storage tracker
    timer = Timer() if profile else None
    prof_dict = {"t_call_solve": 0.0} if profile else {}

    solution_output = _call_solve_with_fallback(
        pde=pde,
        state=state,
        t_end=config.t_end,
        dt=config.dt,
        tracker=tuple(trackers),
        solver_name=solver_name,
        solver_kwargs=solver_kwargs,
        prof=prof_dict,
    )

    if isinstance(solution_output, tuple):
        result, _info = solution_output
    else:
        result = solution_output

    if result is None:
        msg = "solver returned None"
        raise RuntimeError(msg)

    # Type narrowing: assert result is FieldCollection after None check
    assert isinstance(result, FieldCollection), "Expected FieldCollection from solver"

    if profile and timer is not None:
        timer.mark("solve_returned")
        _log_profile(timer=timer, prof=prof_dict)

    return result, storage
