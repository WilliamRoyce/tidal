"""Shared stages for executing one parameter point.

``tidal sweep`` and ``tidal sample`` both answer the same question --
"run the simulation at this point in parameter space and measure it" --
through the same five stages: **resolve params -> probe -> simulate ->
measure -> classify**.  They differ only in the *policy* applied to the
outcome: sweep writes a results row, inference maps the outcome to a
log-likelihood with soft floors.

**A stage used by more than one caller lives HERE.**  If you are about to
copy a stage into a second caller, add it to this module instead.  That
copy is precisely how the sweep and inference paths ended up running two
different tachyonic policies for four months: four of the five stages
were shared and none of them drifted, while the one stage that was
copy-pasted (the pre-flight probe) did -- see GH #454 and the comment it
left behind at the copy site, ``# Pre-flight tachyonic guard -- mirrors
_run_row_inner in _sweep.py``.

Implementations of the simulate/measure stages currently live in
:mod:`tidal.cli._sweep`, and were reached by
:mod:`tidal.inference._likelihood` and
:mod:`tidal.measurement._posthoc_audit` importing private names out of a
CLI module.  This module is the seam that hides that: callers outside
:mod:`tidal.cli` import from here.

**Why they have not simply been moved.**  The obvious next step --- relocate
the ~450 lines of wrappers into this package (GH #480 step 4) --- does not
buy what it looks like it buys.  The wrappers are not the dependency:

* ``simulate_run`` and ``run_inference_step`` call
  ``tidal.cli._simulate._simulate``, the ~3000-line simulation driver.
* ``measure_from_sim_data`` calls eleven private ``_run_*`` measurement
  dispatchers in the ~1000-line :mod:`tidal.cli._measure`.

Moving the wrappers would relocate a large diff across the two hottest
paths in the package and leave both couplings exactly where they are ---
the lazy import would simply be issued from a different file.  The real
work is relocating the *driver* and the *dispatchers* out of
:mod:`tidal.cli`, which is a different and much larger project than
GH #480 step 4 assumed.  Recorded on that issue so it is not re-attempted
on the wrong premise.

Scope: the stages of *running a point*.  Expression evaluation
(``FORMULA_NAMESPACE``, ``safe_formula_eval``) is a separate concern and
is still imported from :mod:`tidal.cli._simulate` directly by the two
callers that evaluate baseline formulae.

:func:`run_point` is that sequence, written once.  Beyond it and the
probe stage re-exported from :mod:`tidal.measurement._stability`, this
module holds no logic of its own --- only :class:`RunStatus`, the shared
vocabulary for what can happen to a point.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from enum import StrEnum
from typing import TYPE_CHECKING, Any

from tidal.measurement._stability import (
    PROBE_METADATA_KEYS,
    probe_for_run,
    probe_metadata,
)

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Callable
    from pathlib import Path

    from tidal.measurement._io import SimulationData
    from tidal.symbolic.json_loader import EquationSystem

__all__ = [
    "PROBE_METADATA_KEYS",
    "PointContext",
    "PointOutcome",
    "RunStatus",
    "measure_from_sim_data",
    "measure_run",
    "parse_bounds",
    "parse_grid_shape",
    "parse_params",
    "probe_for_run",
    "probe_metadata",
    "run_inference_step",
    "run_point",
    "set_single_thread_blas",
    "simulate_run",
]


class RunStatus(StrEnum):
    """Every outcome a parameter point can have, on either entry point.

    ``run_status`` is the column a reader of a sweep ``results.csv`` or a
    chain CSV uses to decide whether a row is usable, and what went wrong
    when it is not.  It had no single definition until GH #480: sweep and
    inference each grew their own vocabulary, sharing only ``success``,
    and three documents described three mutually inconsistent
    taxonomies --- two of them naming tags no code ever emitted.

    A :class:`~enum.StrEnum`, so ``RunStatus.SUCCESS == "success"`` and it
    serializes as the bare string.  Existing readers
    (``row.get("run_status") == "success"``) need no change.

    **Names are not unified across the two paths.**  ``SOLVER_ERROR``
    (sweep) and ``SIMULATION_FAILED`` (inference) mean the same thing
    under different historical names; renaming either would change values
    in newly written CSVs and break comparability with the archive, which
    is a record change and wants its own decision.  This enumeration
    documents the pairing instead of resolving it.
    """

    SUCCESS = "success"
    """The point ran and was measured.  Both paths.

    Note this includes points the stability probe calls tachyonic: since
    GH #454 the probe is a diagnostic and never blocks a run, so growth
    is reported in the ``tachyonic_excess`` / ``n_tachyonic_modes``
    columns rather than by withholding the row."""

    SIMULATION_DIVERGED = "simulation_diverged"
    """The solver raised :exc:`~tidal.solver.SimulationDivergedError`.
    Both paths.

    The failure mode the campaign actually watches for, now that the
    probe does not gate: fields went non-finite or exceeded the norm
    threshold.  Distinct from :attr:`EXCEPTION` on purpose --- before
    GH #480 the inference path swallowed divergence into the bare
    ``except``, so a genuine instability was indistinguishable from a
    ``KeyError`` in measurement code."""

    SIMULATION_FAILED = "simulation_failed"
    """The simulation subprocess exited non-zero.  Inference path.

    Same meaning as :attr:`SOLVER_ERROR` on the sweep path."""

    SOLVER_ERROR = "solver_error"
    """The simulation subprocess exited non-zero.  Sweep path.

    Same meaning as :attr:`SIMULATION_FAILED` on the inference path; the
    two names are historical."""

    KINETIC_ERROR = "kinetic_error"
    """A ``kinetic_coefficient_symbolic`` could not be resolved
    (:exc:`~tidal.solver.KineticEvaluationError`).  Both paths.

    A **configuration** error --- typically a parameter missing from
    ``--param`` --- not a physics result.  It gets its own tag for the
    same reason GH #447 made the exception a ``RuntimeError`` rather than
    a ``ValueError``: so nothing downstream can relabel it as a physics
    verdict.  Before GH #480 the sweep row path caught it in a broad
    ``except RuntimeError`` and recorded it as ``diverged``."""

    MEASUREMENT_ERROR = "measurement_error"
    """The simulation ran but a measurement failed.  Sweep path."""

    METRIC_MISSING = "metric_missing"
    """The requested metric is absent from the simulation output.
    Inference path.  A genuine bug rather than parameter-space signal, so
    it keeps ``-inf`` rather than the soft floor."""

    METRIC_NAN = "metric_nan"
    """The metric evaluated to NaN or Inf.  Inference path.  Soft-floored
    with a distinct tag so post-chain analysis can find these
    specifically."""

    BELOW_NOISE_FLOOR = "below_noise_floor"
    """The simulation returned a finite logL below ``SOFT_FLOOR_LOGL``.
    Inference path.  Observational metadata only --- the value is kept
    verbatim, not clamped (GH #356)."""

    LOGL_MINUS_INF = "logl_minus_inf"
    """The likelihood evaluated to ``-inf``.  Inference path."""

    EXCEPTION = "exception"
    """An unexpected exception.  Both paths.

    The residual bucket: everything not covered by a specific tag above.
    A row landing here is a defect to investigate, not a parameter-space
    feature."""

    # --- historical values -------------------------------------------
    # Nothing emits these any more.  They are retained because archived
    # CSVs contain them and a reader needs to know what they meant.

    DIVERGED = "diverged"
    """**Historical.**  The sweep path's pre-GH #480 catch-all: every
    ``ValueError``/``TypeError``/``KeyError``/``OSError``/``RuntimeError``
    from a run was recorded under this one physics-sounding name,
    configuration errors included.  Split into
    :attr:`SIMULATION_DIVERGED`, :attr:`KINETIC_ERROR`,
    :attr:`MEASUREMENT_ERROR` and :attr:`EXCEPTION`."""

    TACHYONIC_GATED = "tachyonic_gated"
    """**Historical.**  A point the stability probe called tachyonic,
    rejected without being simulated under the ``--gated`` flag.  The flag
    was removed in v0.50.0: rejection on numerical growth is abandoned
    policy, because growth cannot be classified as physics or artifact
    without theory-level analysis (PSALTer, GH #360).  Appears in chains
    and sweeps run with ``--gated`` between 2026-05-10 and v0.50.0."""

    @classmethod
    def live(cls) -> frozenset[RunStatus]:
        """Statuses the current code can emit (i.e. excluding historical)."""
        return frozenset(cls) - {cls.DIVERGED, cls.TACHYONIC_GATED}

    @classmethod
    def from_exception(cls, exc: BaseException) -> RunStatus:
        """Classify an exception raised while running a point.

        One classification shared by every caller that catches broadly,
        so a physics failure, a configuration error and a bug cannot be
        collapsed into one tag again.  The sweep path did exactly that
        until GH #480: a single ``except (ValueError, TypeError, KeyError,
        OSError, RuntimeError, SystemExit)`` recorded all of them as
        ``diverged``.
        """
        from tidal.solver import (
            KineticEvaluationError,
            SimulationDivergedError,
        )

        if isinstance(exc, SimulationDivergedError):
            return cls.SIMULATION_DIVERGED
        if isinstance(exc, KineticEvaluationError):
            return cls.KINETIC_ERROR
        # Deliberately nothing else.  An exception TYPE identifies the stage
        # only when the type is specific to it: a bare OSError is a missing
        # spec file as readily as an unreadable output, and a ValueError is
        # a bad parameter as readily as a failed measurement.  Guessing the
        # stage from the type was an over-claim in the first cut of this
        # helper; :attr:`MEASUREMENT_ERROR` is now set by ``run_point``
        # where it actually knows it is measuring.
        return cls.EXCEPTION


def parse_params(
    param_specs: list[str],
    spec: EquationSystem,
) -> dict[str, float]:
    """Resolve ``--param NAME=VALUE`` strings against a spec's defaults.

    Stage 1 of running a point. See ``tidal.cli._simulate._parse_params``.
    """
    from tidal.cli._simulate import (
        _parse_params,  # pyright: ignore[reportPrivateUsage]
    )

    return _parse_params(param_specs, spec)


def parse_grid_shape(raw: str | None, spatial_dim: int) -> list[int]:
    """Resolve ``--grid-shape``, including the default when it is absent.

    The simulation's own resolution logic
    (``tidal.cli._simulate._parse_grid_shape``), exposed so that anything
    reasoning about the run — notably the stability probe — describes the
    grid the simulation will actually use.  Three separate hardcoded
    fallbacks (256, 64, 256) used to stand in for this and all three
    disagreed with the simulation's 64 (GH #479).
    """
    from tidal.cli._simulate import (
        _parse_grid_shape,  # pyright: ignore[reportPrivateUsage]
    )

    return _parse_grid_shape(raw, spatial_dim)


def parse_bounds(raw: str | None, spatial_dim: int) -> list[tuple[float, float]]:
    """Resolve ``--bounds``, including the default when it is absent.

    Companion to :func:`parse_grid_shape`; see it for why callers must
    not supply their own fallback (GH #479).
    """
    from tidal.cli._simulate import (
        _parse_bounds,  # pyright: ignore[reportPrivateUsage]
    )

    return _parse_bounds(raw, spatial_dim)


def simulate_run(  # noqa: PLR0913
    base_args: Namespace,
    spec_path: Path,
    param_overrides: dict[str, float],
    output_dir: Path,
    grid_shape_override: int | None = None,
    *,
    replicate_seed: int | None = None,
    ic_perturbation: float | None = None,
    spec: EquationSystem | None = None,
) -> tuple[int, float, EquationSystem]:
    """Run one simulation to disk. See ``tidal.cli._sweep.simulate_run``."""
    from tidal.cli._sweep import simulate_run as _impl

    return _impl(
        base_args,
        spec_path,
        param_overrides,
        output_dir,
        grid_shape_override,
        replicate_seed=replicate_seed,
        ic_perturbation=ic_perturbation,
        spec=spec,
    )


def measure_run(  # noqa: PLR0913, PLR0917
    run_dir: Path,
    spec_path: Path,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    spec: EquationSystem | None = None,
) -> dict[str, Any]:
    """Measure a completed run directory. See ``tidal.cli._sweep._measure_run``."""
    from tidal.cli._sweep import (
        _measure_run,  # pyright: ignore[reportPrivateUsage]
    )

    return _measure_run(
        run_dir, spec_path, measurements, source, target, threshold, spec
    )


def measure_from_sim_data(
    data: SimulationData,
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
) -> dict[str, Any]:
    """Measure in-memory simulation output (no disk round-trip).

    See ``tidal.cli._sweep._measure_from_sim_data``.
    """
    from tidal.cli._sweep import (
        _measure_from_sim_data,  # pyright: ignore[reportPrivateUsage]
    )

    return _measure_from_sim_data(data, measurements, source, target, threshold)


def run_inference_step(
    base_args: Namespace,
    spec_path: Path,
    param_overrides: dict[str, float],
    spec: EquationSystem | None = None,
) -> SimulationData:
    """Simulate one point in memory. See ``tidal.cli._sweep.run_inference_step``."""
    from tidal.cli._sweep import run_inference_step as _run_inference_step

    return _run_inference_step(base_args, spec_path, param_overrides, spec)


def set_single_thread_blas() -> None:
    """Pin BLAS to one thread in a worker process.

    See ``tidal.cli._sweep._set_single_thread_blas``.
    """
    from tidal.cli._sweep import (
        _set_single_thread_blas,  # pyright: ignore[reportPrivateUsage]
    )

    _set_single_thread_blas()


# ---------------------------------------------------------------------------
# Running one point
# ---------------------------------------------------------------------------


class _MeasurementStageError(Exception):
    """Internal marker: the failure happened in the MEASURE stage.

    Carries the original exception so the outcome can report it.  Exists
    so the stage is identified by where it was raised rather than guessed
    from the exception's type — a ``ValueError`` is a bad parameter as
    readily as a failed measurement.
    """

    def __init__(self, cause: BaseException) -> None:
        super().__init__(str(cause))
        self.cause = cause


@dataclass(frozen=True)
class PointContext:
    """Everything needed to run one point in parameter space. No policy.

    Carries the differences between the two entry points explicitly, so
    they are visible as data rather than as two divergent code paths:
    ``run_dir`` is a durable artifact for a sweep and a temp directory for
    the inference disk backend; ``replicate_seed`` / ``ic_perturbation``
    are sweep-only; ``sampled_params`` is inference-only.
    """

    spec_path: Path
    base_args: Namespace
    param_overrides: dict[str, float]
    measurements: set[str]
    source: tuple[str, ...] | None
    target: tuple[str, ...] | None
    threshold: float
    run_dir: Path | None = None
    grid_shape_override: int | None = None
    replicate_seed: int | None = None
    ic_perturbation: float | None = None
    #: Sampled (BSM) parameter names, stashed on the spec's metadata so the
    #: modal convolution-block cache can tell which symbols vary per call
    #: from which are geometry-fixed (GH #384 Phase A').  Memory backend only.
    sampled_params: tuple[str, ...] | None = None
    #: Hook for the GH #421 "probe unavailable" marker.  The inference path
    #: passes its warn-once builder so a position-dependent-kinetic spec
    #: records WHY the diagnostic columns are absent; sweep passes nothing.
    unavailable_meta: Callable[[BaseException], dict[str, Any]] | None = None


@dataclass(frozen=True)
class PointOutcome:
    """What happened to a point.  Still no policy — nothing here is a verdict.

    ``status`` says what occurred; it does not say what to do about it.
    Mapping an outcome to a results row or to a log-likelihood is the
    caller's business, and the two callers answer differently: a sweep
    records the row, inference applies a soft floor whose noise is seeded
    on ``theta`` (GH #408) and so cannot be computed here.
    """

    status: RunStatus
    metrics: dict[str, Any] = field(default_factory=dict)
    stability: dict[str, Any] = field(default_factory=dict)
    wall_time_s: float = 0.0
    spec: EquationSystem | None = None
    exit_code: int | None = None
    #: The exception, when one was caught.  Retained rather than
    #: stringified so a caller that wants the original propagation (the
    #: sweep loop builds its error rows from the exception, with the swept
    #: values it holds and ``run_point`` does not) can re-raise it intact.
    exception: BaseException | None = None

    @property
    def ok(self) -> bool:
        return self.status is RunStatus.SUCCESS


def run_point(ctx: PointContext, *, backend: str = "disk") -> PointOutcome:
    """Run one parameter point: probe -> simulate -> measure -> classify.

    The single definition of the sequence both ``tidal sweep`` and
    ``tidal sample`` execute.  They used to hand-roll it separately, which
    is why a policy change to the probe stage reached only one of them for
    four months (GH #454) and why they emitted different metadata schemas
    and different status vocabularies (GH #480).

    **Total**: every path returns an outcome, including the failure ones.
    That is what lets both callers be straight-line mappings instead of
    re-implementing exception handling, and it is why
    :meth:`RunStatus.from_exception` exists.

    **Policy-free**: it never rejects a point, never soft-floors, never
    writes a row.  In particular the stability probe is run for its
    diagnostics and its verdict is *recorded*, never acted on — rejection
    on tachyonic growth is abandoned policy (GH #454).

    Parameters
    ----------
    backend
        ``"disk"`` runs the simulation to ``ctx.run_dir`` and measures from
        it.  ``"memory"`` runs in-process and measures the arrays directly,
        skipping the disk round-trip; it is the inference default, with
        the disk path kept for bisectability (GH #269).
    """
    stability: dict[str, Any] = {}
    try:
        _result, stability = probe_for_run(
            ctx.spec_path,
            ctx.base_args,
            ctx.param_overrides,
            source=ctx.source,
            target=ctx.target,
            measurements=ctx.measurements,
            grid_shape_override=ctx.grid_shape_override,
            unavailable_meta=ctx.unavailable_meta,
        )

        if backend == "memory":
            spec, metrics, wall_time_s = _run_in_memory(ctx)
            exit_code = 0
        else:
            spec, metrics, wall_time_s, exit_code = _run_on_disk(ctx)
            if exit_code != 0:
                return PointOutcome(
                    status=RunStatus.SIMULATION_FAILED,
                    stability=stability,
                    wall_time_s=wall_time_s,
                    spec=spec,
                    exit_code=exit_code,
                )
    except _MeasurementStageError as wrapper:
        return PointOutcome(
            status=RunStatus.MEASUREMENT_ERROR,
            stability=stability,
            exception=wrapper.cause,
        )
    except Exception as exc:  # noqa: BLE001
        return PointOutcome(
            status=RunStatus.from_exception(exc),
            stability=stability,
            exception=exc,
        )

    return PointOutcome(
        status=RunStatus.SUCCESS,
        metrics=metrics,
        stability=stability,
        wall_time_s=wall_time_s,
        spec=spec,
        exit_code=exit_code,
    )


def _run_in_memory(
    ctx: PointContext,
) -> tuple[EquationSystem | None, dict[str, Any], float]:
    """Simulate in-process and measure the arrays directly.

    Raises
    ------
    _MeasurementStageError
        If the measure stage raises, so the outcome can attribute the
        failure to that stage rather than guessing from the type.
    """
    import dataclasses as _dc

    from tidal.symbolic._spec_cache import load_spec_cached

    spec_arg: EquationSystem | None = None
    if ctx.sampled_params is not None:
        raw_spec = load_spec_cached(ctx.spec_path)
        spec_arg = _dc.replace(
            raw_spec,
            metadata={
                **raw_spec.metadata,
                "_inference_sampled_params": tuple(ctx.sampled_params),
            },
        )

    started = time.perf_counter()
    sim_data = run_inference_step(
        ctx.base_args,
        ctx.spec_path,
        ctx.param_overrides,
        spec_arg,
    )
    try:
        metrics = measure_from_sim_data(
            sim_data,
            ctx.measurements,
            ctx.source,
            ctx.target,
            ctx.threshold,
        )
    except Exception as exc:
        raise _MeasurementStageError(exc) from exc
    return sim_data.spec, metrics, time.perf_counter() - started


def _run_on_disk(
    ctx: PointContext,
) -> tuple[EquationSystem | None, dict[str, Any], float, int]:
    """Simulate to ``ctx.run_dir`` and measure from it.

    Raises
    ------
    ValueError
        If ``ctx.run_dir`` is unset; the disk backend has nowhere to write.
    _MeasurementStageError
        If the measure stage raises, so the outcome can attribute the
        failure to that stage rather than guessing from the type.
    """
    if ctx.run_dir is None:
        msg = "run_point(backend='disk') requires PointContext.run_dir"
        raise ValueError(msg)

    exit_code, wall_time_s, spec = simulate_run(
        ctx.base_args,
        ctx.spec_path,
        ctx.param_overrides,
        ctx.run_dir,
        ctx.grid_shape_override,
        replicate_seed=ctx.replicate_seed,
        ic_perturbation=ctx.ic_perturbation,
    )
    if exit_code != 0:
        return spec, {}, wall_time_s, exit_code

    try:
        metrics = measure_run(
            ctx.run_dir,
            ctx.spec_path,
            ctx.measurements,
            ctx.source,
            ctx.target,
            ctx.threshold,
            spec,
        )
    except Exception as exc:
        raise _MeasurementStageError(exc) from exc
    return spec, metrics, wall_time_s, exit_code
