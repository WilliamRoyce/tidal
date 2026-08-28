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
:mod:`tidal.cli._sweep` for historical reasons, and were reached by
:mod:`tidal.inference._likelihood` and
:mod:`tidal.measurement._posthoc_audit` importing private names out of a
CLI module.  This module is the seam that hides that: callers outside
:mod:`tidal.cli` import from here, so the implementations can move into
the library layer without touching them.  Relocating them, and
collapsing ``_run_single`` / ``_evaluate_likelihood`` onto a single
policy-free ``run_point(ctx) -> PointOutcome``, is tracked separately --
those are the two hottest paths in the package and want their own
regression campaign.

Scope: the stages of *running a point*.  Expression evaluation
(``FORMULA_NAMESPACE``, ``safe_formula_eval``) is a separate concern and
is still imported from :mod:`tidal.cli._simulate` directly by the two
callers that evaluate baseline formulae.

Beyond the probe stage it re-exports from
:mod:`tidal.measurement._stability`, this module holds no logic of its
own --- only :class:`RunStatus`, the shared vocabulary for what can
happen to a point.
"""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING, Any

from tidal.measurement._stability import (
    PROBE_METADATA_KEYS,
    probe_for_run,
    probe_metadata,
)

if TYPE_CHECKING:
    from argparse import Namespace
    from pathlib import Path

    from tidal.measurement._io import SimulationData
    from tidal.symbolic.json_loader import EquationSystem

__all__ = [
    "PROBE_METADATA_KEYS",
    "RunStatus",
    "measure_from_sim_data",
    "measure_run",
    "parse_bounds",
    "parse_grid_shape",
    "parse_params",
    "probe_for_run",
    "probe_metadata",
    "run_inference_step",
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
    was removed in v0.49.6: rejection on numerical growth is abandoned
    policy, because growth cannot be classified as physics or artifact
    without theory-level analysis (PSALTer, GH #360).  Appears in chains
    and sweeps run with ``--gated`` between 2026-05-10 and v0.49.6."""

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
        if isinstance(exc, (ValueError, TypeError, KeyError, OSError)):
            return cls.MEASUREMENT_ERROR
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
    """Run one simulation to disk. See ``tidal.cli._sweep._simulate_run``."""
    from tidal.cli._sweep import (
        _simulate_run,  # pyright: ignore[reportPrivateUsage]
    )

    return _simulate_run(
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
