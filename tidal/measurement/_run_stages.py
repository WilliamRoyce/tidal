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

Nothing in this module holds logic of its own beyond the probe stage it
re-exports from :mod:`tidal.measurement._stability`.
"""

from __future__ import annotations

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
