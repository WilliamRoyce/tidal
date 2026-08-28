"""Simulation-based likelihood function for Bayesian inference.

Wraps the shared run stages (:mod:`tidal.measurement._run_stages`) as a
callable ``log_likelihood(theta) -> float``.

Each likelihood evaluation runs one PDE simulation and extracts a scalar
metric (e.g. ``P_max``).

In the v3 architecture (post-2026-05-08 supervisor pivot — see
``docs/V3_ARCHITECTURE.md``) the only ``-inf`` return is for genuine
metric-key-missing bugs.  Numerical failures (sim divergence, NaN, exception)
return a noisy soft floor ``logL = SOFT_FLOOR_LOGL + Normal(0, sigma_explore)``
so PolyChord sees a gradient in the failure region.  The Hwang-Noh
perturbativity gate (``P_max > 0.5 -> -inf``) is removed entirely; large
``P_max`` is faithful linearized-PDE output, not a probability-conservation
violation.  The pre-flight tachyonic probe is run as a metadata measurement
only — nothing acts on its verdict (GH #454; the ``--gated`` escape hatch
was removed in v0.50.0).
"""

from __future__ import annotations

import math
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Sequence


# Soft penalty floor for numerical-failure samples (sim divergence, NaN
# metric, broad exception).  The chain still rejects these in practice
# (``exp(-100) ~ 4e-44`` so they can't out-compete a real sample), but
# PolyChord sees a finite logL and can navigate the failure landscape.
# Per the supervisor's "noise-on-floor" suggestion the actual return value
# is ``SOFT_FLOOR_LOGL + Normal(0, sigma_explore)`` with default
# ``sigma_explore = 1.0`` so the sampler doesn't see a flat plateau where
# it can't tell directions apart.
SOFT_FLOOR_LOGL: float = -100.0

# Once-per-process latch for the GH #421 "stability probe unavailable"
# warning: the condition is a property of the spec, and repeating the
# warning on every likelihood evaluation would flood chain logs.
_posdep_probe_warned = False


def _posdep_probe_unavailable_meta(exc: BaseException) -> dict[str, Any]:
    """Metadata + once-per-process warning for a GH #421-guarded spec.

    The modal-based stability probe raises ``NotImplementedError`` for
    specs with position-dependent kinetic coefficients (GH #421/#427).
    That is a per-spec property, not a per-sample failure: warn ONCE and
    return an explicit ``stability_profile`` marker so chain CSVs record
    WHY the gamma_eff / borderline_stability / n_tachyonic_modes columns
    are absent instead of silently losing them (post-merge review H4 —
    previously this fell into a bare ``except Exception`` at debug level
    and the gate vanished without trace).
    """
    import logging

    global _posdep_probe_warned  # noqa: PLW0603
    if not _posdep_probe_warned:
        _posdep_probe_warned = True
        logging.getLogger("tidal.inference").warning(
            "Pre-flight tachyonic probe unavailable for this spec "
            "(position-dependent kinetic coefficient — GH #421/#427); "
            "sampling continues WITHOUT the stability gate: %s",
            exc,
        )
    return {"stability_profile": "unavailable-posdep-kinetic"}


def _soft_floor_logl(
    sigma_explore: float,
    floor: float = SOFT_FLOOR_LOGL,
    rng: np.random.Generator | None = None,
) -> float:
    """Return the soft penalty floor with optional exploration noise.

    Returns ``floor + Normal(0, sigma_explore)``.  When ``sigma_explore <= 0``
    the noise is disabled and the bare floor is returned (useful for ablation
    tests via ``--soft-floor-noise=0``).  The floor defaults to
    ``SOFT_FLOOR_LOGL`` but can be overridden per-run via
    ``--soft-floor-logl`` when the default -100 contaminates logZ for
    baseline-normalized likelihoods (see issue #375).
    """
    if sigma_explore <= 0.0:
        return floor
    if rng is None:
        rng = np.random.default_rng()
    return floor + float(rng.normal(0.0, sigma_explore))


def _floor_rng(noise_seed: int, theta: Sequence[float]) -> np.random.Generator:
    """Derive the soft-floor generator from the seed and the sample itself.

    Seeding on ``theta`` makes the floor a deterministic function of the
    parameter point, independent of which worker evaluated it or in what
    order.  A call counter cannot do this: ``_likelihood_worker_init``
    gives every ``multiprocessing`` fork a counter starting at zero, so
    the k-th task on each worker would draw the same noise, and
    ``pool.imap`` dispatch order is not reproducible in the first place.
    Storing a ``Generator`` on ``LikelihoodConfig`` fails the same way ---
    the config is pickled into every fork, so all workers would resume
    from one state.

    See issue #408.
    """
    theta_words = np.frombuffer(
        np.asarray(theta, dtype=np.float64).tobytes(),
        dtype=np.uint32,
    )
    return np.random.default_rng(
        np.random.SeedSequence([int(noise_seed), *(int(w) for w in theta_words)])
    )


@dataclass(frozen=True)
class LikelihoodConfig:
    """Configuration for the simulation likelihood.

    Parameters
    ----------
    metric : str
        Name of the metric to extract (e.g. ``"P_max"``, ``"L_mix"``).
    likelihood_type : str
        One of ``"maximize"``, ``"minimize"``, ``"extremize"``,
        ``"gaussian"``, ``"threshold"``.
    target : float
        Target value for gaussian likelihood.
    sigma : float
        Uncertainty for gaussian likelihood.
    min_value : float
        Minimum value for threshold likelihood.
    baseline_formula : str | None
        Baseline formula for ``extremize`` type.  Evaluated per-point
        using the current parameter values + math functions.  Example:
        ``"sin(kappa * B0 * t_end / 2)**2"``.
    """

    metric: str
    likelihood_type: str = "maximize"
    target: float = 0.0
    sigma: float = 1.0
    min_value: float = 0.0
    baseline_formula: str | None = None
    soft_floor_noise_sigma: float = 1.0
    """Standard deviation of the Gaussian noise added to the soft penalty
    floor (sim divergence / NaN / exception).  Default 1.0 gives the
    sampler a finite gradient in the failure region; set to 0 via
    ``--soft-floor-noise 0`` to disable noise for ablation tests."""
    soft_floor_logl: float = SOFT_FLOOR_LOGL
    """Position of the soft penalty floor (default: -100).  When the prior
    contains many tachyonic / failing samples and the natural physics logL
    is near 0 (e.g. baseline-normalized ``maximize`` runs), -100 drags logZ
    to -100 and makes the PolyChord precision criterion unreachable.  Use
    ``--soft-floor-logl -15`` (or similar) to keep logZ in a sensible range
    so that ``precision_criterion × |logZ|`` stays achievable.
    See issue #375."""

    noise_seed: int = 0
    """Seed for the soft-floor exploration noise (default: 0).

    Combined with the sample's own ``theta`` to derive a generator, so a
    given parameter point always receives the same floor value.  This is
    a correctness property as much as a reproducibility one: nested
    sampling assumes the likelihood is a deterministic function of the
    parameters, and an unseeded floor lets a sample's logL change between
    live-point insertion and later re-evaluation.

    Seeding on ``theta`` rather than a call counter is deliberate --- a
    counter restarts at zero in every ``multiprocessing`` worker, so the
    k-th task on every worker would draw identical noise, and ``imap``
    dispatch order is not reproducible anyway.  See issue #408."""

    def __post_init__(self) -> None:
        valid = {"gaussian", "threshold", "maximize", "minimize", "extremize"}
        if self.likelihood_type not in valid:
            msg = f"Unknown likelihood type '{self.likelihood_type}'. Must be one of {sorted(valid)}."
            raise ValueError(msg)
        if self.likelihood_type == "gaussian" and self.sigma <= 0:
            msg = f"Gaussian sigma must be positive, got {self.sigma}"
            raise ValueError(msg)
        if self.likelihood_type == "extremize" and not self.baseline_formula:
            msg = (
                "extremize likelihood requires --baseline-formula "
                '(e.g. "sin(kappa * B0 * t_end / 2)**2")'
            )
            raise ValueError(msg)


def parse_likelihood(
    spec: str,
    *,
    baseline_formula: str | None = None,
    soft_floor_noise_sigma: float = 1.0,
    soft_floor_logl: float = SOFT_FLOOR_LOGL,
    noise_seed: int = 0,
) -> LikelihoodConfig:
    """Parse a CLI likelihood specification string.

    Format: ``METRIC:TYPE[:ARGS]``

    Examples::

        "P_max:maximize"
        "P_max:minimize"
        "P_max:extremize"          # requires baseline_formula
        "P_max:gaussian:0.5:0.1"
        "P_max:threshold:0.01"

    Parameters
    ----------
    spec : str
        Likelihood specification string.
    baseline_formula : str | None
        Baseline formula for ``extremize`` type (passed separately via
        ``--baseline-formula`` CLI flag).
    soft_floor_noise_sigma : float
        Standard deviation of the Gaussian noise added to the soft penalty
        floor (sim divergence / NaN / exception).  Default 1.0; set to 0
        via ``--soft-floor-noise 0`` to disable noise for ablation tests.

    Raises
    ------
    ValueError
        If the specification string is malformed.
    """
    parts = spec.split(":")
    if len(parts) < 2:
        msg = f"Likelihood spec needs METRIC:TYPE, got '{spec}'"
        raise ValueError(msg)

    metric = parts[0]
    ltype = parts[1]

    common = {
        "soft_floor_noise_sigma": soft_floor_noise_sigma,
        "soft_floor_logl": soft_floor_logl,
        "noise_seed": noise_seed,
    }

    if ltype == "gaussian":
        if len(parts) < 4:
            msg = (
                f"Gaussian likelihood needs METRIC:gaussian:TARGET:SIGMA, got '{spec}'"
            )
            raise ValueError(msg)
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="gaussian",
            target=float(parts[2]),
            sigma=float(parts[3]),
            **common,
        )
    if ltype == "threshold":
        if len(parts) < 3:
            msg = f"Threshold likelihood needs METRIC:threshold:MIN_VALUE, got '{spec}'"
            raise ValueError(msg)
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="threshold",
            min_value=float(parts[2]),
            **common,
        )
    if ltype == "maximize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="maximize",
            baseline_formula=baseline_formula,
            **common,
        )
    if ltype == "minimize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="minimize",
            baseline_formula=baseline_formula,
            **common,
        )
    if ltype == "extremize":
        return LikelihoodConfig(
            metric=metric,
            likelihood_type="extremize",
            baseline_formula=baseline_formula,
            **common,
        )

    msg = (
        f"Unknown likelihood type '{ltype}'. "
        "Use: maximize, minimize, extremize, gaussian, threshold."
    )
    raise ValueError(msg)


def compute_log_likelihood(
    metric_value: float,
    config: LikelihoodConfig,
    eval_params: dict[str, float] | None = None,
) -> float:
    """Compute log-likelihood from a metric value and config.

    Parameters
    ----------
    metric_value : float
        The measured value of the metric from simulation.
    config : LikelihoodConfig
        Likelihood configuration.
    eval_params : dict | None
        Current parameter values for formula evaluation (needed by
        ``extremize`` type to compute the per-point baseline).
    """
    if math.isnan(metric_value) or math.isinf(metric_value):
        return -math.inf

    # v3 (post-2026-05-08 supervisor pivot): no upper cap on ``P_max`` at
    # any value.  Large ``P_max`` is faithful linearized-PDE output, not
    # a probability-conservation violation — interaction energies become
    # negative to compensate.  Downstream interpretation handles physical
    # meaning; the chain records every value.  See ``docs/V3_ARCHITECTURE.md``.

    if config.likelihood_type == "gaussian":
        return -0.5 * ((metric_value - config.target) / config.sigma) ** 2

    if config.likelihood_type == "threshold":
        return 0.0 if metric_value >= config.min_value else -math.inf

    # If a baseline formula is provided, use the log-amplification
    # log(A) = log(metric / baseline) as the physics quantity.  This is
    # the Gertsenshtein-normalized conversion (dimensionless, centered
    # on 0 for no-amplification, positive for enhancement, negative for
    # suppression).  Using log-ratio tames the dynamic range so evidence
    # isn't dominated by a few large outliers.
    if config.baseline_formula:
        baseline = _eval_baseline(config.baseline_formula, eval_params)
        if baseline is None or baseline <= 0 or metric_value <= 0:
            return -math.inf
        log_amplification = math.log(metric_value / baseline)
        if config.likelihood_type == "maximize":
            return log_amplification  # find max amplification
        if config.likelihood_type == "minimize":
            return -log_amplification  # find max suppression
        if config.likelihood_type == "extremize":
            return abs(log_amplification)  # find max deviation (both sides)
        # gaussian/threshold fall through to raw-metric branches below

    if config.likelihood_type == "minimize":
        return -metric_value

    if config.likelihood_type == "extremize":
        # extremize requires a baseline (enforced by __post_init__);
        # this branch is unreachable but kept for safety.
        return -math.inf

    # maximize: use metric value directly as log-likelihood
    return metric_value


def _eval_baseline(
    formula: str | None,
    params: dict[str, float] | None,
) -> float | None:
    """Evaluate a baseline formula with the given parameter values.

    Reuses the same ``FORMULA_NAMESPACE`` (sin, cos, sqrt, pi, etc.)
    as the sweep-results derived-columns computation.

    Returns
    -------
    float | None
        The baseline, or ``None`` when no formula is configured.  A
        formula that evaluates to a non-positive number returns it
        unchanged: that is a legitimate parameter-space outcome and the
        caller soft-floors it.

    Both exception types propagate from ``safe_formula_eval``: a name
    nothing supplies raises ``ValueError``, a disallowed construct
    ``TypeError``.  Those are configuration errors, not physics results,
    and are deliberately not swallowed here — see GH #407.
    """
    if formula is None:
        return None

    from tidal.cli._simulate import (
        FORMULA_NAMESPACE,  # pyright: ignore[reportPrivateUsage]
    )

    ns: dict[str, object] = {**FORMULA_NAMESPACE}
    if params:
        ns.update(params)

    from tidal.cli._simulate import safe_formula_eval

    # Deliberately NOT wrapped in a bare except returning None.  Callers
    # treat a None baseline exactly like a baseline of zero — both give
    # -inf, which v3 soft-floors — so swallowing here made a CONFIGURATION
    # error (a name nothing supplies) indistinguishable from a legitimate
    # physics outcome (a baseline that genuinely evaluates to zero), and an
    # inference chain ran to completion on floor values (GH #407).
    #
    # The launch-time check in tidal/cli/_sample.py, _plot_command.py and
    # _analyze.py makes an unresolvable formula unreachable here; this is
    # the backstop, and it is loud.
    return float(safe_formula_eval(formula, ns))  # type: ignore[arg-type]


class SimulationLikelihood:
    """Callable likelihood that runs a simulation per evaluation.

    This wraps the existing sweep infrastructure to produce a single
    scalar log-likelihood from a parameter vector.

    Parameters
    ----------
    base_args : Namespace
        Base simulation arguments (grid, IC, time, etc.).
    spec_path : Path
        Path to the JSON equation specification.
    param_names : list[str]
        Parameter names matching the theta vector ordering.
    measurements : set[str]
        Measurement types to compute.
    source : tuple[str, ...] | None
        Source field(s) for conversion.
    target : tuple[str, ...] | None
        Target field(s) for conversion.
    threshold : float
        Energy conservation threshold.
    likelihood_config : LikelihoodConfig
        How to convert metrics to log-likelihood.
    temp_dir : Path | None
        Base temp directory. If None, uses system temp.
    keep_sims : bool
        If True, keep simulation output directories.
    """

    def __init__(
        self,
        base_args: Namespace,
        spec_path: Path,
        param_names: list[str],
        measurements: set[str],
        source: tuple[str, ...] | None,
        target: tuple[str, ...] | None,
        threshold: float,
        likelihood_config: LikelihoodConfig,
        temp_dir: Path | None = None,
        *,
        keep_sims: bool = False,
    ) -> None:
        self.base_args = base_args
        self.spec_path = spec_path
        self.param_names = param_names
        self.measurements = measurements
        self.source = source
        self.target = target
        self.threshold = threshold
        self.likelihood_config = likelihood_config
        self.temp_dir = temp_dir
        self.keep_sims = keep_sims
        self._call_count = 0
        # Metadata from the most recent __call__: read by the sequential MC
        # path to populate per-sample run_status / tachyonic_excess.  Parallel
        # workers don't share state — they return the metadata from
        # `_likelihood_worker` directly instead.
        self._last_metadata: dict[str, Any] = {}

    def __call__(self, theta: Any) -> float:
        """Evaluate log-likelihood at parameter vector theta.

        Parameters
        ----------
        theta : array-like
            Parameter values in the order of ``param_names``.

        Returns
        -------
        float
            Log-likelihood value. Returns ``-inf`` for failed simulations.
            Per-sample metadata (run_status, tachyonic_excess, etc.) is
            stashed on ``self._last_metadata`` for the caller to read.
        """
        call_idx = self._call_count
        self._call_count += 1
        logl, meta = _evaluate_likelihood(
            theta=theta,
            base_args=self.base_args,
            spec_path=self.spec_path,
            param_names=self.param_names,
            measurements=self.measurements,
            source=self.source,
            target=self.target,
            threshold=self.threshold,
            likelihood_config=self.likelihood_config,
            temp_dir=self.temp_dir,
            keep_sims=self.keep_sims,
            call_index=call_idx,
        )
        self._last_metadata = meta
        return logl


def _evaluate_likelihood(
    *,
    theta: Any,
    base_args: Namespace,
    spec_path: Path,
    param_names: list[str],
    measurements: set[str],
    source: tuple[str, ...] | None,
    target: tuple[str, ...] | None,
    threshold: float,
    likelihood_config: LikelihoodConfig,
    temp_dir: Path | None,
    keep_sims: bool,
    call_index: int,
) -> tuple[float, dict[str, Any]]:
    """Core likelihood evaluation (module-level for pickling).

    This function is the integration point with the sweep infrastructure.

    Returns
    -------
    tuple[float, dict]
        ``(log_likelihood, metadata)`` where ``metadata`` always contains
        a ``run_status`` key with one of:
        ``"success"``, ``"tachyonic_gated"``, ``"simulation_failed"``,
        ``"metric_missing"``, ``"exception"``.  When the rejection is
        tachyonic, ``metadata["tachyonic_excess"]`` carries the maximum
        Re(λ) found (a positive growth rate).  For successful evaluations,
        ``metadata`` may also include the measured metric value.
    """
    from tidal.measurement._run_stages import (
        PointContext,
        RunStatus,
        run_point,
    )

    # Build parameter overrides
    param_overrides = {name: float(theta[i]) for i, name in enumerate(param_names)}

    # Backend selection: in-memory default (fast path), fall back to disk
    # when the user explicitly requests --likelihood-backend=disk.  The
    # disk path is kept for bisectability / rollback (see issue #269).
    backend = getattr(base_args, "likelihood_backend", "memory")

    # The disk backend still needs a temp dir; in-memory does not.
    run_dir: Path | None = None
    if backend == "disk":
        # Create a unique temp directory for this evaluation.
        # Include PID to avoid collisions between multiprocessing workers
        # (each fork gets its own counter starting at 0).
        import os

        base = temp_dir or Path(tempfile.gettempdir())
        run_dir = base / f"inference_run_{os.getpid()}_{call_index:06d}"
        run_dir.mkdir(parents=True, exist_ok=True)

    # Pre-flight stability probe — DIAGNOSTIC ONLY, unconditionally.
    # Records gamma_eff / k_tachyonic / n_tachyonic_modes so chain CSVs
    # support post-hoc filtering.  Nothing acts on the verdict: rejection
    # on tachyonic growth was abandoned as policy (growth cannot be called
    # physics or artifact without theory-level analysis — PSALTer, GH
    # #360) and the ``--gated`` escape hatch was removed in v0.50.0.
    # The probe itself is deliberately conservative:
    # when cond(V) is too high to trust the IC-coupling filter (typical
    # for CDT/PGT), a growing mode counts as tachyonic rather than risk a
    # false negative.  Cost: microseconds per evaluation.
    #
    # The stage lives in tidal.measurement._stability because tidal sweep
    # runs the identical preamble; it used to be copy-pasted here, and the
    # copies drifted into two different policies (GH #454).
    outcome = run_point(
        PointContext(
            spec_path=spec_path,
            base_args=base_args,
            param_overrides=param_overrides,
            measurements=measurements,
            source=source,
            target=target,
            threshold=threshold,
            run_dir=run_dir,
            sampled_params=tuple(param_names) if backend == "memory" else None,
            unavailable_meta=_posdep_probe_unavailable_meta,
        ),
        backend=backend,
    )
    stability_meta = outcome.stability

    # --- policy: outcome -> (logL, metadata) ---------------------------------
    # The soft floor stays HERE, not in run_point: its noise is seeded on
    # theta (GH #408), which makes the likelihood a deterministic function
    # of the parameter point — a correctness property for nested sampling,
    # not a convenience — and run_point has no theta.
    if outcome.status is RunStatus.SIMULATION_FAILED:
        return _soft_floor_logl(
            likelihood_config.soft_floor_noise_sigma,
            likelihood_config.soft_floor_logl,
            _floor_rng(likelihood_config.noise_seed, theta),
        ), {
            "run_status": RunStatus.SIMULATION_FAILED,
            "exit_code": int(outcome.exit_code or 0),
            **stability_meta,
        }

    if outcome.status is RunStatus.METRIC_MISSING:  # pragma: no cover - defensive
        return -math.inf, {
            "run_status": RunStatus.METRIC_MISSING,
            "missing_metric": likelihood_config.metric,
        }

    if not outcome.ok:
        # Divergence, kinetic-evaluation and everything else: all soft-floored
        # identically, so the finer tag cannot move a chain (GH #480).
        return _soft_floor_logl(
            likelihood_config.soft_floor_noise_sigma,
            likelihood_config.soft_floor_logl,
            _floor_rng(likelihood_config.noise_seed, theta),
        ), {"run_status": outcome.status, **stability_meta}

    metrics = outcome.metrics
    spec = outcome.spec

    try:
        # Get the target metric
        metric_value = metrics.get(likelihood_config.metric)
        if metric_value is None:
            # Genuine bug (key missing from sim output) — keep -inf.  Per v3
            # plan, this is distinct from a NaN/Inf metric value (handled
            # below as ``metric_nan`` with the soft floor).
            return -math.inf, {
                "run_status": RunStatus.METRIC_MISSING,
                "missing_metric": likelihood_config.metric,
            }

        # NaN / Inf metric: numerical failure of the sim (overflow at very
        # large growth rates, NaN propagation through the Padé propagator,
        # etc.).  Distinct ``run_status="metric_nan"`` tag (per supervisor
        # request) so post-chain analysis can locate these specifically and
        # inspect what the sim was doing.  Soft floor + noise.
        metric_value_f = float(metric_value)
        if math.isnan(metric_value_f) or math.isinf(metric_value_f):
            return _soft_floor_logl(
                likelihood_config.soft_floor_noise_sigma,
                likelihood_config.soft_floor_logl,
                _floor_rng(likelihood_config.noise_seed, theta),
            ), {
                "run_status": RunStatus.METRIC_NAN,
                likelihood_config.metric: metric_value_f,
                **stability_meta,
            }

        # v3: no Hwang-Noh perturbativity gate.  P_max is recorded faithfully
        # at all values; downstream analysis interprets whether the
        # linearized result is physically meaningful at that parameter point.
        # See ``docs/V3_ARCHITECTURE.md`` and the supervisor meeting notes
        # at ``docs/meetings/2026-05-08_supervisor.md``.

        # Build eval_params for formula-based likelihoods.
        # Merge order (same as simulate_run): spec.metadata["parameters"]
        # defaults -> CLI --param overrides -> swept theta.  Without the
        # first two, formulas like "sin(kappa * B0 * t_end / 2)**2" fail
        # with NameError, _eval_baseline returns None, and every logL
        # collapses to -inf (see #270).
        eval_params: dict[str, float] | None = None
        if likelihood_config.baseline_formula:
            from tidal.measurement._run_stages import parse_params

            eval_params = parse_params(
                list(getattr(base_args, "param", []) or []),
                spec,
            )
            eval_params.update(param_overrides)
            for attr in ("t_end", "dt"):
                val = getattr(base_args, attr, None)
                if val is not None:
                    eval_params[attr] = float(val)

        logl = compute_log_likelihood(
            float(metric_value),
            likelihood_config,
            eval_params,
        )
        success_meta: dict[str, Any] = {
            "run_status": RunStatus.SUCCESS,
            likelihood_config.metric: float(metric_value),
            **stability_meta,
        }
        # In v3 ``compute_log_likelihood`` returns finite values for almost
        # every numeric input — the only ``-inf`` paths left are: NaN/Inf
        # metric (already routed to ``metric_nan`` soft floor upstream),
        # threshold-likelihood below ``min_value``, gaussian/extremize edge
        # cases, and ``baseline <= 0 or metric_value <= 0``.  The first is
        # by construction handled before this point; the others are rare.
        # Route any residual ``-inf`` to the soft floor with a distinct
        # ``logl_minus_inf`` tag so post-chain analysis can locate them.
        if math.isinf(logl) and logl < 0:
            return _soft_floor_logl(
                likelihood_config.soft_floor_noise_sigma,
                likelihood_config.soft_floor_logl,
                _floor_rng(likelihood_config.noise_seed, theta),
            ), {
                **success_meta,
                "run_status": RunStatus.LOGL_MINUS_INF,
            }
        # Distinct tag for samples whose simulation succeeded but logL fell
        # below the soft-floor sentinel — typically P_max ~ 1e-44 or smaller
        # in maximize mode, which is below the simulation's effective
        # numerical noise floor (issue #356).  Value is preserved (no
        # clamping); only the run_status tag changes so post-chain analysis
        # can filter these from "physical" min/max summaries.
        if math.isfinite(logl) and logl < SOFT_FLOOR_LOGL:
            return logl, {**success_meta, "run_status": RunStatus.BELOW_NOISE_FLOOR}
        return logl, success_meta  # noqa: TRY300

    except Exception:  # noqa: BLE001
        # Failures in the POLICY stage only — metric extraction, the
        # likelihood formula.  Everything that can go wrong while probing,
        # simulating or measuring is caught and classified by run_point,
        # which is why the SimulationDivergedError / KineticEvaluationError
        # handlers that used to sit here are gone: those exceptions can no
        # longer reach this frame.
        import logging

        logging.getLogger("tidal.inference").debug(
            "Likelihood evaluation failed at theta=%s",
            theta,
            exc_info=True,
        )
        # v3 soft floor: never -inf for numerical exception.  Sampler can
        # still reject these in practice (exp(-100) ~ 4e-44) but sees a
        # gradient in the failure region.
        return _soft_floor_logl(
            likelihood_config.soft_floor_noise_sigma,
            likelihood_config.soft_floor_logl,
            _floor_rng(likelihood_config.noise_seed, theta),
        ), {"run_status": RunStatus.EXCEPTION}

    finally:
        if run_dir is not None and not keep_sims and run_dir.exists():
            shutil.rmtree(run_dir, ignore_errors=True)


# Module-level wrapper for multiprocessing Pool.map
# Follows the same pattern as _run_single_wrapper in _sweep.py
_LIKELIHOOD_CONFIG: dict[str, Any] = {}


def _likelihood_worker_init(config: dict[str, Any]) -> None:  # pyright: ignore[reportUnusedFunction]
    """Initialize worker process with likelihood config."""
    from tidal.measurement._run_stages import set_single_thread_blas

    set_single_thread_blas()
    _LIKELIHOOD_CONFIG.update(config)


def _likelihood_worker(theta: Any) -> tuple[float, dict[str, Any]]:  # pyright: ignore[reportUnusedFunction]
    """Evaluate likelihood in worker process (picklable module-level function).

    Returns ``(log_likelihood, metadata)``.  The metadata dict carries the
    per-sample ``run_status`` and rejection details (e.g.
    ``tachyonic_excess``) so the parent process can populate per-sample
    columns in the result table.
    """
    import threading

    # Thread-safe call counter
    counter_lock = threading.Lock()
    with counter_lock:
        idx = _LIKELIHOOD_CONFIG.get("_counter", 0)
        _LIKELIHOOD_CONFIG["_counter"] = idx + 1

    return _evaluate_likelihood(
        theta=theta,
        base_args=_LIKELIHOOD_CONFIG["base_args"],
        spec_path=Path(_LIKELIHOOD_CONFIG["spec_path"]),
        param_names=_LIKELIHOOD_CONFIG["param_names"],
        measurements=_LIKELIHOOD_CONFIG["measurements"],
        source=_LIKELIHOOD_CONFIG["source"],
        target=_LIKELIHOOD_CONFIG["target"],
        threshold=_LIKELIHOOD_CONFIG["threshold"],
        likelihood_config=_LIKELIHOOD_CONFIG["likelihood_config"],
        temp_dir=Path(_LIKELIHOOD_CONFIG["temp_dir"])
        if _LIKELIHOOD_CONFIG.get("temp_dir")
        else None,
        keep_sims=_LIKELIHOOD_CONFIG.get("keep_sims", False),
        call_index=idx,
    )
