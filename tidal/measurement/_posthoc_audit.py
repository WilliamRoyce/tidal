"""Post-hoc t-independence audit of inference chains.

Re-runs simulations for a stratified sample of an inference chain at
``2 × t_base`` and tests amplitude-ratio time-stationarity. The
amplification ratio ``A(t) = P_theory(t) / P_GR(t)`` is what should be
approximately stationary for genuine new physics — *not* ``P(t)``,
which Rabi-oscillates. ``A(t)`` grows exponentially when the system
hosts a slow ghost that the pre-flight probe failed to catch (#323's
false-negative regime).

Architecture
------------

The audit complements the pre-flight stability probe; it does not
replace it. The probe runs in microseconds per evaluation and is
campaign-fixed (post-#341) at ``threshold=0.15`` for cross-stage
``D_KL`` comparability — chains run before #341 used ``threshold=0.3``
and have ``stability_profile`` column ``unit-ic-all-k-0.3``. The audit
runs ~3 s per sample and applies a stricter physics gate
(``P_max < 0.5`` Hwang–Noh perturbativity bound) that the probe cannot
encode.

Per-sample classification
-------------------------

Single hard physics gate; everything else is *flag for review* (never
auto-rejected — silent rejection of borderline samples would bias the
posterior in regions where new physics may legitimately modify the
amplification ratio over time).

* **hard_reject** — ``P_max_total > perturbativity_p_max``. Linearized
  description has broken down regardless of mechanism. Auto-zero the
  sample's weight in the reweighted chain.
* **flag_review** — at least one of:

  - ``|log(A_med_2 / A_med_1)| > log(3)`` (exponential proxy: catches
    ``γ_eff_posthoc ≳ 0.11/s`` at ``t_base = 10``).
  - ``P_max_2 / max(P_max_1, ε) > 5`` (independent slow-growth
    signature).
  - ``borderline_stability == 1`` from the original probe.

  Surfaced in ``samples.audit.csv``; weight unchanged in
  ``samples.reweighted.csv``. Human triage decides per stage.
* **pass** — neither gate fires. Weight unchanged.

Output schema
-------------

* ``samples.audit.csv`` — one row per audited sample with full
  classification metadata.
* ``samples.reweighted.csv`` — full chain with the reweighted ``weight``
  column; ``hard_reject`` rows zeroed, others unchanged. Bayesian-
  correct: preserves ``nlive`` accounting; filtering rows would not.
* ``audit_summary.json`` — counts, contamination rate, worst-case
  parameters and ``γ_eff_posthoc``, profile name, t_base used.

References
----------
- Issue #323 (false-negative regime, audit recipe).
- ``docs/tex/troubleshooting.tex`` (post-#323 update).
"""

from __future__ import annotations

import contextlib
import dataclasses
import json
import math
import operator
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

if TYPE_CHECKING:
    from argparse import Namespace
    from collections.abc import Sequence
    from pathlib import Path

    from numpy.typing import NDArray

    from tidal.inference._results import InferenceResult


AuditStatus = Literal["pass", "flag_review", "hard_reject", "sim_failed"]


@dataclasses.dataclass(frozen=True)
class SampleAudit:
    """Audit verdict for one chain sample."""

    index: int
    parameters: dict[str, float]
    p_max_first: float
    p_max_second: float
    p_max_total: float
    a_med_first: float
    a_med_second: float
    a_max_first: float
    a_max_second: float
    a_med_log_ratio: float
    gamma_eff_posthoc: float
    status: AuditStatus
    reason: str
    selection_strategy: str


@dataclasses.dataclass(frozen=True)
class AuditSummary:
    """Aggregate report on an audit run."""

    n_samples_total: int
    n_samples_audited: int
    n_pass: int
    n_flag_review: int
    n_hard_reject: int
    n_sim_failed: int
    contamination_rate: float
    worst_sample_index: int | None
    worst_gamma_eff_posthoc: float
    profile_name: str
    t_base: float
    perturbativity_p_max: float
    selection_strategies_used: tuple[str, ...]


def select_audit_samples(
    result: InferenceResult,
    *,
    n_top: int = 50,
    n_stratified: int = 200,
    n_clusters: int = 5,
    include_map: bool = True,
    rng_seed: int = 0,
) -> dict[int, str]:
    """Choose chain indices to audit.

    Returns a mapping ``{index → strategy_name}``. A sample appearing
    under multiple strategies is recorded under the first strategy
    that selected it (priority order: ``borderline`` > ``top`` > ``map``
    > ``stratified``).

    Strategies
    ----------
    * **borderline** — every sample with ``borderline_stability == 1``
      from the original probe. Highest priority because this is the
      empirical false-negative-prone strip the audit was designed to
      catch.
    * **top** — top ``n_top`` samples by ``log_posterior``. The
      "headline" samples that drive any downstream MAP / corner plot.
    * **map** — the single MAP point if not already covered.
    * **stratified** — k-means on the parameter space; sample
      ``n_stratified // n_clusters`` from each cluster, weighted by
      importance weight (or uniform if no weights). Bias-free coverage
      across multimodal posteriors.
    """
    selected: dict[int, str] = {}

    # 1. Borderline cohort — exhaustive coverage of the false-negative
    #    strip the probe is known to miss (#323).
    metrics = result.metrics or {}
    borderline_arr = metrics.get("borderline_stability")
    if borderline_arr is not None:
        for i, flag in enumerate(borderline_arr):
            try:
                if int(float(flag)) == 1 and i not in selected:
                    selected[i] = "borderline"
            except (ValueError, TypeError):
                continue

    # 2. Top-posterior cohort.
    log_post = result.log_likelihood + result.log_prior
    n_total = len(log_post)
    n_top_eff = min(n_top, n_total)
    if n_top_eff > 0:
        top_idx = np.argsort(log_post)[-n_top_eff:][::-1]
        for i in top_idx:
            i_int = int(i)
            if i_int not in selected:
                selected[i_int] = "top"

    # 3. MAP — single best sample, guaranteed in the audit.
    if include_map and n_total > 0:
        map_idx = int(np.argmax(log_post))
        if map_idx not in selected:
            selected[map_idx] = "map"

    # 4. Stratified random sample by k-means clustering.
    if n_stratified > 0 and n_total > 1:
        rng = np.random.default_rng(rng_seed)
        try:
            from sklearn.cluster import KMeans  # type: ignore[import-untyped]

            k = max(1, min(n_clusters, n_total))
            km = KMeans(n_clusters=k, random_state=rng_seed, n_init=10).fit(  # type: ignore[no-untyped-call]
                result.samples,
            )
            labels = np.asarray(km.labels_, dtype=int)  # type: ignore[reportUnknownArgumentType]
        except ImportError:
            # Fall back to single cluster (uniform random sample).
            labels = np.zeros(n_total, dtype=int)

        per_cluster = max(1, n_stratified // max(1, len(np.unique(labels))))
        for c in np.unique(labels):
            cluster_idx = np.where(labels == c)[0]
            if len(cluster_idx) == 0:
                continue
            weights = (
                result.weights[cluster_idx] if result.weights is not None else None
            )
            if weights is not None and float(np.sum(weights)) > 0:
                p = np.asarray(weights, dtype=float)
                p /= float(np.sum(p))
            else:
                p = None
            n_pick = min(per_cluster, len(cluster_idx))
            picked = rng.choice(cluster_idx, size=n_pick, replace=False, p=p)
            for i in picked:
                i_int = int(i)
                if i_int not in selected:
                    selected[i_int] = "stratified"

    return selected


def _classify_sample(
    p_max_first: float,
    p_max_second: float,
    a_med_first: float,
    a_med_second: float,
    *,
    perturbativity_p_max: float,
    a_log_ratio_threshold: float,
    p_growth_ratio_threshold: float,
    borderline_flag: bool,
) -> tuple[AuditStatus, str]:
    """Apply the hard-gate + flag-review classification.

    See module docstring for the policy. Returns ``(status, reason)``.
    """
    p_max_total = max(p_max_first, p_max_second)
    if not math.isfinite(p_max_total) or p_max_total > perturbativity_p_max:
        return "hard_reject", (
            f"P_max_total={p_max_total:.4g} > {perturbativity_p_max} (Hwang-Noh)"
        )

    reasons: list[str] = []

    if (
        a_med_first > 0
        and a_med_second > 0
        and abs(math.log(a_med_second / a_med_first)) > a_log_ratio_threshold
    ):
        reasons.append(
            f"|log(A_med_2/A_med_1)|={abs(math.log(a_med_second / a_med_first)):.3f} "
            f"> {a_log_ratio_threshold:.3f}",
        )

    p_growth = p_max_second / max(p_max_first, 1e-12)
    if p_growth > p_growth_ratio_threshold:
        reasons.append(
            f"P_max_2/P_max_1={p_growth:.2f} > {p_growth_ratio_threshold:.1f}",
        )

    if borderline_flag:
        reasons.append("probe-borderline")

    if reasons:
        return "flag_review", "; ".join(reasons)

    return "pass", "perturbative; A_med stationary"


# Defaults the simulate code path expects directly (not via getattr).  The
# analyze CLI defines a small subset; the audit re-uses the same simulate
# entry point as the inference likelihood, so any missing attribute would
# raise AttributeError on the first chained sim.  Keep keys aligned with
# _simulate.py's `args.<x>` direct accesses, not the broader getattr set.
def _backfill_simulate_args(audit_args: Namespace) -> None:
    """Fill in simulate-required attrs missing from the analyze namespace.

    The analyze CLI was scoped to chain analysis, so it defines only the
    simulate flags the audit explicitly needs (--ic, --ic-component, etc.).
    The simulate code path direct-accesses a wider set of attributes —
    backfill those with their canonical simulate defaults to avoid
    AttributeError on every sample.
    """
    from tidal.solver._defaults import DEFAULT_ATOL, DEFAULT_RTOL

    defaults: dict[str, Any] = {
        "bc": None,
        # _resolve_scheme returns the value as-is when != "auto", so a
        # missing scheme would propagate None into the solver dispatch.
        # "auto" triggers the standard solver selection used by the
        # campaign chains.
        "scheme": "auto",
        "mode": None,
        "method": None,
        "dt": None,
        # rtol/atol are direct-accessed by _validate_solver_params with
        # `> 0` checks (no None handling), so they need numeric defaults.
        "atol": DEFAULT_ATOL,
        "rtol": DEFAULT_RTOL,
        "max_step": None,
        "ic_width": None,
        "ic_center": None,
        "ic_formula": None,
        "output": None,
        "output_format": None,
        "no_plot": True,
        "quiet": True,
        "resume": None,
        "snapshot": None,
        "t_additional": None,
        "json_path": None,
    }
    for key, default in defaults.items():
        if not hasattr(audit_args, key):
            setattr(audit_args, key, default)


def _evaluate_sample(
    sample_index: int,
    parameters: dict[str, float],
    *,
    base_args: Namespace,
    spec_path: Path,
    source: tuple[str, ...],
    target: tuple[str, ...],
    t_base: float,
    perturbativity_p_max: float,
    baseline_formula: str | None,
    borderline_flag: bool,
    selection_strategy: str,
    a_log_ratio_threshold: float = math.log(3.0),
    p_growth_ratio_threshold: float = 5.0,
) -> SampleAudit:
    """Run one audit simulation at ``2·t_base`` and classify."""
    import copy

    from tidal.measurement._conversion import compute_conversion_probability
    from tidal.measurement._run_stages import run_inference_step

    audit_args = copy.copy(base_args)
    audit_args.t_end = 2.0 * t_base
    # The analyze CLI's --snapshots is an integer count for clarity (the
    # audit's intent is "N snapshots covering both halves"), but simulate
    # interprets --snapshots as a float time interval.  Convert count → dt
    # here so the t-axis is uniform across [0, 2*t_base].
    snapshot_count = getattr(audit_args, "snapshots", None)
    if snapshot_count is not None and snapshot_count > 1:
        audit_args.snapshots = float(audit_args.t_end) / float(snapshot_count - 1)
    else:
        audit_args.snapshots = None
    # The analyze CLI doesn't define every simulate flag, so back-fill any
    # _simulate-required attributes that may be absent on the analyze
    # namespace.  Defaults match _build_parser's "simulate" defaults.
    _backfill_simulate_args(audit_args)

    from tidal.solver._exceptions import SimulationDivergedError

    try:
        sim_data = run_inference_step(audit_args, spec_path, parameters)
    except SimulationDivergedError as exc:
        # The simulate code path's own Padé probe detected non-perturbative
        # divergence at 2*t_base — exactly the contamination this audit
        # exists to surface.  Classify as a hard reject (consistent with
        # the Hwang-Noh gate's intent) rather than as a generic sim failure.
        return SampleAudit(
            index=sample_index,
            parameters=parameters,
            p_max_first=float("nan"),
            p_max_second=float("nan"),
            p_max_total=float("nan"),
            a_med_first=float("nan"),
            a_med_second=float("nan"),
            a_max_first=float("nan"),
            a_max_second=float("nan"),
            a_med_log_ratio=float("nan"),
            gamma_eff_posthoc=float("inf"),
            status="hard_reject",
            reason=f"non-perturbative @ 2*t_base: {exc}",
            selection_strategy=selection_strategy,
        )
    except Exception as exc:  # noqa: BLE001
        return SampleAudit(
            index=sample_index,
            parameters=parameters,
            p_max_first=float("nan"),
            p_max_second=float("nan"),
            p_max_total=float("nan"),
            a_med_first=float("nan"),
            a_med_second=float("nan"),
            a_max_first=float("nan"),
            a_max_second=float("nan"),
            a_med_log_ratio=float("nan"),
            gamma_eff_posthoc=float("nan"),
            status="sim_failed",
            reason=f"{type(exc).__name__}: {exc}",
            selection_strategy=selection_strategy,
        )

    conv = compute_conversion_probability(sim_data, source[0], target[0])
    times = np.asarray(conv.times, dtype=float)
    p_t = np.asarray(conv.probability, dtype=float)

    # Window split. Use t_base as the boundary; first window has at
    # least one snapshot if the snapshot grid covers it.
    mask_first = times <= t_base + 1e-9
    mask_second = times > t_base - 1e-9
    if not np.any(mask_first) or not np.any(mask_second):
        return SampleAudit(
            index=sample_index,
            parameters=parameters,
            p_max_first=float("nan"),
            p_max_second=float("nan"),
            p_max_total=float("nan"),
            a_med_first=float("nan"),
            a_med_second=float("nan"),
            a_max_first=float("nan"),
            a_max_second=float("nan"),
            a_med_log_ratio=float("nan"),
            gamma_eff_posthoc=float("nan"),
            status="sim_failed",
            reason=(
                f"Snapshot grid does not cover both halves "
                f"(t_base={t_base}, snapshots t∈[{times.min()},{times.max()}])"
            ),
            selection_strategy=selection_strategy,
        )

    p_max_first = float(np.max(p_t[mask_first]))
    p_max_second = float(np.max(p_t[mask_second]))
    p_max_total = max(p_max_first, p_max_second)

    # Time-resolved baseline. P_GR(t) ≈ sin²(κB₀·t/2) for the
    # canonical Gertsenshtein channel; the formula is provided by the
    # caller via base_args.baseline_formula or extracted from the
    # likelihood config.
    if baseline_formula is not None:
        p_gr_t = _eval_baseline_at_times(baseline_formula, parameters, times)
    else:
        # No baseline → fall back to A = P (degenerate) so the audit
        # still reports the P-based stats.
        p_gr_t = np.maximum(p_t, 1e-12)

    p_gr_floor = 1e-12
    a_t = p_t / np.maximum(p_gr_t, p_gr_floor)

    a_first = a_t[mask_first]
    a_second = a_t[mask_second]
    a_med_first = float(np.median(a_first)) if a_first.size else 0.0
    a_med_second = float(np.median(a_second)) if a_second.size else 0.0
    a_max_first = float(np.max(a_first)) if a_first.size else 0.0
    a_max_second = float(np.max(a_second)) if a_second.size else 0.0

    a_med_log_ratio = (
        math.log(a_med_second / a_med_first)
        if a_med_first > 0 and a_med_second > 0
        else 0.0
    )
    gamma_eff_posthoc = (
        math.log(p_max_second / max(p_max_first, 1e-12)) / t_base
        if p_max_first > 0
        else 0.0
    )

    status, reason = _classify_sample(
        p_max_first,
        p_max_second,
        a_med_first,
        a_med_second,
        perturbativity_p_max=perturbativity_p_max,
        a_log_ratio_threshold=a_log_ratio_threshold,
        p_growth_ratio_threshold=p_growth_ratio_threshold,
        borderline_flag=borderline_flag,
    )

    return SampleAudit(
        index=sample_index,
        parameters=parameters,
        p_max_first=p_max_first,
        p_max_second=p_max_second,
        p_max_total=p_max_total,
        a_med_first=a_med_first,
        a_med_second=a_med_second,
        a_max_first=a_max_first,
        a_max_second=a_max_second,
        a_med_log_ratio=a_med_log_ratio,
        gamma_eff_posthoc=gamma_eff_posthoc,
        status=status,
        reason=reason,
        selection_strategy=selection_strategy,
    )


def _eval_baseline_at_times(
    formula: str,
    parameters: dict[str, float],
    times: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Evaluate the closed-form Gertsenshtein baseline P_GR(t) elementwise.

    Reuses the same name-space as ``_eval_baseline`` in
    ``tidal.inference._likelihood`` (sin/cos/sqrt/pi/...) and adds the
    free variable ``t_end`` per snapshot. The baseline formula is
    typically ``"sin(kappa*B0*t_end/2)**2"``; substituting each
    snapshot's time gives the time-resolved baseline.
    """
    from tidal.cli._simulate import (
        FORMULA_NAMESPACE,  # pyright: ignore[reportPrivateUsage]
    )

    out = np.zeros_like(times)
    for i, t in enumerate(times):
        ns: dict[str, object] = {**FORMULA_NAMESPACE, **parameters, "t_end": float(t)}
        try:
            from tidal.cli._simulate import safe_formula_eval

            val = float(safe_formula_eval(formula, ns))  # type: ignore[arg-type]
        except Exception:  # noqa: BLE001
            val = float("nan")
        out[i] = val if math.isfinite(val) else float("nan")
    return out


def run_posthoc_audit(
    inference_dir: Path,
    *,
    base_args: Namespace,
    spec_path: Path,
    source: tuple[str, ...],
    target: tuple[str, ...],
    profile_name: str,
    t_base: float,
    baseline_formula: str | None = None,
    output_dir: Path | None = None,
    n_top: int = 50,
    n_stratified: int = 200,
    n_clusters: int = 5,
    perturbativity_p_max: float = 0.5,
    rng_seed: int = 0,
) -> AuditSummary:
    """Audit an inference chain at ``2·t_base`` and write artifacts.

    Loads ``inference.json`` + ``results.csv`` from *inference_dir*,
    selects samples per :func:`select_audit_samples`, re-runs each at
    ``2·t_base``, classifies, and writes ``samples.audit.csv``,
    ``samples.reweighted.csv``, and ``audit_summary.json`` to
    *output_dir* (default: ``inference_dir / "audit"``).

    The original chain is left untouched.
    """
    from tidal.inference._results import InferenceResult

    if output_dir is None:
        output_dir = inference_dir / "audit"
    output_dir.mkdir(parents=True, exist_ok=True)

    result = InferenceResult.from_directory(inference_dir)
    selected = select_audit_samples(
        result,
        n_top=n_top,
        n_stratified=n_stratified,
        n_clusters=n_clusters,
        rng_seed=rng_seed,
    )

    audits: list[SampleAudit] = []
    metrics = result.metrics or {}
    borderline_arr = metrics.get("borderline_stability")
    for index in sorted(selected.keys()):
        params = {
            name: float(result.samples[index, j])
            for j, name in enumerate(result.param_names)
        }
        # Merge in any --param fixed values from base_args so the audit
        # sim runs with the same fixed parameters as the original.
        from tidal.measurement._run_stages import parse_params
        from tidal.symbolic.json_loader import load_equation_system

        spec = load_equation_system(spec_path)
        base_p = parse_params(list(getattr(base_args, "param", []) or []), spec)
        full_params = {**base_p, **params}
        borderline_flag = False
        if borderline_arr is not None:
            try:
                borderline_flag = int(float(borderline_arr[index])) == 1
            except (ValueError, TypeError):
                borderline_flag = False
        audit = _evaluate_sample(
            index,
            full_params,
            base_args=base_args,
            spec_path=spec_path,
            source=source,
            target=target,
            t_base=t_base,
            perturbativity_p_max=perturbativity_p_max,
            baseline_formula=baseline_formula,
            borderline_flag=borderline_flag,
            selection_strategy=selected[index],
        )
        audits.append(audit)

    _write_audit_csv(audits, output_dir / "samples.audit.csv", result.param_names)
    _write_reweighted_csv(
        result,
        audits,
        output_dir / "samples.reweighted.csv",
    )

    summary = _build_summary(
        result,
        audits,
        profile_name=profile_name,
        t_base=t_base,
        perturbativity_p_max=perturbativity_p_max,
        selection_strategies=tuple(sorted(set(selected.values()))),
    )
    (output_dir / "audit_summary.json").write_text(
        json.dumps(dataclasses.asdict(summary), indent=2),
        encoding="utf-8",
    )
    return summary


def _write_audit_csv(
    audits: Sequence[SampleAudit],
    path: Path,
    param_names: Sequence[str],
) -> None:
    import csv

    fields = [
        "sample_index",
        *param_names,
        "p_max_first",
        "p_max_second",
        "p_max_total",
        "a_med_first",
        "a_med_second",
        "a_max_first",
        "a_max_second",
        "a_med_log_ratio",
        "gamma_eff_posthoc",
        "audit_status",
        "audit_reason",
        "selection_strategy",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        for a in audits:
            row: list[Any] = [a.index]
            row.extend(a.parameters.get(name, "") for name in param_names)
            row.extend(
                [
                    a.p_max_first,
                    a.p_max_second,
                    a.p_max_total,
                    a.a_med_first,
                    a.a_med_second,
                    a.a_max_first,
                    a.a_max_second,
                    a.a_med_log_ratio,
                    a.gamma_eff_posthoc,
                    a.status,
                    a.reason,
                    a.selection_strategy,
                ],
            )
            writer.writerow(row)


def _write_reweighted_csv(
    result: InferenceResult,
    audits: Sequence[SampleAudit],
    path: Path,
) -> None:
    """Mirror the original ``results.csv`` with a reweighted ``weight``.

    For ``hard_reject`` audited rows: ``weight *= 0`` (nested-sampling
    accounting preserved). All other rows pass through unchanged. Adds
    an ``audit_status`` column so downstream consumers can find the
    flagged rows without re-loading the audit CSV.
    """
    import csv

    sr = result.to_sweep_results()
    rows = list(sr.rows)
    audit_by_idx = {a.index: a for a in audits}
    for i, row in enumerate(rows):
        a = audit_by_idx.get(i)
        if a is None:
            row["audit_status"] = "not_evaluated"
            continue
        row["audit_status"] = a.status
        if a.status == "hard_reject" and "weight" in row:
            with contextlib.suppress(TypeError, KeyError):
                row["weight"] = 0.0

    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for k in row:
            if k not in seen:
                fieldnames.append(k)
                seen.add(k)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _build_summary(
    result: InferenceResult,
    audits: Sequence[SampleAudit],
    *,
    profile_name: str,
    t_base: float,
    perturbativity_p_max: float,
    selection_strategies: tuple[str, ...],
) -> AuditSummary:
    n_pass = sum(1 for a in audits if a.status == "pass")
    n_flag = sum(1 for a in audits if a.status == "flag_review")
    n_reject = sum(1 for a in audits if a.status == "hard_reject")
    n_failed = sum(1 for a in audits if a.status == "sim_failed")

    finite_gammas = [
        (a.index, a.gamma_eff_posthoc)
        for a in audits
        if math.isfinite(a.gamma_eff_posthoc)
    ]
    if finite_gammas:
        worst_index, worst_gamma = max(finite_gammas, key=operator.itemgetter(1))
    else:
        worst_index, worst_gamma = None, 0.0

    contamination = (n_reject + n_flag) / max(1, len(audits))

    return AuditSummary(
        n_samples_total=int(result.n_samples),
        n_samples_audited=len(audits),
        n_pass=n_pass,
        n_flag_review=n_flag,
        n_hard_reject=n_reject,
        n_sim_failed=n_failed,
        contamination_rate=float(contamination),
        worst_sample_index=worst_index,
        worst_gamma_eff_posthoc=float(worst_gamma),
        profile_name=profile_name,
        t_base=float(t_base),
        perturbativity_p_max=float(perturbativity_p_max),
        selection_strategies_used=selection_strategies,
    )
