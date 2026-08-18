"""Critical field analysis — find minimum field strength for target conversion.

Given sweep results where one swept parameter is a field-strength variable
(e.g. B0, Bpeak), this module finds the **minimum** value of that parameter
at which a conversion metric (typically ``P_final``) first crosses a threshold.

This enables the key comparison for Gertsenshtein characterization:

    amplification = B_min(E-M) / B_min(new theory)

where both are evaluated at the same threshold (full conversion, P ≈ 1).

The output is a reduced ``SweepResults`` with the field-strength parameter
collapsed, containing ``B_min``, ``inv_B_min``, and associated errors as
metric columns.  This plugs directly into the existing ``tidal plot`` heatmap
infrastructure.

References
----------
Boccaletti, D. et al. (1970) Nuovo Cim. 70B, 129-146 (graviton-photon
conversion in a static magnetic field).

Dandoy, V. et al. (2024) arXiv:2406.17853 (gravitational wave constraints
from Gertsenshtein effect).
"""

from __future__ import annotations

import math
import operator
from dataclasses import dataclass
from typing import Any

import numpy as np

from tidal.measurement._sweep_results import SweepResults

__all__ = [
    "CriticalFieldResult",
    "compute_critical_field",
    "compute_reference_threshold",
    "critical_field_to_sweep_results",
]


# ── Crossing quality flags ──────────────────────────────────────────


QUALITY_GOOD = "good"
QUALITY_COARSE = "coarse"
QUALITY_EDGE = "edge"
QUALITY_NONE = "none"


# ── Result dataclass ────────────────────────────────────────────────


@dataclass(frozen=True)
class CriticalFieldResult:
    """Result of critical field extraction.

    Attributes
    ----------
    rows : list[dict[str, Any]]
        One dict per outer-param combination.  Keys include outer param
        values plus ``B_min``, ``inv_B_min``, error columns, and quality.
    field_param : str
        Name of the field-strength parameter that was collapsed.
    threshold : float
        Metric threshold used for the crossing.
    metric : str
        Metric name used for thresholding.
    outer_params : dict[str, list[float]]
        Remaining swept parameters after collapsing *field_param*.
    """

    rows: list[dict[str, Any]]
    field_param: str
    threshold: float
    metric: str
    outer_params: dict[str, list[float]]


# ── Core algorithm ──────────────────────────────────────────────────


def compute_critical_field(
    results: SweepResults,
    field_param: str,
    metric: str = "P_final",
    threshold: float = 0.99,
    *,
    interpolate: bool = True,
) -> CriticalFieldResult:
    """Find minimum field strength for a metric to cross a threshold.

    For each unique combination of the "outer" swept parameters (everything
    except *field_param*), the rows are sorted by *field_param* and scanned
    upward.  The first crossing of ``metric >= threshold`` defines ``B_min``.

    Parameters
    ----------
    results : SweepResults
        Sweep data with *field_param* as one of the swept parameters.
    field_param : str
        The field-strength parameter to threshold on (e.g. ``"B0"``).
    metric : str
        Metric column to compare against *threshold* (default ``"P_final"``).
    threshold : float
        Target value for *metric* (default 0.99 — full conversion).
    interpolate : bool
        If True, linearly interpolate between bracketing grid points for
        sub-grid accuracy.

    Returns
    -------
    CriticalFieldResult

    Raises
    ------
    ValueError
        If *field_param* is not a swept parameter or *metric* is not found.
    """
    if field_param not in results.swept_params:
        msg = (
            f"'{field_param}' is not a swept parameter. "
            f"Available: {list(results.swept_params.keys())}"
        )
        raise ValueError(msg)

    # Check metric exists in at least one row
    if results.rows and not any(metric in row for row in results.rows):
        msg = (
            f"Metric '{metric}' not found in sweep results. "
            f"Available: {results.metric_names}"
        )
        raise ValueError(msg)

    # Identify outer params
    outer_param_names = [p for p in results.swept_params if p != field_param]

    # Group rows by outer param combination
    groups = _group_by_outer(results.rows, outer_param_names)

    # Process each group
    output_rows: list[dict[str, Any]] = []
    outer_param_values: dict[str, set[float]] = {p: set() for p in outer_param_names}

    for outer_key, group_rows in sorted(groups.items()):
        row = _process_group(
            group_rows,
            field_param,
            metric,
            threshold,
            outer_param_names,
            outer_key,
            interpolate=interpolate,
        )
        output_rows.append(row)
        for i, p in enumerate(outer_param_names):
            outer_param_values[p].add(outer_key[i])

    outer_params = {p: sorted(outer_param_values[p]) for p in outer_param_names}

    return CriticalFieldResult(
        rows=output_rows,
        field_param=field_param,
        threshold=threshold,
        metric=metric,
        outer_params=outer_params,
    )


def _group_by_outer(
    rows: list[dict[str, Any]],
    outer_param_names: list[str],
) -> dict[tuple[float, ...], list[dict[str, Any]]]:
    """Group rows by unique outer parameter combinations."""
    groups: dict[tuple[float, ...], list[dict[str, Any]]] = {}
    for row in rows:
        key = tuple(float(row.get(p, float("nan"))) for p in outer_param_names)
        groups.setdefault(key, []).append(row)
    return groups


def _process_group(  # noqa: PLR0913, PLR0914, PLR0917
    group_rows: list[dict[str, Any]],
    field_param: str,
    metric: str,
    threshold: float,
    outer_param_names: list[str],
    outer_key: tuple[float, ...],
    *,
    interpolate: bool,
) -> dict[str, Any]:
    """Extract B_min for a single outer-parameter combination."""
    # Build output row with outer param values
    row: dict[str, Any] = {}
    for i, p in enumerate(outer_param_names):
        row[p] = outer_key[i]

    # Sort by field param, filtering out rows with missing metric
    valid = [
        (float(r[field_param]), float(r[metric]))
        for r in group_rows
        if r.get(field_param) is not None and r.get(metric) is not None
    ]
    if not valid:
        return _fill_nan(row)

    valid.sort(key=operator.itemgetter(0))
    b_vals = np.array([v[0] for v in valid])
    m_vals = np.array([v[1] for v in valid])

    # Find first crossing
    crossing_idx = _find_first_crossing(m_vals, threshold)

    if crossing_idx is None:
        # No crossing found
        return _fill_nan(row, quality=QUALITY_NONE)

    if crossing_idx == 0:
        # Already exceeded at smallest value
        b_min = float(b_vals[0])
        err_grid = b_min / 2.0
        row.update(
            B_min=b_min,
            inv_B_min=1.0 / b_min if b_min > 0 else float("inf"),
            B_min_err=err_grid,
            B_min_err_grid=err_grid,
            B_min_err_metric=float("nan"),
            B_min_err_interp=float("nan"),
            inv_B_min_err=err_grid / b_min**2 if b_min > 0 else float("nan"),
            crossing_quality=QUALITY_EDGE,
        )
        return row

    # Bracketing points
    b_lo, b_hi = float(b_vals[crossing_idx - 1]), float(b_vals[crossing_idx])
    m_lo, m_hi = float(m_vals[crossing_idx - 1]), float(m_vals[crossing_idx])
    delta_b = b_hi - b_lo

    # Linear interpolation
    if interpolate and (m_hi - m_lo) > 0:
        frac = (threshold - m_lo) / (m_hi - m_lo)
        b_min = b_lo + frac * delta_b
    else:
        b_min = b_hi

    # Error: grid spacing
    err_grid = delta_b / 2.0

    # Error: interpolation model (quadratic vs linear)
    err_interp = _interpolation_model_error(
        b_vals,
        m_vals,
        crossing_idx,
        threshold,
        b_min,
    )

    # Combined error (metric error requires replicate data, not available here)
    err_metric = float("nan")
    err_combined = math.sqrt(
        err_grid**2 + (0.0 if math.isnan(err_interp) else err_interp**2),
    )

    # Quality flag
    quality = QUALITY_COARSE if err_grid > 0.1 * abs(b_min) else QUALITY_GOOD

    # Error propagation for 1/B_min: err_{1/B} = err_B / B^2
    inv_b_min = 1.0 / b_min if b_min > 0 else float("inf")
    inv_b_err = err_combined / b_min**2 if b_min > 0 else float("nan")

    row.update(
        B_min=b_min,
        inv_B_min=inv_b_min,
        B_min_err=err_combined,
        B_min_err_grid=err_grid,
        B_min_err_metric=err_metric,
        B_min_err_interp=err_interp,
        inv_B_min_err=inv_b_err,
        crossing_quality=quality,
    )
    return row


def _find_first_crossing(m_vals: np.ndarray, threshold: float) -> int | None:
    """Find the index of the first value >= threshold.

    Returns None if no crossing exists.
    """
    above = np.where(m_vals >= threshold)[0]
    if len(above) == 0:
        return None
    return int(above[0])


def _interpolation_model_error(
    b_vals: np.ndarray,
    m_vals: np.ndarray,
    crossing_idx: int,
    threshold: float,
    b_min_linear: float,
) -> float:
    """Estimate interpolation model error by comparing linear vs quadratic.

    Uses three points around the crossing to fit a quadratic interpolant
    and computes the difference from the linear estimate.

    Returns NaN if fewer than 3 points are available near the crossing.
    """
    min_points_for_quadratic = 3
    n = len(b_vals)
    if n < min_points_for_quadratic or crossing_idx < 1:
        return float("nan")

    # Pick three points: prefer crossing_idx-1, crossing_idx, and one neighbor
    if crossing_idx + 1 < n:
        indices = [crossing_idx - 1, crossing_idx, crossing_idx + 1]
    elif crossing_idx >= min_points_for_quadratic - 1:
        indices = [crossing_idx - 2, crossing_idx - 1, crossing_idx]
    else:
        return float("nan")

    b3 = b_vals[indices].astype(np.float64)
    m3 = m_vals[indices].astype(np.float64)

    # Fit quadratic: m = a*b² + b*b + c
    try:
        coeffs = np.polyfit(b3, m3, 2)
    except (np.linalg.LinAlgError, ValueError):
        return float("nan")

    # Find root of quadratic - threshold = 0
    poly = np.poly1d(coeffs) - threshold
    roots = np.roots(poly)

    # Filter real roots in the bracketing interval
    b_lo, b_hi = float(b_vals[crossing_idx - 1]), float(b_vals[crossing_idx])
    real_roots = [
        float(r.real)
        for r in roots
        if abs(r.imag) < 1e-12  # noqa: PLR2004
        and b_lo <= r.real <= b_hi
    ]

    if not real_roots:
        # Expand to wider interval if no root in bracket
        real_roots = [
            float(r.real)
            for r in roots
            if abs(r.imag) < 1e-12  # noqa: PLR2004
            and b3[0] <= r.real <= b3[-1]
        ]

    if not real_roots:
        return float("nan")

    # Pick the root closest to the linear estimate
    b_min_quad = min(real_roots, key=lambda r: abs(r - b_min_linear))
    return abs(b_min_linear - b_min_quad)


def _fill_nan(
    row: dict[str, Any],
    quality: str = QUALITY_NONE,
) -> dict[str, Any]:
    """Fill B_min columns with NaN for a row with no crossing."""
    row.update(
        B_min=float("nan"),
        inv_B_min=float("nan"),
        B_min_err=float("nan"),
        B_min_err_grid=float("nan"),
        B_min_err_metric=float("nan"),
        B_min_err_interp=float("nan"),
        inv_B_min_err=float("nan"),
        crossing_quality=quality,
    )
    return row


# ── Conversion to SweepResults ──────────────────────────────────────


def critical_field_to_sweep_results(
    result: CriticalFieldResult,
    original: SweepResults,
) -> SweepResults:
    """Convert a :class:`CriticalFieldResult` to a standard :class:`SweepResults`.

    The field-strength parameter is removed from ``swept_params``, and
    ``B_min`` / ``inv_B_min`` become metric columns.  The resulting object
    can be serialized and plotted with existing sweep infrastructure.

    Parameters
    ----------
    result : CriticalFieldResult
        Output of :func:`compute_critical_field`.
    original : SweepResults
        The source sweep data (used for fixed params, sim settings, etc.).

    Returns
    -------
    SweepResults
    """
    metadata: dict[str, Any] = dict(original.metadata)
    metadata["critical_field"] = {
        "field_param": result.field_param,
        "threshold": result.threshold,
        "metric": result.metric,
        "source_spec": original.spec_path,
    }

    return SweepResults(
        swept_params=result.outer_params,
        fixed_params=dict(original.fixed_params),
        sim_settings=dict(original.sim_settings),
        rows=result.rows,
        run_dirs=[],
        spec_path=original.spec_path,
        measurements=["critical_field"],
        source_fields=original.source_fields,
        target_fields=original.target_fields,
        metadata=metadata,
    )


# ── Analytical reference formulas ───────────────────────────────────


_REFERENCE_FORMULAS = {"boccaletti", "uniform"}


def compute_reference_threshold(
    formula: str,
    reference_b: float,
    fixed_params: dict[str, Any],
    sim_settings: dict[str, Any],
) -> float:
    """Compute E-M reference conversion probability from an analytical formula.

    Parameters
    ----------
    formula : str
        ``"boccaletti"`` (localized Gaussian B-field) or ``"uniform"``
        (uniform B, periodic domain).
    reference_b : float
        Reference magnetic field strength.
    fixed_params : dict[str, Any]
        Fixed sweep parameters (must include ``"kappa"``; ``"R"`` for
        Boccaletti).
    sim_settings : dict[str, Any]
        Simulation settings (``"t_end"`` for uniform formula).

    Returns
    -------
    float
        Reference conversion probability P_EM.

    Raises
    ------
    ValueError
        If required parameters are missing or formula is unknown.
    """
    if formula not in _REFERENCE_FORMULAS:
        msg = f"Unknown reference formula '{formula}'. Valid: {_REFERENCE_FORMULAS}"
        raise ValueError(msg)

    kappa = fixed_params.get("kappa")
    if kappa is None:
        msg = f"Fixed parameter 'kappa' required for '{formula}' formula"
        raise ValueError(msg)
    kappa = float(kappa)

    if formula == "boccaletti":
        r_param = fixed_params.get("R")
        if r_param is None:
            msg = "Fixed parameter 'R' required for 'boccaletti' formula"
            raise ValueError(msg)
        r_val = float(r_param)
        # P = sin²(κ/2 · B · R · √(π/2))
        arg = kappa / 2.0 * reference_b * r_val * math.sqrt(math.pi / 2.0)
        return math.sin(arg) ** 2

    # Uniform B, periodic domain
    t_end = sim_settings.get("t_end")
    if t_end is None:
        msg = "Simulation setting 't_end' required for 'uniform' formula"
        raise ValueError(msg)
    t_end = float(t_end)
    # P = sin²(κ · B · t_end / 2)
    arg = kappa * reference_b * t_end / 2.0
    return math.sin(arg) ** 2
