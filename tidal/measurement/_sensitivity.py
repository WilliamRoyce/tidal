"""Sobol and Morris sensitivity analysis for sweep results.

Computes first-order (S1) and total-order (ST) Sobol sensitivity indices,
or Morris screening elementary effects, from existing sweep data.

Requires the optional ``SALib`` package. Import errors produce a clear
message instructing the user to install it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.measurement._sweep_results import SweepResults


def _require_salib() -> None:
    """Check that SALib is importable.

    Raises
    ------
    ImportError
        With installation instructions.
    """
    try:
        import SALib  # noqa: F401, PLC0415
    except ImportError:
        msg = (
            "SALib is required for sensitivity analysis. "
            "Install it with: pip install SALib"
        )
        raise ImportError(msg) from None


@dataclass(frozen=True)
class SensitivityResult:
    """Sensitivity analysis result (Sobol or Morris).

    Attributes
    ----------
    method : str
        ``"sobol"`` or ``"morris"``.
    metric : str
        Which metric was analyzed.
    param_names : list[str]
        Parameter names in analysis order.
    s1 : ndarray or None
        First-order Sobol indices (None for Morris).
    st : ndarray or None
        Total-order Sobol indices (None for Morris).
    s1_conf : ndarray or None
        S1 confidence intervals (None for Morris).
    st_conf : ndarray or None
        ST confidence intervals (None for Morris).
    mu_star : ndarray or None
        Morris |mean elementary effect| (None for Sobol).
    sigma : ndarray or None
        Morris std of elementary effects (None for Sobol).
    """

    method: str
    metric: str
    param_names: list[str]
    s1: NDArray[np.float64] | None
    st: NDArray[np.float64] | None
    s1_conf: NDArray[np.float64] | None
    st_conf: NDArray[np.float64] | None
    mu_star: NDArray[np.float64] | None
    sigma: NDArray[np.float64] | None


def _extract_data(
    results: SweepResults,
    metric: str,
) -> tuple[list[str], NDArray[np.float64], NDArray[np.float64]]:
    """Extract parameter matrix X and metric vector Y from results.

    Returns (param_names, X, Y).

    Raises
    ------
    ValueError
        If no swept parameters in results.
    """
    param_names = list(results.swept_params.keys())
    if not param_names:
        msg = "No swept parameters in results"
        raise ValueError(msg)

    n = len(results.rows)
    x = np.zeros((n, len(param_names)), dtype=np.float64)
    y = np.zeros(n, dtype=np.float64)

    for i, row in enumerate(results.rows):
        for j, name in enumerate(param_names):
            x[i, j] = float(row.get(name, 0.0))
        val = row.get(metric)
        y[i] = float(val) if val is not None else np.nan

    # Remove NaN rows
    valid = ~np.isnan(y)
    return param_names, x[valid], y[valid]


def compute_sobol_indices(
    results: SweepResults,
    metric: str,
    *,
    n_bootstrap: int = 100,
) -> SensitivityResult:
    """Compute Sobol sensitivity indices from sweep results.

    Uses SALib's analyze module with Saltelli-style computation.
    The sweep data should ideally come from a Sobol or Latin Hypercube
    sampling strategy (F2b) for reliable estimates.

    Parameters
    ----------
    results : SweepResults
        Sweep data with 2+ swept parameters.
    metric : str
        Metric column name.
    n_bootstrap : int
        Number of bootstrap resamples for confidence intervals.

    Returns
    -------
    SensitivityResult

    Raises
    ------
    ValueError
        If insufficient data for analysis.
    """
    _require_salib()
    from SALib.analyze import (  # noqa: PLC0415
        sobol as sobol_analyze,  # type: ignore[import-untyped]
    )

    param_names, x, y = _extract_data(results, metric)

    if len(y) < 2 * (len(param_names) + 1):
        msg = (
            f"Sobol analysis needs at least {2 * (len(param_names) + 1)} samples "
            f"for {len(param_names)} parameters, got {len(y)}. "
            f"Use --sweep-strategy sobol --n-samples N for proper sampling."
        )
        raise ValueError(msg)

    # Build SALib problem definition
    bounds = [
        [float(x[:, j].min()), float(x[:, j].max())] for j in range(len(param_names))
    ]
    problem = {
        "num_vars": len(param_names),
        "names": param_names,
        "bounds": bounds,
    }

    si = sobol_analyze.analyze(
        problem, y, calc_second_order=False, num_resamples=n_bootstrap
    )

    return SensitivityResult(
        method="sobol",
        metric=metric,
        param_names=param_names,
        s1=np.asarray(si["S1"], dtype=np.float64),
        st=np.asarray(si["ST"], dtype=np.float64),
        s1_conf=np.asarray(si["S1_conf"], dtype=np.float64),
        st_conf=np.asarray(si["ST_conf"], dtype=np.float64),
        mu_star=None,
        sigma=None,
    )


def compute_morris_screening(
    results: SweepResults,
    metric: str,
) -> SensitivityResult:
    """Compute Morris elementary effects from sweep results.

    Morris screening is faster than Sobol and works well for identifying
    which parameters are important vs unimportant. Best with grid or
    Latin Hypercube sampling.

    Parameters
    ----------
    results : SweepResults
        Sweep data.
    metric : str
        Metric column name.

    Returns
    -------
    SensitivityResult

    Raises
    ------
    ValueError
        If insufficient data.
    """
    _require_salib()
    from SALib.analyze import (  # noqa: PLC0415
        morris as morris_analyze,  # type: ignore[import-untyped]
    )

    param_names, x, y = _extract_data(results, metric)

    if len(y) < len(param_names) + 1:
        msg = (
            f"Morris screening needs at least {len(param_names) + 1} samples "
            f"for {len(param_names)} parameters, got {len(y)}"
        )
        raise ValueError(msg)

    bounds = [
        [float(x[:, j].min()), float(x[:, j].max())] for j in range(len(param_names))
    ]
    problem = {
        "num_vars": len(param_names),
        "names": param_names,
        "bounds": bounds,
    }

    si = morris_analyze.analyze(problem, x, y)

    return SensitivityResult(
        method="morris",
        metric=metric,
        param_names=param_names,
        s1=None,
        st=None,
        s1_conf=None,
        st_conf=None,
        mu_star=np.asarray(si["mu_star"], dtype=np.float64),
        sigma=np.asarray(si["sigma"], dtype=np.float64),
    )


def format_sensitivity_table(result: SensitivityResult) -> str:
    """Format a sensitivity result as a human-readable table."""
    lines = [
        f"{result.method.title()} Sensitivity Analysis: {result.metric}",
        "=" * 60,
    ]

    if result.method == "sobol" and result.s1 is not None and result.st is not None:
        lines.extend(
            (
                f"{'Parameter':<15} {'S1 (main)':>12} {'ST (total)':>12} {'Interaction':>12}",
                "-" * 60,
            )
        )
        for i, name in enumerate(result.param_names):
            s1 = result.s1[i]
            st = result.st[i]
            interaction = (
                "strong"
                if (st - s1) > 0.1  # noqa: PLR2004
                else ("moderate" if (st - s1) > 0.03 else "weak")  # noqa: PLR2004
            )
            s1_str = f"{s1:.3f}"
            st_str = f"{st:.3f}"
            if result.s1_conf is not None:
                s1_str += f" +/- {result.s1_conf[i]:.3f}"
            if result.st_conf is not None:
                st_str += f" +/- {result.st_conf[i]:.3f}"
            lines.append(f"{name:<15} {s1_str:>12} {st_str:>12} {interaction:>12}")

    elif (
        result.method == "morris"
        and result.mu_star is not None
        and result.sigma is not None
    ):
        lines.extend(
            (f"{'Parameter':<15} {'mu* (mean)':>12} {'sigma (std)':>12}", "-" * 60)
        )
        # Sort by mu_star descending
        order = np.argsort(-result.mu_star)
        lines.extend(
            f"{result.param_names[i]:<15} {result.mu_star[i]:>12.4f} {result.sigma[i]:>12.4f}"
            for i in order
        )

    return "\n".join(lines)
