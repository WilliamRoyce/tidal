"""Data model for parameter sweep results.

Provides ``SweepResults``, a frozen dataclass that stores aggregated
metrics from a parameter sweep.  Can be loaded from a sweep output
directory (``sweep.json`` + ``results.json``) or constructed in-memory
from sweep execution.

The CSV/JSON output is designed to be **self-contained and portable** —
every row includes all parameters (swept + fixed), simulation settings,
and measured metrics so the data can be analyzed outside TIDAL.

References
----------
Smith, R.C. (2013) *Uncertainty Quantification: Theory, Implementation,
and Applications*, SIAM. Ch. 3 (sample statistics, SEM, CI).

Efron, B. & Tibshirani, R. (1993) *An Introduction to the Bootstrap*,
Chapman & Hall. Ch. 12-14 (non-parametric bootstrap confidence intervals).
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = ["SweepResults"]

# Default metric candidates for auto-detection (shared across CLI modules)
DEFAULT_METRIC_CANDIDATES: tuple[str, ...] = (
    "P_max",
    "max_energy_error",
    "L_mix",
    "E_total_final",
)


@dataclass  # noqa: PLR0904
class SweepResults:
    """Aggregated results from a parameter sweep.

    Attributes
    ----------
    swept_params : dict[str, list[float]]
        Parameter name -> list of swept values.
    fixed_params : dict[str, float]
        Non-swept parameters (constant across all runs).
    sim_settings : dict[str, Any]
        Simulation settings (grid_shape, t_end, dt, scheme, bc, etc.).
    rows : list[dict[str, Any]]
        One dict per run.  Keys = all param names + metric names +
        ``"wall_time_s"``.  Values = floats or None.
    run_dirs : list[Path]
        Paths to individual run output directories.
    spec_path : str
        Path to the JSON equation specification used.
    measurements : list[str]
        Measurement types requested (e.g. ``["conversion", "mixing"]``).
    source_fields : list[str] | None
        Source field names for conversion measurements.
    target_fields : list[str] | None
        Target field names for conversion measurements.
    metadata : dict[str, Any]
        Provenance metadata (timestamp, version, etc.).
    converge_sizes : list[int] | None
        Grid sizes for convergence mode (None for parameter sweeps).
    """

    swept_params: dict[str, list[float]]
    fixed_params: dict[str, float]
    sim_settings: dict[str, Any]
    rows: list[dict[str, Any]]
    run_dirs: list[Path]
    spec_path: str
    measurements: list[str]
    source_fields: list[str] | None = None
    target_fields: list[str] | None = None
    metadata: dict[str, Any] = field(default_factory=dict[str, Any])
    converge_sizes: list[int] | None = None

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_runs(self) -> int:
        """Number of completed runs."""
        return len(self.rows)

    @property
    def param_names(self) -> list[str]:
        """All parameter names (swept first, then fixed)."""
        return list(self.swept_params.keys()) + list(self.fixed_params.keys())

    @property
    def metric_names(self) -> list[str]:
        """Metric column names (union of all rows, excluding params and settings)."""
        if not self.rows:
            return []
        param_keys = set(self.swept_params.keys()) | set(self.fixed_params.keys())
        setting_keys = set(self.sim_settings.keys())
        seen: set[str] = set()
        result: list[str] = []
        for row in self.rows:
            for k in row:
                if k not in param_keys and k not in setting_keys and k not in seen:
                    result.append(k)
                    seen.add(k)
        return result

    @property
    def is_convergence(self) -> bool:
        """Whether this is a convergence study (vs parameter sweep)."""
        return self.converge_sizes is not None

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    def _with_rows(
        self,
        rows: list[dict[str, Any]],
        run_dirs: list[Path] | None = None,
    ) -> SweepResults:
        """Create a new SweepResults sharing metadata but with different rows."""
        return SweepResults(
            swept_params=self.swept_params,
            fixed_params=self.fixed_params,
            sim_settings=self.sim_settings,
            rows=rows,
            run_dirs=run_dirs if run_dirs is not None else [],
            spec_path=self.spec_path,
            measurements=self.measurements,
            source_fields=self.source_fields,
            target_fields=self.target_fields,
            metadata=self.metadata,
            converge_sizes=self.converge_sizes,
        )

    # ------------------------------------------------------------------
    # Column access
    # ------------------------------------------------------------------

    def column(self, name: str) -> np.ndarray:
        """Extract a column as a numpy array.

        Parameters
        ----------
        name : str
            Column name (parameter name, setting name, or metric name).

        Returns
        -------
        ndarray
            Values for that column across all runs.
        """
        values = [row.get(name) for row in self.rows]
        return np.array(values, dtype=np.float64)

    # ------------------------------------------------------------------
    # Status tracking
    # ------------------------------------------------------------------

    @property
    def failure_rate(self) -> float:
        """Fraction of runs that failed (``run_status != 'success'``).

        Returns 0.0 if there are no rows or no ``run_status`` column.
        """
        if not self.rows:
            return 0.0
        n_failed = sum(
            1 for r in self.rows if r.get("run_status", "success") != "success"
        )
        return n_failed / len(self.rows)

    def successful_rows(self) -> SweepResults:
        """Return only rows with ``run_status == 'success'``."""
        return self._filter_by_status("success", invert=False)

    def failed_rows(self) -> SweepResults:
        """Return only rows with ``run_status != 'success'``."""
        return self._filter_by_status("success", invert=True)

    def _filter_by_status(self, status: str, *, invert: bool) -> SweepResults:
        """Filter rows by run_status value."""
        filtered_rows: list[dict[str, Any]] = []
        filtered_dirs: list[Path] = []
        for i, row in enumerate(self.rows):
            matches = row.get("run_status", "success") == status
            if matches != invert:
                filtered_rows.append(row)
                if i < len(self.run_dirs):
                    filtered_dirs.append(self.run_dirs[i])
        return self._with_rows(filtered_rows, filtered_dirs)

    # ------------------------------------------------------------------
    # Query methods
    # ------------------------------------------------------------------

    def filter(self, **kwargs: float) -> SweepResults:
        """Return a new SweepResults with only rows matching parameter values.

        Uses ``numpy.isclose`` for float comparison.

        Parameters
        ----------
        **kwargs : float
            Parameter name/value pairs to filter on.

        Returns
        -------
        SweepResults
            Filtered results (preserves all metadata).
        """
        filtered_rows: list[dict[str, Any]] = []
        filtered_dirs: list[Path] = []
        for i, row in enumerate(self.rows):
            if all(
                bool(np.isclose(row.get(k, float("nan")), v))  # pyright: ignore[reportUnknownArgumentType]
                for k, v in kwargs.items()
            ):
                filtered_rows.append(row)
                if i < len(self.run_dirs):
                    filtered_dirs.append(self.run_dirs[i])
        return self._with_rows(filtered_rows, filtered_dirs)

    def group_by(self, param: str) -> dict[float, SweepResults]:
        """Group results by distinct values of a parameter.

        Parameters
        ----------
        param : str
            Parameter name to group on.

        Returns
        -------
        dict[float, SweepResults]
            Mapping from parameter value to filtered results.
        """
        groups: dict[float, list[int]] = {}
        for i, row in enumerate(self.rows):
            val = row.get(param)
            if val is None:
                continue
            key = float(val)
            groups.setdefault(key, []).append(i)

        result: dict[float, SweepResults] = {}
        for val, indices in groups.items():
            result[val] = self._with_rows(
                [self.rows[i] for i in indices],
                [self.run_dirs[i] for i in indices if i < len(self.run_dirs)],
            )
        return result

    def best(self, metric: str, *, maximize: bool = True) -> dict[str, Any]:
        """Return the row with the best (max or min) metric value.

        Parameters
        ----------
        metric : str
            Metric column name.
        maximize : bool
            If True, return the row with the highest value; otherwise lowest.

        Returns
        -------
        dict[str, Any]
            The best row.

        Raises
        ------
        KeyError
            If *metric* is not found in any row.
        ValueError
            If all metric values are None.
        """
        best_val: float | None = None
        best_row: dict[str, Any] | None = None
        found_key = False
        for row in self.rows:
            if metric not in row:
                continue
            found_key = True
            val = row[metric]
            if val is None:
                continue
            fval = float(val)
            if (
                best_val is None
                or (maximize and fval > best_val)
                or (not maximize and fval < best_val)
            ):
                best_val = fval
                best_row = row
        if not found_key:
            msg = f"Metric '{metric}' not found in any row"
            raise KeyError(msg)
        if best_row is None:
            msg = f"All values for metric '{metric}' are None"
            raise ValueError(msg)
        return best_row

    def summary(self) -> str:
        """Return a human-readable summary of the sweep results.

        Returns
        -------
        str
            Formatted multi-line summary.
        """
        header = f"SweepResults summary ({self.n_runs} runs)"
        lines: list[str] = [header, "=" * len(header)]

        if self.swept_params:
            lines.append("\nSwept parameters:")
            for name, vals in self.swept_params.items():
                lines.append(
                    f"  {name}: {min(vals):.6g} to {max(vals):.6g} ({len(vals)} values)"
                )

        if self.fixed_params:
            lines.append("\nFixed parameters:")
            for name, val in self.fixed_params.items():
                lines.append(f"  {name}: {val:.6g}")

        metrics = self.metric_names
        if metrics and self.rows:
            lines.append("\nMetrics:")
            for m in metrics:
                vals = [row.get(m) for row in self.rows if row.get(m) is not None]
                if not vals:
                    lines.append(f"  {m}: no data")
                    continue
                arr = np.array(vals, dtype=np.float64)
                if len(arr) == 1:
                    lines.append(f"  {m}: {arr[0]:.6g} (1 value)")
                else:
                    lines.append(
                        f"  {m}: mean={np.mean(arr):.6g} \u00b1 {np.std(arr):.6g}, "
                        f"min={np.min(arr):.6g}, max={np.max(arr):.6g}"
                    )

        return "\n".join(lines)

    # ------------------------------------------------------------------
    # Ensemble / replicate support
    # ------------------------------------------------------------------

    @property
    def has_replicates(self) -> bool:
        """Whether results contain replicate runs."""
        return any("replicate" in row for row in self.rows)

    @property
    def n_replicates(self) -> int:
        """Number of replicates (max replicate index + 1, or 1 if none)."""
        if not self.has_replicates:
            return 1
        max_rep = max(row.get("replicate", 0) for row in self.rows)
        return int(max_rep) + 1

    def group_by_point(self) -> dict[tuple[float, ...], list[dict[str, Any]]]:
        """Group rows by parameter-point identity, ignoring replicate/seed.

        Uses ``{param}_nominal`` columns if present (param noise),
        otherwise the swept parameter values themselves.

        Returns
        -------
        dict[tuple[float, ...], list[dict[str, Any]]]
            Mapping from parameter-point tuple to list of replicate rows.
        """
        param_names = list(self.swept_params.keys())
        groups: dict[tuple[float, ...], list[dict[str, Any]]] = {}
        for row in self.rows:
            key_vals: list[float] = []
            for pname in param_names:
                # Use nominal value if available (param noise), else realized
                nominal_key = f"{pname}_nominal"
                val = row.get(nominal_key, row.get(pname))
                key_vals.append(float(val) if val is not None else float("nan"))
            key = tuple(key_vals)
            groups.setdefault(key, []).append(row)
        return groups

    def aggregate(
        self,
        metrics: list[str] | None = None,
    ) -> SweepResults:
        """Aggregate replicate rows into one row per parameter point.

        Computes per-metric statistics: mean, std, sem (std/sqrt(N)),
        95% CI bounds, median, 25th/75th percentiles, and count.

        The plain ``{metric}`` column holds the mean for backward
        compatibility with existing plot code.

        Parameters
        ----------
        metrics : list[str] or None
            Metrics to aggregate. If None, uses all metric_names.

        Returns
        -------
        SweepResults
            New instance with one row per parameter point and
            ``{metric}_mean``, ``{metric}_std``, etc. columns.
        """
        if not self.has_replicates:
            return self

        # Columns that are never aggregated (non-numeric or metadata)
        skip = {
            "replicate",
            "seed",
            "run_status",
            "error_message",
            "solver_exit_code",
            "error",
            "energy_conserved",
            "conservation_error",
        }
        if metrics is None:
            metrics = [
                m
                for m in self.metric_names
                if m not in skip and not m.endswith("_nominal")
            ]

        from scipy.stats import (  # noqa: PLC0415  # type: ignore[reportMissingTypeStubs]
            t as t_dist,
        )

        param_names = list(self.swept_params.keys())
        groups = self.group_by_point()
        agg_rows: list[dict[str, Any]] = []

        for point_key, rep_rows in groups.items():
            row: dict[str, Any] = dict(zip(param_names, point_key, strict=True))

            # Fixed params and sim settings from first row
            row.update(self.fixed_params)
            row.update(
                {k: rep_rows[0].get(k) for k in self.sim_settings if k in rep_rows[0]}
            )

            # Aggregate each metric
            for m in metrics:
                vals = [
                    r[m]
                    for r in rep_rows
                    if r.get(m) is not None and r.get("run_status") == "success"
                ]
                n = len(vals)
                row[f"{m}_n"] = n
                if n == 0:
                    row[m] = None
                    for suffix in (
                        "_mean",
                        "_std",
                        "_sem",
                        "_ci95_lo",
                        "_ci95_hi",
                        "_median",
                        "_q25",
                        "_q75",
                    ):
                        row[f"{m}{suffix}"] = None
                    continue

                arr = np.array(vals, dtype=np.float64)
                mean = float(np.mean(arr))
                row[m] = mean  # backward compat
                row[f"{m}_mean"] = mean

                if n == 1:
                    row[f"{m}_std"] = 0.0
                    row[f"{m}_sem"] = 0.0
                    row[f"{m}_ci95_lo"] = mean
                    row[f"{m}_ci95_hi"] = mean
                    row[f"{m}_median"] = mean
                    row[f"{m}_q25"] = mean
                    row[f"{m}_q75"] = mean
                else:
                    std = float(np.std(arr, ddof=1))
                    sem = std / np.sqrt(n)
                    row[f"{m}_std"] = std
                    row[f"{m}_sem"] = sem

                    # 95% CI via t-distribution
                    t_crit = float(t_dist.ppf(0.975, df=n - 1))
                    row[f"{m}_ci95_lo"] = mean - t_crit * sem
                    row[f"{m}_ci95_hi"] = mean + t_crit * sem

                    row[f"{m}_median"] = float(np.median(arr))
                    row[f"{m}_q25"] = float(np.percentile(arr, 25))
                    row[f"{m}_q75"] = float(np.percentile(arr, 75))

            agg_rows.append(row)

        return self._with_rows(agg_rows)

    def bootstrap_ci(
        self,
        metric: str,
        *,
        n_bootstrap: int = 10_000,
        alpha: float = 0.05,
        seed: int = 0,
    ) -> SweepResults:
        """Compute non-parametric bootstrap confidence intervals.

        For each parameter point, resamples the replicate metric values
        with replacement *n_bootstrap* times, computing the mean of each
        resample. The CI bounds are the alpha/2 and 1-alpha/2 percentiles
        of the bootstrap distribution.

        More robust than t-distribution CI for non-Gaussian distributions
        or very small sample sizes.

        Parameters
        ----------
        metric : str
            Metric column to bootstrap.
        n_bootstrap : int
            Number of bootstrap resamples (default 10,000).
        alpha : float
            Significance level (default 0.05 for 95% CI).
        seed : int
            RNG seed for reproducibility.

        Returns
        -------
        SweepResults
            New instance with one row per parameter point, containing
            ``{metric}_boot_lo`` and ``{metric}_boot_hi`` columns.

        References
        ----------
        Efron, B. & Tibshirani, R. (1993) *An Introduction to the
        Bootstrap*, Chapman & Hall. Ch. 12-14.
        """
        if not self.has_replicates:
            return self

        rng = np.random.default_rng(seed)
        param_names = list(self.swept_params.keys())
        groups = self.group_by_point()
        boot_rows: list[dict[str, Any]] = []

        lo_pct = 100 * alpha / 2
        hi_pct = 100 * (1 - alpha / 2)

        for point_key, rep_rows in groups.items():
            row: dict[str, Any] = dict(zip(param_names, point_key, strict=True))
            row.update(self.fixed_params)
            row.update(
                {k: rep_rows[0].get(k) for k in self.sim_settings if k in rep_rows[0]}
            )

            vals = [
                r[metric]
                for r in rep_rows
                if r.get(metric) is not None and r.get("run_status") == "success"
            ]
            n = len(vals)
            if n == 0:
                row[f"{metric}_boot_lo"] = None
                row[f"{metric}_boot_hi"] = None
            elif n == 1:
                row[f"{metric}_boot_lo"] = vals[0]
                row[f"{metric}_boot_hi"] = vals[0]
            else:
                arr = np.array(vals, dtype=np.float64)
                # Bootstrap: resample with replacement, compute mean
                boot_means = np.array(
                    [
                        float(np.mean(rng.choice(arr, size=n, replace=True)))
                        for _ in range(n_bootstrap)
                    ],
                    dtype=np.float64,
                )
                row[f"{metric}_boot_lo"] = float(np.percentile(boot_means, lo_pct))
                row[f"{metric}_boot_hi"] = float(np.percentile(boot_means, hi_pct))

            # Copy the plain metric as mean for convenience
            if vals:
                row[metric] = float(np.mean(vals))
            else:
                row[metric] = None

            boot_rows.append(row)

        return self._with_rows(boot_rows)

    # ------------------------------------------------------------------
    # Derived metrics
    # ------------------------------------------------------------------

    def add_derived_columns(
        self,
        baseline_formula: str = "sin(kappa * B0 * t_end / 2)**2",
        *,
        metric: str = "P_max",
        validity_threshold: float = 0.1,
    ) -> None:
        """Add amplification-related derived columns in place.

        Computes for each row:
        - ``P_EM``: baseline conversion probability from *baseline_formula*
        - ``C0``: conversion coefficient ``P / B0**2``
        - ``A``: amplification factor ``P / P_EM``
        - ``log10_A``: ``log10(A)`` (useful for divergent colormaps)
        - ``P_valid``: ``True`` if ``0 < P < validity_threshold``

        Parameters
        ----------
        baseline_formula : str
            Python expression for the GR baseline probability.  May
            reference any fixed or swept parameter name plus standard
            math functions (``sin``, ``cos``, ``sqrt``, ``pi``, etc.)
            and simulation settings (``t_end``, ``grid_shape``).
        metric : str
            Column name for the measured conversion probability.
        validity_threshold : float
            Maximum *P* for the linear regime (default 0.1).
        """
        import math as _math  # noqa: PLC0415

        from tidal.cli._simulate import (  # noqa: PLC0415
            FORMULA_NAMESPACE,  # pyright: ignore[reportPrivateUsage]
        )

        for row in self.rows:
            p_val = row.get(metric)
            status = row.get("run_status", "success")

            # Build namespace from fixed params, sim settings, and this row
            ns: dict[str, object] = {**FORMULA_NAMESPACE}
            ns.update(self.fixed_params)
            for k, v in self.sim_settings.items():
                with contextlib.suppress(TypeError, ValueError):
                    ns[k] = float(v)
            # Swept params from row override fixed
            for pname in self.swept_params:
                val = row.get(pname)
                if val is not None:
                    with contextlib.suppress(TypeError, ValueError):
                        ns[pname] = float(val)

            # Evaluate baseline
            try:
                p_em = float(
                    eval(baseline_formula, {"__builtins__": {}}, ns)  # noqa: S307
                )
            except Exception:  # noqa: BLE001
                p_em = _math.nan

            # Compute derived values
            try:
                p = float(p_val)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                p = _math.nan

            valid = status == "success" and 0 < p < validity_threshold
            b0 = float(ns.get("B0", 0) or 0)  # type: ignore[arg-type]

            row["P_EM"] = p_em if _math.isfinite(p_em) else None
            row["C0"] = p / b0**2 if b0 > 0 and _math.isfinite(p) else None
            row["A"] = p / p_em if p_em > 0 and valid and _math.isfinite(p) else None
            row["log10_A"] = (
                _math.log10(p / p_em)
                if p_em > 0 and valid and p > 0 and _math.isfinite(p)
                else None
            )
            row["P_valid"] = valid

    def add_paired_baseline(
        self,
        baseline: SweepResults,
        *,
        metric: str = "P_max",
        match_tolerance: float = 1e-6,
    ) -> None:
        """Add amplification relative to a paired baseline sweep.

        Matches each signal row with the baseline row whose parameter
        values are closest (within *match_tolerance*).  Adds columns:

        - ``P_baseline``: matched baseline metric value
        - ``A_paired``: ``signal_metric / P_baseline``
        - ``log10_A_paired``: ``log10(A_paired)``

        The baseline sweep must cover the same parameter grid for all
        shared swept parameters.  Unmatched rows get ``None``.

        Parameters
        ----------
        baseline : SweepResults
            Baseline sweep (typically at deltam=0).
        metric : str
            Column name for the conversion probability.
        match_tolerance : float
            Maximum absolute difference for parameter matching.
        """
        import math as _math  # noqa: PLC0415

        # Identify shared swept parameters between signal and baseline
        base_params = set(baseline.swept_params.keys())
        sig_params = set(self.swept_params.keys())
        shared = base_params & sig_params

        # Also include baseline fixed params that are signal swept params
        # (e.g., baseline has deltam=0 fixed, signal sweeps deltam)
        # → don't match on those
        match_params = sorted(shared) if shared else sorted(base_params)

        # Build lookup: tuple of rounded param values → baseline P_max
        base_lookup: dict[tuple[float, ...], float] = {}
        for row in baseline.rows:
            if row.get("run_status") != "success":
                continue
            key = tuple(
                round(float(row.get(p, 0)), 6)
                for p in match_params  # type: ignore[arg-type]
            )
            val = row.get(metric)
            if val is not None:
                base_lookup[key] = float(val)

        n_matched = 0
        for row in self.rows:
            key = tuple(
                round(float(row.get(p, 0)), 6)
                for p in match_params  # type: ignore[arg-type]
            )
            p_base = base_lookup.get(key)

            if p_base is None:
                # Try fuzzy match within tolerance
                for bkey, bval in base_lookup.items():
                    if all(abs(a - b) <= match_tolerance for a, b in zip(key, bkey, strict=False)):
                        p_base = bval
                        break

            p_val = row.get(metric)
            try:
                p = float(p_val)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                p = _math.nan

            if p_base is not None and p_base > 0 and _math.isfinite(p):
                a = p / p_base
                row["P_baseline"] = p_base
                row["A_paired"] = a
                row["log10_A_paired"] = _math.log10(a) if a > 0 else None
                n_matched += 1
            else:
                row["P_baseline"] = p_base
                row["A_paired"] = None
                row["log10_A_paired"] = None

        import logging  # noqa: PLC0415

        logging.getLogger(__name__).info(
            "Paired baseline: %d/%d rows matched (%d baseline points, matching on %s)",
            n_matched,
            len(self.rows),
            len(base_lookup),
            match_params,
        )

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def _csv_columns(self) -> list[str]:
        """Ordered list of CSV column names."""
        cols: list[str] = []
        cols.extend(self.swept_params.keys())
        cols.extend(self.fixed_params.keys())
        cols.extend(self.sim_settings.keys())
        cols.extend(self.metric_names)
        return cols

    def to_csv(self, path: Path | None = None) -> str:
        """Write results to CSV.

        Parameters
        ----------
        path : Path or None
            If given, write to this file.  Always returns the CSV string.

        Returns
        -------
        str
            CSV content.
        """
        columns = self._csv_columns()
        buf = io.StringIO()
        writer = csv.DictWriter(buf, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for row in self.rows:
            writer.writerow(row)
        content = buf.getvalue()
        if path is not None:
            path.write_text(content)
        return content

    def to_json(self, path: Path | None = None) -> str:
        """Write results to JSON.

        Parameters
        ----------
        path : Path or None
            If given, write to this file.  Always returns the JSON string.

        Returns
        -------
        str
            JSON content.
        """
        data = {
            "columns": self._csv_columns(),
            "rows": self.rows,
        }
        content = json.dumps(data, indent=2, default=_json_default)
        if path is not None:
            path.write_text(content)
        return content

    def save_sweep_json(self, path: Path) -> None:
        """Write ``sweep.json`` provenance file.

        Parameters
        ----------
        path : Path
            Output file path.
        """
        data: dict[str, Any] = {
            "version": 1,
            "spec_path": self.spec_path,
            "swept_parameters": self.swept_params,
            "fixed_parameters": self.fixed_params,
            "sim_settings": self.sim_settings,
            "measurements": self.measurements,
        }
        if self.source_fields is not None:
            data["source"] = self.source_fields
        if self.target_fields is not None:
            data["target"] = self.target_fields
        if self.converge_sizes is not None:
            data["converge_sizes"] = self.converge_sizes
        data["run_dirs"] = [str(d) for d in self.run_dirs]
        data.update(self.metadata)
        data["completed_runs"] = self.n_runs
        data["total_runs"] = self.metadata.get("total_runs", self.n_runs)

        path.write_text(
            json.dumps(data, indent=2, default=_json_default), encoding="utf-8"
        )

    # ------------------------------------------------------------------
    # Loading
    # ------------------------------------------------------------------

    @classmethod
    def from_directory(cls, path: Path | str) -> SweepResults:
        """Load sweep results from a directory.

        Expects ``sweep.json`` and ``results.json`` in *path*.

        Parameters
        ----------
        path : Path or str
            Sweep output directory.

        Returns
        -------
        SweepResults

        Raises
        ------
        FileNotFoundError
            If ``sweep.json`` or ``results.json`` is missing.
        """
        dirp = Path(path)
        sweep_file = dirp / "sweep.json"
        results_file = dirp / "results.json"

        if not sweep_file.exists():
            msg = f"No sweep.json in {dirp}"
            raise FileNotFoundError(msg)
        if not results_file.exists():
            msg = f"No results.json in {dirp}"
            raise FileNotFoundError(msg)

        sweep_data = json.loads(sweep_file.read_text())
        results_data = json.loads(results_file.read_text())

        # Use stored run_dirs if available, else infer (backward compat)
        stored_dirs = sweep_data.get("run_dirs")
        run_dirs: list[Path] = []
        if stored_dirs:
            run_dirs = [Path(d) for d in stored_dirs]
        else:
            for row in results_data.get("rows", []):
                run_dir = _infer_run_dir(dirp, sweep_data, row)
                run_dirs.append(run_dir)

        return cls(
            swept_params=sweep_data.get("swept_parameters", {}),
            fixed_params=sweep_data.get("fixed_parameters", {}),
            sim_settings=sweep_data.get("sim_settings", {}),
            rows=results_data.get("rows", []),
            run_dirs=run_dirs,
            spec_path=sweep_data.get("spec_path", ""),
            measurements=sweep_data.get("measurements", []),
            source_fields=sweep_data.get("source"),
            target_fields=sweep_data.get("target"),
            metadata={
                k: v
                for k, v in sweep_data.items()
                if k
                not in {
                    "version",
                    "spec_path",
                    "swept_parameters",
                    "fixed_parameters",
                    "sim_settings",
                    "measurements",
                    "source",
                    "target",
                    "converge_sizes",
                    "completed_runs",
                    "total_runs",
                    "run_dirs",
                }
            },
            converge_sizes=sweep_data.get("converge_sizes"),
        )


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------


def _json_default(obj: object) -> object:
    """Serialize numpy types and Paths for JSON output.

    Raises
    ------
    TypeError
        If *obj* is not a recognized serializable type.
    """
    if isinstance(obj, np.integer):
        return int(obj)  # type: ignore[call-overload]  # numpy integer subtype
    if isinstance(obj, np.floating):
        return float(obj)  # type: ignore[call-overload]  # numpy floating subtype
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    msg = f"Object of type {type(obj).__name__} is not JSON serializable"
    raise TypeError(msg)


def _infer_run_dir(
    sweep_dir: Path,
    sweep_data: dict[str, Any],
    row: dict[str, Any],
) -> Path:
    """Infer the run subdirectory from a results row."""
    swept = sweep_data.get("swept_parameters", {})
    converge = sweep_data.get("converge_sizes")

    if converge is not None:
        # Convergence mode: subdirs named by grid size
        gs = row.get("grid_shape", "")
        return sweep_dir / f"N_{gs}"

    # Parameter sweep: subdir named from swept param values
    parts: list[str] = []
    for param_name in swept:
        val = row.get(param_name, "")
        parts.append(f"{param_name}_{val}")
    subdir_name = "_".join(parts) if parts else "run_0"
    return sweep_dir / subdir_name
