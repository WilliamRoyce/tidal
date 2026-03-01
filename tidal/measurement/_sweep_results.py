"""Data model for parameter sweep results.

Provides ``SweepResults``, a frozen dataclass that stores aggregated
metrics from a parameter sweep.  Can be loaded from a sweep output
directory (``sweep.json`` + ``results.json``) or constructed in-memory
from sweep execution.

The CSV/JSON output is designed to be **self-contained and portable** —
every row includes all parameters (swept + fixed), simulation settings,
and measured metrics so the data can be analyzed outside TIDAL.
"""

from __future__ import annotations

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


@dataclass
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
    metadata: dict[str, Any] = field(default_factory=dict)
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
            if all(np.isclose(row.get(k, float("nan")), v) for k, v in kwargs.items()):
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
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
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
