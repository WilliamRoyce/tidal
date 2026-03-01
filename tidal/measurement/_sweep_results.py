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
    # Serialization
    # ------------------------------------------------------------------

    def _csv_columns(self) -> list[str]:
        """Ordered list of CSV column names."""
        cols: list[str] = []
        # Swept params first
        cols.extend(self.swept_params.keys())
        # Fixed params
        cols.extend(self.fixed_params.keys())
        # Simulation settings
        cols.extend(self.sim_settings.keys())
        # Metrics (union of all rows, preserving insertion order)
        if self.rows:
            seen = set(cols)
            for row in self.rows:
                for k in row:
                    if k not in seen:
                        cols.append(k)
                        seen.add(k)
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

        path.write_text(json.dumps(data, indent=2, default=_json_default), encoding="utf-8")

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
