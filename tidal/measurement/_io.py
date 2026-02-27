"""Uniform data abstraction for simulation output.

Provides ``SimulationData``, a frozen dataclass that stores the full time
history of field and momentum arrays.  Can be constructed from:

- A solver result dict (IDA/leapfrog output)
- A snapshot directory written by :class:`~tidal.measurement.SnapshotWriter`
  (memory-mapped, O(1) RAM)
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from math import prod
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

from tidal.solver.state import StateLayout

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from tidal.solver._types import SolverResult
    from tidal.solver.grid import GridInfo
    from tidal.symbolic.json_loader import EquationSystem


@dataclass(frozen=True)
class SimulationData:
    """Full time-history of a simulation, ready for measurement.

    Both *fields* and *momenta* store arrays of shape
    ``(n_snapshots, *grid_shape)`` — one spatial snapshot per recorded time.
    Constraint fields (``time_derivative_order == 0``) have no momentum entry.

    Attributes
    ----------
    times : ndarray, shape (n_snapshots,)
        Snapshot times.
    fields : dict[str, ndarray]
        Mapping ``field_name → (n_snapshots, *grid_shape)`` arrays.
    momenta : dict[str, ndarray]
        Mapping ``field_name → (n_snapshots, *grid_shape)`` arrays.
        Only present for 2nd-order (wave) fields.
    grid_spacing : tuple[float, ...]
        Cell size per spatial axis, e.g. ``(dx, dy)``.
    grid_bounds : tuple[tuple[float, float], ...]
        Domain bounds per spatial axis.
    periodic : tuple[bool, ...]
        Whether each spatial axis is periodic.
    spec : EquationSystem
        The equation specification (fields, equations, matrices).
    parameters : dict[str, float]
        Resolved parameter values used in the simulation.
    bc_types : tuple[str, ...] or None
        Per-axis boundary condition type (e.g. ``("periodic", "neumann")``).
        ``None`` for legacy data where BC info was not recorded.
        Used by the energy module for BC-aware gradient computation.
    """

    times: NDArray[np.float64]
    fields: dict[str, NDArray[np.float64]]
    momenta: dict[str, NDArray[np.float64]]
    grid_spacing: tuple[float, ...]
    grid_bounds: tuple[tuple[float, float], ...]
    periodic: tuple[bool, ...]
    spec: EquationSystem
    parameters: dict[str, float]
    bc_types: tuple[str, ...] | None = None
    dt: float | None = None

    # ------------------------------------------------------------------
    # Derived helpers
    # ------------------------------------------------------------------

    @property
    def n_snapshots(self) -> int:
        """Number of time snapshots."""
        return len(self.times)

    @property
    def volume_element(self) -> float:
        """Uniform cell volume ``dx * dy * ...``."""
        return float(prod(self.grid_spacing))

    @property
    def dynamical_fields(self) -> tuple[str, ...]:
        """Field names with ``time_derivative_order >= 2`` (have momenta)."""
        return tuple(
            eq.field_name
            for eq in self.spec.equations
            if eq.time_derivative_order >= 2  # noqa: PLR2004
        )

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_result(
        cls,
        result: SolverResult,
        spec: EquationSystem,
        grid_info: GridInfo,
        parameters: dict[str, float] | None = None,
        dt: float | None = None,
    ) -> SimulationData:
        """Build from a solver result dict (IDA/leapfrog output).

        This is the native-path constructor — no py-pde types involved.
        Directly slices the flat state vector using ``StateLayout``.

        Parameters
        ----------
        result : dict
            Solver output with keys ``"t"`` (1D times array) and
            ``"y"`` (2D array of shape ``(n_snapshots, total_flat_size)``).
        spec : EquationSystem
            Equation specification.
        grid_info : GridInfo
            Spatial grid descriptor.
        parameters : dict, optional
            Resolved parameter values.
        dt : float, optional
            Time-step size used by the solver (for conservation diagnostics).

        Raises
        ------
        ValueError
            If *result* has no snapshots or flat vector size doesn't match
            the layout.
        """
        if not result["success"]:
            msg = f"Cannot build SimulationData from failed solver result: {result['message']}"
            raise ValueError(msg)

        times = np.asarray(result["t"], dtype=np.float64)
        y_all = np.asarray(result["y"], dtype=np.float64)

        if len(times) == 0:
            msg = "Solver result contains no snapshots"
            raise ValueError(msg)

        layout = StateLayout.from_spec(spec, grid_info.num_points)

        if y_all.ndim == 1:
            y_all = y_all.reshape(1, -1)

        if y_all.shape[1] != layout.total_size:
            msg = (
                f"Flat vector size {y_all.shape[1]} doesn't match "
                f"layout.total_size {layout.total_size}"
            )
            raise ValueError(msg)

        if y_all.shape[0] != len(times):
            msg = (
                f"Snapshot count mismatch: {len(times)} time points but "
                f"{y_all.shape[0]} state vectors in result['y']"
            )
            raise ValueError(msg)

        n_pts = grid_info.num_points
        shape = grid_info.shape

        fields: dict[str, NDArray[np.float64]] = {}
        momenta: dict[str, NDArray[np.float64]] = {}

        for i, slot in enumerate(layout.slots):
            start = i * n_pts
            end = start + n_pts
            # Slice all snapshots at once: (n_snapshots, *grid_shape)
            arr = y_all[:, start:end].reshape(-1, *shape)

            if slot.kind == "velocity":
                momenta[slot.field_name] = arr
            else:
                fields[slot.name] = arr

        return cls(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=grid_info.dx,
            grid_bounds=grid_info.bounds,
            periodic=grid_info.periodic,
            spec=spec,
            parameters=parameters or {},
            bc_types=grid_info.bc_types,
            dt=dt,
        )

    @classmethod
    def from_directory(  # noqa: C901, PLR0912, PLR0914, PLR0915
        cls,
        path: Path | str,
        spec: EquationSystem,
    ) -> SimulationData:
        """Load from a snapshot directory with memory-mapped arrays (O(1) RAM).

        The directory must contain ``metadata.json`` and per-field ``.npy``
        files written by :class:`~tidal.measurement.SnapshotWriter`.  Arrays
        are opened as read-only memory maps — only the pages actually accessed
        by measurement functions are loaded into RAM.

        Parameters
        ----------
        path : Path or str
            Path to the snapshot directory.
        spec : EquationSystem
            JSON-derived equation specification.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist or is not a directory.
        ValueError
            If ``metadata.json`` is missing or corrupt.
        """
        p = Path(path)
        if not p.is_dir():
            msg = f"Snapshot directory not found: {p}"
            raise FileNotFoundError(msg)

        metadata_path = p / "metadata.json"
        metadata: dict[str, object] = {}
        if metadata_path.exists():
            try:
                raw = metadata_path.read_text()
                metadata = json.loads(raw)
            except (json.JSONDecodeError, UnicodeDecodeError):
                # Corrupt metadata (e.g. truncated write) — fall through
                metadata = {}
        if not metadata:
            # Crash recovery: infer snapshot count from times.npy size
            times_path = p / "times.npy"
            if not times_path.exists():
                msg = f"Snapshot directory {p} missing both metadata.json and times.npy"
                raise ValueError(msg)
            times_arr = np.load(str(times_path), mmap_mode="r")
            # Find actual count: last non-zero time (or all if first is 0.0)
            n_recovered = _infer_snapshot_count(times_arr)
            metadata = {
                "n_snapshots": n_recovered,
                "fields": list(spec.component_names),
                "momenta": [
                    eq.field_name
                    for eq in spec.equations
                    if eq.time_derivative_order >= 2  # noqa: PLR2004
                ],
            }

        n = int(metadata["n_snapshots"])  # type: ignore[arg-type]

        # Load times
        times_npy = p / "times.npy"
        if not times_npy.exists():
            msg = f"Snapshot directory {p} missing times.npy"
            raise ValueError(msg)
        times_full = np.load(str(times_npy), mmap_mode="r")
        if n > len(times_full):
            msg = (
                f"Metadata claims {n} snapshots but times.npy "
                f"has only {len(times_full)} entries"
            )
            raise ValueError(msg)
        times = times_full[:n]

        # Load fields (memory-mapped)
        fields: dict[str, NDArray[np.float64]] = {}
        for name in spec.component_names:
            npy = p / f"{name}.npy"
            if npy.exists():
                fields[name] = np.load(str(npy), mmap_mode="r")[:n]

        # Load momenta (memory-mapped)
        momenta: dict[str, NDArray[np.float64]] = {}
        momentum_names = cast("list[str]", metadata.get("momenta", []))
        for name in momentum_names:
            npy = p / f"v_{name}.npy"
            if npy.exists():
                momenta[name] = np.load(str(npy), mmap_mode="r")[:n]

        # Grid metadata — from metadata.json or spec defaults
        grid_spacing: tuple[float, ...]
        if "grid_spacing" in metadata:
            raw_spacing = cast("list[float]", metadata["grid_spacing"])
            grid_spacing = tuple(float(v) for v in raw_spacing)
        else:
            grid_spacing = (1.0,) * spec.spatial_dimension

        grid_bounds: tuple[tuple[float, float], ...]
        if "grid_bounds" in metadata:
            raw_bounds = cast("list[list[float]]", metadata["grid_bounds"])
            grid_bounds = tuple((float(b[0]), float(b[1])) for b in raw_bounds)
        else:
            raw_shape = cast(
                "list[int]",
                metadata.get("grid_shape", [1] * spec.spatial_dimension),
            )
            grid_bounds = tuple(
                (0.0, float(s) * grid_spacing[i]) for i, s in enumerate(raw_shape)
            )

        periodic: tuple[bool, ...]
        if "periodic" in metadata:
            raw_periodic = cast("list[bool]", metadata["periodic"])
            periodic = tuple(bool(v) for v in raw_periodic)
        else:
            periodic = (False,) * spec.spatial_dimension

        # Parameters
        raw_params = cast("dict[str, float]", metadata.get("parameters", {}))
        parameters: dict[str, float] = {str(k): float(v) for k, v in raw_params.items()}

        # BC types (version 2+; None for legacy data)
        bc_types: tuple[str, ...] | None = None
        if "bc_types" in metadata:
            raw_bc = cast("list[str]", metadata["bc_types"])
            bc_types = tuple(str(v) for v in raw_bc)

        # Solver time-step (for conservation diagnostics); None for legacy
        dt_val: float | None = None
        if "dt" in metadata:
            dt_val = float(metadata["dt"])  # type: ignore[arg-type]

        return cls(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=grid_spacing,
            grid_bounds=grid_bounds,
            periodic=periodic,
            spec=spec,
            parameters=parameters,
            bc_types=bc_types,
            dt=dt_val,
        )

    def save(self, path: Path | str) -> Path:
        """Save to a snapshot directory (metadata.json + .npy files).

        This is the inverse of ``load()`` / ``from_directory()``.
        Overwrites any existing files in the directory.

        Parameters
        ----------
        path : Path or str
            Directory to write into (created if it does not exist).

        Returns
        -------
        Path
            The directory that was written.
        """
        p = Path(path)
        p.mkdir(parents=True, exist_ok=True)

        np.save(str(p / "times.npy"), np.asarray(self.times))
        for name, arr in self.fields.items():
            np.save(str(p / f"{name}.npy"), np.asarray(arr))
        for name, arr in self.momenta.items():
            np.save(str(p / f"v_{name}.npy"), np.asarray(arr))

        # Infer grid shape from the first field array
        first_field = next(iter(self.fields.values()))
        grid_shape = list(first_field.shape[1:])  # strip snapshot dim

        metadata: dict[str, object] = {
            "version": 1,
            "n_snapshots": self.n_snapshots,
            "grid_shape": grid_shape,
            "grid_spacing": list(self.grid_spacing),
            "grid_bounds": [list(b) for b in self.grid_bounds],
            "periodic": list(self.periodic),
            "parameters": self.parameters,
            "fields": list(self.fields.keys()),
            "momenta": list(self.momenta.keys()),
            "dtype": "float64",
        }
        if self.bc_types is not None:
            metadata["bc_types"] = list(self.bc_types)
        if self.dt is not None:
            metadata["dt"] = self.dt
        (p / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        return p

    @classmethod
    def load(
        cls,
        path: Path | str,
        spec: EquationSystem,
    ) -> SimulationData:
        """Load from a snapshot directory (memory-mapped, O(1) RAM).

        Parameters
        ----------
        path : Path or str
            Path to snapshot directory.
        spec : EquationSystem
            JSON-derived equation specification.

        Raises
        ------
        ValueError
            If *path* is not a directory.
        """
        p = Path(path)
        if not p.is_dir():
            msg = (
                f"Expected a snapshot directory, got file '{p}'. "
                f"NPZ format is no longer supported — "
                f"use 'tidal simulate --output <directory>'"
            )
            raise ValueError(msg)
        return cls.from_directory(p, spec)


def _infer_snapshot_count(times_arr: NDArray[np.float64]) -> int:
    """Infer actual snapshot count from a times array (crash recovery).

    When ``metadata.json`` is missing (writer wasn't closed), we look at
    the pre-allocated ``times.npy`` and find how many entries were actually
    written.  Unwritten entries are zero (from memmap pre-allocation).

    Strategy: Find the last index where ``times[i] > 0`` or ``i == 0``
    (the initial time t=0 is legitimately zero).
    """
    n = len(times_arr)
    if n == 0:
        return 0
    # The first entry (t=0) is always valid.  After that, unwritten
    # entries are 0.0 from memmap pre-allocation.  Find the last
    # non-zero entry.
    for i in range(n - 1, 0, -1):
        t = float(times_arr[i])
        if t != 0.0 and np.isfinite(t):
            return i + 1
    # Only t=0 was written
    return 1
