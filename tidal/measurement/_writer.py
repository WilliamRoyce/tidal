"""Streaming snapshot writer for disk-backed simulation storage.

Writes simulation snapshots incrementally to a directory of ``.npy`` files
using memory-mapped I/O, consuming O(1) memory regardless of snapshot count.

The directory layout is framework-agnostic — any solver (py-pde, Julia, etc.)
can write the same format, and :meth:`SimulationData.from_directory` reads it
with lazy memory-mapped access.

Directory structure::

    output_dir/
      metadata.json          # Grid info, parameters, field list, snapshot count
      times.npy              # shape (n_snapshots,), float64
      phi_0.npy              # shape (n_snapshots, *grid_shape), float64
      v_phi_0.npy            # shape (n_snapshots, *grid_shape), float64
      ...
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Self

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from tidal.symbolic.json_loader import EquationSystem

# Current format version.  Increment when the metadata schema changes.
_FORMAT_VERSION = 1

_DEFAULT_DTYPE: np.dtype[np.float64] = np.dtype(np.float64)

# Flush mmaps every N snapshots for crash resilience.
# Higher values improve I/O performance; lower values reduce data loss on crash.
_DEFAULT_FLUSH_INTERVAL = 10


def _open_memmap(
    path: Path,
    shape: tuple[int, ...],
    dtype: np.dtype[np.float64] = _DEFAULT_DTYPE,
) -> np.memmap[tuple[int, ...], np.dtype[np.float64]]:
    """Create a writable memory-mapped ``.npy`` file with the given shape."""
    # np.lib.format.open_memmap creates a proper .npy file (with header)
    # that can later be loaded with np.load(..., mmap_mode="r").
    return np.lib.format.open_memmap(  # type: ignore[no-any-return]
        str(path),
        mode="w+",
        dtype=dtype,
        shape=shape,
    )


class SnapshotWriter:
    """Stream simulation snapshots to disk with O(1) memory.

    Pre-allocates one ``.npy`` file per field/velocity (plus ``times.npy``)
    using ``numpy.memmap``, then writes each snapshot in-place.  The exact
    number of snapshots must be known at construction time — compute it as
    ``int(t_end / snapshot_interval) + 1``.

    Use as a context manager for automatic ``close()``::

        with SnapshotWriter(output_dir, ...) as writer:
            for t, fields, velocities in simulation:
                writer.append(t, fields, velocities)

    Parameters
    ----------
    output_dir : Path
        Directory to create.  Must not already exist.
    field_names : list[str]
        Names of field arrays to store (e.g. ``["phi_0", "chi_0"]``).
    velocity_names : list[str]
        Names of velocity arrays to store (e.g. ``["phi_0", "chi_0"]``).
        Stored as ``v_{name}.npy``.
    grid_shape : tuple[int, ...]
        Spatial grid shape (e.g. ``(96, 96)``).
    n_snapshots : int
        Exact number of snapshots to write.
    grid_spacing : tuple[float, ...]
        Cell size per spatial axis.
    grid_bounds : tuple[tuple[float, float], ...]
        Domain bounds per spatial axis.
    periodic : tuple[bool, ...]
        Whether each spatial axis is periodic.
    parameters : dict[str, float] or None
        Resolved parameter values.
    spec_path : Path or None
        Path to the JSON spec file (for auto-discovery by ``tidal measure``).
    """

    def __init__(  # noqa: PLR0913, PLR0917
        self,
        output_dir: Path,
        field_names: list[str],
        velocity_names: list[str],
        grid_shape: tuple[int, ...],
        n_snapshots: int,
        grid_spacing: tuple[float, ...],
        grid_bounds: tuple[tuple[float, float], ...],
        periodic: tuple[bool, ...],
        parameters: dict[str, float] | None = None,
        spec_path: Path | None = None,
        flush_interval: int = _DEFAULT_FLUSH_INTERVAL,
        bc_types: tuple[str, ...] | None = None,
        dt: float | None = None,
    ) -> None:
        if n_snapshots < 1:
            msg = f"n_snapshots must be >= 1, got {n_snapshots}"
            raise ValueError(msg)
        if not field_names:
            msg = "field_names must be non-empty"
            raise ValueError(msg)

        self._output_dir = Path(output_dir)
        if self._output_dir.exists() and not self._output_dir.is_dir():
            msg = f"Output path exists but is not a directory: {self._output_dir}"
            raise ValueError(msg)
        self._output_dir.mkdir(parents=True, exist_ok=True)

        # Remove stale .npy / metadata files from a previous run
        existing_npy = list(self._output_dir.glob("*.npy"))
        stale_meta = self._output_dir / "metadata.json"
        if existing_npy or stale_meta.exists():
            for f in existing_npy:
                f.unlink()
            if stale_meta.exists():
                stale_meta.unlink()

        self._n_snapshots = n_snapshots
        self._grid_shape = grid_shape
        self._field_names = list(field_names)
        self._velocity_names = list(velocity_names)
        self._grid_spacing = grid_spacing
        self._grid_bounds = grid_bounds
        self._periodic = periodic
        self._parameters = parameters or {}
        self._spec_path = spec_path
        self._bc_types = bc_types
        self._dt = dt
        self._count = 0
        self._closed = False
        self._flush_interval = max(flush_interval, 1)

        # Pre-allocate times
        self._times_mmap = _open_memmap(
            self._output_dir / "times.npy",
            shape=(n_snapshots,),
        )

        # Pre-allocate per-field .npy files
        snapshot_shape = (n_snapshots, *grid_shape)
        self._field_mmaps: dict[
            str, np.memmap[tuple[int, ...], np.dtype[np.float64]]
        ] = {}
        for name in field_names:
            self._field_mmaps[name] = _open_memmap(
                self._output_dir / f"{name}.npy",
                shape=snapshot_shape,
            )

        # Pre-allocate per-velocity .npy files
        self._velocity_mmaps: dict[
            str, np.memmap[tuple[int, ...], np.dtype[np.float64]]
        ] = {}
        for name in velocity_names:
            self._velocity_mmaps[name] = _open_memmap(
                self._output_dir / f"v_{name}.npy",
                shape=snapshot_shape,
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def append(  # noqa: C901
        self,
        t: float,
        fields: dict[str, NDArray[np.float64]],
        velocities: dict[str, NDArray[np.float64]],
    ) -> None:
        """Write one snapshot at the next time index.

        Parameters
        ----------
        t : float
            Simulation time for this snapshot.
        fields : dict[str, ndarray]
            Mapping ``field_name → spatial_array`` for this snapshot.
        velocities : dict[str, ndarray]
            Mapping ``field_name → spatial_array`` for velocities.

        Raises
        ------
        ValueError
            If the writer is closed, the snapshot count is exceeded,
            the time is non-finite or non-monotonic, or a field/velocity
            array has the wrong shape.
        """
        if self._closed:
            msg = "SnapshotWriter is already closed"
            raise ValueError(msg)
        if self._count >= self._n_snapshots:
            msg = (
                f"Cannot write snapshot {self._count}: "
                f"pre-allocated for {self._n_snapshots} snapshots"
            )
            raise ValueError(msg)

        # Validate time
        if not np.isfinite(t):
            msg = f"Time must be finite, got {t}"
            raise ValueError(msg)
        if self._count > 0:
            prev_t = float(self._times_mmap[self._count - 1])
            if t < prev_t:
                msg = (
                    f"Times must be non-decreasing: "
                    f"snapshot {self._count - 1} at t={prev_t}, "
                    f"snapshot {self._count} at t={t}"
                )
                raise ValueError(msg)

        idx = self._count
        self._times_mmap[idx] = t

        for name in self._field_names:
            if name not in fields:
                msg = f"Missing field '{name}' in snapshot {idx}"
                raise ValueError(msg)
            arr = fields[name]
            if arr.shape != self._grid_shape:
                msg = (
                    f"Field '{name}' has shape {arr.shape}, expected {self._grid_shape}"
                )
                raise ValueError(msg)
            self._field_mmaps[name][idx] = arr

        for name in self._velocity_names:
            if name not in velocities:
                msg = f"Missing velocity '{name}' in snapshot {idx}"
                raise ValueError(msg)
            arr = velocities[name]
            if arr.shape != self._grid_shape:
                msg = (
                    f"Velocity '{name}' has shape {arr.shape}, "
                    f"expected {self._grid_shape}"
                )
                raise ValueError(msg)
            self._velocity_mmaps[name][idx] = arr

        self._count += 1

        # Flush periodically for crash resilience (every _flush_interval snapshots)
        if self._count % self._flush_interval == 0:
            self._flush_mmaps()

    @property
    def count(self) -> int:
        """Number of snapshots written so far."""
        return self._count

    @property
    def n_snapshots(self) -> int:
        """Total number of snapshots pre-allocated."""
        return self._n_snapshots

    @property
    def output_dir(self) -> Path:
        """Directory where snapshot files are written."""
        return self._output_dir

    def close(self) -> None:
        """Flush all mmaps and write ``metadata.json``.

        It is safe to call ``close()`` multiple times.
        """
        if self._closed:
            return

        if self._count != self._n_snapshots:
            warnings.warn(
                f"SnapshotWriter closed with {self._count}/{self._n_snapshots} "
                f"snapshots written",
                stacklevel=2,
            )

        # Flush and release mmaps
        self._flush_mmaps()
        del self._times_mmap
        for mmap in self._field_mmaps.values():
            del mmap
        self._field_mmaps.clear()
        for mmap in self._velocity_mmaps.values():
            del mmap
        self._velocity_mmaps.clear()

        # Write metadata
        metadata: dict[str, Any] = {
            "version": _FORMAT_VERSION,
            "n_snapshots": self._count,
            "grid_spacing": list(self._grid_spacing),
            "grid_bounds": [list(b) for b in self._grid_bounds],
            "grid_shape": list(self._grid_shape),
            "periodic": list(self._periodic),
            "parameters": self._parameters,
            "fields": self._field_names,
            "velocities": self._velocity_names,
            "momenta": self._velocity_names,  # backward compat
            "dtype": "float64",
        }
        if self._bc_types is not None:
            metadata["bc_types"] = list(self._bc_types)
        if self._dt is not None:
            metadata["dt"] = self._dt
        if self._spec_path is not None:
            metadata["spec_path"] = str(self._spec_path)

        metadata_path = self._output_dir / "metadata.json"
        # Atomic write: temp file + rename to avoid corrupt JSON on crash
        temp_path = self._output_dir / "metadata.json.tmp"
        temp_path.write_text(json.dumps(metadata, indent=2) + "\n")
        temp_path.replace(metadata_path)

        self._closed = True

    # ------------------------------------------------------------------
    # Context manager
    # ------------------------------------------------------------------

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _flush_mmaps(self) -> None:
        """Flush all memory-mapped files to disk."""
        self._times_mmap.flush()
        for mmap in self._field_mmaps.values():
            mmap.flush()
        for mmap in self._velocity_mmaps.values():
            mmap.flush()


def compute_snapshot_count(t_end: float, snapshot_interval: float) -> int:
    """Compute the exact number of snapshots for a simulation.

    Parameters
    ----------
    t_end : float
        Total simulation time.
    snapshot_interval : float
        Time between snapshots.

    Returns
    -------
    int
        Exact snapshot count (includes the initial state at t=0).

    Raises
    ------
    ValueError
        If *t_end* or *snapshot_interval* are non-positive.
    """
    if t_end <= 0:
        msg = f"t_end must be positive, got {t_end}"
        raise ValueError(msg)
    if snapshot_interval <= 0:
        msg = f"snapshot_interval must be positive, got {snapshot_interval}"
        raise ValueError(msg)
    return int(t_end / snapshot_interval) + 1


def _field_names_from_spec(
    spec: EquationSystem,
    *,
    exclude_constraints: bool = False,
) -> tuple[list[str], list[str]]:
    """Extract field and velocity names from an equation spec.

    Parameters
    ----------
    spec : EquationSystem
        The equation specification.
    exclude_constraints : bool
        If True, omit fields with ``time_derivative_order == 0``.

    Returns
    -------
    field_names : list[str]
        Field names to store.
    velocity_names : list[str]
        Velocity names to store (fields with ``time_derivative_order >= 2``).
    """
    field_names: list[str] = []
    velocity_names: list[str] = []
    for eq in spec.equations:
        if exclude_constraints and eq.time_derivative_order == 0:
            continue
        field_names.append(eq.field_name)
        if eq.time_derivative_order >= 2:  # noqa: PLR2004
            velocity_names.append(eq.field_name)
    return field_names, velocity_names


def create_snapshot_callback(  # noqa: PLR0913, PLR0917
    output_dir: Path | str,
    spec: EquationSystem,
    grid: Any,  # noqa: ANN401
    t_end: float,
    snapshot_interval: float,
    parameters: dict[str, float] | None = None,
    spec_path: Path | None = None,
) -> tuple[SnapshotWriter, Callable[[Any, float], None]]:
    """Create a SnapshotWriter + callback for streaming snapshots to disk.

    The returned *callback* accepts ``(state, time)`` where *state* is a
    ``pde.FieldCollection``.  Wrap it in ``CallbackTracker(callback, ...)``
    and pass to ``pde.solve()``.  Call ``writer.close()`` after the solve.

    Parameters
    ----------
    output_dir : Path or str
        Directory to write snapshot files into.
    spec : EquationSystem
        Equation system (provides field/velocity names and state layout).
    grid : CartesianGrid
        py-pde grid (provides shape, spacing, bounds, periodicity).
    t_end : float
        Total simulation time.
    snapshot_interval : float
        Time between snapshots.
    parameters : dict or None
        Resolved parameter values to store in metadata.
    spec_path : Path or None
        Path to the JSON spec file (for auto-discovery by ``tidal measure``).

    Returns
    -------
    writer : SnapshotWriter
        Call ``writer.close()`` after the solve finishes.
    callback : callable
        Pass to ``CallbackTracker(callback, interrupts=snapshot_interval)``.
    """
    out = Path(output_dir)
    field_names, velocity_names = _field_names_from_spec(spec)
    n_snapshots = compute_snapshot_count(t_end, snapshot_interval)

    grid_bounds_raw = grid.bounds
    spacing = tuple(
        float((b[1] - b[0]) / s)
        for b, s in zip(grid_bounds_raw, grid.shape, strict=True)
    )
    bounds = tuple((float(b[0]), float(b[1])) for b in grid_bounds_raw)
    periodic_flags = tuple(bool(p) for p in grid.periodic)

    writer = SnapshotWriter(
        output_dir=out,
        field_names=field_names,
        velocity_names=velocity_names,
        grid_shape=tuple(grid.shape),
        n_snapshots=n_snapshots,
        grid_spacing=spacing,
        grid_bounds=bounds,
        periodic=periodic_flags,
        parameters=parameters,
        spec_path=Path(spec_path) if spec_path else None,
    )

    # Build slot maps from spec.state_layout
    field_slots: dict[str, int] = {}
    velocity_slots: dict[str, int] = {}
    for idx, (name, slot_type) in enumerate(spec.state_layout):
        if slot_type == "field":
            field_slots[name] = idx
        elif slot_type == "velocity":
            velocity_slots[name] = idx

    field_set = set(field_names)
    velocity_set = set(velocity_names)

    def _on_snapshot(state_view: Any, time: float) -> None:  # noqa: ANN401
        fields = {
            name: np.asarray(state_view[slot].data, dtype=np.float64)
            for name, slot in field_slots.items()
            if name in field_set
        }
        vels = {
            name: np.asarray(state_view[slot].data, dtype=np.float64)
            for name, slot in velocity_slots.items()
            if name in velocity_set
        }
        writer.append(time, fields, vels)

    return writer, _on_snapshot
