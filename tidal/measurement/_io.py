"""Uniform data abstraction for simulation output.

Provides ``SimulationData``, a frozen dataclass that stores the full time
history of field and momentum arrays.  Can be constructed from either a live
``MemoryStorage`` object (straight from the solver) or from an NPZ file
previously saved by ``tidal simulate --output``.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import prod
from pathlib import Path
from typing import TYPE_CHECKING, cast

import numpy as np

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from pde import CartesianGrid, FieldCollection, MemoryStorage

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
    """

    times: NDArray[np.float64]
    fields: dict[str, NDArray[np.float64]]
    momenta: dict[str, NDArray[np.float64]]
    grid_spacing: tuple[float, ...]
    grid_bounds: tuple[tuple[float, float], ...]
    periodic: tuple[bool, ...]
    spec: EquationSystem
    parameters: dict[str, float]

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
    def from_storage(
        cls,
        storage: MemoryStorage,
        spec: EquationSystem,
        grid: CartesianGrid,
        parameters: dict[str, float] | None = None,
    ) -> SimulationData:
        """Build from a live ``MemoryStorage`` produced by the solver.

        Iterates through **all** snapshots (not just the final state) so
        that transient dynamics such as Rabi oscillations are captured.

        Parameters
        ----------
        storage : MemoryStorage
            Time-indexed snapshot collection from the solver.
        spec : EquationSystem
            JSON-derived equation specification.
        grid : CartesianGrid
            Spatial grid used for the simulation.
        parameters : dict, optional
            Resolved parameter values (e.g. ``{"m2": 1.0, "g": 0.5}``).

        Raises
        ------
        ValueError
            If *storage* is empty or snapshot field count does not match
            ``spec.state_size``.
        """
        n = len(storage)
        if n == 0:
            msg = "MemoryStorage is empty — no snapshots to extract"
            raise ValueError(msg)

        # Validate first snapshot size
        first_snapshot = cast("FieldCollection", storage[0])
        if len(first_snapshot) != spec.state_size:
            msg = (
                f"Snapshot has {len(first_snapshot)} slots but spec.state_size "
                f"is {spec.state_size}"
            )
            raise ValueError(msg)

        # Build slot maps: field_name → slot_idx, momentum_name → slot_idx
        field_slot_map: dict[str, int] = {}
        momentum_slot_map: dict[str, int] = {}
        for idx, (name, slot_type) in enumerate(spec.state_layout):
            if slot_type == "field":
                field_slot_map[name] = idx
            else:  # "momentum"
                momentum_slot_map[name] = idx

        # Extract full time history
        times = np.array(storage.times, dtype=np.float64)

        fields: dict[str, NDArray[np.float64]] = {}
        for name, slot_idx in field_slot_map.items():
            fields[name] = np.stack(
                [
                    np.asarray(cast("FieldCollection", storage[t])[slot_idx].data, dtype=np.float64)
                    for t in range(n)
                ]
            )

        momenta: dict[str, NDArray[np.float64]] = {}
        for name, slot_idx in momentum_slot_map.items():
            momenta[name] = np.stack(
                [
                    np.asarray(cast("FieldCollection", storage[t])[slot_idx].data, dtype=np.float64)
                    for t in range(n)
                ]
            )

        # Grid metadata
        spacing = tuple(
            float((b[1] - b[0]) / s)
            for b, s in zip(grid.axes_bounds, grid.shape, strict=True)
        )
        bounds = tuple((float(b[0]), float(b[1])) for b in grid.axes_bounds)
        periodic_flags = tuple(bool(p) for p in grid.periodic)

        return cls(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=spacing,
            grid_bounds=bounds,
            periodic=periodic_flags,
            spec=spec,
            parameters=parameters or {},
        )

    @classmethod
    def from_npz(  # noqa: C901, PLR0913, PLR0917
        cls,
        path: Path | str,
        spec: EquationSystem,
        grid_spacing: tuple[float, ...],
        grid_bounds: tuple[tuple[float, float], ...],
        periodic: tuple[bool, ...],
        parameters: dict[str, float] | None = None,
    ) -> SimulationData:
        """Load from an NPZ file saved by ``tidal simulate --output``.

        Parameters
        ----------
        path : Path or str
            Path to the ``.npz`` file.
        spec : EquationSystem
            JSON-derived equation specification.
        grid_spacing : tuple[float, ...]
            Cell sizes per spatial axis.
        grid_bounds : tuple[tuple[float, float], ...]
            Domain bounds per spatial axis.
        periodic : tuple[bool, ...]
            Per-axis periodicity flags.
        parameters : dict, optional
            Resolved parameter values.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the NPZ is missing ``"times"`` or expected field keys.
        """
        p = Path(path)
        if not p.exists():
            msg = f"NPZ file not found: {p}"
            raise FileNotFoundError(msg)

        data = np.load(str(p))

        if "times" not in data:
            msg = f"NPZ file {p} missing required 'times' key"
            raise ValueError(msg)
        times = np.asarray(data["times"], dtype=np.float64)
        n = len(times)

        # Reconstruct field arrays
        fields: dict[str, NDArray[np.float64]] = {}
        for name in spec.component_names:
            slices: list[NDArray[np.float64]] = []
            for t_idx in range(n):
                key = f"{name}_t{t_idx}"
                if key not in data:
                    msg = f"NPZ file {p} missing expected key '{key}'"
                    raise ValueError(msg)
                slices.append(np.asarray(data[key], dtype=np.float64))
            fields[name] = np.stack(slices)

        # Reconstruct momentum arrays (may not be present in older NPZs)
        momenta: dict[str, NDArray[np.float64]] = {}
        for eq in spec.equations:
            if eq.time_derivative_order < 2:  # noqa: PLR2004
                continue
            name = eq.field_name
            key0 = f"pi_{name}_t0"
            if key0 not in data:
                # Older NPZ format without momenta — skip gracefully
                continue
            slices_m: list[NDArray[np.float64]] = []
            for t_idx in range(n):
                key = f"pi_{name}_t{t_idx}"
                if key not in data:
                    msg = f"NPZ file {p} missing expected momentum key '{key}'"
                    raise ValueError(msg)
                slices_m.append(np.asarray(data[key], dtype=np.float64))
            momenta[name] = np.stack(slices_m)

        return cls(
            times=times,
            fields=fields,
            momenta=momenta,
            grid_spacing=grid_spacing,
            grid_bounds=grid_bounds,
            periodic=periodic,
            spec=spec,
            parameters=parameters or {},
        )

    @classmethod
    def from_npz_auto(
        cls,
        path: Path | str,
        spec: EquationSystem,
    ) -> SimulationData:
        """Load from an NPZ file, reading grid metadata and parameters from the file.

        This is a convenience wrapper around :meth:`from_npz` that reads
        ``grid_spacing``, ``grid_bounds``, ``grid_periodic``, and resolved
        parameters directly from the NPZ — no manual specification needed.

        Parameters
        ----------
        path : Path or str
            Path to the ``.npz`` file (saved by ``tidal simulate --output``).
        spec : EquationSystem
            JSON-derived equation specification.

        Raises
        ------
        FileNotFoundError
            If *path* does not exist.
        ValueError
            If the NPZ is missing grid metadata keys.
        """
        p = Path(path)
        if not p.exists():
            msg = f"NPZ file not found: {p}"
            raise FileNotFoundError(msg)

        data = np.load(str(p))

        for key in ("grid_spacing", "grid_bounds", "grid_periodic"):
            if key not in data:
                msg = (
                    f"NPZ file {p} missing '{key}' metadata — "
                    f"use from_npz() with explicit grid parameters instead"
                )
                raise ValueError(msg)

        grid_spacing = tuple(float(v) for v in data["grid_spacing"])
        grid_bounds = tuple(
            (float(row[0]), float(row[1])) for row in data["grid_bounds"]
        )
        periodic = tuple(bool(v) for v in data["grid_periodic"])

        # Reconstruct parameters if saved
        parameters: dict[str, float] = {}
        if "_param_names" in data and "_param_values" in data:
            names = list(data["_param_names"])
            values = list(data["_param_values"])
            parameters = dict(zip(names, (float(v) for v in values), strict=True))

        return cls.from_npz(
            path=p,
            spec=spec,
            grid_spacing=grid_spacing,
            grid_bounds=grid_bounds,
            periodic=periodic,
            parameters=parameters,
        )
