"""Typed field container for TIDAL solvers.

``FieldSet`` owns named field data on a grid, backed by a single contiguous
flat numpy array.  Named access returns zero-copy views (numpy slices).

This consolidates the state packing/unpacking logic previously duplicated
across ida.py, leapfrog.py, _simulate.py, _io.py, and _writer.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from typing_extensions import override

if TYPE_CHECKING:
    from tidal.solver.state import StateLayout


class FieldSet:
    """Typed container owning named field data on a grid.

    Backed by a single contiguous flat ``np.ndarray``.  Named access
    returns zero-copy views (numpy slices) into this array.

    Parameters
    ----------
    layout : StateLayout
        Describes the slot structure (field/momentum names and ordering).
    grid_shape : tuple[int, ...]
        Shape of the spatial grid (e.g. ``(64,)`` or ``(32, 32)``).
    data : np.ndarray, optional
        Flat data array of length ``layout.total_size``.
        If ``None``, initializes to zeros.
    """

    __slots__ = ("_data", "_grid_shape", "_layout", "_n", "_name_to_idx")

    def __init__(
        self,
        layout: StateLayout,
        grid_shape: tuple[int, ...],
        data: np.ndarray | None = None,
    ) -> None:
        self._layout = layout
        self._grid_shape = grid_shape
        self._n = layout.num_points

        if data is not None:
            if data.shape != (layout.total_size,):
                msg = (
                    f"Expected flat array of length {layout.total_size}, "
                    f"got shape {data.shape}"
                )
                raise ValueError(msg)
            self._data = data
        else:
            self._data = np.zeros(layout.total_size)

        # Build unified name → slot index mapping
        self._name_to_idx: dict[str, int] = {}
        for i, slot in enumerate(layout.slots):
            self._name_to_idx[slot.name] = i

    # ---- Named access (zero-copy views) ----

    def __getitem__(self, name: str) -> np.ndarray:
        """Return grid-shaped view of named slot.

        Raises
        ------
        KeyError
            If *name* is not a valid slot name.
        """
        idx = self._name_to_idx.get(name)
        if idx is None:
            msg = f"Unknown slot '{name}'. Valid slots: {sorted(self._name_to_idx)}"
            raise KeyError(msg)
        start = idx * self._n
        return self._data[start : start + self._n].reshape(self._grid_shape)

    def __setitem__(self, name: str, value: np.ndarray) -> None:
        """Write grid-shaped data into named slot.

        Raises
        ------
        KeyError
            If *name* is not a valid slot name.
        """
        idx = self._name_to_idx.get(name)
        if idx is None:
            msg = f"Unknown slot '{name}'. Valid slots: {sorted(self._name_to_idx)}"
            raise KeyError(msg)
        start = idx * self._n
        self._data[start : start + self._n] = np.asarray(value).ravel()

    def __contains__(self, name: object) -> bool:
        """Check if *name* is a valid slot name."""
        return name in self._name_to_idx

    def __len__(self) -> int:
        """Return the number of slots."""
        return len(self._name_to_idx)

    @override
    def __repr__(self) -> str:
        slots = ", ".join(self._name_to_idx)
        return f"FieldSet(slots=[{slots}], grid_shape={self._grid_shape})"

    # ---- Flat array access ----

    @property
    def flat(self) -> np.ndarray:
        """Underlying flat vector (no copy). Used by IDA/leapfrog."""
        return self._data

    # ---- Metadata ----

    @property
    def layout(self) -> StateLayout:
        """The state layout this FieldSet wraps."""
        return self._layout

    @property
    def grid_shape(self) -> tuple[int, ...]:
        """Spatial grid shape."""
        return self._grid_shape

    @property
    def field_names(self) -> tuple[str, ...]:
        """Ordered field slot names (e.g. ``("phi_0", "A_0")``).

        Excludes momentum slots.
        """
        return tuple(
            slot.name for slot in self._layout.slots if slot.kind != "momentum"
        )

    @property
    def momentum_names(self) -> tuple[str, ...]:
        """Ordered momentum slot names (e.g. ``("pi_phi_0",)``).

        Only present for second-order equations.
        """
        return tuple(
            slot.name for slot in self._layout.slots if slot.kind == "momentum"
        )

    @property
    def slot_names(self) -> tuple[str, ...]:
        """All slot names in order."""
        return tuple(slot.name for slot in self._layout.slots)

    # ---- Dict-like views (zero-copy) ----

    def fields_dict(self) -> dict[str, np.ndarray]:
        """Return dict of field name → grid-shaped view (zero-copy)."""
        return {name: self[name] for name in self.field_names}

    def momenta_dict(self) -> dict[str, np.ndarray]:
        """Return dict of momentum name → grid-shaped view (zero-copy)."""
        return {name: self[name] for name in self.momentum_names}

    def as_dict(self) -> dict[str, np.ndarray]:
        """Return dict of all slot names → grid-shaped views (zero-copy)."""
        return {name: self[name] for name in self.slot_names}

    # ---- Constructors ----

    def copy(self) -> FieldSet:
        """Deep copy (new flat array)."""
        return FieldSet(self._layout, self._grid_shape, self._data.copy())

    @classmethod
    def zeros(cls, layout: StateLayout, grid_shape: tuple[int, ...]) -> FieldSet:
        """Create a FieldSet initialized to zero."""
        return cls(layout, grid_shape)

    @classmethod
    def from_flat(
        cls,
        layout: StateLayout,
        grid_shape: tuple[int, ...],
        flat: np.ndarray,
    ) -> FieldSet:
        """Wrap an existing flat array (no copy).

        The caller must ensure the array is not unexpectedly mutated
        elsewhere.
        """
        return cls(layout, grid_shape, flat)

    @classmethod
    def from_dict(
        cls,
        layout: StateLayout,
        grid_shape: tuple[int, ...],
        slot_data: dict[str, np.ndarray],
    ) -> FieldSet:
        """Pack named arrays into a FieldSet.

        Missing slots default to zero.  Extra keys not in the layout
        are silently ignored.

        Parameters
        ----------
        layout : StateLayout
            State layout descriptor.
        grid_shape : tuple[int, ...]
            Spatial grid shape.
        slot_data : dict[str, np.ndarray]
            Mapping from slot name → grid-shaped array.
        """
        fs = cls(layout, grid_shape)
        for name, arr in slot_data.items():
            if name in fs:
                fs[name] = arr
        return fs
