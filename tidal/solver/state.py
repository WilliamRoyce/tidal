"""State vector layout and flat ↔ field-dict conversions for TIDAL solvers.

IDA and leapfrog operate on flat numpy arrays.  This module provides:

- ``StateLayout``: describes which slots map to which fields/velocities
- ``state_to_flat`` / ``flat_to_fields``: convert between flat vectors and
  named field dictionaries
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import TYPE_CHECKING

import numpy as np

from tidal.solver._defaults import SECOND_ORDER

if TYPE_CHECKING:
    from tidal.symbolic.json_loader import EquationSystem


# ---------------------------------------------------------------------------
# Slot descriptor
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SlotInfo:
    """Metadata for a single slot in the flat state vector.

    Each slot holds ``grid.num_points`` contiguous values in C-order.

    Attributes
    ----------
    name : str
        Field or velocity name (e.g. ``"phi_0"`` or ``"v_phi_0"``).
    field_name : str
        The physical field this slot belongs to (always the field name,
        even for velocity slots).
    kind : str
        One of ``"field"``, ``"velocity"``, ``"constraint"``.
    time_order : int
        Original LHS time-derivative order (0, 1, or 2).
    dynamical_index : int | None
        Index into the dynamical-field subset (for K matrix lookup).
        None for constraints and first-order fields.
    """

    name: str
    field_name: str
    kind: str
    time_order: int
    dynamical_index: int | None = None


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StateLayout:
    """Describes the mapping between flat state vector and named fields.

    Built from an ``EquationSystem`` spec.  The flat vector is partitioned
    into contiguous blocks of ``num_points`` elements, one per slot.

    Attributes
    ----------
    slots : tuple[SlotInfo, ...]
        Ordered slot descriptors.
    num_points : int
        Grid points per slot (``grid.num_points``).
    field_slot_map : dict[str, int]
        Maps field name → slot index.
    velocity_slot_map : dict[str, int]
        Maps field name → velocity slot index (second-order fields only).
    dynamical_fields : tuple[str, ...]
        Names of dynamical fields (time_order >= 2), in order.
    """

    slots: tuple[SlotInfo, ...]
    num_points: int
    field_slot_map: dict[str, int]
    velocity_slot_map: dict[str, int]
    dynamical_fields: tuple[str, ...]

    @property
    def momentum_slot_map(self) -> dict[str, int]:
        """Alias for ``velocity_slot_map`` (transition aid)."""
        return self.velocity_slot_map

    @classmethod
    def from_spec(cls, spec: EquationSystem, num_points: int) -> StateLayout:
        """Build layout from an equation system specification.

        Follows the same ordering as py-pde FieldCollection:
        for each equation, emit a field slot; if second-order, also emit
        a velocity slot immediately after.
        """
        slots: list[SlotInfo] = []
        field_slot_map: dict[str, int] = {}
        velocity_slot_map: dict[str, int] = {}
        dynamical_fields: list[str] = []
        dyn_idx = 0

        for eq in spec.equations:
            name = eq.field_name
            order = eq.time_derivative_order

            kind = "constraint" if order == 0 else "field"

            d_idx = None
            if order >= SECOND_ORDER:
                d_idx = dyn_idx
                dynamical_fields.append(name)
                dyn_idx += 1

            field_slot_map[name] = len(slots)
            slots.append(
                SlotInfo(
                    name=name,
                    field_name=name,
                    kind=kind,
                    time_order=order,
                    dynamical_index=d_idx,
                ),
            )

            if order >= SECOND_ORDER:
                vel_name = f"v_{name}"
                velocity_slot_map[name] = len(slots)
                slots.append(
                    SlotInfo(
                        name=vel_name,
                        field_name=name,
                        kind="velocity",
                        time_order=order,
                        dynamical_index=d_idx,
                    ),
                )

        return cls(
            slots=tuple(slots),
            num_points=num_points,
            field_slot_map=field_slot_map,
            velocity_slot_map=velocity_slot_map,
            dynamical_fields=tuple(dynamical_fields),
        )

    @cached_property
    def slot_name_to_idx(self) -> dict[str, int]:
        """Map from slot name to slot index. Cached on frozen dataclass."""
        return {slot.name: i for i, slot in enumerate(self.slots)}

    @cached_property
    def total_size(self) -> int:
        """Total flat vector length (num_slots * num_points)."""
        return len(self.slots) * self.num_points

    @cached_property
    def num_slots(self) -> int:
        """Number of slots in the state vector."""
        return len(self.slots)

    def slot_slice(self, slot_idx: int) -> slice:
        """Return the flat-array slice for a given slot index."""
        n = self.num_points
        return slice(slot_idx * n, (slot_idx + 1) * n)

    # ---- Pre-computed slot groups (branch-free hot-path iteration) ----

    @cached_property
    def velocity_slot_groups(self) -> tuple[tuple[int, slice, str], ...]:
        """Pre-computed ``(slot_idx, flat_slice, field_name)`` for velocity slots."""
        return tuple(
            (i, self.slot_slice(i), s.field_name)
            for i, s in enumerate(self.slots)
            if s.kind == "velocity"
        )

    @cached_property
    def dynamical_field_slot_groups(self) -> tuple[tuple[int, slice, int], ...]:
        """Pre-computed ``(slot_idx, flat_slice, vel_slot_idx)`` for 2nd-order fields."""
        return tuple(
            (i, self.slot_slice(i), self.velocity_slot_map[s.field_name])
            for i, s in enumerate(self.slots)
            if s.kind == "field" and s.time_order >= SECOND_ORDER
        )

    @cached_property
    def drift_slot_pairs(self) -> tuple[tuple[slice, slice], ...]:
        """Pre-computed ``(field_slice, vel_slice)`` for zero-copy drift.

        Allows ``y[field_slice] += dt * y[vel_slice]`` without copying
        velocity data into a separate buffer.
        """
        return tuple(
            (self.slot_slice(i), self.slot_slice(self.velocity_slot_map[s.field_name]))
            for i, s in enumerate(self.slots)
            if s.kind == "field" and s.time_order >= SECOND_ORDER
        )

    @cached_property
    def first_order_slot_groups(self) -> tuple[tuple[int, slice, str], ...]:
        """Pre-computed ``(slot_idx, flat_slice, field_name)`` for 1st-order fields."""
        return tuple(
            (i, self.slot_slice(i), s.field_name)
            for i, s in enumerate(self.slots)
            if s.kind == "field" and s.time_order == 1
        )

    @cached_property
    def constraint_slot_groups(self) -> tuple[tuple[int, slice, str], ...]:
        """Pre-computed ``(slot_idx, flat_slice, field_name)`` for constraint slots."""
        return tuple(
            (i, self.slot_slice(i), s.field_name)
            for i, s in enumerate(self.slots)
            if s.time_order == 0
        )

    @cached_property
    def algebraic_indices(self) -> list[int]:
        """Flat indices of algebraic (constraint) variables for IDA.

        IDA needs to know which entries in the state vector are algebraic
        (not differential) so it can handle them appropriately.
        """
        indices: list[int] = []
        n = self.num_points
        for i, slot in enumerate(self.slots):
            if slot.time_order == 0:
                indices.extend(range(i * n, (i + 1) * n))
        return indices


# ---------------------------------------------------------------------------
# Flat ↔ field-dict conversions
# ---------------------------------------------------------------------------


def state_to_flat(
    fields: dict[str, np.ndarray],
    layout: StateLayout,
) -> np.ndarray:
    """Pack named field arrays into a single flat vector.

    Parameters
    ----------
    fields : dict[str, np.ndarray]
        Mapping from slot name (e.g. ``"phi_0"``, ``"v_phi_0"``) to
        grid-shaped arrays.
    layout : StateLayout
        State layout descriptor.

    Returns
    -------
    np.ndarray
        Flat vector of length ``layout.total_size``.

    """
    y = np.zeros(layout.total_size)
    n = layout.num_points
    for i, slot in enumerate(layout.slots):
        y[i * n : (i + 1) * n] = fields[slot.name].ravel()
    return y


def flat_to_fields(
    y: np.ndarray,
    layout: StateLayout,
    shape: tuple[int, ...],
) -> dict[str, np.ndarray]:
    """Unpack a flat vector into named field arrays.

    Parameters
    ----------
    y : np.ndarray
        Flat state vector.
    layout : StateLayout
        State layout descriptor.
    shape : tuple[int, ...]
        Grid shape to reshape each slot's data into.

    Returns
    -------
    dict[str, np.ndarray]
        Mapping from slot name to grid-shaped array.
    """
    fields: dict[str, np.ndarray] = {}
    n = layout.num_points
    for i, slot in enumerate(layout.slots):
        fields[slot.name] = y[i * n : (i + 1) * n].reshape(shape)
    return fields
