"""Tests for tidal.solver.state — state layout and flat ↔ dict conversions."""

from __future__ import annotations

from typing import Any

import numpy as np

from tidal.solver.grid import GridInfo
from tidal.solver.state import StateLayout, flat_to_fields, state_to_flat
from tidal.symbolic.json_loader import EquationSystem


def _make_spec(equations: list[dict[str, Any]], dim: int = 2) -> EquationSystem:
    """Build a minimal EquationSystem from equation dicts."""
    fields = [{"name": eq["field"], "index": i} for i, eq in enumerate(equations)]
    data: dict[str, Any] = {
        "spacetime": {"dimension": dim, "signature": [-1, 1]},
        "fields": fields,
        "equations": equations,
    }
    return EquationSystem.from_dict(data)


def _wave_eq(name: str = "phi_0") -> dict[str, Any]:
    return {
        "field": name,
        "lhs": {"expression": f"d2_t({name})", "order": {"time": 2}},
        "rhs": {
            "type": "linear_combination",
            "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": name}],
        },
    }


def _constraint_eq(name: str = "A_0") -> dict[str, Any]:
    return {
        "field": name,
        "lhs": {"expression": name, "order": {"time": 0, "space": 0}},
        "rhs": {
            "type": "linear_combination",
            "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": name}],
        },
    }


def _first_order_eq(name: str = "psi_0") -> dict[str, Any]:
    return {
        "field": name,
        "lhs": {"expression": f"d_t({name})", "order": {"time": 1}},
        "rhs": {
            "type": "linear_combination",
            "terms": [{"coefficient": 1.0, "operator": "laplacian", "field": name}],
        },
    }


class TestStateLayout:
    def test_single_wave(self) -> None:
        spec = _make_spec([_wave_eq()])
        layout = StateLayout.from_spec(spec, num_points=16)
        assert layout.num_slots == 2  # field + momentum
        assert layout.total_size == 32
        assert layout.field_slot_map == {"phi_0": 0}
        assert layout.momentum_slot_map == {"phi_0": 1}
        assert layout.dynamical_fields == ("phi_0",)
        assert layout.algebraic_indices == []

    def test_constraint_plus_wave(self) -> None:
        spec = _make_spec([_constraint_eq("A_0"), _wave_eq("A_1")])
        layout = StateLayout.from_spec(spec, num_points=8)
        assert layout.num_slots == 3  # constraint + field + momentum
        assert layout.field_slot_map == {"A_0": 0, "A_1": 1}
        assert layout.momentum_slot_map == {"A_1": 2}
        assert layout.dynamical_fields == ("A_1",)
        assert layout.algebraic_indices == list(range(8))  # constraint slot

    def test_mixed_orders(self) -> None:
        spec = _make_spec(
            [
                _constraint_eq("A_0"),
                _first_order_eq("psi_0"),
                _wave_eq("phi_0"),
            ]
        )
        layout = StateLayout.from_spec(spec, num_points=4)
        assert layout.num_slots == 4  # constraint + first_order + field + velocity
        assert layout.slots[0].kind == "constraint"
        assert layout.slots[0].time_order == 0
        assert layout.slots[1].kind == "field"
        assert layout.slots[1].time_order == 1
        assert layout.slots[2].kind == "field"
        assert layout.slots[2].time_order == 2
        assert layout.slots[3].kind == "velocity"

    def test_two_coupled_waves(self) -> None:
        spec = _make_spec([_wave_eq("phi_0"), _wave_eq("chi_0")])
        layout = StateLayout.from_spec(spec, num_points=4)
        assert layout.num_slots == 4
        assert layout.dynamical_fields == ("phi_0", "chi_0")
        assert layout.slots[0].dynamical_index == 0
        assert layout.slots[2].dynamical_index == 1


class TestFlatConversions:
    def test_round_trip(self) -> None:
        spec = _make_spec([_wave_eq()])
        grid = GridInfo(bounds=((0, 1),), shape=(8,), periodic=(True,))
        layout = StateLayout.from_spec(spec, grid.num_points)

        fields = {
            "phi_0": np.sin(np.linspace(0, 2 * np.pi, 8)),
            "v_phi_0": np.zeros(8),
        }
        flat = state_to_flat(fields, layout)
        assert flat.shape == (16,)

        recovered = flat_to_fields(flat, layout, grid.shape)
        np.testing.assert_allclose(recovered["phi_0"], fields["phi_0"])
        np.testing.assert_allclose(recovered["v_phi_0"], fields["v_phi_0"])

    def test_2d_round_trip(self) -> None:
        spec = _make_spec([_wave_eq()], dim=3)
        grid = GridInfo(bounds=((0, 1), (0, 1)), shape=(4, 4), periodic=(True, True))
        layout = StateLayout.from_spec(spec, grid.num_points)

        xs, ys = grid.coord_arrays()
        fields = {
            "phi_0": np.sin(xs) * np.cos(ys),
            "v_phi_0": np.zeros((4, 4)),
        }
        flat = state_to_flat(fields, layout)
        assert flat.shape == (32,)

        recovered = flat_to_fields(flat, layout, grid.shape)
        np.testing.assert_allclose(recovered["phi_0"], fields["phi_0"])
