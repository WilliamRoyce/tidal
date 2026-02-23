"""Tests for tidal.solver.fields — FieldSet typed container."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tidal.solver.fields import FieldSet
from tidal.solver.state import StateLayout
from tidal.symbolic.json_loader import EquationSystem

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_kg_spec() -> EquationSystem:
    """Klein-Gordon 1D: second-order → 2 slots (phi_0, pi_phi_0)."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 2, "signature": [-1, 1]},
        "fields": [{"name": "phi_0", "index": 0}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            }
        ],
    }
    return EquationSystem.from_dict(data)


def _make_coupled_spec() -> EquationSystem:
    """Coupled 2-field system: phi_0 (2nd order) + A_0 (constraint)."""
    data: dict[str, Any] = {
        "spacetime": {"dimension": 3, "signature": [-1, 1, 1]},
        "fields": [
            {"name": "A_0", "index": 0},
            {"name": "phi_0", "index": 1},
        ],
        "equations": [
            {
                "field": "A_0",
                "lhs": {"expression": "0", "order": {"time": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "A_0"},
                    ],
                },
            },
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi_0"},
                    ],
                },
            },
        ],
    }
    return EquationSystem.from_dict(data)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestFieldSetConstruction:
    def test_zeros(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 32)
        fs = FieldSet.zeros(layout, (32,))

        assert fs.grid_shape == (32,)
        assert fs.flat.shape == (64,)  # 2 slots x 32 points
        np.testing.assert_array_equal(fs.flat, 0.0)

    def test_from_flat_no_copy(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 16)
        flat = np.arange(32, dtype=np.float64)
        fs = FieldSet.from_flat(layout, (16,), flat)

        assert fs.flat is flat  # same object, no copy

    def test_from_dict(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        phi = np.ones(8)
        pi = np.full(8, 2.0)
        fs = FieldSet.from_dict(layout, (8,), {"phi_0": phi, "pi_phi_0": pi})

        np.testing.assert_array_equal(fs["phi_0"], 1.0)
        np.testing.assert_array_equal(fs["pi_phi_0"], 2.0)

    def test_from_dict_partial(self) -> None:
        """Missing slots should default to zero."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.from_dict(layout, (8,), {"phi_0": np.ones(8)})

        np.testing.assert_array_equal(fs["phi_0"], 1.0)
        np.testing.assert_array_equal(fs["pi_phi_0"], 0.0)

    def test_from_dict_extra_keys_ignored(self) -> None:
        """Extra keys not in layout should be silently ignored."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.from_dict(
            layout, (8,), {"phi_0": np.ones(8), "nonexistent": np.zeros(8)}
        )
        np.testing.assert_array_equal(fs["phi_0"], 1.0)

    def test_wrong_flat_size_raises(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        with pytest.raises(ValueError, match="length 16"):
            FieldSet(layout, (8,), data=np.zeros(10))


class TestFieldSetAccess:
    def test_getitem_returns_view(self) -> None:
        """Modifying the view should change the underlying flat array."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        view = fs["phi_0"]
        view[:] = 42.0

        # Flat array should reflect the change
        np.testing.assert_array_equal(fs.flat[:8], 42.0)

    def test_setitem(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        fs["pi_phi_0"] = np.full(8, 2.72)
        np.testing.assert_array_equal(fs["pi_phi_0"], 2.72)

    def test_contains(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        assert "phi_0" in fs
        assert "pi_phi_0" in fs
        assert "nonexistent" not in fs

    def test_unknown_slot_raises_keyerror(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        with pytest.raises(KeyError, match="nonexistent"):
            _ = fs["nonexistent"]

    def test_setitem_unknown_raises_keyerror(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        with pytest.raises(KeyError, match="bad"):
            fs["bad"] = np.zeros(8)

    def test_2d_grid_shape(self) -> None:
        """FieldSet should reshape to 2D grid."""
        spec = _make_coupled_spec()
        n = 4 * 4  # 16 points
        layout = StateLayout.from_spec(spec, n)
        fs = FieldSet.zeros(layout, (4, 4))

        # A_0 is constraint (1 slot), phi_0 is 2nd order (2 slots: field + mom)
        assert fs["A_0"].shape == (4, 4)
        assert fs["phi_0"].shape == (4, 4)
        assert fs["pi_phi_0"].shape == (4, 4)


class TestFieldSetMetadata:
    def test_field_names(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))

        # field_names includes fields and constraints, but NOT momenta
        assert fs.field_names == ("A_0", "phi_0")

    def test_momentum_names(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))

        assert fs.momentum_names == ("pi_phi_0",)

    def test_slot_names(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))

        assert fs.slot_names == ("A_0", "phi_0", "pi_phi_0")

    def test_len(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))
        assert len(fs) == 2  # phi_0 + pi_phi_0

    def test_repr(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))
        r = repr(fs)
        assert "phi_0" in r
        assert "pi_phi_0" in r
        assert "(8,)" in r


class TestFieldSetDicts:
    def test_fields_dict(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))
        fs["phi_0"] = np.ones(4)

        fd = fs.fields_dict()
        assert set(fd.keys()) == {"A_0", "phi_0"}
        np.testing.assert_array_equal(fd["phi_0"], 1.0)

    def test_momenta_dict(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))
        fs["pi_phi_0"] = np.full(4, 5.0)

        md = fs.momenta_dict()
        assert set(md.keys()) == {"pi_phi_0"}
        np.testing.assert_array_equal(md["pi_phi_0"], 5.0)

    def test_as_dict(self) -> None:
        spec = _make_coupled_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet.zeros(layout, (4,))

        d = fs.as_dict()
        assert set(d.keys()) == {"A_0", "phi_0", "pi_phi_0"}

    def test_dict_views_are_zero_copy(self) -> None:
        """Dict values should be views, not copies."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))

        fd = fs.fields_dict()
        fd["phi_0"][:] = 99.0
        np.testing.assert_array_equal(fs["phi_0"], 99.0)


class TestFieldSetCopy:
    def test_copy_is_independent(self) -> None:
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet.zeros(layout, (8,))
        fs["phi_0"] = np.ones(8)

        fs2 = fs.copy()
        fs2["phi_0"] = np.full(8, 42.0)

        # Original unchanged
        np.testing.assert_array_equal(fs["phi_0"], 1.0)
        np.testing.assert_array_equal(fs2["phi_0"], 42.0)


class TestFieldSetRoundTrip:
    def test_from_dict_to_dict(self) -> None:
        """Round-trip: from_dict → fields_dict + momenta_dict → from_dict."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 16)
        original = {
            "phi_0": np.sin(np.linspace(0, 2 * np.pi, 16)),
            "pi_phi_0": np.cos(np.linspace(0, 2 * np.pi, 16)),
        }
        fs1 = FieldSet.from_dict(layout, (16,), original)

        # Extract and reconstruct
        combined = {**fs1.fields_dict(), **fs1.momenta_dict()}
        fs2 = FieldSet.from_dict(layout, (16,), combined)

        np.testing.assert_array_equal(fs1.flat, fs2.flat)

    def test_flat_roundtrip(self) -> None:
        """from_flat → modify via __setitem__ → flat reflects changes."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        flat = np.zeros(16)
        fs = FieldSet.from_flat(layout, (8,), flat)

        fs["phi_0"] = np.arange(8, dtype=np.float64)
        np.testing.assert_array_equal(flat[:8], np.arange(8))


class TestFieldSetDiagnostics:
    def test_max_norm(self) -> None:
        """max_norm returns the largest absolute value across all slots."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet(layout, (8,))
        fs["phi_0"] = np.array([1.0, -3.0, 2.0, 0.0, 0.5, -0.5, 1.5, -2.5])
        fs["pi_phi_0"] = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0])
        assert fs.max_norm() == pytest.approx(3.0)

    def test_max_norm_all_zeros(self) -> None:
        """max_norm returns 0.0 when all fields are zero."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 4)
        fs = FieldSet(layout, (4,))
        assert fs.max_norm() == pytest.approx(0.0)

    def test_check_finite_true(self) -> None:
        """check_finite returns True when all values are finite."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet(layout, (8,))
        fs["phi_0"] = np.linspace(-1.0, 1.0, 8)
        assert fs.check_finite() is True

    def test_check_finite_false_nan(self) -> None:
        """check_finite returns False when a NaN is present."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet(layout, (8,))
        arr = np.linspace(-1.0, 1.0, 8)
        arr[3] = float("nan")
        fs["phi_0"] = arr
        assert fs.check_finite() is False

    def test_check_finite_false_inf(self) -> None:
        """check_finite returns False when an Inf is present."""
        spec = _make_kg_spec()
        layout = StateLayout.from_spec(spec, 8)
        fs = FieldSet(layout, (8,))
        arr = np.linspace(-1.0, 1.0, 8)
        arr[0] = float("inf")
        fs["phi_0"] = arr
        assert fs.check_finite() is False
