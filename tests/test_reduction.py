"""Tests for plane-wave dimensional reduction (tidal.symbolic.reduction)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from tidal.symbolic.reduction import (
    _build_operator_remap,
    _check_killed_coord_refs,
    _find_degenerate_constraint,
    _find_gradient_zero_constraint,
    _remap_coord_string,
    _remap_operator,
    _substitute_in_terms,
    _transform_hamiltonian_term,
    _transform_term,
    eliminate_degenerate_constraints,
    reduce_spec,
)

# ---------------------------------------------------------------------------
# Fixtures: minimal 3+1D spec for unit tests
# ---------------------------------------------------------------------------

_EXAMPLES_DATA = Path(__file__).resolve().parent.parent / "examples" / "data"


def _make_3d_kg_spec() -> dict[str, Any]:
    """Minimal flat 3+1D Klein-Gordon spec for unit tests."""
    return {
        "metadata": {"source": "test", "parameters": {"m2": 1.0}},
        "spacetime": {
            "dimension": 4,
            "signature": [-1, 1, 1, 1],
            "coordinates": ["t", "x", "y", "z"],
        },
        "fields": [
            {"name": "phi_0", "index": 0, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": -1.0,
                            "operator": "identity",
                            "field": "phi_0",
                            "coefficient_symbolic": "-m2",
                        },
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi_0"},
                        {"coefficient": 1.0, "operator": "laplacian_y", "field": "phi_0"},
                        {"coefficient": 1.0, "operator": "laplacian_z", "field": "phi_0"},
                    ],
                },
            }
        ],
        "coupling": {"mass_matrix_symbolic": [["-m2"]]},
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "identity"},
                    "factor_b": {"field": "phi_0", "operator": "identity"},
                    "coefficient_symbolic": "m2/2",
                },
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "gradient_x"},
                    "factor_b": {"field": "phi_0", "operator": "gradient_x"},
                },
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "gradient_y"},
                    "factor_b": {"field": "phi_0", "operator": "gradient_y"},
                },
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "gradient_z"},
                    "factor_b": {"field": "phi_0", "operator": "gradient_z"},
                },
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                    "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                },
            ],
        },
    }


def _make_2field_3d_spec() -> dict[str, Any]:
    """Two-field 3+1D spec where one field has only killed-axis operators."""
    spec = _make_3d_kg_spec()
    spec["fields"].append({"name": "psi_0", "index": 1, "is_dynamical": True})
    # psi_0 equation only has transverse operators (should be eliminated)
    spec["equations"].append({
        "field": "psi_0",
        "lhs": {"expression": "d2_t(psi_0)", "order": {"time": 2, "space": 0}},
        "rhs": {
            "type": "linear_combination",
            "terms": [
                {"coefficient": 1.0, "operator": "laplacian_x", "field": "psi_0"},
                {"coefficient": 1.0, "operator": "laplacian_y", "field": "psi_0"},
            ],
        },
    })
    spec["coupling"]["mass_matrix_symbolic"].append(["-m2"])
    spec["coupling"]["mass_matrix_symbolic"][0].append("0")
    return spec


def _make_curved_spec_surviving() -> dict[str, Any]:
    """3+1D spec with coefficient_symbolic on the propagation axis only."""
    return {
        "metadata": {"source": "test", "parameters": {}},
        "spacetime": {
            "dimension": 4,
            "signature": [-1, 1, 1, 1],
            "coordinates": ["t", "x", "y", "z"],
        },
        "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "gradient_z",
                            "field": "phi_0",
                            "coefficient_symbolic": "2/z[]",
                            "coordinate_dependent": ["z"],
                        },
                        {"coefficient": 1.0, "operator": "laplacian_z", "field": "phi_0"},
                    ],
                },
            }
        ],
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "gradient_z"},
                    "factor_b": {"field": "phi_0", "operator": "gradient_z"},
                },
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "phi_0", "operator": "time_derivative"},
                    "factor_b": {"field": "phi_0", "operator": "time_derivative"},
                },
            ],
            "volume_element": "z[]^2",
        },
    }


def _make_curved_spec_killed_ref() -> dict[str, Any]:
    """Spec where a surviving term references a killed coordinate (incompatible)."""
    return {
        "metadata": {"source": "test", "parameters": {}},
        "spacetime": {
            "dimension": 3,
            "signature": [-1, 1, 1],
            "coordinates": ["t", "x", "y"],
        },
        "fields": [{"name": "phi_0", "index": 0, "is_dynamical": True}],
        "equations": [
            {
                "field": "phi_0",
                "lhs": {"expression": "d2_t(phi_0)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {
                            "coefficient": 1.0,
                            "operator": "laplacian_x",
                            "field": "phi_0",
                            "coefficient_symbolic": "y[]^2/(x[]^2 + y[]^2)",
                            "coordinate_dependent": ["x", "y"],
                        },
                    ],
                },
            }
        ],
        "canonical": {"hamiltonian_terms": []},
    }


# ---------------------------------------------------------------------------
# Unit tests: operator remap
# ---------------------------------------------------------------------------


class TestBuildOperatorRemap:
    """Tests for _build_operator_remap."""

    def test_3d_propagation_z(self) -> None:
        """Propagation along z: laplacian_z->laplacian_x, killed x/y."""
        remap = _build_operator_remap("z", ["x", "y", "z"])
        assert remap["laplacian_z"] == "laplacian_x"
        assert remap["gradient_z"] == "gradient_x"
        assert remap["laplacian_x"] is None
        assert remap["laplacian_y"] is None
        assert remap["gradient_x"] is None
        assert remap["gradient_y"] is None

    def test_3d_propagation_x(self) -> None:
        """Propagation along x: laplacian_x->laplacian_x, killed y/z."""
        remap = _build_operator_remap("x", ["x", "y", "z"])
        assert remap["laplacian_x"] == "laplacian_x"
        assert remap["gradient_x"] == "gradient_x"
        assert remap["laplacian_y"] is None
        assert remap["laplacian_z"] is None

    def test_2d_propagation_x(self) -> None:
        """2+1D propagation along x: killed y."""
        remap = _build_operator_remap("x", ["x", "y"])
        assert remap["laplacian_x"] == "laplacian_x"
        assert remap["laplacian_y"] is None

    def test_cross_derivatives_all_removed(self) -> None:
        """All cross derivatives are removed."""
        remap = _build_operator_remap("z", ["x", "y", "z"])
        assert remap["cross_derivative_xy"] is None
        assert remap["cross_derivative_xz"] is None
        assert remap["cross_derivative_yz"] is None

    def test_identity_preserved(self) -> None:
        """Identity and first_derivative_t are preserved."""
        remap = _build_operator_remap("z", ["x", "y", "z"])
        assert remap["identity"] == "identity"
        assert remap["first_derivative_t"] == "first_derivative_t"
        assert remap["time_derivative"] == "time_derivative"
        assert remap["laplacian"] == "laplacian"
        assert remap["biharmonic"] == "biharmonic"


# ---------------------------------------------------------------------------
# Unit tests: coordinate string remapping
# ---------------------------------------------------------------------------


class TestRemapCoordString:
    """Tests for _remap_coord_string."""

    def test_z_to_x(self) -> None:
        """z[] should be remapped to x[]."""
        assert _remap_coord_string("2/z[]", "z") == "2/x[]"

    def test_x_unchanged(self) -> None:
        """If propagation axis is already x, no change."""
        assert _remap_coord_string("2/x[]", "x") == "2/x[]"

    def test_no_coord_refs(self) -> None:
        """No coordinate references → unchanged."""
        assert _remap_coord_string("m2/2", "z") == "m2/2"

    def test_multiple_refs(self) -> None:
        """Multiple z[] references all remapped."""
        assert _remap_coord_string("z[]^2 + z[]", "z") == "x[]^2 + x[]"


# ---------------------------------------------------------------------------
# Unit tests: killed coordinate checking
# ---------------------------------------------------------------------------


class TestCheckKilledCoordRefs:
    """Tests for _check_killed_coord_refs."""

    def test_no_killed_refs_ok(self) -> None:
        """No killed coordinate refs → no error."""
        _check_killed_coord_refs("2/z[]", {"x", "y"}, "test")

    def test_killed_ref_raises(self) -> None:
        """Reference to killed coordinate → ValueError."""
        with pytest.raises(ValueError, match="killed coordinate 'y"):
            _check_killed_coord_refs("y[]^2/(x[]^2 + y[]^2)", {"y"}, "test")

    def test_no_refs_at_all_ok(self) -> None:
        """No coordinate references at all → OK."""
        _check_killed_coord_refs("m2", {"x", "y"}, "test")


# ---------------------------------------------------------------------------
# Unit tests: operator remap lookup
# ---------------------------------------------------------------------------


class TestRemapOperator:
    """Tests for _remap_operator."""

    def test_known_remap(self) -> None:
        """Known operator is remapped."""
        remap = {"laplacian_z": "laplacian_x"}
        assert _remap_operator("laplacian_z", remap) == "laplacian_x"

    def test_known_none(self) -> None:
        """Killed operator maps to None."""
        remap: dict[str, str | None] = {"laplacian_x": None}
        assert _remap_operator("laplacian_x", remap) is None

    def test_unknown_axis_op(self) -> None:
        """Unknown axis-specific operator defaults to None."""
        assert _remap_operator("gradient_w", {}) is None

    def test_unknown_non_axis_preserved(self) -> None:
        """Unknown non-axis operator is preserved."""
        assert _remap_operator("my_custom_op", {}) == "my_custom_op"


# ---------------------------------------------------------------------------
# Unit tests: term transformation
# ---------------------------------------------------------------------------


class TestTransformTerm:
    """Tests for _transform_term."""

    def test_killed_operator_returns_none(self) -> None:
        """Term with killed operator returns None."""
        remap: dict[str, str | None] = {"laplacian_x": None}
        term: dict[str, Any] = {
            "coefficient": 1.0,
            "operator": "laplacian_x",
            "field": "phi_0",
        }
        assert _transform_term(term, remap, "z", {"x", "y"}) is None

    def test_surviving_operator_remapped(self) -> None:
        """Surviving operator is remapped correctly."""
        remap: dict[str, str | None] = {"laplacian_z": "laplacian_x"}
        term: dict[str, Any] = {
            "coefficient": 1.0,
            "operator": "laplacian_z",
            "field": "phi_0",
        }
        result = _transform_term(term, remap, "z", {"x", "y"})
        assert result is not None
        assert result["operator"] == "laplacian_x"

    def test_coefficient_symbolic_remapped(self) -> None:
        """coefficient_symbolic has coordinate references remapped."""
        remap: dict[str, str | None] = {"gradient_z": "gradient_x"}
        term: dict[str, Any] = {
            "coefficient": 1.0,
            "operator": "gradient_z",
            "field": "phi_0",
            "coefficient_symbolic": "2/z[]",
            "coordinate_dependent": ["z"],
        }
        result = _transform_term(term, remap, "z", {"x", "y"})
        assert result is not None
        assert result["coefficient_symbolic"] == "2/x[]"
        assert result["coordinate_dependent"] == ["x"]

    def test_killed_coord_in_surviving_term_raises(self) -> None:
        """Surviving term with killed coordinate_dependent entry raises."""
        remap: dict[str, str | None] = {"laplacian_x": "laplacian_x"}
        term: dict[str, Any] = {
            "coefficient": 1.0,
            "operator": "laplacian_x",
            "field": "phi_0",
            "coordinate_dependent": ["x", "y"],
        }
        with pytest.raises(ValueError, match="killed coordinate"):
            _transform_term(term, remap, "x", {"y"})


# ---------------------------------------------------------------------------
# Unit tests: hamiltonian term transformation
# ---------------------------------------------------------------------------


class TestTransformHamiltonianTerm:
    """Tests for _transform_hamiltonian_term."""

    def test_killed_factor_returns_none(self) -> None:
        """If either factor has a killed operator, return None."""
        remap: dict[str, str | None] = {"gradient_y": None, "gradient_z": "gradient_x"}
        term: dict[str, Any] = {
            "coefficient": 0.5,
            "factor_a": {"field": "phi_0", "operator": "gradient_y"},
            "factor_b": {"field": "phi_0", "operator": "gradient_y"},
        }
        assert _transform_hamiltonian_term(term, remap, "z", {"x", "y"}) is None

    def test_surviving_factors_remapped(self) -> None:
        """Surviving factor operators are remapped."""
        remap: dict[str, str | None] = {"gradient_z": "gradient_x"}
        term: dict[str, Any] = {
            "coefficient": 0.5,
            "factor_a": {"field": "phi_0", "operator": "gradient_z"},
            "factor_b": {"field": "phi_0", "operator": "gradient_z"},
        }
        result = _transform_hamiltonian_term(term, remap, "z", {"x", "y"})
        assert result is not None
        assert result["factor_a"]["operator"] == "gradient_x"
        assert result["factor_b"]["operator"] == "gradient_x"


# ---------------------------------------------------------------------------
# Integration tests: reduce_spec on synthetic specs
# ---------------------------------------------------------------------------


class TestReduceSpec:
    """Tests for reduce_spec on synthetic specs."""

    def test_flat_3d_to_1d_z(self) -> None:
        """Flat 3+1D KG reduced along z produces correct 1+1D spec."""
        spec = _make_3d_kg_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        assert result["spacetime"]["dimension"] == 2
        assert result["spacetime"]["signature"] == [-1, 1]
        assert result["spacetime"]["coordinates"] == ["t", "x"]

    def test_flat_3d_to_1d_x(self) -> None:
        """Flat 3+1D KG reduced along x also works."""
        spec = _make_3d_kg_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "x"})

        assert result["spacetime"]["dimension"] == 2
        assert result["spacetime"]["coordinates"] == ["t", "x"]

    def test_surviving_equation_terms(self) -> None:
        """Only identity and laplacian_z survive when reducing along z."""
        spec = _make_3d_kg_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        terms = result["equations"][0]["rhs"]["terms"]
        operators = [t["operator"] for t in terms]
        assert "identity" in operators
        assert "laplacian_x" in operators
        assert len(operators) == 2
        assert "laplacian_y" not in operators
        assert "laplacian_z" not in operators

    def test_hamiltonian_terms_filtered(self) -> None:
        """Killed gradient operators removed from hamiltonian_terms."""
        spec = _make_3d_kg_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        ham_terms = result["canonical"]["hamiltonian_terms"]
        # Should keep: identity*identity, gradient_x (from z), time_derivative
        ops = [
            (h["factor_a"]["operator"], h["factor_b"]["operator"])
            for h in ham_terms
        ]
        assert ("identity", "identity") in ops
        assert ("gradient_x", "gradient_x") in ops
        assert ("time_derivative", "time_derivative") in ops
        # Should remove: gradient_x*gradient_x (original x), gradient_y*gradient_y
        assert len(ham_terms) == 3

    def test_field_elimination(self) -> None:
        """Field with all-killed operators is eliminated."""
        spec = _make_2field_3d_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        field_names = [f["name"] for f in result["fields"]]
        assert "phi_0" in field_names
        assert "psi_0" not in field_names

    def test_eliminated_fields_in_provenance(self) -> None:
        """Provenance records eliminated fields."""
        spec = _make_2field_3d_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        reduction = result["metadata"]["reduction"]
        assert reduction["type"] == "plane_wave"
        assert reduction["original_dimension"] == 4
        assert reduction["propagation_axis"] == "z"
        assert "psi_0" in reduction["eliminated_fields"]

    def test_field_reindexing(self) -> None:
        """Surviving fields are reindexed from 0."""
        spec = _make_2field_3d_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        for i, field in enumerate(result["fields"]):
            assert field["index"] == i

    def test_curved_surviving_coord_remapped(self) -> None:
        """Curved spec with propagation-axis-only coefficients succeeds."""
        spec = _make_curved_spec_surviving()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        terms = result["equations"][0]["rhs"]["terms"]
        gradient_term = next(t for t in terms if t["operator"] == "gradient_x")
        assert gradient_term["coefficient_symbolic"] == "2/x[]"
        assert gradient_term["coordinate_dependent"] == ["x"]

    def test_curved_volume_element_remapped(self) -> None:
        """Volume element depending on propagation axis only is remapped."""
        spec = _make_curved_spec_surviving()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        assert result["canonical"]["volume_element"] == "x[]^2"

    def test_volume_element_trivial_removed(self) -> None:
        """Volume element of '1' is removed from output."""
        spec = _make_3d_kg_spec()
        spec["canonical"]["volume_element"] = "1"
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        assert "volume_element" not in result["canonical"]

    def test_killed_coord_in_coefficient_raises(self) -> None:
        """Incompatible reduction: surviving coeff references killed coord."""
        spec = _make_curved_spec_killed_ref()
        with pytest.raises(ValueError, match="killed coordinate"):
            reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "x"})

    def test_killed_coord_in_volume_element_raises(self) -> None:
        """Volume element referencing killed coord raises ValueError."""
        spec = _make_curved_spec_surviving()
        spec["canonical"]["volume_element"] = "z[]^2 * y[]"
        with pytest.raises(ValueError, match="killed coordinate"):
            reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

    def test_coupling_matrix_reduced(self) -> None:
        """Coupling matrices are filtered to surviving fields."""
        spec = _make_2field_3d_spec()
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        matrix = result["coupling"]["mass_matrix_symbolic"]
        assert len(matrix) == 1
        assert len(matrix[0]) == 1

    def test_deep_copy_input_unchanged(self) -> None:
        """Input spec_data is not mutated."""
        spec = _make_3d_kg_spec()
        original_dim = spec["spacetime"]["dimension"]
        reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})
        assert spec["spacetime"]["dimension"] == original_dim


# ---------------------------------------------------------------------------
# Integration tests: reduce_spec on real JSON files
# ---------------------------------------------------------------------------


class TestReduceRealSpecs:
    """Tests against real JSON specs from examples/data/."""

    @staticmethod
    def _load_spec(name: str) -> dict[str, Any]:
        path = _EXAMPLES_DATA / name
        if not path.exists():
            pytest.skip(f"{name} not found")
        with path.open() as f:
            return json.load(f)

    def test_klein_gordon_3d_along_z(self) -> None:
        """KG 3+1D reduced along z → 1+1D with 1 field, wave equation."""
        spec = self._load_spec("klein_gordon_3d.json")
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        assert result["spacetime"]["dimension"] == 2
        assert result["spacetime"]["coordinates"] == ["t", "x"]
        assert len(result["fields"]) == 1
        assert len(result["equations"]) == 1

        operators = [t["operator"] for t in result["equations"][0]["rhs"]["terms"]]
        assert "identity" in operators
        assert "laplacian_x" in operators
        assert "laplacian_y" not in operators
        assert "laplacian_z" not in operators

    def test_klein_gordon_3d_along_x(self) -> None:
        """KG 3+1D reduced along x → same result (isotropic)."""
        spec = self._load_spec("klein_gordon_3d.json")
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "x"})

        assert result["spacetime"]["dimension"] == 2
        assert len(result["equations"]) == 1
        operators = [t["operator"] for t in result["equations"][0]["rhs"]["terms"]]
        assert "laplacian_x" in operators

    def test_cylindrical_kg_along_z(self) -> None:
        """Cylindrical KG reduced along z → flat 1D.

        The existing cylindrical_kg.json has volume_element="Abs[x[]]" which
        references the killed coordinate x(=rho). In a real derive, the Wolfram
        pipeline would produce a factored volume element of 1. Here we strip
        the volume element to simulate that, testing equation transformation.
        """
        spec = self._load_spec("cylindrical_kg.json")
        # Strip vol element that would be factored by Wolfram pipeline
        spec.get("canonical", {}).pop("volume_element", None)
        result = reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

        assert result["spacetime"]["dimension"] == 2
        assert result["spacetime"]["coordinates"] == ["t", "x"]
        assert len(result["equations"]) == 1

        # Cylindrical z-axis is flat: only identity and laplacian_x survive
        terms = result["equations"][0]["rhs"]["terms"]
        operators = [t["operator"] for t in terms]
        assert "identity" in operators
        assert "laplacian_x" in operators

    def test_cylindrical_kg_unfactored_vol_errors(self) -> None:
        """Existing cylindrical_kg.json vol element references killed coords."""
        spec = self._load_spec("cylindrical_kg.json")
        with pytest.raises(ValueError, match="killed coordinate"):
            reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "z"})

    def test_sphere_kg_incompatible(self) -> None:
        """Sphere KG (stereographic) is incompatible — non-separable metric."""
        spec = self._load_spec("sphere_kg.json")
        with pytest.raises(ValueError, match="killed coordinate"):
            reduce_spec(spec, {"type": "plane_wave", "propagation_axis": "x"})


# ---------------------------------------------------------------------------
# WLS generation tests (dry-run)
# ---------------------------------------------------------------------------


class TestWlsReduction:
    """Tests for plane-wave reduction WLS generation."""

    @staticmethod
    def _generate(config: dict[str, object]) -> str:
        """Build TOML config dict, convert to WLS string."""
        from tidal.cli._derive import generate_wls

        base: dict[str, object] = {
            "theory": {"name": "Test"},
            "spacetime": {"dimension": 4, "metric": "minkowski"},
            "fields": [{"name": "phi", "type": "scalar"}],
            "lagrangian": {"expression": "-1/2 (d phi)^2 - m2/2 phi^2"},
            "output": {"path": "out.json"},
        }
        base.update(config)
        return generate_wls(base)

    def test_reduction_before_el(self) -> None:
        """Reduction code appears before E-L derivation when [reduction] present."""
        wls = self._generate({
            "reduction": {"type": "plane_wave", "propagation_axis": "z"},
        })
        reduction_idx = wls.find("Plane-wave reduction")
        el_idx = wls.find("Euler-Lagrange equations")
        assert reduction_idx != -1, "Plane-wave reduction block not found"
        assert el_idx != -1, "Euler-Lagrange equations block not found"
        assert reduction_idx < el_idx

    def test_no_reduction_no_code(self) -> None:
        """Without [reduction], no reduction code in WLS."""
        wls = self._generate({})
        assert "Plane-wave reduction" not in wls

    def test_field_elimination_code_present(self) -> None:
        """Field elimination code is generated when reduction is active."""
        wls = self._generate({
            "reduction": {"type": "plane_wave", "propagation_axis": "z"},
        })
        assert "zeroFieldNames" in wls
        assert "Eliminating zero fields" in wls

    def test_factored_volume_element_for_reduction(self) -> None:
        """Factored volume element code generated when reduction is active."""
        wls = self._generate({
            "reduction": {"type": "plane_wave", "propagation_axis": "z"},
            "spacetime": {
                "dimension": 4,
                "metric": "diagonal",
                "diagonal": [-1, "x[]^2", "x[]^2", 1],
            },
        })
        assert "FreeQ" in wls or "killedVars" in wls

    def test_derivative_slots_correct_for_z(self) -> None:
        """z-propagation in 3+1D kills slots 2 (x) and 3 (y)."""
        wls = self._generate({
            "reduction": {"type": "plane_wave", "propagation_axis": "z"},
        })
        # Slots 2 and 3 should be zeroed (x=slot 2, y=slot 3 in t,x,y,z)
        assert "{ords}[[2]]" in wls
        assert "{ords}[[3]]" in wls

    def test_derivative_slots_correct_for_x(self) -> None:
        """x-propagation in 3+1D kills slots 3 (y) and 4 (z)."""
        wls = self._generate({
            "reduction": {"type": "plane_wave", "propagation_axis": "x"},
        })
        assert "{ords}[[3]]" in wls
        assert "{ords}[[4]]" in wls


# ---------------------------------------------------------------------------
# TOML validation tests
# ---------------------------------------------------------------------------


class TestValidateReduction:
    """Tests for [reduction] TOML validation in _validate_config."""

    @staticmethod
    def _validate(config: dict[str, object]) -> None:
        from tidal.cli._derive import _validate_reduction

        _validate_reduction(config)

    def test_valid_plane_wave(self) -> None:
        """Valid plane_wave reduction config accepted."""
        self._validate({
            "spacetime": {"dimension": 4},
            "reduction": {"type": "plane_wave", "propagation_axis": "z"},
        })

    def test_no_reduction_ok(self) -> None:
        """No [reduction] section is fine."""
        self._validate({"spacetime": {"dimension": 4}})

    def test_invalid_type(self) -> None:
        """Non-plane_wave type rejected."""
        with pytest.raises(ValueError, match="plane_wave"):
            self._validate({
                "spacetime": {"dimension": 4},
                "reduction": {"type": "spherical"},
            })

    def test_missing_propagation_axis(self) -> None:
        """Missing propagation_axis rejected."""
        with pytest.raises(ValueError, match="propagation_axis"):
            self._validate({
                "spacetime": {"dimension": 4},
                "reduction": {"type": "plane_wave"},
            })

    def test_invalid_propagation_axis(self) -> None:
        """Invalid propagation_axis rejected."""
        with pytest.raises(ValueError, match="valid spatial"):
            self._validate({
                "spacetime": {"dimension": 4},
                "reduction": {"type": "plane_wave", "propagation_axis": "w"},
            })

    def test_1d_cannot_reduce(self) -> None:
        """Cannot reduce a 1+1D theory."""
        with pytest.raises(ValueError, match="cannot reduce"):
            self._validate({
                "spacetime": {"dimension": 2},
                "reduction": {"type": "plane_wave", "propagation_axis": "x"},
            })

    def test_2d_valid(self) -> None:
        """2+1D can be reduced along x."""
        self._validate({
            "spacetime": {"dimension": 3},
            "reduction": {"type": "plane_wave", "propagation_axis": "x"},
        })


# ---------------------------------------------------------------------------
# Tests for degenerate algebraic constraint elimination
# ---------------------------------------------------------------------------


def _make_traceless_constraint_spec() -> dict[str, Any]:
    """Minimal spec with a degenerate traceless constraint.

    Mimics TT gauge: h_9 = h_4 + h_7 + h_9 → 0 = h_4 + h_7 → h_7 = -h_4.
    h_4 and h_7 have cross-Laplacians (the unstable pattern).
    """
    return {
        "metadata": {},
        "spacetime": {"dimension": 2, "signature": [-1, 1], "coordinates": ["t", "x"]},
        "fields": [
            {"name": "h_4", "index": 0, "is_dynamical": True},
            {"name": "h_7", "index": 1, "is_dynamical": True},
            {"name": "h_9", "index": 2, "is_dynamical": True},
        ],
        "equations": [
            {
                "field": "h_4",
                "lhs": {"expression": "d2_t(h_4)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": -1.0, "operator": "laplacian_x", "field": "h_7"},
                        {"coefficient": 0.3, "operator": "gradient_x", "field": "a_2"},
                    ],
                },
            },
            {
                "field": "h_7",
                "lhs": {"expression": "d2_t(h_7)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 0.5, "operator": "identity", "field": "h_4"},
                        {"coefficient": -0.5, "operator": "identity", "field": "h_7"},
                        {"coefficient": -1.0, "operator": "laplacian_x", "field": "h_4"},
                    ],
                },
            },
            {
                "field": "h_9",
                "lhs": {"expression": "h_9", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "identity", "field": "h_4"},
                        {"coefficient": 1.0, "operator": "identity", "field": "h_7"},
                        {"coefficient": 1.0, "operator": "identity", "field": "h_9"},
                    ],
                },
            },
            {
                "field": "a_2",
                "lhs": {"expression": "d2_t(a_2)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": -0.3, "operator": "gradient_x", "field": "h_4"},
                        {"coefficient": 1.0, "operator": "laplacian_x", "field": "a_2"},
                    ],
                },
            },
        ],
        "coupling": {
            "mass_matrix_symbolic": [
                [None, None, None],
                [None, "-0.5", None],
                [None, None, None],
            ],
        },
        "canonical": {"hamiltonian_terms": []},
    }


class TestFindDegenerateConstraint:
    """Tests for _find_degenerate_constraint."""

    def test_finds_traceless(self) -> None:
        """Detects h_9 = h_4 + h_7 + h_9 as degenerate."""
        spec = _make_traceless_constraint_spec()
        result = _find_degenerate_constraint(spec["equations"])
        assert result is not None
        constraint_eq_field, dep_field, substitution = result
        assert constraint_eq_field == "h_9"
        # Should eliminate one of h_4/h_7
        assert dep_field in ("h_4", "h_7")
        # Substitution: dep = -1 * other
        assert len(substitution) == 1
        other = next(iter(substitution))
        assert other in ("h_4", "h_7")
        assert other != dep_field
        assert substitution[other] == pytest.approx(-1.0)

    def test_no_degenerate_in_normal_spec(self) -> None:
        """Normal algebraic equations (self-term != 1) are not flagged."""
        equations = [
            {
                "field": "psi",
                "lhs": {"expression": "psi", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 0.5, "operator": "identity", "field": "psi"},
                        {"coefficient": 1.0, "operator": "identity", "field": "phi"},
                    ],
                },
            },
        ]
        assert _find_degenerate_constraint(equations) is None

    def test_ignores_dynamical_equations(self) -> None:
        """Dynamical (time_order=2) equations are skipped."""
        equations = [
            {
                "field": "phi",
                "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "identity", "field": "phi"},
                    ],
                },
            },
        ]
        assert _find_degenerate_constraint(equations) is None

    def test_ignores_non_identity(self) -> None:
        """Constraints with derivative operators are not degenerate (for now)."""
        equations = [
            {
                "field": "psi",
                "lhs": {"expression": "psi", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "gradient_x", "field": "phi"},
                        {"coefficient": 1.0, "operator": "identity", "field": "psi"},
                    ],
                },
            },
        ]
        assert _find_degenerate_constraint(equations) is None


class TestSubstituteInTerms:
    """Tests for _substitute_in_terms."""

    def test_simple_substitution(self) -> None:
        """h_7 → -h_4 in a Laplacian term."""
        terms = [
            {"coefficient": -1.0, "operator": "laplacian_x", "field": "h_7"},
            {"coefficient": 0.3, "operator": "gradient_x", "field": "a_2"},
        ]
        result = _substitute_in_terms(terms, "h_7", {"h_4": -1.0})
        assert len(result) == 2
        # -1.0 * laplacian_x(h_7) → -1.0 * (-1.0) * laplacian_x(h_4) = +1.0
        lap = next(t for t in result if t["operator"] == "laplacian_x")
        assert lap["field"] == "h_4"
        assert lap["coefficient"] == pytest.approx(1.0)

    def test_merge_terms(self) -> None:
        """Terms on same (op, field) are merged."""
        terms = [
            {"coefficient": 2.0, "operator": "identity", "field": "h_4"},
            {"coefficient": 3.0, "operator": "identity", "field": "h_7"},
        ]
        # h_7 → -h_4: 3.0 * identity(h_7) → 3.0 * (-1) * identity(h_4) = -3.0
        # Merged with existing 2.0 * identity(h_4) → -1.0
        result = _substitute_in_terms(terms, "h_7", {"h_4": -1.0})
        assert len(result) == 1
        assert result[0]["field"] == "h_4"
        assert result[0]["coefficient"] == pytest.approx(-1.0)

    def test_cancel_to_zero(self) -> None:
        """Terms that cancel to zero are removed."""
        terms = [
            {"coefficient": 1.0, "operator": "gradient_x", "field": "h_4"},
            {"coefficient": 1.0, "operator": "gradient_x", "field": "h_7"},
        ]
        # h_7 → -h_4: grad(h_4) + grad(h_7) → grad(h_4) - grad(h_4) = 0
        result = _substitute_in_terms(terms, "h_7", {"h_4": -1.0})
        assert len(result) == 0

    def test_zero_substitution(self) -> None:
        """dep_field = 0 removes all terms referencing it."""
        terms = [
            {"coefficient": 1.0, "operator": "identity", "field": "h_7"},
            {"coefficient": 2.0, "operator": "laplacian_x", "field": "phi"},
        ]
        result = _substitute_in_terms(terms, "h_7", {})
        assert len(result) == 1
        assert result[0]["field"] == "phi"


class TestEliminateDegenerateConstraints:
    """Tests for the full elimination pipeline."""

    def test_traceless_elimination(self) -> None:
        """TT gauge traceless constraint is properly eliminated."""
        spec = _make_traceless_constraint_spec()
        result = eliminate_degenerate_constraints(spec)

        # Should have eliminated h_7 (or h_4) and h_9
        elim = result["metadata"]["constraint_elimination"]["eliminated_fields"]
        assert len(elim) >= 2  # dep_field + constraint_eq_field

        field_names = {f["name"] for f in result["fields"]}
        assert "h_9" not in field_names
        # One of h_4/h_7 survives, the other is eliminated
        assert "h_4" in field_names or "h_7" in field_names

        eq_fields = {eq["field"] for eq in result["equations"]}
        assert "h_9" not in eq_fields

    def test_cross_laplacian_resolved(self) -> None:
        """After elimination, h_4 equation gets self-Laplacian (no cross)."""
        spec = _make_traceless_constraint_spec()
        result = eliminate_degenerate_constraints(spec)

        # Find the surviving graviton equation
        surviving = next(
            f["name"] for f in result["fields"]
            if f["name"] in ("h_4", "h_7")
        )
        eq = next(eq for eq in result["equations"] if eq["field"] == surviving)

        # Should have laplacian_x on itself (self-Laplacian), not cross
        lap_terms = [
            t for t in eq["rhs"]["terms"]
            if t["operator"] == "laplacian_x"
        ]
        assert len(lap_terms) == 1
        assert lap_terms[0]["field"] == surviving
        assert lap_terms[0]["coefficient"] == pytest.approx(1.0)

    def test_trivial_equations_removed(self) -> None:
        """Equations that become 0=0 after substitution are removed."""
        spec = _make_traceless_constraint_spec()
        # Add a subsidiary constraint that references only h_4 and h_7
        spec["fields"].append(
            {"name": "h_tz", "index": 3, "is_dynamical": True},
        )
        spec["equations"].append({
            "field": "h_tz",
            "lhs": {"expression": "h_tz", "order": {"time": 0, "space": 0}},
            "rhs": {
                "type": "linear_combination",
                "terms": [
                    {"coefficient": -1.0, "operator": "gradient_x", "field": "h_4"},
                    {"coefficient": -1.0, "operator": "gradient_x", "field": "h_7"},
                ],
            },
        })
        result = eliminate_degenerate_constraints(spec)

        # h_tz equation should be gone (gradient(h_4) + gradient(h_7) = 0)
        eq_fields = {eq["field"] for eq in result["equations"]}
        assert "h_tz" not in eq_fields

    def test_no_degenerate_passthrough(self) -> None:
        """Spec without degenerate constraints is returned unchanged."""
        spec = {
            "metadata": {},
            "fields": [{"name": "phi", "index": 0, "is_dynamical": True}],
            "equations": [
                {
                    "field": "phi",
                    "lhs": {"expression": "d2_t(phi)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "phi"},
                        ],
                    },
                },
            ],
            "coupling": {},
            "canonical": {"hamiltonian_terms": []},
        }
        result = eliminate_degenerate_constraints(spec)
        assert result is spec  # Same object (no copy needed)

    def test_field_reindexing(self) -> None:
        """Surviving fields are reindexed 0, 1, 2, ..."""
        spec = _make_traceless_constraint_spec()
        result = eliminate_degenerate_constraints(spec)
        for idx, field in enumerate(result["fields"]):
            assert field["index"] == idx

    def test_symbolic_coefficients_preserved(self) -> None:
        """Symbolic coefficients survive substitution and merging."""
        spec = {
            "metadata": {},
            "fields": [
                {"name": "f1", "index": 0, "is_dynamical": True},
                {"name": "f2", "index": 1, "is_dynamical": True},
                {"name": "trace", "index": 2, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "f1",
                    "lhs": {"expression": "d2_t(f1)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 0.5,
                                "operator": "identity",
                                "field": "f2",
                                "coefficient_symbolic": "alpha/2",
                            },
                            {
                                "coefficient": 1.0,
                                "operator": "laplacian_x",
                                "field": "f1",
                            },
                        ],
                    },
                },
                {
                    "field": "f2",
                    "lhs": {"expression": "d2_t(f2)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {
                                "coefficient": 0.5,
                                "operator": "identity",
                                "field": "f1",
                                "coefficient_symbolic": "alpha/2",
                            },
                            {
                                "coefficient": -0.5,
                                "operator": "identity",
                                "field": "f2",
                                "coefficient_symbolic": "-alpha/2",
                            },
                            {
                                "coefficient": -1.0,
                                "operator": "laplacian_x",
                                "field": "f1",
                            },
                        ],
                    },
                },
                {
                    "field": "trace",
                    "lhs": {"expression": "trace", "order": {"time": 0, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "identity", "field": "f1"},
                            {"coefficient": 1.0, "operator": "identity", "field": "f2"},
                            {"coefficient": 1.0, "operator": "identity", "field": "trace"},
                        ],
                    },
                },
            ],
            "coupling": {},
            "canonical": {"hamiltonian_terms": []},
        }

        result = eliminate_degenerate_constraints(spec)

        # After elimination: f1 → -f2, so f2's equation should have merged
        # mass terms with symbolic: "-(alpha/2) + (-alpha/2)" = "-alpha"
        surviving = next(
            f["name"] for f in result["fields"]
            if f["name"] in ("f1", "f2")
        )
        eq = next(eq for eq in result["equations"] if eq["field"] == surviving)
        mass_term = next(
            t for t in eq["rhs"]["terms"]
            if t["operator"] == "identity"
        )
        # Must have symbolic coefficient
        assert "coefficient_symbolic" in mass_term, (
            f"Symbolic coefficient lost during elimination for "
            f"{mass_term['coefficient']}*identity({mass_term['field']})"
        )
        assert mass_term["coefficient"] == pytest.approx(-1.0)

    def test_on_real_gertsenshtein_json(self) -> None:
        """Integration test: constraint elimination on actual Gertsenshtein JSON."""
        json_path = _EXAMPLES_DATA / "gertsenshtein.json"
        if not json_path.exists():
            pytest.skip("gertsenshtein.json not found")

        spec = json.loads(json_path.read_text(encoding="utf-8"))
        result = eliminate_degenerate_constraints(spec)

        elim = result["metadata"]["constraint_elimination"]["eliminated_fields"]
        assert len(elim) >= 2  # h_7 + h_9 at minimum

        field_names = {f["name"] for f in result["fields"]}
        assert "h_9" not in field_names

        # The surviving system should have self-Laplacians, no cross-Laplacians
        for eq in result["equations"]:
            for term in eq["rhs"]["terms"]:
                if term["operator"] == "laplacian_x":
                    assert term["field"] == eq["field"], (
                        f"Cross-Laplacian found: {eq['field']} has "
                        f"laplacian_x({term['field']})"
                    )


# ---------------------------------------------------------------------------
# Gradient-zero constraint elimination
# ---------------------------------------------------------------------------


class TestFindGradientZeroConstraint:
    """Tests for _find_gradient_zero_constraint()."""

    def test_detects_gradient_self_constraint(self) -> None:
        """Single gradient_z self-term → field constrained to zero."""
        equations = [
            {
                "field": "h_6",
                "lhs": {"expression": "h_6", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) == "h_6"

    def test_detects_laplacian_self_constraint(self) -> None:
        """Single laplacian self-term → field constrained to zero."""
        equations = [
            {
                "field": "phi",
                "lhs": {"expression": "phi", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "laplacian", "field": "phi"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) == "phi"

    def test_ignores_identity_constraint(self) -> None:
        """Identity self-term is NOT a gradient-zero constraint."""
        equations = [
            {
                "field": "h_9",
                "lhs": {"expression": "h_9", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "identity", "field": "h_9"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) is None

    def test_ignores_dynamical_equations(self) -> None:
        """time_order > 0 equations are not constraints."""
        equations = [
            {
                "field": "h_6",
                "lhs": {"expression": "d2_t(h_6)", "order": {"time": 2, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) is None

    def test_ignores_multi_term_constraints(self) -> None:
        """Multi-field divergence constraint is NOT a single-term gradient-zero."""
        equations = [
            {
                "field": "h_6",
                "lhs": {"expression": "h_6", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "gradient_x", "field": "h_4"},
                        {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) is None

    def test_ignores_cross_field_gradient(self) -> None:
        """Gradient of a DIFFERENT field is not a self-constraint."""
        equations = [
            {
                "field": "h_t",
                "lhs": {"expression": "h_t", "order": {"time": 0, "space": 0}},
                "rhs": {
                    "type": "linear_combination",
                    "terms": [
                        {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                    ],
                },
            },
        ]
        assert _find_gradient_zero_constraint(equations) is None


class TestGradientZeroElimination:
    """Integration tests for gradient-zero elimination in eliminate_degenerate_constraints."""

    def test_gradient_zero_elimination(self) -> None:
        """Field with gradient_z(f)=0 constraint is eliminated and substituted with zero."""
        spec = {
            "metadata": {},
            "fields": [
                {"name": "h_5", "index": 0, "is_dynamical": True},
                {"name": "h_6", "index": 1, "is_dynamical": True},
                {"name": "h_7", "index": 2, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "h_5",
                    "lhs": {"expression": "d2_t(h_5)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_5"},
                        ],
                    },
                },
                {
                    "field": "h_6",
                    "lhs": {"expression": "h_6", "order": {"time": 0, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                        ],
                    },
                },
                {
                    "field": "h_7",
                    "lhs": {"expression": "d2_t(h_7)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_7"},
                        ],
                    },
                },
            ],
            "coupling": {},
            "canonical": {"hamiltonian_terms": []},
        }

        result = eliminate_degenerate_constraints(spec)

        field_names = {f["name"] for f in result["fields"]}
        assert field_names == {"h_5", "h_7"}, f"Expected h_5, h_7 but got {field_names}"
        assert "h_6" in result["metadata"]["constraint_elimination"]["eliminated_fields"]

    def test_gradient_zero_removes_decoupled_terms(self) -> None:
        """Zero-field references are removed from other equations."""
        spec = {
            "metadata": {},
            "fields": [
                {"name": "h_5", "index": 0, "is_dynamical": True},
                {"name": "h_6", "index": 1, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "h_5",
                    "lhs": {"expression": "d2_t(h_5)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_5"},
                            {"coefficient": 0.5, "operator": "identity", "field": "h_6"},
                        ],
                    },
                },
                {
                    "field": "h_6",
                    "lhs": {"expression": "h_6", "order": {"time": 0, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                        ],
                    },
                },
            ],
            "coupling": {},
            "canonical": {"hamiltonian_terms": []},
        }

        result = eliminate_degenerate_constraints(spec)

        # h_6 eliminated, h_5 survives without the h_6 coupling term
        field_names = {f["name"] for f in result["fields"]}
        assert field_names == {"h_5"}
        eq = result["equations"][0]
        assert len(eq["rhs"]["terms"]) == 1
        assert eq["rhs"]["terms"][0]["field"] == "h_5"

    def test_cascade_degenerate_then_gradient_zero(self) -> None:
        """Degenerate elimination first, then gradient-zero cleans up remaining.

        Mimics GW 1D structure: traceless constraint h_9 = h_4 + h_7 + h_9
        (3-field, eliminates h_4 as h_4 = -h_7), then gradient-zero eliminates h_6.
        """
        spec = {
            "metadata": {},
            "fields": [
                {"name": "h_4", "index": 0, "is_dynamical": True},
                {"name": "h_5", "index": 1, "is_dynamical": True},
                {"name": "h_6", "index": 2, "is_dynamical": True},
                {"name": "h_7", "index": 3, "is_dynamical": True},
                {"name": "h_9", "index": 4, "is_dynamical": True},
            ],
            "equations": [
                {
                    "field": "h_4",
                    "lhs": {"expression": "d2_t(h_4)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_4"},
                        ],
                    },
                },
                {
                    "field": "h_5",
                    "lhs": {"expression": "d2_t(h_5)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_5"},
                        ],
                    },
                },
                {
                    "field": "h_6",
                    "lhs": {"expression": "h_6", "order": {"time": 0, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "gradient_z", "field": "h_6"},
                        ],
                    },
                },
                {
                    "field": "h_7",
                    "lhs": {"expression": "d2_t(h_7)", "order": {"time": 2, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "laplacian_x", "field": "h_7"},
                        ],
                    },
                },
                {
                    "field": "h_9",
                    "lhs": {"expression": "h_9", "order": {"time": 0, "space": 0}},
                    "rhs": {
                        "type": "linear_combination",
                        "terms": [
                            {"coefficient": 1.0, "operator": "identity", "field": "h_4"},
                            {"coefficient": 1.0, "operator": "identity", "field": "h_7"},
                            {"coefficient": 1.0, "operator": "identity", "field": "h_9"},
                        ],
                    },
                },
            ],
            "coupling": {},
            "canonical": {"hamiltonian_terms": []},
        }

        result = eliminate_degenerate_constraints(spec)

        # h_4 eliminated by degenerate (h_4 = -h_7 from traceless)
        # h_6 eliminated by gradient-zero
        # h_9 eliminated (constraint equation field)
        field_names = {f["name"] for f in result["fields"]}
        assert field_names == {"h_5", "h_7"}, f"Expected h_5, h_7 but got {field_names}"
        elim = result["metadata"]["constraint_elimination"]["eliminated_fields"]
        assert "h_4" in elim
        assert "h_6" in elim
