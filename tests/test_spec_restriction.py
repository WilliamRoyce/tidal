"""GH #468 route 3: dependency closure + observable-sector spec restriction.

``EquationSystem.dependency_closure`` is the ONE definition of "which
fields the evolution of a seed set can ever read"; ``restrict_spec_dict``
turns a closure into a first-class spec JSON (fields outside the closure
ABSENT, never present-but-wrong; Hamiltonian terms touching them dropped
and counted). These tests are the regression contract:

* unit tests pin the closure semantics on a hand-built chain spec;
* restriction round-trips through ``EquationSystem.from_dict`` (the
  loader's own validation is the consistency check), re-indexes fields,
  filters Hamiltonian terms, records provenance, and refuses non-closed
  sets;
* the corpus pin fixes the E.cal fact the whole route rests on:
  closure({h_5, a_1}) == {h_5, a_1} while closure({a_2}) pulls the
  implicit-dynamical sector (the documented counter-case).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from tidal.symbolic.json_loader import EquationSystem, restrict_spec_dict

REPO_ROOT = Path(__file__).resolve().parent.parent
ECAL = REPO_ROOT / "examples/data/gertsenshtein_ungauged_e_dual_gaussian.json"


def _eq(name: str, order: int, terms: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "field": name,
        "lhs": {"expression": name, "order": {"time": order, "space": 0}},
        "rhs": {"type": "linear_combination", "terms": terms},
    }


def _term(coeff: float, op: str, field: str) -> dict[str, Any]:
    return {"coefficient": coeff, "operator": op, "field": field}


def _chain_spec_dict() -> dict[str, Any]:
    """A reads v_b, b reads c, c is closed, d reads a (d is sourced BY a)."""
    return {
        "spacetime": {"dimension": 2, "coordinates": ["t", "x"]},
        "fields": [
            {"name": "a", "index": 0},
            {"name": "b", "index": 1},
            {"name": "c", "index": 2},
            {"name": "d", "index": 3},
        ],
        "equations": [
            _eq("a", 2, [_term(-1.0, "laplacian", "a"), _term(0.5, "identity", "v_b")]),
            _eq("b", 2, [_term(-1.0, "laplacian", "b"), _term(0.1, "identity", "c")]),
            _eq("c", 2, [_term(-1.0, "laplacian", "c")]),
            _eq("d", 2, [_term(-1.0, "laplacian", "d"), _term(0.3, "identity", "a")]),
        ],
        "canonical": {
            "hamiltonian_terms": [
                {
                    "coefficient": 0.5,
                    "factor_a": {"field": "v_a", "operator": "identity"},
                    "factor_b": {"field": "v_a", "operator": "identity"},
                },
                {
                    "coefficient": 0.3,
                    "factor_a": {"field": "d", "operator": "identity"},
                    "factor_b": {"field": "a", "operator": "identity"},
                },
            ]
        },
        "metadata": {"parameters": {}},
    }


class TestDependencyClosure:
    def test_chain_semantics(self) -> None:
        spec = EquationSystem.from_dict(_chain_spec_dict())
        assert spec.dependency_closure({"c"}) == {"c"}
        assert spec.dependency_closure({"b"}) == {"b", "c"}
        # velocity reference v_b counts as b
        assert spec.dependency_closure({"a"}) == {"a", "b", "c"}
        # d is sourced BY a's sector, but reads a → its closure is everything
        assert spec.dependency_closure({"d"}) == {"a", "b", "c", "d"}
        assert spec.dependency_closure(["c", "d"]) == {"a", "b", "c", "d"}

    def test_unknown_seed_raises(self) -> None:
        spec = EquationSystem.from_dict(_chain_spec_dict())
        with pytest.raises(ValueError, match="unknown field"):
            spec.dependency_closure({"zeta"})


class TestRestrictSpecDict:
    def test_round_trip_and_provenance(self) -> None:
        data = _chain_spec_dict()
        restricted, record = restrict_spec_dict(
            data, {"a", "b", "c"}, parent_spec="chain.json", seeds=["a"], reason="test"
        )
        spec = EquationSystem.from_dict(restricted)  # loader re-validation
        assert spec.component_names == ("a", "b", "c")
        assert [f["index"] for f in restricted["fields"]] == [0, 1, 2]
        assert [eq["field"] for eq in restricted["equations"]] == ["a", "b", "c"]
        # the d·a Hamiltonian term is dropped and counted; v_a·v_a survives
        assert len(restricted["canonical"]["hamiltonian_terms"]) == 1
        assert record.dropped_hamiltonian_terms == 1
        assert record.evolved == ("a", "b", "c")
        assert record.omitted == ("d",)
        assert record.seeds == ("a",)
        assert restricted["metadata"]["restriction"]["parent_spec"] == "chain.json"
        assert spec.metadata["restriction"]["omitted"] == ["d"]
        # the parent dict is untouched
        assert len(data["fields"]) == 4

    def test_not_closed_refuses(self) -> None:
        with pytest.raises(ValueError, match="not closed"):
            restrict_spec_dict(_chain_spec_dict(), {"a", "b"})  # b reads c

    def test_unknown_field_refuses(self) -> None:
        with pytest.raises(ValueError, match="unknown field"):
            restrict_spec_dict(_chain_spec_dict(), {"a", "zeta"})

    def test_full_set_is_identity_on_structure(self) -> None:
        data = _chain_spec_dict()
        restricted, record = restrict_spec_dict(data, {"a", "b", "c", "d"})
        assert record.omitted == ()
        assert record.dropped_hamiltonian_terms == 0
        assert EquationSystem.from_dict(restricted).component_names == (
            "a",
            "b",
            "c",
            "d",
        )


class TestECalCorpusPin:
    """The E.cal facts route 3 rests on (localized implicit-dynamical class)."""

    @pytest.fixture
    def ecal(self) -> tuple[dict[str, Any], EquationSystem]:
        if not ECAL.exists():
            pytest.skip("E.cal spec not present")
        data = json.loads(ECAL.read_text())
        return data, EquationSystem.from_dict(data)

    def test_observable_sector_is_exactly_closed(self, ecal) -> None:
        _data, spec = ecal
        assert spec.second_order_sector.promoted  # the class that refuses
        assert spec.dependency_closure({"h_5", "a_1"}) == {"h_5", "a_1"}

    def test_counter_case_pulls_the_implicit_dynamical_sector(self, ecal) -> None:
        _data, spec = ecal
        closure = spec.dependency_closure({"a_2"})
        assert closure & spec.second_order_sector.promoted

    def test_restricted_ecal_has_no_implicit_dynamical_sector(self, ecal) -> None:
        data, spec = ecal
        restricted, record = restrict_spec_dict(
            data, spec.dependency_closure({"h_5", "a_1"})
        )
        sub = EquationSystem.from_dict(restricted)
        assert sub.component_names == ("a_1", "h_5") or set(sub.component_names) == {
            "a_1",
            "h_5",
        }
        assert not sub.second_order_sector.promoted
        assert record.dropped_hamiltonian_terms > 0
        assert set(record.omitted) == set(spec.component_names) - {"h_5", "a_1"}
