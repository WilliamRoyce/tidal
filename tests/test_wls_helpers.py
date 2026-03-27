"""Tests for tidal.cli._wls_helpers — Wolfram Language code generation helpers."""

from __future__ import annotations

from tidal.cli._wls_helpers import (
    wl_bg_rule_entry,
    wl_component_value,
    wl_diag_matrix,
    wl_flatten_list,
    wl_index,
    wl_list,
    wl_rule,
    wl_rules_block,
    wl_skip_tuples,
    wl_zero_component,
)


class TestWlList:
    def test_empty(self) -> None:
        assert wl_list() == "{}"

    def test_single(self) -> None:
        assert wl_list("a") == "{a}"

    def test_multiple(self) -> None:
        assert wl_list("a", "b", "c") == "{a, b, c}"

    def test_nested(self) -> None:
        inner = wl_list("1", "2")
        assert wl_list(inner, "x") == "{{1, 2}, x}"

    def test_with_expressions(self) -> None:
        assert wl_list("x -> 0", "y :> 1") == "{x -> 0, y :> 1}"


class TestWlIndex:
    def test_covariant(self) -> None:
        assert wl_index(0, "geChart") == "{0, -geChart}"

    def test_contravariant(self) -> None:
        assert wl_index(2, "geChart", contra=True) == "{2, geChart}"

    def test_slot_3(self) -> None:
        assert wl_index(3, "ch") == "{3, -ch}"


class TestWlRule:
    def test_normal_rule(self) -> None:
        assert wl_rule("x", "0") == "x -> 0"

    def test_delayed_rule(self) -> None:
        assert wl_rule("f[x_]", "x^2", delayed=True) == "f[x_] :> x^2"


class TestWlRulesBlock:
    def test_single_rule(self) -> None:
        r = wl_rule("a", "0")
        assert wl_rules_block([r]) == "{a -> 0}"

    def test_multiple_rules(self) -> None:
        r1 = wl_rule("a[___]", "0", delayed=True)
        r2 = wl_rule("b[___]", "0", delayed=True)
        assert wl_rules_block([r1, r2]) == "{a[___] :> 0, b[___] :> 0}"


class TestWlComponentValue:
    def test_rank1(self) -> None:
        result = wl_component_value("geH", [(0, "ch")], "0")
        assert result == "ComponentValue[geH[{0, -ch}], 0];"

    def test_rank2(self) -> None:
        result = wl_component_value("geH", [(0, "ch"), (1, "ch")], "0")
        assert result == "ComponentValue[geH[{0, -ch}, {1, -ch}], 0];"


class TestWlSkipTuples:
    def test_dim2(self) -> None:
        result = wl_skip_tuples(2)
        assert result == '"SkipTuples" -> {{0, 0}, {0, 1}}'

    def test_dim4(self) -> None:
        result = wl_skip_tuples(4)
        assert result == '"SkipTuples" -> {{0, 0}, {0, 1}, {0, 2}, {0, 3}}'


class TestWlBgRuleEntry:
    def test_with_components(self) -> None:
        result = wl_bg_rule_entry("geAbar", "0, 0, -B0*z[], 0", "0, 0, B0*z[], 0")
        assert result == "{geAbar, {0, 0, -B0*z[], 0}, {0, 0, B0*z[], 0}}"

    def test_empty_components(self) -> None:
        result = wl_bg_rule_entry("geAbar", "", "")
        assert result == "{geAbar, {}, {}}"


class TestWlFlattenList:
    def test_basic(self) -> None:
        assert wl_flatten_list("a, b") == "Flatten[{a, b}]"


class TestWlDiagMatrix:
    def test_basic(self) -> None:
        assert wl_diag_matrix("-1, 1, 1, 1") == "DiagonalMatrix[{-1, 1, 1, 1}]"


class TestWlZeroComponent:
    def test_produces_lines(self) -> None:
        lines = wl_zero_component("geH0", "h", "temporal")
        assert len(lines) == 4
        assert "fieldEquations" in lines[1]
        assert "geH0[args___] :> 0" in lines[1]
        assert "Derivative[ders__][geH0][args___] :> 0" in lines[1]
