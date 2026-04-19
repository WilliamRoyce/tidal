"""Tests for tidal.cli._wls_helpers — Wolfram Language code generation helpers."""

from __future__ import annotations

from tidal.cli._wls_helpers import (
    validate_wls_brackets,
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
        # Default (apply_to_lagrangian=True): emits the fieldEquations rule
        # plus an If[ValueQ[lagComp], ...] block for Phase C / Phase D
        # consistency.
        lines = wl_zero_component("geH0", "h", "temporal")
        assert any("fieldEquations" in L for L in lines)
        assert any("geH0[args___] :> 0" in L for L in lines)
        assert any("Derivative[ders__][geH0][args___] :> 0" in L for L in lines)
        assert any("ValueQ[lagComp]" in L for L in lines)

    def test_apply_to_lagrangian_false(self) -> None:
        """When apply_to_lagrangian=False, no lagComp block is emitted."""
        lines = wl_zero_component("geH0", "h", "temporal", apply_to_lagrangian=False)
        assert not any("ValueQ[lagComp]" in L for L in lines)
        assert any("fieldEquations" in L for L in lines)


class TestValidateWlsBrackets:
    def test_balanced(self) -> None:
        assert validate_wls_brackets("f[x, {a, b}]") == []

    def test_nested_balanced(self) -> None:
        assert validate_wls_brackets("Module[{x = 0}, Do[Print[x], {x, 10}]]") == []

    def test_unmatched_close(self) -> None:
        errs = validate_wls_brackets("f[x]]")
        assert len(errs) == 1
        assert "unmatched" in errs[0]

    def test_unclosed_open(self) -> None:
        errs = validate_wls_brackets("f[x, {a, b}")
        assert len(errs) == 1
        assert "unclosed" in errs[0]

    def test_mismatched_pair(self) -> None:
        errs = validate_wls_brackets("f[x}")
        assert len(errs) == 1
        assert "closes" in errs[0]

    def test_string_literals_skipped(self) -> None:
        assert validate_wls_brackets('Print["unbalanced {"]') == []

    def test_comments_skipped(self) -> None:
        assert validate_wls_brackets("(* unbalanced [ *)") == []

    def test_nested_comments(self) -> None:
        assert validate_wls_brackets("(* outer (* inner *) *)") == []

    def test_escaped_bracket_in_string(self) -> None:
        assert validate_wls_brackets(r'x = "\[Alpha]"') == []

    def test_character_escape_skipped(self) -> None:
        # \[Alpha] should not be counted as an unmatched [
        assert validate_wls_brackets("x = \\[Alpha]") == []

    def test_real_wolfram_code(self) -> None:
        code = """
Module[{lagComp = 0, fieldEquations = {}},
  Do[
    Module[{termComp},
      termComp = DecomposeScalarExpression[lagTerms[[k]], chart, {geH, gea}];
      lagComp += termComp;
    ],
    {k, Length[lagTerms]}
  ];
  Print["Done: ", lagComp];
]
"""
        assert validate_wls_brackets(code) == []

    def test_empty_script(self) -> None:
        assert validate_wls_brackets("") == []

    def test_generated_script_valid(self) -> None:
        """Full integration: generate a real .wls script and validate brackets."""
        import tomllib
        from pathlib import Path

        from tidal.cli._derive import generate_wls

        toml_path = Path("examples/coupled_scalars/theory.toml")
        if not toml_path.exists():
            return  # skip if not available
        with toml_path.open("rb") as f:
            config = tomllib.load(f)
        script = generate_wls(config, None, config_dir=toml_path.parent)
        errs = validate_wls_brackets(script)
        assert errs == [], f"Bracket errors in coupled_scalars .wls: {errs}"
