"""Helper functions for generating and validating Wolfram Language code strings.

Eliminates fragile triple-brace f-string patterns (``f"{{{x}}}"`` → ``wl_list(x)``)
when constructing .wls scripts in ``_derive.py``.  All functions return plain
strings — no Wolfram execution, no side effects.

Also provides :func:`validate_wls_brackets` for pre-execution syntax checking
of generated scripts (#191).
"""

from __future__ import annotations

from itertools import starmap


def wl_list(*items: str) -> str:
    """Wolfram list literal: ``wl_list("a", "b") → '{a, b}'``."""
    return "{" + ", ".join(items) + "}"


def wl_index(slot: int, chart: str, *, contra: bool = False) -> str:
    """Wolfram basis index: ``wl_index(0, "ch") → '{0, -ch}'``.

    Parameters
    ----------
    slot
        Integer slot index (0-based).
    chart
        Chart symbol name (e.g. ``"geChart"``).
    contra
        If True, omit the minus sign (contravariant index).
    """
    sign = "" if contra else "-"
    return "{" + f"{slot}, {sign}{chart}" + "}"


def wl_rule(lhs: str, rhs: str, *, delayed: bool = False) -> str:
    """Wolfram rule: ``wl_rule("x", "0") → 'x -> 0'``.

    Use ``delayed=True`` for ``:>`` (``RuleDelayed``).
    """
    op = ":>" if delayed else "->"
    return f"{lhs} {op} {rhs}"


def wl_rules_block(rules: list[str]) -> str:
    """Wrap rules in Wolfram ``{rule1, rule2, ...}`` list."""
    return wl_list(*rules)


def wl_component_value(head: str, idx_pairs: list[tuple[int, str]], value: str) -> str:
    """Full ``ComponentValue[tensor[indices], value];`` statement.

    Parameters
    ----------
    head
        Tensor head symbol (e.g. ``"geH"``).
    idx_pairs
        List of ``(slot, chart)`` pairs.  Each produces ``{slot, -chart}``.
    value
        RHS value (e.g. ``"0"``).
    """
    idx_str = ", ".join(starmap(wl_index, idx_pairs))
    return f"ComponentValue[{head}[{idx_str}], {value}];"


def wl_skip_tuples(dim: int) -> str:
    """SkipTuples option for TT gauge: ``"SkipTuples" -> {{0,0},{0,1},...}``."""
    tuples = ", ".join(wl_list(str(0), str(mu)) for mu in range(dim))
    return f'"SkipTuples" -> {wl_list(tuples)}'


def wl_bg_rule_entry(head: str, comps: str, contra: str) -> str:
    """BackgroundFieldRules entry: ``{head, {comps}, {contra}}``."""
    return wl_list(
        head, wl_list(comps) if comps else "{}", wl_list(contra) if contra else "{}"
    )


def wl_flatten_list(items_str: str) -> str:
    """Wrap in Flatten[{...}]: ``wl_flatten_list("a, b") → 'Flatten[{a, b}]'``."""
    return f"Flatten[{{{items_str}}}]"


def wl_diag_matrix(entries: str) -> str:
    """``DiagonalMatrix[{entries}]``."""
    return f"DiagonalMatrix[{{{entries}}}]"


# ---------------------------------------------------------------------------
# Bracket validation (#191)
# ---------------------------------------------------------------------------


def _skip_wl_string(script: str, i: int) -> int:
    """Advance past a Wolfram string literal starting at ``script[i] == '"'``."""
    n = len(script)
    i += 1  # skip opening quote
    while i < n and script[i] != '"':
        if script[i] == "\\" and i + 1 < n:
            i += 2  # skip escaped char
        else:
            i += 1
    return i + 1  # skip closing quote


def _skip_wl_comment(script: str, i: int) -> int:
    """Advance past a (possibly nested) Wolfram comment starting at ``(*``."""
    n = len(script)
    depth = 1
    i += 2  # skip opening (*
    while i + 1 < n and depth > 0:
        if script[i] == "(" and script[i + 1] == "*":
            depth += 1
            i += 2
        elif script[i] == "*" and script[i + 1] == ")":
            depth -= 1
            i += 2
        else:
            i += 1
    return i


def _skip_wl_char_escape(script: str, i: int) -> int:
    r"""Advance past a ``\[Name]`` character escape (e.g. ``\[Alpha]``)."""
    n = len(script)
    i += 2  # skip \[
    while i < n and script[i] != "]":
        i += 1
    return i + 1  # skip closing ]


def validate_wls_brackets(script: str) -> list[str]:
    r"""Check balanced ``[]``, ``{}``, ``()`` in a generated Wolfram script.

    Skips string literals (``"..."``), Wolfram comments (``(* ... *)``),
    and ``\[...]`` character escapes (e.g. ``\[Alpha]``).

    Returns a list of human-readable error strings.  Empty list means valid.
    """
    errors: list[str] = []
    stack: list[tuple[str, int]] = []
    pairs = {"}": "{", "]": "[", ")": "("}
    i = 0
    n = len(script)

    while i < n:
        c = script[i]

        if c == '"':
            i = _skip_wl_string(script, i)
            continue
        if c == "(" and i + 1 < n and script[i + 1] == "*":
            i = _skip_wl_comment(script, i)
            continue
        if c == "\\" and i + 1 < n and script[i + 1] == "[":
            i = _skip_wl_char_escape(script, i)
            continue

        if c in "{[(":
            stack.append((c, i))
        elif c in "}])":
            if not stack:
                line = script[:i].count("\n") + 1
                errors.append(f"line {line}: unmatched closing '{c}'")
            else:
                open_char, open_pos = stack.pop()
                if pairs[c] != open_char:
                    open_line = script[:open_pos].count("\n") + 1
                    close_line = script[:i].count("\n") + 1
                    errors.append(
                        f"line {close_line}: '{c}' closes '{open_char}' "
                        f"opened at line {open_line}"
                    )
        i += 1

    for open_char, pos in stack:
        line = script[:pos].count("\n") + 1
        errors.append(f"line {line}: unclosed '{open_char}'")

    return errors


def wl_zero_component(comp_name: str, field_name: str, gauge_type: str) -> list[str]:
    """Generate WLS to substitute a component and its derivatives with zero.

    Returns a list of WLS code lines (same contract as ``_wls_*`` functions).
    """
    return [
        f"(* {gauge_type.capitalize()} gauge: {field_name} component {comp_name} = 0 *)",
        f"fieldEquations = fieldEquations /. "
        f"{wl_list(wl_rule(f'{comp_name}[args___]', '0', delayed=True), wl_rule(f'Derivative[ders__][{comp_name}][args___]', '0', delayed=True))};",
        f'Print["Applied {gauge_type} gauge: {comp_name} = 0"];',
        "",
    ]
