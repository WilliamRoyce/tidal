r"""Unit tests for scripts/listings/wrap_long_lines.py.

Covers the six cases specified in Fix 1 of the Appendix E plan:

(i)   short lines pass through unchanged;
(ii)  a depth-0 ` + ` is split;
(iii) a ` + ` inside `\\frac{a + b}{c}` is NOT split;
(iv)  a ` + ` inside `\\tensor{T}{_a + b}` is NOT split;
(v)   the LHS `&=` alignment column is preserved on the first line;
(vi)  continuation lines start with `  &\\quad`.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

# Load the wrap module from its repo-relative path. `scripts/listings/`
# is not on `sys.path`, so importlib does the work explicitly.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_WRAP_PATH = _REPO_ROOT / "scripts" / "listings" / "wrap_long_lines.py"
_spec = importlib.util.spec_from_file_location("wrap_long_lines", _WRAP_PATH)
assert _spec is not None
assert _spec.loader is not None
wrap_long_lines = importlib.util.module_from_spec(_spec)
sys.modules["wrap_long_lines"] = wrap_long_lines
_spec.loader.exec_module(wrap_long_lines)


def test_short_line_passes_through_unchanged() -> None:
    """(i) Short lines are returned verbatim, including trailing newline."""
    line = "\\mathcal{L} &= -\\tfrac{1}{4} F^2 \\\\\n"
    assert wrap_long_lines.wrap_line(line, width=120) == line


def test_depth_zero_plus_minus_is_split() -> None:
    """(ii) ` + ` and ` - ` at depth 0 produce a multi-line wrap."""
    # Long enough that the 120-char threshold is exceeded.
    body = " + ".join(
        f"\\alpha_{{{i}}} \\partial_z^2 \\mathcal{{T}}_{{{i}}}" for i in range(12)
    )
    line = f"\\partial_t^{{2}} \\mathcal{{A}}_{{0}} &= {body} \\\\\n"
    assert len(line) > 120, "test setup: line must exceed the threshold"
    out = wrap_long_lines.wrap_line(line, width=120)
    # More than one physical line in the output
    assert out.count("\n") >= 2, f"expected wrap; got:\n{out}"


def test_plus_inside_frac_is_not_split() -> None:
    r"""(iii) ` + ` inside `\\frac{a + b}{c}` does not become a split point."""
    rhs = "\\frac{\\alpha + \\beta}{\\gamma}"
    terms = wrap_long_lines._split_top_level_terms(rhs)
    assert terms == [("+", rhs)], (
        f"expected single term, got {terms!r} (split inside \\frac{{}}{{}} is forbidden)"
    )


def test_plus_inside_tensor_indices_is_not_split() -> None:
    r"""(iv) ` + ` inside `\\tensor{T}{_a + b}` does not become a split point."""
    rhs = "\\tensor{T}{_a + b} \\tensor{g}{^a^c}"
    terms = wrap_long_lines._split_top_level_terms(rhs)
    # One " + " sits inside `{...}`, the other " " is between tensor groups —
    # neither is a top-level operator, so the result is a single term.
    assert terms == [("+", rhs)], (
        f"expected single term, got {terms!r} (split inside \\tensor{{}}{{}} is forbidden)"
    )


def test_lhs_equals_alignment_preserved_on_first_line() -> None:
    """(v) The first emitted line still has the `LHS &= …` alignment column."""
    body = " + ".join(f"\\alpha_{{{i}}} X_{{{i}}}" for i in range(20))
    line = f"\\partial_t^{{2}} A &= {body} \\\\\n"
    out = wrap_long_lines.wrap_line(line, width=80)
    first_line = out.splitlines()[0]
    assert " &= " in first_line, f"first line lost the `&=` column: {first_line!r}"
    assert first_line.lstrip().startswith("\\partial_t^{2} A"), first_line


def test_continuation_lines_use_amp_quad() -> None:
    r"""(vi) Continuation lines start with `  &\\quad`."""
    body = " + ".join(f"\\alpha_{{{i}}} X_{{{i}}}" for i in range(20))
    line = f"\\partial_t^{{2}} A &= {body} \\\\\n"
    out = wrap_long_lines.wrap_line(line, width=80)
    continuations = out.splitlines()[1:]
    assert continuations, "expected at least one continuation line"
    for cont in continuations:
        assert cont.startswith("  &\\quad "), f"bad continuation prefix: {cont!r}"


def test_trailing_backslashes_round_trip() -> None:
    r"""Every emitted line ends with ` \\\\` (the align row terminator)."""
    body = " + ".join(f"\\alpha_{{{i}}} X_{{{i}}}" for i in range(20))
    line = f"\\partial_t^{{2}} A &= {body} \\\\\n"
    out = wrap_long_lines.wrap_line(line, width=80).rstrip("\n")
    for emitted in out.split("\n"):
        assert emitted.rstrip().endswith("\\\\"), (
            f"emitted line missing trailing align terminator: {emitted!r}"
        )


def test_leading_minus_is_consumed_into_first_term_sign() -> None:
    """An RHS that starts with `-` produces a negative first term, not a stray sign."""
    rhs = "-B_{0} \\partial_z \\mathcal{H}_{2} + \\partial_z^2 \\mathcal{A}_{0}"
    terms = wrap_long_lines._split_top_level_terms(rhs)
    assert terms == [
        ("-", "B_{0} \\partial_z \\mathcal{H}_{2}"),
        ("+", "\\partial_z^2 \\mathcal{A}_{0}"),
    ], terms


# ---------------------------------------------------------------------------
# Issue #371 regression guard: no \dot{\mathcal{X}_{N}} or \dot{\mathcal{X}^{N}}
# in generated Appendix-E listings. The `accents` package (loaded in
# manuscript/macros.tex for \accentset) leaks \mathcal font scope into a
# sub/super-script that sits *inside* a \dot{...} argument, emitting
# "Missing character: There is no <digit> in font stix-mathcal" warnings.
# tidal/symbolic/latex.py::_wrap_dot lifts trailing _{...} / ^{...} groups
# outside the dot; this check fails fast if a regenerated listing
# reintroduces the pattern.
# ---------------------------------------------------------------------------

import re  # noqa: E402

_LISTINGS_DIR = _REPO_ROOT / "manuscript" / "sections" / "appendices" / "listings"
_DOT_SUPSUB_RE = re.compile(r"\\dot\{[^{}]*(?:\{[^{}]*\}[^{}]*)*[_^]\{[^{}]+\}\}")


def test_no_dot_wraps_subscript_or_superscript() -> None:
    r"""Guard against ``\dot{...}`` wrapping a trailing ``_{N}`` / ``^{N}``.

    The accents package leaks ``\mathcal`` font scope into a subscript that
    sits inside a ``\dot{...}`` argument (issue #371). Listings must lift
    such subscripts outside the dot.
    """
    if not _LISTINGS_DIR.exists():
        return  # Listings dir absent in some checkouts; nothing to guard.
    offenders: list[tuple[str, int, str]] = []
    for path in sorted(_LISTINGS_DIR.glob("eom_*_full.tex")):
        for lineno, line in enumerate(path.read_text().splitlines(), start=1):
            if _DOT_SUPSUB_RE.search(line):
                offenders.append((path.name, lineno, line.strip()))
    assert not offenders, (
        "Appendix-E listings contain \\dot{...} wrapping a trailing "
        "_{...} or ^{...} (issue #371). Regenerate via "
        "`bash scripts/listings/render_app_e_symbolic.sh`. Offenders:\n"
        + "\n".join(f"  {n}:{ln}: {body}" for n, ln, body in offenders[:10])
    )
