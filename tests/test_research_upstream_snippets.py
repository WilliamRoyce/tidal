"""The snippet checker for the upstream reports (``scripts/research/perturbations/upstream``).

It needs a Wolfram kernel to run a report, which CI does not have, so what is tested here is
the part that decides what gets checked and how: which trailing comments are exact
expectations and which are descriptions for a person, and that every snippet line reaches the
kernel on its own line, unsplit -- the filed xMAG #2 snippet was wrong precisely because a name
was read on the same line as the ``Needs`` that defines it, and a checker that split that line
would have hidden the error.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CHECKER = (
    REPO / "scripts" / "research" / "perturbations" / "upstream" / "check_snippets.py"
)
SAMPLE = REPO / "tests" / "data" / "upstream_snippets" / "sample.md"


def _load_checker():
    spec = importlib.util.spec_from_file_location("check_snippets", CHECKER)
    assert spec
    assert spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["check_snippets"] = module
    spec.loader.exec_module(module)
    return module


check_snippets = _load_checker()


@pytest.mark.parametrize(
    ("comment", "expected", "exact"),
    [
        ("1", "1", True),
        ("q[a, b]*T[-a]", "q[a, b]*T[-a]", True),
        ("-K[-b, a, -a], unreduced", "-K[-b, a, -a]", True),
        ("still unreduced", "still unreduced", False),
    ],
)
def test_exact_and_descriptive_comments(comment: str, expected: str, exact: bool):
    assert check_snippets.expectation(comment) == (expected, exact)


def test_program_keeps_each_line_whole_and_runs_setup_first():
    program, checks = check_snippets.build_program(SAMPLE)
    lines = program.splitlines()
    assert lines[0] == "DefManifold[M, 4, {a, b, c, d}];"
    # the Needs and the name stay on one line, so the kernel parses them as a reader's would
    assert lines[1].startswith(
        'Print["SNIPPET_OUTPUT 0 ", ToString[Needs["xAct`xTensor`"]; $RiemannSign'
    )
    assert [c.exact for c in checks] == [True, True, True, False]
    assert [c.source for c in checks] == ["sample.md block 1"] * 2 + [
        "sample.md block 2"
    ] * 2
    assert lines[-1] == "f[]"
