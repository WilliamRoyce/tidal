"""The pre-flight lint for the R-C Wolfram lane scripts (GH #591).

Mathematica returns a wrong call unevaluated instead of failing, so the defects this lint looks
for surface as plausible wrong answers rather than errors -- one of them produced a claim that
had to be retracted. Two properties are tested: every committed lane script is clean, and the
lint is seen to fail on a deliberately broken file for each defect it claims to catch (a check
that has never been watched failing is not yet a check).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LINT_PATH = REPO / "scripts" / "research" / "perturbations" / "wl_lint.py"
FIXTURES = REPO / "tests" / "data" / "wl_lint"


def _load_lint():
    spec = importlib.util.spec_from_file_location("wl_lint", LINT_PATH)
    assert spec
    assert spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["wl_lint"] = module
    spec.loader.exec_module(module)
    return module


wl_lint = _load_lint()


def test_every_lane_script_is_clean():
    targets = wl_lint.default_targets()
    assert len(targets) >= 18, "the lane's scripts were not found"
    defects = [str(d) for path in targets for d in wl_lint.lint(path)]
    assert defects == []


@pytest.mark.parametrize(
    ("fixture", "code"),
    [
        ("arity.wls", "ARITY"),
        ("balance.wls", "BALANCE"),
        ("multiline.wls", "MULTILINE"),
        ("needs.wls", "NEEDS"),
        ("shadow.wls", "SHADOW"),
    ],
)
def test_each_defect_is_caught(fixture: str, code: str):
    codes = {d.code for d in wl_lint.lint(FIXTURES / fixture)}
    assert code in codes, f"{fixture} should report {code}, reported {sorted(codes)}"


def test_correct_forms_are_not_flagged():
    """Optional arguments (typed or not), a parenthesized multi-line definition, a Needs on its
    own line and a Needs nested inside Check are all correct and must pass.
    """
    assert [str(d) for d in wl_lint.lint(FIXTURES / "good.wls")] == []


def test_command_line_exit_status(capsys):
    assert wl_lint.main([str(FIXTURES / "good.wls")]) == 0
    assert wl_lint.main([str(FIXTURES / "arity.wls")]) == 1
    assert "ARITY" in capsys.readouterr().out
