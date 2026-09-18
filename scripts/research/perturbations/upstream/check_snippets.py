#!/usr/bin/env python3
"""Run every code snippet in the upstream reports the way a reader would, and check its outputs.

Two reports filed from this lane contained a snippet that did not do what its comment said
(xMAG #2, item 1: a name read on the same line as the ``Needs`` that defines it). A report is
not ready until every snippet in it has been run as written and every ``(* expected *)``
comment compared with what actually came back. This script does that.

Markdown reports (``*.md``): every fenced ``wolfram`` block in a file is run, in order, in ONE
fresh kernel -- as a reader working down the report would. An HTML comment directly before a
block, ``<!-- setup: ... -->``, holds lines the report assumes but does not show; they are run
first and are invisible on GitHub. Each line ending in a comment is checked:

- an *exact* comment (``(* g *)``, ``(* -1 *)``, ``(* -ContorsionCDT[-b, a, -a], unreduced *)``
  -- the part before the first comma, when it is not an English phrase) must equal the
  line's output in ``InputForm``, ignoring spaces;
- a *descriptive* comment (``(* still unreduced *)``) is printed beside the actual output for a
  person to read, and never counts as a pass.

Attached scripts (``*.wls``) are run with ``wolframscript -file``. Their ``control`` checks must
print ``ok`` and their ``reported`` checks ``reproduced`` -- no ``MISMATCH``, no
``NOT REPRODUCED`` -- and they must end with ``DONE``.

Refuses to start if a Wolfram kernel is already running (one engine license). Exit status 1 on
any mismatch or failure. Usage: ``python3 check_snippets.py [FILE ...]`` (default: every report
in this directory).
"""

from __future__ import annotations

import os
import re
import subprocess  # noqa: S404 -- running a Wolfram kernel is this script's purpose
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
FENCE = re.compile(
    r"(?:<!--\s*setup:(?P<setup>.*?)-->\s*\n)?```wolfram\n(?P<body>.*?)```", re.DOTALL
)
TRAILING_COMMENT = re.compile(r"^(?P<code>.*?)\s*\(\*\s*(?P<comment>.*?)\s*\*\)\s*$")
ENGLISH = re.compile(r"[A-Za-z]{2,} [A-Za-z]{2,}")
MARK = "SNIPPET_OUTPUT"


@dataclass
class Check:
    """One commented line: where it is, what it says to expect, whether that is exact."""

    source: str
    line: str
    expected: str
    exact: bool


def top_level_head(comment: str) -> str:
    """The comment up to its first comma outside any brackets: ``q[a, b]*T[-a]`` stays whole,
    ``-ContorsionCDT[-b, a, -a], unreduced`` loses only its trailing note.
    """
    depth = 0
    for i, c in enumerate(comment):
        if c in "[({":
            depth += 1
        elif c in "])}":
            depth -= 1
        elif c == "," and depth == 0:
            return comment[:i].strip()
    return comment.strip()


def expectation(comment: str) -> tuple[str, bool]:
    """The expected output in a trailing comment, and whether it is exact or descriptive."""
    head = top_level_head(comment)
    if ENGLISH.search(head):
        return comment, False
    return head, True


def kernel_busy() -> bool:
    for name in ("wolframscript", "WolframKernel"):
        if (
            subprocess.run(
                ["pgrep", "-x", name], check=False, capture_output=True
            ).returncode
            == 0
        ):
            return True
    return False


def run_code(code: str, timeout: int = 900) -> str:
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    result = subprocess.run(
        ["wolframscript", "-code", code],
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
    )
    return result.stdout + result.stderr


def build_program(path: Path) -> tuple[str, list[Check]]:
    """Turn a report's snippets into one program; commented lines print their output."""
    text = path.read_text()
    lines: list[str] = []
    checks: list[Check] = []
    for block_number, match in enumerate(FENCE.finditer(text), start=1):
        if match.group("setup"):
            lines.extend(
                s for s in match.group("setup").strip().splitlines() if s.strip()
            )
        for raw in match.group("body").splitlines():
            m = TRAILING_COMMENT.match(raw)
            if not m or not m.group("code").strip():
                lines.append(raw)
                continue
            code = m.group("code").rstrip().rstrip(";")
            expected, exact = expectation(m.group("comment"))
            index = len(checks)
            checks.append(
                Check(
                    f"{path.name} block {block_number}", code.strip(), expected, exact
                )
            )
            # The line is kept whole, so it is parsed exactly as the reader's would be.
            lines.append(f'Print["{MARK} {index} ", ToString[{code}, InputForm]]')
    return "\n".join(lines), checks


def check_markdown(path: Path) -> bool:
    program, checks = build_program(path)
    if not checks:
        print(f"{path.name}: no commented lines to check")
        return True
    output = run_code(program)
    actual: dict[int, str] = {}
    for line in output.splitlines():
        if line.startswith(MARK + " "):
            _, index, value = line.split(" ", 2)
            actual[int(index)] = value
    ok = True
    for index, check in enumerate(checks):
        got = actual.get(index, "<no output: the line failed or never ran>")
        if check.exact:
            match = got.replace(" ", "") == check.expected.replace(" ", "")
            ok &= match
            status = "ok" if match else "MISMATCH"
            print(
                f"{status:8} {check.source}: {check.line}\n         expected {check.expected!r}, got {got!r}"
            )
        else:
            print(
                f"{'READ':8} {check.source}: {check.line}\n         says {check.expected!r}; got {got!r}"
            )
    return ok


def check_script(path: Path) -> bool:
    env = {**os.environ, "QT_QPA_PLATFORM": "offscreen"}
    result = subprocess.run(
        ["wolframscript", "-file", str(path)],
        check=False,
        capture_output=True,
        text=True,
        timeout=900,
        env=env,
    )
    output = result.stdout + result.stderr
    reported = [ln for ln in output.splitlines() if ln.startswith(("CHECK", "DONE"))]
    for line in reported:
        print(f"{path.name}: {line}")
    bad = any("MISMATCH" in ln or "NOT REPRODUCED" in ln for ln in reported)
    return not bad and any(ln.startswith("DONE") for ln in reported)


def main(argv: list[str]) -> int:
    if kernel_busy():
        print("check_snippets: a Wolfram kernel is running; refusing (single license)")
        return 3
    targets = [Path(a) for a in argv] or sorted(
        [*HERE.glob("*.md"), *HERE.glob("*.wls")]
    )
    targets = [
        t for t in targets if t.name != "README.md" and not t.name.startswith("claims")
    ]
    results = {}
    for target in targets:
        print(f"=== {target.name} ===")
        results[target.name] = (
            check_script(target) if target.suffix == ".wls" else check_markdown(target)
        )
    failed = [name for name, ok in results.items() if not ok]
    print(
        "check_snippets: "
        + ("all as expected" if not failed else f"FAILED: {', '.join(failed)}")
    )
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
