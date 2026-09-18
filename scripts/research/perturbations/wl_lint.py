#!/usr/bin/env python3
"""Pre-flight check for the Wolfram scripts in this directory -- run before any kernel starts.

Mathematica rarely reports the mistakes this lane has actually made. A call that matches no
definition comes back unevaluated; a definition that spans lines can be cut short at the first
line break; a comment containing the closing delimiter truncates the rest of the file; a name
read on the same line as the ``Needs`` that defines it is created in the wrong context. Each
of those produced a plausible-looking wrong answer rather than an error (see the memo's
section 8 and GH #591). This script reads the source, without a kernel, and refuses to let
such a file reach one.

Checks, each with the defect code it reports:

- ``BALANCE``   brackets, braces, parentheses, strings and comments are balanced.
- ``ARITY``     every call of a helper defined in the file -- or in the shared harness
                ``RCSetup.wl`` / ``RCSetupCore.wl`` beside it -- passes an argument count one
                of its definitions accepts (found f7's four-argument calls to a three-argument
                helper, the source of the retracted "silent empty first order").
- ``MULTILINE`` a top-level line that starts with a binary operator after a line that already
                completed a statement: the continuation is silently dropped.
- ``NEEDS``     a top-level ``Needs[...]`` followed on the same line by another statement,
                which is parsed before the package is loaded.
- ``SHADOW``    a bare ``$Version``. Every xAct package exports its own, so once one is loaded
                the bare name is that package's version, not the kernel's: measured, it is
                xPand's after xPand loads and xMAG's after xMAG loads. The harness printed
                xPand's version under the label ``WOLFRAM_VERSION`` for that reason. Write
                ``System`$Version`` or the package's qualified name.

Usage: ``python3 wl_lint.py FILE...`` (exit 1 if any defect is found), or with no arguments
to check every ``*.wls`` and ``*.wl`` under ``wolfram/`` and every ``*.wls`` under ``upstream/``.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

HERE = Path(__file__).resolve().parent
HARNESS = ("RCSetup.wl", "RCSetupCore.wl")
OPEN, CLOSE = "[({", "])}"
PAIR = dict(zip(CLOSE, OPEN, strict=True))


@dataclass(frozen=True)
class Defect:
    """One finding: the file, the 1-based line, a defect code, and what is wrong."""

    path: Path
    line: int
    code: str
    message: str

    def __str__(self) -> str:
        return f"{self.path}:{self.line}: {self.code}: {self.message}"


# ------------------------------------------------------------------------------------------
# Lexing: blank out comments and replace string contents, keeping line numbers intact
# ------------------------------------------------------------------------------------------


def neutralize(text: str) -> tuple[str, list[tuple[int, str]]]:
    """Return the source with comments blanked and strings reduced to ``"S"``.

    Newlines are preserved so line numbers survive. Also returns balance problems found while
    lexing (unterminated comment or string) as ``(line, message)`` pairs.
    """
    out: list[str] = []
    problems: list[tuple[int, str]] = []
    i, n, line = 0, len(text), 1
    while i < n:
        if text.startswith("(*", i):
            depth, start_line = 1, line
            i += 2
            while i < n and depth:
                if text.startswith("(*", i):
                    depth += 1
                    i += 2
                elif text.startswith("*)", i):
                    depth -= 1
                    i += 2
                else:
                    if text[i] == "\n":
                        line += 1
                        out.append("\n")
                    else:
                        out.append(" ")
                    i += 1
            if depth:
                problems.append((start_line, "unterminated comment"))
            continue
        c = text[i]
        if c == '"':
            j, start_line, newlines = i + 1, line, 0
            while j < n and text[j] != '"':
                if text[j] == "\n":
                    newlines += 1
                j += 2 if text[j] == "\\" else 1
            if j >= n:
                problems.append((start_line, "unterminated string"))
            out.append('"S"' + "\n" * newlines)
            line += newlines
            i = j + 1
            continue
        if c == "\n":
            line += 1
        out.append(c)
        i += 1
    return "".join(out), problems


def check_balance(src: str) -> list[tuple[int, str]]:
    """Report unmatched or mismatched brackets in comment- and string-free source."""
    stack: list[tuple[str, int]] = []
    problems: list[tuple[int, str]] = []
    line = 1
    for c in src:
        if c == "\n":
            line += 1
        elif c in OPEN:
            stack.append((c, line))
        elif c in CLOSE:
            if not stack:
                problems.append((line, f"unmatched {c}"))
            elif stack[-1][0] != PAIR[c]:
                opener, opened = stack.pop()
                problems.append((line, f"{c} closes {opener} opened on line {opened}"))
            else:
                stack.pop()
    problems.extend((opened, f"{opener} never closed") for opener, opened in stack[-5:])
    return problems


# ------------------------------------------------------------------------------------------
# Arity: definitions and call sites of locally defined helpers
# ------------------------------------------------------------------------------------------

NAME_CALL = re.compile(r"(?<![\w`$@])([A-Za-z$][\w$]*)\[")
DEFINITION_OPERATOR = re.compile(r"\s*(:=|=(?![=.!]))")


def split_args(src: str, start: int) -> tuple[list[str], int]:
    """Split the arguments of a call whose ``[`` ends just before ``start``.

    Returns the top-level argument strings and the index of the closing ``]``.
    """
    depth, i, current, args = 1, start, [], []
    while i < len(src):
        c = src[i]
        if c in OPEN:
            depth += 1
        elif c in CLOSE:
            depth -= 1
            if depth == 0:
                break
        if depth == 1 and c == ",":
            args.append("".join(current))
            current = []
        else:
            current.append(c)
        i += 1
    args.append("".join(current))
    if len(args) == 1 and not args[0].strip():
        args = []
    return args, i


@dataclass(frozen=True)
class Signature:
    """The argument counts one definition accepts."""

    low: int
    high: int
    sequence: bool

    def accepts(self, count: int) -> bool:
        return self.sequence or self.low <= count <= self.high


def definitions(src: str) -> dict[str, list[Signature]]:
    """Collect ``f[pattern, ...] :=`` and ``f[pattern, ...] =`` definitions by name."""
    found: dict[str, list[Signature]] = {}
    for match in NAME_CALL.finditer(src):
        args, end = split_args(src, match.end())
        if not DEFINITION_OPERATOR.match(src, end + 1):
            continue
        if not args or not all("_" in a for a in args):
            continue
        sequence = any("__" in a for a in args)
        # an optional argument is x_:default, x_Head:default or x_. (not x_ :> or x_ :=)
        optional = sum(1 for a in args if re.search(r"_[A-Za-z`$]*\s*:(?![=>])|_\.", a))
        found.setdefault(match.group(1), []).append(
            Signature(len(args) - optional, len(args), sequence)
        )
    return found


def check_arity(src: str, known: dict[str, list[Signature]]) -> list[tuple[int, str]]:
    """Report calls of known helpers whose argument count no definition accepts."""
    problems: list[tuple[int, str]] = []
    for name, signatures in known.items():
        for match in re.finditer(r"(?<![\w`$@])" + re.escape(name) + r"\[", src):
            args, end = split_args(src, match.end())
            if DEFINITION_OPERATOR.match(src, end + 1) and all("_" in a for a in args):
                continue
            if not any(s.accepts(len(args)) for s in signatures):
                accepted = sorted({(s.low, s.high) for s in signatures})
                problems.append(
                    (
                        src.count("\n", 0, match.start()) + 1,
                        f"{name} called with {len(args)} argument(s); definitions accept {accepted}",
                    )
                )
    return problems


# ------------------------------------------------------------------------------------------
# Line-level checks: multi-line continuations and same-line Needs
# ------------------------------------------------------------------------------------------

CONTINUATION = re.compile(r"^\s*(\+|-\s|\*(?!\))|/(?![/.@;])|&&|\|\||<>|->|:>)")
OPEN_ENDING = re.compile(r"(\+|-|\*|/|&&|\|\||<>|->|:>|:=|=|,|\[|\(|\{|@|/@|//)\s*$")


def top_level_depths(src: str) -> list[int]:
    """Bracket depth at the start of each line of comment- and string-free source."""
    depths, depth = [0], 0
    for c in src:
        if c in OPEN:
            depth += 1
        elif c in CLOSE:
            depth = max(depth - 1, 0)
        elif c == "\n":
            depths.append(depth)
    return depths


def check_lines(src: str) -> list[tuple[int, str, str]]:
    """Report dropped continuations (MULTILINE) and same-line ``Needs`` (NEEDS)."""
    lines = src.split("\n")
    depths = top_level_depths(src)
    problems: list[tuple[int, str, str]] = []
    previous = ""
    for number, (text, depth) in enumerate(zip(lines, depths, strict=False), start=1):
        stripped = text.strip()
        if (
            depth == 0
            and CONTINUATION.match(text)
            and previous
            and not OPEN_ENDING.search(previous)
        ):
            problems.append(
                (
                    number,
                    "MULTILINE",
                    "a top-level line starts with an operator after a complete statement; "
                    "the continuation is dropped -- wrap the whole right-hand side in ( ... )",
                )
            )
        needs = re.match(r"^\s*Needs\[", text)
        if depth == 0 and needs:
            _, end = split_args(text, needs.end())
            rest = text[end + 1 :]
            if rest.strip().lstrip(";").strip():
                problems.append(
                    (
                        number,
                        "NEEDS",
                        "Needs[...] shares its line with another statement, which is parsed "
                        "before the package loads -- put Needs on a line of its own",
                    )
                )
        if stripped:
            previous = stripped
    return problems


SHADOWED = re.compile(r"(?<![\w`$])\$Version\b")


def check_shadowed(src: str) -> list[tuple[int, str]]:
    """Report bare names that a loaded package silently redefines (``$Version``)."""
    return [
        (
            src.count("\n", 0, m.start()) + 1,
            "bare $Version is the last-loaded xAct package's version, not the kernel's -- "
            "write System`$Version or the package's qualified name",
        )
        for m in SHADOWED.finditer(src)
    ]


# ------------------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------------------


def harness_definitions(directory: Path) -> dict[str, list[Signature]]:
    """Definitions from the shared harness beside ``directory``, used by every script."""
    known: dict[str, list[Signature]] = {}
    for name in HARNESS:
        path = directory / name
        if path.exists():
            for key, value in definitions(neutralize(path.read_text())[0]).items():
                known.setdefault(key, []).extend(value)
    return known


def lint(path: Path) -> list[Defect]:
    """All defects in one Wolfram source file."""
    text = path.read_text()
    src, lex_problems = neutralize(text)
    defects = [Defect(path, line, "BALANCE", msg) for line, msg in lex_problems]
    defects += [Defect(path, line, "BALANCE", msg) for line, msg in check_balance(src)]
    known = harness_definitions(path.parent)
    for key, value in definitions(src).items():
        known.setdefault(key, []).extend(value)
    defects += [
        Defect(path, line, "ARITY", msg) for line, msg in check_arity(src, known)
    ]
    defects += [Defect(path, line, code, msg) for line, code, msg in check_lines(src)]
    defects += [Defect(path, line, "SHADOW", msg) for line, msg in check_shadowed(src)]
    return sorted(defects, key=lambda d: (str(d.path), d.line, d.code))


def default_targets() -> list[Path]:
    """Every Wolfram script in this lane, and every script attached to an upstream report."""
    wolfram = HERE / "wolfram"
    return sorted(
        [
            *wolfram.glob("*.wls"),
            *wolfram.glob("*.wl"),
            *(HERE / "upstream").glob("*.wls"),
        ]
    )


def main(argv: list[str]) -> int:
    targets = [Path(a) for a in argv] or default_targets()
    defects = [d for path in targets for d in lint(path)]
    for defect in defects:
        print(defect)
    print(f"wl_lint: {len(targets)} file(s), {len(defects)} defect(s)")
    return 1 if defects else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
