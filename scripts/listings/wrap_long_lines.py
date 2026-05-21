r"""Wrap long generated-LaTeX equation lines onto multiple aligned lines.

This script is a thin post-processor for the output of `tidal inspect
--latex --latex-format align` (after the `\\begin{align}` / `\\end{align}`
wrapper has been stripped by the driver). It rewrites any line of the
form

    LHS &= term_1 +/- term_2 +/- ... +/- term_N \\\\

into a multi-line aligned chain

    LHS &= term_1 +/- term_2 +/- ... \\\\
        &\\quad +/- term_k +/- ... \\\\
        ...

whenever the original line exceeds a configurable character threshold.
The split is taken on the top-level ` + ` / ` - ` boundaries only;
operators inside `\\frac{a + b}{c}`, `\\tensor{T}{_a + b}`, or any
other braced grouping are left untouched.

Designed for inclusion via `\\input{...}` inside an `aligned` math
environment housed in a `figure*` float (see Appendix E plan).
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterable


def _split_top_level_terms(rhs: str) -> list[tuple[str, str]]:
    r"""Walk `rhs` and split it into (sign, term-body) pairs.

    The walk tracks nesting depth via `{` / `}` and `[` / `]`. At depth 0,
    every ` + ` and ` - ` is a candidate split point; at non-zero depth
    those tokens belong to a grouped construct and are kept inside the
    current term.

    The first term gets sign `+` unless the RHS begins with a `-`. Backslash
    commands (e.g. `\\frac`, `\\tensor`) are treated as opaque text — only
    `{...}` and `[...]` change the depth counter.
    """
    terms: list[tuple[str, str]] = []
    depth = 0
    last_split = 0
    i = 0
    # If the RHS starts with a leading `-`, fold it into the first term sign.
    first_sign = "+"
    rhs_stripped = rhs.lstrip()
    leading_offset = len(rhs) - len(rhs_stripped)
    if rhs_stripped.startswith("-"):
        first_sign = "-"
        last_split = leading_offset + 1  # skip the `-` itself
    elif rhs_stripped.startswith("+"):
        last_split = leading_offset + 1
    while i < len(rhs):
        ch = rhs[i]
        if ch in {"{", "["}:
            depth += 1
        elif ch in {"}", "]"}:
            depth -= 1
        elif depth == 0 and i + 1 < len(rhs) and ch == " " and rhs[i + 1] in "+-":
            # Top-level term boundary. The generator format is asymmetric:
            #   ` + <term>` always has a space after the plus, but
            #   ` -<term>` has *no* space after the minus (e.g. ` -\frac{…}`
            #   or ` -2 \, \zeta_1`).
            # Accept both forms; for ` - ` (with trailing space) also skip
            # the trailing space when advancing past the operator.
            sign_ch = rhs[i + 1]
            has_trailing_space = i + 2 < len(rhs) and rhs[i + 2] == " "
            if sign_ch == "+" and not has_trailing_space:
                # `+` without trailing space is unusual; treat conservatively
                # by NOT splitting here. (Skip past the space and continue.)
                i += 1
                continue
            term_body = rhs[last_split:i].strip()
            terms.append((first_sign, term_body))
            first_sign = sign_ch
            skip = 3 if has_trailing_space else 2
            last_split = i + skip
            i += skip
            continue
        i += 1
    # Tail term (whatever remains after the last split point).
    tail = rhs[last_split:].strip()
    if tail:
        terms.append((first_sign, tail))
    return terms


# Match either `&= RHS` (equations of motion, Lagrangians) or
# `&\supset RHS` (Hamiltonian densities restricted to the self-GW + self-EM
# sector; cf. tidal/symbolic/latex.py::hamiltonian_to_latex). The operator
# is captured so the emission preserves whichever the source line used.
_LHS_RE = re.compile(
    r"^(?P<lhs>.*?)\s*&(?P<op>=|\\supset)\s*(?P<rhs>.*?)(?P<trailer>\s*\\\\\s*)?$"
)


def _join_terms(signed_terms: Iterable[tuple[str, str]]) -> str:
    """Render `[(+, t1), (-, t2), ...]` as the original infix string."""
    out: list[str] = []
    for idx, (sign, body) in enumerate(signed_terms):
        if idx == 0:
            if sign == "-":
                out.append(f"-{body}")
            else:
                out.append(body)
        else:
            out.append(f" {sign} {body}")
    return "".join(out)


def wrap_line(line: str, width: int) -> str:
    """Return a (possibly multi-line) replacement for `line`.

    Short lines and lines that do not match the `LHS &= RHS` shape pass
    through unchanged.
    """
    if len(line.rstrip("\n")) <= width:
        return line
    m = _LHS_RE.match(line.rstrip("\n"))
    if not m:
        return line
    lhs = m.group("lhs")
    op = m.group("op")  # "=" or r"\supset"
    rhs = m.group("rhs")
    trailer = m.group("trailer") or ""
    terms = _split_top_level_terms(rhs)
    if len(terms) <= 1:
        # Single term that's too long — nothing we can split on cleanly.
        return line
    # Greedy fill: pack terms onto the first line up to ~width chars,
    # then start a continuation; repeat. Always emit at least one term
    # per physical line so we make progress even when a single term is
    # near the width threshold.
    chunks: list[list[tuple[str, str]]] = [[]]
    # First-line budget: `LHS &<op> ` (LHS + 4 for "&= " or LHS + len("&\supset ")
    # for the Hamiltonian \supset form).
    chunks_first_prefix = f"{lhs} &{op} "
    chunk_lengths: list[int] = [len(chunks_first_prefix)]
    for sign, body in terms:
        candidate = (
            f" {sign} {body}" if chunks[-1] else (f"-{body}" if sign == "-" else body)
        )
        if chunks[-1] and chunk_lengths[-1] + len(candidate) > width:
            chunks.append([])
            chunk_lengths.append(len("  &\\quad ") + len(f"{sign} {body}"))
        else:
            chunk_lengths[-1] += len(candidate)
        chunks[-1].append((sign, body))
    # Emit. Preserve the source operator (= for L/EOMs, \supset for H).
    lines: list[str] = []
    first = chunks[0]
    lines.append(f"{lhs} &{op} {_join_terms(first)} \\\\")
    for chunk in chunks[1:]:
        # Each continuation line starts with the sign of the first term.
        head_sign, head_body = chunk[0]
        rest = chunk[1:]
        head = f"  &\\quad {head_sign} {head_body}"
        _join_terms([("+", "")][:0] + [(s, b) for s, b in rest])
        # `_join_terms` prefixes the first item with no operator; for
        # continuation lines we want ` <sign> <body>` for every item after
        # the head — re-render manually:
        cont_parts = [head]
        for s, b in rest:
            cont_parts.append(f" {s} {b}")
        lines.append("".join(cont_parts) + " \\\\")
    # Restore the original trailing `\\\\` only on the *last* emitted line;
    # we have already appended ` \\\\` to every line above, so drop the
    # extra one added at the end if the source had no trailer.
    if not trailer.strip():
        # No trailing backslashes in the input — strip from our output.
        lines[-1] = lines[-1].removesuffix(" \\\\")
    return "\n".join(lines) + "\n"


def process_stream(infile, outfile, width: int) -> None:
    for raw in infile:
        outfile.write(wrap_line(raw, width))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--width",
        type=int,
        default=120,
        help="character threshold above which a line is split (default: 120)",
    )
    parser.add_argument(
        "input",
        nargs="?",
        type=argparse.FileType("r"),
        default=sys.stdin,
        help="input .tex file (default: stdin)",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="output .tex file (default: stdout)",
    )
    args = parser.parse_args(argv)
    process_stream(args.input, args.output, args.width)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
