"""Read-only BibTeX field extractor for the citation audit.

This module NEVER writes to the manuscript. It performs simple, tolerant text
extraction of entries and their fields so they can be compared against canonical
records fetched from INSPIRE-HEP / Crossref. All mutation of ``references.bib``
happens elsewhere via the Edit tool (exact string replacement) — see
``scripts/bibaudit/README.md`` and the audit plan.

The extractor handles both delimiter styles present in this project:
  * local entries use ``field = {value}`` with nested braces / LaTeX accents
  * INSPIRE bibtex exports use ``field = "value"``
and bare numeric values (``year = 1962``). It also tracks the ``% === ... ===``
section header each entry falls under, which drives audit batching.
"""

from __future__ import annotations

import itertools
import re
import sys
from dataclasses import dataclass
from dataclasses import field as dc_field
from pathlib import Path

# --- entry / field tokenizer ---------------------------------------------------

_ENTRY_START = re.compile(r"@(?P<type>[A-Za-z]+)\s*\{\s*(?P<key>[^,\s]+)\s*,")
_RULE = re.compile(r"^%\s*=+\s*$")
_BOILERPLATE = re.compile(
    r"all entries verified|^bibtex from|^cited|^consolidated", re.IGNORECASE
)


def _find_sections(text: str) -> list[tuple[int, str]]:
    """Locate ``% === / % Title / % ===`` banner blocks.

    Returns a sorted list of ``(char_offset, title)`` where the section becomes
    effective at ``char_offset`` (just after the closing rule line). The title is
    the joined comment lines between the two rules, minus boilerplate. Only tight
    rule pairs (<= 6 lines apart) count, so per-entry prose comments and the
    ``% --- arXiv:... ---`` markers between distant rules are ignored.
    """
    lines = text.splitlines(keepends=True)
    offsets, off = [], 0
    for ln in lines:
        offsets.append(off)
        off += len(ln)
    rules = [i for i, ln in enumerate(lines) if _RULE.match(ln)]
    out: list[tuple[int, str]] = []
    for a, b in itertools.pairwise(rules):
        if b - a > 6:
            continue
        title_lines = []
        for k in range(a + 1, b):
            s = lines[k].strip().lstrip("%").strip()
            if not s or _BOILERPLATE.search(s):
                continue
            title_lines.append(s)
        if not title_lines:
            continue
        title = " ".join(title_lines).strip().rstrip(",")
        out.append((offsets[b] + len(lines[b]), title))
    return out


def _section_for(sections: list[tuple[int, str]], pos: int) -> str:
    sec = ""
    for o, title in sections:
        if o <= pos:
            sec = title
        else:
            break
    return sec


@dataclass
class Entry:
    """One parsed BibTeX entry. ``raw`` is the verbatim source span."""

    type: str
    key: str
    fields: dict[str, str] = dc_field(default_factory=dict)
    raw: str = ""
    start_line: int = 0
    section: str = ""


def _read_delimited(text: str, i: int) -> tuple[str, int]:
    """Read a field value starting at index ``i`` (first non-space after ``=``).

    Returns ``(value, next_index)``. Supports ``{...}`` (brace-depth balanced),
    ``"..."`` (may contain braces), and bare tokens up to ``,`` / ``}`` / newline.
    """
    n = len(text)
    ch = text[i]
    if ch == "{":
        depth = 0
        start = i
        while i < n:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    return text[start + 1 : i], i + 1
            i += 1
        return text[start + 1 :], n  # unbalanced; return rest
    if ch == '"':
        depth = 0
        start = i
        i += 1
        while i < n:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            elif c == '"' and depth == 0:
                return text[start + 1 : i], i + 1
            i += 1
        return text[start + 1 :], n
    # bare token
    start = i
    while i < n and text[i] not in ",}\n":
        i += 1
    return text[start:i].strip(), i


def _parse_entry_body(body: str) -> dict[str, str]:
    """Parse ``name = value`` pairs from an entry body (text after ``key,``)."""
    fields: dict[str, str] = {}
    i, n = 0, len(body)
    while i < n:
        m = re.compile(r"([A-Za-z][A-Za-z0-9_-]*)\s*=\s*").match(body, i)
        if not m:
            i += 1
            continue
        name = m.group(1).lower()
        j = m.end()
        while j < n and body[j] in " \t\r\n":
            j += 1
        if j >= n:
            break
        value, j = _read_delimited(body, j)
        fields[name] = value.strip()
        # skip trailing separator
        while j < n and body[j] in " \t\r\n,":
            j += 1
        i = j
    return fields


def parse_bibtex(text: str, track_sections: bool = True) -> list[Entry]:
    """Parse all entries from a BibTeX string.

    Works for both ``references.bib`` and a single fetched INSPIRE/Crossref entry.
    ``@string``/``@preamble``/``@comment`` are skipped.
    """
    entries: list[Entry] = []
    sections = _find_sections(text) if track_sections else []
    line_starts = [0]
    for mline in re.finditer(r"\n", text):
        line_starts.append(mline.end())

    def line_of(pos: int) -> int:
        lo, hi = 0, len(line_starts) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if line_starts[mid] <= pos:
                lo = mid
            else:
                hi = mid - 1
        return lo + 1

    pos = 0
    n = len(text)
    while pos < n:
        m = _ENTRY_START.search(text, pos)
        if not m:
            break
        section = _section_for(sections, m.start()) if track_sections else ""
        etype = m.group("type").lower()
        if etype in {"string", "preamble", "comment"}:
            pos = m.end()
            continue
        key = m.group("key")
        # find matching closing brace for the whole entry
        depth = 0
        i = text.index("{", m.start())
        start = m.start()
        while i < n:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        end = i + 1
        raw = text[start:end]
        body = text[m.end() : i]  # between "key," and final "}"
        ent = Entry(
            type=etype,
            key=key,
            fields=_parse_entry_body(body),
            raw=raw,
            start_line=line_of(start),
            section=section if track_sections else "",
        )
        entries.append(ent)
        pos = end
    return entries


def load_bib(path: str | Path) -> list[Entry]:
    return parse_bibtex(Path(path).read_text(encoding="utf-8"), track_sections=True)


# --- self-check / extraction validation ---------------------------------------

# Fields we require for a normal @article to consider extraction "clean".
_ARTICLE_REQUIRED = {"author", "title", "year"}


def _validate(path: str) -> int:
    entries = load_bib(path)
    print(f"parsed {len(entries)} entries from {path}")
    by_type: dict[str, int] = {}
    for e in entries:
        by_type[e.type] = by_type.get(e.type, 0) + 1
    print("by type:", dict(sorted(by_type.items())))

    problems = []
    seen_keys: dict[str, int] = {}
    for e in entries:
        seen_keys[e.key] = seen_keys.get(e.key, 0) + 1
        if not e.fields:
            problems.append((e.key, "NO FIELDS extracted"))
            continue
        if e.type == "article":
            missing = _ARTICLE_REQUIRED - set(e.fields)
            if missing:
                problems.append((e.key, f"missing {sorted(missing)}"))
        # sanity: a value that still contains an unbalanced brace count
        for fname, fval in e.fields.items():
            if fval.count("{") != fval.count("}"):
                problems.append((e.key, f"unbalanced braces in {fname!r}"))

    dupes = {k: c for k, c in seen_keys.items() if c > 1}
    if dupes:
        print(f"DUPLICATE keys: {dupes}")

    if problems:
        print(f"\n{len(problems)} extraction problem(s):")
        for key, msg in problems:
            print(f"  - {key}: {msg}")
        return 1
    print(
        "\nOK: all entries extracted cleanly (required fields present, braces balanced)."
    )
    return 0


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "manuscript/references.bib"
    raise SystemExit(_validate(target))
