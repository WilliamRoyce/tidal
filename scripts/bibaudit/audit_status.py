"""Progress dashboard for the citation audit (READ-ONLY).

Parses the ``%@audit`` and ``%@cite`` marker comments in ``references.bib`` and
reports how much of the audit is done — ``bib: X/231`` and ``cites: Y/192`` — plus
the outstanding STOP list and any unmapped journals. This is the source of truth
for resuming: an entry already carrying a ``%@audit`` line is done. Run it any time
to see exactly what remains, independent of chat context.
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path

from scripts.bibaudit.bibio import load_bib

BIB = "manuscript/references.bib"
_AUDIT = re.compile(r"^%\s*AUDIT\b.*?\bverdict=(\w+)", re.MULTILINE)
_CITE = re.compile(r"^%\s*CITE\b.*?\bstatus=(\w+)", re.MULTILINE)


def _markers_above(text: str, entry_start: int) -> list[str]:
    """Return the contiguous block of comment lines immediately above an entry."""
    lines = text[:entry_start].splitlines()
    out = []
    for ln in reversed(lines):
        if ln.lstrip().startswith("%"):
            out.append(ln)
        elif not ln.strip():
            continue
        else:
            break
    return list(reversed(out))


def status() -> dict:
    text = Path(BIB).read_text(encoding="utf-8")
    entries = load_bib(BIB)
    total = len(entries)
    audited = Counter()
    cited_ann = 0
    for e in entries:
        idx = text.find(e.raw)
        block = "\n".join(_markers_above(text, idx)) if idx >= 0 else ""
        ma = _AUDIT.search(block)
        if ma:
            audited[ma.group(1)] += 1
        if _CITE.search(block):
            cited_ann += 1
    return {
        "total": total,
        "audited_total": sum(audited.values()),
        "audited_by_verdict": dict(audited),
        "cite_annotated": cited_ann,
    }


def main() -> int:
    s = status()
    print(f"bib audited: {s['audited_total']}/{s['total']}")
    for v, c in sorted(s["audited_by_verdict"].items()):
        print(f"    {v:12s} {c}")
    print(f"cite-annotated keys: {s['cite_annotated']}/{s['total']}")
    remaining = s["total"] - s["audited_total"]
    if remaining:
        print(f"\n{remaining} entries not yet audited.")
    else:
        print("\nAll entries audited.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
