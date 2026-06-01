r"""Task 2 worklist builder: every ``\\cite`` call-site, grouped by key (READ-ONLY).

Scans ``manuscript/sections/**/*.tex`` for ``\\cite{...}`` calls, captures the
enclosing sentence as the "claim" the citation backs, and groups by bib key.
Cross-references the offline ``literature/<eprint>/`` cache and
``literature_content_notes.md`` so the layered content check can start cheap.

Writes only ``cache/cites.json``. Does not touch the manuscript.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from scripts.bibaudit.bibio import load_bib

ROOT = Path("manuscript")
SECTIONS = ROOT / "sections"
LIT = Path("literature")
NOTES = ROOT / "planning" / "literature_content_notes.md"
CACHE = Path(__file__).with_name("cache")

_CITE = re.compile(r"\\cite\{([^}]*)\}")


def _sentence_around(text: str, start: int, end: int) -> str:
    """Return the sentence containing the cite span [start, end)."""
    left = max(
        text.rfind(". ", 0, start),
        text.rfind("\n\n", 0, start),
        text.rfind("}\n", 0, start),
    )
    left = left + 1 if left > 0 else 0
    right_dot = text.find(". ", end)
    right_par = text.find("\n\n", end)
    cands = [x for x in (right_dot, right_par) if x != -1]
    right = min(cands) + 1 if cands else len(text)
    snippet = text[left:right].strip()
    snippet = re.sub(r"\s+", " ", snippet)
    return snippet[:400]


def extract() -> dict:
    keys_eprint = {
        e.key: e.fields.get("eprint", "").strip()
        for e in load_bib("manuscript/references.bib")
    }
    cached_eprints = {p.name for p in LIT.iterdir()} if LIT.exists() else set()
    notes_text = NOTES.read_text() if NOTES.exists() else ""

    sites: list[dict] = []
    for tex in sorted(SECTIONS.rglob("*.tex")):
        text = tex.read_text(encoding="utf-8")
        # strip full-line comments so commented-out cites are not counted
        for m in _CITE.finditer(text):
            # ignore if the cite is on a commented line
            line_start = text.rfind("\n", 0, m.start()) + 1
            if text[line_start : m.start()].lstrip().startswith("%"):
                continue
            keys = [k.strip() for k in m.group(1).split(",") if k.strip()]
            sentence = _sentence_around(text, m.start(), m.end())
            line_no = text.count("\n", 0, m.start()) + 1
            rel = tex.relative_to(ROOT).as_posix()
            for k in keys:
                sites.append(
                    {"key": k, "file": rel, "line": line_no, "claim": sentence}
                )

    by_key: dict[str, dict] = {}
    for s in sites:
        rec = by_key.setdefault(
            s["key"],
            {
                "key": s["key"],
                "sites": [],
                "eprint": "",
                "offline": False,
                "in_notes": False,
            },
        )
        rec["sites"].append({"file": s["file"], "line": s["line"], "claim": s["claim"]})
    for k, rec in by_key.items():
        ep = keys_eprint.get(k, "")
        rec["eprint"] = ep
        rec["offline"] = bool(ep) and ep in cached_eprints
        rec["in_notes"] = (bool(ep) and ep in notes_text) or k in notes_text

    out = {
        "n_sites": len(sites),
        "n_keys": len(by_key),
        "keys": sorted(by_key.values(), key=lambda r: -len(r["sites"])),
    }
    CACHE.mkdir(exist_ok=True)
    (CACHE / "cites.json").write_text(json.dumps(out, indent=2))
    return out


def main() -> int:
    out = extract()
    print(f"call-sites: {out['n_sites']}   distinct keys: {out['n_keys']}")
    offline = sum(1 for r in out["keys"] if r["offline"])
    in_notes = sum(1 for r in out["keys"] if r["in_notes"])
    print(f"keys with offline literature/ TeX: {offline}")
    print(f"keys mentioned in literature_content_notes.md: {in_notes}")
    print("\n=== top-cited keys ===")
    for r in out["keys"][:15]:
        flags = ("O" if r["offline"] else "-") + ("N" if r["in_notes"] else "-")
        print(f"  {len(r['sites']):3d}x [{flags}] {r['key']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
