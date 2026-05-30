"""Task 1 driver: classify every bib entry against its canonical record.

READ-ONLY w.r.t. the manuscript. Produces:
  1. a human report on stdout, grouped by verdict;
  2. ``cache/report.json`` with full per-entry detail (verdict, source URL,
     canonical fields, substantive/minor/enrich diffs) and *proposed* exact
     ``(old_line -> new_line)`` edits + the ``%@audit`` marker line.

Nothing here edits ``references.bib``. Application of the proposed edits happens
separately, one at a time, through the Edit tool (exact string replacement that
fails loudly on any mismatch). compare.py only suggests; Edit verifies.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

from scripts.bibaudit import normalize as nrm
from scripts.bibaudit.bibio import Entry, load_bib
from scripts.bibaudit.fetch_canonical import fetch_for_entry

CACHE = Path(__file__).with_name("cache")
TODAY = "2026-05-30"


def _field_line(raw: str, field: str) -> tuple[str, str, str] | None:
    """Find the verbatim source line for ``field`` in an entry.

    Returns (full_line, indent, trailing) where trailing is ',' or ''. None if
    the field is not present on its own line.
    """
    pat = re.compile(
        rf"^([ \t]*){re.escape(field)}([ \t]*=[ \t]*)(\{{.*\}}|\".*\"|[^,\n]+?)([ \t]*,?)[ \t]*$",
        re.MULTILINE,
    )
    m = pat.search(raw)
    if not m:
        return None
    return m.group(0), m.group(1), m.group(4).strip()


def _replace_value_line(raw: str, field: str, new_value: str) -> tuple[str, str] | None:
    """Build (old_line, new_line) replacing ``field``'s value with ``new_value``.

    Preserves indentation, the ``= `` spacing, brace delimiters and trailing comma.
    """
    info = _field_line(raw, field)
    if not info:
        return None
    old_line, indent, trailing = info
    m = re.match(rf"^[ \t]*{re.escape(field)}([ \t]*=[ \t]*)", old_line)
    sep = m.group(1) if m else " = "
    comma = "," if trailing == "," else ""
    new_line = f"{indent}{field}{sep}{{{new_value}}}{comma}"
    return old_line, new_line


def _marker(canon: dict, verdict: str) -> str:
    # NB: BibTeX does NOT treat % as a comment and scans for '@' to start entries,
    # so marker lines must contain NO '@' (no '%@audit' prefix, no mailto e-mail).
    src = canon.get("url", "manual").replace("?mailto=wr286@cam.ac.uk", "")
    return f"% AUDIT v=1 date={TODAY} src={src} verdict={verdict}"


def build_edits(entry: Entry, canon: dict, cls: dict) -> dict:
    """Compute proposed exact edits for one entry (suggestions; Edit verifies)."""
    edits: list[dict] = []
    cfields = canon.get("fields", {})

    # journal abbreviation (MINOR 'abbreviate') and pages dash etc.
    for note in cls["minor"]:
        if note.startswith("journal:") and "(abbreviate)" in note:
            # target = the INSPIRE-abbreviated token for the LOCAL journal name
            abbrev = nrm.norm_journal(entry.fields.get("journal", ""))
            rep = _replace_value_line(entry.raw, "journal", abbrev) if abbrev else None
            if rep:
                edits.append({"kind": "journal", "old": rep[0], "new": rep[1]})
        elif note.startswith("pages:") and cfields.get("pages"):
            rep = _replace_value_line(
                entry.raw, "pages", nrm.norm_pages(cfields["pages"])
            )
            if rep:
                edits.append({"kind": "pages", "old": rep[0], "new": rep[1]})

    # enrichment: fields canonical has that local lacks
    for f, val in cls["enrich"].items():
        edits.append({"kind": f"enrich:{f}", "field": f, "value": val})

    return {
        "marker": _marker(canon, cls["verdict"]),
        "edits": edits,
    }


def run(refresh: bool = False) -> dict:
    entries = load_bib("manuscript/references.bib")
    report = []
    for i, e in enumerate(entries, 1):
        canon = fetch_for_entry(e, refresh=refresh)
        if not canon.get("found"):
            verdict = canon.get("reason", "NOTFOUND")
            cls = {
                "verdict": verdict,
                "substantive": [],
                "minor": [],
                "enrich": {},
                "journal_needs_map": False,
            }
        else:
            cls = nrm.classify(e.fields, canon.get("fields", {}))
        rec = {
            "key": e.key,
            "type": e.type,
            "section": e.section,
            "line": e.start_line,
            "verdict": cls["verdict"],
            "service": canon.get("service"),
            "via": canon.get("via"),
            "src": canon.get("url", "manual"),
            "texkey": canon.get("texkey"),
            "substantive": cls["substantive"],
            "minor": cls["minor"],
            "enrich": cls["enrich"],
            "journal_needs_map": cls["journal_needs_map"],
            "proposed": build_edits(e, canon, cls)
            if canon.get("found")
            else {"marker": _marker(canon, cls["verdict"]), "edits": []},
        }
        report.append(rec)
        print(f"[{i:3d}/{len(entries)}] {e.key:42s} {cls['verdict']}")
    CACHE.mkdir(exist_ok=True)
    (CACHE / "report.json").write_text(json.dumps(report, indent=2))
    return {"report": report}


def summarize(report: list[dict]) -> None:
    from collections import Counter

    counts = Counter(r["verdict"] for r in report)
    print("\n=== verdict summary ===")
    for v in nrm.SEVERITY:
        if counts.get(v):
            print(f"  {v:12s} {counts[v]}")
    print(f"  {'TOTAL':12s} {len(report)}")

    subs = [r for r in report if r["verdict"] == "SUBSTANTIVE"]
    if subs:
        print(f"\n=== SUBSTANTIVE / STOP review ({len(subs)}) ===")
        for r in subs:
            print(f"  - {r['key']} ({r['src']})")
            for s in r["substantive"]:
                print(f"      {s}")

    unmapped = [r for r in report if r["journal_needs_map"]]
    if unmapped:
        print(f"\n=== unmapped journals ({len(unmapped)}) ===")
        for r in unmapped:
            for note in r["minor"]:
                if "UNMAPPED" in note:
                    print(f"  - {r['key']}: {note}")

    nf = [r for r in report if r["verdict"] in {"NOTFOUND", "NOID"}]
    if nf:
        print(f"\n=== NOTFOUND / NOID ({len(nf)}) ===")
        for r in nf:
            print(f"  - {r['key']} [{r['verdict']}] ({r['type']})")

    pending = [r for r in report if r["verdict"] == "PENDING"]
    if pending:
        print(
            f"\n=== PENDING — transient fetch failure, NOT cached ({len(pending)}) ==="
        )
        for r in pending:
            print(f"  - {r['key']}")
        print("  >> re-run `python -m scripts.bibaudit.compare` to retry these.")


def _has_marker_above(text: str, idx: int) -> bool:
    """True if a %@audit marker is in the contiguous comment block above ``idx``."""
    for ln in reversed(text[:idx].rstrip("\n").splitlines()):
        s = ln.lstrip()
        if s.startswith("%"):
            if s.startswith("%@audit"):
                return True
        elif s == "":
            continue
        else:
            break
    return False


def _apply_entry(raw: str, rec: dict) -> tuple[str | None, str | None]:
    """Return (new_entry_text, error). Applies line edits + enrichment to ``raw``."""
    new = raw
    for ed in rec["proposed"]["edits"]:
        if ed.get("old") is not None and ed.get("new") is not None:
            if new.count(ed["old"]) < 1:
                return None, f"line not found: {ed['kind']}"
            new = new.replace(ed["old"], ed["new"], 1)
    enrich = [ed for ed in rec["proposed"]["edits"] if ed["kind"].startswith("enrich:")]
    if enrich:
        nl = new.find("\n")  # insert right after the "@type{key," opening line
        if nl < 0:
            return None, "no newline after opening"
        ins = "".join(f"  {ed['field']:<13} = {{{ed['value']}}},\n" for ed in enrich)
        new = new[: nl + 1] + ins + new[nl + 1 :]
    return new, None


def apply_all(report: list[dict], bib_path: str, *, dry_run: bool) -> dict:
    """Apply markers + journal/pages edits + enrichment via exact block replacement.

    Each entry's verbatim block is located (must be unique), edits applied to it,
    the %@audit marker prepended, and the block swapped in. Entries that already
    carry a %@audit marker are skipped (idempotent). Aborts a single entry on any
    mismatch rather than guessing.
    """
    text = Path(bib_path).read_text(encoding="utf-8")
    entries = {e.key: e for e in load_bib(bib_path)}
    applied, skipped, failed = 0, 0, []
    for rec in report:
        e = entries.get(rec["key"])
        if not e:
            continue
        raw = e.raw
        n = text.count(raw)
        if n != 1:
            failed.append((rec["key"], f"block not unique (count={n})"))
            continue
        if _has_marker_above(text, text.find(raw)):
            skipped += 1
            continue
        new_entry, err = _apply_entry(raw, rec)
        if err:
            failed.append((rec["key"], err))
            continue
        block = rec["proposed"]["marker"] + "\n" + new_entry
        text = text.replace(raw, block, 1)
        applied += 1
    if not dry_run:
        Path(bib_path).write_text(text, encoding="utf-8")
    print(
        f"\n{'DRY-RUN ' if dry_run else ''}apply: {applied} applied, {skipped} skipped "
        f"(already marked), {len(failed)} failed"
    )
    for key, why in failed:
        print(f"  FAIL {key}: {why}")
    return {"applied": applied, "skipped": skipped, "failed": failed}


def main(argv: list[str]) -> int:
    refresh = "--refresh" in argv
    if "--summary-only" in argv and (CACHE / "report.json").exists():
        report = json.loads((CACHE / "report.json").read_text())
    else:
        report = run(refresh=refresh)["report"]
    if "--apply" in argv or "--dry-run" in argv:
        apply_all(report, "manuscript/references.bib", dry_run=("--apply" not in argv))
    else:
        summarize(report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
