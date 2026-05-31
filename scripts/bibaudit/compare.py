"""Task 1 driver: classify every bib entry against its canonical record.

READ-ONLY w.r.t. the manuscript. Produces:
  1. a human report on stdout, grouped by verdict;
  2. ``cache/report.json`` with full per-entry detail (verdict, source URL,
     canonical fields, substantive/minor/enrich diffs) and *proposed* exact
     ``(old_line -> new_line)`` edits + the ``% AUDIT`` marker line.

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
    # so marker lines must contain NO '@' (no '% AUDIT' prefix, no mailto e-mail).
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
    """True if a ``% AUDIT`` marker is in the contiguous comment block above ``idx``."""
    for ln in reversed(text[:idx].rstrip("\n").splitlines()):
        s = ln.lstrip()
        if s.startswith("%"):
            if s.replace(" ", "").upper().startswith("%AUDIT"):
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
    the % AUDIT marker prepended, and the block swapped in. Entries that already
    carry a % AUDIT marker are skipped (idempotent). Aborts a single entry on any
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


def _sanitize(s: str) -> str:
    """Strip BibTeX-hazardous '@' and newlines from comment text."""
    return (s or "").replace("@", " at ").replace("\n", " ").replace("  ", " ").strip()


def cite_final_status(rec: dict) -> str:
    """The recheck-aware final verdict (adversarial recheck overrides the first pass)."""
    rc = rec.get("recheck") or {}
    return rc.get("status") or rec.get("status", "?")


def _cite_block(rec: dict) -> str:
    """Build a BibTeX-safe ``% CITE`` block (no '@') from a Task-2 verdict.

    Uses the recheck-aware FINAL status; notes the first-pass verdict when the
    adversarial recheck changed it (``was=...``). Prefers the recheck's evidence/
    depth/source when present.
    """
    rc = rec.get("recheck") or {}
    final = cite_final_status(rec)
    first = rec.get("status", "?")
    sites = rc.get("sites") or rec.get("sites") or []
    used = ", ".join(
        f"{(s.get('file', '') or '').split('/')[-1]} L{s.get('line', '')} {s.get('role', '')}"
        for s in sites[:8]
    )
    depth = rc.get("depth_reached") or rec.get("depth_reached", "")
    src = rc.get("source_checked") or rec.get("source_checked", "")
    conf = rc.get("confidence") or rec.get("confidence", "")
    evidence = rc.get("evidence_quote") or rec.get("evidence_quote", "")
    head = f"% CITE v=1 status={final} depth={_sanitize(depth)} conf={conf}"
    if rc and final != first:
        head += f" was={first}"
    lines = [head, f"%! purpose: {_sanitize(rec.get('purpose', ''))[:200]}"]
    if used:
        lines.append(f"%! used: {_sanitize(used)[:300]}")
    if src:
        lines.append(f"%! checked: {_sanitize(src)[:160]} [{TODAY}]")
    if final != "SUPPORTS" and evidence:
        lines.append(f"%! evidence: {_sanitize(evidence)[:220]}")
    return "\n".join(lines) + "\n"


def _strip_cite_block(text: str, key: str) -> str:
    """Remove an existing ``% CITE`` block (the % CITE line + its %! lines, which sit
    above the % AUDIT line) for ``key``, leaving % AUDIT and the entry intact.
    """
    pat = re.compile(
        r"% CITE[^\n]*\n(?:%![^\n]*\n)*(?=% AUDIT[^\n]*\n(?:%![^\n]*\n)*@\w+\{"
        + re.escape(key)
        + r",)"
    )
    return pat.sub("", text, count=1)


def apply_cites(
    verify: list[dict], bib_path: str, *, dry_run: bool, force: bool = False
) -> dict:
    """Prepend a ``% CITE`` block above each entry's contiguous comment block (above
    the existing ``% AUDIT`` line). With ``force``, an existing ``% CITE`` block is
    stripped and rewritten (refresh); otherwise such entries are skipped. Exact-span
    replacement; aborts an entry on any non-uniqueness.
    """
    text = Path(bib_path).read_text(encoding="utf-8")
    entries = {e.key: e for e in load_bib(bib_path)}
    applied, skipped, refreshed, failed = 0, 0, 0, []
    for rec in verify:
        key = rec.get("key")
        e = entries.get(key)
        if not e:
            failed.append((key, "key not in bib"))
            continue
        raw = e.raw
        if text.count(raw) != 1:
            failed.append((key, f"block not unique (count={text.count(raw)})"))
            continue
        idx = text.find(raw)
        pre = text[:idx].rstrip("\n").split("\n") if idx > 0 else []
        nmark = 0
        for ln in reversed(pre):
            if ln.lstrip().startswith("%"):
                nmark += 1
            else:
                break
        block_lines = pre[len(pre) - nmark :] if nmark else []
        if any("% CITE" in ln for ln in block_lines):
            if not force:
                skipped += 1
                continue
            text = _strip_cite_block(text, key)
            refreshed += 1
            # recompute the (now CITE-free) comment span
            idx = text.find(raw)
            pre = text[:idx].rstrip("\n").split("\n") if idx > 0 else []
            nmark = 0
            for ln in reversed(pre):
                if ln.lstrip().startswith("%"):
                    nmark += 1
                else:
                    break
            block_lines = pre[len(pre) - nmark :] if nmark else []
        old = ("\n".join(block_lines) + "\n" if nmark else "") + raw
        if text.count(old) != 1:
            failed.append((key, "span not unique"))
            continue
        text = text.replace(old, _cite_block(rec) + old, 1)
        applied += 1
    if not dry_run:
        Path(bib_path).write_text(text, encoding="utf-8")
    print(
        f"\n{'DRY-RUN ' if dry_run else ''}cites-apply: {applied} applied "
        f"({refreshed} refreshed), {skipped} skipped, {len(failed)} failed"
    )
    for key, why in failed:
        print(f"  FAIL {key}: {why}")
    return {
        "applied": applied,
        "skipped": skipped,
        "refreshed": refreshed,
        "failed": failed,
    }


def cites_summary(verify: list[dict]) -> None:
    """Print the Task-2 STOP list on the recheck-FINAL verdicts: MISMATCH / WEAK /
    PARTIAL / UNVERIFIABLE groups.
    """
    from collections import Counter

    tally = Counter(cite_final_status(r) for r in verify)
    print("\n=== citation-support tally (recheck-final) ===")
    for s in ["SUPPORTS", "PARTIAL", "WEAK", "MISMATCH", "UNVERIFIABLE"]:
        if tally.get(s):
            print(f"  {s:13s} {tally[s]}")
    for grp in ["MISMATCH", "WEAK", "PARTIAL", "UNVERIFIABLE"]:
        rows = [r for r in verify if cite_final_status(r) == grp]
        if rows:
            print(f"\n=== {grp} ({len(rows)}) ===")
            for r in rows:
                first = r.get("status")
                tag = f" [was {first}]" if first != grp else ""
                print(f"  - {r.get('key')}{tag}: {(r.get('purpose') or '')[:70]}")


def main(argv: list[str]) -> int:
    refresh = "--refresh" in argv
    if (
        "--cites-apply" in argv
        or "--cites-dry-run" in argv
        or "--cites-summary" in argv
    ):
        data = json.loads((CACHE / "cite_verify.json").read_text())
        recs = data.get("results", data) if isinstance(data, dict) else data
        if "--cites-summary" in argv:
            cites_summary(recs)
        else:
            apply_cites(
                recs,
                "manuscript/references.bib",
                dry_run=("--cites-apply" not in argv),
                force=("--force" in argv),
            )
        return 0
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
