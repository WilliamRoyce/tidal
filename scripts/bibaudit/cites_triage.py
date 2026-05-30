"""Task 2 Stage 1 — triage the citation worklist (READ-ONLY, no network).

Joins ``cache/cites.json`` (per-key call-sites + claim sentences, from
``extract_cites.py``) with each entry's identifiers in ``references.bib`` and the
offline ``literature/`` cache, then tags every cited key with:

  * a content-source **tier** (drives how deep verification can go, and whether the
    serial abstract pre-fetch in Stage 2 needs to fetch it):
      A  offline TeX in literature/<eprint>/        -> deepest local read
      B  arXiv eprint, not offline                  -> abstract pre-fetch
      C  DOI only (no eprint)                        -> Crossref/INSPIRE abstract
      D  ISBN only (book)                            -> metadata only
      E  no identifier                              -> bib fields only
  * a per-site **citation role** heuristic: EXISTENCE/CONTEXT vs ATTRIBUTION vs
    SPECIFIC (a formula/number/result is attributed) — SPECIFIC/ATTRIBUTION sites
    are the ones worth escalating to a deep read.

Writes ``cache/cites_worklist.json``. Pure read; the Stage-3 verification agents
consume this.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

from scripts.bibaudit.bibio import load_bib

CACHE = Path(__file__).with_name("cache")
LIT = Path("literature")

# claim-sentence cues for the per-site role heuristic
_ATTRIB = re.compile(
    r"\b(first|originally|introduced|derived|proposed|due to|following|"
    r"convention of|formul\w+ of|result of|theorem of|as shown (?:by|in))\b",
    re.IGNORECASE,
)
_SPECIFIC = re.compile(
    r"(=|\\sin|\\kappa|\\propto|probabilit|amplitud|\d+\s*%|"
    r"factor|equation|\\eqref|\\cref|formula|coefficient|cross[- ]section|"
    r"\\frac|order|eigenvalue|spectrum|bound|constraint)",
    re.IGNORECASE,
)
_EXISTENCE = re.compile(
    r"\b(see also|e\.g\.|cf\.|review|overview|among others|such as|"
    r"detect\w*|observ\w*|measur\w*|survey|catalog)\b",
    re.IGNORECASE,
)


def _role(claim: str) -> str:
    """Heuristic citation role for one call-site's claim sentence."""
    if _ATTRIB.search(claim):
        return "ATTRIBUTION"
    if _SPECIFIC.search(claim):
        return "SPECIFIC"
    return "EXISTENCE"  # default: context/landmark cite -> cheap pass suffices


def _tier(fields: dict, offline_ids: set[str]) -> tuple[str, str]:
    ep = fields.get("eprint", "").strip()
    doi = fields.get("doi", "").strip()
    isbn = fields.get("isbn", "").strip()
    if ep and ep in offline_ids:
        return "A", ep
    if ep:
        return "B", ep
    if doi:
        return "C", doi
    if isbn:
        return "D", isbn
    return "E", ""


def triage() -> dict:
    cites = json.loads((CACHE / "cites.json").read_text())
    by_key = {r["key"]: r for r in cites["keys"]}
    entries = {e.key: e for e in load_bib("manuscript/references.bib")}
    offline = {p.name for p in LIT.iterdir()} if LIT.exists() else set()

    work = []
    for key, rec in by_key.items():
        e = entries.get(key)
        fields = e.fields if e else {}
        tier, ident = _tier(fields, offline)
        sites = []
        roles = set()
        for s in rec["sites"]:
            role = _role(s.get("claim", ""))
            roles.add(role)
            sites.append({**s, "role": role})
        # the key's required depth = the most demanding role across its sites
        needs_deep = bool(roles & {"SPECIFIC", "ATTRIBUTION"})
        work.append(
            {
                "key": key,
                "tier": tier,
                "ident": ident,
                "eprint": fields.get("eprint", ""),
                "doi": fields.get("doi", ""),
                "n_sites": len(sites),
                "needs_deep": needs_deep,
                "in_notes": rec.get("in_notes", False),
                "sites": sites,
            }
        )

    # batches: tier first (A local-deep, B/C abstract, D/E shallow), then by site count
    work.sort(key=lambda r: (r["tier"], -r["n_sites"]))
    out = {
        "n_keys": len(work),
        "by_tier": {t: sum(1 for w in work if w["tier"] == t) for t in "ABCDE"},
        "needs_deep": sum(1 for w in work if w["needs_deep"]),
        "keys": work,
    }
    CACHE.mkdir(exist_ok=True)
    (CACHE / "cites_worklist.json").write_text(json.dumps(out, indent=2))
    return out


def main() -> int:
    out = triage()
    print(
        f"keys: {out['n_keys']} | needs-deep (SPECIFIC/ATTRIBUTION): {out['needs_deep']}"
    )
    print("by content-source tier:")
    labels = {
        "A": "offline TeX",
        "B": "arXiv (prefetch)",
        "C": "DOI (prefetch)",
        "D": "book/ISBN",
        "E": "no-id",
    }
    for t in "ABCDE":
        n = out["by_tier"][t]
        deep = sum(1 for w in out["keys"] if w["tier"] == t and w["needs_deep"])
        print(f"  {t} {labels[t]:18s}: {n:3d} keys ({deep} need deep read)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
