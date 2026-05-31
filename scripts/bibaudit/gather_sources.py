"""Stage B of the deep-completion: gather the BEST available source for every flagged
citation (MISMATCH/WEAK/PARTIAL/UNVERIFIABLE), so the Stage-C re-verify agents read
locally (no concurrent network). Writes per-key ``cache/flagged_src/<key>.json`` with the
consolidated source + local full-text paths, and downloads open-access PDFs into
``literature/<id>/oa.txt`` (pdftotext).

Sources, in order of richness: local full TeX (literature/) > open-access PDF text >
INSPIRE-by-DOI abstract (HEP) > Semantic Scholar abstract+TLDR > Google Books preview.
READ-ONLY w.r.t. the manuscript; writes only to cache/ and the gitignored literature/.
"""

from __future__ import annotations

import json
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

import requests

from scripts.bibaudit.bibio import load_bib
from scripts.bibaudit.compare import cite_final_status

CACHE = Path(__file__).with_name("cache")
LIT = Path("literature")
UA = "tidal-bibaudit/1.0 (mailto:wr286@cam.ac.uk)"
FLAGGED = {"MISMATCH", "WEAK", "PARTIAL", "UNVERIFIABLE"}


def _get_json(url, *, data=None, timeout=40):
    headers = {"User-Agent": UA}
    if data is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(data).encode()
    try:
        r = (
            requests.post(url, data=data, headers=headers, timeout=timeout)
            if data
            else requests.get(url, headers=headers, timeout=timeout)
        )
        if r.status_code == 200:
            return r.json()
    except Exception:
        return None
    return None


def s2_batch(ids):
    return (
        _get_json(
            "https://api.semanticscholar.org/graph/v1/paper/batch?fields=title,abstract,tldr,openAccessPdf",
            data={"ids": ids},
        )
        or []
    )


def inspire_abstract(doi):
    j = _get_json(
        "https://inspirehep.net/api/literature?q=doi:"
        + urllib.parse.quote(doi)
        + "&fields=titles,abstracts"
    )
    hits = (j or {}).get("hits", {}).get("hits", [])
    if hits:
        return (hits[0]["metadata"].get("abstracts") or [{}])[0].get("value", "")
    return ""


def googlebooks_preview(isbn):
    j = _get_json(
        "https://www.googleapis.com/books/v1/volumes?q=isbn:" + urllib.parse.quote(isbn)
    )
    items = (j or {}).get("items", [])
    if items:
        vi = items[0].get("volumeInfo", {})
        return (vi.get("description", "") or "")[:1500]
    return ""


def download_oapdf(url, safe):
    dest = LIT / safe
    txt = dest / "oa.txt"
    if txt.exists() and txt.stat().st_size > 500:
        return str(txt)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        raw = urllib.request.urlopen(req, timeout=60).read()
    except Exception:
        return ""
    if not raw[:5].startswith(b"%PDF"):
        return ""  # HTML landing page, not a PDF
    dest.mkdir(parents=True, exist_ok=True)
    pdf = dest / "oa.pdf"
    pdf.write_bytes(raw)
    try:
        subprocess.run(
            ["pdftotext", "-q", str(pdf), str(txt)], timeout=120, check=False
        )
    except Exception:
        return ""
    return str(txt) if txt.exists() and txt.stat().st_size > 500 else ""


def main() -> int:
    v = json.loads((CACHE / "cite_verify.json").read_text())["results"]
    bib = {e.key: e for e in load_bib("manuscript/references.bib")}
    flagged = [r for r in v if cite_final_status(r) in FLAGGED]
    out_dir = CACHE / "flagged_src"
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Semantic Scholar batch (one call)
    ids, keymap = [], {}
    for r in flagged:
        e = bib.get(r["key"])
        if not e:
            continue
        ep = e.fields.get("eprint", "").strip()
        doi = e.fields.get("doi", "").strip()
        sid = f"ARXIV:{ep}" if ep else (f"DOI:{doi}" if doi else None)
        if sid:
            ids.append(sid)
            keymap[sid] = r["key"]
    s2 = dict(zip(ids, s2_batch(ids), strict=False))

    for i, r in enumerate(flagged, 1):
        key = r["key"]
        e = bib.get(key)
        if not e:
            continue
        safe = key.replace("/", "_").replace(":", "_")
        ep = e.fields.get("eprint", "").strip()
        doi = e.fields.get("doi", "").strip()
        isbn = e.fields.get("isbn", "").strip()
        sid = f"ARXIV:{ep}" if ep else (f"DOI:{doi}" if doi else None)
        x = s2.get(sid) or {}
        rec = {
            "key": key,
            "final_status": cite_final_status(r),
            "first_status": r.get("status"),
            "sites": [
                {k2: s.get(k2) for k2 in ("file", "line", "role", "claim")}
                for s in (r.get("sites") or [])
            ],
            "abstract": (x.get("abstract") or ""),
            "tldr": (x.get("tldr") or {}).get("text") if x.get("tldr") else "",
            "oapdf_url": (x.get("openAccessPdf") or {}).get("url", "")
            if x.get("openAccessPdf")
            else "",
            "fulltext_paths": [],
        }
        # local full TeX (arXiv source fetched earlier, possibly)
        for d in {ep, ep.replace("/", "_")}:
            if d:
                p = LIT / d
                if p.exists():
                    rec["fulltext_paths"] += [str(q) for q in p.glob("*.tex")]
        # INSPIRE abstract if no S2 abstract and HEP DOI
        if not rec["abstract"] and doi:
            rec["abstract"] = inspire_abstract(doi)
            time.sleep(2.0)
        # Google Books preview for books
        if isbn and not rec["abstract"] and not rec["fulltext_paths"]:
            rec["gbooks"] = googlebooks_preview(isbn)
            time.sleep(1.5)
        # download OA-PDF -> text (richest if no local TeX)
        if rec["oapdf_url"] and not rec["fulltext_paths"]:
            txt = download_oapdf(rec["oapdf_url"], safe)
            if txt:
                rec["fulltext_paths"].append(txt)
            time.sleep(1.0)
        (out_dir / f"{safe}.json").write_text(json.dumps(rec, indent=2))
        has = (
            "TeX/PDF"
            if rec["fulltext_paths"]
            else (
                "abs"
                if rec["abstract"]
                else (
                    "tldr"
                    if rec["tldr"]
                    else ("gbooks" if rec.get("gbooks") else "NONE")
                )
            )
        )
        print(f"[{i:2d}/{len(flagged)}] {key:30s} {rec['final_status']:12s} src={has}")
    print(f"\nwrote {len(flagged)} flagged-source records to cache/flagged_src/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
