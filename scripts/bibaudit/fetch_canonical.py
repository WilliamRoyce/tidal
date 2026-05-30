"""Canonical-record fetcher for the citation audit (READ-ONLY w.r.t. manuscript).

Pulls the authoritative record for each bib entry from, in order of preference:
  1. INSPIRE-HEP  (arXiv eprint, then DOI)   -- physics literature, abbreviated journals
  2. Crossref     (DOI)                       -- non-HEP journals / books with a DOI
  3. Google Books (ISBN)                      -- books

Every raw response is cached under ``scripts/bibaudit/cache/<service>/`` so the run
is idempotent and resumable: re-runs read the cache and never re-hit the network
unless ``--refresh`` is passed. Writes ONLY to the cache directory.

Network etiquette: 1 s throttle between live calls, exponential backoff on 429/503,
and a ``User-Agent`` carrying the project's academic contact email so Crossref
routes us through its faster "polite pool".
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
from pathlib import Path

import requests

from scripts.bibaudit.bibio import load_bib, parse_bibtex

CACHE = Path(__file__).with_name("cache")
USER_AGENT = "tidal-bibaudit/1.0 (mailto:wr286@cam.ac.uk)"
# INSPIRE tolerates ~6s spacing (faster bursts trigger a 429 penalty); Crossref's
# polite pool is far more generous. Override per-run via BIBAUDIT_THROTTLE_S.
_THROTTLE_S = float(os.environ.get("BIBAUDIT_THROTTLE_S", "6"))
_last_call = [0.0]


def _throttle() -> None:
    dt = time.monotonic() - _last_call[0]
    if dt < _THROTTLE_S:
        time.sleep(_THROTTLE_S - dt)
    _last_call[0] = time.monotonic()


def _get(
    url: str,
    *,
    accept: str | None = None,
    tries: int | None = None,
    wait429: float | None = None,
) -> requests.Response | None:
    """GET with configurable rate-limit handling.

    On 429/503 we wait ``wait429`` seconds (honoring Retry-After if larger) and
    retry up to ``tries`` times. Defaults come from the environment so a whole run
    can be made fail-fast: ``BIBAUDIT_TRIES`` (default 6) and ``BIBAUDIT_WAIT429``
    (default 60). With ``BIBAUDIT_TRIES=1 BIBAUDIT_WAIT429=0`` every rate-limited
    entry gives up immediately (-> PENDING, never cached) so the run never stalls;
    re-running then mops up the transients (cached successes are skipped).
    """
    if tries is None:
        tries = int(os.environ.get("BIBAUDIT_TRIES", "6"))
    if wait429 is None:
        wait429 = float(os.environ.get("BIBAUDIT_WAIT429", "60"))
    headers = {"User-Agent": USER_AGENT}
    if accept:
        headers["Accept"] = accept
    for attempt in range(tries):
        _throttle()
        try:
            r = requests.get(url, headers=headers, timeout=30)
        except requests.RequestException as exc:
            print(
                f"    ! request error ({exc}); retry {attempt + 1}/{tries}",
                file=sys.stderr,
            )
            if wait429 <= 0:
                return None
            time.sleep(20.0)
            continue
        if r.status_code in {429, 503}:
            if wait429 <= 0:
                return None  # fail-fast: let the caller mark PENDING
            ra = r.headers.get("Retry-After")
            wait = max(float(ra) if (ra and ra.isdigit()) else 0.0, wait429)
            print(
                f"    ! {r.status_code}; wait {wait}s (attempt {attempt + 1}/{tries})",
                file=sys.stderr,
            )
            time.sleep(wait)
            continue
        return r
    return None


def _cache_path(service: str, ident: str) -> Path:
    safe = ident.replace("/", "_").replace(":", "_")
    return CACHE / service / f"{safe}.json"


def _load_cache(service: str, ident: str) -> dict | None:
    p = _cache_path(service, ident)
    if p.exists():
        return json.loads(p.read_text())
    return None


def _save_cache(service: str, ident: str, payload: dict) -> None:
    p = _cache_path(service, ident)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=2))


# --- per-service fetchers ------------------------------------------------------


def inspire(query: str, ident: str, *, refresh: bool = False) -> dict | None:
    """Fetch one INSPIRE record as bibtex. ``query`` is the ``q=`` value."""
    cached = None if refresh else _load_cache("inspire", ident)
    if cached is not None:
        return cached
    url = (
        "https://inspirehep.net/api/literature?q="
        + urllib.parse.quote(query)
        + "&format=bibtex"
    )
    r = _get(url)  # rate-limit policy from BIBAUDIT_TRIES / BIBAUDIT_WAIT429
    if r is None or r.status_code != 200:
        # transient (network / rate-limit exhaustion): do NOT cache, so a re-run retries
        return {
            "found": False,
            "transient": True,
            "url": url,
            "status": getattr(r, "status_code", None),
        }
    if not r.text.strip().startswith("@"):
        # genuine: INSPIRE has no bibtex record for this query
        payload = {"found": False, "url": url, "status": 200}
        _save_cache("inspire", ident, payload)
        return payload
    entries = parse_bibtex(r.text, track_sections=False)
    payload = {
        "found": bool(entries),
        "url": url,
        "bibtex": r.text.strip(),
        "fields": entries[0].fields if entries else {},
        "texkey": entries[0].key if entries else None,
    }
    _save_cache("inspire", ident, payload)
    return payload


def crossref(doi: str, *, refresh: bool = False) -> dict | None:
    cached = None if refresh else _load_cache("crossref", doi)
    if cached is not None:
        return cached
    # `mailto` query param is Crossref's documented "polite pool" signal -> faster,
    # more reliable than the anonymous pool (which 429s under shared load).
    url = (
        "https://api.crossref.org/works/"
        + urllib.parse.quote(doi, safe="")
        + "?mailto=wr286@cam.ac.uk"
    )
    r = _get(url, accept="application/json")
    if r is None or r.status_code not in {200, 404}:
        return {
            "found": False,
            "transient": True,
            "url": url,
            "status": getattr(r, "status_code", None),
        }
    if r.status_code == 404:
        payload = {"found": False, "url": url, "status": 404}
        _save_cache("crossref", doi, payload)
        return payload
    msg = r.json().get("message", {})
    authors = " and ".join(
        f"{a.get('family', '')}, {a.get('given', '')}".strip(", ")
        for a in msg.get("author", [])
    )
    fields = {
        "title": (msg.get("title") or [""])[0],
        "author": authors,
        "year": str((msg.get("issued", {}).get("date-parts") or [[None]])[0][0] or ""),
        "journal": (msg.get("container-title") or [""])[0],
        "volume": msg.get("volume", ""),
        "pages": msg.get("page", ""),
        "doi": msg.get("DOI", ""),
    }
    payload = {
        "found": True,
        "url": url,
        "fields": {k: v for k, v in fields.items() if v},
    }
    _save_cache("crossref", doi, payload)
    return payload


def googlebooks(isbn: str, *, refresh: bool = False) -> dict | None:
    cached = None if refresh else _load_cache("googlebooks", isbn)
    if cached is not None:
        return cached
    url = "https://www.googleapis.com/books/v1/volumes?q=isbn:" + urllib.parse.quote(
        isbn
    )
    r = _get(url, accept="application/json")
    if r is None or r.status_code != 200:
        return {
            "found": False,
            "transient": True,
            "url": url,
            "status": getattr(r, "status_code", None),
        }
    if not r.json().get("items"):
        payload = {"found": False, "url": url, "status": 200}
        _save_cache("googlebooks", isbn, payload)
        return payload
    vi = r.json()["items"][0].get("volumeInfo", {})
    fields = {
        "title": vi.get("title", "")
        + (f": {vi['subtitle']}" if vi.get("subtitle") else ""),
        "author": " and ".join(vi.get("authors", [])),
        "year": (vi.get("publishedDate", "")[:4]),
        "publisher": vi.get("publisher", ""),
    }
    payload = {
        "found": True,
        "url": url,
        "fields": {k: v for k, v in fields.items() if v},
    }
    _save_cache("googlebooks", isbn, payload)
    return payload


# --- book resolution (title+author search; no usable id on the entry) ----------


def crossref_search(title: str, author: str, *, rows: int = 3) -> list[dict]:
    """Crossref bibliographic search restricted to book-like types -> candidates."""
    q = urllib.parse.quote(f"{title} {author}".strip())
    url = (
        f"https://api.crossref.org/works?query.bibliographic={q}"
        "&filter=type:book,type:monograph,type:reference-book,type:edited-book"
        f"&rows={rows}&mailto=wr286@cam.ac.uk"
    )
    r = _get(url, accept="application/json")
    if r is None or r.status_code != 200:
        return []
    out = []
    for m in r.json().get("message", {}).get("items", []):
        out.append(
            {
                "doi": m.get("DOI", ""),
                "title": (m.get("title") or [""])[0],
                "author": ", ".join(a.get("family", "") for a in m.get("author", [])),
                "year": str(
                    (m.get("issued", {}).get("date-parts") or [[None]])[0][0] or ""
                ),
                "publisher": m.get("publisher", ""),
                "isbn": (m.get("ISBN") or [""])[0],
                "type": m.get("type", ""),
            }
        )
    return out


def googlebooks_search(title: str, author: str) -> list[dict]:
    """Google Books title+author search -> candidates (for ISBN / publisher / year)."""
    q = urllib.parse.quote(f"intitle:{title} inauthor:{author}")
    url = f"https://www.googleapis.com/books/v1/volumes?q={q}&maxResults=3"
    r = _get(url, accept="application/json")
    if r is None or r.status_code != 200:
        return []
    out = []
    for it in r.json().get("items", []):
        vi = it.get("volumeInfo", {})
        isbns = {
            x.get("type"): x.get("identifier")
            for x in vi.get("industryIdentifiers", [])
        }
        out.append(
            {
                "title": vi.get("title", "")
                + (f": {vi['subtitle']}" if vi.get("subtitle") else ""),
                "author": ", ".join(vi.get("authors", [])),
                "year": (vi.get("publishedDate", "") or "")[:4],
                "publisher": vi.get("publisher", ""),
                "isbn13": isbns.get("ISBN_13", ""),
                "isbn10": isbns.get("ISBN_10", ""),
            }
        )
    return out


def resolve_books(*, refresh: bool = False) -> dict:
    """For every book-like entry, gather Crossref + Google Books candidates by
    title+author so a human/agent can verify + enrich. Writes cache/books/<key>.json
    and prints a per-book candidate report. Does NOT edit the bib.
    """
    entries = load_bib("manuscript/references.bib")
    booklike = {"book", "inbook", "incollection"}
    out = {}
    for e in entries:
        if e.type not in booklike:
            continue
        title = e.fields.get("title", "")
        author = e.fields.get("author", e.fields.get("editor", ""))
        lead = author.split(" and ")[0].split(",")[0].strip()
        cr = crossref_search(title[:80], lead)
        gb = googlebooks_search(title[:60], lead)
        rec = {
            "key": e.key,
            "local": {
                k: e.fields.get(k, "")
                for k in (
                    "title",
                    "author",
                    "editor",
                    "year",
                    "publisher",
                    "isbn",
                    "doi",
                )
            },
            "crossref": cr,
            "googlebooks": gb,
        }
        out[e.key] = rec
        cdoi = cr[0]["doi"] if cr else "-"
        gisbn = (gb[0]["isbn13"] or gb[0]["isbn10"]) if gb else "-"
        print(
            f"\n=== {e.key} ({e.type}) — local doi={e.fields.get('doi', '-')} isbn={e.fields.get('isbn', '-')}"
        )
        print(
            f"    crossref[0]: doi={cdoi} | {cr[0]['title'][:60] if cr else '-'} | {cr[0]['author'][:40] if cr else ''} ({cr[0]['year'] if cr else ''})"
        )
        print(
            f"    gbooks[0]:   isbn={gisbn} | {gb[0]['title'][:60] if gb else '-'} | {gb[0]['author'][:40] if gb else ''} ({gb[0]['year'] if gb else ''})"
        )
    (CACHE / "books").mkdir(parents=True, exist_ok=True)
    (CACHE / "books" / "candidates.json").write_text(json.dumps(out, indent=2))
    print(f"\n{len(out)} book-like entries; candidates -> cache/books/candidates.json")
    return out


# --- dispatch ------------------------------------------------------------------


def fetch_for_entry(entry, *, refresh: bool = False) -> dict:
    """Pick the best identifier and fetch. Returns the canonical payload + source.

    Source order is chosen for reliability + coverage:
      1. reuse an already-cached INSPIRE-by-arXiv record (don't waste prior fetches);
      2. DOI -> Crossref (generous polite pool, authoritative for published data,
         and covers the non-HEP journals INSPIRE does not index);
      3. eprint -> INSPIRE (the only source for arXiv-only entries);
      4. ISBN -> Google Books.
    INSPIRE rate-limits aggressively, so we only hit it live for arXiv-only entries.
    If every attempt failed only transiently, reason='PENDING' (not cached) so a
    re-run retries it.
    """
    f = entry.fields
    ep = f.get("eprint", "").strip()
    doi = f.get("doi", "").strip()
    isbn = f.get("isbn", "").strip().replace("-", "")
    transient = False

    # 1. reuse cached INSPIRE-by-arXiv (authoritative + already fetched)
    if ep and not refresh:
        cached = _load_cache("inspire", f"arxiv_{ep}")
        if cached and cached.get("found"):
            return {**cached, "service": "inspire", "via": "arxiv"}

    # 2. DOI -> Crossref
    if doi:
        p = crossref(doi, refresh=refresh)
        if p and p.get("found"):
            return {**p, "service": "crossref", "via": "doi"}
        transient = transient or bool(p and p.get("transient"))

    # 3. eprint -> INSPIRE (arXiv-only entries; also retried for DOI'd ones not on Crossref)
    if ep:
        p = inspire(f"arxiv:{ep}", f"arxiv_{ep}", refresh=refresh)
        if p and p.get("found"):
            return {**p, "service": "inspire", "via": "arxiv"}
        transient = transient or bool(p and p.get("transient"))

    # 4. ISBN -> Google Books
    if isbn:
        p = googlebooks(isbn, refresh=refresh)
        if p and p.get("found"):
            return {**p, "service": "googlebooks", "via": "isbn"}
        transient = transient or bool(p and p.get("transient"))

    has_id = any((ep, doi, isbn))
    reason = "PENDING" if transient else ("NOTFOUND" if has_id else "NOID")
    return {"found": False, "service": None, "reason": reason}


def main(argv: list[str]) -> int:
    refresh = "--refresh" in argv
    if "--mode" in argv and argv[argv.index("--mode") + 1] == "books":
        resolve_books(refresh=refresh)
        return 0
    bib = "manuscript/references.bib"
    entries = load_bib(bib)
    found = notfound = noid = 0
    for i, e in enumerate(entries, 1):
        p = fetch_for_entry(e, refresh=refresh)
        if p.get("found"):
            found += 1
            tag = f"{p['service']}/{p['via']}"
        elif p.get("reason") == "NOID":
            noid += 1
            tag = "NOID"
        else:
            notfound += 1
            tag = "NOTFOUND"
        print(f"[{i:3d}/{len(entries)}] {e.key:42s} {tag}")
    print(f"\nfound={found} notfound={notfound} noid={noid} total={len(entries)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
