r"""Regenerate the ``literature/README.md`` inventory table from what is on disk.

``literature/`` holds arXiv TeX sources kept locally so papers can be read without
repeated web fetches (see ``.claude/rules/literature.md``). Its README used to be
hand-maintained and drifted badly -- it documented 15 of 124 directories and claimed
there were 20. This module rebuilds the table from the directories themselves:

  * the directory name gives the arXiv id (``astro-ph_9603033`` -> ``astro-ph/9603033``);
  * the main TeX is the one file containing ``\begin{document}``;
  * title and authors come from the arXiv API, batched and cached.

Only the region between the ``<!-- BEGIN generated index -->`` / ``<!-- END ... -->``
markers is rewritten; the surrounding prose is left alone. ``docs/references.md`` remains
the canonical, curated index -- this file is the on-disk inventory of what is downloaded.

Usage::

    uv run python -m scripts.bibaudit.index_literature            # rewrite the table
    uv run python -m scripts.bibaudit.index_literature --check    # exit 1 if stale
    uv run python -m scripts.bibaudit.index_literature --refresh  # ignore the cache
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from html import unescape
from pathlib import Path

CACHE = Path(__file__).with_name("cache") / "arxiv_meta"
USER_AGENT = "tidal-bibaudit/1.0 (mailto:wr286@cam.ac.uk)"
API = "https://export.arxiv.org/api/query"

BEGIN_MARKER = "<!-- BEGIN generated index -->"
END_MARKER = "<!-- END generated index -->"

#: arXiv accepts long ``id_list`` batches; 60 keeps each response comfortably small.
BATCH = 60
#: Courtesy pause between live API calls.
THROTTLE_S = 3.0

_NEW_STYLE = re.compile(r"^\d{4}\.\d{4,5}$")
_OLD_STYLE = re.compile(r"^(?P<archive>[a-z-]+(?:\.[A-Z]{2})?)_(?P<num>\d{7})$")

#: Directories holding a paper that predates arXiv or is stored as a publisher PDF.
#: Resolved once via Crossref and pinned here so the generator stays a single-service
#: tool; keyed by directory name.
MANUAL: dict[str, tuple[str, str, str]] = {
    "BF02721794": (
        "Hamiltonian structure of the theory of gravity with R+T^2 type of Lagrangian",
        "Blagojevic & Nikolic (1983)",
        "https://doi.org/10.1007/BF02721794",
    ),
    "PhysRep258.1": (
        "Metric-affine gauge theory of gravity: field equations, Noether identities, "
        "world spinors, and breaking of dilation invariance",
        "Hehl, McCrea, Mielke & Ne'eman (1995)",
        "https://doi.org/10.1016/0370-1573(94)00111-F",
    ),
    "PhysRevD.28.2455": (
        "Hamiltonian dynamics of Poincare gauge theory: General structure in the time gauge",
        "Blagojevic & Nikolic (1983)",
        "https://doi.org/10.1103/PhysRevD.28.2455",
    ),
}


@dataclass(frozen=True)
class Row:
    """One inventory line: a directory, what it holds, and where it came from."""

    directory: str
    arxiv_id: str | None
    link: str | None
    main_tex: str | None
    title: str
    authors: str


# --- directory -> arXiv id ----------------------------------------------------


def dir_to_arxiv_id(name: str) -> str | None:
    """Map a ``literature/`` directory name back to its arXiv id.

    Old-style ids are stored with ``/`` replaced by ``_`` (the convention in
    ``fetch_fulltext.py``), so ``gr-qc_0305049`` -> ``gr-qc/0305049``. Returns
    ``None`` for directories that do not name an arXiv paper.
    """
    if _NEW_STYLE.match(name):
        return name
    m = _OLD_STYLE.match(name)
    if m:
        return f"{m['archive']}/{m['num']}"
    return None


# --- main TeX detection -------------------------------------------------------


def find_main_tex(directory: Path) -> str | None:
    """Return the file name of the directory's main TeX, or ``None`` if it has none.

    The main file is the one that opens the document body. When a source ships
    several such files (rare -- usually a stray copy of a template), prefer
    ``main.tex`` and otherwise the largest, so the choice is deterministic.
    """
    tex = sorted(directory.glob("*.tex"))
    if not tex:
        return None
    bodies = [
        p
        for p in tex
        if r"\begin{document}" in p.read_text(encoding="utf-8", errors="replace")
    ]
    if not bodies:
        bodies = tex
    for p in bodies:
        if p.name == "main.tex":
            return p.name
    return max(bodies, key=lambda p: p.stat().st_size).name


# --- arXiv metadata -----------------------------------------------------------


def _parse_entries(xml: str) -> dict[str, dict[str, str]]:
    """Pull ``{id: {title, authors}}`` out of an Atom response from the arXiv API."""
    out: dict[str, dict[str, str]] = {}
    for entry in re.findall(r"<entry>(.*?)</entry>", xml, re.DOTALL):
        m_id = re.search(r"<id>.*?abs/(.*?)(?:v\d+)?</id>", entry, re.DOTALL)
        m_title = re.search(r"<title>(.*?)</title>", entry, re.DOTALL)
        if not m_id or not m_title:
            continue
        names = [
            unescape(n).strip()
            for n in re.findall(r"<name>(.*?)</name>", entry, re.DOTALL)
        ]
        out[m_id.group(1)] = {
            "title": " ".join(unescape(m_title.group(1)).split()),
            "authors": _format_authors(names),
        }
    return out


def _format_authors(names: list[str]) -> str:
    """Condense an author list to the surname form used in the table.

    Collaboration entries ("LIGO Scientific Collaboration") are kept whole -- their
    last word is not a surname.
    """
    surnames = [
        n if "Collaboration" in n else (n.split()[-1] if n.split() else n)
        for n in names
    ]
    if not surnames:
        return "-"
    if len(surnames) == 1:
        return surnames[0]
    if len(surnames) == 2:
        return f"{surnames[0]} & {surnames[1]}"
    return f"{surnames[0]} et al."


def _cache_path(arxiv_id: str) -> Path:
    return CACHE / f"{arxiv_id.replace('/', '_')}.json"


def fetch_metadata(
    ids: list[str], *, refresh: bool = False
) -> dict[str, dict[str, str]]:
    """Resolve titles and authors for ``ids``, caching one JSON file per paper.

    Cached ids are never re-requested, so a second run needs no network at all.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    meta: dict[str, dict[str, str]] = {}
    missing: list[str] = []
    for aid in ids:
        path = _cache_path(aid)
        if path.exists() and not refresh:
            meta[aid] = json.loads(path.read_text(encoding="utf-8"))
        else:
            missing.append(aid)

    for start in range(0, len(missing), BATCH):
        batch = missing[start : start + BATCH]
        query = urllib.parse.urlencode(
            {"id_list": ",".join(batch), "max_results": len(batch)}
        )
        req = urllib.request.Request(
            f"{API}?{query}", headers={"User-Agent": USER_AGENT}
        )
        print(f"  arXiv API: {len(batch)} ids...", file=sys.stderr)
        with urllib.request.urlopen(req, timeout=90) as resp:
            found = _parse_entries(resp.read().decode("utf-8", errors="replace"))
        for aid in batch:
            record = found.get(aid)
            if record is None:
                print(f"  WARNING: arXiv returned nothing for {aid}", file=sys.stderr)
                continue
            _cache_path(aid).write_text(json.dumps(record, indent=2), encoding="utf-8")
            meta[aid] = record
        if start + BATCH < len(missing):
            time.sleep(THROTTLE_S)
    return meta


# --- rows and rendering -------------------------------------------------------


def build_rows(literature: Path, *, refresh: bool = False) -> list[Row]:
    """Inventory every subdirectory of ``literature`` as a table row."""
    dirs = sorted(p for p in literature.iterdir() if p.is_dir())
    ids = [aid for p in dirs if (aid := dir_to_arxiv_id(p.name))]
    meta = fetch_metadata(ids, refresh=refresh)

    rows: list[Row] = []
    for path in dirs:
        aid = dir_to_arxiv_id(path.name)
        if aid:
            record = meta.get(aid, {})
            title = record.get("title", "(not resolved on arXiv)")
            authors = record.get("authors", "-")
            link = f"https://arxiv.org/abs/{aid}"
        elif path.name in MANUAL:
            title, authors, link = MANUAL[path.name]
        else:
            title, authors, link = "(unidentified)", "-", None
        rows.append(
            Row(
                directory=path.name,
                arxiv_id=aid,
                link=link,
                main_tex=find_main_tex(path),
                title=title,
                authors=authors,
            )
        )
    return rows


def _escape(text: str) -> str:
    """Make a cell safe to drop into a markdown table."""
    return text.replace("|", r"\|")


def render_table(rows: list[Row]) -> str:
    """Render the inventory rows as the README's markdown table."""
    lines = [
        "| Directory | Source | Main TeX | Title | Authors |",
        "| --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        if r.arxiv_id and r.link:
            source = f"[{r.arxiv_id}]({r.link})"
        elif r.link:
            source = f"[{r.directory}]({r.link})"
        else:
            source = "--"
        tex = f"`{r.main_tex}`" if r.main_tex else "-- (PDF only)"
        lines.append(
            f"| `{r.directory}/` | {source} | {tex} | {_escape(r.title)} | {_escape(r.authors)} |"
        )
    return "\n".join(lines)


def splice(readme: str, table: str) -> str:
    """Replace the generated region of the README, leaving the prose untouched."""
    if BEGIN_MARKER not in readme or END_MARKER not in readme:
        msg = (
            f"README is missing the {BEGIN_MARKER} / {END_MARKER} markers; "
            "add them around the table before running the generator."
        )
        raise SystemExit(msg)
    head, rest = readme.split(BEGIN_MARKER, 1)
    _, tail = rest.split(END_MARKER, 1)
    return f"{head}{BEGIN_MARKER}\n\n{table}\n\n{END_MARKER}{tail}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--literature-dir",
        type=Path,
        default=Path("literature"),
        help="directory holding the local sources (default: ./literature)",
    )
    parser.add_argument(
        "--readme",
        type=Path,
        default=None,
        help="README to rewrite (default: <dir>/README.md)",
    )
    parser.add_argument(
        "--check", action="store_true", help="exit 1 if the README is out of date"
    )
    parser.add_argument(
        "--refresh", action="store_true", help="re-fetch cached metadata"
    )
    args = parser.parse_args(argv)

    literature: Path = args.literature_dir
    readme: Path = args.readme or literature / "README.md"
    if not literature.is_dir():
        print(f"error: no such directory: {literature}", file=sys.stderr)
        return 2
    if not readme.is_file():
        print(f"error: no such file: {readme}", file=sys.stderr)
        return 2

    rows = build_rows(literature, refresh=args.refresh)
    updated = splice(readme.read_text(encoding="utf-8"), render_table(rows))

    if args.check:
        if updated != readme.read_text(encoding="utf-8"):
            print(f"{readme} is out of date; re-run without --check", file=sys.stderr)
            return 1
        print(f"{readme} is up to date ({len(rows)} entries)")
        return 0

    readme.write_text(updated, encoding="utf-8")
    no_tex = sum(1 for r in rows if r.main_tex is None)
    print(f"wrote {readme}: {len(rows)} entries ({no_tex} without TeX)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
