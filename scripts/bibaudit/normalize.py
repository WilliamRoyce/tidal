r"""Field normalization + entry classification for the citation audit.

Pure functions, no I/O beyond loading ``journal_map.json``. Used by ``compare.py``
to decide, for each entry, whether the local fields match the canonical INSPIRE /
Crossref record and — if not — whether the difference is cosmetic (MINOR),
fillable (ENRICHED), or a real attribution problem (SUBSTANTIVE).

The normalization is the false-positive defense: it must treat "Physical Review D"
== "Phys. Rev. D", "Poincar\\'e" == "Poincaré" == "Poincare", "84--85" == "84-85",
and "{Exact}" == "Exact" as equal, while still catching a swapped author or a
wrong title/year (the documented mis-attribution class).
"""

from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

_MAP_PATH = Path(__file__).with_name("journal_map.json")
_JOURNAL_MAP: dict[str, str] = json.loads(_MAP_PATH.read_text())["map"]

# --- generic text helpers ------------------------------------------------------

# common LaTeX accent commands -> plain letter is handled by stripping the command
# wrapper and then NFKD-folding; e.g. {\'e} / \'{e} / \'e -> e
_ACCENT_CMD = re.compile(r"\\[`'^\"~=.uvHcdbra]+\s*\{?([A-Za-z])\}?")
_BRACES = re.compile(r"[{}]")
_WS = re.compile(r"\s+")
_TEX_DASH = re.compile(r"\s*-{2,3}\s*")  # -- or ---
_LATEX_CMD = re.compile(r"\\[a-zA-Z@]+")
_NONPRINT = re.compile(r"[\ufffd\x00-\x1f]")  # mojibake replacement char / control

# text-mode LaTeX commands that Crossref/INSPIRE emit; map to a plain equivalent
# BEFORE accent-stripping so they are not mangled into stray letters.
_TEXT_CMDS = {
    r"\textendash": "-",
    r"\textemdash": "-",
    r"\textquoteright": "'",
    r"\textquoteleft": "'",
    r"\textquotedblright": '"',
    r"\textquotedblleft": '"',
    r"\textbackslash": " ",
    r"\&": "&",
}


# curly quotes, thin/no-break spaces, unicode dashes -> plain ASCII equivalents
# (Crossref emits these; they otherwise create spurious author/title diffs)
_PUNCT_FOLD = str.maketrans(
    {
        "’": "'",
        "‘": "'",
        "“": '"',
        "”": '"',
        " ": " ",
        " ": " ",
        "–": "-",
        "—": "-",
        "−": "-",
    }
)
_XML_TAG = re.compile(r"<[^>]+>")
_XML_ENT = re.compile(r"&[a-zA-Z]+;|&#\d+;")


def _strip_latex(s: str) -> str:
    s = s.translate(_PUNCT_FOLD)
    s = _XML_TAG.sub(" ", s)  # strip MathML/HTML tags Crossref leaves in titles
    s = _XML_ENT.sub(" ", s)
    for cmd, rep in _TEXT_CMDS.items():
        s = s.replace(cmd, rep)
    s = _ACCENT_CMD.sub(r"\1", s)  # \'e -> e
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = _LATEX_CMD.sub(" ", s)  # any remaining \commands
    return s.replace("\\", " ")


def _strip_accents(s: str) -> str:
    return _strip_latex(s)


def debrace(s: str) -> str:
    return _BRACES.sub("", s)


def collapse_ws(s: str) -> str:
    return _WS.sub(" ", s).strip()


def fold(s: str) -> str:
    """Aggressive fold for token comparison: LaTeX, accents, braces, case, ws."""
    return collapse_ws(_NONPRINT.sub(" ", _strip_latex(debrace(s)))).casefold()


# --- per-field normalizers -----------------------------------------------------


def norm_journal(s: str) -> str | None:
    """Return the canonical token for a journal name, or None if unmapped."""
    if not s:
        return None
    key = collapse_ws(debrace(s)).casefold().rstrip(".")
    # try exact, then with trailing period restored
    for variant in (key, key + "."):
        if variant in _JOURNAL_MAP:
            return _JOURNAL_MAP[variant]
    return None


def norm_title(s: str) -> str:
    """Token set of the title for content comparison (order-insensitive)."""
    folded = fold(s).rstrip(".")
    # drop punctuation that varies between sources
    folded = re.sub(r"[^\w\s]", " ", folded)
    return " ".join(sorted(t for t in folded.split() if t))


def norm_title_seq(s: str) -> str:
    """Title as an ordered string (accents/braces/case-folded). For display diff."""
    return collapse_ws(re.sub(r"[^\w\s]", " ", fold(s))).rstrip()


def _surname(author: str) -> str:
    """Extract a comparable surname token from one author entry."""
    a = _strip_accents(debrace(author)).strip()
    if "," in a:
        sur = a.split(",", 1)[0]
    else:
        parts = a.split()
        sur = parts[-1] if parts else a
    sur = collapse_ws(sur).casefold()
    for ch in "-. '":
        sur = sur.replace(ch, "")
    return sur


def split_authors(s: str) -> list[str]:
    # bibtex author lists are separated by " and "
    return [a.strip() for a in re.split(r"\s+and\s+", s.strip()) if a.strip()]


def author_surnames(s: str) -> list[str]:
    out = []
    for a in split_authors(s):
        if a.lower() in {"others", "et al.", "et al", "collaboration"}:
            out.append("others")
        else:
            out.append(_surname(a))
    return out


def norm_pages(s: str) -> str:
    s = s.strip()
    s = _TEX_DASH.sub("-", s)
    s = re.sub(r"\s*[–—]\s*", "-", s)  # en/em dash -> hyphen
    return s.strip()


def norm_arxiv(s: str) -> str:
    s = s.strip().replace("arXiv:", "").replace("arxiv:", "")
    s = re.sub(r"v\d+$", "", s)  # strip version suffix
    return s.casefold()


def norm_doi(s: str) -> str:
    s = s.strip()
    s = re.sub(r"^https?://(dx\.)?doi\.org/", "", s, flags=re.IGNORECASE)
    s = re.sub(r"^doi:\s*", "", s, flags=re.IGNORECASE)
    return s.casefold()


def norm_year(s: str) -> str:
    m = re.search(r"\d{4}", s)
    return m.group(0) if m else s.strip()


# --- classification ------------------------------------------------------------

# verdict severity ordering for aggregation
SEVERITY = ["MATCH", "MINOR", "ENRICHED", "NOID", "NOTFOUND", "PENDING", "SUBSTANTIVE"]


def classify(local: dict[str, str], canon: dict[str, str]) -> dict:
    """Compare a local entry's fields against a canonical record.

    Returns a dict with: verdict, substantive (list of field-level mismatch
    descriptions), minor (list), enrich (dict of fields canonical has that local
    lacks/abbreviates and could be filled), journal_needs_map (bool).
    ``canon`` is expected to already be parsed into a field dict.
    """
    substantive: list[str] = []
    minor: list[str] = []
    enrich: dict[str, str] = {}
    journal_needs_map = False

    # --- author (mis-attribution detection) ---
    # The reliable mis-attribution signal is the LEAD author. Crossref routinely
    # reorders compound surnames, truncates author lists, and prefixes a stray
    # " and ", so positional/length comparison false-positives heavily. Rule:
    #   * lead surnames match (or one side's lead is contained in the other) AND
    #     the named-surname sets overlap -> same paper (MINOR if not identical);
    #   * otherwise the entry likely points at the wrong work -> SUBSTANTIVE.
    if local.get("author") and canon.get("author"):
        lnames = [x for x in author_surnames(local["author"]) if x and x != "others"]
        cnames = [x for x in author_surnames(canon["author"]) if x and x != "others"]
        if lnames and cnames:
            lead_l, lead_c = lnames[0], cnames[0]
            lead_ok = lead_l == lead_c or lead_l in lead_c or lead_c in lead_l
            overlap = len(set(lnames) & set(cnames)) / min(len(lnames), len(cnames))
            same_paper = lead_ok or overlap >= 0.5
            if not same_paper:
                substantive.append(
                    f"author: local={local['author']!r} canonical={canon['author']!r}"
                )
            elif set(lnames) != set(cnames):
                minor.append(
                    f"author: list differs (order/truncation) "
                    f"local={local['author']!r} canonical={canon['author']!r}"
                )

    # --- title (token-overlap, tolerant of truncation / encoding noise) ---
    if local.get("title") and canon.get("title"):
        lt = {t for t in norm_title(local["title"]).split() if len(t) > 1}
        ct = {t for t in norm_title(canon["title"]).split() if len(t) > 1}
        if lt and ct:
            inter = lt & ct
            jacc = len(inter) / len(lt | ct)
            subset = lt <= ct or ct <= lt
            if lt == ct:
                if norm_title_seq(local["title"]) != norm_title_seq(canon["title"]):
                    minor.append("title: brace/case/punctuation differs")
            elif subset or jacc >= 0.6:
                # truncated source (e.g. Crossref 'SUNDIALS'), 'Part I' suffix, or
                # an encoding-mangled token: cosmetic, not a mis-attribution
                minor.append(
                    f"title: minor diff (jaccard={jacc:.2f}) "
                    f"local={local['title']!r} canonical={canon['title']!r}"
                )
            else:
                substantive.append(
                    f"title: local={local['title']!r} canonical={canon['title']!r}"
                )

    # --- year (preprint vs publication year differs by <=1 routinely -> MINOR) ---
    if local.get("year") and canon.get("year"):
        ly, cy = norm_year(local["year"]), norm_year(canon["year"])
        if ly != cy:
            try:
                delta = abs(int(ly) - int(cy))
            except ValueError:
                delta = 99
            note = f"year: local={ly!r} canonical={cy!r}"
            (minor if delta <= 1 else substantive).append(note)

    # --- journal ---
    # Journal standardization is a LOCAL transform: map the local name to its
    # INSPIRE-abbreviated token via journal_map (so it works the same whether the
    # canonical came from INSPIRE, Crossref, or nothing). The canonical journal is
    # used only as a secondary identity cross-check.
    lj = local.get("journal", "")
    cj = canon.get("journal", "")
    if not lj and cj:
        enrich["journal"] = norm_journal(cj) or cj
    elif lj:
        lt = norm_journal(lj)
        if lt is None:
            journal_needs_map = True
            minor.append(f"journal: {lj!r} UNMAPPED (extend journal_map)")
        else:
            ct = norm_journal(cj) if cj else None
            if ct is not None and ct != lt:
                minor.append(
                    f"journal: local maps {lt!r} but canonical={cj!r} ({ct}) (check)"
                )
            if lj.strip() != lt:
                minor.append(f"journal: {lj!r} -> {lt!r} (abbreviate)")

    # --- volume / number / pages ---
    # NB: volume/number conventions legitimately differ across sources (e.g. JHEP /
    # JCAP list volume = year with the issue in `number`, while INSPIRE puts the
    # issue in `volume`). So these are MINOR notes for review, never SUBSTANTIVE —
    # the reliable mis-attribution signals are author, title and year.
    for f in ("volume", "number", "pages"):
        lv, cv = local.get(f, ""), canon.get(f, "")
        if not cv:
            continue
        if not lv:
            enrich[f] = cv
            continue
        if f == "pages":
            if norm_pages(lv) != norm_pages(cv):
                if norm_pages(lv).replace("-", "") == norm_pages(cv).replace("-", ""):
                    minor.append(f"pages: {lv!r} -> {cv!r} (dash)")
                else:
                    minor.append(f"pages: local={lv!r} canonical={cv!r}")
        elif lv.strip() != cv.strip():
            minor.append(f"{f}: local={lv!r} canonical={cv!r} (convention?)")

    # --- identifiers (enrich if missing) ---
    # A differing DOI is usually an alias (e.g. a JSTOR vs AMS DOI for one paper)
    # rather than a mis-attribution -> MINOR. Genuine wrong-DOI entries are caught
    # by the author/title checks above (the DOI resolves to a different work).
    for f, normf in (("doi", norm_doi), ("eprint", norm_arxiv)):
        lv, cv = local.get(f, ""), canon.get(f, "")
        if cv and not lv:
            enrich[f] = cv
        elif lv and cv and normf(lv) != normf(cv):
            minor.append(f"{f}: local={lv!r} canonical={cv!r} (differs)")
        elif lv and cv and lv.strip() != cv.strip():
            minor.append(f"{f}: {lv!r} -> {cv!r} (format)")

    # --- verdict ---
    if substantive:
        verdict = "SUBSTANTIVE"
    elif enrich:
        verdict = "ENRICHED"
    elif minor:
        verdict = "MINOR"
    else:
        verdict = "MATCH"

    return {
        "verdict": verdict,
        "substantive": substantive,
        "minor": minor,
        "enrich": enrich,
        "journal_needs_map": journal_needs_map,
    }
