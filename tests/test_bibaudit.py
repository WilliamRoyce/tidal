"""Offline tests for the citation-audit tooling (scripts/bibaudit/).

No network: validates the BibTeX field extractor and the field-comparison
classifier against the real manuscript bibliography. These guard the read-only
analysis layer; all actual edits to references.bib go through the Edit tool.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from scripts.bibaudit import normalize as nrm
from scripts.bibaudit.bibio import load_bib, parse_bibtex

BIB = "manuscript/references.bib"

pytestmark = pytest.mark.skipif(
    not Path(BIB).exists(), reason="manuscript bibliography not present"
)


# --- extractor ----------------------------------------------------------------


def test_all_entries_extract_cleanly():
    entries = load_bib(BIB)
    assert len(entries) == 231
    for e in entries:
        assert e.fields, f"{e.key}: no fields extracted"
        for name, val in e.fields.items():
            assert val.count("{") == val.count("}"), f"{e.key}: unbalanced {name}"


def test_article_required_fields_present():
    for e in load_bib(BIB):
        if e.type == "article":
            assert {"author", "title", "year"} <= set(e.fields), e.key


def test_no_duplicate_keys():
    keys = [e.key for e in load_bib(BIB)]
    assert len(keys) == len(set(keys))


def test_sections_assigned():
    # every entry should fall under a detected % === section === banner
    entries = load_bib(BIB)
    assert all(e.section for e in entries)


def test_parses_both_delimiter_styles():
    # local uses {} ; INSPIRE exports use ""
    text = '@article{a, title = {Braced {Nested} Title}, author = "Doe, J."}'
    (e,) = parse_bibtex(text, track_sections=False)
    assert e.fields["title"] == "Braced {Nested} Title"
    assert e.fields["author"] == "Doe, J."


# --- normalizers --------------------------------------------------------------


def test_journal_abbreviation_map():
    assert nrm.norm_journal("Journal of High Energy Physics") == "JHEP"
    assert nrm.norm_journal("Physical Review D") == "Phys. Rev. D"
    assert nrm.norm_journal("Phys. Rev. D") == "Phys. Rev. D"
    assert nrm.norm_journal("Totally Made Up Journal") is None


def test_accent_and_brace_folding():
    assert nrm.fold("Poincar\\'e") == nrm.fold("Poincaré") == "poincare"
    assert nrm.norm_pages("84--85") == nrm.norm_pages("84-85") == "84-85"


def test_author_surnames():
    assert nrm.author_surnames("Ejlli, Damian") == ["ejlli"]
    assert nrm.author_surnames("Raffelt, G. and Stodolsky, L.") == [
        "raffelt",
        "stodolsky",
    ]


# --- classifier ---------------------------------------------------------------

_CANON = {
    "author": "Ejlli, Damian",
    "title": "Graviton-photon mixing. Exact solution in a constant magnetic field",
    "journal": "JHEP",
    "volume": "06",
    "pages": "029",
    "year": "2020",
    "doi": "10.1007/JHEP06(2020)029",
    "eprint": "2004.02714",
}


# Fixed local fixture (the *original* Ejlli entry, full journal name + Springer
# volume convention). Independent of the live bib, which the audit mutates.
_LOCAL = {
    "author": "Ejlli, Damian",
    "title": "Graviton--photon mixing. {Exact} solution in a constant magnetic field",
    "journal": "Journal of High Energy Physics",
    "volume": "2020",
    "number": "06",
    "pages": "029",
    "year": "2020",
    "doi": "10.1007/JHEP06(2020)029",
    "eprint": "2004.02714",
}


def test_classify_match_identical():
    assert nrm.classify(dict(_CANON), dict(_CANON))["verdict"] == "MATCH"


def test_classify_minor_journal_abbrev():
    c = nrm.classify(dict(_LOCAL), _CANON)
    assert c["verdict"] == "MINOR"
    assert any("abbreviate" in m for m in c["minor"])
    assert not c["substantive"]  # volume convention diff must NOT be substantive


def test_classify_substantive_on_swapped_author():
    c = nrm.classify(dict(_LOCAL), dict(_CANON, author="Domcke, Valerie"))
    assert c["verdict"] == "SUBSTANTIVE"
    assert any("author" in s for s in c["substantive"])


def test_classify_enriches_missing_field():
    stub = {k: v for k, v in _LOCAL.items() if k != "journal"}
    c = nrm.classify(stub, _CANON)
    assert "journal" in c["enrich"]
