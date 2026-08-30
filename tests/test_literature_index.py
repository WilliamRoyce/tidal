"""Offline tests for the ``literature/README.md`` inventory generator.

No network: the arXiv metadata stage is exercised against a seeded cache directory,
and ``urlopen`` is replaced by a fixture that fails if it is ever reached. Kept in its
own file (rather than ``test_bibaudit.py``) because these tests do not depend on the
manuscript bibliography, which that module skips on.
"""

from __future__ import annotations

import json
from typing import TYPE_CHECKING

import pytest

from scripts.bibaudit import index_literature as idx

if TYPE_CHECKING:
    from pathlib import Path


def _seed_cache(tmp_path: Path, records: dict[str, dict[str, str]]) -> Path:
    """Build a metadata cache directory pre-filled with ``{id: record}``."""
    cache = tmp_path / "arxiv_meta"
    cache.mkdir()
    for aid, rec in records.items():
        (cache / f"{aid.replace('/', '_')}.json").write_text(json.dumps(rec))
    return cache


# --- directory -> arXiv id ----------------------------------------------------


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("2606.30785", "2606.30785"),
        ("0803.1967", "0803.1967"),
        ("astro-ph_9603033", "astro-ph/9603033"),
        ("gr-qc_0305049", "gr-qc/0305049"),
        ("hep-th_0103093", "hep-th/0103093"),
        # Pre-arXiv / publisher-PDF directories are not arXiv ids.
        ("BF02721794", None),
        ("PhysRep258.1", None),
        ("PhysRevD.28.2455", None),
        ("README.md", None),
    ],
)
def test_dir_to_arxiv_id(name: str, expected: str | None) -> None:
    assert idx.dir_to_arxiv_id(name) == expected


# --- main TeX detection -------------------------------------------------------


def test_find_main_tex_picks_the_document_body(tmp_path: Path) -> None:
    (tmp_path / "macros.tex").write_text(r"\newcommand{\x}{y}")
    (tmp_path / "paper.tex").write_text(r"\begin{document}hello\end{document}")
    assert idx.find_main_tex(tmp_path) == "paper.tex"


def test_find_main_tex_prefers_main_over_template_leftovers(tmp_path: Path) -> None:
    # 2205.13962 ships four REVTeX sample files alongside the real manuscript, and
    # some are large -- the name must win over size.
    (tmp_path / "main.tex").write_text(r"\begin{document}real")
    (tmp_path / "apssamp.tex").write_text(r"\begin{document}" + "x" * 5000)
    assert idx.find_main_tex(tmp_path) == "main.tex"


def test_find_main_tex_falls_back_to_largest_body(tmp_path: Path) -> None:
    # No main.tex: the real manuscript is the biggest file with a document body.
    (tmp_path / "ms.tex").write_text(r"\begin{document}" + "real" * 2000)
    (tmp_path / "sorsamp.tex").write_text(r"\begin{document}template")
    assert idx.find_main_tex(tmp_path) == "ms.tex"


def test_find_main_tex_none_when_pdf_only(tmp_path: Path) -> None:
    (tmp_path / "main.pdf").write_bytes(b"%PDF-1.4")
    assert idx.find_main_tex(tmp_path) is None


# --- metadata formatting ------------------------------------------------------


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (["Eiichiro Komatsu"], "Komatsu"),
        (["Chung-Pei Ma", "Edmund Bertschinger"], "Ma & Bertschinger"),
        (["S. Blanes", "F. Casas", "J. A. Oteo", "J. Ros"], "Blanes et al."),
        # A collaboration's last word is not a surname.
        (["The LIGO Scientific Collaboration"], "The LIGO Scientific Collaboration"),
        ([], "-"),
    ],
)
def test_format_authors(names: list[str], expected: str) -> None:
    assert idx._format_authors(names) == expected


def test_parse_entries_reads_the_atom_response() -> None:
    xml = """
    <feed><entry>
      <id>http://arxiv.org/abs/2606.30785v2</id>
      <title>Numerical polology:
      towards next-generation model-building</title>
      <author><name>Will Barker</name></author>
      <author><name>Carlo Marzo</name></author>
    </entry></feed>
    """
    parsed = idx._parse_entries(xml)
    # Version suffix stripped, wrapped title reflowed onto one line.
    assert parsed["2606.30785"]["title"] == (
        "Numerical polology: towards next-generation model-building"
    )
    assert parsed["2606.30785"]["authors"] == "Barker & Marzo"


# --- rendering ----------------------------------------------------------------


def test_render_table_escapes_pipes() -> None:
    row = idx.Row(
        directory="1234.5678",
        arxiv_id="1234.5678",
        link="https://arxiv.org/abs/1234.5678",
        main_tex="main.tex",
        title="A | B",
        authors="Doe",
    )
    body = idx.render_table([row]).splitlines()[-1]
    assert r"A \| B" in body
    # The escaped pipe must not open a sixth column.
    assert body.count("|") - body.count(r"\|") == 6


def test_render_table_marks_pdf_only_entries() -> None:
    row = idx.Row(
        directory="PhysRep258.1",
        arxiv_id=None,
        link="https://doi.org/10.1016/0370-1573(94)00111-F",
        main_tex=None,
        title="Metric-affine gauge theory of gravity",
        authors="Hehl et al.",
    )
    assert "(PDF only)" in idx.render_table([row])


# --- end to end (offline) -----------------------------------------------------


def test_build_rows_is_offline_when_cache_is_seeded(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    lit = tmp_path / "literature"
    (lit / "2606.30785").mkdir(parents=True)
    (lit / "2606.30785" / "Manuscript.tex").write_text(r"\begin{document}x")
    (lit / "PhysRep258.1").mkdir()
    (lit / "PhysRep258.1" / "paper.pdf").write_bytes(b"%PDF")

    monkeypatch.setattr(
        idx,
        "CACHE",
        _seed_cache(
            tmp_path,
            {"2606.30785": {"title": "Numerical polology", "authors": "Barker et al."}},
        ),
    )

    def _no_network(*_args: object, **_kwargs: object) -> None:
        msg = "hit the network despite a seeded cache"
        raise AssertionError(msg)

    monkeypatch.setattr(idx.urllib.request, "urlopen", _no_network)

    rows = {r.directory: r for r in idx.build_rows(lit)}
    assert rows["2606.30785"].title == "Numerical polology"
    assert rows["2606.30785"].main_tex == "Manuscript.tex"
    # Pre-arXiv entries come from the pinned MANUAL table, not the API.
    assert rows["PhysRep258.1"].main_tex is None
    assert "Metric-affine" in rows["PhysRep258.1"].title


def test_splice_replaces_only_the_generated_region() -> None:
    doc = f"prose\n\n{idx.BEGIN_MARKER}\n\nold\n\n{idx.END_MARKER}\n\ntail\n"
    out = idx.splice(doc, "new")
    assert "old" not in out
    assert out.startswith("prose")
    assert out.endswith("tail\n")
    assert idx.splice(out, "new") == out  # idempotent


def test_splice_refuses_a_readme_without_markers() -> None:
    with pytest.raises(SystemExit, match="markers"):
        idx.splice("no markers here", "table")
