"""Each figure script produces a PDF when invoked from a clean tempdir.

Marked `publication_slow` because matplotlib + PDF backend is slower than
the table/data tests. Run with:

    uv run pytest -m "publication and publication_slow" tests/publication/
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from typing import TYPE_CHECKING

import pytest
import yaml

from .conftest import MANIFEST_PATH, REPO_ROOT, iter_artefacts

if TYPE_CHECKING:
    from pathlib import Path

pytestmark = [
    pytest.mark.publication,
    pytest.mark.publication_slow,
    pytest.mark.skipif(
        shutil.which("latex") is None,
        reason="latex binary not installed; figure scripts use matplotlib usetex=True",
    ),
]


@pytest.mark.parametrize(
    ("appendix", "name", "entry"),
    [
        pytest.param(*x, id=f"{x[0]}/{x[1]}")
        for x in iter_artefacts(
            yaml.safe_load(MANIFEST_PATH.read_text()), kinds=("figure",)
        )
    ],
)
def test_figure_script_runs(
    appendix: str,
    name: str,
    entry: dict,
    tmp_path: Path,
) -> None:
    if entry.get("status") == "stub":
        pytest.skip(f"entry marked stub: {entry.get('figure_script')}")

    canon = entry["canonical_data"]
    paths = canon if isinstance(canon, list) else [canon]
    for canon_one in paths:
        if not (REPO_ROOT / canon_one).exists():
            pytest.skip(f"canonical data not yet produced: {canon_one}")

    script = REPO_ROOT / entry["figure_script"]
    out = tmp_path / "figure.pdf"

    result = subprocess.run(
        [sys.executable, str(script), "--out", str(out)],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, f"{script} failed:\n{result.stderr}"
    assert out.exists(), f"{script} produced no PDF"
    assert out.stat().st_size > 0, f"{script} produced an empty PDF"
    assert out.read_bytes()[:4] == b"%PDF", f"{script} output is not a PDF"
