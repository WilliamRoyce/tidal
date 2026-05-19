"""Each table script renders to a non-empty .tex with a tabular body."""

from __future__ import annotations

import subprocess
import sys
from typing import TYPE_CHECKING

import pytest
import yaml

from .conftest import MANIFEST_PATH, REPO_ROOT, iter_artefacts

if TYPE_CHECKING:
    from pathlib import Path

pytestmark = pytest.mark.publication


@pytest.mark.parametrize(
    ("appendix", "name", "entry"),
    [
        pytest.param(*x, id=f"{x[0]}/{x[1]}")
        for x in iter_artefacts(
            yaml.safe_load(MANIFEST_PATH.read_text()), kinds=("table",)
        )
    ],
)
def test_table_script_runs(
    appendix: str,
    name: str,
    entry: dict,
    tmp_path: Path,
) -> None:
    canon = REPO_ROOT / entry["canonical_data"]
    if not canon.exists():
        pytest.skip(f"canonical data not yet produced: {entry['canonical_data']}")

    script = REPO_ROOT / entry["table_script"]
    out = tmp_path / "table.tex"

    result = subprocess.run(
        [sys.executable, str(script), "--out", str(out)],
        check=False,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.returncode == 0, f"{script} failed:\n{result.stderr}"
    assert out.exists(), f"{script} produced no output"
    assert out.stat().st_size > 0, f"{script} produced an empty output"
    body = out.read_text()
    assert "\\begin{tabular}" in body, f"{script} output is not a tabular: {body[:200]}"
    assert "\\end{tabular}" in body
