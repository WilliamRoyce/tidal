"""Static manifest checks: every entry references real files on disk."""

from __future__ import annotations

import pytest
import yaml

from .conftest import MANIFEST_PATH, REPO_ROOT, iter_artefacts

pytestmark = pytest.mark.publication


def test_manifest_loads() -> None:
    with MANIFEST_PATH.open() as fh:
        data = yaml.safe_load(fh)
    assert isinstance(data, dict)
    assert data, "manifest is empty"


@pytest.mark.parametrize(
    ("appendix", "name", "entry"),
    [
        pytest.param(*x, id=f"{x[0]}/{x[1]}")
        for x in iter_artefacts(yaml.safe_load(MANIFEST_PATH.read_text()))
    ],
)
def test_artefact_scripts_exist(appendix: str, name: str, entry: dict) -> None:  # noqa: ARG001
    if entry.get("status") == "stub":
        pytest.skip("entry marked stub")
    kind = entry["kind"]
    if kind == "figure":
        assert (REPO_ROOT / entry["figure_script"]).is_file(), entry["figure_script"]
    elif kind == "table":
        assert (REPO_ROOT / entry["table_script"]).is_file(), entry["table_script"]
    elif kind == "static_tikz":
        assert (REPO_ROOT / entry["source"]).is_file(), entry["source"]
    else:
        pytest.fail(f"unknown kind: {kind}")


@pytest.mark.parametrize(
    ("appendix", "name", "entry"),
    [
        pytest.param(*x, id=f"{x[0]}/{x[1]}")
        for x in iter_artefacts(
            yaml.safe_load(MANIFEST_PATH.read_text()), kinds=("figure", "table")
        )
    ],
)
def test_benchmark_scripts_exist(appendix: str, name: str, entry: dict) -> None:  # noqa: ARG001
    if entry.get("status") == "stub":
        pytest.skip("entry marked stub")
    bs = entry.get("benchmark_script")
    if bs is None:
        pytest.skip("no benchmark_script for this entry")
    assert (REPO_ROOT / bs).is_file(), bs
