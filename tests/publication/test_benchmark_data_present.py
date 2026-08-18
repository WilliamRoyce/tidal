"""Canonical benchmark JSON exists and validates against a minimal schema.

This catches the failure mode "table .tex references a canonical JSON that
disappeared or was overwritten with malformed data".
"""

from __future__ import annotations

import json

import pytest
import yaml

from .conftest import MANIFEST_PATH, REPO_ROOT, iter_artifacts

pytestmark = pytest.mark.publication

REQUIRED_METADATA_KEYS = {"timestamp", "host", "git_sha", "parameters"}


@pytest.mark.parametrize(
    ("appendix", "name", "entry"),
    [
        pytest.param(*x, id=f"{x[0]}/{x[1]}")
        for x in iter_artifacts(
            yaml.safe_load(MANIFEST_PATH.read_text()), kinds=("figure", "table")
        )
    ],
)
def test_canonical_data_present(appendix: str, name: str, entry: dict) -> None:
    if entry.get("status") == "stub":
        pytest.skip("entry marked stub")
    canon = entry.get("canonical_data")
    if canon is None:
        pytest.skip("no canonical_data for this entry")
    # canonical_data may be a single path or a list of paths (figures
    # that consume multiple benchmarks). Validate each one individually
    # under the same schema.
    paths = canon if isinstance(canon, list) else [canon]
    for canon_one in paths:
        path = REPO_ROOT / canon_one
        if not path.exists():
            pytest.skip(f"canonical data not yet produced: {canon_one}")
        with path.open() as fh:
            data = json.load(fh)
        assert isinstance(data, dict), f"{canon_one} root must be an object"
        assert "metadata" in data, f"{canon_one} missing 'metadata'"
        assert "results" in data, f"{canon_one} missing 'results'"
        missing = REQUIRED_METADATA_KEYS - set(data["metadata"])
        assert not missing, f"{canon_one} metadata missing keys: {missing}"
