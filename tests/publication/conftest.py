"""Shared fixtures for publication-artefact tests.

These tests validate the publication-figure infrastructure:
    - the manifest is parseable
    - canonical benchmark JSON exists where the manifest expects it
    - figure/table scripts can be re-run and emit their declared outputs

They are skipped by default in the main test lane (`uv run pytest tests/`).
Run explicitly with `uv run pytest -m publication tests/publication/`.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO_ROOT / "scripts" / "publication_manifest.yaml"


@pytest.fixture(scope="session")
def manifest() -> dict:
    with MANIFEST_PATH.open() as fh:
        return yaml.safe_load(fh)


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


def iter_artefacts(
    manifest: dict, *, kinds: tuple[str, ...] = ()
) -> list[tuple[str, str, dict]]:
    out: list[tuple[str, str, dict]] = []
    for appendix, entries in manifest.items():
        for name, entry in entries.items():
            if kinds and entry.get("kind") not in kinds:
                continue
            out.append((appendix, name, entry))
    return out
