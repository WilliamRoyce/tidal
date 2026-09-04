"""Repository hygiene: no environment-specific absolute paths in tracked files.

TIDAL is distributed and installed by other people, so a committed path naming
one machine, one user, or one clone location works only for its author.  The
2026-08 repository rename exposed how far that had spread: ``devcontainer.json``
hardcoded ``/workspaces/torsion-gertsenshtein``, so cloning into a folder of any
other name produced a failing container build.

This module makes the ``CLAUDE.md`` convention self-enforcing rather than
advisory.  Derive roots instead of stating them -- see
``scripts/hpc_shuttle.sh`` for the shell idiom.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent

# Directories whose contents are checked.  ``manuscript/`` and ``research/`` are
# excluded on purpose: they hold dated provenance records (bib-audit
# "%! checked:" markers, captured tracebacks) that state where a file was read at
# a point in time.  Rewriting those would falsify the record.
CHECKED_PREFIXES = (
    "tidal/",
    "tidalcosmo/",
    "tests_cosmo/",
    "scripts/",
    "tests/",
    "examples/",
    ".claude/",
    ".devcontainer/",
    "docs/",
    ".github/",
)

# Patterns that indicate a path tied to one machine or clone location.
PATTERNS = {
    "container/clone path": re.compile(r"/workspaces/"),
    "user home path": re.compile(r"/home/[A-Za-z0-9_.-]+/"),
    "Claude project slug": re.compile(r"-workspaces-[A-Za-z0-9-]+"),
}

# Deliberate exceptions, each with the reason it cannot be derived.
ALLOWLIST: dict[str, str] = {
    # Mount targets and remoteEnv must match the image's user literally; the
    # devcontainer spec does not interpolate arbitrary variables there.
    ".devcontainer/devcontainer.json": "mount targets/remoteEnv must be literal",
    # Permission patterns are matched literally -- they cannot take variables.
    ".claude/settings.json": "permission patterns do not interpolate",
    # Documents the literal devcontainer.json mount block verbatim.
    ".devcontainer/docs/WOLFRAM_GUIDE.md": "quotes the literal mount config",
    # This file names the patterns it forbids.
    "tests/test_repo_hygiene.py": "defines the patterns themselves",
}


def _tracked_files() -> list[str]:
    """Return git-tracked paths, or skip when git is unavailable (sdist/wheel)."""
    git = shutil.which("git")
    if git is None:
        pytest.skip("git unavailable -- cannot enumerate tracked files")
    try:
        out = subprocess.run(
            [git, "-C", str(REPO_ROOT), "ls-files"],
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        pytest.skip("git unavailable -- cannot enumerate tracked files")
    return [line for line in out.splitlines() if line]


def _candidate_files() -> list[str]:
    return [
        path
        for path in _tracked_files()
        if path.startswith(CHECKED_PREFIXES) or "/" not in path
        if path not in ALLOWLIST
    ]


def test_no_environment_specific_paths() -> None:
    """No tracked file may hardcode a machine-, user- or clone-specific path."""
    violations: list[str] = []

    for rel in _candidate_files():
        full = REPO_ROOT / rel
        try:
            text = full.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue  # binary or unreadable -- nothing to check

        for lineno, line in enumerate(text.splitlines(), start=1):
            for label, pattern in PATTERNS.items():
                match = pattern.search(line)
                if match:
                    violations.append(f"{rel}:{lineno}: {label} -- {match.group(0)!r}")
                    break

    assert not violations, (
        "Environment-specific absolute paths found in tracked files.\n"
        "TIDAL must work from any clone location for any user. Derive the root "
        'instead: REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" in '
        "shell, relative paths in devcontainer lifecycle commands, "
        "$CLAUDE_PROJECT_DIR in hooks. See CLAUDE.md.\n\n" + "\n".join(violations)
    )


def test_allowlist_entries_still_exist() -> None:
    """An allowlist entry for a deleted file hides future regressions."""
    missing = [p for p in ALLOWLIST if not (REPO_ROOT / p).exists()]
    assert not missing, f"Allowlisted files no longer exist: {missing}"
