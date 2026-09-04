"""The strangler-fig boundary, made self-enforcing.

``tidalcosmo`` is built beside legacy ``tidal``, not on top of it: the programme's
hardest structural rule is that **new code never imports old code** (D3, see
``docs/COSMOLOGY_PROGRAM.md`` and ``docs/cosmology/repo_reshape.md`` §8).  Useful
legacy capabilities are *ported* -- redesigned and moved, with their docstrings and
issue references -- never reached into.

``tidalcosmo/README.md`` stated this as already enforced.  It was not: no test
checked it, and ``pyrightconfig.json``/coverage/``testpaths`` all skipped the tree,
so the WS1 gate ("suite green, no old-code imports, pyright clean") would have
passed vacuously on an unchecked package.  This module closes that gap before the
first ``.py`` file lands, which is the cheapest moment to do it.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
NEW_PACKAGE = REPO_ROOT / "tidalcosmo"

# ``import tidal`` / ``from tidal import ...`` / ``from tidal.x import ...`` but NOT
# ``import tidalcosmo`` -- the word boundary is what separates the two packages.
LEGACY_IMPORT = re.compile(r"^\s*(?:from|import)\s+tidal\b(?!cosmo)", re.MULTILINE)


def test_new_package_never_imports_legacy() -> None:
    """No module under ``tidalcosmo/`` may import from legacy ``tidal``."""
    violations: list[str] = []

    for path in sorted(NEW_PACKAGE.rglob("*.py")):
        text = path.read_text(encoding="utf-8")
        for match in LEGACY_IMPORT.finditer(text):
            lineno = text.count("\n", 0, match.start()) + 1
            rel = path.relative_to(REPO_ROOT)
            violations.append(f"{rel}:{lineno}: {match.group(0).strip()}")

    assert not violations, (
        "tidalcosmo must never import legacy tidal -- port the capability instead "
        "(docs/cosmology/repo_reshape.md section 8):\n  " + "\n  ".join(violations)
    )
