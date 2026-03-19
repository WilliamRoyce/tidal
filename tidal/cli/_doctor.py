"""``tidal doctor`` — Environment health check.

Checks Python version, key dependencies, Wolfram Engine availability,
xAct libraries, and example specifications. Inspired by ``brew doctor``
and ``flutter doctor``.
"""

from __future__ import annotations

import shutil
import subprocess  # noqa: S404
import sys
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from argparse import Namespace


def _check_python() -> tuple[bool, str]:
    v = sys.version_info
    label = f"Python {v.major}.{v.minor}.{v.micro}"
    if v >= (3, 10):
        return True, label
    return False, f"{label} (requires >= 3.10)"


def _check_package(name: str, import_name: str | None = None) -> tuple[bool, str]:
    """Check if a package is importable and return its version."""
    mod = import_name or name
    try:
        pkg = __import__(mod)
    except ImportError:
        return False, f"{name} not installed"
    else:
        ver = getattr(pkg, "__version__", "unknown")
        return True, f"{name} {ver}"


def _check_sundials() -> tuple[bool, str]:
    """Check scikit-sundae (SUNDIALS wrapper)."""
    try:
        import sksundae
    except ImportError:
        return False, "scikit-sundae not installed — IDA/CVODE solvers unavailable"
    else:
        ver = getattr(sksundae, "__version__", "unknown")
        return True, f"scikit-sundae {ver} (SUNDIALS IDA/CVODE)"


def _check_wolframscript() -> tuple[bool, str]:
    """Check wolframscript availability."""
    wscript = shutil.which("wolframscript")
    if wscript is None:
        return False, "wolframscript not found — 'tidal derive' will not work"
    try:
        result = subprocess.run(  # noqa: S603
            [wscript, "-code", "$VersionNumber"],
            check=False,
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False, "wolframscript found but timed out"
    else:
        if result.returncode == 0:
            ver = result.stdout.strip()
            return True, f"Wolfram Engine {ver}"
        return False, "wolframscript found but not functional (check license)"


def _check_xact() -> tuple[bool, str]:
    """Check xAct Mathematica library installation."""
    # xAct installs to ~/.Mathematica/Applications/xAct/
    xact_paths = [
        Path.home() / ".Mathematica" / "Applications" / "xAct",
        Path.home() / ".WolframEngine" / "Applications" / "xAct",
        Path("/usr/share/Mathematica/Applications/xAct"),
    ]
    for p in xact_paths:
        if p.is_dir():
            return True, f"xAct libraries ({p})"
    return False, "xAct libraries not found — 'tidal derive' needs xAct/xTensor"


def _check_examples() -> tuple[bool, str]:
    """Check for example JSON specifications."""
    data_dir = Path("examples/data")
    if not data_dir.is_dir():
        return False, "examples/data/ directory not found"
    specs = list(data_dir.glob("*.json"))
    if specs:
        return True, f"{len(specs)} example JSON specs in examples/data/"
    return False, "no JSON specs found in examples/data/"


def doctor_command(_args: Namespace) -> int:
    """Run environment health checks."""
    from tidal.cli._console import success as _success
    from tidal.cli._console import warn as _warn

    checks = [
        _check_python(),
        _check_package("numpy"),
        _check_package("scipy"),
        _check_package("matplotlib"),
        _check_sundials(),
        _check_package("tqdm"),
        _check_wolframscript(),
        _check_xact(),
        _check_examples(),
    ]

    all_ok = True
    for ok, msg in checks:
        if ok:
            _success(msg)
        else:
            _warn(msg)
            all_ok = False

    if all_ok:
        print("\nAll checks passed.")
    else:
        print("\nSome checks had warnings — see above.")

    return 0
