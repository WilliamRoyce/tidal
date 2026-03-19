"""Colored console output for the TIDAL CLI.

Mirrors the ``[INFO]``/``[WARN]``/``[ERROR]``/``[OK]`` pattern used in
the shell scripts (``scripts/verify-wolfram-setup.sh`` etc.).

**Accessibility principle:** Prefixes are ALWAYS present regardless of
colour capability.  Colour is an additive enhancement, never the sole
carrier of meaning.  Structural separators (``====``, ``---``) remain.
Output must be fully readable without colour.

All messages go to *stderr* so that structured data on stdout (``--json``,
piped output) is never contaminated.
"""

from __future__ import annotations

import os
import sys

# ── ANSI colour codes ────────────────────────────────────────────

_RESET = "\033[0m"
_BOLD = "\033[1m"
_DIM = "\033[2m"

_RED = "\033[0;31m"
_GREEN = "\033[0;32m"
_YELLOW = "\033[1;33m"
_BLUE = "\033[0;34m"
_CYAN = "\033[0;36m"

# ── Colour support detection ────────────────────────────────────
# Replicates the logic in tidal/banner.py but checks stderr (where
# all console output is directed).


def _supports_color() -> bool:
    """Return True when stderr supports ANSI colour codes."""
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not hasattr(sys.stderr, "isatty"):
        return False
    return sys.stderr.isatty()


# Module-level state set once by ``configure()`` in main().
_use_color: bool | None = None  # lazy-init on first call
_verbose: bool = False
_quiet: bool = False


def configure(*, verbose: bool = False, quiet: bool = False) -> None:
    """Set verbosity levels.  Called once from ``main()``."""
    global _verbose, _quiet  # noqa: PLW0603
    _verbose = verbose
    _quiet = quiet


def _color() -> bool:
    """Lazy colour detection (cached after first call)."""
    global _use_color  # noqa: PLW0603
    if _use_color is None:
        _use_color = _supports_color()
    return _use_color


def _styled(prefix: str, color: str, msg: str) -> str:
    """Format a prefixed message, optionally with colour."""
    if _color():
        return f"{color}{_BOLD}{prefix}{_RESET} {msg}"
    return f"{prefix} {msg}"


# ── Public API ───────────────────────────────────────────────────


def info(msg: str) -> None:
    """Print ``[INFO] msg`` in blue to stderr."""
    print(_styled("[INFO]", _BLUE, msg), file=sys.stderr)


def warn(msg: str) -> None:
    """Print ``[WARN] msg`` in yellow to stderr."""
    print(_styled("[WARN]", _YELLOW, msg), file=sys.stderr)


def error(msg: str) -> None:
    """Print ``[ERROR] msg`` in red to stderr."""
    print(_styled("[ERROR]", _RED, msg), file=sys.stderr)


def success(msg: str) -> None:
    """Print ``[OK] msg`` in green to stderr."""
    print(_styled("[OK]", _GREEN, msg), file=sys.stderr)


def debug(msg: str) -> None:
    """Print ``[DEBUG] msg`` in dim to stderr — only when ``--verbose``."""
    if not _verbose:
        return
    if _color():
        print(f"{_DIM}[DEBUG] {msg}{_RESET}", file=sys.stderr)
    else:
        print(f"[DEBUG] {msg}", file=sys.stderr)


def log(msg: str) -> None:
    """Print plain *msg* to stderr.  Suppressed by ``--quiet``."""
    if _quiet:
        return
    print(msg, file=sys.stderr)


def error_with_hint(msg: str, hints: list[str]) -> None:
    """Print an error followed by actionable hints."""
    error(msg)
    for h in hints:
        if _color():
            print(f"  {_DIM}Hint: {h}{_RESET}", file=sys.stderr)
        else:
            print(f"  Hint: {h}", file=sys.stderr)


def header(title: str) -> str:
    """Return a section header with ``====`` separators.

    Colour is additive — the ``====`` bars remain for non-colour contexts.
    """
    sep = "=" * 64
    if _color():
        return f"{_CYAN}{_BOLD}{sep}\n{title}\n{sep}{_RESET}"
    return f"{sep}\n{title}\n{sep}"


def pass_fail(label: str, *, passed: bool) -> str:
    """Return ``PASS``/``FAIL`` with colour when available."""
    tag = "PASS" if passed else "FAIL"
    if _color():
        c = _GREEN if passed else _RED
        return f"{label}: {c}{_BOLD}{tag}{_RESET}"
    return f"{label}: {tag}"
