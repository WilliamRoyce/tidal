"""Colored console output for the TIDAL CLI.

Mirrors the ``[INFO]``/``[WARN]``/``[ERROR]``/``[OK]`` pattern used in
the shell scripts (``scripts/verify-wolfram-setup.sh`` etc.).

**Accessibility principle:** Prefixes are ALWAYS present regardless of
color capability.  Color is an additive enhancement, never the sole
carrier of meaning.  Structural separators (``====``, ``---``) remain.
Output must be fully readable without color.

All messages go to *stderr* so that structured data on stdout (``--json``,
piped output) is never contaminated.
"""

from __future__ import annotations

import os
import sys

# ── ANSI color codes ────────────────────────────────────────────

_RESET = "\033[0m"
_BOLD = "\033[1m"
_DIM = "\033[2m"

_RED = "\033[0;31m"
_GREEN = "\033[0;32m"
_YELLOW = "\033[1;33m"
_BLUE = "\033[0;34m"
_CYAN = "\033[0;36m"

# ── Color support detection ────────────────────────────────────
# Replicates the logic in tidal/banner.py but checks stderr (where
# all console output is directed).


def _supports_color() -> bool:
    """Return True when stderr supports ANSI color codes."""
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
    """Lazy color detection (cached after first call)."""
    global _use_color  # noqa: PLW0603
    if _use_color is None:
        _use_color = _supports_color()
    return _use_color


def _styled(prefix: str, color: str, msg: str) -> str:
    """Format a prefixed message, optionally with color."""
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

    Color is additive — the ``====`` bars remain for non-color contexts.
    """
    sep = "=" * 64
    if _color():
        return f"{_CYAN}{_BOLD}{sep}\n{title}\n{sep}{_RESET}"
    return f"{sep}\n{title}\n{sep}"


def key_value(key: str, val: str, *, width: int = 20) -> str:
    """Return an aligned key-value pair.  Key is colored when available."""
    if _color():
        return f"  {_CYAN}{key:<{width}}{_RESET} {val}"
    return f"  {key:<{width}} {val}"


def table(
    headers: list[str],
    rows: list[list[str]],
    *,
    indent: int = 2,
) -> str:
    """Return a simple aligned table with header underline.

    Color is additive — the table is readable without color.
    """
    # Compute column widths
    n_cols = len(headers)
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row[:n_cols]):
            widths[i] = max(widths[i], len(cell))

    pad = " " * indent
    # Header row
    hdr = pad + "  ".join(h.ljust(widths[i]) for i, h in enumerate(headers))
    sep = pad + "  ".join("-" * widths[i] for i in range(n_cols))
    # Data rows
    data_lines: list[str] = []
    for row in rows:
        cells: list[str] = [
            (row[i] if i < len(row) else "").ljust(widths[i]) for i in range(n_cols)
        ]
        data_lines.append(pad + "  ".join(cells))

    if _color():
        hdr = f"{_CYAN}{_BOLD}{hdr}{_RESET}"
        sep = f"{_DIM}{sep}{_RESET}"

    return "\n".join([hdr, sep, *data_lines])


def pass_fail(label: str, *, passed: bool) -> str:
    """Return ``PASS``/``FAIL`` with color when available."""
    tag = "PASS" if passed else "FAIL"
    if _color():
        c = _GREEN if passed else _RED
        return f"{label}: {c}{_BOLD}{tag}{_RESET}"
    return f"{label}: {tag}"
