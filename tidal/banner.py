"""
TIDAL CLI Banner — nabla + ripple rings with colour gradients.

Algorithmically rasterized from the TIDAL SVG logo using half-block
characters (▀ ▄ █) for universal terminal compatibility.

Usage:
    from tidal.banner import print_banner
    print_banner()                      # default ocean theme
    print_banner(theme="sunset")        # warm gradient
    print_banner(theme="plain")         # no colour (safe for pipes / CI)

Themes: "ocean", "ice", "sunset", "neon", "plain"
"""

from __future__ import annotations

import os
import sys
from typing import TextIO

# ────────────────────────────────────────────────────────────
#  Pre-composed layouts (no dynamic padding — immune to
#  variable-width Unicode rendering across terminals)
# ────────────────────────────────────────────────────────────

# Compact: icon left, TIDAL text right (pre-composed, ~65 display cols)
_LAYOUT_COMPACT = r"""
     ▗▄▄▄▄▄▄▄▄▖
      ▜█▛▀▀▜█▛         ████████╗██╗██████╗  █████╗ ██╗
       ▜█▖▗█▛          ╚══██╔══╝██║██╔══██╗██╔══██╗██║
     ▗  ▝██▘  ▖           ██║   ██║██║  ██║███████║██║
 ▗▞▘█▙   ▝▘   ▟█▝▚▖       ██║   ██║██║  ██║██╔══██║██║
 █▙ ▀▀████████▀▀ ▟█       ██║   ██║██████╔╝██║  ██║███████╗
 ▝▜█▙▄▄▄▄▄▄▄▄▄▄▟█▛▘       ╚═╝   ╚═╝╚═════╝ ╚═╝  ╚═╝╚══════╝
    ▝▀▀▀▀▀▀▀▀▀▀▘
""".strip("\n")

# Stacked: icon on top, TIDAL text below
_LAYOUT_STACKED = r"""
         ▗▄▄▄▄▄▄▄▄▖
          ▜█▛▀▀▜█▛
           ▜█▖▗█▛
         ▗  ▝██▘  ▖
     ▗▞▘█▙   ▝▘   ▟█▝▚▖
     █▙ ▀▀████████▀▀ ▟█
     ▝▜█▙▄▄▄▄▄▄▄▄▄▄▟█▛▘
        ▝▀▀▀▀▀▀▀▀▀▀▘

████████╗██╗██████╗  █████╗ ██╗
╚══██╔══╝██║██╔══██╗██╔══██╗██║
   ██║   ██║██║  ██║███████║██║
   ██║   ██║██║  ██║██╔══██║██║
   ██║   ██║██████╔╝██║  ██║███████╗
   ╚═╝   ╚═╝╚═════╝ ╚═╝  ╚═╝╚══════╝
""".strip("\n")

# Minimal: icon + small text
_LAYOUT_MINIMAL = r"""
     ▗▄▄▄▄▄▄▄▄▖
      ▜█▛▀▀▜█▛
       ▜█▖▗█▛
     ▗  ▝██▘  ▖
 ▗▞▘█▙   ▝▘   ▟█▝▚▖
 █▙ ▀▀████████▀▀ ▟█
 ▝▜█▙▄▄▄▄▄▄▄▄▄▄▟█▛▘
    ▝▀▀▀▀▀▀▀▀▀▀▘

      T I D A L
""".strip("\n")


def _build_compact() -> list[str]:
    """Icon left, text right — pre-composed."""
    return _LAYOUT_COMPACT.split("\n")


def _build_stacked() -> list[str]:
    """Icon on top, text below — pre-composed."""
    return _LAYOUT_STACKED.split("\n")


def _build_minimal() -> list[str]:
    """Icon + small text — pre-composed."""
    return _LAYOUT_MINIMAL.split("\n")


# ────────────────────────────────────────────────────────────
#  Colour themes
# ────────────────────────────────────────────────────────────

THEMES: dict[str, list[tuple[int, int, int]]] = {
    "ocean": [
        (100, 220, 255),
        (50, 180, 240),
        (0, 140, 220),
        (0, 100, 200),
        (0, 60, 170),
        (10, 30, 130),
    ],
    "ice": [
        (200, 230, 255),
        (140, 200, 255),
        (80, 170, 240),
        (40, 140, 220),
        (20, 100, 200),
        (10, 60, 160),
    ],
    "sunset": [
        (255, 100, 100),
        (255, 140, 70),
        (255, 180, 50),
        (255, 210, 60),
        (255, 230, 90),
        (255, 245, 130),
    ],
    "neon": [
        (0, 255, 200),
        (0, 220, 255),
        (80, 160, 255),
        (140, 100, 255),
        (200, 60, 255),
        (255, 40, 200),
    ],
    "plain": [],
}


# ────────────────────────────────────────────────────────────
#  Rendering
# ────────────────────────────────────────────────────────────


def _supports_color() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    if os.environ.get("FORCE_COLOR"):
        return True
    if not hasattr(sys.stdout, "isatty"):
        return False
    return sys.stdout.isatty()


def _lerp(
    a: tuple[int, int, int], b: tuple[int, int, int], t: float
) -> tuple[int, int, int]:
    return (
        int(a[0] + (b[0] - a[0]) * t),
        int(a[1] + (b[1] - a[1]) * t),
        int(a[2] + (b[2] - a[2]) * t),
    )


def _gradient_row(
    idx: int, total: int, stops: list[tuple[int, int, int]]
) -> tuple[int, int, int]:
    if not stops or total <= 1:
        return (200, 200, 200)
    t = idx / (total - 1)
    seg = t * (len(stops) - 1)
    lo = int(seg)
    hi = min(lo + 1, len(stops) - 1)
    return _lerp(stops[lo], stops[hi], seg - lo)


_RESET = "\033[0m"


def _colorize(lines: list[str], stops: list[tuple[int, int, int]]) -> str:
    n = len(lines)
    return "\n".join(
        f"\033[38;2;{r};{g};{b}m{line}{_RESET}"
        for i, line in enumerate(lines)
        for r, g, b in [_gradient_row(i, n, stops)]
    )


# ────────────────────────────────────────────────────────────
#  Public API
# ────────────────────────────────────────────────────────────

_MIN_COMPACT_COLS = 68
_MIN_STACKED_COLS = 38


def get_banner(theme: str = "ocean", layout: str = "auto") -> str:
    """Return the TIDAL ASCII banner as a string.

    Parameters
    ----------
    theme : str
        One of "ocean", "ice", "sunset", "neon", "plain".
    layout : str
        "compact"  — icon left, text right  (needs ~52 cols)
        "stacked"  — icon on top, text below (needs ~38 cols)
        "minimal"  — small icon + small text (needs ~20 cols)
        "auto"     — pick based on terminal width
    """
    if layout == "auto":
        try:
            cols = os.get_terminal_size().columns
        except (OSError, ValueError):
            cols = 80
        if cols >= _MIN_COMPACT_COLS:
            layout = "compact"
        elif cols >= _MIN_STACKED_COLS:
            layout = "stacked"
        else:
            layout = "minimal"

    builders = {
        "compact": _build_compact,
        "stacked": _build_stacked,
        "minimal": _build_minimal,
    }
    lines = builders.get(layout, _build_compact)()

    stops = THEMES.get(theme, THEMES["ocean"])
    if stops and _supports_color():
        return _colorize(lines, stops)
    return "\n".join(lines)


def print_banner(
    theme: str = "ocean", layout: str = "auto", file: TextIO | None = None
) -> None:
    """Print the TIDAL banner to *file* (default: sys.stderr)."""
    if file is None:
        file = sys.stderr
    print(get_banner(theme=theme, layout=layout), file=file)


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Preview TIDAL CLI banner")
    p.add_argument("--theme", choices=list(THEMES.keys()), default="ocean")
    p.add_argument(
        "--layout", choices=["compact", "stacked", "minimal", "auto"], default="auto"
    )
    a = p.parse_args()
    print_banner(theme=a.theme, layout=a.layout, file=sys.stdout)
