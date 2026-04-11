#!/usr/bin/env python3
r"""Analyze Phase 1 results: 1D parameter slices.

Thin wrapper around ``tidal plot campaign``. The full analysis logic
has been integrated into the CLI plotting infrastructure.

Usage:
  python analyze_phase1.py [PHASE1_DIR]

Equivalent CLI command:
  tidal plot campaign PHASE1_DIR --type campaign --metric log10_A \\
      --scatter --annotate-extremes --log-y --classify "active:log10_A>0.01"
"""

from __future__ import annotations

import subprocess  # noqa: S404
import sys
from pathlib import Path

PHASE1_DIR = (
    Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_sweeps_v2/phase1")  # noqa: S108
)


def main() -> None:
    cmd = [
        "uv",
        "run",
        "tidal",
        "plot",
        PHASE1_DIR.as_posix(),
        "--type",
        "campaign",
        "--metric",
        "log10_A",
        "--layout",
        "3x3",
        "--scatter",
        "--annotate-extremes",
        "--log-y",
        "--classify",
        "active:log10_A>0.01",
        "--output",
        str(PHASE1_DIR / "phase1_overview.png"),
    ]
    print(f"Running: {' '.join(str(c) for c in cmd)}")
    subprocess.run(cmd, check=True)  # noqa: S603


if __name__ == "__main__":
    main()
