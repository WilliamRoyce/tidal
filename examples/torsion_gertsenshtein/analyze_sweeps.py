#!/usr/bin/env python3
r"""Cross-phase sweep analysis and stability classification.

Thin wrapper around ``tidal plot`` and ``tidal analyze``. The full
analysis logic has been integrated into the CLI plotting infrastructure.

Usage:
  python analyze_sweeps.py [DATA_DIR]

Equivalent CLI commands:
  tidal plot PHASE2_DIR --type sweep-divergence --output divergence.png
  tidal plot PHASE3_DIR --type sweep-scatter --metric log10_A \\
      --divergent-center 0 --cmap RdYlBu_r --log-diagonal
"""

from __future__ import annotations

import subprocess  # noqa: S404
import sys
from pathlib import Path

DATA_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/tidal_sweeps_v2")  # noqa: S108


def main() -> None:
    phase3 = DATA_DIR / "phase3" / "lhs1000"
    if phase3.exists():
        cmd = [
            "uv",
            "run",
            "tidal",
            "plot",
            phase3.as_posix(),
            "--type",
            "sweep-scatter",
            "--metric",
            "log10_A",
            "--divergent-center",
            "0",
            "--cmap",
            "RdYlBu_r",
            "--log-diagonal",
            "--params",
            "delta1,zeta2,zeta3",
            "--output",
            str(DATA_DIR / "phase3_scatter.png"),
        ]
        print(f"Running: {' '.join(str(c) for c in cmd)}")
        subprocess.run(cmd, check=True)  # noqa: S603
    else:
        print(f"Phase 3 data not found at {phase3}")


if __name__ == "__main__":
    main()
