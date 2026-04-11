#!/usr/bin/env python3
r"""Analyze dark photon torsion sweep results using C₀ = P/B₀² metric.

Thin wrapper around ``tidal plot``. The full C₀ analysis logic has been
integrated into the CLI via derived metrics (``--metric C0`` or ``--metric A``).

Usage:
    python analyze_sweep.py <sweep_dir> [--output <png_path>]

Equivalent CLI commands:
  # 1D sweep with C0 metric:
  tidal plot SWEEP_DIR --type sweep --metric C0 --log-y \\
      --hline "625.0:C_EM (GR):green" --output analysis_C0.png

  # With arctan axes for unbounded parameter range:
  tidal plot SWEEP_DIR --type sweep --metric A --arctan-axes --output analysis_A.png

Refs:
    Hwang & Noh (2310.04150) — linearized regime validity
    Domcke & Garcia-Cely (2301.02072) — Gertsenshtein thin-magnet formula
"""

from __future__ import annotations

import argparse
import subprocess  # noqa: S404


def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze dark photon sweep")
    parser.add_argument("sweep_dir", help="Sweep output directory")
    parser.add_argument("--output", default=None, help="Output PNG path")
    args = parser.parse_args()

    output = args.output or f"{args.sweep_dir}/analysis_C0.png"
    cmd = [
        "uv",
        "run",
        "tidal",
        "plot",
        args.sweep_dir,
        "--type",
        "sweep",
        "--metric",
        "C0",
        "--log-y",
        "--output",
        output,
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)  # noqa: S603


if __name__ == "__main__":
    main()
