"""Render the sparse-CSC vs dense ``expm_multiply`` wall-time table.

Reads:   benchmark_results/canonical/sparse_csc_vs_dense.json
Writes:  manuscript/tables/sparse_csc_benchmark.tex
Serves:  manuscript/sections/appendices/numerical.tex:490 (tab:SparseCSCBenchmark)

Status: stub — see scripts/benchmarks/sparse_csc_vs_dense.py.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = (
    REPO_ROOT / "benchmark_results" / "canonical" / "sparse_csc_vs_dense.json"
)
DEFAULT_OUT = REPO_ROOT / "manuscript" / "tables" / "sparse_csc_benchmark.tex"


def _render_placeholder() -> str:
    return (
        "% Placeholder body: scripts/benchmarks/sparse_csc_vs_dense.py is a stub.\n"
        "\\begin{tabular}{ccccc}\n"
        "  \\hline\\hline\n"
        "  $N$ & density (\\%) & dense (s) & sparse (s) & speedup \\\\\n"
        "  \\hline\n"
        "  \\multicolumn{5}{c}{\\itshape benchmark pending} \\\\\n"
        "  \\hline\\hline\n"
        "\\end{tabular}\n"
    )


def _render_from_data(data: dict) -> str:
    rows: list[str] = [
        f"  ${r['n']}$ & ${r['density_pct']:.2f}$ & "
        f"${r['dense_s']:.3g}$ & ${r['sparse_s']:.3g}$ & "
        f"${r['speedup']:.1f}\\times$ \\\\"
        for r in data["results"]
    ]
    body = "\n".join(rows)
    return (
        "\\begin{tabular}{ccccc}\n"
        "  \\hline\\hline\n"
        "  $N$ & density (\\%) & dense (s) & sparse (s) & speedup \\\\\n"
        "  \\hline\n"
        f"{body}\n"
        "  \\hline\\hline\n"
        "\\end{tabular}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    if args.data.exists():
        with args.data.open() as fh:
            data = json.load(fh)
        body = _render_from_data(data)
    else:
        body = _render_placeholder()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(body)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
