r"""Render the Padé-vs-eigendecomposition wall-time table.

Reads:   benchmark_results/canonical/pade_vs_eig.json
Writes:  manuscript/tables/pade_benchmark.tex
Serves:  manuscript/sections/appendices/numerical.tex:348 (tab:PadeBenchmark)

Status: stub. Renders a ``\\begin{tabular}`` body. When the benchmark
script populates canonical/pade_vs_eig.json with the schema documented in
``scripts/benchmarks/pade_vs_eig.py``, this script reads it and emits one
row per result entry. Until then, it emits a TODO placeholder so the
manifest round-trip stays consistent.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "pade_vs_eig.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "tables" / "pade_benchmark.tex"


def _render_placeholder() -> str:
    return (
        "% Placeholder body: scripts/benchmarks/pade_vs_eig.py is a stub.\n"
        "% Implement the benchmark, then re-run this script.\n"
        "\\begin{tabular}{lcccc}\n"
        "  \\hline\\hline\n"
        "  Theory & $N$ & Eig.\\ (s) & Pad\\'e (s) & Ratio \\\\\n"
        "  \\hline\n"
        "  \\multicolumn{5}{c}{\\itshape benchmark pending} \\\\\n"
        "  \\hline\\hline\n"
        "\\end{tabular}\n"
    )


def _render_from_data(data: dict) -> str:
    rows: list[str] = [
        f"  \\lstinline!{r['theory']}! & "
        f"${r['n']}$ & "
        f"${r['eig_s']:.4g}$ & "
        f"${r['pade_s']:.4g}$ & "
        f"${r['ratio']:.2f}$ \\\\"
        for r in data["results"]
    ]
    body = "\n".join(rows)
    return (
        "\\begin{tabular}{lcccc}\n"
        "  \\hline\\hline\n"
        "  Theory & $N$ & Eig.\\ (s) & Pad\\'e (s) & Ratio \\\\\n"
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
