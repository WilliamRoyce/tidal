"""Render the Nyquist energy-drift |dE/E| table.

Reads:   benchmark_results/canonical/nyquist_energy.json
Writes:  manuscript/tables/nyquist_energy.tex
Serves:  manuscript/sections/appendices/numerical.tex:525 (tab:NyquistEnergy)

Status: stub — see scripts/benchmarks/nyquist_energy.py.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA = REPO_ROOT / "benchmark_results" / "canonical" / "nyquist_energy.json"
DEFAULT_OUT = REPO_ROOT / "manuscript" / "tables" / "nyquist_energy.tex"


def _render_placeholder() -> str:
    return (
        "% Placeholder body: scripts/benchmarks/nyquist_energy.py is a stub.\n"
        "\\begin{tabular}{lcc}\n"
        "  \\hline\\hline\n"
        "  Theory & $N$ & $|dE/E|$ \\\\\n"
        "  \\hline\n"
        "  \\multicolumn{3}{c}{\\itshape benchmark pending} \\\\\n"
        "  \\hline\\hline\n"
        "\\end{tabular}\n"
    )


def _render_from_data(data: dict) -> str:
    rows: list[str] = [
        f"  \\lstinline!{r['theory']}! & ${r['n']}$ & ${r['abs_dE_over_E']:.2g}$ \\\\"
        for r in data["results"]
    ]
    body = "\n".join(rows)
    return (
        "\\begin{tabular}{lcc}\n"
        "  \\hline\\hline\n"
        "  Theory & $N$ & $|dE/E|$ \\\\\n"
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
