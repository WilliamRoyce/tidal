# Table-generation scripts for the manuscript

Each script reads a canonical benchmark JSON in
`benchmark_results/canonical/` and emits one `.tex` file to
`manuscript/tables/`. The appendix sources the resulting tabular body via
`\input{tables/<name>}`.

## Conventions

1. **One script per table.** Filename `tab_<short_name>.py` mirrors the
   `scripts/figures/figN_*.py` pattern.
2. **No measurements.** Tables only render. If the JSON is stale, rerun
   the corresponding benchmark in `scripts/benchmarks/`.
3. **Output is a `tabular` body** (not a full `table` float). The appendix
   provides `\begin{table}`, `\caption{}`, `\label{}`, and any
   surrounding prose — only the rows are auto-generated.
4. **booktabs style.** Use `\toprule`, `\midrule`, `\bottomrule` to match
   existing manuscript tables.
5. **Stable formatting.** Re-rendering with unchanged canonical data must
   produce a byte-identical `.tex` (no timestamps, no host metadata in
   the rendered output — those belong only in the JSON).

## What's here

| Script | Reads | Writes | Appendix line |
|---|---|---|---|
| `tab_pade_benchmark.py` | `benchmark_results/canonical/pade_vs_eig.json` | `manuscript/tables/pade_benchmark.tex` | C:348 |
| `tab_sparse_csc.py` | `benchmark_results/canonical/sparse_csc_vs_dense.json` | `manuscript/tables/sparse_csc_benchmark.tex` | C:490 |
| `tab_nyquist_energy.py` | `benchmark_results/canonical/nyquist_energy.json` | `manuscript/tables/nyquist_energy.tex` | C:525 |

## Regenerate a single table

    uv run python scripts/tables/tab_pade_benchmark.py

## Regenerate all tables

    make publication-tables
