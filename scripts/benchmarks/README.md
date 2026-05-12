# Benchmark scripts for publication artefacts

Each script in this directory measures one quantity used by an App C figure
or table in the manuscript. The convention mirrors `scripts/figures/`:

1. **One script per benchmark.** Filename describes the quantity measured
   (e.g. `fd_convergence.py`, `pade_vs_eig.py`).
2. **Each script reads no manuscript state and writes exactly one JSON
   file** to `benchmark_results/canonical/<name>.json` by default. The
   `--out PATH` flag overrides the destination for per-machine reruns.
3. **Self-contained.** Re-runnable from a fresh checkout via
   `uv run python scripts/benchmarks/<name>.py`.
4. **Metadata.** Every output JSON records: `timestamp`, `host`,
   `python_version`, `numpy_version`, `scipy_version`, `git_sha`, and the
   `parameters` swept.
5. **Optional-dependency guards.** If a backend is unavailable (e.g.
   SUNDIALS missing), the script records `skipped: <reason>` for that row
   rather than failing.

## What's here

| Script | Produces | Consumed by |
|---|---|---|
| `fd_convergence.py` | L₂ error vs N for FD orders 2/4/6 + spectral | `scripts/figures/figC2_fd_convergence.py` |
| `pade_vs_eig.py` | Padé vs eigendecomposition wall-times across theories | `scripts/tables/tab_pade_benchmark.py` |
| `sparse_csc_vs_dense.py` | Sparse-CSC vs dense `expm_multiply` timings | `scripts/tables/tab_sparse_csc.py` |
| `nyquist_energy.py` | |dE/E| sweep at N ∈ {128…1024} | `scripts/tables/tab_nyquist_energy.py` |
| `solver_convergence.py` | Modal solver convergence vs Boccaletti formula | `scripts/figures/figC1_solver_convergence.py` |

## Run a single benchmark

    uv run python scripts/benchmarks/fd_convergence.py

## Run all benchmarks (per-machine; writes to `benchmark_results/<host>-<date>/`)

    bash scripts/run_benchmarks.sh

## Promote results to canonical (commit-ready) data

    bash scripts/run_benchmarks.sh --canonical

The `--canonical` flag requires a clean git state and writes directly to
`benchmark_results/canonical/`. Reserve for "I am promoting these numbers
to the manuscript" — the figures/tables read from canonical data so this
gesture is what changes the appendix.

## Why separate benchmarks from figures

Benchmark reruns can be slow (minutes). Figure/table regeneration must be
fast (seconds) so reviewers regenerating from a fresh checkout get the
same artefacts without paying the benchmark cost. The committed canonical
JSON is the single source of truth; benchmark scripts exist to reproduce
those numbers, not to be re-run on every figure build.
