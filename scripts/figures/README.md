# Figure scripts for the MSci report

Each script in this directory generates exactly one PDF figure for the report.
The convention is:

1. One script per figure. Filename `figN_short_description.{py,tex}` where `N`
   is the figure number in the report (`fig1`, `fig2`, `fig3` for main-body;
   `figA1`, `figC1`, `figD1`, `figD2`, `figD3` for appendix figures).
2. Each script reads campaign data from `hpc_results/`, regression fixtures
   in `examples/`, or other documented sources, and writes one PDF to
   `manuscript/figures/`. It does nothing else (no logging, no side effects,
   no interactive output).
3. Each script is self-contained: a reader can run it from a fresh checkout
   to regenerate the figure exactly as it appears in the report.
4. The architecture diagram (figA1) is TikZ rather than matplotlib, so
   `figA1_architecture_tikz.tex` is a `.tex` file rather than `.py`.
5. Scripts are stubs at Phase 2 scaffold time. Each per-figure drafting
   session implements one stub.

Run a single figure:

    uv run python scripts/figures/fig1_boccaletti_validation.py

Each stub records its data source in the module docstring so the
implementation session can begin without re-deciding which campaign
feeds which figure.
