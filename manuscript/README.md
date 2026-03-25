# Manuscript

This directory contains the masters project report/article. All contents except
this README and `.gitignore` are git-ignored to prevent the manuscript from
appearing in the public repository (plagiarism detection protection).

## Structure

```
main.tex            Master document (\input's all sections)
macros.tex          Packages + macros (supervisor's template + TIDAL notation)
references.bib      Bibliography (seeded from docs/tex/references.bib)
sections/
  abstract.tex      Abstract
  introduction.tex  Introduction
  theory.tex        Theory
  methods.tex       Methods
  results.tex       Results
  discussion.tex    Concluding remarks
  acknowledgements.tex
  appendices.tex    Appendix
figures/            Manuscript-specific figures
```

## Building

```bash
latexmk -pdf -pvc -interaction=nonstopmode main.tex
```

Requires `revtex4-2` (available in TeX Live: `texlive-publishers`).

## Repo resources

- **`docs/tex/`** — LaTeX-formatted equations and derivations from the project
  documentation. Useful as a source to copy content into the manuscript.
- **`docs/figures/`** — TikZ architecture and pipeline diagrams with a shared
  style file (`tidal-tikz-styles.sty`).
- **`docs/tex/references.bib`** — Canonical BibTeX references (manuscript copy
  was seeded from this).
