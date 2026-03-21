# Documentation Directory

This directory contains living documentation in LaTeX format (`docs/tex/`) for easy inclusion in the project's Overleaf report, plus project management files in Markdown.

## LaTeX Documentation (`docs/tex/`)

All technical documentation lives in `docs/tex/` as LaTeX fragments. Each file is self-contained (no `\documentclass`), starts with `\section{Title}\label{sec:slug}`, and uses macros from `preamble.tex`. To compile any fragment standalone:

```latex
\documentclass[11pt,a4paper]{article}
\input{preamble}
\begin{document}
\input{fragment_name}
\bibliographystyle{unsrt}
\bibliography{references}
\end{document}
```

### Infrastructure

| File | Purpose |
| ---- | ------- |
| `preamble.tex` | Shared packages and macros (amsmath, physics, tensor, listings, booktabs, siunitx) |
| `references.bib` | BibTeX database (Gertsenshtein, torsion, numerical methods, xAct) |

### Physics

| File | Source | Content |
| ---- | ------ | ------- |
| `gertsenshtein.tex` | gertsenshtein.md | Gertsenshtein effect: physics background, validation targets |
| `gertsenshtein_formula.tex` | gertsenshtein_formula.md | Conversion formula derivation, literature comparison |
| `gertsenshtein_localized.tex` | gertsenshtein_localized.md | Boccaletti formula, localized B-field scattering |
| `critical_field.tex` | critical_field.md | Critical field analysis, amplification factor |
| `torsion.tex` | torsion.md | Poincare gauge theory, torsion implementation |
| `chern_simons.tex` | chern-simons-notes.md | Chern-Simons 2+1D implementation |

### Architecture

| File | Source | Content |
| ---- | ------ | ------- |
| `architecture.tex` | architecture/README.md | Pipeline overview, module roles, TikZ figure refs |
| `json_schema.tex` | JSON_SCHEMA_GUIDE.md | Complete JSON specification reference |
| `solver_migration.tex` | solver_migration.md | py-pde to SUNDIALS migration |
| `modal_solver.tex` | modal_solver.md | Fourier modal solver |
| `solver_optimizations.tex` | solver_optimizations.md | FD stencils, Yoshida, spectral phases |
| `adaptive_timestepping.tex` | adaptive_timestepping.md | Tolerance-controlled solvers |
| `kinetic_matrix.tex` | kinetic_matrix_alternatives.md | Non-diagonal kinetic matrix handling |

### Features

| File | Source | Content |
| ---- | ------ | ------- |
| `background_fields.tex` | background_fields.md | Position-dependent coefficients |
| `constraint_fields.tex` | constraint_fields.md | Mixed time-derivative orders, DAE handling |
| `gauge_fixing.tex` | gauge_fixing.md | Per-field gauge presets |
| `multi_field_perturbation.tex` | multi_field_perturbation.md | Multi-field linearization (xPert) |

### Operational & User-Facing

| File | Source | Content |
| ---- | ------ | ------- |
| `troubleshooting.tex` | troubleshooting.md | Error encyclopedia |
| `cli_reference.tex` | source/cli.md | CLI subcommand reference |
| `pipeline.tex` | source/pipeline.md | Two-stage data flow |
| `examples.tex` | source/examples.md | Working examples catalog |
| `derivation_performance.tex` | derivation_performance.md | Wolfram bottleneck analysis |
| `adr_disk_storage.tex` | adr-disk-storage.md | ADR: mmap NumPy storage |
| `volume_element_fix.tex` | volume-element-fix.md | sqrt|g| volume element fix |

## TikZ Figures (`docs/figures/`)

17 standalone TikZ diagrams (pipeline, solvers, constraints, etc.) with shared styles in `tidal-tikz-styles.sty`. Each compiles independently with `\documentclass[border=10pt]{standalone}`.

## Project Management (Markdown)

| File | Purpose |
| ---- | ------- |
| `ROADMAP.md` | Feature roadmap |
| `NEXT_PHASES.md` | Implementation phases A-I |
| `COMMUNITY.md` | Support channels |
| `references.md` | Curated bibliography (browsable) |
| `next-features.md` | Sweep framework features |
| `torsion_implementation_checklist.md` | PGT implementation tracking |

## Sphinx API Docs (`docs/source/`)

Auto-generated API documentation via Sphinx (`.rst` files). Build with `make html` from `docs/`.

## Maintenance

- **Update immediately** when solving non-trivial bugs
- **Add patterns** after implementing new features
- New `.tex` files: follow `gertsenshtein_formula.tex` as template
- New example-specific docs: use `chern_simons.tex` as template

---

Last updated: 2026-03-21 (Migrated from Markdown to LaTeX)
