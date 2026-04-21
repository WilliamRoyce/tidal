# Documentation Directory

This directory contains living documentation in LaTeX format (`docs/tex/`) for easy inclusion in the project's Overleaf report, plus project management files in Markdown.

## LaTeX Documentation (`docs/tex/`)

All technical documentation lives in `docs/tex/` as LaTeX fragments. Each file is self-contained (no `\documentclass`), starts with `\section{Title}\label{sec:slug}`, and uses macros from `preamble.tex`. These `.tex` files are the primary documentation — update them directly. To compile any fragment standalone:

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

| File | Location | Content |
| ---- | -------- | ------- |
| `gertsenshtein.tex` | `docs/tex/gertsenshtein.tex` | Gertsenshtein effect: physics background, validation targets |
| `gertsenshtein_formula.tex` | `docs/tex/gertsenshtein_formula.tex` | Conversion formula derivation, literature comparison |
| `gertsenshtein_localized.tex` | `docs/tex/gertsenshtein_localized.tex` | Boccaletti formula, localized B-field scattering |
| `background_validity.tex` | `docs/tex/background_validity.tex` | Background validity, B₀→0 argument, EFT structure, sweep methodology |
| `critical_field.tex` | `docs/tex/critical_field.tex` | Critical field analysis, amplification factor |
| `torsion.tex` | `docs/tex/torsion.tex` | Poincare gauge theory, torsion implementation |
| `chern_simons.tex` | `docs/tex/chern_simons.tex` | Chern-Simons 2+1D implementation |

### Architecture

| File | Location | Content |
| ---- | -------- | ------- |
| `architecture.tex` | `docs/tex/architecture.tex` | Pipeline overview, module roles, component E-L, Ostrogradsky |
| `perturbative_reduction.tex` | `docs/tex/perturbative_reduction.tex` | v6 iterative order reduction: Pass 0 / Pass 1, Parker–Simon + FKY validity, closed-form Duhamel kernel, constraint-field Schur recovery, EH Power-normalisation and matter-only CD precompute gate (issue #271) |
| `perturbative_reduction_design.tex` | `docs/tex/perturbative_reduction_design.tex` | Engineer-facing implementation specification: algorithm pseudocode, module layout, gate helpers, regression matrix |
| `json_schema.tex` | `docs/tex/json_schema.tex` | Complete JSON specification reference |
| `solver_migration.tex` | `docs/tex/solver_migration.tex` | py-pde to SUNDIALS migration |
| `modal_solver.tex` | `docs/tex/modal_solver.tex` | Fourier modal solver |
| `solver_optimizations.tex` | `docs/tex/solver_optimizations.tex` | FD stencils, Yoshida, spectral, component E-L |
| `adaptive_timestepping.tex` | `docs/tex/adaptive_timestepping.tex` | Tolerance-controlled solvers |
| `kinetic_matrix.tex` | `docs/tex/kinetic_matrix.tex` | Non-diagonal kinetic matrix handling |

### Features

| File | Location | Content |
| ---- | -------- | ------- |
| `background_fields.tex` | `docs/tex/background_fields.tex` | Position-dependent coefficients |
| `constraint_fields.tex` | `docs/tex/constraint_fields.tex` | Mixed time-derivative orders, DAE handling |
| `gauge_fixing.tex` | `docs/tex/gauge_fixing.tex` | Per-field gauge presets |
| `multi_field_perturbation.tex` | `docs/tex/multi_field_perturbation.tex` | Multi-field linearization (xPert) |

### Operational & User-Facing

| File | Location | Content |
| ---- | -------- | ------- |
| `inference.tex` | `docs/tex/inference.tex` | Bayesian inference: priors, constraints, MC and nested sampling, posterior analysis |
| `troubleshooting.tex` | `docs/tex/troubleshooting.tex` | Error encyclopedia |
| `cli_reference.tex` | `docs/tex/cli_reference.tex` | CLI subcommand reference |
| `pipeline.tex` | `docs/tex/pipeline.tex` | Two-stage data flow |
| `examples.tex` | `docs/tex/examples.tex` | Working examples catalog |
| `derivation_performance.tex` | `docs/tex/derivation_performance.tex` | Wolfram bottleneck analysis, component E-L timings |
| `adr_disk_storage.tex` | `docs/tex/adr_disk_storage.tex` | ADR: mmap NumPy storage |
| `volume_element_fix.tex` | `docs/tex/volume_element_fix.tex` | sqrt|g| volume element fix |

## TikZ Figures (`docs/figures/`)

18 standalone TikZ diagrams (pipeline, solvers, constraints, etc.) with shared styles in `tidal-tikz-styles.sty`. Each compiles independently with `\documentclass[border=10pt]{standalone}`.

## Project Management (Markdown)

| File | Purpose |
| ---- | ------- |
| `ROADMAP.md` | Feature roadmap |
| `NEXT_PHASES.md` | Implementation phases A-I |
| `COMMUNITY.md` | Support channels |
| `references.md` | Curated bibliography (browsable) |
| `next-features.md` | Sweep framework features |
| `torsion_implementation_checklist.md` | PGT implementation tracking |

## Research (`research/`)

Systematic enumeration of the most general quadratic PGT+EM Lagrangian using xAct/xTras.

| File | Content |
| ---- | ------- |
| `general_quadratic_lagrangian.tex` | Complete enumeration: 35 core couplings + derivative extensions |
| `general_quadratic_lagrangian.wls` | xTras `MakeContractionAnsatz` enumeration script |
| `make_ansatz.wls` | Core quadratic ansatz generation |
| `classify_sectors.wls` | Ghost/parity/mixing classification |
| `check_constraints.wls` | DDI analysis and projective invariance |
| `enumeration_physical.json` | Sector classification with physics metadata |
| `enumeration_classified.json` | Full classification with ghost analysis |
| `enumeration_results.json` | Term counts by interaction type |

## Sphinx API Docs (`docs/source/`)

Auto-generated API documentation via Sphinx (`.rst` files). Build with `make html` from `docs/`.

## Maintenance

- **Update immediately** when solving non-trivial bugs
- **Add patterns** after implementing new features
- New `.tex` files: follow `gertsenshtein_formula.tex` as template
- New example-specific docs: use `chern_simons.tex` as template

---

Last updated: 2026-04-02
