# Documentation Directory

This directory contains three kinds of document:

- **`docs/cosmology/`** — the **active program**: the Cobaya-extension direction (umbrella #488). Start at `COSMOLOGY_PROGRAM.md`.
- **`docs/tex/`** — living technical documentation in LaTeX, for inclusion in the project's Overleaf report. It describes the **legacy** `tidal/` pipeline, which still exists and still runs; it is accurate, not stale.
- **Markdown files** — project management: roadmaps, trackers, campaign records. Most of the roadmap material is **parked** under the 2026-08-29 pivot; each carries an amendment banner saying so.

## Cosmology Program (`docs/cosmology/`) — ACTIVE

The current direction (pivot 2026-08-29): a **Cobaya extension** evolving a new sector's
linear perturbations as **spectators** on a CAMB LCDM background, producing CMB observables
and real likelihoods. Umbrella **#488**.

**Read `COSMOLOGY_PROGRAM.md` first** — it carries the decisions register, the observable
ladder, the workstreams and the current state. The rest are the research records behind it,
produced by eight sequential handoff sessions (H1–H8).

| File | Content |
| ---- | ------- |
| `../COSMOLOGY_PROGRAM.md` | **The operational record.** Decisions D1–D9, observable ladder O0–O4, workstreams WS0–WS6, verification gates. Start here |
| `primer.md` | Physics orientation: what a Boltzmann code does, what a spectator field is, why this route works |
| `spectator_route.md` | H7 — the three conditions a spectator sector must satisfy, and what breaks if one fails |
| `torc_pipeline_audit.md` | H1 — audit of TorC (arXiv:2507.09228), the predecessor this program is framed against. Settles O1's scope (R1) and how the CAMB patch is made (R2) |
| `observable_ladder.md` | H2 — per-rung feasibility: capabilities, known-answer targets, costs, risks. Establishes that O2 and O3 are **different numerical problems** |
| `magnetic_field_background.md` | H2 — the primordial B-field model O3 needs, and the spectator condition no paper enforces on it |
| `solver_design.md` | H3 — two engines over one shared core; the integration-target decision (option iii); the O2 ladder, matrix-WKB design, benchmark protocol |
| `repo_reshape.md` | H4 — the strangler-fig migration: `tidalcosmo/` beside legacy `tidal/`, the capability ledger, milestones M0–M7, CI and oracle strategy |
| `spectrum_design.md` | H6 — the two-stage spectrum architecture; Schur-complement kinetic-matrix criterion (arXiv:2506.02111) as the primary algorithm |
| `stage1_engineering_plan.md` | H8 — Stage-1 engineering **study**: six live-source findings correcting H6 (#521, #522, #523) |
| `birefringence_notes.md` | O4 foundation: rotation formulae, the per-operator frequency-scaling problem, likelihood options |
| `planning_session_record.md` | Archived planning transcript. **Not authoritative** — kept for provenance only |
| `handoffs/H1–H6, H8.md` | The prompts each session was given. All complete; each carries a status header. Kept as the record of *what was asked* |

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
| `torsion.tex` | `docs/tex/torsion.tex` | Poincare gauge theory, torsion implementation |
| `pgt_stability_priors.tex` | `docs/tex/pgt_stability_priors.tex` | D2 stability priors with TIDAL/Blagojević/Barker cross-check |
| `chern_simons.tex` | `docs/tex/chern_simons.tex` | Chern-Simons 2+1D implementation |
| `amplification_mechanism.tex` | `docs/tex/amplification_mechanism.tex` | Light-mediator amplification at the stability boundary |
| `dark_photon_torsion.tex` | `docs/tex/dark_photon_torsion.tex` | Dark-photon torsion analogue |
| `gertsenshtein_plasma.tex` | `docs/tex/gertsenshtein_plasma.tex` | Plasma-mass Gertsenshtein conversion |

### Architecture

| File | Location | Content |
| ---- | -------- | ------- |
| `architecture.tex` | `docs/tex/architecture.tex` | Pipeline overview, module roles, component E-L, Ostrogradsky |
| `perturbative_reduction.tex` | `docs/tex/perturbative_reduction.tex` | v6 iterative order reduction: Pass 0 / Pass 1, Parker–Simon + FKY validity, closed-form Duhamel kernel, constraint-field Schur recovery, EH Power-normalization and matter-only CD precompute gate (issue #271) |
| `perturbative_reduction_design.tex` | `docs/tex/perturbative_reduction_design.tex` | Engineer-facing implementation specification: algorithm pseudocode, module layout, gate helpers, regression matrix |
| `json_schema.tex` | `docs/tex/json_schema.tex` | Complete JSON specification reference |
| `solver_migration.tex` | `docs/tex/solver_migration.tex` | py-pde to SUNDIALS migration |
| `modal_solver.tex` | `docs/tex/modal_solver.tex` | Fourier modal solver |
| `pencil_engine.tex` | `docs/tex/pencil_engine.tex` | Matrix-pencil engine: deflation contract, implicit-dynamical sector, gauge quotient, rejected alternatives, pin + certify and observable-sector closure (GH #457–#474) |
| `solver_optimizations.tex` | `docs/tex/solver_optimizations.tex` | FD stencils, Yoshida, spectral, component E-L |
| `adaptive_timestepping.tex` | `docs/tex/adaptive_timestepping.tex` | Tolerance-controlled solvers |
| `kinetic_matrix.tex` | `docs/tex/kinetic_matrix.tex` | Non-diagonal kinetic matrix handling |
| `modal_jax.tex` | `docs/tex/modal_jax.tex` | JAX modal backend |
| `stability_probe.tex` | `docs/tex/stability_probe.tex` | Pre-flight tachyonic stability probe |
| `plane_wave_ic.tex` | `docs/tex/plane_wave_ic.tex` | Plane-wave IC mode-snapping |
| `perturbative_reduction_constraint_barrier.tex` | `docs/tex/perturbative_reduction_constraint_barrier.tex` | Constraint-promotion barrier analysis |

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
| `PERTURBATIVE_REDUCTION_IMPLEMENTATION.md` | v6 perturbative reduction implementation log |
| `MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md` | Reverted template-cache experiment record |
| `AMPLIFICATION_INVESTIGATION.md` | Amplification-mechanism investigation log |
| `PHASE_6_COMPARISON.md` | Phase 6 comparison notes |

## Campaign & Results Records (Markdown)

These are internal bookkeeping for the survey campaign — the thesis
(`manuscript/`, frozen archive) is the authoritative record of results, and
`docs/RESULTS_AMENDMENTS.md` is the living record of what has since been
superseded.

| File | Purpose |
| ---- | ------- |
| `RESULTS_AMENDMENTS.md` | **Living corrected-results record**: claim-by-claim status of the archived thesis numbers (start here) |
| `dkl_recompute_report.md` | GH #420 marginal D_KL recompute: old-vs-new evidence for all publication chains |
| `campaign_plan.md` | Campaign stages and success thresholds (see its amendment pointer) |
| `V3_ARCHITECTURE.md` | v3 inference architecture reference (see its amendment banner) |
| `V3_PHASE_TRACKER.md` | v3 phase status tracker |
| `V3_2_DESIGN_INVESTIGATION.md` | v3.2 design investigation |
| `V3_D1_REPLAY_NOTES.md` | D1 replay notes |
| `PHASE_E_TRACKER.md` | Phase E (localized geometry) tracker |
| `PHASE_E_ATLAS_TRACKER.md` | Phase E cubed-sphere atlas tracker |
| `lagrangian_depruning_audit.md` | Lagrangian de-pruning audit |
| `V3_PHASE_C_REFERENCE.md` | v3 Phase C reference |
| `V3_PHASE_D_DESIGN.md` | v3 Phase D design |
| `V3_PHASE_E_DESIGN.md` | v3 Phase E design |
| `V3_PHASE_E_PROTOTYPE.md` | v3 Phase E prototype record |
| `hpc_workflow.md` | **Archived** CSD3 HPC workflow — inactive under D4 (no HPC without explicit permission); reactivate when authorized |

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

Last updated: 2026-09-04
