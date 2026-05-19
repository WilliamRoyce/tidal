# MSci Report — Phase 0 & Phase 1 Checklist

Live tracking document. Tick items as completed. See `report_plan.md` for full context.

---

## Pre-execution (done first, before downloading any papers)

- [x] Write `manuscript/planning/report_plan.md` — repo-local archive of the full plan
- [x] Write `manuscript/planning/report_checklist.md` — this file
- [x] Write `manuscript/planning/literature_content_notes.md` — per-paper notes from reading all 23 locally-available papers
- [x] Update `docs/MEMORY.md` with pointer to planning files
- [x] Commit: `chore(report): add report plan archive and Phase 0/1 checklist` (582685e)

---

## Phase 0 — Fetch references and freeze style intel

### Tier 1 — Barker software-package papers (highest priority)

- [x] Download `literature/2406.09500/` — PSALTer original (125 .tex files)
- [x] Download `literature/2506.02111/` — PSALTer v2 parity-violating (26 .tex files)
- [x] Download `literature/2512.25007/` — Hamilcar (103 .tex files)

### Tier 2 — BHL joint corpus on PGT/torsion (highest priority for main body)

- [x] Verify `literature/2406.12826/` — *Every PGT is conformal* (4 .tex files, already local)
- [x] Download `literature/2309.14783/` — *Manifestly covariant variational principle* (1 .tex, single-file gz)
- [x] Download `literature/2303.11094/` — *Gravitational confinement and galactic rotation curves* (5 .tex files)
- [x] Download `literature/2101.02645/` — *Nonlinear Hamiltonian analysis, quadratic torsion I* (4 .tex files)
- [x] Download `literature/2006.03581/` — *PGT cosmology → Horndeski* (2 .tex files)
- [x] Download `literature/2003.02690/` — *H₀ tension with emergent dark radiation* (14 .tex files)

### Tier 3 — Hobson + Lasenby joint, close-topic

- [x] Download `literature/2008.09053/` — *Fresh perspective on gauging the conformal group* (1 .tex, single-file gz)
- [x] Download `literature/2005.02228/` — *Ghost and tachyon free Weyl gauge theories* (1 .tex file)

### Tier 4 — Barker solo voice

- [x] Download `literature/2205.13534/` — solo-author multipliers paper (2 .tex files)
- [x] Download `literature/2311.11790/` — solo-author Yang-Mills gravity (2 .tex files)

### Tier 5 — Additional Barker papers (recent metric-affine and PGT work)

- [x] Download `literature/2510.17094/` — Gertsenshtein effect on curved spacetime (1 .tex, DIRECT relevance)
- [x] Download `literature/2507.09228/` — TorC torsion condensation cosmology (Handley+Barker, 2 .tex)
- [x] Verify `literature/2510.08201/` — R² Hamiltonian analysis (already local, 1 .tex)
- [x] Download `literature/2506.21662/` — Infrared foundations I: rank-3 field theories (5 .tex)
- [x] Download `literature/2507.05349/` — Infrared foundations II: torsion-like + ghost analysis (5 .tex)
- [x] Download `literature/2505.23894/` — Can metric-affine gravity be saved? (2 .tex)
- [x] Download `literature/2402.07641/` — Particle spectra of Palatini/metric-affine (2 .tex)
- [x] Download `literature/2402.14917/` — Metric-affine gravity from extended projective symmetry (4 .tex)

### Tier 6 — Key collaborator / field papers

- [x] Download `literature/2406.11956/` — Weyl-invariant Einstein-Cartan gravity (Karananas et al., 1 .tex)
- [x] Download `literature/2506.17017/` — Cosmology of Cubic PGT (Bahamonde et al., 1 .tex)

### Already-in-`literature/` papers that must be studied for CONTENT (not yet analysed)

These are local — no download needed, but must be read for §1, §2, §4, §5:
- [x] Read and note `literature/1812.02675/` — Lin, Hobson, Lasenby: ghost-free PGT classification (§1.2, §5.2)
- [x] Read and note `literature/1804.05556/` — Blagojević & Cvetković PGT Hamiltonian structure (§5.2)
- [x] Read and note `literature/2301.02072/` — Simple derivation of Gertsenshtein effect (§1.1, §4)
- [x] Read and note `literature/2310.04150/` — EM field definitions in graviton-photon conversion (§2.3)
- [x] Read and note `literature/2004.02714/` — Exact graviton-photon mixing in B-field (Ejlli, §4)
- [x] Read and note `literature/2105.04565/` — Dark Photon Limits Handbook (Caputo et al., §4.2)
- [x] Read and note `literature/hep-th_0103093/` — Shapiro: Physical aspects of torsion (§2 background)
- [x] Read and note `literature/gr-qc_0001010/` — Hehl & Obukhov: How does EM couple to gravity? (§2)
- [x] Read and note `literature/gr-qc_0307063/` — Itin & Hehl: Maxwell + quadratic torsion (§2, §5)

### Pre-arXiv must-cite papers (bibliography entries only)

These need entries in `references.bib` — no TeX to download:
- [x] `gertsenshtein1962wave` — already present in references.bib
- [x] `boccaletti1970conversion` — already present in references.bib
- [x] `raffelt1988mixing` — already present in references.bib
- [x] `kibble1961lorentz` — added (J. Math. Phys. 2, 212 (1961))
- [x] `hehl1976spin` — already present in references.bib
- [x] `sezgin1980ghostfree` — added (PRD 21, 3269 (1980))
- [x] `cillis1996graviton` — added (PRD 54, 4757 (1996))

### Style intel artefacts

- [x] Create `manuscript/style_intel/` directory
- [x] Write `manuscript/style_intel/style_analysis.md` — 117 lines, covers all 10 patterns (a–j), analytical prose
- [x] Write `manuscript/style_intel/genre_conventions.md` — covers (a–g) + BHL physics-calculation sections (h–k) + cross-paper evolution notes
- [x] Write `manuscript/style_intel/genre_conventions_theory.md` — comprehensive theoretical paper genre (18 papers across 3 clusters): 6 paper sub-types, title/abstract/intro/theory/results/discussion/appendix conventions, Gertsenshtein-specific patterns, agent-use summary
- [x] Create `manuscript/STYLE_GUIDE.md` — actionable checklist (a–j), per-section templates §1–§5, self-review checklist

### Commit

- [x] Commit: `chore(report): freeze style analysis and supervisor reference corpus for MSci draft`

---

## Phase 1 — Outline approval, notation set, bibliography

### Supervisor sign-off

- [ ] Send outline (from `report_plan.md` § "Proposed report outline") to Barker for sign-off
- [ ] Receive sign-off (or incorporate feedback and re-send)

### Bibliography hygiene

- [x] Add `barker2024psalter` (2406.09500) to `manuscript/references.bib`
- [x] Add `barker2025psalter2` (2506.02111) to `manuscript/references.bib`
- [x] Add `barker2025hamilcar` (2512.25007) to `manuscript/references.bib`
- [x] Add software comparators: NRPy+, OGRe, Cadabra2 (×2 entries); Dedalus already present; xAct/xPert already present
- [x] Add dark-photon / hidden-photon review (`caputo2021darkphoton`, PRD 104, 095029 (2021))
- [x] Add ghost-instability PGT literature: `sezgin1980ghostfree` (PRD 21, 3269); `lin2019pgt` (Lin-Hobson-Lasenby 2019); `blagojevic2018hamiltonian` — covers the key no-ghost constraints
- [x] Add GW-detection context anchors: `amaro2017lisa` (LISA); `nanograv2023fifteen` (NANOGrav 15-yr); `abbott2023gwtc3` (GWTC-3 / LIGO O4 context)
- [x] Verify entry count ≥ 55 — **81 entries** total

### Notation survey and macros

- [x] Read Tier 1–2 papers focusing on: metric signature, vierbein/spin-connection notation, torsion irreducible decomposition labels, Gertsenshtein coupling symbols
- [x] Record standard choices + source traceback in `manuscript/macros.tex` comments (see §"Agreed notation set" block added to macros.tex)
- [x] Agree macro set for all report symbols: metric (-,+,+,+); indices μ,ν,ρ (spacetime) / a,b,c (frame); e^a_μ vierbein; ω^{ab}_μ spin connection; T_μ torsion trace; S^μ axial; q_{μνρ} tensor; κ²=16πG; P=sin²(κB₀D/2)
- [x] Check for conflicts with `docs/tex/preamble.tex` — no conflicts; both use `\varkappa` for κ, same torsion/contorsion names

### Commit

- [ ] Commit: `chore(report): Phase 1 — bibliography 81 entries, notation agreed in macros.tex, literature notes`

---

## Phase 0/1 complete when

- All Tier 1–4 papers present in `literature/`
- `manuscript/STYLE_GUIDE.md`, `manuscript/style_intel/style_analysis.md`, `manuscript/style_intel/genre_conventions.md` all exist
- `references.bib` ≥ 55 entries
- `manuscript/macros.tex` has agreed notation set with literature-traceable comments
- Outline sign-off received from Barker

**Phase 2 (drafting) is a separate future session — do not begin writing report prose until Phase 0/1 are complete.**

---

## Phase 2 — Manuscript scaffolding (no prose drafting)

Phase 2 produces *scaffolding with content plans*: section/subsection headers
with semantic labels, named-paragraph topic stubs, equation placeholders with
labels, figure placeholders with caption stubs, and cross-references that
resolve. **No prose for any section** — each per-section drafting session is
a separate later effort. See plan archived as
`/home/vscode/.claude/plans/we-are-planning-to-federated-papert.md`.

### Bibliography hygiene additions

- [x] Add `hwangnoh2023graviton` (2310.04150)
- [x] Add `hwangnoh2024nonlinear` (2405.11786)
- [x] Add `higham2008functions` (Higham, *Functions of Matrices*, SIAM 2008)
- [x] Add `molervanloan2003nineteen` (Moler & Van Loan, SIAM Review 2003)
- [x] Verify `palessandro2024gertsenshtein` already present (it is)
- [x] Verify entry count ≥ 84 — **85 entries** total

### Figure infrastructure

- [x] Create `scripts/figures/` directory
- [x] Create `scripts/figures/README.md` (figure-script convention)
- [x] Stub `fig1_boccaletti_validation.py` (§4.1)
- [x] Stub `fig2_propagating_torsion_map.py` (§4.4 headline)
- [x] Stub `fig3_trace_ghost_growth.py` (§4.5 conditional)
- [x] Stub `figA1_architecture_tikz.tex` (App A diagram)
- [x] Stub `figC1_solver_convergence.py` (App C; cross-ref from App D)
- [x] Stub `figD1_stageA_heatmap.py` (App D)
- [x] Stub `figD2_stageB_sweep.py` (App D)
- [x] Stub `figD3_channel_comparison.py` (App D)

### Section scaffolding

- [x] Scaffold `manuscript/sections/abstract.tex` (4 placeholder slots)
- [x] Scaffold `manuscript/sections/introduction.tex` (4 BHL four-stage paragraphs)
- [x] Scaffold `manuscript/sections/theory.tex` (conventions + 4 subsections, equation labels)
- [x] Create `manuscript/sections/computational_approach.tex` (5 paragraphs, no subsections)
- [x] Update `manuscript/main.tex` to `\input` the new computational-approach file
- [x] Scaffold `manuscript/sections/results.tex` (5 subsections, 3 figure placeholders, ghost-condition equation label)
- [x] Scaffold `manuscript/sections/discussion.tex` (3 discussion subsections + Conclusion)
- [x] Scaffold `manuscript/sections/appendices.tex` (5 appendices A–E with subsections, named paragraphs, equation labels, figure placeholders)
- [x] Keep `manuscript/sections/acknowledgements.tex` unchanged (Barker boilerplate is correct as-is)

### Verification

- [ ] `latexmk -pdf manuscript/main.tex` compiles cleanly
- [ ] `\cref{}` resolution clean (no "undefined reference" warnings)
- [ ] `chktex` clean on all scaffolding files
- [ ] `\cite{}` resolution clean (no "Citation undefined" warnings)
- [ ] Forbidden-words grep on `manuscript/sections/*.tex` returns no matches
- [ ] `python -m py_compile scripts/figures/*.py` succeeds
- [ ] Compiled PDF table of contents matches the plan outline

### Commit

- [x] Commit Phase 2: scaffolding, figure stubs, bibliography hygiene (944b447)

---

## Phase 2.1 — Structural review and restructure (post initial scaffolding)

A three-agent structural review identified that the initial Phase 2 scaffolding tied §4 too closely to specific projected results, while several HPC campaigns are still in progress and the eventual headline finding cannot be confidently predicted. This phase implements the review's recommendations with a survey-shape principle that does not commit to specific findings.

### Style-guide refinements

- [x] Update `STYLE_GUIDE.md §(c)` to specify bare-CamelCase labels for all types (equations, sections, figures, appendices) — matches BHL corpus convention; remove `eq:`, `sec:`, `fig:`, `app:` prefixes. Cross-refs to labels in `docs/tex/` keep their existing prefixes.

### Section-level restructure

- [x] `introduction.tex`: rename paragraph titles to content-bearing micro-headings (no "Frontier tension"-style stage labels). Split Stage 4 into "In this report" + "Structure of this report" (PSALTer §I "Aims"+"Structure" split). Add bridging logic between the §1.2 ghost-restriction motivation and §4 sector coverage. Word allocation 800→900.
- [x] `theory.tex`: switch to bare CamelCase labels. Add Levi-Civita ε convention to Conventions paragraph. Upgrade hx↔ax framing from "structural observation" to "structural identity at linear order under the admissibility constraint of gr-qc/0001010". Fold §2.4 (Channel decomposition) into §2.3 (Gertsenshtein kernel and channel structure) — combined ~400 w subsection. Remove `\paragraph*{}` labels (BHL doesn't label them).
- [x] `computational_approach.tex`: shrink from 600 w to 400 w. Add new `Relation to other automation tools` paragraph distinguishing TIDAL from PSALTer/Hamilcar (per Agent 3 review I1).
- [x] `results.tex`: restructure around survey shape — Validation / Sector Survey (table-driven) / Channel Structure (mechanism explaining nulls) / Headline Finding (TBD slot for the strongest result emerging from in-progress campaigns). Add `Survey perimeter` opening paragraph. Add `Coverage rationale` paragraph (per Agent 3 review C4). Remove §4.5 ghost diagnosis subsection (moved to App D Validation Pathological). Word allocation 1700→1900 to absorb the survey table and headline slot.
- [x] `discussion.tex`: collapse the 4 subsections to flat named paragraphs (per Agent 2 review S3 — matches BHL letter convention). Connect path 1 explicitly to non-minimal couplings (Itin-Hehl, Karananas et al.). Use `bahamonde2025cubic` bibkey for path 3. Drop "or equivalent" hedges.
- [x] `appendices.tex`: switch to bare CamelCase labels. Move trace-channel ghost diagnosis to App D `Pathological cases` paragraph (was §4.5). Add new `Stability guard for plane-wave initial conditions` paragraph in App D documenting the engineering work that was previously invisible (per Agent 3 review I7). Rename App E equation-block subsections to be theory-agnostic (`SymbolicReference`, `SymbolicMechanism`, `SymbolicHeadline`) so they accommodate whichever theories prove most informative at drafting time.

### Word allocation (final)

- §1 Introduction: ~900 w (was 800)
- §2 Theory: ~1100 w (was 1200; folded §2.4 into §2.3)
- §3 Computational approach: ~400 w (was 600)
- §4 Results: ~1900 w (was 1700)
- §5+§6 Discussion+Conclusion: ~700 w (unchanged)
- Total: 5000 w main body

### Verification

- [x] `pdflatex` compiles cleanly to a 4-page PDF showing every planned heading.
- [x] All section/subsection headings appear in the compiled TOC; cross-refs resolve.
- [x] Forbidden-words grep on all section files returns no matches.
- [x] `python -m py_compile scripts/figures/*.py` succeeds.
- [x] No undefined references or citation warnings.

### Commit

- [x] Commit Phase 2.1 restructure (f6e0efd)

---

## Phase 2.2 — Appendix split and BHL-flavoured renaming

A check of the cached BHL/Barker corpus showed that appendix titles are uniformly content-bearing (compare PSALTer v2's "Construction of operators", "Moore--Penrose pseudoinversion", "Chequer-Hermicity"; Hamilcar's "Gauss--Codazzi and the Ricci scalar"; BHL long paper 2101.02645's "Irreducible decomposition of the fields", "Nonlinear Poisson brackets"). Generic stage labels like "Symbolic stage" / "Numerical stage" do not appear anywhere in the corpus. Splitting the single `appendices.tex` into per-appendix files is also a forward-looking decision given the anticipated content volume.

### Per-appendix file split

- [x] Create `manuscript/sections/appendices/` subdirectory
- [x] `appendices/architecture.tex` — App A "TIDAL pipeline architecture"
- [x] `appendices/symbolic.tex` — App B "Symbolic computation of the equations of motion" (was "Symbolic stage")
- [x] `appendices/numerical.tex` — App C "Numerical evolution of the linearised system" (was "Numerical stage")
- [x] `appendices/validation.tex` — App D "Validation against analytic limits and reference models" (was "Validation suite")
- [x] `appendices/symbolic_outputs.tex` — App E "Symbolic-pipeline outputs for the surveyed theories"
- [x] `appendices.tex` rewritten as a five-line wrapper that `\input`s each per-appendix file

### Label renames (cross-refs updated everywhere)

- [x] `\label{SymbolicStage}` → `\label{SymbolicComputation}` (App B)
- [x] `\label{NumericalStage}` → `\label{NumericalEvolution}` (App C)
- [x] `computational_approach.tex` cross-refs updated: `\cref{NumericalStage}` → `\cref{NumericalEvolution}`
- Other labels (`Architecture`, `Validation`, `SymbolicOutputs`, `ModalSolver`, `AnalyticalJacobian`, `OtherBackends`, etc.) unchanged.

### Verification

- [x] `pdflatex` compiles cleanly to a 4-page PDF.
- [x] Compiled PDF shows "Appendix A: TIDAL pipeline architecture", "Appendix B: Symbolic computation of the equations of motion", "Appendix C: Numerical evolution of the linearised system", "Appendix D: Validation against analytic limits and reference models", "Appendix E: Symbolic-pipeline outputs for the surveyed theories".
- [x] No undefined references; cross-refs from main body to appendices resolve.

### Commit

- [ ] Commit Phase 2.2 appendix split
