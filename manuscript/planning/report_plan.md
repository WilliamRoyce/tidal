# Plan: TIDAL MSci Project Report — Context-Building & Writing Strategy

## Context

The student (William Royce) is preparing the Cambridge **MSci project report** for **TIDAL** — a symbolic→numerical pipeline (Wolfram/xAct → Python PDE solvers) used to study the Gertsenshtein effect (gravitational↔EM wave conversion) in Poincaré gauge theory with torsion. Primary supervisor: **Dr Will Barker** (Charles University Prague, formerly Cambridge KICC); also Prof Mike Hobson and Prof Anthony Lasenby. The repo is `/workspaces/torsion-gertsenshtein` (branch `hpc/pgt-survey`).

**Hard format constraints**:
- **Main body: 5 000 words maximum.** This excludes the abstract, figures and captions, bibliography, and appendices.
- **Appendices: unlimited.** This is the lever that lets us preserve a real software-package contribution without spending main-body word budget on implementation detail.
- **Template**: Barker's supplied graduate-supervision report template (revtex-style scientific-paper format), intentionally chosen because the student plans to develop the report into a publication. So formatting/binding decisions are already locked in by the template.
- **Current manuscript state**: `manuscript/main.tex` and the `manuscript/sections/*.tex` files are *unchanged from Barker's template*. None of the existing prose, label conventions, or PACS-code examples reflect anything the student has written; they are placeholders to be overwritten. Treat them as scaffolding only.

**The report has dual character**: physics results (torsion-Gertsenshtein channel structure, ghost analysis, three constructive paths) AND a software-package announcement (TIDAL: 23 k LOC, 1 700+ tests). The central numerical contribution is the **Fourier modal solver** — the key solver innovation whose development, validation and deployment drove the campaign results. The other backends (IDA/CVODE/leapfrog/scipy) provide context and fallback but are less central now; App C focuses on the modal solver and the analytical Jacobian, with the others documented only briefly. With only 5 000 main-body words, the strategy is to keep the **main body focused on theory, results and discussion**, and push **TIDAL's full software documentation into appendices**, referenced from the main text. Implementation detail does not steal word budget from physics.

This plan describes (a) what context-gathering is still needed before writing, (b) the proposed report outline and writing order under the 5 000-word constraint, (c) the supporting artefacts (style analysis, figures plan, bibliography updates), and (d) the open questions the student must resolve before drafting begins.

---

## Key finding from this planning round

The three arXiv IDs the user flagged (**2506.02111**, **2406.09500**, **2512.25007**) are **all Will Barker's own software-package papers**:

| arXiv ID | Title | Year | Software | Authorship |
|----------|-------|------|----------|-----------|
| **2406.09500** | PSALTer: Particle Spectrum for Any Tensor Lagrangian | 2024 | PSALTer (Wolfram/xAct) | Barker, Marzo, Rigouzzo |
| **2506.02111** | The particle spectra of parity-violating theories (PSALTer v2) | 2025 | PSALTer upgrade | Barker, Karananas, Tu |
| **2512.25007** | Fast Poisson brackets and constraint algebras in canonical gravity | 2025 | **Hamilcar** (Wolfram/xAct) | **Barker, sole author** |

The user provided these as personal exemplars (not as a directive from Barker), but they are the closest-topic genre matches for TIDAL's software dimension and so dominate the §3 imitation list.

**The selection criterion for style imitation is topic-relevance to TIDAL, not authorship.** All three supervisors are senior, experienced academics whose writing is worth learning from — but their portfolios are heterogeneous. Hobson's CMB / Bayesian-cosmology work and Lasenby's geometric-algebra particle-physics work are excellent in their own right but use a different register from a torsion / gauge-gravity / Hamiltonian-analysis report. Pulling stylistic cues wholesale from off-topic papers will produce a register clash. The right rule is:

- **Mine heavily** any paper on PGT / torsion / gauge gravity / Hamiltonian analysis / ghost-tachyon analysis / gauge-gravity cosmology — regardless of whether Barker, Hobson, Lasenby, or any combination is on it.
- **Mine selectively** for narrow specialised patterns when the report calls for them: Hobson's Bayesian-evidence framing (e.g., 2102.11511) only if/where TIDAL's inference framework needs that language; Lasenby's GA notation conventions (e.g., 1912.05960) only if GA-style derivations enter the theory section.
- **Do not import** stylistic patterns from clearly off-topic supervisor work (pure CMB cosmology, geometric-algebra in octonions, etc.) — different topic, different register, different audience.

The Barker–Hobson–Lasenby joint corpus on PGT/torsion is the headline match because every one of those papers is *both* in the supervisory voice *and* on a topic adjacent to TIDAL. Their joint papers therefore dominate the §1, §2, §4, §5 imitation list. Solo papers (any author) earn a place on the imitation list if and only if they're on-topic.

### Confirmed Barker–Hobson–Lasenby joint corpus (highest-leverage style templates)

| arXiv ID | Year | Title | Authorship | Topic relevance |
|----------|------|-------|------------|-----------------|
| **2406.12826** | 2024 | Every Poincaré gauge theory is conformal | Barker + Hobson + Lasenby + Lin + Wei | **Direct** — vector torsion, conformal embedding |
| **2309.14783** | 2023 | Manifestly covariant variational principle for gauge theories of gravity | Barker + Hobson + Lasenby | **Direct** — Weyl/eWGT variational calculus |
| **2303.11094** | 2023 | Does gravitational confinement sustain flat galactic rotation curves? | Barker + Hobson + Lasenby | Adjacent — PN constraints on PGT |
| **2101.02645** | 2021 | Nonlinear Hamiltonian analysis of new quadratic torsion theories I | Barker + Lasenby + Hobson + Handley | **Direct** — Hamiltonian analysis, particle spectra |
| **2006.03581** | 2020 | Mapping Poincaré gauge cosmology to Horndeski theory for emergent dark energy | Barker + Lasenby + Hobson + Handley | **Direct** — PGT cosmology |
| **2003.02690** | 2020 | Addressing the H₀ tension with emergent dark radiation in unitary gravity | Barker + Lasenby + Hobson + Handley | Adjacent — PGT cosmology |
| **2008.09053** | 2021 | Fresh perspective on gauging the conformal group (Hobson + Lasenby joint, no Barker) | Hobson + Lasenby | Adjacent — conformal gauge theory |
| **2005.02228** | 2020 | Ghost and tachyon free Weyl gauge theories | Hobson + Lasenby + co-authors | **Direct** — ghost analysis in gauge gravity |

None of these eight joint papers are currently in `literature/`. **Fetching their TeX sources is the highest-priority context action** (Phase 0 below), alongside the three Barker package papers.

### Additional Barker papers (Phase 0 Tier 5 — recent work on metric-affine and PGT)

| arXiv ID | Year | Title | Topic relevance |
|----------|------|-------|-----------------|
| **2510.17094** | 2025 | Gertsenshtein effect on curved spacetime with geometric optics | **Direct** — on our exact topic |
| **2507.09228** | 2025 | Alleviating H₀ tension with Torsion Condensation (Legner, Handley, Barker) | **Direct** — propagating torsion cosmology |
| **2510.08201** | 2025 | Spectrum of pure R² gravity: full Hamiltonian analysis (Barker et al.) | **Direct** — Hamiltonian spectral analysis |
| **2506.21662** | 2025 | Infrared foundations I: catalogue of rank-three field theories | Adjacent — systematic field theory catalogue |
| **2507.05349** | 2025 | Infrared foundations II: all torsion-like theories, ghost-tachyon-free cases | **Direct** — torsion-like theories with ghost analysis |
| **2505.23894** | 2025 | Can metric-affine gravity be saved? (Barker, Marzo, Santoni) | Adjacent — ghost pathologies in metric-affine |
| **2402.07641** | 2024 | Particle spectra of Palatini/metric-affine theories (Barker, Marzo) | **Direct** — ghost/spectrum methodology |
| **2402.14917** | 2024 | Consistent particle physics in metric-affine gravity (Barker, Zell) | Adjacent — extended projective symmetry → ghost control |

### Key collaborator / field papers (Phase 0 Tier 6)

| arXiv ID | Year | Title | Topic relevance |
|----------|------|-------|-----------------|
| **2406.11956** | 2024 | Weyl-invariant Einstein-Cartan gravity (Karananas, Shaposhnikov, Zell) | Adjacent — conformal protection in EC/PGT |
| **2506.17017** | 2025 | Cosmology of Cubic PGT (Bahamonde, Iosifidis et al.) | **Direct** — cubic PGT cosmology, constructive path §5.2 |

### Papers already in `literature/` that must be studied (not yet analysed)

These were present before Phase 0 but not included in the style analysis. They are primary references for CONTENT (§1, §2, §4, §5) even though they are not supervisor papers:

| ID | Title / Authors | Section relevance |
|----|----------------|------------------|
| **1812.02675** | Nair, Park, Yoon — ghost-free PGT systematic classification (450 cases) | §1.2, §5.2 |
| **1804.05556** | Blagojević & Cvetković — General PGT Hamiltonian structure | §5.2 |
| **2301.02072** | Palessandro & Rothman — Simple derivation of Gertsenshtein effect | §1.1, §4 |
| **2310.04150** | Hwang & Noh — EM field definitions in graviton-photon conversion | §2.3 |
| **2004.02714** | Ejlli — Exact graviton-photon mixing in uniform B-field | §4 |
| **2405.11786** | Hwang & Noh — Graviton-photon in Euler-Heisenberg regime | §5 context |
| **2405.01407** | Palessandro — Gertsenshtein as probe of quantum gravity | §1, §5 |
| **2105.04565** | Caputo et al. — Dark Photon Limits: a Handbook | §4.2 |
| **2507.16609** | Domcke et al. — GW scattering on magnetic fields (2025) | §5 context |
| **hep-th_0103093** | Shapiro — Physical aspects of space-time torsion (comprehensive review) | §2 background |
| **gr-qc_0001010** | Hehl & Obukhov — How does EM couple to gravity? | §2 framework |
| **gr-qc_0307063** | Itin & Hehl — Maxwell nonminimally coupled to quadratic torsion | §2, §5 path 1 |
| **gr-qc_0305049** | Rubilar & Obukhov — Torsion nonminimally coupled to EM | §2 |
| **2510.08201** | Barker et al. — Spectrum of pure R² gravity (already local) | §5.2 Hamiltonian |

### Pre-arXiv foundational papers (bib entries only — no TeX available)

These must be in `references.bib` but cannot be downloaded as TeX sources. They are essential citations:

| Cite key | Reference | Section use |
|----------|-----------|-------------|
| `gertsenshtein1962` | M.E. Gertsenshtein, Sov. Phys. JETP **14**, 84 (1962) | §1.1 — original discovery |
| `boccaletti1970` | Boccaletti et al., Nuovo Cim. B (1970) | §4.1 — the formula we validate to 0.04% |
| `raffelt1988` | Raffelt & Stodolsky, PRD **37**, 1237 (1988) | §2.3, §4.2 — foundational photon mixing |
| `kibble1961` | T.W.B. Kibble, J. Math. Phys. **2**, 212 (1961) | §1.2 — gauge origin of PGT |
| `hehl1976` | Hehl, von der Heyde, Kerlick, Nester, Rev. Mod. Phys. **48**, 393 (1976) | §2.1 — definitive PGT review |
| `sezgin1980` | Sezgin & Van Nieuwenhuizen, PRD **21**, 3269 (1980) | §1.2 — first systematic ghost analysis |
| `cillis1996` | Cillis & Harari, PRD **54**, 4757 (1996) | §4 — cosmological conversion + plasma |

---

## What we already have (do not re-do)

- **Barker corpus inventory** (25 papers identified via arXiv author search): six local TeX sources + nineteen catalogued by abstract. Style synthesis (the user has confirmed this hits the target) covers his abstract pattern, intro scaffold ("modern frontier", 4-stage), section structure (long sections with subsections, not many short), equations sandwiched in prose with semantic labels (`\label{TorsionTransform}`), first-person plural with strategic hedging, footnote erudition, ironic/classical references, "Catch-22" framing for paradoxes. Solo-author papers (2206.00658, 2205.13534, 2311.11790, **2512.25007**) expose his unmediated voice most cleanly — these are the closest-fit imitation targets.
- **Software-paper genre survey** (16 exemplars): xAct, xPert, xTras, OGRe, FieldsX, Cadabra2, EinsteinPy, NRPy+, Dedalus, GRChombo, Einstein Toolkit, emcee, dynesty, PolyChord, Bilby, SciPy. Universal conventions extracted: name-first titles, problem-centric framing, dedicated architecture/validation/examples/conclusion sections, multi-method validation (analytic limits + reference codes + convergence + physical tests), parametric robustness over raw speed, modest first-release scope.
- **Existing material that may help (each with a clear judgement on its actual usefulness)**:
  - `docs/tex/*.tex` — architecture, gertsenshtein, modal_solver, solver_migration, derivation_performance, kinetic_matrix, multi_field_perturbation, background_fields, stability_probe, plane_wave_ic. **Reference material only** — these docs are maintained primarily as Claude context and are written in a continuous-update style that is not appropriate for direct transfer. Every sentence in the appendices must be freshly composed. Read the relevant supervisor and package papers to understand correct framing and what to include; then use the docs as a factual reference to check technical details (equations, parameter names, performance numbers). Never paste from them verbatim.
  - `CAMPAIGN.md`, `AMPLIFICATION_INVESTIGATION.md`, `docs/MEMORY.md` — research log with all numerical findings (Boccaletti 0.04 %; Stages A/B closed NULL; 276-run dark-photon sweep NULL; hx↔ax torsion-independence; |δ₁|<0.005 stability window; trace-channel ghost diagnosis). **Primary source for §4 results.**
  - `phd-application-context.md` — **NOT a seed for this report.** Was written for PhD CV updates, target audience is admissions committees, content is not technically calibrated for an MSci report.
  - `docs/talks/practice_2026_04/script.md` — **not a content source.** Authored for a non-expert audience; it lacks the technical depth and focus that the report requires. For motivating the introduction, study how the supervisor's papers open instead.
- **Manuscript skeleton**: `manuscript/main.tex` (revtex-based), `manuscript/sections/{abstract,introduction,theory,results,discussion,acknowledgements,appendices}.tex`, `manuscript/macros.tex`, `manuscript/references.bib` (43 entries — strong on Gertsenshtein/PGT/numerics, missing the three Barker software papers and a handful of dark-photon and ghost-instability refs), `manuscript/figures/` (empty except `.gitkeep`). **All `.tex` files are unmodified Barker-template content** — no original writing yet. We will overwrite them; existing prose, label conventions, PACS-code examples, and section names are placeholders, not constraints.

---

## The plan, by phase

### Phase 0 — Fetch missing references and freeze style intel  *(1 day)*

**Pre-execution: archive this plan and create a checklist in the repo before any other action.**
- Copy this plan to `manuscript/planning/report_plan.md` in the repo so it survives container rebuilds and session losses.
- Create `manuscript/planning/report_checklist.md` — a flat, actionable checklist version of Phase 0 and Phase 1 items only (one checkbox per concrete action: download paper X, create file Y, add bib entry Z). This becomes the live tracking document for Phase 0/1 execution and replaces re-reading this plan during each work session.
- Update Claude memory (`MEMORY.md`) with a pointer to both files.
- Commit: `chore(report): add report plan archive and Phase 0/1 checklist`.

1. **Tier 1 — Barker software-package papers** (highest priority for §3 of the report). Download TeX sources to `literature/`:
   - `literature/2406.09500/` — PSALTer original (47 pp, Wolfram/xAct, github.com/wevbarker/PSALTer)
   - `literature/2506.02111/` — PSALTer v2 (58 pp, parity-violating extension)
   - `literature/2512.25007/` — Hamilcar (44 pp, **sole author** — closest-fit voice model)
2. **Tier 2 — Barker–Hobson–Lasenby joint corpus on PGT/torsion** (highest priority for §1, §2, §4, §5 — same topic family as TIDAL):
   - `literature/2406.12826/` — *Every Poincaré gauge theory is conformal* (already local — verify TeX completeness)
   - `literature/2309.14783/` — *Manifestly covariant variational principle for gauge theories of gravity*
   - `literature/2303.11094/` — *Does gravitational confinement sustain flat galactic rotation curves?*
   - `literature/2101.02645/` — *Nonlinear Hamiltonian analysis of new quadratic torsion theories Part I*
   - `literature/2006.03581/` — *Mapping Poincaré gauge cosmology to Horndeski theory*
   - `literature/2003.02690/` — *Addressing H₀ tension with emergent dark radiation in unitary gravity*
3. **Tier 3 — Hobson + Lasenby joint (no Barker), close-topic** (supporting):
   - `literature/2008.09053/` — *Fresh perspective on gauging the conformal group*
   - `literature/2005.02228/` — *Ghost and tachyon free Weyl gauge theories*
4. **Tier 4 — additional Barker solo voice**:
   - `literature/2205.13534/` — solo-author multipliers paper (raw voice, 2022)
   - `literature/2311.11790/` — solo-author Yang-Mills gravity (3 pp brief)
5. **Skip wholesale stylistic imitation, but do not skip reading**: solo-Hobson CMB / Bayesian-cosmology papers, solo-Lasenby geometric-algebra-in-particle-physics work. They are excellent papers worth reading for general academic-writing practice, but their topic divergence means imitating their register would produce a register clash. Pull narrow patterns from them (Hobson's prior-sensitivity language; Lasenby's GA notation) only where the report's content specifically calls for those patterns.
6. **Distil and save style understanding (analytical, not extractive) into the repo**:
   - **Goal**: build *our own* understanding of how Barker / the Cambridge gauge-gravity group writes — abstract structure, intro scaffold, section organisation, equation handling, paragraph rhythm, transition devices, voice — so we can write *original prose in the same register*. The point is internalisation and emulation, **not** to maintain a pile of quotes for copy-paste reuse. Plagiarising the supervisors' own sentences would be self-defeating.
   - `manuscript/STYLE_GUIDE.md` — concrete checklists derived from the analysis (section opening rule, equation rule, label rule, voice rule, figures rule, self-review checklist). See "Style guide artefact" section below for contents.
   - `manuscript/style_intel/style_analysis.md` — analytical write-up: for each pattern (e.g., "the four-stage Barker introduction"), describe the pattern in our own words, with one or two short illustrative excerpts (≤ 1 sentence each, attributed) only where needed to anchor the analysis. The document's purpose is to teach the student the patterns; once the patterns are internalised, the document is no longer consulted line-by-line during writing.
   - `manuscript/style_intel/genre_conventions.md` — analytical synthesis of physics-software-paper conventions across the 16-paper genre sample, in our own words. Used to plan the appendices on TIDAL.
7. Commit as one prep commit: `chore(report): freeze style analysis and supervisor reference corpus for MSci draft`.

### Phase 1 — Outline approval, notation set, bibliography hygiene  *(0.5 day)*

1. **Get supervisor sign-off on the outline** (Section "Proposed report outline" below) before drafting any section. With only 5 000 main-body words, every word matters — outline drift mid-draft is expensive.
2. **Bibliography**: extend `manuscript/references.bib` with:
   - The three Barker software papers (`barker2024psalter`, `barker2025psalter2`, `barker2025hamilcar`).
   - Software comparators cited in App A (NRPy+, Dedalus, OGRe, FieldsX, Cadabra2, xAct/xPert if not present).
   - Dark-photon / hidden-photon review (e.g. Caputo–Millar–O'Hare–Vitagliano 2021 *"Dark Photon Limits: a Cookbook"*).
   - Ghost-instability literature in PGT (Sezgin–Van Nieuwenhuizen, Beltrán-Jiménez–Heisenberg–Koivisto, Nair–Park–Yoon).
   - GW-detection context anchors (LISA, NANOGrav, LIGO O4 — one each).
   - **Target: 55–60 entries total** after additions.
3. **Agree the notation set in `manuscript/macros.tex` before drafting begins.** *What this means in practice*: list every symbol the report will use (metric signature, indices, torsion components, contortion, the irreducible decomposition, Gertsenshtein constants, plasma quantities) and define a macro for each. Once that set is agreed, do not introduce new macros mid-draft — because adding `\Tirr` halfway through means going back to retrofit earlier sections that wrote out the long form. The "freeze" is about not paying that retrofit cost twice. New macros are still allowed during drafting if genuinely needed; the rule is just "add it and propagate to all already-drafted sections in the same commit", not "never add another macro".
4. **Survey the notation conventions in the BHL joint corpus and Barker's solo papers.** The goal is not to invent symbols but to match what the relevant literature uses — the same discipline applied to parameter names in `research/`. Read through the Tier 1–2 papers (Phase 0) with a focus on: metric signature, vierbein/spin-connection notation, torsion irreducible decomposition labels, Gertsenshtein coupling symbols. Record the standard choices in `manuscript/macros.tex` comments so the macro set is explicitly traceable to its source. Where multiple conventions exist in the literature, pick the one that appears most frequently in the papers being cited.
5. **Notation alignment with `docs/tex/preamble.tex`**: the technical docs and the report should agree on symbol choices. Where they currently disagree, the report notation (sourced from the literature survey above) wins; update the docs to match.
6. Commit: `chore(report): outline approved, notation agreed, bibliography extended`.

### Phase 2 — Drafting  *(dedicated future session — NOT in scope of this plan's execution)*

This phase is reserved for a dedicated future writing session with its own plan. The current plan ends at Phase 1. No writing of report prose happens until Phase 0 and Phase 1 are complete and the student has internalised the style analysis.

**Planned writing order** (to be followed in the future session): Appendices A–E (TIDAL detail) → §2 Theory → §4 Results → §3 Computational approach → §1 Introduction → §5 Discussion → §6 Conclusion → Abstract.

Justification for the order: Appendices A–E are the most factual sections and their content is best understood after having read the Tier 1 software-paper exemplars (Phase 0). Writing them first builds momentum, pins notation, and makes every later main-body cross-reference concrete. `docs/tex/` files serve as a factual reference for technical details — not prose to be copied. Every sentence must be freshly written after consulting the relevant literature to understand correct framing. §2 describes the physics; §4 reuses §2 labels; §3 is a 600-word précis with cross-references; §1 last-but-one because the hook depends on §4's final shape; Abstract last (≤ 180 words).

Draft slightly long, then cut: aim for 5 500–6 000 main-body words in the first draft, then trim to 5 000. Cutting is easier than padding and improves clarity.

### Phase 3 — Figures  *(parallel with Phase 2 §4 / App D)*

See "Figures plan" below. Main body gets ~3 figures (validation, headline result, ghost diagnosis); appendices absorb the rest (architecture diagram, sweep heatmaps, performance benchmarks, channel comparisons). **Existing figures already produced** (architecture diagrams, algorithm flowcharts, channel comparisons) are available and should be reviewed — but the report is not designed around them. If an existing figure serves the purpose at publication quality it is used; if not, it is regenerated or replaced. The report design determines what figures are needed, not the reverse.

**Reproducibility requirement**: every figure that appears in the final report must be generated by a dedicated, self-contained script under `scripts/figures/` (e.g. `scripts/figures/fig1_boccaletti_validation.py`). Each script reads campaign data, produces a PDF, and does nothing else. Scripts must be committed alongside the figures so that any figure can be regenerated from raw data by a reader. This is in addition to any existing `tidal plot` convenience wrappers. All figures land in `manuscript/figures/`. The architecture diagram is TikZ — `Notes.svgz` and `docs/TIDAL_Logo_TikZ_Figure.svg` show TikZ already in use.

### Phase 4 — Supervisor iteration  *(weave across Phases 2–3, ~6 weeks total)*

| Week | Send to | Content | Form |
|------|---------|---------|------|
| 1 | Barker only | App A + App C drafts | Compilable `.tex` on a `feedback/wk1` branch |
| 2 | Barker + Hobson | App B + App D + App E + §2 (theory) | Same |
| 3 | Barker + Lasenby | §4 results + App F/G (Lasenby for ghost-channel framing) | Same |
| 4 | All three | Full draft incl. §1, §3, §5, §6, abstract | Compiled PDF + `.tex` |
| 5 | All three | Revisions, figures finalised, word-count audit on main body | Same |
| 6 | — | Final polish, bibliography sweep, submission | Final PDF |

Pre-agree a 7-day turnaround norm in week 1; if missed, push without blocking.

### Phase 5 — Software release alongside the report  *(submission day)*

If supervisors approve (open question 5 below): tag a Zenodo DOI for the repo at submission, update `CITATION.cff`, cut a GitHub release `v1.0-msci`. Appendix H records exact commit hash + TIDAL version + HPC campaign IDs (`hpc_results/2873…`) keyed to each figure.

---

## Proposed report outline

**Strategy under the 5 000-word constraint**: keep the main body **physics-first** (introduction, theory, results, discussion) and **push almost all software detail into appendices**. Reference the appendices liberally from the main text — readers who care about implementation can follow the cross-references; word budget stays free for physics. This mirrors how Barker's HiGGS paper (2206.00658) handles deep technical material via appendices while keeping the main narrative driven by the physics result.

**Body word allocation: 5 000 words main + abstract + appendices (unlimited).** Targets below sum to 5 000; tolerate ±5 % per section in drafting and rebalance at the end.

**Template-priority hierarchy** (selector is topic-relevance to TIDAL, then authorship/voice match):

1. **For App A–E (TIDAL pipeline detail) — Barker's three software papers**: PSALTer (2406.09500) for capability/scope and "Examples with Code" structure; PSALTer v2 (2506.02111) for evolution / extension presentation; Hamilcar (2512.25007) for sole-author voice and validation case-study (`§2.3 pure GR at two loops` is the gold standard).
2. **For §1, §2, §3, §4, §5 (main body) — Barker–Hobson–Lasenby joint corpus on PGT/torsion**: 2406.12826, 2309.14783, 2101.02645, 2006.03581, 2510.08201, 2407.09598, 2206.00658.
3. **For supplementary on-topic patterns**: Hobson + Lasenby joint on gauge gravity (2008.09053, 2005.02228); solo-Lasenby on-topic gauge-gravity work (gr-qc/0405033 Doran-Lasenby-Gull) if GA derivations enter §2.
4. **For specialised narrow patterns only** — pulled from off-topic but methodologically-related supervisor work where the content specifically calls for that pattern. Do not import register wholesale.

### Main body (5 000 words)

| # | Section | Words | Primary content seed | Imitation template |
|---|---------|-------|-----------------------|--------------------|
| — | Abstract (separate, ≤ 250 words) | ~180 | written last | 2206.00658 abstract; 2406.12826 abstract |
| 1 | Introduction | 800 | supervisor PGT intros (2406.12826, 2101.02645) for shaping + relevant literature | 2206.00658 + 2510.08201 + 2406.12826 introductions |
| 1.1 | Gertsenshtein & the GW–EM bridge | 250 | relevant literature; docs/tex as factual check only | — |
| 1.2 | Why torsion: PGT, ghosts, the open question | 300 | own synthesis from `literature_critical_analysis.md` and the BHL joint corpus | 2406.12826 intro paradox; 2101.02645 §I |
| 1.3 | This work and its contribution | 250 | — | 2510.08201 §I closing |
| 2 | Theoretical framework | 1 200 | BHL joint corpus for physics content and framing; docs/tex as factual check | 2206.00658 §II, 2510.08201 §II–III, 2406.12826 §II |
| 2.1 | PGT, contortion, irreducible torsion | 350 | 2406.12826 §II and adjacent BHL papers | 2406.12826 §II |
| 2.2 | Linearisation around flat + B₀ background | 350 | 2510.08201 §II and related | 2510.08201 §II |
| 2.3 | The Gertsenshtein conversion kernel | 250 | literature; docs/tex as factual check | — |
| 2.4 | Channel decomposition: hx↔ax vs trace sector | 250 | own analysis (referenced to App F) | — |
| 3 | Computational approach (overview only) | 600 | docs/tex/architecture.tex (facts only — prose freshly written) | 2206.00658 §III opening; 2406.09500 §I-B |
|   | Single section, no subsections. States the symbolic→numerical pipeline at architecture level; leads with the Fourier modal solver as the key numerical innovation (auto-selection logic, eigendecomposition, machine-precision accuracy), mentions the other backends briefly as fallbacks, declares the validation strategy and the HPC scale, then defers all detail to App A–E with explicit `\cref{app:tidal-architecture}` etc. references. | | | |
| 4 | Results | 1 700 | `CAMPAIGN.md`, `AMPLIFICATION_INVESTIGATION.md`, `hpc_results/.../`, regression fixtures | 2406.12826 §V; 2407.09598 results; 2510.08201 §IV |
| 4.1 | Einstein–Maxwell baseline (validation, Boccaletti 0.04 %) | 250 | `examples/einstein_maxwell_baseline/` | — |
| 4.2 | Stage A null: T1 dark-photon-plasma | 350 | `CAMPAIGN.md` Stage A | — |
| 4.3 | Stage B null: T2 Einstein–Cartan, 276-run sweep | 350 | `AMPLIFICATION_INVESTIGATION.md` | — |
| 4.4 | Propagating PGT: hx↔ax torsion-independence + |δ₁| stability window | 400 | `hpc_results/.../d22_shapiro_sup_intr/` | — |
| 4.5 | Trace-channel ghost diagnosis | 350 | same parameter point, longer integration | — |
| 5 | Discussion | 600 | own synthesis | 2206.00658 §V, 2510.08201 §VI, 2406.12826 §VI |
| 5.1 | What the null constrains, what hx↔ax independence reveals | 250 | — | — |
| 5.2 | Three constructive paths beyond torsion-independence | 250 | `literature_critical_analysis.md` | 2406.12826 §VI |
| 5.3 | Limitations and threats to validity | 100 | — | — |
| 6 | Conclusion | 100 | — | 2510.08201 closing paragraph |

Sum: **5 000 words** main body.

### Appendices (unlimited; this is where TIDAL lives)

| Label | Title | Approx. words | Source |
|-------|-------|---------------|--------|
| **A** | TIDAL pipeline architecture | ~1 800 | Barker software papers for framing; docs/tex as factual reference |
| **B** | Symbolic stage: xAct, xPert, Euler-Lagrange in Wolfram | ~1 600 | Barker software papers; docs/tex as factual reference |
| **C** | Numerical stage: Fourier modal solver (primary), analytical Jacobian, and solver backend overview | ~2 000 | Barker software papers + NRPy+/Dedalus genre conventions; docs/tex as factual reference |
| **D** | Validation suite: analytic limits, convergence, ghost-stability self-consistency | ~1 200 | docs/tex/stability_probe.tex, background_validity.tex, regression-test fixtures |
| **E** | HPC campaign infrastructure and reproducibility | ~1 000 | `CAMPAIGN.md`, `hpc_results/campaigns/` |
| **F** | Theoretical derivation: hx↔ax torsion-independence at quadratic order | ~600 | own derivation |
| **G** | Trace-channel ghost analysis | ~400 | own derivation + numerical diagnostic |
| **H** | Software release manifest (Zenodo DOI, commit hash, campaign IDs keyed to figures) | ~200 | `CITATION.cff`, repo HEAD |

Appendices total: ~8 800 words. Appendices may be longer if needed; the constraint is reader patience, not the regulations.

**Cross-referencing pattern in the main text**: at every point where a main-body sentence would otherwise need to spend words on implementation detail, drop a forward reference instead. Examples:
- *"We solve the linearised system using the Fourier modal backend (App C.3), which auto-selects on flat metrics with periodic boundary conditions."* — one sentence, defers everything.
- *"The pipeline reproduces Boccaletti's analytic conversion kernel to 0.04 % at the validation point (App D.1, Fig. 2)."* — saves the convergence-study real estate for the appendix.
- *"All numerical results below come from runs on CSD3 sapphire/icelake (App E)."* — a single declarative reference covers all of §4.

---

## Style guide artefact

Create `manuscript/STYLE_GUIDE.md` with concrete checklists derived from the style analysis. The checklists encode patterns Barker (and the BHL joint corpus) actually uses; they do not contain his sentences.

- **Section opening rule**: every (sub)section starts with one *what* sentence + one *why* sentence before any equation, citation, or list. (Pattern observed across the BHL corpus.)
- **Equation rule**: every numbered equation has prose immediately before AND after; no orphan displays. Punctuate equations as part of the sentence (revtex default).
- **Label rule**: semantic labels only — `\label{TorsionTransform}`-style descriptive CamelCase, never `\label{eq:1}`. For section labels, follow whatever convention the BHL joint papers actually use — verify by grepping the local TeX sources in Phase 0 (the existing `sec:Theory` in our current `sections/*.tex` is template placeholder, not Barker's actual practice). The point is consistency with his published convention, not with the template scaffolding.
- **Voice rule**: first-person plural ("we"), active for our own work; passive sparingly, mainly for problem diagnosis. Strategic hedging ("we conjecture", "appears to") allowed where uncertain.
- **Footnote rule**: footnotes for substantive asides only; citations always inline with `\cite{}`. Ironic / classical references allowed but rare — capped at one per major section.
- **Numbering**: equations numbered by section (revtex default).
- **Citation density**: introduction ≥ 1 cite per paragraph; theory/results ≥ 1 cite per subsection; cite every named tool on first mention.
- **Forbidden words list**: "novel" (let the result speak), "in this paper" (use "we" + verb), "methods" as a section title (prefer "computational approach" or similar).
- **Figures rule**: every caption is two sentences — what the panel shows + what to take from it.
- **Self-review checklist**: macros consistent, all `\cref{}` resolves, no TODO comments, line numbering off for shared drafts, hyperref colours matched to the revtex template.

The companion `manuscript/style_intel/style_analysis.md` (described in Phase 0) is an *analytical* document — it describes the patterns in our own words, with brief illustrative excerpts only where needed to anchor a point. It is not a quote bank for reuse. The output of using the style guide is original prose in the same register.

---

## Figures plan

Figures don't count against the 5 000-word cap, so use them as load-bearing communication where prose would have to spend words. Place a *headline* figure in each main-body section that has one; route detail figures into the appendices.

### Main-body figures (3)

| # | Title | Section | Data source | Tooling |
|---|-------|---------|-------------|---------|
| 1 | Boccaletti validation: P_conv vs t, sim vs analytic, with residual panel | §4.1 | `examples/einstein_maxwell_baseline/` regression fixtures | dedicated `scripts/figures/fig1_boccaletti_validation.py` |
| 2 | Headline result — propagating-torsion parameter map: A vs δ₁ stability window, with hx↔ax channel torsion-independence visible at A = 1 | §4.4 | `hpc_results/28739692/d22_shapiro_sup_intr/` | dedicated `scripts/figures/fig2_propagating_torsion_map.py` |
| 3 | Trace-channel ghost growth: log\|amplitude\| vs t, log scale | §4.5 | longer integration of the same parameter point | dedicated `scripts/figures/fig3_trace_ghost_growth.py` |

### Appendix figures (5)

| # | Title | Appendix | Data source | Tooling |
|---|-------|----------|-------------|---------|
| A.1 | TIDAL architecture diagram | App A | schematic | TikZ |
| C.1 | Solver performance: wall-clock per backend, ± analytical Jacobian | App C | benchmark data | dedicated script |
| D.1 | Stage A heatmap: P_max in (m_A, ε) plane | App D | `hpc_results/campaigns/d20_bahamonde_amp_rerun_v2/` | dedicated script |
| D.2 | Stage B sweep: 276-run scatter, P_max vs sin²(κB₀t/2) | App D | `hpc_results/campaigns/d21_barker_sup/` | dedicated script |
| D.3 | Channel comparison: hx↔ax vs trace, time series | App D | same as headline figure | dedicated script |

**All figures must have a dedicated, self-contained script in `scripts/figures/`** — each reads raw campaign data and writes a PDF only. This is mandatory for publication reproducibility; it is not optional even when `tidal plot` produces the same output. Existing figures (architecture diagrams, flowcharts) are reviewed during Phase 3 and either adopted directly or regenerated to match the report's needs.

---

## Risk register

| Risk | Mitigation |
|------|------------|
| **5 000-word cap pressure**: physics + theory + results + discussion is genuinely a lot to fit in 5 000 words. | Aggressive use of appendices for software detail (App A–E absorb everything implementation-related). Cross-reference, don't re-explain. Draft slightly long, then cut. The cut is itself an editorial discipline that improves clarity. |
| **Null-result framing**: Stages A/B closed NULL; dark-photon NULL; propagating-torsion A = 1.0. Reads as failure unless reframed. | Frame nulls as **constraints on theory space** (Barker's framing in 2406.12826: "we exclude X"). Lead §4 with the Boccaletti positive validation; present each null as "we constrain channel X to amplification < Y". The hx↔ax torsion-independence is the **central positive analytic finding** — make it the headline, not the nulls. |
| **Appendix bloat eclipsing main body**: with unlimited appendices it's tempting to keep adding. | Each appendix has a stated purpose; if a paragraph doesn't serve the purpose, cut. Appendices total ~8 000–10 000 words is plenty. |
| **Scope creep into PhD territory**: the three constructive paths are tempting to half-derive. | §5.2 lists them and costs them; **no derivations in main body**. Appendix F sketches one path's first step only if needed. |
| **Supervisor feedback latency** (three supervisors). | Send appendix chunks to Barker only first (fast, low-stakes — the technical-detail material is easiest to review). Reserve full-draft circulation for week 4. Pre-agree a 7-day turnaround in week 1. |
| **Reproducibility audit at examination**. | Tag a Zenodo DOI at submission. Appendix H records exact commit hash, TIDAL version, HPC campaign IDs keyed to each figure. |
| **Notation drift** between `docs/tex/*.tex` and the report. | Agree the macro set in `manuscript/macros.tex` before Phase 2; whenever a new macro is needed during drafting, propagate to all already-drafted sections in the same commit (and update `docs/tex/preamble.tex` so the technical docs stay aligned). |
| **Plagiarism risk from heavy supervisor mining**. | The style files are *analytical* (patterns described in our own words), not extractive. Every sentence in the report is original prose. Imitate the *register*, not the wording. |

---

## Open questions for the student to resolve before Phase 2

1. **Submission deadline** — drives the drafting schedule's start date.
2. **Title**: current placeholder is template default. Lock now — title shapes §1's framing. Working suggestion: *"Torsion in the Gertsenshtein effect: a symbolic-numerical survey of Poincaré gauge theory"*. Alt: *"TIDAL: a symbolic-numerical pipeline for the Gertsenshtein effect in Poincaré gauge theory"* (software-name-first, follows Barker's PSALTer/Hamilcar convention).
3. **Author affiliation block**: Cavendish Astrophysics, University of Cambridge. *(Closed.)*
4. **Supervisor draft preference**: incremental (plan above assumes this — appendices first, main body in chunks) or single full-draft submission?
5. **Software release at submission**: does Barker want TIDAL released alongside the report (Zenodo DOI, GitHub release `v1.0-msci`, `CITATION.cff` updated)? The infrastructure exists; decision is yes/no plus timing (release with submission vs after examiners report).
6. **arXiv classification**: at submission time, replace template's PACS codes with current arXiv categories (likely `gr-qc + hep-th + astro-ph.CO`). Not a planning question; just an execution detail.
7. **Examiner**: known but the report is written independently of the examiner — no stylistic slanting. *(Closed.)*
8. **Verbatim reuse from `docs/tex/*.tex`**: none. The docs were written without the benefit of the style research we are now doing and are primarily Claude context, not report prose. All report sentences are freshly composed. *(Closed.)*

---

## Critical files (paths to remember during execution)

- `/workspaces/torsion-gertsenshtein/manuscript/planning/report_plan.md` — THIS FILE (repo-local archive)
- `/workspaces/torsion-gertsenshtein/manuscript/planning/report_checklist.md` — live Phase 0/1 checklist
- `/workspaces/torsion-gertsenshtein/manuscript/main.tex` — document root (Barker template)
- `/workspaces/torsion-gertsenshtein/manuscript/sections/{abstract,introduction,theory,results,discussion,acknowledgements,appendices}.tex` — section files; all currently template placeholders, to be overwritten with original prose
- `/workspaces/torsion-gertsenshtein/manuscript/macros.tex` — notation set; agree contents in Phase 1
- `/workspaces/torsion-gertsenshtein/manuscript/references.bib` — 43 entries; extend in Phase 1
- `/workspaces/torsion-gertsenshtein/manuscript/figures/` — empty; populate in Phase 3
- `/workspaces/torsion-gertsenshtein/manuscript/STYLE_GUIDE.md` — to be created in Phase 0
- `/workspaces/torsion-gertsenshtein/manuscript/style_intel/style_analysis.md` — analytical write-up of supervisor style patterns; to be created in Phase 0
- `/workspaces/torsion-gertsenshtein/manuscript/style_intel/genre_conventions.md` — analytical write-up of physics-software-paper genre conventions; to be created in Phase 0
- `/workspaces/torsion-gertsenshtein/docs/tex/{architecture,gertsenshtein,modal_solver,solver_migration,kinetic_matrix,multi_field_perturbation,background_fields,stability_probe,derivation_performance,solver_optimizations}.tex` — factual reference material (not prose sources)
- `/workspaces/torsion-gertsenshtein/CAMPAIGN.md`, `AMPLIFICATION_INVESTIGATION.md`, `docs/MEMORY.md` — research log with all numerical findings
- `/workspaces/torsion-gertsenshtein/literature/{2406.09500,2506.02111,2512.25007,2309.14783,2303.11094,2101.02645,2006.03581,2003.02690,2008.09053,2005.02228}/` — to be downloaded in Phase 0
- `/workspaces/torsion-gertsenshtein/hpc_results/campaigns/` — figure data sources

---

## Verification (how to know the plan is working)

- **End of Phase 0**: all Tier 1–4 papers in `literature/`; `manuscript/STYLE_GUIDE.md` and `manuscript/style_intel/style_analysis.md` exist and the student can describe Barker's intro scaffold from memory.
- **End of Phase 1**: outline approved by Barker (email or chat); `references.bib` ≥ 55 entries; `manuscript/macros.tex` has an agreed notation set committed.
- **End of each section draft**: `latexmk -pdf main.tex` compiles cleanly; `\cref{}` resolves; no `TODO` comments; main-body word count tracked against the 5 000 cap; runs through `chktex` without warnings.
- **End of Phase 4 week 4**: full draft compiles to PDF; main body within 5 000 words; appendices coherent and cross-referenced; circulated to all three supervisors.
- **Submission day**: Appendix H contains commit hash matching repo HEAD; Zenodo DOI minted (if open question 5 = yes); `CITATION.cff` updated; main-body word count verified ≤ 5 000.
