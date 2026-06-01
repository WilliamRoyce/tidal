# Title — drafting spec (TIDAL MSci thesis)

Binding stylistic and structural specification for the thesis Title. Distilled from a `\title{}` survey of 38 papers in `literature/` (Barker solo/first-author, BHL joint, Hobson+Lasenby joint, Gertsenshtein primary literature, wider PGT/non-minimal-EM-torsion), the existing [abstract_drafting_spec.md](abstract_drafting_spec.md) and [conclusion_drafting_spec.md](conclusion_drafting_spec.md), the supervisor-voice analysis in [style_analysis.md](../style_intel/style_analysis.md), and the user's framing direction this session (physics-led headline; TIDAL is not the leading noun phrase; this is an MSci physics report imitating real journal-paper style, not a software-announcement). This file is the brainstorm session's reference for length, structure, terminology, punctuation, capitalisation, and the placement of the package name. It does **not** propose final titles — that is the next session.

The two working titles already in [report_plan.md](report_plan.md) §"Open questions" #2 are the first inputs the next session must evaluate against the rules below.

## §1 Binding rule (set by the user this session)

The title leads with the **physics finding** — the coupling-level Gertsenshtein survey across the quadratic Poincaré gauge theory operator space — and not with the package TIDAL. The MSci report imitates the *style* of journal papers in the Barker corpus; it does **not** imitate the *genre* of a software-announcement paper. The package name may appear as a trailing qualifier or be omitted entirely; it must not be the title's leading noun phrase. This rule supersedes the second working title in `report_plan.md` ("TIDAL: a symbolic-numerical pipeline…"), which is non-preferred for that reason.

The headline must also be consistent with the §4 framing in [abstract_drafting_spec.md](abstract_drafting_spec.md): a coupling-level Bayes-factor survey of the quadratic-PGT operator space, with named couplings ($\delta_1$, $\chi_1$, …) admitting or forbidding ghost-free Gertsenshtein amplification. The title does **not** quote plane-wave amplification magnitudes and does **not** flag what the work does not cover (no localised-geometry promise, no follow-up deferment).

## §2 Corpus surveyed

38 papers, all local. Title strings extracted directly from `\title{}` in the TeX sources; multi-line/`\boldmath`/`\Large` formatting stripped to the printed string. Barker solo/first-author papers (Tier 1, weight 3): 2206.00658, 2205.13534, 2311.11790, 2406.09500, 2506.02111, 2512.25007, 2402.07641, 2402.14917, 2505.23894, 2510.08201, 2407.09598, 2506.21662, 2507.05349. BHL joint corpus on PGT/torsion (Tier 2, weight 2): 2406.12826, 2309.14783, 2303.11094, 2101.02645, 2006.03581, 2003.02690. Hobson+Lasenby joint (Tier 3, weight 1): 2008.09053, 2005.02228. Gertsenshtein primary literature (Tier 4, weight 1): 2301.02072, 2310.04150, 2510.17094 (Tomomatsu, Suyama, Gondolo — **not Barker**), 2406.17853, 2507.16609, 2004.02714, 2405.01407, 2405.11786, 2405.08865. Wider PGT / non-minimal EM-torsion / genre cross-reference (Tier 5, weight 1): 1812.02675, 1804.05556, 2506.17017, 2410.01355, gr-qc/0001010, gr-qc/0307063, gr-qc/0305049, hep-th/0103093, 2105.04565, 2406.11956.

## §3 Length

Word-count distribution across the corpus: minimum 6 (`Can metric-affine gravity be saved?`), maximum 22 (`Systematic study of background cosmology in unitary Poincaré gauge theories with application to emergent dark radiation and $H_0$ tension`, 2003.02690), median ~10, modal band **8–12 words**. Barker solo/first-author median ~10 (e.g. `Geometric multipliers and partial teleparallelism in Poincaré gauge theory` is 9; `Spectrum of pure R² gravity: full Hamiltonian analysis` is 8; `Every Poincaré gauge theory is conformal: a compelling case for dynamical vector torsion` is 12). **Target for the thesis: 9–14 words.** Below 8 the title cannot carry both the operator-space scope and the coupling-level finding; above 14 it starts to read as the abstract's first sentence.

## §4 Structural archetypes observed

Each archetype is named, illustrated with corpus IDs (one example each, attributed), and graded for fitness against §1.

1. **Imperative-stance declarative** — a position is stated as the title (`Every Poincaré gauge theory is conformal: a compelling case for dynamical vector torsion`, 2406.12826; `Does gravitational confinement sustain flat galactic rotation curves without dark matter?`, 2303.11094; `Can metric-affine gravity be saved?`, 2505.23894). **Best fit** for a coupling-level structural finding. The colon-qualifier (`a compelling case for…`) is the Barker-signature subtitle device — it narrows the headline into a specific claim. Question form is admissible only when the body actually answers it.
2. **Topic-then-Hamiltonian / classification register** — `Topic: full Hamiltonian analysis` / `Topic: a systematic approach` (2510.08201, 1812.02675, 2005.02228). **Strong fit** for the thesis if the headline phrasing leans on the survey/classification character of §4. The colon role is `statement → method qualifier`, not `name → expansion`.
3. **Topic-in-context** — `Gertsenshtein effect on the spacetime curved by background magnetic field with geometric optics` (2510.17094 by Tomomatsu et al., the closest *topical* analogue in the corpus though not a Barker paper); `Cosmology of Cubic Poincaré Gauge gravity` (2506.17017); `Graviton-photon mixing. Exact solution in a constant magnetic field` (2004.02714). **Good fit** but tends to lose the survey character — usable as a subtitle pattern, less so as a headline.
4. **Mechanism-in-medium** — `Graviton-photon X in Y` (2310.04150, 2405.11786). Standard Gertsenshtein-community register, but **does not advertise the PGT-survey contribution** of this work. Admissible only with a colon-qualifier that restores it.
5. **Tool-against-problem** — `Supercomputers against strong coupling in gravity with curvature and torsion` (2206.00658, HiGGS). Barker's flagship demonstration that a software paper can — and in his solo voice does — keep the package name out of the title and lead with the physics purpose. **The single most relevant Barker precedent for the thesis.**
6. **Software-colon** — `Name: long-form expansion` (PSALTer, 2406.09500; PSALTer v2, 2506.02111). **De-prioritised by §1.** Note that even Barker's other solo software paper (Hamilcar, 2512.25007) eschews this archetype: its title is `Fast Poisson brackets and constraint algebras in canonical gravity` — the package name appears nowhere. The corpus therefore offers two Barker precedents (HiGGS, Hamilcar) where the software is announced *without* its name in the title, against one named exception (the PSALTer family). For an MSci report, the supervisor-aligned default is the HiGGS/Hamilcar pattern.

## §5 Punctuation and capitalisation

- **Colon** — used in roughly half the Barker corpus, in two distinct roles: `name → expansion` (PSALTer, Infrared foundations I/II) and `statement → method qualifier` (R² spectrum, BHL conformal, 1812.02675, 2005.02228, Weyl-invariant EC). Only the second role is admissible for the thesis; the first is the software-colon archetype.
- **Question mark** — three corpus examples (2303.11094, 2505.23894, gr-qc/0001010). Admissible only if the body delivers an unambiguous answer; for a survey thesis whose §4 verdict is differentiated across couplings, the question form risks under-committing. Default to the declarative.
- **Period / full stop mid-title** — one example (2004.02714, `Graviton-photon mixing. Exact solution…`). Not in BHL register. Avoid.
- **Dashes, em-dashes, semicolons, exclamation marks, parentheses** — **zero** corpus examples in any of the 38 titles. Forbidden for the thesis.
- **Capitalisation** — BHL/Barker corpus is overwhelmingly **sentence case** (capital on the first word and on proper nouns only): `Every Poincaré gauge theory is conformal…`, `Spectrum of pure R² gravity…`, `Geometric multipliers and partial teleparallelism…`. Title case appears only in the older or non-Barker Gertsenshtein subcorpus (`A Simple Derivation of the Gertsenshtein Effect`, `Graviton-Photon Oscillations as a Probe of Quantum Gravity`). **Use sentence case.**
- **Accents and special characters** — `Poincaré` always carries the accent (every BHL instance: 2406.12826, 2205.13534, 2101.02645, 2006.03581, 2003.02690, 2506.17017, 1812.02675, 1804.05556). $R^2$, $H_0$, etc. appear inline in math mode. No bold, no `\Large` in the printed string.

## §6 Canonical terminology (exact spelling and hyphenation)

The following noun phrases recur and are canonical; the thesis title must match the corpus orthography exactly.

- **Poincaré gauge theory** — always accented, never abbreviated to `PGT` *in a title* (the abbreviation is fine inline in prose). Cf. 2406.12826, 2205.13534, 1812.02675, 2506.17017.
- **quadratic torsion** / **quadratic Poincaré gauge theory** — 2101.02645, gr-qc/0307063. Use `quadratic` rather than `four-parameter`, `second-order`, etc.
- **ghost and tachyon free** — no hyphen between `tachyon` and `free` in the BHL convention (1812.02675, 2005.02228). For the thesis the cleaner attributive form `ghost-free` (single hyphen, no `tachyon`) is admissible and consistent with §4 phrasing in [abstract_drafting_spec.md](abstract_drafting_spec.md); the longer `ghost and tachyon free` is the more conservative match to the corpus.
- **graviton-photon** — hyphenated when used as a compound modifier (2310.04150, 2004.02714, 2405.01407, 2405.11786, 2407.09598). Capitalised only sentence-initial.
- **Gertsenshtein effect** — `Gertsenshtein` capitalised as a proper noun; `effect` lowercase. No possessive form (no `Gertsenshtein's`). Cf. 2301.02072, 2510.17094.
- **Einstein-Cartan**, **Einstein-Maxwell**, **Euler-Heisenberg**, **metric-affine** — all hyphenated, both letters of each proper-noun pair capitalised. Cf. 2406.11956, 2405.11786, 2402.07641.
- **operator space** / **theory space** — both attested in the BHL register (`across the operator space`, `restricted to a sub-region of theory space`). Use either; do not invent `Lagrangian space`.

## §7 What the corpus does not put in titles

Forbidden patterns for the thesis title, by corpus evidence:

- **Verbs of presentation** — no `we present`, `a study of`, `an investigation into`. Zero corpus examples.
- **`novel` / `new` as adjectives** — one corpus exception (2101.02645 `new quadratic torsion theories`); otherwise absent. Use only if the noun it modifies is genuinely a new construction, not a restatement. Default: omit.
- **Software-engineering metrics** — no `23k LOC`, `1700 tests`, `parallel`, `GPU-accelerated`. Zero corpus examples.
- **Parenthetical scope-limits** — no `(linearised regime)`, `(uniform B-field)`. Zero examples.
- **Forward deferment** — no `…and beyond`, `…and prospects`, `Part I` (used by Barker only when a published Part II already exists, cf. 2101.02645 / 2003.02690; 2405.08865 uses it likewise. Do not use for the thesis.)
- **Naming what the work does not cover** — mirrors the abstract-spec exclusion (no localised-geometry caveat, no plasma caveat in the title).
- **Boccaletti or baseline-validation framing** — mirrors §4 and the abstract spec. The baseline check is not the headline.

## §8 Mood and register signatures specific to Barker

Patterns observable in the Tier 1 corpus that the thesis title is allowed to imitate:

- **The "compelling case for…" subtitle device** (2406.12826). A colon-qualifier that narrows a broad headline statement into a specific argued position.
- **Framing the paper around a position** (`Every…`, `Does…`, `Can metric-affine gravity be saved?`). The title takes a stance rather than describing a topic.
- **Colon-qualifier that narrows scope** (`Spectrum of pure R² gravity: full Hamiltonian analysis`; `Ghost and tachyon free Weyl gauge theories: a systematic approach`). Lifts the title from a topic to a delivered survey/analysis.
- **`Catalogue of…` / `Spectrum of…` framing** (2506.21662, 2507.05349, 2510.08201). When the work itself is a classification, the title says so directly. Admissible for the thesis if the headline noun phrase calls the survey what it is (a coupling-level catalogue).

## §9 Placement of TIDAL in a physics-led title

Three admissible options, in supervisor-aligned order:

1. **Absent entirely** (HiGGS precedent, 2206.00658; Hamilcar precedent, 2512.25007). The title carries the physics; TIDAL is named only at first mention in §3 and the abstract. **This is the default.**
2. **Trailing qualifier in the subtitle** (`…: a coupling-level survey with TIDAL`). Reads as a tool acknowledgement rather than a software announcement. Admissible if the colon-qualifier first does its scope-narrowing job; the package name is the last word of the title, not the first.
3. **Acronym-tag after a colon** (`…: a TIDAL survey`). Admissible only when the colon-clause also names a physics scope; the bare `…: TIDAL` is the software-colon archetype and is excluded by §1.

The two working titles in `report_plan.md` map as follows: the first ("Torsion in the Gertsenshtein effect: a symbolic-numerical survey of Poincaré gauge theory") is an Option-1 title but does not yet advertise the coupling-level finding (still a "topic-in-context" headline); the second ("TIDAL: a symbolic-numerical pipeline for the Gertsenshtein effect in Poincaré gauge theory") is software-colon and excluded.

## §10 Decision matrix for the next session

Three headline choices × structural archetypes; rows are ordered by §1 preference. The skeletons are prototypes — the brainstorm session fills in the named couplings.

| Headline | Archetype (§4) | Skeleton |
|----------|---------------|----------|
| **Coupling-level survey finding (recommended)** | Imperative-stance declarative (1) | `[Verb-stance about the operator space]: [colon-qualifier naming the ghost-free amplifying couplings]` — e.g. cousin of 2406.12826 |
| Coupling-level survey finding (recommended) | Topic-then-Hamiltonian / classification (2) | `[Survey noun phrase over the operator space]: a coupling-level classification` — cousin of 2510.08201 / 1812.02675 |
| Coupling-level survey finding (recommended) | `Catalogue of…` / `Spectrum of…` (§8) | `Catalogue of Gertsenshtein-active couplings in quadratic Poincaré gauge theory` — cousin of 2506.21662 / 2507.05349 |
| Methodological / Hamiltonian-analysis contribution | Topic-then-Hamiltonian (2) | `Gertsenshtein conversion in Poincaré gauge theory: full Hamiltonian survey` — cousin of 2510.08201; admissible if §4 framing leans methodological |
| Methodological / Hamiltonian-analysis contribution | Tool-against-problem (5) | `[Tool noun] against [physics obstacle] in Gertsenshtein conversion under torsion` — cousin of 2206.00658; admissible because the tool is *implicit*, not named |
| TIDAL software announcement (**non-preferred per §1**) | Software-colon (6) | `TIDAL: …` — included for completeness only; excluded by the binding rule unless §1 is renegotiated |

Place the named couplings ($\delta_1$, $\chi_1$, and any others load-bearing at submission) only in the colon-qualifier slot, never in the headline noun phrase; the headline must remain readable to a non-expert.

## §11 How to use this spec in the next session

Read this file end-to-end, then for each of the two working titles in [report_plan.md](report_plan.md) §"Open questions" #2 grade them against §1, §3, §4, §5, §6, §7, §9, and the matrix in §10. Reject any that fail §1. Generate 5–8 candidate titles, each compatible with at least one **physics-led** row of the matrix; constrain length to 9–14 words; constrain capitalisation to sentence case; constrain punctuation to at most one colon (or one question mark) and nothing else from §5's forbidden list. Trim to a shortlist of three; circulate to Barker. Final title locks in Phase 1 of [report_plan.md](report_plan.md).
