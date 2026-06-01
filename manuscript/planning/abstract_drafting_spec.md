# Abstract — drafting spec (TIDAL MSci thesis)

Binding stylistic and structural specification for the thesis Abstract. Distilled from a thorough read of the BHL / Barker / Gertsenshtein corpus, the existing `manuscript/style_intel/genre_conventions_theory.md` §(b) and `manuscript/style_intel/style_analysis.md` §(a), and the supervisor's direction this session (coupling-level structural framing; Boccaletti is *not* the headline; plane-wave amplification magnitudes are not the abstract's job; the abstract does not flag what is *not* in the work). This file is the drafter's reference for length, shape, register, content, and exclusions. It does not duplicate `STYLE_GUIDE.md` or `genre_conventions_theory.md` §(b) — cross-reference those for general rules. Companion to `conclusion_drafting_spec.md`.

The current `manuscript/sections/abstract.tex` scaffold is treated as outdated and is not binding. Slot allocations and the content checklist below supersede the scaffold's TODO list (in particular, the scaffold's instruction to "lead with Boccaletti validation" and to centre the closing on the §4 nulls is incorrect for the present headline).

## Context

The thesis main body now runs ~4 900 w (intro 899, theory 1 060, computational approach 242, results 1 770, discussion 1 031 once filled — plus appendices). The abstract is the most-read section of the document: examiners, supervisors, and any future arXiv reader meet the work here first. Per `manuscript/planning/report_plan.md`, the abstract is written last, after all main-body sections are locked.

The §4 headline has shifted from "mostly null with one positive structural finding" to a Bayes-factor sector survey of the quadratic Poincaré gauge theory operator space. Several Lagrangian couplings — the $\delta_1$ non-minimal EM–Riemann–Cartan curvature coupling, the $\chi_1$ derivative coupling, and others under finalisation — admit ghost-free Gertsenshtein amplification relative to the Einstein–Maxwell baseline; the parity-odd minimal sub-sector suppresses; controls separating propagating from non-propagating torsion isolate which structural ingredient drives each verdict. The plane-wave geometry used for the survey is *not* the physically realistic configuration (localised dual-Gaussian + wavepacket campaigns are deferred to follow-up work); the abstract therefore frames its findings at the *coupling level* — which Lagrangian couplings admit / forbid amplification — rather than quoting plane-wave amplification magnitudes that would over-claim physical reach.

## What is the Abstract, and what is it not?

- The abstract is **distinct from** §1 Stage 4 ("In this report"). Stage 4 announces the contribution to a reader who has already accepted the framing; the abstract delivers the contribution to a reader who has not yet decided whether to read further.
- The abstract is **not** a re-summary of §5 / §6. §5 carries implications and constructive paths; §6 closes the document. The abstract states only what *was done* and *what was found*, in present-perfect register.
- The abstract is **not** a software-paper announcement. TIDAL is named, but the headline is the coupling-level structural finding, not the package itself. PSALTer (2406.09500) is the closest software-first abstract in the corpus — it is the secondary template here, not the primary one. The primary template is the Barker–Hobson–Lasenby physics-paper abstract (2406.12826, 2101.02645, 1812.02675) where the *finding* leads and the *instrument* serves it.
- The abstract is **not** a forward pointer or deferment. No "will be explored in future work"; no companion-paper deferment; no mention of what the present work does not cover.
- The abstract is **not** a Boccaletti-validation announcement. That the pipeline reproduces the Einstein–Maxwell baseline is taken for granted by virtue of comparing to it; elevating it to the abstract would frame an expected baseline check as a thesis result. (Mirrors the conclusion-spec exclusion of Boccaletti from §6.)

## Corpus surveyed

Thirteen abstracts harvested end-to-end from the local TeX corpus and the wider Gertsenshtein literature: **2406.12826, 2303.11094, 2101.02645, 2510.17094, 2406.17853, 2205.13534, 2311.11790, 2406.09500 (PSALTer), 2512.25007 (Hamilcar), 2206.00658 (HiGGS), 2510.08201 (R² spectrum), 2301.02072 (Palessandro–Rothman), 1812.02675 (Nair–Park–Yoon).** Word counts span 50 (2311.11790 letter) to 285 (2303.11094); the modal range across the eleven non-letter abstracts is 120–200 w. Patterns reported below are synthesised across all thirteen; no single paper is the anchor. Cross-paper synthesis is also drawn from `genre_conventions_theory.md` §(b) and `style_analysis.md` §(a).

## Length

Target: **180 w**, tolerance **150–220 w**. Justification: the corpus modal band for full-paper abstracts is 150–200 w; `genre_conventions_theory.md` §(b) prescribes 150–200; the thesis hybrid character (physics survey + named software package + four content slots) pushes toward the upper end of that window. Below 150 w one slot must be cut, and the obvious candidate (Slot 4) carries the closing-line register. Above 220 w the abstract begins to compete with §1 Stage 1 for context-setting word budget. 180 w is the smallest figure that hits all four slots without crowding any one of them.

The abstract is a **single substantive prose paragraph**. No bullet lists, no enumerations, no headings, no `\paragraph*{}` micro-structure. This is unanimous across the corpus.

## Heading

`\begin{abstract} ... \end{abstract}` exactly as the revtex template provides. No `\section*{Abstract}` wrapper. The abstract environment is the only structural marker required.

## Structural shape (4 slots, one paragraph)

The corpus shows a strict four-part contract for full-paper physics abstracts. The TIDAL abstract follows it verbatim. All word targets are approximate; tolerate ±10 w per slot.

1. **Phenomenon (≈ 30 w, 1–2 sentences).** Open directly on the physics. The Gertsenshtein effect — gravitational-to-electromagnetic conversion of waves in a background magnetic field — is named and contextualised within the high-frequency gravitational-wave detection programme. Do *not* begin "We present…" or "In this paper…". Do *not* begin with historiography ("Since Gertsenshtein 1962…"). Open on the active research tension: Einstein–Maxwell is the standard baseline; whether and how torsion modifies the channel across the Poincaré gauge theory operator space is unknown. One sentence is usually enough; two if the channel and the open question need separating.

2. **Tool + capability (≈ 40 w, 2 sentences).** Announce TIDAL by name, expanded once (Tensor Integration and Derivation for Any Lagrangian). State precisely what it does: derives the linearised equations of motion directly from any quadratic gauge-gravity Lagrangian and propagates them numerically in a chosen background. Frame automation as the *physics-scale requirement* (operator space too large for manual derivation) per `genre_conventions_theory.md` §(b) and the PSALTer / HiGGS / Hamilcar precedent — *not* as a software-engineering achievement. One sentence on the tool, one on the capability that makes the survey possible.

3. **Coupling-level structural findings (≈ 70 w, 3–4 sentences).** The substantive scientific content of the abstract. Frame at the *coupling level* — i.e. by named Lagrangian coupling, not by "theory class" or "sector" or "campaign-grouping name". The campaign labels in `results.tex` (Yang–Mills–PGT base, Barker extension, Shapiro extension, $\chi$-closure, $\xi$-kinetic closure, etc.) are *mostly arbitrary groupings* used to organise the sweep and are not physics-level natural classes; the abstract must speak about the specific *couplings within those groupings* that were found to carry structure (candidates from §4: the $\delta_1$ non-minimal EM–Riemann–Cartan coupling, the $\chi_1$ derivative coupling, and any others that crystallise as load-bearing at submission). The corpus precedent for coupling-level enumeration in an abstract is the Nair–Park–Yoon classification register (1812.02675), which closes by enumerating viable cases by their operator content rather than by external naming labels. The TIDAL abstract mirrors that register: name the *couplings* that introduce Gertsenshtein amplification, name those that suppress, and assert that the amplifying configurations are ghost-free at linear order (per `genre_conventions_theory.md` §(g) Pattern 4). **Do not quote plane-wave amplification magnitudes** (e.g. "~450:1 in T1", per-sector log B values) — the plane-wave geometry is not the physically realistic showcase and the magnitudes would over-claim physical reach.

4. **Closing (≈ 40 w, 1–2 sentences).** Definitive statement of what the survey *delivers*, framed in the visionary-infrastructure or structural-finding register (`genre_conventions_theory.md` §(h) admissible registers; same as §6 conclusion). The closing sentence is forward-direction in the sense that it asserts what the work delivers and what TIDAL enables across the theory space — *not* in the sense of pointing to follow-up settings, localised-geometry extensions, or any other work-not-yet-done. The abstract does not name what is not in the work. Acceptable patterns: "TIDAL identifies [the coupling content] as the locus of Gertsenshtein amplification, and makes the coupling-level survey tractable for any quadratic gauge-gravity Lagrangian." (Visionary-infrastructure register; do not copy verbatim, paraphrase.) Or: "The amplifying coupling content is restricted to [named couplings]; these configurations are ghost-free at linear order." (Structural-finding register.) Both are corpus-attested; mixing them in a single sentence is not.

## Closing-line register — what to imitate, what to refuse

`style_analysis.md` §(a) records that the corpus's abstract closer is "a definitive statement or emphatic characterisation, not a promise." The conclusion-spec's catalogue of admissible / inadmissible registers transfers to the abstract with one adjustment (software-paper precedent enters the admissible set, because the corpus's software-paper abstracts close with capability-statements that *are* definitive).

**Admissible:**

- *Visionary-infrastructure statement* (HiGGS / PSALTer / 2406.12826 register). Asserts what the survey or tool changes about what is dismissable from the theory space. Strongest analogue for the TIDAL abstract because automation is exactly what makes the coupling-level survey possible.
- *Structural-finding statement* (1812.02675 / 2510.17094 register). Closes on the precise coupling content of an exact identity or constraint: "The amplifying coupling content is restricted to …" Admissible if the §4 headline crystallises as a clean named set.
- *Capability-with-scope-limit* (PSALTer register). "The initial release allows for parity-preserving operators …" — definitive, forward-direction, no aspiration. Admissible because it ties the tool announcement to a concrete scope statement; not admissible if it deteriorates into "future versions will …" wording, and not admissible if the scope-limit clause crosses into naming what the present work does not cover (e.g. localised geometry).

**Inadmissible (corpus-attested but not for this abstract):**

- *Boccaletti-validation closing* — would frame the baseline check as the thesis result. Excluded for the same reason §6 excludes it (supervisor direction this session; mirrors the conclusion-spec exclusion).
- *Plane-wave amplification magnitude closing* — over-claims physical reach. The geometry is acknowledged not to be physically realistic; the localised-geometry magnitudes are deferred to follow-up work. Quoting "~450:1 in T1" or equivalent Bayes-factor magnitudes in the abstract would invite a referee to ask for the physical conversion probability under realistic background fields, which is a question this thesis explicitly does not answer.
- *Companion-paper deferment* (2205.13534 register: "consider their tree-level phenomena", "to be analysed in subsequent work"). There is no companion paper.
- *Aspirational closing-line wording* ("we hope that …", "it would be interesting to …", "future work should explore …", "this will be explored in subsequent work"). Forbidden as a *register*, not as content.
- *Naming what the present work does not cover.* The localised-geometry next-step content belongs in §5 / §6, not here.

Draft the closing sentence first, in isolation, before drafting Slots 1–3 — and revise it until it carries one of the three admissible registers without hope/expectation verbs and without quoting plane-wave magnitudes.

## Voice, citations, equations, cross-references

- **First-person plural** throughout ("we present", "we identify", "we constrain"), per `genre_conventions_theory.md` §(j) and the BHL solo-paper precedent (2205.13534, 2510.08201) where solo authors still write "we".
- **Zero equations.** The full corpus places zero equations in the abstract environment (the closest precedent is the Boccaletti-formula display in gr-qc/0307063, a single-calculation paper unlike TIDAL). Display formulae belong in §2 or §4.
- **Zero figures and zero `\cref{}`.** Abstracts in this corpus never cross-reference the body. A reader who reads only the abstract should grasp the contribution without flipping pages.
- **Zero in-text citations.** Corpus-attested across all thirteen abstracts surveyed; none carry a `\cite{}` in the abstract environment. The bibliographic anchoring belongs in §1.
- **TIDAL is named and acronym-expanded once.** "TIDAL (Tensor Integration and Derivation for Any Lagrangian)" — PSALTer / Hamilcar / HiGGS precedent. Subsequent mentions in the abstract use "TIDAL" only.
- **"PGT" is not used unexpanded.** Corpus check across 2406.12826, 2205.13534, 1812.02675, 2206.00658 confirms two admissible patterns: either expand on first use ("Poincaré gauge theory (PGT)…" — 2406.12826 register; subsequent mentions in the abstract then use "PGT") or use the full form "Poincaré gauge theory" throughout with no acronym at all (2205.13534, 1812.02675, 2206.00658 register). The bare unexpanded "PGT" is corpus-absent in abstracts. Pick one of the two admissible patterns at draft time; default to expand-on-first-use because the abstract refers to the theory multiple times.
- **Forbidden words** from `STYLE_GUIDE.md` §g apply with extra force here: no "it can be seen that", "it should be noted that", "the results indicate", "novel", "state-of-the-art", "obviously", "in this paper" (use "we" + verb instead).

## Mandatory content checklist for the drafter

The drafter must, at writing time, ensure the abstract contains each of the following — none more, none fewer:

- One-sentence statement of the Gertsenshtein phenomenon in Poincaré-gauge-theory context, with the open question (the coupling-level structure of the channel across the quadratic gauge-gravity action) named explicitly.
- One- or two-sentence announcement of TIDAL by name and acronym expansion, with the symbolic-to-numerical capability stated as a physics-scale requirement (operator space too large for manual derivation), not a software achievement.
- The structural / coupling-level findings of the §4 survey, framed by named Lagrangian coupling (not by campaign-grouping name): which specific couplings introduce ghost-free Gertsenshtein amplification, which suppress, and the ghost-free status of the amplifying configurations stated explicitly. **No** plane-wave amplification magnitudes; **no** Boccaletti agreement percentage; **no** campaign-grouping labels masquerading as physics-level classes.
- A closing sentence in the visionary-infrastructure, structural-finding, or capability-with-scope-limit register, with an active subject and a concrete verb. No aspirational locutions. No forward-pointer deferment. No mention of next settings the present work does not cover.

## Forbidden in the abstract

- Any equation, figure, or `\cref{}` to a §2/§4/appendix object.
- Any in-text `\cite{}`. The full corpus has zero in-abstract citations across the thirteen papers surveyed.
- Any framing that elevates the Boccaletti validation to a headline finding, or that announces the validation precision as a thesis result.
- Any quoted plane-wave amplification magnitude (per-sector log B values, "~450:1 in T1", Bayes-factor numerical anchors). The plane-wave geometry is not the abstract's setting for physical-magnitude claims.
- Any companion-paper / sequel-paper / "in future work we will" framing.
- Passive-aspirational locutions: "we hope", "future work should explore", "it would be interesting to", "this will be explored in subsequent work".
- Any explicit naming of what the present work does not cover (the localised-geometry deferral, follow-up campaigns, the EFT-perspective interpretation, etc.). That content belongs in §5 / §6.
- Bullet lists, numbered enumerations, sub-paragraph headings, or `\paragraph*{}` markers.
- "TIDAL is a Python package that …" software-paper framing as the *headline*. TIDAL is the instrument; the headline is the coupling-level finding.
- Bare unexpanded "PGT". Either expand on first use ("Poincaré gauge theory (PGT)") or use the full form throughout.
- More than 220 w by `texcount` on the abstract block alone.

## Open decisions for the drafter

- **The exact set of named couplings** to enumerate in Slot 3 depends on the final state of §4 at submission. The right granularity is the *coupling*, not the campaign-grouping name. Candidates emerging from `results.tex` and the §4 discussion as load-bearing structural findings: the $\delta_1$ non-minimal EM–Riemann–Cartan coupling, the $\chi_1$ derivative coupling, and whichever other Lagrangian couplings crystallise as the actual locus of the amplification verdict once the discussion section is filled. Draft Slot 3 with placeholder `[COUPLINGS: amplifying = …; suppressing = …]` and resolve the names from §4 / §5 at submission. Two to three named couplings is the maximum the slot can carry without crowding; pick the most informative subset, and refuse the temptation to list the campaign-grouping labels (Yang–Mills–PGT extensions, $\chi$-closure, $\xi$-kinetic closure, etc.) which are organisational artefacts of the sweep rather than physics-level findings.
- **Whether to include a structural side-note about torsion decoupling at minimal coupling.** §4 supplies one candidate compact structural observation in this register — that minimal-coupling configurations decouple torsion from the Gertsenshtein channel, isolating the amplification effect to the non-minimal couplings. This is the *only* structural identity worth considering for the abstract; the hx↔ax / kinetic-matrix coincidence and similar in-text observations are not abstract material. The decoupling side-note is optional: include it only if it reads as a tight one-clause structural framing that explains *why* the amplifying couplings sit where they do, and only if doing so does not push Slot 3 over 80 w.
- **Closing-sentence register**: visionary-infrastructure vs structural-finding vs capability-with-scope-limit. Pick at draft time. If the §4 headline crystallises as a clean coupling-content identity (a small named set of couplings carrying the amplification verdict), the structural-finding register is preferred. Otherwise the visionary-infrastructure register, asserting that TIDAL makes the coupling-level survey tractable across any quadratic gauge-gravity Lagrangian, is the safest default.
- **Naming any next setting (localised geometry, etc.) is not an open decision — it is excluded.** The abstract does not flag what is *not* in the present work. The PSALTer / Hamilcar "initial release allows for …" pattern is corpus-attested but inadmissible here: it crosses into "what we did not do" framing, which the supervisor has ruled out for this abstract. The BHL physics-letter convention (2406.12826, 2311.11790, 1812.02675 — none of them name a next setting in the abstract) is the binding precedent. The closing sentence states what was done in present-perfect register, full stop; the "natural next setting" content belongs in §5 / §6, not here.

## Quality bar before the abstract is considered done

- The abstract is 150–220 w by `texcount` on the abstract block, with 180 w as the working target.
- Zero in-text citations, zero `\cref{}`, zero equations, zero figures, zero `\paragraph*{}` markers, zero bullet lists.
- TIDAL is named once with acronym expansion. "Poincaré gauge theory" is either expanded on first use as "Poincaré gauge theory (PGT)" or used in full throughout; never appears as bare unexpanded "PGT".
- The Boccaletti validation does not appear as a headline result and ideally does not appear at all (admissible as a single sub-clause inside Slot 2 only if the drafter judges the validation provenance load-bearing; default is to leave it out).
- No quoted plane-wave amplification magnitude.
- The closing sentence ends on an active forward verb in one of the three admissible registers; it does not contain "hope", "expect", "future work", "companion paper", "this will be explored", "in conclusion", "to summarise".
- The closing sentence does not name any setting the present work does not cover.
- A reader who reads the abstract cold (no prior pages) can state in one sentence: which Lagrangian couplings of the quadratic gauge-gravity action TIDAL identifies as admitting Gertsenshtein amplification, and which suppress.
- Self-review against `STYLE_GUIDE.md` §k checklist passes the items relevant to a single-paragraph abstract.
- `latexmk -pdf main.tex` still compiles cleanly; main-body `texcount` remains ≤ 5 000 w (the abstract is excluded from the main-body cap per `report_plan.md`, but verify the build is unaffected).

## Verification

- Open `manuscript/main.pdf` after compile and confirm the abstract reads cleanly as a single paragraph with no orphan widow on the next page.
- Read the abstract aloud once: a closing sentence with the wrong register is much easier to hear than to see.
- Diff the abstract content against the §1 Stage 4 hook in `introduction.tex` and the §4 headline table in `results.tex`; the named couplings and amplifying/suppressing verdicts must match. If they do not, the §4 table is canonical — adjust the abstract, not the table.
