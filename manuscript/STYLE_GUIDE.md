# MSci Report Style Guide

Actionable rules for drafting, derived from analysis of the Barker / BHL corpus.
Check each rule before sharing any draft section.

---

## (a) Section Opening Rule

Every section and every subsection must open with a complete declarative sentence that states the section's purpose or finding. Do not begin with a definition, an equation, or a citation. Do not begin with "In this section we will…" — state the claim directly. Named paragraphs (`\paragraph*{Topic}`) within sections may begin with a topic-stating sentence or an assumption statement.

---

## (b) Equation Rule

Every displayed equation must be preceded or followed by at least one complete prose sentence that interprets it (not just names it). The sentence *before* an equation states what it means; the sentence *after* draws a consequence. Never drop an equation into the text without surrounding prose. Never refer to equations as "the above" or "the following"; always use `\cref{label}`.

---

## (c) Label Rule

Every equation, section, subsection, figure, and appendix that will be cited — even once — must have a label. Labels must be semantic, not serial.

**Use bare CamelCase descriptive names**: `\label{QuadraticAction}`, `\label{TorsionTransform}`, `\label{GertsenshteinFormula}`, `\label{BoccalettiValidation}`. **Do not** use type prefixes like `eq:`, `sec:`, `fig:`, `app:`. The bare-CamelCase convention is the BHL corpus standard (verified across 2406.12826 and 2101.02645 — every label is bare). `cleveref` auto-detects the reference type from the LaTeX context, so prefixes are unnecessary; the bare names also read as a flat conceptual index of the paper rather than a typed namespace.

**Choose names that cannot clash across types**. If a section and an equation share a topic, give them distinct names (e.g. `\label{Linearisation}` for the subsection and `\label{KineticMatrix}` for the equation it produces — never two labels with the same bare name).

**Do not label `\paragraph*{...}` headers** unless they are explicit cross-reference targets. Named paragraphs are organisational devices, not citeable objects in the BHL convention.

**Cross-references to labels in `docs/tex/`**: those files use their own conventions (some keep the `sec:` prefix). Cross-references from manuscript scaffolding to docs labels keep whatever prefix the docs use — for example, `\cref{sec:kappa-conventions}` resolves to a label in `docs/tex/gertsenshtein_formula.tex` that uses the prefix.

---

## (d) Voice Rule

Use first-person plural throughout, including in solo-authored work: "we show", "we find", "we confirm". Use "we" even for results that are entirely your own computation. Never use "it can be seen that", "it should be noted that", or "the results indicate". State findings actively: "this confirms X" not "X appears to be confirmed".

Hedge only when genuinely uncertain. When hedging, state the finding first, then the caveat: "vector torsion appears viable — pending a full nonlinear analysis." Do not front-load with caveats.

---

## (e) Footnote Rule

Footnotes are for substantive asides that experts need but which interrupt the main argument. Use them for: technical caveats that a specialist reader may want but a general reader can skip; ironic or classical cultural observations that would be a tonal mismatch in the body; precise distinctions that would derail the paragraph but are genuinely important.

Do not use footnotes for: pure citations (embed those in the body), definitions of standard quantities, or content that belongs in the main text but is hard to place. A footnote that exceeds four lines should be moved to a parenthetical sentence in the body, or cut.

---

## (f) Citation Density Rule

- **Introductions**: Dense citation clusters are appropriate when surveying prior work. Cite five to fifteen sources in a single sentence for a broad claim about the field. After establishing the prior-art inventory, reduce citations to one to three per claim.
- **Body sections**: Cite only the foundational or directly relevant work for each claim. One to three citations per claim is the norm. Do not cite support for your own derivations.
- **Equations**: Citations belong in the prose surrounding an equation, not inline with the label.
- **Classical canon**: Always cite the original paper for foundational results (Einstein, Maxwell, Proca, Klein–Gordon), not a textbook.
- **Software dependencies**: Cite each dependency with its own reference (xAct, xTensor, Mathematica, SUNDIALS, numpy). "Standard tools" without references is not acceptable.

---

## (g) Forbidden Words List

Remove these from every draft:

| Forbidden | Replace with |
|-----------|-------------|
| "it can be seen that" | state the result directly |
| "it should be noted that" | state the observation directly |
| "the results indicate" | "we find" / "this shows" |
| "in order to" | "to" |
| "due to the fact that" | "because" |
| "at this point in time" | "now" / "currently" |
| "a large number of" | "many" |
| "it is clear that" | (delete the phrase; just state the claim) |
| "obviously" | (delete; if it is obvious, state it plainly) |
| "interesting" (as standalone adjective) | describe *why* it is interesting |
| "novel" (in abstract or intro) | describe specifically what is new |
| "state of the art" | cite the specific frontier work |

---

## (h) Figures and Caption Rule

Every figure must be cited in the main text before it appears, using `\cref{fig:label}`. Captions must be self-contained: a reader who sees only the figure and its caption must be able to understand what is being shown and what the key result is. Captions should: (1) say what the figure shows, (2) state the parameter regime or calculation conditions, (3) identify the key result visible in the figure. Captions do not end with "see text for details".

Label axes with physical quantities and units. Do not use variable names alone (not "$\xi$" but "coupling constant $\xi$"). If colour coding is used, include a legend or explain colour assignments in the caption.

---

## (i) Per-Section Templates

Use these when drafting each section. Each entry names the closest-match template paper and the specific structural moves to imitate.

---

### §1 Introduction

**Template papers:** 2206.00658 §I (HiGGS — software-physics hybrid intro); 2510.08201 §I (R² — concise, paradox-driven); 2406.12826 §I (conformal PGT — letter format, dense).

**Structural moves:**
1. Open with a living research tension, not a historical overview. The tension is: GW detectors are probing new physics, but the gauge-gravity sector they probe is theoretically vast and poorly constrained.
2. State the deep problem in one paragraph: the Gertsenshtein channel structure in PGT is unknown, and computing it manually is infeasible for more than one model at a time.
3. Survey prior art: Gertsenshtein/Boccaletti (mechanism), BHL (PGT landscape), PSALTer/Hamilcar (spectrum tools), TIDAL (this work). Dense citation cluster here (8–15 refs in 2–3 sentences is normal in BHL intros). Identify the gap: no systematic channel decomposition for PGT + EM has been done.
4. Hook: "We present TIDAL, a symbolic–numerical pipeline that automates this calculation for any quadratic PGT+EM Lagrangian. We apply it to..." One paragraph.
5. Close with one sentence per remaining section (the structural map).

**Sentence-level pattern for opening:** Begin with the phenomenon, not the history. Compare: ✗ "The Gertsenshtein effect was discovered by..." vs ✓ "When gravitational and electromagnetic waves propagate through a background magnetic field, they mix — a phenomenon with direct consequences for high-frequency GW detection."

---

### §2 Theory

**Template papers:** 2406.12826 §II (PGT framework); 2101.02645 §II (quadratic action); 2510.17094 §II–III (linearisation around background); 2510.08201 §II (action + ADM).

**Structural moves for §2.1 (PGT framework):**
1. One sentence recalling the Riemannian baseline (to fix notation, not to teach GR).
2. Two sentences on the PGT generalisation: connection promoted to independent field; field strengths are torsion and curvature.
3. Paired equations: curvature tensor, then torsion tensor, with identical syntactic structure.
4. General quadratic action with coupling constants `\Bet{1}`, `\Bet{2}`, `\Bet{3}` (torsion²), etc. State: this is the most general action; special cases correspond to known theories.
5. Irreducible decomposition in a named paragraph or subsection after the action.

**Structural moves for §2.2 (linearisation):**
1. State the background: flat metric $\eta_{\mu\nu}$, uniform static $B_0$ in the $z$-direction.
2. Write the perturbation ansatz: $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, $A_\mu = A_\mu^{(0)} + \delta A_\mu$. State the order of truncation.
3. State the order-counting rule immediately after the ansatz: "We work to linear order in $h_{\mu\nu}$ and $\delta A_\mu$, and to first order in $B_0$." (Pattern from 2510.17094 §II.)
4. Derive kinetic matrix — credit the algebra with "after expanding the Lagrangian to quadratic order in the perturbations, one finds..." without displaying intermediate steps.

**Do not** open §2.1 with a definition of a Lie group or a review of Yang-Mills. The BHL template assumes the reader knows gauge theory; §2.1 fixes the PGT-specific notation.

---

### §3 Computational Approach

**Template papers (revised 2026-05-30 after wider corpus audit).** The genuine comparators are physics-first papers that USE software, not standalone software papers.  Verified examples with direct quotations:

- **2101.02645** (Barker+Lasenby+Hobson+Handley, uses HiGGS) — HiGGS is mentioned briefly in §I; implementation lives in §II–§III (~3000 w inline because that paper has no comprehensive appendix).
- **2206.00658:168** (HiGGS software paper, introductory-sentence pattern) — *"The core of the implementation is (version 1.0.0 of~\cite{supp2}) the \emph{Hamiltonian Gauge Gravity Surveyor} (HiGGS), a Mathematica package grounded in the popular open-source xAct tensor manipulation suite~\cite{Martin-Garcia:2007bqa}."*
- **2406.12826:331** (BHL letter using PSALTer) — *"The particle spectra are obtained using the Particle Spectrum for Any Tensor Lagrangian (PSALTer) software~\cite{…}."*  Single sentence in methodology.
- **2510.08201:466** (R² gravity using Cadabra) — single-sentence footnote: *"…is greatly facilitated by the use of Cadabra~\cite{…}."*

**Single section, no subsections.** Target: ~80 words, 3 named paragraphs of 1–2 sentences each.

**Why so brief (corpus + word-budget argument).** Software-first paper conventions (HiGGS abstract, PSALTer §I "Previous advances", PSALTer v2 "computational barriers") prescribe heavy motivation, but those are *standalone software papers* where the tool IS the contribution.  In physics-first papers that USE software (the comparators above), tool motivation is essentially absent in the main body — and our ~8000 w App A–J carries the implementation in full, so §3 does not need to.  The MSci 5000 w cap also forces brevity: §1 (800) + §2 (1200) + §4 (1700) + §5 (600) already consume ~4300 w; §3 at ~80 w fits the remaining envelope.  The 2101.02645 inline-methodology pattern (~3000 w in main body) is structurally unavailable to us.

**Voice (verified empirically across the corpus).** Passive structural is the dominant voice in tool-mention paragraphs in physics-first papers (2406.12826:331, 2510.08201:466 are both passive; 2206.00658:168 is declarative-structural).  Active first-person plural ("we [verb]") is acceptable when claiming methodological agency (2206.00658:150, 2510.08201:214).  Mixed voice within a section is corpus-attested (2206.00658 itself mixes within a paragraph).  For §3: ¶1 may open active ("We obtain…") or passive ("The equations are obtained…"); ¶2 and ¶3 typically open passive.  Both are valid; do not over-prescribe.

**Acronym introduction order (verified BHL convention).** Use **Full Name (Acronym)** on first mention, NOT Acronym (Full Name).  Per 2206.00658:168: *"the \emph{Hamiltonian Gauge Gravity Surveyor} (HiGGS), a Mathematica package…"*  Italicise the full name with `\emph{}`; the acronym is plain.  If a parenthetical needs to carry both the acronym and a `\cref{}` cluster, separate them with a **comma**, not a semicolon (no surveyed BHL parenthetical uses an internal semicolon for this).

**Structural moves (3 paragraphs):**

1. **¶1 — TIDAL (~32 w, 1 sentence).**  Single sentence in the HiGGS line-168 pattern: name the tool with **Full Name in `\emph{}` followed by acronym in parens**, foundation cite for `\xAct{}~\cite{martingarcia2007xact}` and `\Mathematica{}~\cite{mathematica}`.  Optionally include the App A/B/C cref cluster in the same parenthetical, separated from the acronym by a comma.
2. **¶2 — Numerical evolution (~28 w, 2 sentences).**  Passive structural opening: *"The linearised equations are evolved with \TIDAL{}'s Fourier modal solver (\cref{ModalSolver})."*  Second sentence on the technical method: per-Fourier-mode matrix exponentials computed by **Padé scaling-and-squaring** via `scipy.linalg.expm`, cite `higham2009scaling`.
3. **¶3 — Parameter inference (~18 w, 1 sentence with semicolon).**  Passive: *"Parameter exploration is by nested sampling with PolyChord~\cite{Handley:2015fda,Handley:2015vkr}; the inference architecture is documented in \cref{InferenceArchitecture}."*  Both PolyChord papers cited per the PolyChordLite README requirement.

**Do NOT include in §3** (each item explicitly vetoed during the four-round revision investigation; preserving these vetoes is the purpose of this section):

- **Architectural breakdown of the symbolic→numerical pipeline** — App A's job.  §3 just names TIDAL and cross-refs.
- **"Why automation is necessary" sentence** — corpus shows tier-2 physics papers omit it.  Methodological-scope motivation belongs in §1 Stage 2 (template: 2101.02645's *"singular kinetic Hessian incentivises the Hamiltonian approach"*), not §3.
- **Fallback solver enumeration (IDA, CVODE, leapfrog, scipy)** — all §4 results used the modal solver; defer to App C.
- **"Auto-selected on flat-metric, periodic, time-independent systems" eligibility qualifier** — App C eligibility conditions, not main-paper material.  Say the modal solver is *used*, not that the dispatcher *chose* it.
- **"Machine-precision evolution over the long integration windows"** — TECHNICALLY MISLEADING.  The modal solver is O(1) in `t_end` (single matrix-exponential application per mode), not stepped.  Use "single application of its matrix exponential" instead.  App C `numerical.tex:275-277` is explicit.
- **"Per-mode eigendecomposition" claim** — TECHNICALLY INCORRECT.  TIDAL explicitly abandons eigendecomposition because `cond(V) > 1e13` in the multi-field PGT regime; uses Padé scaling-and-squaring instead.  See App C `numerical.tex:315-344`.
- **Validation paragraph** — belongs in §4.1 (Boccaletti calibration anchors §4.1 directly).
- **"tachyon-permissive" inference label, or enumeration of compactified priors / soft-floor likelihood / `run_status` tags** — defer to App J (the inference-architecture appendix) and §4 (where the inference outputs need the language).
- **Closing "Detail deferred" paragraph** — none of the BHL comparators contain a standalone "deferred to appendix X" closing.  Spending words to announce non-content is the wrong move.
- **"see X" connector before crefs** — BHL writes bare inline anchors (verified across 2406.12826, 2101.02645).
- **PSALTer / Hamilcar / HiGGS peer citations** — BHL convention cites the foundation (xAct) as the only software dependency, never sibling tools as ecosystem positioning.

**Note on §(i) Appendix Rule cref-anchor requirement.** The codified §(i) Appendix Rule below says *"every main-text cross-reference to an appendix includes a one-clause explanation of what is deferred"*.  The empirical corpus does NOT support this rule for brief tool-mention contexts — 2206.00658:168 has no anchor on `~\cite{Martin-Garcia:2007bqa}`; 2406.12826:331 cites PSALTer without anchor.  The rule should be read as **"optional anchor for brief tool-mention crefs; required anchor for substantive main-text deferrals to results subsections or appendix derivations"**.  Bare `\cref{}` in §3 parentheticals is corpus-attested and acceptable.

---

### §4 Results

**Template papers:** 2406.12826 §V (conformal PGT results); 2303.11094 §III–IV (null + constraint framing); 2510.17094 §IV (competing-effects result).

**Opening:** One sentence mapping the four subsections. "§4.1 validates TIDAL against the Boccaletti formula; §4.2–4.3 present the null results for the T1 and T2 sectors; §4.4–4.5 present the propagating-torsion findings."

**§4.1 Boccaletti validation:** Lead with the formula. Show fig. 1 (simulation vs. analytic). State: "TIDAL reproduces the Boccaletti conversion probability to 0.04%." This is the positive result that validates everything downstream.

**§4.2–4.3 Null results:** Use the Dandoy framing (2406.17853): "we constrain the [T1/T2] sector conversion amplification to $A < X$ across the full parameter plane (Fig. 2)." Never write "we did not find amplification". Explain the physical mechanism in one sentence: "the trace coupling vanishes at linear order because..." Each subsection opens directly on the result — do not repeat the §4 roadmap sentence inside each subsection.

**§4.4 Propagating torsion:** Use the 2510.17094 plane-wave-cancellation framing: state the two sectors (hx and ax), show they enter the kinetic matrix identically, state the consequence (torsion-independence), cite App F for the derivation. One paragraph, one figure.

**§4.5 Ghost diagnosis:** Use the 2510.08201 ghost-reporting pattern: state the no-ghost condition first, show it is violated in the trace sector, state the consequence. Do not call this a "failure" — call it a "diagnosis" that "identifies the trace channel as ghost-unstable at the parameter point..."

---

### §5 Discussion

**Template papers:** 2406.12826 conclusions; 2303.11094 §V; 2510.08201 §V.

**Three-move structure (see genre_conventions_theory.md §h):**
1. Findings restatement: plain language, no equations, ~150 words.
2. Implications: what the null + torsion-independence mean for PGT as a physical model; the "three constructive paths" enumerated here (not in §4).
3. Specific next steps: name one concrete theory (the non-minimal PGT with torsion–curvature cross-terms) and one concrete observable (non-zero channel amplification if the cross-term coefficient enters the kinetic matrix). This last paragraph should be completeable in one sentence per path.

**Do not** put the three constructive paths in §4. They are speculative implications, not results. §4 reports what was found; §5 interprets it.

---

---

## (i) Appendix Structure Rule

Every appendix carries a semantic title that states its content (not "Appendix A: Technical details" but "Appendix A: TIDAL pipeline architecture"). Each appendix opens with a single declarative sentence stating its purpose and what it contains — this is the appendix's equivalent of the section opening rule (a).

Appendix subsections are numbered A.1, A.2, etc. Each subsection opens with a declarative sentence. Appendix figures are numbered in sequence with main-body figures (not restarted at A.1).

**Cross-referencing convention in the main text**: at every point where a main-body sentence would otherwise spend words on implementation detail, replace the detail with a one-clause forward reference: "We solve the linearised system using the Fourier modal backend (\cref{app:modal}), which auto-selects on flat metrics with periodic boundary conditions." The clause after the `\cref{}` is mandatory — it tells the reader why they might want to follow the cross-reference.

**Per-appendix structural guidance for TIDAL**:
- **App A (architecture)**: begin with a 2–3 sentence overview of the pipeline stages (symbolic → JSON → numerical), then a schematic diagram (TikZ), then prose expanding each stage. ~1800 words + 1 diagram.
- **App B (symbolic stage)**: begin with the xAct/xPert environment and the Euler–Lagrange derivation procedure. Show representative code only if it serves as user documentation; avoid listing internal routines. ~1600 words.
- **App C (numerical stage)**: begin with the Fourier modal solver (the primary backend) — describe the eigendecomposition algorithm and auto-selection logic. The other backends (IDA, CVODE, leapfrog, scipy) each get one named paragraph. ~2000 words.
- **App D (validation suite)**: four named paragraphs following the Barker validation hierarchy (analytic limits → reference models → pathological cases → scale-up). Each named paragraph has a figure; caption supplies the quantitative reading; prose gives the interpretation. ~1200 words.
- **App E (HPC infrastructure)**: one paragraph on the CSD3 campaign setup, one on reproducibility (commit hash, TIDAL version, campaign IDs keyed to figures). A table of campaign parameters is appropriate. ~1000 words.

---

## (k) Self-Review Checklist Before Sharing a Draft

Before circulating any draft section, confirm all of the following:

- [ ] Introduction has four identifiable stages: Frontier, Deep Problem, Prior Art, Hook
- [ ] §2 has a "Conventions" named paragraph declaring signature, index conventions, and $\kappa^2 = 16\pi G$
- [ ] §2 states explicitly whether working in Palatini or metric formulation
- [ ] §2.2 linearisation has: background declared, perturbation ansatz, order-counting rule immediately after ansatz
- [ ] Every section and named paragraph opens with a declarative sentence (no equation-first openings)
- [ ] Every displayed equation has surrounding prose (one sentence before or after minimum)
- [ ] All cited equations use `\cref{}` (no bare equation numbers, no "Eq. (N)")
- [ ] All equation labels are semantic, not serial
- [ ] Coupling constants are introduced in the action equation first, then glossed in prose
- [ ] No forbidden words from the list above appear in the draft
- [ ] Every figure is cited in text before it appears; caption is self-contained; caption names the parameter regime
- [ ] All hedges are in the form "finding first, caveat second"
- [ ] Every footnote passes the "could this go in the body?" test
- [ ] Citation density is high in the prior-art inventory and sparse in the derivation sections
- [ ] Null results are framed as precision findings using "we constrain", "we exclude", "we confirm X to within Y", not "we did not find"
- [ ] Abstract ends with a definitive claim, scope statement, or emphatic result — not a promise or a limitation
- [ ] Each appendix opens with a declarative sentence; subsections are numbered A.1, A.2 etc.
- [ ] Every main-text cross-reference to an appendix includes a one-clause explanation of what is deferred
