# Genre Conventions: The BHL Theoretical Physics Paper

This document synthesises the conventions of the theoretical physics paper genre as practised by the Barker–Hobson–Lasenby (BHL) group, drawing on comprehensive analysis of eighteen papers across three thematic clusters:

**Cluster I — BHL joint papers on PGT/torsion**: 2406.12826, 2101.02645, 2309.14783, 2303.11094, 2008.09053, 2005.02228

**Cluster II — Barker solo / Barker-led recent**: 2205.13534, 2311.11790, 2510.08201, 2407.09598, 2006.03581, 2003.02690

**Cluster III — Gertsenshtein-effect and torsion-EM papers**: 2301.02072, 2406.17853, gr-qc/0307063, hep-th/0103093, gr-qc/0001010, 1812.02675

The companion file `genre_conventions_software.md` covers the software-paper genre (PSALTer, Hamilcar). The present file covers theoretical physics content — the conventions relevant to §1, §2, §4, and §5 of the TIDAL report.

---

## Paper Sub-Types in the BHL Corpus

Not all theoretical papers share the same structure. Six recognisable sub-types appear across the corpus, each with distinct genre expectations:

| Sub-type | Archetype | Key signature |
|----------|-----------|---------------|
| **BHL PGT letter** | 2406.12826, 2311.11790 | ≤10 pp, one theorem, no standalone methods section |
| **BHL long PGT paper** | 2101.02645, 2309.14783 | 20–35 pp, subsection-heavy theory, extended appendices |
| **BHL cosmology paper** | 2006.03581, 2003.02690 | Observational context up front, cosmological evolution equations |
| **Barker solo paper** | 2205.13534, 2510.08201 | Single author, direct voice, more explicit roadmap in intro |
| **Gertsenshtein derivation paper** | gr-qc/0307063, hep-th/0103093 | Calculation-centric, sparse prose, derivation is the result |
| **Observational constraint paper** | 2406.17853, 2303.11094, 1812.02675 | Observational/phenomenological framing, quantitative bounds |

The TIDAL report most closely resembles **a hybrid of the BHL PGT letter and the Barker solo paper**: a single author, a single central calculation result, a compact main body, and software contribution pushed into appendices.

---

## (a) Title Conventions

BHL theoretical papers use nominative-technical titles that encode the central claim rather than announcing a method. Three patterns dominate:

1. **Claim-as-title** (the most common): the title states the finding or identification. *"Every Poincaré gauge theory is conformal"* (2406.12826); *"Does gravitational confinement sustain flat galactic rotation curves?"* (2303.11094). The claim structure implies a proof or derivation follows; the reader knows exactly what the paper will deliver before opening it.

2. **Capability-as-title**: the title announces what the paper enables or establishes. *"Manifestly covariant variational principle for gauge theories of gravity"* (2309.14783); *"Ghost and tachyon free Weyl gauge theories"* (2005.02228). These titles are used when the contribution is a new theoretical framework rather than a specific physical result.

3. **Question-as-title**: a genuine question is posed and answered. *"Does gravitational confinement sustain flat galactic rotation curves?"* (2303.11094). Questions are used only when the answer is non-obvious enough to justify the uncertainty framing — BHL does not use question titles for results that are clearly affirmative or clearly refuted.

**For TIDAL**: a claim-as-title or capability-as-title is appropriate. The current working title (*"Torsion in the Gertsenshtein effect: a symbolic–numerical survey of Poincaré gauge theory"*) follows the capability convention with a colon-subtitle structure. Avoid question titles for null results — "Does torsion amplify..." would frame the main finding as a failure.

---

## (b) Abstract Conventions

BHL abstracts follow a strict four-part structure in 150–200 words:

1. **Phenomenon** (1–2 sentences): state the physical phenomenon or theoretical context without preamble. No "we study" — open directly on the physics. "When gravitational and electromagnetic waves propagate through a background magnetic field..." BHL abstracts never begin with the words "We" or "In this paper."

2. **Contribution** (2–3 sentences): state precisely what the paper computes or establishes. Use "we show", "we find", "we prove", "we identify". State the main result in quantitative or structural terms if possible — not "we investigate X" but "we find that X vanishes at linear order."

3. **Scope statement or secondary finding** (1–2 sentences): delimit what is included and what is not. BHL uses this slot to prevent over-reading: "the analysis is restricted to the linearised regime"; "we consider only the parity-preserving sector." Alternatively, place a secondary finding here if there is one.

4. **Implications or forward pointer** (1 sentence): close with the significance or an explicit pointer to future work. This sentence is forward-looking and positive. Do not end the abstract with a caveat or a limitation. The PSALTer abstract closes: "The initial release allows for parity-preserving operators..." — a scope statement that reads positively.

**Barker solo abstracts** (2205.13534, 2510.08201) tend to be denser in the contribution slot — more results per sentence, tighter scoping, less hedging. Joint abstracts (2101.02645, 2006.03581) expand the context slightly to ensure all co-authors' perspectives are represented.

**Gertsenshtein-paper abstracts** (gr-qc/0307063, 1812.02675) place the physical formula (the conversion probability) in the abstract itself — this is unusual and specific to calculation papers where the result *is* a formula. If the TIDAL report's main numerical result is a confirmed null (A = 1.0 to quantified precision), state that number in the abstract.

---

## (c) Introduction Architecture: The BHL Four-Stage Scaffold

Every BHL theoretical paper introduction follows a four-stage scaffold. This is not accidental — it is the most consistent cross-paper structural pattern in the corpus.

### Stage 1 — The Frontier Tension (~1 paragraph)

Open with the phenomenon that motivates the paper, stated as an active research tension. The dominant BHL pattern avoids warm-up historiography ("since Einstein introduced...") and opens directly on the current state of the field. However, note that the opening style varies across the corpus: some papers (2406.12826) begin with a mini-historical framing ("General relativity offers a remarkably successful description...") before immediately pivoting to the tension, while others (2101.02645, 2303.11094) open on the tension directly. The key criterion is that the *tension* is established within the first two sentences — not that historical background is entirely absent.

Examples of the tension statement (regardless of how the paragraph opens):
- 2406.12826: the conformal structure of PGT has been disputed for decades — but the resolution requires examining torsion's role under dilations.
- 2101.02645: the Hamiltonian analysis of Poincaré gauge theories is central to ghost and tachyon elimination — but the full nonlinear analysis has not been completed.
- 2303.11094: flat galactic rotation curves are conventionally attributed to dark matter — but modified gravity is a live alternative.

The tension is always stated in one rhetorical move: X (the phenomenon) is important because Y (the stakes), but Z (the gap) prevents us from knowing W (the specific unknown).

### Stage 2 — The Deep Problem (~1–2 paragraphs)

After the frontier tension, narrow to the specific problem the paper addresses. The deep problem is the intersection of the large tension (Stage 1) and the specific theoretical challenge that makes the problem hard. This is where technical vocabulary is first introduced — but only the vocabulary needed to state the problem.

BHL uses the second paragraph for two purposes simultaneously: (a) establishing the depth of the problem (why it is not easily solved), and (b) introducing the notation or framework in which the solution will be expressed. The two are not separated.

### Stage 3 — Prior Art Survey (~2–3 paragraphs)

The prior art is surveyed with high citation density. BHL introductions contain 5–10 citations in a single sentence when establishing that a result is known across the field. The prior art section has a specific rhetorical goal: establish what *has* been done and *where the gap lies*. The gap is always stated explicitly — not implied. The last sentence of the prior art survey typically begins "However..." or "But..." and names the specific gap that the paper fills.

Citation placement in BHL prior art: citations always appear in parentheses at the end of a sentence or clause, never at the start. Never "As Jones showed [1], ..." — instead "the conversion probability for the Gertsenshtein effect is known analytically [Jones 1962; Boccaletti 1970; ...]."

### Stage 4 — The Hook (~1 paragraph)

The hook announces the paper's contribution and maps the remaining sections. BHL hooks are specific and direct: "We present X, which does Y. We apply it to Z and find W." No vague capability claims — the hook states the actual contribution. The structural map closes the introduction: "In §2 we... In §3 we... In §4 we..."

**For the Barker solo sub-type**, the introduction runs 5–7 stages rather than 4, with an additional paragraph contrasting the present approach with the closest prior method. This makes solo papers slightly longer in the introduction (~1000 words vs ~600–800 for joint letters). The TIDAL report should use the 4-stage structure (joint-paper style) to stay within the 800-word introduction budget.

---

## (d) How BHL Papers Introduce Gauge Theory Frameworks

The theory section of a BHL gauge-gravity paper — whether on PGT, Weyl gauge theory, or Einstein–Cartan theory — follows a strict pedagogical sequence. This sequence is so consistent across the corpus (2406.12826, 2101.02645, 2309.14783, 2005.02228, 2008.09053) that deviating from it signals a pedagogical or structural error.

**Step (0): Index and notation conventions.** Before any tensor definitions, state the index convention: lowercase Greek for coordinate/holonomic indices ($\mu, \nu, \rho, \sigma$), uppercase Roman for Lorentz/frame indices ($A, B, C$) if the vierbein formalism is used. BHL papers declare this in the first named paragraph of the theory section ("Conventions" or as part of the Riemannian baseline). For the TIDAL report, this convention paragraph should also declare the signature $(+,-,-,-)$ and the Planck mass definition $\kappa^2 = 16\pi G$.

**Step (a): Recall the Riemannian baseline (1–2 sentences).** Even in papers where the reader knows GR, the Riemannian baseline is stated. Its purpose is not pedagogy but notation: "We work with a metric $g_{\mu\nu}$ of signature $(+,-,-,-)$; the Levi-Civita connection is $\mathring{\Gamma}^\sigma{}_{\mu\nu}$." This sentence is purely notational. Note: some BHL papers work in the metric (second-order) formulation, others in the Palatini (first-order, metric + connection independent) formulation. State the choice explicitly here: "We adopt the Palatini formalism in which both the metric and connection are independent fields."

**Step (b): State the generalisation (1–2 sentences).** "In Poincaré gauge theory, the connection $\Gamma^\sigma{}_{\mu\nu}$ is promoted to an independent field." State the physical meaning of the generalisation and name its gauge symmetry. No more than two sentences. If the vierbein formalism is used, note it here: "The connection is most naturally expressed via the vierbein $e^A{}_\mu$ and spin connection $\omega^{AB}{}_\mu$."

**Step (c): Define the two field strengths (paired adjacent equations).** The curvature tensor $R^\sigma{}_{\mu\nu\rho}$ and the torsion tensor $T^\sigma{}_{\mu\nu}$ are defined in paired adjacent equations with nearly identical syntactic structure — the reader sees them as a system, not as independent definitions. "The field strengths of the gauge connection are the curvature (i) and torsion (ii)..." followed by the two equations.

**Step (d): State the most general quadratic action.** Write the action with named coupling constants. State explicitly that this is the most general quadratic action and that special cases correspond to known theories. Do not derive which special cases — just assert the claim. The coupling constants are given generic labels ($\alpha_n$, $\beta_n$, $\gamma_n$ in BHL notation) with brief physical glosses. After the action equation, one prose sentence glosses the coupling constant groups: "The $\alpha_n$ couple the three irreducible torsion-squared invariants; the $\beta_n$ the curvature-squared invariants."

**Step (e): Introduce the irreducible decomposition (in a named paragraph or subsection).** The irreducible torsion decomposition (traceless tensor ${}^{(1)}T$, trace vector ${}^{(2)}T$, axial pseudovector ${}^{(3)}T$ following Hehl 1976) is always introduced *after* the action, never before. The reason: the decomposition is used to write the action in canonical form; introducing it before the action inverts the logical order.

**What BHL does NOT do:**
- Open the theory section with a review of Yang-Mills gauge theory or Lie group theory. The reader is assumed to know gauge theory.
- Define the Lorentz group or the Poincaré algebra from scratch. One sentence on the gauge structure suffices.
- Introduce index conventions mid-derivation. All conventions are declared in the notation paragraph at the start of §1.

---

## (e) Equation Handling Conventions

BHL papers handle equations with consistent conventions that are worth recording explicitly:

**The sandwich rule**: every displayed equation is embedded in a sentence — the sentence before states what the equation expresses, the sentence after draws a consequence or names a quantity. "The Lagrangian for Einstein–Cartan theory reads [equation], where the coupling constants $\alpha_n$ are real parameters." Orphan equations (displayed without surrounding prose) do not appear in BHL papers.

**Introducing equations with verbs**: BHL uses a small vocabulary of standard verbs to introduce equations, each carrying distinct semantic weight:
- "reads" — for a Lagrangian, equation of motion, or defining relation (the equation *is* the object)
- "is given by" — for a derived or computed quantity
- "is defined by" — for a notation convention
- "takes the form" — when an equation appears after a transformation or substitution
- "reduces to" — when a general expression specialises to a special case

Avoid "can be written as" (passive voice, implies arbitrariness) and "is simply" (patronising).

**Paired definitions**: when two related objects are defined (curvature and torsion; kinetic and potential terms; tensor and axial sectors), they appear in adjacent equations with the same syntactic template. The reader perceives them as a dual system.

**Crediting algebra without displaying it**: BHL does not display intermediate algebra steps. After setting up a calculation, the result appears with a phrase such as "after expanding to quadratic order in perturbations, one finds..." or "a direct calculation yields...". This pattern appears in every paper in the corpus. The purpose is to credit the effort without spending page budget on it.

**Closure sentences after a calculation**: after presenting a derived result, BHL writes a closure sentence that draws the interpretive payoff. This is either teleological ("this is the linearised equation of motion that we will solve in §4") or interpretive ("the coupling constant $\alpha_3$ is therefore the only free parameter in the conversion amplitude"). Never leave a displayed result hanging.

**Label discipline**: every equation that will be cited — even once — has a semantic label: `\label{QuadraticAction}`, `\label{TorsionTransform}`, not `\label{eq:17}`. BHL papers are consistent on this; equations without labels are either unreferenced or inline.

---

## (f) Theory Section Sub-Structure

The theory section of a BHL PGT paper is almost always divided into three sub-movements, regardless of whether they are labelled as subsections:

1. **The field content and action** (what fields, what symmetry, what quadratic Lagrangian). This is Steps (a)–(d) from §(d) above.

2. **The decomposition and degrees of freedom** (irreducible torsion, linearisation, perturbation ansatz). This is Step (e) plus the perturbative expansion around a background.

3. **The observable or kinetic matrix** (what the equations of motion look like in the chosen basis; what the propagating modes are). This is the bridge to the results section.

In long papers (2101.02645, 2309.14783), these three movements are labelled subsections. In letters (2406.12826, 2311.11790), they are named paragraphs within a single section or unlabelled prose blocks.

**Linearisation sub-movement** (present in all papers with perturbation calculations):
1. State the background explicitly: flat metric $\eta_{\mu\nu}$, uniform static magnetic field $B_0$ in a named direction.
2. Write the perturbation ansatz: $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, $A_\mu = A_\mu^{(0)} + \delta A_\mu$. Name the expansion parameter.
3. State the order-counting rule immediately after the ansatz: "We work to linear order in $h_{\mu\nu}$ and $\delta A_\mu$, and to first order in $B_0$." This pattern appears verbatim in multiple BHL papers on linearisation (2510.17094, 2510.08201) and is a genre requirement, not a stylistic choice.
4. State what the quadratic Lagrangian produces: "after expanding to quadratic order, the kinetic matrix is..."

---

## (g) Results Section Conventions

BHL results sections vary more than theory sections, but six consistent patterns appear across the corpus:

**Pattern 1: Lead with positive validation.** Even when the main scientific contribution is a null result or a constraint, the first result presented is always a known result reproduced by the new method. Papers that introduce a new calculational framework (2406.12826, 2101.02645, 2510.08201) always begin by showing it recovers known spectra before applying it to new theories. This is not optional: it establishes credibility for the novel results that follow.

**Pattern 2: Null-as-precision-finding.** When a calculation yields a null (no amplification, no mixing, no ghost in a given sector), BHL frames it as a quantitative constraint rather than an absence: "we constrain the conversion probability to $P < 10^{-5}$ across the full parameter space" rather than "we found no conversion". The Dandoy framing (2406.17853) is the closest published precedent for TIDAL's null results: converting upper bounds from absence-of-signal into quantitative exclusion regions is the genre's standard move.

**Pattern 3: Competing-effects or cancellation framing.** When a mechanism that might produce an effect is present but neutralised by a compensating term, BHL states both the mechanism and the cancellation explicitly. The positive statement ("vector torsion couples to the photon") comes before the cancellation ("but the coupling vanishes at linear order because..."). Plane-wave cancellation arguments (2510.17094) follow this structure precisely.

**Pattern 4: Ghost and tachyon reporting.** The standard BHL ghost-reporting pattern (observed in 2005.02228, 2510.08201, 2510.17094) is: (i) state the no-ghost condition as a closed inequality on the coupling constants; (ii) show whether the condition is satisfied at the model point under consideration; (iii) state the consequence. Do not call ghost instability a "failure" — call it a "diagnosis" that "identifies the [sector] as ghost-unstable at [parameter point]."

**Pattern 5: Section-internal roadmapping.** Long results sections (2101.02645, 2406.12826 §IV–V) open with one sentence mapping their subsections: "§4.1 validates the spectrum against known results; §4.2–4.3 present the new propagating-torsion sectors; §4.4 gives the ghost-free conditions." This is identical to the Introduction's structural map but at section rather than paper level.

**Pattern 6: Figures carry the primary result load in parameter scans.** Prose accompanying a parameter scan figure identifies the key qualitative feature visible in the figure (the boundary of the stable region, the torsion-independent line) without repeating numbers that are already readable in the plot. The caption supplies the quantitative reading; the prose supplies the interpretation.

---

## (h) Discussion and Conclusion Conventions

The discussion section of a BHL physics paper has a recognisable three-move structure. This structure is distinct from the results section — it does not repeat results, it interprets them.

### Move 1: Restatement of Findings (~1 paragraph)

Plain language, no equations, no hedging. "We have shown that...", "we found that...", "our main finding is...". Results that were presented with hedges in the results section are stated confidently here — the hedge appeared to signal epistemic status during the calculation; the conclusion confirms the finding. Do not retrospectively hedge results that were reported with confidence in §4.

### Move 2: Implications for the Broader Programme (~1–2 paragraphs)

What the findings imply for the broader theoretical programme, not just the paper's own results. This is where the "three constructive paths" framing belongs in the TIDAL report — not in §4. The standard BHL move is: "If our null result is robust, it implies X for the theory space. However, three avenues remain open: (i) non-minimal coupling, (ii) ghost-free kinetics in the [sector], (iii) cubic extensions." Enumerate the avenues, do not develop them.

### Move 3: Specific Next Steps (~1 paragraph)

Concrete, specific, actionable. Not "future work should explore..." but "the natural next step is to apply TIDAL to the full non-minimal PGT Lagrangian (App. B), which includes the torsion–curvature cross-term excluded by construction from the present campaign." One theory, one observable, one calculational step. This last paragraph is the invitation to readers who want to build on the work.

**BHL discussion anti-patterns** (seen in non-BHL comparison papers but absent from the corpus):
- Do not enumerate all limitations in a bullet list. Acknowledge limitations inline, in the context of the finding they limit.
- Do not end with "we hope that..." (passive aspiration). End with "the next step is..." (active direction).
- Do not put derivations or equations in the discussion. Equations belong in §2 or the appendices.

---

## (i) Appendix Conventions

BHL papers use appendices extensively and purposefully. Three functions appear consistently:

**Function 1: Proofs and extended derivations.** Material that is essential for the result's validity but would interrupt the main narrative. The main text states the result and cites the appendix; the appendix gives the derivation. This is the primary use across all long BHL papers (2101.02645, 2309.14783).

**Function 2: Tabular primary data.** Classification tables, complete spectra, enumeration of cases. When a paper's contribution is partly an enumeration (all ghost-free theories, all PGT coupling-constant combinations), the table goes in an appendix with a brief prose interpretation in the main text.

**Function 3: Computational infrastructure and reproducibility.** Code availability statements, numerical methods, software versions. This function appears only in post-2020 papers (2406.12826, 2510.08201, PSALTer, Hamilcar). It is a genre evolution — pre-2020 BHL papers rarely had computational appendices; post-2020 papers treat reproducibility as a first-class obligation.

**Appendix titles and cross-referencing**: appendices carry semantic titles that state their content, not labels: "Appendix A: Irreducible decomposition of the torsion tensor" rather than "Appendix A: Technical details." Cross-references from the main text use `\cref{app:label}` with a brief one-clause description of what is being deferred: "...as we show in App. A, the irreducible decomposition under the Lorentz group is..."

---

## (j) Voice, Hedging, and Rhetorical Modes

### Voice

BHL theoretical papers use first-person plural throughout, including in solo-authored work (2205.13534, 2510.08201). "We show", "we find", "we confirm" — even when the result is entirely the author's own computation. The passive voice appears only in three contexts: (a) historical attribution ("the Gertsenshtein effect was predicted in..."), (b) problem diagnosis ("the trace sector is ghost-unstable"), (c) standard mathematical descriptions ("the action is invariant under...").

Never use "it can be seen that", "it should be noted that", "the results indicate", "it is clear that", or "obviously". These are forbidden in STYLE_GUIDE.md and are consistently absent from the BHL corpus.

### Hedging

BHL papers hedge only when genuinely uncertain. The pattern is: finding first, hedge second. "The trace-channel amplification is suppressed — pending a full nonlinear analysis." Front-loading with caveats ("although we cannot rule out... we find...") is not BHL practice. 

Strategic hedges appear in the discussion, not in the results section. A result is stated in §4 with the confidence the calculation supports; uncertainty about whether the result generalises is addressed in §5.

### Footnote culture

BHL papers use footnotes for substantive technical asides that specialists need but which would interrupt the narrative. The footnote culture is visible in 2101.02645 (extended notes on subtleties in the Hamiltonian analysis), 2406.12826 (technical caveats on the conformal identification), and Barker's solo papers (ironical or classical cultural references). Footnotes in BHL papers are typically 1–4 lines. A footnote that runs to 8+ lines signals content that belongs in the body or an appendix.

---

## (k) Conventions Specific to Gertsenshtein and Torsion-EM Papers

The Gertsenshtein-effect literature (gr-qc/0307063, hep-th/0103093, 1812.02675, 2301.02072, 2406.17853) has additional genre conventions that do not appear in the wider PGT corpus:

**The conversion probability formula is always displayed prominently.** In Gertsenshtein papers, $P = \sin^2(\kappa B_0 D / 2)$ (or its generalisation) is displayed as a named, labelled equation in the theory or results section — not buried in prose. This formula is the paper's central technical object. Papers that extend or modify the formula display both the original and the generalisation for comparison.

**The B₀-scaling argument typically appears.** Papers computing conversion probabilities usually establish the $B_0$ scaling of the result, because $B_0$ is the experimental control parameter that readers comparing with observations need. This is not a hard genre requirement — in some papers (2406.17853) the scaling emerges as part of the solution rather than being foregrounded — but when writing §4, including the $B_0$ scaling explicitly (e.g., "$P \propto B_0^2$") aids comparison with observational bounds. If the TIDAL null result is B₀-independent (as expected for the torsion sectors), state this explicitly.

**Wave decomposition as standard technique**: torsion-EM papers decompose gravitational waves into polarisation modes (cross, plus, axial, trace) before computing coupling. The decomposition is standard enough that it is not derived — the modes are named and the decomposition is stated as a known fact with a citation.

**Two Gertsenshtein-paper sub-genres require different prose styles**:

*Pre-2015 calculation papers* (gr-qc/0307063, hep-th/0103093): sparse prose, dense equations, derivation is the result. These papers use compact equation chaining within paragraphs — in a long derivation, sequential equations are introduced with minimal transition prose ("expanding...," "diagonalising..."). This style is not appropriate for the TIDAL report, which is written in the post-2020 narrative register.

*Post-2020 constraint papers* (1812.02675, 2301.02072, 2406.17853): these papers have absorbed the full BHL narrative style — introduction with tension, theory with named conventions, results with figures, conclusion with outlook. The TIDAL report should follow this sub-genre. The post-2020 papers connect to LISA, NANOGrav, and LIGO sensitivity curves and cite at least one observational instrument in the introduction and one sensitivity curve in the results. The report should anchor results to detection-relevant field strengths, even if the tidal result is a null.

**Observational instrument citation as genre norm** (post-2020 papers only): cite at least one future detector (LISA, DECIGO, BBO) and one current bound (LIGO O3, PPTA, NANOGrav) by the end of §1. This positions the calculation within a detectability programme, even for a theory-side null result.

---

## (l) Cross-Paper Evolution Notes

**From pre-2015 to post-2020 BHL papers**: the clearest evolution is in computational transparency. Pre-2015 BHL papers (Barker's earliest with Lasenby) state results without computational infrastructure claims. Post-2020 papers (2406.12826, PSALTer, Hamilcar) include explicit statements about software availability, reproducibility, and computational cost. This is not a change in physics conventions — it is a response to community expectations about open science.

**From early Gertsenshtein papers to recent**: the original Gertsenshtein papers (1962, Boccaletti 1970) are calculation papers with minimal prose. By 2019–2024 (1812.02675, 2301.02072, 2406.17853), the genre has absorbed the BHL narrative style: introduction with tension, theory with named conventions, results with figures, conclusion with outlook. A modern torsion-Gertsenshtein paper reads like a BHL letter, not like a 1960s calculation note.

**From solo-Barker to joint-BHL voice**: Barker's solo papers (2205.13534, 2311.11790) use a more compressed, ironic voice with occasional footnotes of pure cultural content. Joint papers (2101.02645, 2406.12826) smooth out the ironism and adopt a more neutral collaborative register. The TIDAL report is solo-authored but intended for publication as a joint paper — write in the collaborative register (joint-paper voice) while keeping the directness of Barker's solo style.

**Software as first-class contributor**: the Hamilcar and PSALTer papers (2512.25007, 2406.09500) represent a genuine genre innovation: software tools are co-contributors to a physics result, not merely implementation details. This framing is appropriate for TIDAL. The report should assert TIDAL as the enabling instrument, not merely a technical aid. The PGT channel results are not independently derivable without automation — this is the paper's central claim, and it should be stated plainly.

---

## Summary for Agent Use

### §1 Introduction (~800 words, four stages)

**Stage 1 (~1 paragraph)**: Open on the active research tension — gravitational wave detectors probe the gauge-gravity sector, but the torsion channel structure in PGT is unknown. The tension may begin with a brief framing reference to GR before pivoting to the gap. Do not begin "We present..." or "In this paper..." Open on the physics.

**Stage 2 (~1–2 paragraphs)**: Narrow to the specific problem: the Gertsenshtein kinetic matrix for PGT+EM has not been computed across the theory space; manual computation is infeasible for more than one model at a time. Name what makes the problem hard (scale of theory space, entanglement of coupling constants and spectrum).

**Stage 3 (~2 paragraphs)**: Survey prior art with 5–10 citations per broad claim. Cover: original Gertsenshtein/Boccaletti mechanism [cite]; BHL PGT landscape (2406.12826, 2101.02645); PSALTer/Hamilcar as spectral tools (2406.09500, 2512.25007); prior Gertsenshtein-torsion papers (1812.02675, 2301.02072, 2406.17853). Close with "However, no systematic Gertsenshtein channel decomposition for PGT+EM has been performed."

**Stage 4 (~1 paragraph + structural map)**: "We present TIDAL, a symbolic–numerical pipeline that automates this calculation for any quadratic PGT+EM Lagrangian. We apply it to [specific sectors] and find [headline results: null/torsion-independence/ghost]." Then: "§2 establishes the theoretical framework; §3 describes the computational approach; §4 presents the results; §5 discusses implications."

Also cite ≥1 future detector (LISA, DECIGO) and ≥1 current bound (LIGO O3, NANOGrav) to anchor the calculation within the detectability programme.

### §2 Theory (~1200 words)

**First named paragraph ("Conventions")**: declare metric signature $(+,-,-,-)$, index conventions (Greek = coordinate, Roman = Lorentz if vierbein used), $\kappa^2 = 16\pi G$, irreducible torsion labelling (Hehl 1976).

**§2.1 PGT framework**:
1. One sentence: Riemannian baseline (notation only)
2. State whether working in Palatini (first-order) or metric (second-order) formulation — do this explicitly
3. One sentence: promotion of connection to independent field
4. Paired equations: curvature $R^\sigma{}_{\mu\nu\rho}$ then torsion $T^\sigma{}_{\mu\nu}$, identical syntactic structure
5. General quadratic action — label it `\label{QuadraticAction}` — with coupling constants $\alpha_n$, $\beta_n$, $\gamma_n$. Follow with: "The $\alpha_n$ couple the three irreducible torsion-squared invariants."
6. Irreducible decomposition in a named paragraph after the action (never before)

**§2.2 Linearisation**:
1. State the background: flat metric $\eta_{\mu\nu}$, uniform static $B_0$ in the $z$-direction
2. Write the perturbation ansatz explicitly: $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, $A_\mu = A_\mu^{(0)} + \delta A_\mu$
3. *Immediately* state the order-counting rule: "We work to linear order in $h_{\mu\nu}$ and $\delta A_\mu$, and to first order in $B_0$."
4. Derive kinetic matrix: "After expanding the Lagrangian to quadratic order, one finds..." — do not display intermediate algebra

### §4 Results (~1700 words, four subsections)

**§4.1 Boccaletti validation**: Lead with the formula $P = \sin^2(\kappa B_0 D / 2)$. Show figure. State: "TIDAL reproduces the Boccaletti conversion probability to 0.04%." This is the positive result that validates §4.2–4.5.

**§4.2–4.3 Null results**: Framing template: "we constrain the [T1/T2] sector conversion amplification to $A < X$ across the full [$(m_A, \varepsilon)$] parameter plane." Explain the physical mechanism in one sentence: "the [coupling term] vanishes at linear order because [reason]." Do not write "we did not find amplification." Cite the Dandoy precision-bound framing (2406.17853) as the published precedent.

**§4.4 Propagating torsion**: Competing-effects framing. Two sentences: "The tensor (hx) and axial (ax) torsion sectors enter the kinetic matrix identically at quadratic order. Consequently, the conversion probability is torsion-independent in this sector — deviating from the Einstein–Maxwell prediction by less than [bound] across the full stability window $|\delta_1| < 0.005$."

**§4.5 Ghost diagnosis**: Three sentences: (1) state the no-ghost condition as a closed inequality; (2) show it is violated at the model point; (3) state the consequence. Do NOT call this a "failure" — use "diagnosis." Template: "Ghost freedom in the trace sector requires [condition]. At the parameter point [values], this condition is violated. We diagnose the trace channel as ghost-unstable, confirming the linear PGT ghost analysis of [cite]."

### §5 Discussion (~600 words)

**Move 1 (~150 words)**: Restatement — plain language, no equations. "We have shown that... We found that... Our main finding is..." State each result from §4 confidently. Do not hedge retrospectively.

**Move 2 (~250 words)**: Implications. "If the null result is robust, it implies that [physical consequence]. However, three avenues remain open: (i) non-minimal coupling (torsion–curvature cross-terms); (ii) ghost-free kinetics in the propagating-torsion sector; (iii) cubic-order extensions." These are speculative implications — do not develop them quantitatively. They belong here, not in §4.

**Move 3 (~100 words)**: Specific next step. Name one concrete theory (the non-minimal PGT Lagrangian with torsion–curvature cross-term, already in App. B), one observable (non-zero channel amplification if the cross-term coefficient enters the kinetic matrix), and one calculational tool (TIDAL applied to the full non-minimal Lagrangian). Close: "The natural next step is to apply TIDAL to [specific Lagrangian], which [concrete prediction]."
- Do NOT put the three constructive paths in §4; they are implications, not results
