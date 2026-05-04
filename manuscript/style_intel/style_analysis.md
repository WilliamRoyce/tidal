# Style Intelligence: Writing Physics Papers in the Barker / BHL Register

This document is a teaching resource for writing original prose in the register of Will Barker and the Barker–Hobson–Lasenby (BHL) Cambridge group. It is based on close reading of papers from the corpus: 2406.09500, 2506.02111, 2512.25007, 2406.12826, 2303.11094, 2101.02645, 2006.03581, 2003.02690, 2206.00658, 2311.11790, 2309.14783, 2008.09053, 2005.02228, 2510.08201, 2407.09598, 2205.13534; and from the Gertsenshtein literature for results-section conventions: 2510.17094 (Tomomatsu et al., not Barker-authored), 2406.17853 (Dandoy et al.), 2301.02072 (Domcke et al.). Each pattern below is described analytically. One short attributed phrase illustrates the pattern where it anchors the analysis; further quotation is avoided.

---

## (a) Abstract Structure and Length

BHL abstracts are substantive narrative paragraphs, typically 150–220 words for journal-length papers and 80–130 words for letters. They do not use bullet lists. The structure is: one sentence establishing the physical context or problem, one or two sentences stating what the paper presents or proves, and a closing sentence on scope, limitations, or promise of future work. Crucially, the abstract commits to a finding — even a null result is stated boldly. The gravitational confinement paper (2303.11094) opens its abstract with "The short answer is *probably no*", immediately signalling the paper's conclusion. Software papers add one sentence identifying the tool's computational advantage (parallelisation, automation of a previously manual procedure) and one on scope limitations ("The initial release allows for parity-preserving operators…").

**Abstract closing sentence convention**: the final sentence is a definitive statement or emphatic characterisation, not a promise. The conformal PGT letter (2406.12826) closes with a phrase equivalent to "vector torsion is back on the menu" — a rhetorical flourish that is nonetheless a firm claim. The closing sentence may be colloquial in register as long as it is factually definitive. The abstract never ends with "this will be explored in future work" or a limitation.

**What to avoid:** An abstract that spends more than a third of its words on background rather than on what the paper does or finds.

---

## (b) The Barker Introduction Scaffold: A Four-Stage Architecture

Barker's introductions, whether solo or collaborative, follow a recognisable four-stage structure that can be described as: Modern Frontier → Deep Problem → Inventory of Prior Work → Hook into the Paper.

**Stage 1 — Modern Frontier.** The introduction opens by placing the reader at an active research boundary, not at the beginning of history. The opening paragraph of the HiGGS/supercomputers paper (2206.00658) begins at the landscape of competing modified gravity models, and immediately characterises it as "problematic by dint of their heterogeneity". There is no warm-up paragraph on Newton or Einstein. The frontier is stated as a living tension, not a fait accompli.

**Stage 2 — Deep Problem.** The frontier gives way to a specific technical or conceptual pathology that makes progress hard. This is not just "X is difficult" but "X is difficult *because* of Y, and Y has a structural origin." In the HiGGS paper this is the strong coupling problem: not merely that it is bad, but that the Hamiltonian analysis which would reveal it is itself too expensive to carry out manually. The problem is thus nested — you cannot diagnose the disease without a tool that does not yet exist. In the PGT cosmology paper (2003.02690) the deep problem is that the parameter space is vast and the unitarity conditions (which restrict it to viable theories) are entangled with the cosmological predictions in ways that are not transparent.

**Stage 3 — Inventory of Prior Art.** After establishing the problem, Barker surveys what has been tried. This survey is dense with citations but not exhaustive: it identifies the key landmarks and the specific gap that remains. The inventory is structured to show that prior work either (i) made progress on a special case that misses the general problem, or (ii) developed the tools at the wrong level of abstraction, or (iii) required heroic manual labour that does not scale. The PSALTer paper (2406.09500) notes that "several such applications can be found in the literature, which are particularly impressive because they appear to have been computed manually" — a phrase that simultaneously credits predecessors and underlines why automation is needed.

**Stage 4 — Hook.** The fourth stage introduces the paper's contribution, often beginning with a bold "In this paper…" or "We present…" paragraph. It does not hedge. The contribution is framed as closing the gap identified in Stage 2. In letters, Stage 4 merges with a brief structural paragraph ("In \cref{X} we … In \cref{Y} we … Conclusions follow…"). In full-length papers, structural and conventions paragraphs follow separately.

This four-stage pattern is so consistent across solo Barker and BHL joint papers that learning to recognise it is more important than any other single style lesson. If the reader of a draft cannot immediately identify which paragraph belongs to which stage, the introduction needs structural revision.

---

## (c) Section Structure Conventions

BHL papers prefer a small number of long sections over many short ones. A typical 30–40 page paper might have three to five top-level sections: Introduction, followed by one or two technical sections, followed by Conclusions, followed by appendices. Within sections, the work is structured through subsections and then — most distinctively — through *named paragraphs* using the `\paragraph*{Topic}` command. Each named paragraph begins with a bold label that functions as a micro-heading: "Particle spectra", "The algorithm", "Demand in gravity". This paragraph-naming convention allows very fine-grained organisation without proliferating section numbers.

**Three-level hierarchy**: long BHL papers (2101.02645, 2309.14783) use a consistent three-level structure: section → subsection → named paragraph. The reader navigates at the section level, reads at the subsection level, and absorbs detail at the named-paragraph level. Named paragraphs within a subsection are typically 4–6 sentences each, long enough to make a point and no longer. This creates the characteristic visual rhythm of BHL writing: a sequence of named moves, each self-contained, separated by bold headings, with no undifferentiated prose block running more than eight to ten sentences without a paragraph break.

Letters (2406.12826, 2311.11790) push this to the extreme: they have no subsections at all, only named paragraphs. Everything from the introduction through the conclusions is delivered via paragraph headers. This is only viable because the headers do real organisational work; without them the paper would be unnavigable.

**Paragraph length norm**: BHL paragraphs — whether named or unnamed — consistently run 3–8 sentences. A one-sentence paragraph signals a transition, not a finding. A paragraph exceeding 10–12 lines signals that two separate points have been merged; split it. For the MSci report, aim for 5–7 sentences per named paragraph in the theory sections and 4–5 in the results section.

The practical implication for MSci writing: use subsections to delineate major logical breaks, name your paragraphs within each subsection, and keep each named paragraph to a single point.

---

## (d) Equation Handling: Embedding in Prose and Label Naming

Equations in BHL papers are deeply embedded in the prose. The standard pattern is: a complete grammatical sentence ending with a colon or a comma, followed by the displayed equation, followed immediately by a sentence that interprets or uses the equation. An equation is almost never dropped in without either a preceding or a following prose sentence. The prose sentence before an equation typically says what the equation *means*, not just what it *is*; the prose sentence after it draws an immediate consequence or interprets one of its terms.

Cross-references use `\cref{}` throughout, never "Eq." or "(N)". This ensures that labels are load-bearing: every equation that matters gets a label, and labels are named semantically rather than mnemonically. The 2406.12826 paper uses `\label{TorsionTransform}` for the equation showing how torsion transforms under dilations — not `\label{eq17}`. The label `\label{BasicPositionAction}` appears across multiple papers for the quadratic free action, `\label{MasslessSpectrum}` for the massless propagator analysis, `\label{TorsionTransform}` for gauge-covariance equations. This naming discipline means the labels double as a conceptual index of the paper.

Equations in BHL papers are numbered fully (not "unnumbered for display-only"), and grouped equations use `subequations` environments with `\label{ECTorsion}`, `\label{ECCurvature}` labelling the individual lines. Where multiple equations in a `subequations` block need to be cited together, both labels exist and are combined with `\cref{ECTorsion,ECCurvature}`.

---

## (e) Voice and Hedging Patterns

The default voice is first-person plural ("we present", "we show", "we find") even in solo Barker papers. This is not royal we; it is the conventional physics plural that includes the reader in the intellectual journey. Single-author papers use "we" consistently, not "I". The voice is assertive rather than tentative: the group does not write "it may be possible that" or "results suggest" when a calculation has been done; they write "we show" or "we confirm".

**Finding-first, caveat-second pattern**: when presenting a surprising or counterintuitive result, state the result first, then the caveat. "Our result means that every PGT is conformal and — after a two-decade hiatus — vector torsion is back on the menu." The caveat (the hiatus, the prior negative results) comes second as context, not as a qualification that weakens the claim. This pattern is pervasive and important: front-loading with caveats ("although one might expect...", "while it appears that...") is a different register — not BHL. The caveat serves the finding; the finding never serves the caveat.

Strategic hedging appears in two contexts only. First, when the paper's result is explicitly negative or incomplete: the gravitational confinement abstract uses "probably no" with italics, and the body uses "we are disappointed to find" — the hedge here is honest and carries rhetorical weight because it is so unusual in the register. Second, when claiming generality: "there are no convincing theoretical grounds for excluding them" is a softer framing than "they should be included", preserving intellectual integrity while signalling a position.

**Prohibited constructions:** "it can be seen that", "it should be noted that", "the results indicate". These passive hedges dilute the directness of the BHL register.

---

## (f) Footnote Use: Substantive Asides Only

Footnotes in BHL papers are not parenthetical clarifications or literature nods. They are substantive asides that would interrupt the flow of the main argument but are genuinely interesting in their own right. The HiGGS paper (2206.00658) uses a long footnote at a key result to explain the "catch-22" — that the linear unitarity conditions and the nonlinear ghost charges are entangled in a way that traps the theory — where putting this in the main body would distort the rhetorical structure. Another long footnote (ibid.) discusses the "hic sunt dracones" attitude of the community toward PGT; this has an ironic or classical flavour that belongs in a footnote but would be tonal misfire in the body. **This "ironic footnote" pattern** — where classical literary or historical allusions appear as footnotes because their tone does not fit the body's register — is a distinctive Barker style marker. It signals erudition without interrupting the argument. Aim for at most one such footnote per major section.

In the PSALTer paper (2406.09500), footnotes clarify limitations of Mathematica's built-in pseudoinverse routines and explain when `Method->"Easy"` suffices. These are user-manual style footnotes appropriate to a software paper.

The Hamilcar paper (2512.25007) uses a footnote to note that Hamilcar "completely deprecates" HiGGS, Barker's own earlier package — this is a quiet status update that does not need to be in the main text.

Three rules distilled from the corpus: (1) never put a citation-only thought in a footnote if it can be folded into a sentence in the body; (2) use footnotes for technical caveats that experts need but novices can skip; (3) a footnote that would require more than four lines is probably a parenthetical that belongs back in the body.

---

## (g) Transition Device Patterns

Transitions between named paragraphs are almost never explicit connectors ("furthermore", "in addition", "on the other hand"). Instead, the last sentence of a paragraph announces the next move, and the next paragraph's opening sentence executes it. The PGT conformal paper (2406.12826) transitions from the scale-invariant toy model to the gravitational application with a "Recap and extension" paragraph that uses a numbered list to collect the three prior observations, ending: "We will now connect these three observations by showing that eWGT is the unique scale-invariant embedding of PGT." The next paragraph then delivers exactly that. The connector is the sentence "We will now connect…", which is both a transition and a commitment.

**Subsection opening convention**: when beginning a new subsection, BHL papers use a *recall-and-contrast* opening — a brief reference to the standard or baseline framework before introducing the new content. "Recall that in Einstein's theory, the covariant derivative acting on a vector is simply…" (2101.02645) — then the torsion generalisation follows. This is distinct from within-paragraph transitions; it is a pedagogical re-anchor that reminds the reader where they are in the logical chain.

Within paragraphs, the BHL group uses explicit logical markers when the argument has a two-step structure: "Firstly…", "Secondly…", "Nor shall we…", "Evidently…". These are not filler words; they appear only when the logical structure genuinely has parts that need labelling.

The most distinctive transition device is the **recap-after-complexity**: after a dense calculation or a long derivation, the text stops and issues a plain-English summary. In the gravitational confinement paper this appears after the complex perturbation analysis: the reader is explicitly told what was just shown, and what question this opens. This is both pedagogically effective and rhetorically disciplined — it enforces clarity on the writer as well as the reader.

---

## (h) Framing Nulls and Constraints as Positive Results

The BHL group has a consistent rhetoric for null results and constraint findings. A parameter space that has been ruled out is not "eliminated"; it is *characterised* or *enumerated*. The statement "only spin-zero scalar torsion is allowed to propagate" is not framed as a defeat for PGT model-builders but as a precision result that defines the viable corner of theory space. Similarly, the stability conditions derived from the PSO algorithm are written as "no-ghost condition" and "no-tachyon condition" — positive statements about what the theory is, not what it lacks.

When the paper's own result is null — as in the gravitational confinement paper — the abstract leads with "the short answer is probably no" and the body provides the detailed demonstration. The framing is: "we worked hard to find the effect and did not find it; here is what we did find." The null is thus positioned as a positive contribution to the knowledge base, not as a failure.

For constraints on parameters: the standard phrase is "we exclude" rather than "we are unable to find" or "we cannot allow". Exclusion is an active, positive choice.

---

## (i) "Catch-22" Framing for Paradoxes

A recurring rhetorical device in BHL papers is the identification of structural paradoxes in the theory landscape as "Catch-22" situations. The name is used explicitly in the HiGGS paper: the strong coupling problem in PGT is called a "catch-22" because the unitarity conditions required by the linear theory demand that the nonlinearly activated modes contribute negative energies. You cannot make the linear theory healthy without creating a nonlinear disease; you cannot repair the nonlinear disease without compromising the linear health. The footnote names this explicitly, explaining that the historical objection to strongly coupled PGT modes is their "ghostly character, as inferred by inspection of the signs of squared momenta in the Hamiltonian."

This Catch-22 framing signals that the paradox is structural and not resolvable by parameter tuning — which elevates the difficulty level and motivates the paper's proposed resolution. The device also has the practical effect of giving the referee a memorable name for the core tension, which aids recall during review.

---

## (j) Citation Density Norms

Citation density in BHL papers is high but disciplined. Introductions carry dense citation clusters where the literature survey appears (five to fifteen citations in a single sentence for a broad claim), but the body text is much sparser. A technical paragraph introducing a calculation typically cites the foundational work (one to three references) and no more. Equations themselves are rarely cited inline; citations appear in the surrounding prose.

The PSALTer introduction (2406.09500) demonstrates the pattern clearly: a sentence cataloguing manual computations in the literature carries thirteen citations, but the sentences describing the new algorithm have none — the algorithm is the paper's own contribution. The convention is: cite heavily where establishing that a problem exists or a method was developed elsewhere; cite sparsely where presenting your own work.

Self-citation follows a distinct pattern: prior papers by the same group are cited without fanfare, as factual antecedents. "Previous iterations of PSALTer were used in [refs]" is a matter-of-fact acknowledgement, not a promotional sentence. Software papers consistently cite the underlying dependency stack (xAct, Mathematica, specific sub-packages) with precise references rather than vague "standard tools" language.

One unusual feature: the BHL group sometimes cites classical papers (Maxwell 1865, Einstein 1915, Proca 1936, Klein 1926, Gordon 1926) using their original bibliographic form even in a 2024 software paper. This is not pedantry; it signals that the examples are canonical physics tests, not novelties.

---

## (k) Theory and Calculation Section Conventions

BHL physics papers follow a consistent internal structure within their theory sections, distinct from both the introduction scaffold and the software-paper architecture section.

**Notation-declaration paragraph** (often the first named paragraph in §2): BHL papers always include a notation-declaration block near the start of the theory section, after the introduction. It specifies: metric signature, index ranges (e.g., Greek holonomic/coordinate indices vs. Roman Lorentz/frame indices), differential geometry conventions, and any non-standard abbreviations. The 2101.02645 paper has a dedicated "Conventions" named paragraph at the very start of §2. For the TIDAL report, this paragraph should declare: $(+{-}{-}{-})$ signature, $\kappa^2 = 16\pi G$, uppercase Roman for Lorentz indices, lowercase Greek for coordinate indices, the irreducible torsion decomposition labelling (Hehl 1976), and the notation for the background field $B_0$. This paragraph should be the first named paragraph in §2.

**Coupling constant introduction convention**: BHL papers introduce coupling constants in the action equation first — a complete set of coupling constants ($\alpha_1, \alpha_2, \alpha_3$ for the torsion-squared sector; $\beta_1, \ldots, \beta_6$ for the curvature-squared sector) appears in the general quadratic action. The constants are subscripted sequentially by sector, not by appearance in the paper. After the action equation, a brief prose sentence glosses each group: "The $\alpha_n$ couple the three irreducible torsion-squared invariants; the $\beta_n$ couple the curvature-squared terms." In long papers (2101.02645), a table then lists all coupling constants with their physical interpretations. The constants are never renamed or relabelled in later subsections — use them consistently from first introduction.

**Subsection opening moves.** Every subsection opens with a complete declarative sentence that situates the calculation: either (i) a contrastive recall-and-reference to the standard framework — "Recall that in Einstein's theory, the covariant derivative acting on a vector is simply..." (2101.02645) — or (ii) a direct statement of the first computational step — "The pure R² theory... is defined by its action," (2510.08201). Neither form begins with motivation; the motivation was supplied by the introduction. The subsection opens by *doing*.

**Calculation-level motivation sentence**: immediately before beginning a new calculation *within* a subsection, include a one-sentence statement of why this particular step is needed: "To determine the propagating modes, we linearise the action to quadratic order in the perturbation." This differs from the broader introduction-level motivation — it is local justification for the next algebraic move, ensuring the reader is never confused about why a calculation is being done. This sentence consistently uses the construction "To [purpose], we [method]."

**Equation introduction verbs.** The prose sentence immediately before a displayed equation almost always uses a structural framing verb: "is defined by", "is given by", "reads", "takes the form", "reduces to". The verb choice signals the equation's role: "reads" for Lagrangians and equations of motion; "is given by" for derived results; "is defined by" for notational conventions; "takes the form" after a transformation; "reduces to" for special cases. The vocabulary also includes "is expressed as" and "can be written" in slightly longer derivations. The meaning sentence after the equation uses "where", "here", "this encodes", "this ensures that".

**Paired definitions.** When two tensors are related (curvature and torsion, metric and connection), BHL papers consistently define them in adjacent equations of nearly identical syntactic structure. Torsion always appears after curvature in this pairing — "These are the Riemann-Cartan curvature and torsion tensors" (2101.02645) — reinforcing the parallelism of the field-strength structure. Never define torsion in isolation; always pair it with the curvature it accompanies.

**Crediting effort without showing algebra.** When a calculation is long, BHL papers credit the work with a single honest phrase and give only the result: "After some work, the ADM decomposition of the Ricci scalar follows..." (2510.08201). The reader is told that algebra happened; the algebra does not appear unless it is itself the contribution. This is the correct pattern for §2 of the report — state what was done, give the result, move on.

**Derivation closure sentences.** Derivations end with one of two patterns: (a) a teleological sentence marking completion of a phase — "The canonical action is then written in the standard form" (2510.08201) — or (b) an interpretive payoff that anchors the result to physics — "Not only do we recover... but we see how local scale invariance additionally points to..." (2406.12826). Pattern (a) is used when the next subsection continues the same chain; pattern (b) is used at major section boundaries.

**Linearisation around a background.** The Tomomatsu et al. paper (2510.17094) and 2510.08201 both introduce linearisation with the same move: state the background field (uniform magnetic field, flat metric), write the metric ansatz $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$, then immediately state the order-counting rule ("accurate to $\mathcal{O}(B^2)$"). The order-counting statement is always present and always immediately follows the ansatz. For §2.2 of the report, introduce the linearisation ansatz, state the order of accuracy, and then state what is assumed about the background (uniform $B_0$, flat metric) before proceeding to derive anything.

---

## (l) Results Section Conventions in Physics Papers

The results section of a BHL physics paper (or a Gertsenshtein-field paper) has conventions distinct from both the introduction and the theory section. These patterns are drawn from 2510.17094, 2406.17853, 2407.09598, and 2303.11094.

**Lead with the positive validation before the null.** Every results section in this corpus opens with a calculation that works — a case where the expected physics is confirmed. Tomomatsu et al. (2510.17094) presents the plane-wave solution (which gives exact cancellation, a clean analytic result) before the spherical-wave solution (which gives the non-trivial competition). Dandoy et al. (2406.17853) opens the results with the Gertsenshtein conversion probability formula (already known) before applying it to derive the new constraint. The principle: establish that the tool reproduces known physics before using it to generate new physics.

**Null results as precision findings.** When a calculation yields zero or no excess, the BHL convention — and the Gertsenshtein field convention — is to state the null as a positive characterisation of theory space. Dandoy et al. (2406.17853): "we set stringent constraints on the gravitational-wave strain $h_c$, strengthening current astrophysical bounds by ∼1–2 orders of magnitude." The 276-run null sweep in §4.3 of the report should be presented as "we confirm that the trace-channel amplification satisfies $A < \epsilon$ across the full (parameter) plane" — a precision statement, not a failure report.

**Verb lexicon for constraint and null results** (active verbs only):
- "we constrain [X] to [bound]"
- "we exclude [sector] at [level]"
- "we set bounds on [parameter]"
- "we confirm [X] to [precision]"
- "we place an upper limit of [N] on [quantity]"

Avoid: "we did not find", "there is no signal", "we are unable to detect", "the results indicate no". These passive or negative phrasings dilute the finding. The BHL convention is that even a null result is an active characterisation of theory space.

**Quantify with explicit comparisons.** Results in this field are almost always stated as ratios or improvements relative to prior work or reference calculations. Dandoy et al. quantifies every bound against the prior state of the art. Tomomatsu et al. quantifies the competing effects (focusing vs. conversion) as a ratio. For the report, every numerical result in §4 must include: (a) the value, (b) the reference value it is being compared to (Boccaletti formula, prior PGT computation, or theoretical prediction), and (c) the precision or relative deviation.

**Competing effects and cancellations.** When two effects partially cancel (torsion-independence in the hx/ax sector is an exact decoupling, analogous to the plane-wave cancellation in 2510.17094), the framing pattern is: state the two competing mechanisms, show they cancel to a given order, and then name the surviving effect. The sentence structure is: "The [first effect] and the [second effect] precisely cancel, leaving [the surviving physics]." For §4.4, present hx↔ax torsion-independence as: "the torsion coupling in the tensor sector and in the axial sector enter the kinetic matrix with equal and opposite sign at quadratic order, leaving the conversion probability torsion-independent in this sector." This follows the 2510.17094 plane-wave cancellation framing exactly.

**Ghost instability results.** The R² paper (2510.08201) and the PSALTer papers set the convention for reporting ghost diagnoses: state the ghost condition first as a no-ghost requirement (positive definite residue), then show the condition is violated (residue changes sign), then state the consequence (the mode is a ghost). Never lead with "there is a ghost" — lead with the condition that would have to be satisfied and show it is not. This is the pattern for §4.5.

**Section-internal roadmapping.** Results sections in 2510.08201 and 2407.09598 open with a brief forward map of the subsections: "In Sec. X we derive...; In Sec. Y we show...; In Sec. Z we demonstrate..." This is distinct from the global roadmap in the introduction. For §4 of the report, open with one sentence mapping the four subsections (validation → Stage A null → Stage B null → propagating PGT) before diving in.

---

## Summary: The Barker Register in Five Sentences

The Barker register is assertive, technically dense, structurally transparent, and historically anchored. Every section has a named logical purpose. Equations are embedded in interpretive prose, never dropped unattended. Null results and constraints are framed as precision findings. The introduction moves through four identifiable stages: frontier, deep problem, prior art, hook. A reader who has absorbed these patterns can produce prose in this register without conscious effort.
