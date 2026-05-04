# Style Intelligence: Writing Physics Papers in the Barker / BHL Register

This document is a teaching resource for writing original prose in the register of Will Barker and the Barker–Hobson–Lasenby (BHL) Cambridge group. It is based on close reading of ten papers from the corpus (2406.09500, 2506.02111, 2512.25007, 2406.12826, 2303.11094, 2303.11094, 2101.02645, 2006.03581, 2003.02690, 2206.00658, 2311.11790, 2309.14783, 2008.09053, 2005.02228). Each pattern below is described analytically. One short attributed phrase illustrates the pattern where it anchors the analysis; further quotation is avoided.

---

## (a) Abstract Structure and Length

BHL abstracts are substantive narrative paragraphs, typically 150–220 words for journal-length papers and 80–130 words for letters. They do not use bullet lists. The structure is: one sentence establishing the physical context or problem, one or two sentences stating what the paper presents or proves, and a closing sentence on scope, limitations, or promise of future work. Crucially, the abstract commits to a finding — even a null result is stated boldly. The gravitational confinement paper (2303.11094) opens its abstract with "The short answer is *probably no*", immediately signalling the paper's conclusion. Software papers add one sentence identifying the tool's computational advantage (parallelisation, automation of a previously manual procedure) and one on scope limitations ("The initial release allows for parity-preserving operators…"). The abstract never ends with a promise; it ends with a result or a scope statement.

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

BHL papers prefer a small number of long sections over many short ones. A typical 30–40 page paper might have three to five top-level sections: Introduction, followed by one or two technical sections, followed by Conclusions, followed by appendices. Within sections, the work is structured through subsections and then — most distinctively — through *named paragraphs* using the `\paragraph*{Topic}` command. Each named paragraph begins with a bold label that functions as a micro-heading: "Particle spectra", "The algorithm", "Demand in gravity". This paragraph-naming convention allows very fine-grained organisation without proliferating section numbers. The result is that a long section reads as a sequence of named moves, each self-contained in three to ten sentences, rather than as undifferentiated prose.

Letters (2406.12826, 2311.11790) push this to the extreme: they have no subsections at all, only named paragraphs. Everything from the introduction through the conclusions is delivered via paragraph headers. This is only viable because the headers do real organisational work; without them the paper would be unnavigable.

The practical implication for MSci writing: use subsections to delineate major logical breaks, but name your paragraphs within each subsection so the reader always knows where they are.

---

## (d) Equation Handling: Embedding in Prose and Label Naming

Equations in BHL papers are deeply embedded in the prose. The standard pattern is: a complete grammatical sentence ending with a colon or a comma, followed by the displayed equation, followed immediately by a sentence that interprets or uses the equation. An equation is almost never dropped in without either a preceding or a following prose sentence. The prose sentence before an equation typically says what the equation *means*, not just what it *is*; the prose sentence after it draws an immediate consequence or interprets one of its terms.

Cross-references use `\cref{}` throughout, never "Eq." or "(N)". This ensures that labels are load-bearing: every equation that matters gets a label, and labels are named semantically rather than mnemonically. The 2406.12826 paper uses `\label{TorsionTransform}` for the equation showing how torsion transforms under dilations — not `\label{eq17}`. The label `\label{BasicPositionAction}` appears across multiple papers for the quadratic free action, `\label{MasslessSpectrum}` for the massless propagator analysis, `\label{TorsionTransform}` for gauge-covariance equations. This naming discipline means the labels double as a conceptual index of the paper.

Equations in BHL papers are numbered fully (not "unnumbered for display-only"), and grouped equations use `subequations` environments with `\label{ECTorsion}`, `\label{ECCurvature}` labelling the individual lines. Where multiple equations in a `subequations` block need to be cited together, both labels exist and are combined with `\cref{ECTorsion,ECCurvature}`.

---

## (e) Voice and Hedging Patterns

The default voice is first-person plural ("we present", "we show", "we find") even in solo Barker papers. This is not royal we; it is the conventional physics plural that includes the reader in the intellectual journey. Single-author papers use "we" consistently, not "I". The voice is assertive rather than tentative: the group does not write "it may be possible that" or "results suggest" when a calculation has been done; they write "we show" or "we confirm".

Strategic hedging appears in two contexts only. First, when the paper's result is explicitly negative or incomplete: the gravitational confinement abstract uses "probably no" with italics, and the body uses "we are disappointed to find" — the hedge here is honest and carries rhetorical weight because it is so unusual in the register. Second, when claiming generality: "there are no convincing theoretical grounds for excluding them" is a softer framing than "they should be included", preserving intellectual integrity while signalling a position.

The hedging instrument of choice when presenting a finding that is surprising or counterintuitive is to state the result first, then state the caveat: "Our result means that every PGT is conformal and — after a two-decade hiatus — vector torsion is back on the menu." The caveat (the hiatus, the prior negative results) comes second, as context for the claim, not as a qualification that weakens it.

**Prohibited constructions:** "it can be seen that", "it should be noted that", "the results indicate". These passive hedges dilute the directness of the BHL register.

---

## (f) Footnote Use: Substantive Asides Only

Footnotes in BHL papers are not parenthetical clarifications or literature nods. They are substantive asides that would interrupt the flow of the main argument but are genuinely interesting in their own right. The HiGGS paper (2206.00658) uses a long footnote at a key result to explain the "catch-22" — that the linear unitarity conditions and the nonlinear ghost charges are entangled in a way that traps the theory — where putting this in the main body would distort the rhetorical structure. Another long footnote (ibid.) discusses the "hic sunt dracones" attitude of the community toward PGT; this has an ironic or classical flavour that belongs in a footnote but would be tonal misfire in the body.

In the PSALTer paper (2406.09500), footnotes clarify limitations of Mathematica's built-in pseudoinverse routines and explain when `Method->"Easy"` suffices. These are user-manual style footnotes appropriate to a software paper.

The Hamilcar paper (2512.25007) uses a footnote to note that Hamilcar "completely deprecates" HiGGS, Barker's own earlier package — this is a quiet status update that does not need to be in the main text.

Three rules distilled from the corpus: (1) never put a citation-only thought in a footnote if it can be folded into a sentence in the body; (2) use footnotes for technical caveats that experts need but novices can skip; (3) a footnote that would require more than four lines is probably a parenthetical that belongs back in the body.

---

## (g) Transition Device Patterns

Transitions between named paragraphs are almost never explicit connectors ("furthermore", "in addition", "on the other hand"). Instead, the last sentence of a paragraph announces the next move, and the next paragraph's opening sentence executes it. The PGT conformal paper (2406.12826) transitions from the scale-invariant toy model to the gravitational application with a "Recap and extension" paragraph that uses a numbered list to collect the three prior observations, ending: "We will now connect these three observations by showing that eWGT is the unique scale-invariant embedding of PGT." The next paragraph then delivers exactly that. The connector is the sentence "We will now connect…", which is both a transition and a commitment.

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

## Summary: The Barker Register in Five Sentences

The Barker register is assertive, technically dense, structurally transparent, and historically anchored. Every section has a named logical purpose. Equations are embedded in interpretive prose, never dropped unattended. Null results and constraints are framed as precision findings. The introduction moves through four identifiable stages: frontier, deep problem, prior art, hook. A reader who has absorbed these patterns can produce prose in this register without conscious effort.
