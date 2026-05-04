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

Every equation that will be cited — even once — must have a label. Labels must be semantic, not serial. Use `\label{QuadraticAction}`, `\label{TorsionTransform}`, `\label{GertsenshteinProbability}`, not `\label{eq17}`. Subsections should also carry semantic labels: `\label{sec:LinearisedPGT}`, not `\label{sec:2.3}`. Labels are the paper's conceptual index; name them accordingly.

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

## (i) Self-Review Checklist Before Sharing a Draft

Before circulating any draft section, confirm all of the following:

- [ ] Introduction has four identifiable stages: Frontier, Deep Problem, Prior Art, Hook
- [ ] Every section and named paragraph opens with a declarative sentence (no equation-first openings)
- [ ] Every displayed equation has surrounding prose (one sentence before or after minimum)
- [ ] All cited equations use `\cref{}` (no bare equation numbers, no "Eq. (N)")
- [ ] All equation labels are semantic, not serial
- [ ] No forbidden words from the list above appear in the draft
- [ ] Every figure is cited in text before it appears; caption is self-contained
- [ ] All hedges are in the form "finding first, caveat second"
- [ ] Every footnote passes the "could this go in the body?" test
- [ ] Citation density is high in the prior-art inventory and sparse in the derivation sections
- [ ] Null results are framed as precision findings ("we exclude", "we show that X is not present"), not as failures
- [ ] Abstract ends with a result or a scope statement, not a promise
