# §6 Conclusion — drafting spec (TIDAL MSci thesis)

Binding stylistic and structural specification for the §6 Conclusion section. Distilled from a thorough read of the BHL/Barker corpus, the Gertsenshtein-effect literature, and the existing prescriptive files under `manuscript/` (`STYLE_GUIDE.md`, `style_intel/`, `planning/`). This file is the drafter's reference for length, shape, register, content, and exclusions. It does not duplicate `STYLE_GUIDE.md` or `style_intel/genre_conventions_theory.md` — cross-reference those for general rules.

## Context

The thesis main body currently ends with a short `\paragraph*{Closing remarks}` inside the §4 Discussion subsection (the "B fold" of 2026-05-25). The user has decided to **restore a separate §6 Conclusion section** in line with the original `report_plan.md` order (Appendices → §2 → §4 → §3 → §1 → §5 → **§6** → Abstract).

The thesis is **terminal**: it is a Cambridge Part III MSci dissertation, not a journal letter, and has no companion paper to defer to. The conclusion must close the document standing alone. A later journal version may exist but is irrelevant to how §6 is written today.

## What is the §6 Conclusion, and what is it not?

- §6 is **distinct** from the in-results "Closing remarks" paragraph in `discussion.tex`. The closing-remarks paragraph (~70 w, inside §4 Discussion) is the local sign-off of the Discussion-and-outlook block. §6 sits above all four §4 subsections and the §5 Discussion, and is what an examiner who jumps to the last page reads as a self-contained statement of what the thesis delivered.
- §6 is **not** a re-summary of §5. §5 carries the three constructive paths and the Move 1/2/3 structure in full. §6 names them in one sentence at most.
- §6 is **not** a software roadmap. Plans for TIDAL's future development belong in App E / §3.
- §6 is **not** a deferment to follow-up work. The thesis is terminal; no "companion paper" / "in the sequel we will" framing is admissible. The corpus has examples of that register (2205.13534 closes by deferring to a companion paper) — they are not the right anchor here.

## Corpus surveyed

Eleven BHL/Barker/Gertsenshtein-adjacent papers were read end-to-end for their closing sections: **2406.12826, 2303.11094, 2101.02645, 2510.17094, 2406.17853, 2205.13534, 2311.11790, 2406.09500 (PSALTer), 2512.25007 (Hamilcar), 2206.00658 (HiGGS), 2510.08201 (R² spectrum).** Patterns reported below are synthesised across all eleven; no single paper is used as the anchor.

## Length

The corpus shows a clear bimodal split:

- **Letter-length closings**: 250 w (2311.11790), 357 w (Hamilcar), 550 w (2510.17094), 638 w (2510.08201 Discussion).
- **Long-paper closings**: 900–2100 w (2406.12826, 2205.13534, PSALTer, HiGGS, 2101.02645, 2303.11094, 2406.17853).

For a thesis whose §5 already carries ~600 w of Move 1/2/3 + Scope, the right scale for a *separate* §6 is the letter-length register: **≈ 250–350 w**, single section, no subsections, 3–4 paragraphs of 60–100 w each. This matches 2311.11790 and Hamilcar — the two corpus papers whose role most resembles TIDAL §6 (a brief closing on top of an already-discursive earlier section). It is also the smallest feasible footprint against the 5000 w main-body cap.

If the 5000 w cap proves binding, the §5 closing-remarks paragraph (~70 w) should be removed in favour of §6 — never the reverse.

## Heading

`\section{Conclusions}` with `\label{Conclusions}`. Use the plural ("Conclusions") — the corpus is unanimous on this when the heading is a `\section`. Reserve "Closing remarks" / `\paragraph*{Closing remarks}` for the in-results sign-off only.

## Structural shape (3–4 paragraphs, in order)

The corpus shows two reusable structural skeletons. The TIDAL §6 should fuse them as follows.

1. **Opening: content-first restatement of the frontier and what TIDAL made tractable (≈ 80 w).** Echo §1's frontier tension in one sentence; in a second sentence state what the thesis built / did about it. The corpus's content-first opening (2406.12826: "Poincaré gauge theory strongly motivates spacetime torsion, but…") is preferred over the software-paper "In this paper we have presented X" opener — the latter reads as a software README, not a thesis closing. Frame TIDAL as an *instrument* that made a previously-intractable survey possible, not as the subject of the paper.

2. **Headline finding in plain language (≈ 80 w).** The Boccaletti reproduction is *not* the headline — it is an expected baseline check on a known result and should not be elevated to a thesis result. Restate instead the actual scientific finding(s) of the §4 survey: the structural pattern observed in the kinetic matrix across the surveyed PGT classes, the precision constraint on channel amplification ("we constrain A < X across …"), the ghost diagnosis where applicable. Use the Dandoy precision-bound register ("we constrain", "we exclude"), not absence language ("we did not find"). One claim per sentence; no hedges that §4 itself did not supply.

3. **One-sentence index of the open theoretical directions (≈ 50 w).** Name the constructive theoretical paths that §5 develops in full (e.g. non-minimal couplings, ghost-free propagating sectors, higher-order extensions) in a single sentence, with no in-line citations and no development. §6 acknowledges them and moves on — the actual *primary* next directions live in paragraph 4 below.

4. **The forward-direction closing (≈ 60 w; ends with the closing sentence of the thesis).** This is the operational next step for the TIDAL programme, distinct from the speculative theoretical paths in paragraph 3. Two concrete extensions are already identifiable from the present work and should be named here:

   (a) **Repeating the §4 sector survey on the localised geometry setup** — i.e. re-running every theory class already covered under the plane-wave background under the localised (Phase E, dual-Gaussian + wavepacket) geometry, to test whether the structural patterns observed under the plane-wave background survive a physically realistic localised configuration. This is the highest-priority extension because it directly stress-tests the headline result.

   (b) **Atlas-based direction sampling across the PGT parameter space** — the inference appendix's atlas/joint-prior machinery has the capacity to identify "interesting directions" (regions where the kinematic identity breaks, where amplification re-emerges, or where ghost diagnoses flip) inside the full coupling-constant space, going beyond the one-coupling-at-a-time sweeps of §4. Frame this as the natural way to locate non-trivial channel structure that a coarse sweep would miss.

   The closing sentence of the thesis fuses these two extensions into a single active-direction statement — e.g. "The natural extension is to re-run the present survey on the localised-geometry setup and to apply atlas-based direction sampling across the full PGT coupling space, which together will localise any non-trivial channel structure the present plane-wave sweep cannot resolve." (Paraphrase; do not copy verbatim.) These are TIDAL-internal extensions — the same instrument, applied to richer geometries and a higher-dimensional sampling space — not new theoretical paths. That makes them the right closing content for a thesis whose central claim is that automation makes such extensions tractable.

The 3-paragraph variant is achieved by folding moves 3 and 4 into a single paragraph. Pick the variant that hits the word target most cleanly at draft time.

## Closing-line register — what to imitate, what to refuse

The corpus offers four reusable closing-line registers. Two are admissible for the thesis, two are not.

**Admissible:**

- *Visionary infrastructure statement* (HiGGS / PSALTer / 2406.12826 register). Asserts that the survey/tool changes what is dismissable from the theory space. The strongest analogue for TIDAL because automation is exactly what makes the previously-intractable survey possible. Example register (do not copy verbatim, paraphrase to fit): "A future is then suggested in which [old concern] is no longer valid grounds upon which to dismiss [class of theories]."
- *Structural-finding statement* (2510.17094 register). Closes on the precise scope of an exact identity or constraint observed in the data: "The independence from [parameter] is specific to …" Admissible if the §4 headline turns out to be a clean structural identity (e.g. the hx↔ax kinetic-matrix coincidence as a theorem-shaped result).

**Inadmissible (corpus-attested but not for this thesis):**

- *Companion-paper deferment* (2205.13534 register: "In the companion paper we present …"). There is no companion paper. Do not use this register.
- *Aspirational closing-line wording* ("we hope that …", "it would be interesting to …", "future work should explore …"). It is the *wording register* that is forbidden, not the topic. Talking about future work, extensions, and the next steps for the TIDAL programme is not only fine but mandatory in §6 — and is exactly what paragraph 4 above is for. What the BHL corpus and `STYLE_GUIDE.md §g` / `style_intel/genre_conventions_theory.md §h` rule out is the *passive-aspirational locution*: sentences whose main verb is "hope", "expect", "wish", or whose subject is the impersonal "future work" / "it would be interesting". The same content stated with an active subject and a concrete verb is correct: "The natural extension is to re-run the survey on the localised-geometry setup" is admissible; "We hope that future work will re-run the survey" is not. Verify each future-facing sentence by reading its main verb: if the verb is "hope", "expect", "wish", or the sentence has no concrete agent, rewrite it.

The drafter should write the closing sentence first, in isolation, before drafting paragraphs 1–3 — and revise it until it carries the visionary-infrastructure or structural-finding register without hope/expectation verbs.

## Voice, citations, equations, cross-references

- **First-person plural** throughout ("we have shown", "we constrain", "we present"), even though the thesis is solo-authored. This is the BHL collaborative register and is mandatory per `style_intel/genre_conventions_theory.md §j` and `report_plan.md`.
- **Zero equations.** The entire surveyed corpus places zero equations in `\section{Conclusions}` (one exception: 2406.17853, a constraint paper, displays a single BBN inequality — not relevant to TIDAL). Match the default.
- **Zero figures and zero `\cref{}` to figures.** A reader who jumps to §6 should not need to flip back. Refer to results by name ("the surveyed sector nulls", "the kinetic-matrix coincidence"), not by figure number.
- **Citations: 0–3.** Letter-length closings sit at 0–3 (Hamilcar 1; 2311.11790 2–3; 2510.17094 0). Long closings cite 10–20 — not the target here. Use citations only to anchor the next-step claim. Re-citing Boccaletti, TIDAL-internal references, or anything already cited in §1/§4 is not admissible — those debts are paid earlier.
- **Forbidden words** from `STYLE_GUIDE.md §g` apply with extra force in §6 because it is the most-read section of the thesis: no "it can be seen that", "it should be noted that", "the results indicate", "novel", "state of the art", "obviously", "in conclusion", "to summarise".

## Mandatory content checklist for the drafter

The drafter must, at writing time, ensure §6 contains each of the following:

- One-sentence echo of the §1 frontier tension (the Gertsenshtein channel in PGT was unknown across the theory space and manual analysis was infeasible).
- One-sentence statement of what TIDAL is and what it made tractable, framed as an instrument-of-the-thesis (not as the subject).
- The actual §4 scientific finding restated in plain language, using precision-bound framing — **excluding** the Boccaletti reproduction, which is an expected baseline check rather than a thesis result.
- A single-sentence index of the constructive theoretical paths developed in §5, no citations, no development.
- A concrete forward-direction paragraph naming the two TIDAL-programme extensions: (a) re-running the §4 survey under the localised-geometry setup; (b) atlas-based direction sampling across the PGT coupling space (per the inference appendix).
- A final sentence in the visionary-infrastructure or structural-finding register, with an active forward verb and no aspirational language.

## Forbidden in §6

- Any equation, figure, or `\cref` to a §2/§4/appendix object.
- Any retrospective hedge of a result that §4 stated with confidence.
- Any framing that elevates the Boccaletti validation to a headline finding.
- Any companion-paper / sequel-paper / "in future work we will" framing.
- Passive-aspirational locutions: "we hope", "future work should explore", "it would be interesting to", "in conclusion", "to summarise". (Future-work *content* is mandatory — see paragraph 4 of the structural shape; only the aspirational wording register is excluded.)
- Bullet lists, numbered enumerations, or sub-paragraph headings.
- Subsections; the only valid structure at 250–350 w is a flat 3–4 paragraph block.
- Software roadmap material — that belongs in App E and §3.
- New citations not already in `references.bib`. The default for §6 is **zero citations**; if any are needed they must be drawn from keys already present and already cited earlier in the thesis. Do not introduce §6-only references.

## Open decisions for the drafter

- **The exact §4 headline** to restate in paragraph 2 depends on the final state of the §4.2–4.5 results at draft time. Write paragraph 2 with the placeholder "[HEADLINE: …]" and fill it from §4 last.
- **Whether to cite anything at all in §6.** The default is zero citations — the corpus tolerates this (2510.17094 closes with zero). If a citation is needed, it should anchor either the localised-geometry extension (paragraph 4a) or the atlas-based sampling extension (paragraph 4b) to a concrete prior reference already present in `references.bib`. Mirror §5 attribution choices; do not invent new attribution in §6.
- **Three- vs four-paragraph variant.** Decide by word count at draft time. The three-paragraph variant folds the open-directions index into the next-step paragraph; the four-paragraph variant keeps them separate. Both are corpus-attested.
- **Final closing-sentence register** (visionary infrastructure vs structural-finding). Pick at draft time based on whether the §4 headline is best described as an exact identity (structural-finding) or as a methodological capability claim (visionary). Either is admissible; mixing the two within a single sentence is not.

## Quality bar before §6 is considered done

- §6 is 250–350 w by `texcount`.
- §6 contains zero equations, zero figures, zero `\cref{}` references.
- §6 cites no more than three references, all already in `references.bib`.
- §6 contains no item from the `STYLE_GUIDE.md §g` forbidden-words list.
- The closing sentence ends on an active forward verb; it does not contain "hope", "expect", "future work", "companion paper", "in conclusion", "to summarise".
- The Boccaletti validation does not appear as a headline result.
- A reader who reads §6 cold (no prior pages) understands what the thesis claims and what the natural next investigation is.
- Self-review against `STYLE_GUIDE.md §k` checklist passes the items relevant to a closing section.
- `latexmk -pdf main.tex` still compiles; main-body `texcount` ≤ 5000 w with §6 included.
