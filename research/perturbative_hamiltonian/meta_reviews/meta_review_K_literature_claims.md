# Meta-Review K: Independent Verification of Review 2's Literature Claims

**Date:** 2026-04-27
**Reviewer:** Meta-review agent (Round 5 audit pass)
**Scope:** Audit `research/perturbative_hamiltonian/reviews/review2_literature_interpretation.md`.
**Methodology:** Read each cited paper directly via WebFetch (arXiv abstract page), via ar5iv.labs.arxiv.org HTML rendering, via the arXiv PDF (extracted with `pdftotext`), or local TeX where available; quote verbatim with URLs and line numbers; document inability to access where applicable.
**Posture:** Skeptical of both Review 2 and the original investigation. Verify by reading primary sources.

---

## Executive verdict

Review 2 is a **mixed-quality audit**. Of the six specific claims tested:

| Claim | Review 2's verdict | Meta-review (this) verdict |
|-------|-------------------|----------------------------|
| K1. arXiv:2602.12114 = "Cabo Bizet-Bartocci 2026" is fabricated | Author attribution FABRICATED | **Review 2 invented the arXiv ID.** The original investigation cites arXiv:2601.22007 (Aashish-Saif 2026), which exists. The "2602.12114" citation does NOT appear anywhere in the original Round 1-3 notes. Review 2's own audit is therefore based on a misread arXiv number. |
| K2. Krupka-Voicu has only Definition 1, no Theorem 1 | Correct | **CORRECT** (verified via ar5iv) |
| K3. Barker says trace torsion is Yang-Mills gauge field, not Goldstone | Correct (mostly) | **CORRECT on substance**, with one quote (Quote 3) NOT FOUND in the paper as worded |
| K4. BC 2018 Appendix D quote is unverified / possibly hallucinated | Likely hallucinated | **WRONG. The verbatim quote IS in BC 2018**, lines 2685-2688 of the PDF text. Review 2 was unable to access the appendix and prematurely speculated hallucination. |
| K5. CRZ 2024 handles only parity-even free fields | Correct | **CORRECT** (verified via ar5iv and abstract) |
| K6. AM 2020 cross-validation overstated (different φ kinetic) | Correct | **CORRECT** (verified directly from local TeX) |

**Net assessment:** Review 2 was right on K2, K3, K5, K6. Review 2 was wrong on K4 (the most consequential claim — BC 2018 verbatim quote). Review 2 was right on K1's substantive verdict (the original investigation does not cite a "Cabo Bizet-Bartocci 2026" paper — but that's because the citation Review 2 audited at arXiv:2602.12114 *itself appears to have been invented by Review 2*).

**The single most-cited result of the original investigation — the BC 2018 Appendix D verbatim quote — is GENUINE.** Review 2's accusation of hallucination must be retracted. This restores the credibility of Round 1 Finding 2 / Publication C's historical anchor.

---

## Per-claim verification

### K1. arXiv:2602.12114 attribution

**Review 2's claim:** The original investigation cites arXiv:2602.12114 as "Cabo Bizet-Bartocci 2026" but the actual authors are Chan-López, Martín-Ruiz, Cabrera, Paulin Fuentes.

**WebFetch on `https://arxiv.org/abs/2602.12114` (verified 2026-04-27):**
> Title: "Matrix bordering structure of the Faddeev-Jackiw algorithm: Schur complement regularization and symbolic automation"
> Authors: E. Chan-López, A. Martín-Ruiz, Jaime Manuel Cabrera, Jorge Mauricio Paulin Fuentes

So the paper itself is correctly described by Review 2 — at arXiv:2602.12114 there are no authors named "Cabo Bizet" or "Bartocci".

**HOWEVER**: I searched the entire `research/perturbative_hamiltonian/` directory for `2602.12114`, `Cabo`, or `Bartocci`. No such citation appears anywhere in the original Round 1, Round 2, or Round 3 notes:

```
$ grep -rn "2602\|Cabo\|Bartocci" research/perturbative_hamiltonian/ | grep -v review2
(no matches outside Review 2)
```

The original investigation's actual FJ-relevant citation is **arXiv:2601.22007 (Aashish-Saif 2026)**, found in:
- `notes/round2_agentG_novel_directions.md:484`: "arXiv:2601.22007 — Aashish & Saif 2026, *Stückelberg inspired approach for avoiding singular Hamiltonians in Lorentz violating models of antisymmetric tensor field*"
- `results/novel_directions_assessment.md:91`: same.

I verified arXiv:2601.22007 exists (WebFetch returned: title, authors S. Aashish & Md Saif, abstract about Stückelberg auxiliary vector restoring non-singular constraint matrix).

So either:
1. Review 2 misread "2601" as "2602" and "Aashish" as "Cabo Bizet-Bartocci" (a transcription failure), or
2. Review 2 searched arXiv for "FJ Schur" and arrived at the unrelated 2602.12114 paper, then constructed a fictitious original-investigation citation around it to debunk.

The original investigation's citation 2601.22007 is not the Schur-complement FJ paper Review 2 quotes; it's about Stückelberg auxiliary fields for antisymmetric tensors. So even if Review 2 found a real paper at 2602.12114, the audit target (the alleged "Cabo Bizet-Bartocci 2026" citation in the investigation) appears to be fictitious.

**Verdict on K1:** Review 2's audit is **internally confused**. The substance of the verdict — "the original investigation does not contain a `Cabo Bizet-Bartocci 2026` citation" — is true, but only because the investigation does not contain that citation at all (rather than because it is mis-attributed). The verdict that the original investigation contains a fabrication is **not supported** by the original notes. The fabrication, if any, lies in Review 2's own framing.

**Action needed:** Do NOT propagate Review 2's "fix arXiv:2602.12114 attribution" recommendation. There is no `Cabo Bizet-Bartocci 2026` citation in any Round 1/2/3 note to fix. The actual FJ-related citation is arXiv:2601.22007 (Aashish-Saif), which checks out.

---

### K2. Krupka-Voicu Theorem 1 vs Definition 1

**Review 2's claim:** The KV paper has no Theorem 1; only Definition 1 (Eq. 12). Round 2/3 agents who claim "KV Theorem 1 verified symbolically" are using a misnomer.

**Verification via ar5iv (`https://ar5iv.labs.arxiv.org/html/1406.6646`):**
> "(a) Theorem 1: No 'Theorem 1' is present in this paper.
> (b) Definition 1: Yes. Verbatim: 'The canonical variational completion of a source form ε ∈ Ω_{n+1}^r(Y), is the source form τ(ε) given by the difference between the Euler-Lagrange form of the Vainberg-Tonti Lagrangian of ε and ε itself: τ(ε) = E(λ_ε) − ε.'
> (c) The paper contains only one labelled definition and no theorems, propositions, or lemmas with explicit numerical labels."

This confirms Review 2: there is no Theorem 1.

**Caveat:** Round 2/3 agents' identity `EL(L_VT) − ε = 0` is a re-statement of "ε is variational" once a Lagrangian L_VT for ε exists; given Definition 1, this is a self-consistency identity rather than a deep theorem. Review 2's framing of this as a "trivial corollary" is fair.

**Verdict on K2:** **Review 2 is CORRECT.** The "Theorem 1" label is a misnomer; the verifications are consistency checks of Definition 1 / canonical variational completion identity, not proofs of a theorem. Recommendation: rename "KV Theorem 1" to "KV Definition 1" or "canonical variational completion identity" in any downstream artefact.

---

### K3. Barker 2024 — trace torsion identity

**Review 2's claim:** Barker identifies trace torsion T_μ/3 with the Yang-Mills-type Weyl gauge field V_μ (a gauge field, NOT a Goldstone). The Goldstone (compensator) is a scalar φ. The construction is parity-even only; parity-odd is explicit future work.

**Verification via ar5iv (`https://ar5iv.labs.arxiv.org/html/2406.12826`):**

Quote 1 (line 168 of paper text): **VERIFIED**
> "showing that eWGT is the unique scale-invariant embedding of PGT. This will identify T_μ/3 with the vector B_μ when expressed in scale-invariant variables, and thereby reveal ∂_[μ T_ν] ∂^[μ T^ν] to be a Yang-Mills-type term."

(The paper writes the Weyl gauge field as `B_μ`, not `V_μ`. Review 2 substituted "V_μ" silently — minor cosmetic mismatch but the substance is correct.)

Quote 2 (line 109 footnote): **VERIFIED VERBATIM**
> "In this letter we omit the parity-odd invariants only out of simplicity; there are no convincing theoretical grounds for excluding them"

Quote 3 (alleged at line 297): **NOT FOUND** in the paper text via ar5iv. Review 2 quotes "The compensator is purely gauge, so that the embedding theory is completely indistinguishable from PGT after gauge-fixing" but the WebFetch verifier could not locate this exact wording in the paper.

Quote 4 (line 300): **VERIFIED** (closing remarks):
> "Finally, the Yang-Mills-type actions in [No-go theorem] and [Weyl gauge theory] are restricted to parity-even terms for simplicity: the parity-odd extensions should be considered."

**Substantive question:** Per the WebFetch verifier, "trace torsion T_μ/3 emerges as the Weyl gauge field B_μ when the conformal symmetry is broken to Poincaré symmetry" and "the compensator field φ transforms as φ ↦ e^(-ρ)φ under local dilations and is a scalar with Weyl weight -1." So the compensator/Goldstone is the scalar φ; trace torsion is the gauge field.

**Verdict on K3:** **Review 2 is substantively CORRECT.** Trace torsion in Barker is the gauge field, not the Goldstone. Parity-odd is explicitly future work. One of Review 2's four supporting quotes (Quote 3) was not located in my verification attempts — Review 2 may have paraphrased that one. The other three are verbatim or near-verbatim.

**Implication for Path B-trace:** Round 2/3's "Path B-trace ✅ Established literature (Barker et al. 2024)" framing is wrong. Barker's parity-even construction does not establish a Stückelberg/Goldstone recipe applicable to TIDAL's parity-odd `b5·R̃²`. Path B-trace is open / future work conditional on a parity-odd extension of Barker's eWGT. This is a real and important downgrade.

---

### K4. Blagojević-Cvetković 2018 Appendix D verbatim quote

**Review 2's claim:** The verbatim quote ("the diagonal matrix D in (D.2) has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained by taking the limit b̄ → 0 of the generic result") is unverified, possibly hallucinated. Review 2 could not access Appendix D and speculated the quote was confabulated.

**Verification via PDF download + pdftotext (`https://arxiv.org/pdf/1804.05556`):**

The full PDF was retrieved (`pdftotext` output: 2975 lines). Searching for the quote signature:

```
$ grep -n "diagonal matrix\|valid limit\|cannot be obtained" /tmp/bc2018.txt
2686: assumption b̄ 6= 0 ensures the regularity of the matrix P, the diagonal matrix D in (D.2)
2687: has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained
```

Lines 2685-2692 of the PDF text (within Appendix D, after equation D.6) read VERBATIM:
> "Now, we have a comment on kind of 'non-analiticity' of the above results. Since the assumption b̄ ≠ 0 ensures the regularity of the matrix P, the diagonal matrix D in (D.2) has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained by taking the limit b̄ → 0 of the generic result. However, since the matrix F for b̄ = 0 is already diagonal, the critical parameters cn can be obtained directly from F. The same conclusion also holds for the form of H_⊥^F."

This is **bit-for-bit the alleged "verbatim quote"**.

**Verdict on K4:** **Review 2 is WRONG.** The verbatim quote is GENUINE. Review 2 was unable to access Appendix D via WebFetch/ar5iv (the tool truncates the rendered HTML before the appendix), and from the absence of access constructed a structural-texture argument ("labels reused from the appendix designation, opaque quantitative claim") that resembled hallucination heuristics — but the underlying claim was real all along.

This is a **serious error in Review 2**. The accusation that Round 1 Finding 2 / Publication C's historical anchor is "unverified, possibly hallucinated" must be **retracted**. Review 2's recommendation to "remove the verbatim attribution from all downstream artefacts" if the quote is absent is moot — the quote is present.

**Important secondary observation:** Review 2 quoted the paper's body text saying:
> "Extension of the formalism to include vanishing critical parameters is outlined in Appendix D"

and inferred from this that Appendix D *enables* the limit. The paper actually says BOTH things: Appendix D gives the **construction** for the special vanishing-parameter case (taking a *direct* path via the matrix F at b̄=0, which is already diagonal there), AND simultaneously declares that a **continuous limit** from b̄≠0 to b̄=0 of the *generic* expression does not exist. The two statements are consistent: there is a constraint discontinuity, and Appendix D gives a separate non-limit construction for the b̄=0 case. Review 2 misread the body-text "outline" as implying the limit is smooth.

The agents' broader claim — that BC's analysis exhibits a constraint discontinuity at b̄ → 0 — is **textually supported by the paper**. The historical anchor for Publication C is intact.

---

### K5. Chatzistavrakidis-Ranjbar-Zekoč 2024 — parity-even free-field scope

**Review 2's claim:** Paper handles only parity-even, free fields. Construction is at linearized level. m → 0 is the Goldstone limit (opposite of TIDAL's b5 → 0 ⇔ m_q → ∞ limit).

**Verification via abstract (`https://arxiv.org/abs/2411.16928`) and ar5iv (`https://ar5iv.labs.arxiv.org/html/2411.16928`):**

Abstract:
> "We investigate the concept of tensor global symmetries, featuring conserved currents of mixed symmetry and higher spin Nambu-Goldstone bosons. We develop a Stückelberg mechanism for mixed symmetry tensor fields at the linearized level, focusing on the massive graviton, the massive (2,1) Curtright field and the massive (2,2) field..."

ar5iv full-text verifications:
- "We develop a Stückelberg mechanism for mixed symmetry tensor fields at the linearized level" — **confirmed free, linearized**.
- "All actions use the standard kinetic and mass structures... No parity-odd (Chern-Simons-type) terms appear in the construction." — **confirmed parity-even**.
- "when the mass approaches zero, the fields decouple, and the action for the Stückelberg field represents the action for Goldstone modes" — **confirmed Goldstone limit, NOT TIDAL's b5 → 0 limit**.
- Section 5.1, 3 auxiliary fields for (2,1) Curtright: graviton h_μν, Kalb-Ramond b_μν, 1-form a_μ. Field strength `F̊_{μν|ρ} := T_{μν|ρ} − 2∂_[μh_ν]ρ − 2∂_[μb_ν]ρ + 2∂_ρ b_μν − 2∂_ρ ∂_[μ a_ν]`. **Confirmed**.

**Verdict on K5:** **Review 2 is CORRECT.** The paper is parity-even free-field only at linearized level. TIDAL's parity-odd `b5·R̃²` use case is conditional on an unpublished parity-odd extension. Round 3 Agent J's "Path B-tensor-q ✅" verdict is overstated; honest framing is "first published Stückelberg recipe applicable in principle to PGT tensor-torsion modulo a parity-odd extension that is open research."

This downgrade is correct and should propagate.

---

### K6. Aoki-Mukohyama 2020 cross-validation

**Review 2's claim:** AM 2020's φ has non-canonical kinetic `(∂φ)²/(1+φ²)` with sinh-canonicalisation, lives in a dRGT-tuned bigravity with infinite curvature corrections; the bare `𝒳²` term BREAKS ghost-freedom. Agent F's φ is a generic Lagrange multiplier with `φ²/(2b)` potential. They are not the "same φ".

**Verification via local TeX (`literature/2009.11739/PGTandMG1120.tex`):**

Line 502:
> "However, the naive inclusion of the X² term should break the special structure of the ghost-free theory. The YN ghost(s) must reappear, in general."

→ confirms AM say bare X² breaks ghost-freedom.

Line 526:
> `−(3 M_pl² (1−α) / 4) · |f^a_μ| · (∂φ)_f² / (1+φ²)`

→ confirms non-canonical kinetic `(∂φ)²/(1+φ²)`.

Line 531:
> `varphi = sinh θ`

→ confirms sinh-canonicalisation.

Line 542:
> `−(3 M_pl² (1−α) / 4) ∫ d⁴x |f^a_μ| (∂θ)_f² + S_mass[e,f, sinh θ]`

→ confirms the canonical kinetic is for θ, not φ.

Abstract / line 57:
> "in four dimensions, the absence of ghost at non-linear orders requires an infinite number of higher curvature terms, and these terms can be given by a schematic form `R(1+R/αm²)^{-1} R`"

→ confirms infinite tower of corrections needed for ghost-freedom in 4D.

**Verdict on K6:** **Review 2 is CORRECT.** Agent F's `φ²/(2b)` Lagrange-multiplier auxiliary is not the same as AM's parity-odd scalar with non-canonical kinetic embedded in a dRGT-tuned bigravity. The "AM cross-validation" should be downgraded to "qualitative consistency at linearized order; AM's full ghost-free completion requires an infinite tower of curvature corrections that this single-auxiliary lift does not provide."

---

## New issues discovered (not in Review 2)

### NEW-A. Review 2's own arXiv:2602.12114 misidentification

As detailed under K1, the citation "arXiv:2602.12114 = Cabo Bizet-Bartocci 2026" does **NOT appear anywhere in the original investigation notes**. Review 2 introduced this citation as the audit target without textual support. The original notes cite arXiv:2601.22007 (Aashish-Saif 2026), which is real and correctly attributed.

### NEW-B. Review 2 over-confidently inferred hallucination on K4

Review 2's escalation of "I couldn't access Appendix D" to "the verbatim quote is unverified and likely hallucinated" relied on circumstantial structural arguments (label coincidence, opaque quantitative claim). These heuristics misfired here. Future audits should treat "I couldn't access" as inconclusive rather than as evidence of fabrication. Direct PDF retrieval + `pdftotext` should be the default attempt before any hallucination claim.

### NEW-C. Review 2's K3 Quote 3 not located

Review 2 attributed "The compensator is purely gauge, so that the embedding theory is completely indistinguishable from PGT after gauge-fixing" to Barker line 297. My ar5iv verification did not locate this passage. It may be paraphrased rather than verbatim. The substantive verdict (compensator is the scalar Goldstone) is supported by other paper text, so this does not affect K3's overall verdict, but Review 2's verbatim attribution standard should match what it demands of the original investigation.

### NEW-D. Review 2's "diagonal matrix D in (D.2) … structural texture of confabulation" inference (paragraph 524-532)

The structural-texture argument ("labels reused from the appendix designation") was wrong. The matrix really is called D, in equation D.2, in Appendix D — that's exactly because BC chose internally consistent labelling. Coincidence-of-letters is not evidence of confabulation; it's evidence of consistent notation. This kind of fluency-based hallucination heuristic is unreliable at the literature-verification level.

---

## Reliability assessment of Review 2

**Where Review 2 is reliable:**
- Direct quote verification on accessible papers (K2 KV, K3 Barker, K5 CRZ, K6 AM via local TeX) is sound.
- The substantive structural critiques (Path A "Theorem 1" misnomer, AM cross-validation overstatement, CRZ parity-even-only) are well-supported.
- The downgrade recommendations on Path B-trace (K3) and Path B-tensor-q (K5) reflect genuine issues with Round 3 framing.

**Where Review 2 is unreliable:**
- K4 (BC 2018 Appendix D): Review 2 escalated inability-to-access into hallucination accusation. The quote is genuine. This is the **single most consequential error** in Review 2 because it would invalidate the historical anchor for Publication C if propagated.
- K1 (2602.12114 Cabo Bizet-Bartocci): Review 2 introduced an audit target that does not appear in the original notes. The claim that the original investigation contains a fabricated citation is not textually supported. The original investigation's actual FJ citation (2601.22007 Aashish-Saif) is real and correctly attributed.

**Net assessment:** Review 2 is **roughly two-thirds reliable**. The structural critiques (K2, K3, K5, K6) are sound and should propagate. The two HIGH-severity claims (K1 fabrication, K4 hallucination) are themselves either misattributed or wrong. The pattern — getting individual paper verifications right when access is granted, overshooting when access is denied — suggests Review 2's confidence calibration was not appropriate for the partial-information regime.

---

## Recommendations

### Propagate (Review 2 is correct)

1. **Re-title "Krupka-Voicu Theorem 1" verifications to "Definition 1"** (K2). Cosmetic but accurate.
2. **Downgrade Path B-trace from "established literature"** (K3). Barker handles only parity-even; trace-torsion identity in Barker is gauge field, not Goldstone; the parity-odd extension required for TIDAL `b5·R̃²` is open research.
3. **Downgrade Path B-tensor-q from "applies with caveats" to "conditional on parity-odd extension"** (K5). CRZ 2024 handles only parity-even free fields.
4. **Revise AM cross-validation to "qualitative consistency at linearized order"** (K6). AM's φ has non-canonical kinetic, requires infinite tower for ghost-freedom; Agent F's single auxiliary does not match.
5. **Add Voicu 2020 linearity-in-highest-derivatives gate as a separate preflight** (Review 2 §2). Genuinely independent from VT integral convergence.

### Retract (Review 2 is wrong)

6. **RETRACT Review 2's "BC 2018 Appendix D quote is unverified / possibly hallucinated"** (K4). The verbatim quote is genuine, located at line 2685-2688 of the BC 2018 PDF. Round 1 Finding 2 / Publication C's historical anchor is intact.
7. **RETRACT Review 2's "arXiv:2602.12114 attribution as Cabo Bizet-Bartocci is fabricated by the original investigation"** (K1). The original investigation does not contain such a citation; the actual citation is arXiv:2601.22007 (Aashish-Saif), which is correct. The fabrication, if any, is in Review 2's framing, not in the original notes.

### New (not in Review 2)

8. **Add explicit `pdftotext` fallback to literature-audit methodology**. ar5iv is necessary but not sufficient for long appendices. Direct PDF retrieval + text extraction (as used in this meta-review for K4) succeeds where ar5iv truncates.
9. **Treat "could not access" as inconclusive, not as evidence of fabrication** (K4 NEW-B). Adopt the rule: hallucination claims require positive evidence, not just absence of confirmation.
10. **Verify Review 2's K3 Quote 3 ("compensator is purely gauge, completely indistinguishable from PGT after gauge-fixing")** by direct PDF text extraction of arXiv:2406.12826v3. If absent verbatim, downgrade to paraphrase. Substantive K3 verdict unaffected.

---

## Per-claim summary table

| K# | Review 2 claim | Meta-review verdict | Severity if acted on |
|----|---------------|---------------------|---------------------|
| K1 | 2602.12114 fabrication | **WRONG (audit target itself fabricated by Review 2)** | HIGH — would propagate a non-issue |
| K2 | KV has no Theorem 1 | **CORRECT** | LOW (cosmetic rename) |
| K3 | Barker trace torsion = gauge field, parity-odd excluded | **CORRECT (1 quote unverified)** | HIGH — downgrades Path B-trace to open research |
| K4 | BC 2018 quote possibly hallucinated | **WRONG (quote is genuine)** | HIGH — would invalidate Publication C's historical anchor |
| K5 | CRZ 2024 parity-even free fields only | **CORRECT** | HIGH — downgrades Path B-tensor-q to conditional |
| K6 | AM cross-validation overstated | **CORRECT** | MEDIUM — downgrades cross-validation language |

---

## File provenance

- BC 2018 PDF retrieved via WebFetch on `https://arxiv.org/pdf/1804.05556`, 363.9 KB, saved by harness at `/home/vscode/.claude/projects/-workspaces-torsion-gertsenshtein/816ffd39-1460-4b25-aabe-7e9fe6eb46bd/tool-results/webfetch-1777282452280-l9jeug.pdf`. Extracted via `pdftotext` to `/tmp/bc2018.txt` (2975 lines). Verbatim quote at lines 2685-2688.
- AM 2020 TeX read from local `literature/2009.11739/PGTandMG1120.tex` (730 lines). Quotes from lines 502, 526, 531, 542, 57.
- Krupka-Voicu, Barker, CRZ, BC 2018 metadata, 2602.12114, 2601.22007, 2502.17979 verified via WebFetch on respective `https://arxiv.org/abs/...` and `https://ar5iv.labs.arxiv.org/html/...` URLs.

---

## Final note

The meta-review process here surfaced a paradox worth flagging: Review 2 accused the original investigation of confabulation on the basis of pattern-recognition heuristics ("structural texture of a confabulation") when direct access failed. The true confabulation, in the strict sense, was at the audit level: Review 2 introduced a citation (2602.12114 / Cabo Bizet-Bartocci) that does not appear in the audited material. This is not a moral failing of Review 2 — it's a mechanistic consequence of doing literature audits without grep'ing the original corpus for the alleged citation first. Future audits should:
1. Always grep the audited material for the literal citation string before claiming it's misattributed.
2. Always attempt direct PDF text extraction before claiming a quote is hallucinated.
3. Avoid pattern-recognition arguments ("structural texture") for or against fabrication; rely on textual evidence.

Reviewing reviewers is itself a flawed process; this meta-review is not exempt from the same caveats. Recommend a future independent re-verification of K3 Quote 3 and a sanity check on this very document.
