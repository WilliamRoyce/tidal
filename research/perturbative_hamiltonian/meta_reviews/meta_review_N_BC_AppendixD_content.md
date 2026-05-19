# Meta-Review N — Blagojević-Cvetković 2018 Appendix D: Actual Content

**Date:** 2026-04-27
**Investigator:** Claude (Opus 4.7, harness session)
**Source:** arXiv:1804.05556v3 ("General Poincaré gauge theory: Hamiltonian structure and particle spectrum"), Phys. Rev. D 98, 024014 (2018)
**Working copy:** `/tmp/bc2018/paper.txt` (extracted via `pdftotext -layout`)

---

## Section 1 — What I could access

| Channel | Result |
|---|---|
| `https://arxiv.org/abs/1804.05556` | Abstract page, metadata only — Appendix D NOT visible |
| `https://ar5iv.labs.arxiv.org/html/1804.05556` | HTML render — content truncated mid-Sec 4 by the WebFetch harness, Appendix D body not delivered |
| `https://arxiv.org/pdf/1804.05556` (curl + pdftotext) | **FULL TEXT** — complete Appendix D extracted verbatim, lines 1613-1684 of `paper.txt` |
| Inspire-HEP recid 1667893 | 58 forward citations, no paper found that explicitly elaborates Appendix D's method |

I have the full verbatim text of Appendix D. Every quote in this review is a literal transcription from `paper.txt`. The PDF confirmed.

---

## Section 2 — Appendix D's purpose and structure

Appendix D is titled **"General construction of H⊥"** (NOT "vanishing critical parameters" — that's its application). It runs 30 lines of body + four equations (D.1-D.6).

**Stated purpose** (verbatim, lines 1614-1616):

> "In this appendix, we discuss the general structure of H⊥, including the case when some of the critical parameters vanish. In a simplified but self-evident notation, the relations that define critical parameters have the following typical form (see Sec. 3)..."

**Structure:**

- **D.1** — defines the typical "if-constraint ⇄ velocity" relation `φ = F V` with 2×2 matrix `F = [[a, b̄], [c̄, d]]`, with `c̄ = κ b̄` proportional to a parity-odd parameter.
- **D.2** — gives the diagonalisation `D = P⁻¹ F P`, where `P = [[−b̄, −b̄], [a−c1, a−c2]]`, valid only when `det P = b̄(c2−c1) ≠ 0`.
- **D.3** — change of basis: `φ' = D V'` with `φ' = P⁻¹ φ`, `V' = P⁻¹ V`.
- **D.4** — generic dynamical-Hamiltonian piece `H⊥ᶠ = φᵀ Q V`.
- **D.5** — explicit form of `H⊥ᶠ` when **only one** critical parameter vanishes (`c2 = 0`).
- **D.6** — corresponding no-ghost condition `det F = ad − b̄c̄ = 0, σ c1 > 0`.

The appendix then has a closing paragraph (lines 1678-1684) titled informally as a "comment on kind of 'non-analiticity'", which is the source of the disputed quote.

**Three cases distinguished** (lines 1654-1670, paraphrasing the structure but reproducing the substance):

1. `c1, c2 ≠ 0` — generic case, recovers the Sec. 4 result.
2. `c1 = c2 = 0` — "rather trivial: both if-constraints φ'n become true constraints that appear in the total Hamiltonian, but H⊥ᶠ = 0".
3. **Only one** critical parameter vanishes (`c2 = 0`, `det F = 0`) — yields the explicit formula (D.5) `H⊥ᶠ = (1/c1)(b̄²q1 + d²q2)(φ'1)² + φ'1(b̄²q1 − adq2)V'2`. The 1/c1 factor is the "typical dependence on the critical parameters", and the V'2 term "can be absorbed into the total Hamiltonian, see [14, 15, 23]". An extra constraint φ'2 then "requires to complete the whole consistency procedure".

---

## Section 3 — The b̄ → 0 limit treatment

This is the closing paragraph, **verbatim** (lines 1678-1684):

> "Now, we have a comment on kind of 'non-analiticity' of the above results. Since the assumption b̄ ≠ 0 ensures the regularity of the matrix P, the diagonal matrix D in (D.2) has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained by taking the limit b̄ → 0 of the generic result. However, since the matrix F for b̄ = 0 is already diagonal, the critical parameters cn can be obtained directly from F. The same conclusion also holds for the form of H⊥ᶠ."

**This is a fully constructive resolution.** The text says: yes, the diagonalisation procedure D.2 itself is singular at b̄ → 0, but **you do not need it** — at b̄ = 0 the matrix F is *already diagonal*, so you read the critical parameters off directly. The same reasoning extends to H⊥ᶠ.

A further confirmation appears earlier in the body (Sec. 3.5, lines 513-516):

> "The generic set of the critical parameters c±(F), F = A, B0, B1, B2, is defined provided the parity odd parameters in F do not vanish, see Appendix D. Hence, the limit of the final expressions c±(F) when these parameters tend to zero is not well defined. However, since in that case F is already diagonal, one can identify c± directly from F."

So the b̄ → 0 case is **not** an unsolved open problem in BC 2018. It is a discontinuous limit of the *parametrisation*, not of the *physics* — and the physics for b̄ = 0 is recovered by stepping outside the b̄ ≠ 0 parametrisation and reading off the (already-diagonal) F directly.

---

## Section 4 — Verbatim-quote verification

The earlier round's "no valid limit / cannot be obtained" quote is **literally in the paper, word for word**. From the paper (lines 1679-1681):

> "the diagonal matrix D in (D.2) has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained by taking the limit b̄ → 0 of the generic result."

The quote is real. **However**, the surrounding context inverts its meaning. Quoting two more sentences (the immediate follow-up, lines 1681-1684):

> "However, since the matrix F for b̄ = 0 is already diagonal, the critical parameters cn can be obtained directly from F. The same conclusion also holds for the form of H⊥ᶠ."

The original investigation's framing reproduced the negative half ("no valid limit") and dropped the constructive resolution that immediately follows ("can be obtained directly from F"). The "no valid limit" sentence is a methodological *caveat* about a particular parametric procedure, not a no-go theorem.

This is the major framing error Review 2 anticipated. It is real and it materially changes what BC 2018 says about the singular case.

---

## Section 5 — Verdict: is this a missed solution?

### **PARTIAL — see qualifications below.**

#### What BC 2018 Appendix D *does* solve (constructively)

For **their specific class of theories** — quadratic Poincaré gauge theory of the form (2.6), with curvature coefficients `b1...b6, b̄1...b̄6` and torsion coefficients `a1...a3, ā1...ā3` — Appendix D supplies:

- A complete Hamiltonian-construction recipe in the generic case (all `cn ≠ 0`).
- An explicit formula (D.5) for `H⊥ᶠ` when one critical parameter vanishes (`det F = 0`), with consistent treatment of the resulting extra true-constraint and the absorbed-into-total-Hamiltonian piece.
- A consistent prescription for the `b̄ = 0` (parity-even) edge: drop the diagonalisation, work with the already-diagonal F directly.
- Closing remark: "An extra constraint φ'2 requires to complete the whole consistency procedure" — they acknowledge that secondary/tertiary consistency must still be carried out, but the kinetic piece is handled.

#### What BC 2018 Appendix D *does not* solve (and is not relevant to TIDAL)

The framework of BC 2018 is **purely 2nd-order in time**. Their `H⊥` is built from the standard Dirac-ADM canonical Hamiltonian for a Lagrangian quadratic in field strengths `T, R`. The "critical parameters" are eigenvalues of small (2×2) coefficient matrices `A, B0, B1, B2` that mix two if-constraints in a single irreducible spin sector. When such a 2×2 matrix becomes degenerate, an algebraic-constraint mode is upgraded to a primary-true-constraint mode. **No higher-derivative kinetic structure ever appears.**

In the TIDAL theory `L = (1/κ²)R̃ + α₁I₁ + α₂I₂ + α₃I₃ + b5·R̃² − ¼F²`, the **`b5·R̃²` term is not a curvature-component coefficient in the BC sense.** R̃ here is the *full Riemann-Cartan Ricci scalar*, so `b5·R̃²` is a quadratic-in-curvature **scalar** that, when expanded in the connection, generates **fourth-order time derivatives** for the constraint-promoted fields h_4, h_7, h_9. At b5 = 0 those fields are algebraic-constraint slots (Lagrange multipliers); at b5 ≠ 0 they propagate via 4th-order kinetic terms.

That is an **Ostrogradsky-type rank-jump**, not the BC-type 2×2-mixing rank-jump. BC's eigenvalue `c± = (1/2)(tr A ± √((tr A)² − 4 det A))` is the eigenvalue of a 2nd-order kinetic block; it has no analogue of the Ostrogradsky doubling that the b5·R̃² term induces.

#### Why TIDAL's b5 is *not* one of BC's critical parameters

Concretely, BC's `b5` is the coefficient of the irreducible curvature component `(5)Rijkl` in Lagrangian (2.6); their text never treats `bn → 0` as a singular limit (the Lagrangian itself stays smooth in `bn`). The "critical-parameter" small quantities are the *eigenvalues* `c± = c±(B0)`, `c± = c±(B1)`, `c± = c±(B2)`, `c± = c±(A)`, which are nonlinear functions of `(b2, b3, b4, b5, b6, b̄2, b̄3, b̄5, b̄6)`. For TIDAL's specific parameter choice the relevant BC critical parameters all sit at generic, nonzero values; the singular thing in TIDAL is the *rank-jump in time-derivative order*, which happens **above** the BC framework, not inside it.

#### Concretely: does BC Appendix D give us a recipe for the b5·R̃² rank-jump?

**No.** Three reasons:

1. **Wrong derivative order.** BC's analysis is built around the standard Dirac-ADM Legendre transform of a Lagrangian containing only first time derivatives of the fields. The Ostrogradsky structure of `R̃²` (where R̃ contains time derivatives of the connection, so `R̃²` contains squared second time derivatives) is outside the framework that produces (D.1).
2. **Wrong constraint topology.** BC's "if-constraint" is `φ_n = 0` becoming a primary constraint when one critical parameter `cn` vanishes. TIDAL's b5 = 0 limit involves *primary constraints disappearing as you move off b5 = 0* (h_4, h_7, h_9 stop being Lagrange multipliers and become propagating fields). The two situations are topologically inverse: BC adds a constraint when c → 0, TIDAL removes a constraint when b5 ≠ 0.
3. **Wrong dimension of the singular matrix.** BC's diagonalisation problem is a 2×2 mixing in each spin sector. The rank jump in TIDAL's perturbative-Hamiltonian reduction is a `det(M) ∝ b5^N` with `N = 6` in Round 1's three-promoted-field analysis — not a 2×2 problem with a single zero eigenvalue, but a multi-block determinant degenerate to non-trivial order in b5.

#### Verdict

**The "25-year-unsolved" framing for the *constraint-promotion barrier in higher-derivative torsion-Gertsenshtein PGT* stands.** BC 2018 Appendix D is a constructive recipe, but for a **different class of theories** (2nd-order quadratic PGT with parity-odd 2×2 mixing). The recipe does not transfer to the b5·R̃²-induced Ostrogradsky rank-jump that is the actual TIDAL blocker.

What was wrong in the original investigation:

- The verbatim "no valid limit" quote was real but **stripped of its constructive follow-up sentence**. As a description of BC's intent, the original framing was misleading.
- The investigation appears to have assumed BC's Appendix D was the closest published analogue, then transferred the negative half of its closing remark into a "25-year-old unsolved problem" claim.

What was right:

- BC's recipe genuinely *does not* solve TIDAL's specific problem. The two singular-limit problems are structurally different, even though both involve a parameter going to zero in a PGT Hamiltonian.
- The claim that there is no published Hamiltonian recipe **for TIDAL's b5·R̃² Ostrogradsky rank-jump** is, to the best of my forward-citation searching, supported. None of the 58 papers citing BC 2018 (per InspireHEP) appears to extend Appendix D to higher-derivative R̃² PGT.

---

## Section 6 — Concrete next steps

### 6.1 Documentation correction (immediate)

In any TIDAL document or memory file that cites BC 2018 Appendix D as evidence of a "25-year-old unsolved problem" or "no valid limit" no-go:

- **Add the constructive follow-up sentence**: "However, since the matrix F for b̄ = 0 is already diagonal, the critical parameters cn can be obtained directly from F."
- **Reframe the citation**: BC 2018 Appendix D solves the *parity-violating 2×2-mixing* singular limit constructively; it does not address the Ostrogradsky rank-jump arising from `b5·R̃²` quadratic curvature.
- **Tighten the "unsolved" claim**: the unsolved problem is *higher-derivative-induced rank-jump in PGT canonical structure*, which is genuinely outside BC 2018's scope. The 2018 paper is not evidence for or against a solution; it is silent on this case.

Files to update:
- `docs/tex/perturbative_reduction_constraint_barrier.tex` — primary doc, will reference BC 2018
- `research/perturbative_hamiltonian/notes/round1_synthesis.md` — likely cites BC
- Any meta-review or plan that propagated the "no-go" framing

### 6.2 Re-check whether BC's case-(3) formula has structural lessons for TIDAL

BC's case (3) (eq. D.5) treats *one* zero eigenvalue: it produces (a) a true constraint, (b) a 1/c1 singular kinetic term, (c) an undetermined velocity that absorbs into the total Hamiltonian.

The TIDAL Round 1 result `det(M) ∝ b5⁶` resembles case (2) (`c1 = c2 = 0`, "trivial: H⊥ᶠ = 0") in that the entire constraint matrix collapses to zero, not just one eigenvalue. Worth checking explicitly whether the TIDAL block matrix at b5 = 0 is in BC's case (2) regime — if so, BC says this case is "trivial: both if-constraints φ'n become true constraints that appear in the total Hamiltonian, but H⊥ᶠ = 0", which would parallel the TIDAL claim that h_4/h_7/h_9 are Lagrange multipliers at b5 = 0 (no kinetic piece).

This is not a *solution* — it is a *consistency check* that BC's classification at least recognises the limit shape we see.

### 6.3 Forward-citation deep dive (optional, lower priority)

The InspireHEP search returned 58 forward citations of BC 2018. I sampled a few; none addresses higher-derivative `R̃²` Hamiltonian directly. If we want to be exhaustive:

- Glavan, Zlosnik & Lin 2024 "Hamiltonian analysis of metric-affine-R²" (arXiv:2311.17459) — different Lagrangian (metric-affine, not torsion-only PGT) but methodologically close. **Worth reading next.**
- Nikjoo & Zlosnik 2024 "Hamiltonian formulation of gravity as a spontaneously-broken gauge theory of the Lorentz group" — degree-of-freedom-vs-parameter discontinuity discussion.
- Barker et al. 2024 (arXiv:2406.12826) — already in TIDAL's reference list; conformal-embedding lift for trace torsion sector.

Estimated effort: half a day per paper. Recommended only if Path A (Krupka-Voicu) and Path B (sector-by-sector, currently active) both stall.

### 6.4 What this finding does *not* change

- The Round 1 `det(M) ∝ b5⁶` result stands.
- The current Path A (Vainberg-Tonti / Helmholtz residue) and Path B (sector-by-sector Bopp-Podolsky) directions remain the project's best paths.
- The supervisor-facing framing "constraint-promotion barrier in higher-derivative PGT lacks a published Hamiltonian recipe" is correct **as long as we are precise about what kind of barrier**.

What this finding *does* change:

- The "BC 2018 said it's impossible" rhetorical move must be retired. BC 2018 said no such thing about higher-derivative PGT; BC 2018 *did* say their parametric diagonalisation has no smooth b̄ → 0 limit, but they immediately gave a workaround. The original investigation mis-cited.
- We should remove BC 2018 Appendix D from the "blockers / negative results" list. It belongs in the "neighbouring-but-distinct case studies" list.

---

## Appendix to this meta-review — Appendix D's full text

For the record, here is the literal Appendix D as it appears in arXiv:1804.05556v3 (lines 1613-1684 of `paper.txt`, reformatted for readability but preserving all sentences and equation labels):

> **D — General construction of H⊥**
>
> In this appendix, we discuss the general structure of H⊥, including the case when some of the critical parameters vanish. In a simplified but self-evident notation, the relations that define critical parameters have the following typical form (see Sec. 3):
>
> φ = F V                                       (D.1)
>
> where φ = (φ1, φ2)ᵀ, F = [[a, b̄], [c̄, d]], V = (V1, V2)ᵀ. Here, φ represents the if-constraints, V are the corresponding velocities, and F is the matrix with eigenvalues c1, c2. Since F is chosen to represent A, B0, B1 or B2, the parameter c̄ is proportional to b̄, c̄ = κ b̄. If b̄ = 0, the matrix F is already diagonal, and the construction of H⊥ is quite simple. When b̄ ≠ 0, which is typical for the parity-violating PG, the matrix F needs first to be diagonalized. The diagonal form D of F is constructed as
>
> D = P⁻¹ F P,   P = [[−b̄, −b̄], [a−c1, a−c2]],   D = diag(c1, c2)        (D.2)
>
> where P is invertible provided det P = b̄(c2 − c1) ≠ 0, and
>
> P⁻¹ = (1/det P) [[a − c2, b̄], [−a + c1, −b̄]].
>
> Left multiplication of (D.1) by P⁻¹ yields
>
> φ' = D V'                                     (D.3a)
>
> where φ' := P⁻¹ φ and V' := P⁻¹ V, or equivalently,
>
> φ'1 = c1 V'1,    φ'2 = c2 V'2.                (D.3b)
>
> To construct the related F-part of H⊥, note that its typical form reads
>
> H⊥ᶠ = φᵀ Q V ≡ φ'ᵀ (Pᵀ Q P) V',   Q = diag(q1, q2)     (D.4)
>
> see Sec. 4. Further discussion depends on the specific values of c1 and c2.
>
> (1) When c1, c2 ≠ 0, Eq. (D.1) implies V = F⁻¹ φ, and H⊥ᶠ = φᵀ Q F⁻¹ φ coincides with the result found in Sec. 4. (2) The case c1 = c2 = 0 is rather trivial: both if-constraints φ'n become true constraints that appear in the total Hamiltonian, but H⊥ᶠ = 0. (3) Finally, when only one critical parameter vanishes (which requires det F = 0), say c2 = 0, then φ'2 = 0 (a new constraint), V'2 remains undetermined and φ'1 = c1 V'1. Hence, (D.4) implies that
>
> H⊥ᶠ = (1/c1)(b̄² q1 + d² q2)(φ'1)² + φ'1(b̄² q1 − a d q2) V'2.       (D.5)
>
> The result can be also expressed in terms of the original if-constraints φn by noting that φ'2 = 0 implies φ'1 = −φ1/b̄. The factor 1/c1 in the first term shows a typical dependence on the critical parameters, known from PG+, whereas the second term, linear in the undetermined velocity V'2, can be absorbed into the total Hamiltonian, see [14, 15, 23]. The presence of an extra constraint φ'2 requires to complete the whole consistency procedure.
>
> In the context of the weak field approximation, the form of H⊥ᶠ in (D.5) determines the no-ghost conditions for the case (3):
>
> det F = a d − b̄ c̄ = 0,   σ c1 > 0                    (D.6)
>
> where σ is the sign of (b̄² q1 + d² q2) and c1 = a + d.
>
> Now, we have a comment on kind of "non-analiticity" of the above results. Since the assumption b̄ ≠ 0 ensures the regularity of the matrix P, the diagonal matrix D in (D.2) has no valid limit for b̄ → 0. Hence, the expressions for cn when b̄ = 0 cannot be obtained by taking the limit b̄ → 0 of the generic result. However, since the matrix F for b̄ = 0 is already diagonal, the critical parameters cn can be obtained directly from F. The same conclusion also holds for the form of H⊥ᶠ.

End of Appendix D.
