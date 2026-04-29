# Review 2 — Literature Interpretation Audit

**Date:** 2026-04-26
**Reviewer:** Critical-review agent (Round 4 review pass)
**Scope:** Audit Round 1–3 agents' interpretations of the cited literature.
**Methodology:** Read each paper directly via WebFetch / ar5iv / local TeX where
available; compare verbatim paper claims against agent claims; flag
misattributions, over-summarisation, and unstated hypotheses.
**Posture:** Skeptical. Sources of authority are the papers themselves, not
the agent summaries.

---

## Executive verdict

Of the ten papers reviewed, **four are interpreted essentially correctly**
(Lyakhovich 2021, Abakumova-Lyakhovich 2021, Aoki-Mukohyama 2020 — but with
caveats — and Hinterbichler-Saravani 2015), **three are partially over-stated**
(Krupka-Voicu 2015, Voicu 2020, Chatzistavrakidis-Ranjbar-Zekoč 2024), **two
are seriously misinterpreted** (Barker et al. 2024 and Blagojević-Cvetković
2018), and **one is fabricated as cited** (the alleged "Cabo Bizet-Bartocci
2026" paper at arXiv:2602.12114).

The strongest individual finding: **Path A is built on a paper that the agents
do not appear to have read carefully.** Krupka-Voicu 2015 has *no Theorem 1*;
the central object is *Definition 1* (Eq. 12). The agents repeatedly cite
"Krupka-Voicu Theorem 1 verified symbolically" — but there is no such
theorem to verify. The thing that *is* verified (`E_λ(L_VT) - ε = 0`)
is a *consistency check* of the *definition* of canonical variational
completion, not a proof that anything has been "completed". The Helmholtz
test δE = 0 is a *separate* property that must hold *if* the source form is
to be the EL of *some* Lagrangian; KV's Definition 1 produces a completion
*regardless of whether δE = 0*, but in that case `E_λ(L_VT) - ε ≠ 0` and the
result is not a variational completion of ε but of `ε - τ(ε)`. Path A's
"both preflights pass" framing conflates these two facts.

The two seriously misinterpreted papers have direct downstream consequences:

- **Barker 2024** does NOT support Agent C's claim that "trace torsion is a
  Stückelberg/Goldstone for broken Weyl invariance". Barker explicitly
  identifies the trace torsion `T_μ/3` as an honest Yang-Mills-type vector
  (a *gauge field*, not a Goldstone), and a footnote explicitly excludes
  parity-odd terms. Path B-trace is therefore not "established literature"
  — the parity-odd `R̃²` case is named but unsolved by Barker.
- **Blagojević-Cvetković 2018** Appendix D was claimed to contain a verbatim
  quote about "the diagonal matrix D in (D.2) has no valid limit for b̄ → 0".
  I could not access Appendix D in full via WebFetch, but the paper's
  explicit framing of Appendix D ("Extension of the formalism to include
  vanishing critical parameters is outlined in Appendix D") suggests
  Appendix D *enables* the limit rather than declaring it impossible. The
  alleged verbatim quote is unverified and likely paraphrased or
  hallucinated.

The arXiv:2602.12114 citation as "Cabo Bizet-Bartocci 2026" is **factually
wrong**: the paper is by Chan-López, Martín-Ruiz, Cabrera, and Paulin Fuentes,
titled *Matrix bordering structure of the Faddeev-Jackiw algorithm: Schur
complement regularization and symbolic automation*. The author attribution
is fabricated; this needs correcting in any downstream artefact that cites it.

The remainder of this document goes paper-by-paper.

---

## 1. Krupka-Voicu 2015 (arXiv:1406.6646)

**Agent claim** (Round 2, Agent E; reused by Round 3 Agent I):
> "Krupka-Voicu Theorem 1 verified symbolically for all three dynamical EOMs:
> EL_y_a(L_VT) − ε_y_a = 0 exactly."

**What the paper actually says:**

The KV paper does not contain a labelled "Theorem 1". The central object
is **Definition 1** (Eq. 12), which defines the *canonical variational
completion* of a source form `ε` as the difference

```
τ(ε) = E(λ_ε) - ε
```

where `λ_ε` is the Vainberg-Tonti Lagrangian (Eq. 11) and `E` is the
Euler-Lagrange operator. The completion `τ(ε)` is the *correction* needed to
make the system variational — NOT a guarantee that `ε` itself is variational.

In other words: KV always *constructs* a completion, but the completion
*equals zero* (i.e. `E(λ_ε) = ε`) only when `ε` is already variational —
i.e. when the Helmholtz residue δE = 0 vanishes.

The Round 2/3 agent verifications are computing `E(λ_ε) - ε` and finding it
= 0. This is not a verification of "KV Theorem 1" — it is a re-statement
of the fact that the parent source form `ε` is variational (which Round 1
Agent B already established via δE = 0). KV's *content* is the *definition*
of the canonical completion and the procedure for computing `λ_ε`; the
agents have been verifying the trivial corollary that "if ε is variational,
its canonical completion is the variational form itself".

**More importantly**: KV 2015 makes **no statement at all** about
convergence of the VT integral. The paper assumes convergence as part
of Definition 1. Quoting WebFetch on the ar5iv-rendered KV paper:
> "(c) Convergence/regularity conditions: **None are stated.** The
> Vainberg-Tonti Lagrangian (Eq. 11) uses
> `ℒ_ε = y^σ ∫₀¹ ε_σ(...,uy_J^σ) du`, but no polynomial degree or growth
> assumptions are given for convergence."

So the agents' framing of "Path A's preflight gates 1 and 2: Helmholtz residue
+ VT convergence" is a defensible *de facto* requirement, but it is *not*
something KV proves or even discusses. Convergence is a separate concern,
analyzed in Voicu 2020 — but not by KV 2015.

**Agent claim about higher-derivative source forms:**
KV does discuss this briefly. Quoting WebFetch:
> "Generally speaking, the Vainberg-Tonti Lagrangian … of order r … are of
> order 2r. Still, under certain conditions … the Vainberg-Tonti Lagrangian
> is actually equivalent to a Lagrangian of order r."

The "certain conditions" are unspecified in the WebFetch summary. This is a
non-trivial caveat: for `b5·R̃²` PGT, the source form `ε_y_a` from Round 3
Agent I is order 6 in jets, so the VT Lagrangian is naively order 12 unless
the "certain conditions" hold to bring it down to order 6. The agents do
not check whether these conditions hold — they simply observe that the EL
match works modulo IBP, which is consistent with but does not prove the
order-r reduction.

**Round-4 Verdict:**
Agent E's headline claim that VT converges polynomially for T2/T3/T4 is
operationally correct (the integrand is degree-1 polynomial in fibres), and
the symbolic check `EL(L_VT) = ε` is genuinely informative. But the framing
as "Krupka-Voicu Theorem 1 verified" is sloppy: there is no Theorem 1, and
the verified identity is the existence-of-completion definition, not a
hypothesis-to-conclusion theorem. The paper's stated *content* — Definition
1 of canonical variational completion — gives no convergence guarantee and
no order-reduction guarantee.

**Concrete language fix for any publication built on this:** "We compute
the canonical variational completion (Krupka-Voicu 2015, Definition 1) and
verify by direct calculation that `E(λ_ε) - ε = 0`, confirming the source
form is variational and that `λ_ε` is a Lagrangian for it." Drop "Theorem 1".

---

## 2. Voicu 2020 (arXiv:2009.05459)

**Agent claim** (Round 2, Agent E):
> "VT integral diverges for 4D-GB because of negative fibre homogeneity.
> TIDAL avoids this because PS-reduced ε_σ are degree-1 polynomial."

**What the paper actually says:**

WebFetch on the ar5iv version confirms that the divergence is driven by the
homogeneity degree of the source form under metric rescaling — but the
characterization is *not* "negative fibre homogeneity". Quoting:
> "the resulting improper integral diverges. This is, e.g., the case of PDE
> systems … that are homogeneous of negative degree smaller or equal to −1"

So the threshold is "homogeneity degree ≤ −1", not "negative homogeneity".
Degrees in (−1, 0) would still be fine. This is a quantitative refinement,
not a paradigm difference, but it matters: a source form of degree 0 (e.g.
linearised Einstein tensor) is on the convergence boundary, not inside the
divergent regime. PGT b5·R̃² source forms scale as degree +1 because the
parent Lagrangian is quadratic — Agent E's framing is correct on the
*outcome* but the criterion is the +1 / −1 dichotomy, not "polynomial vs
rational".

**Critical: another failure mode the agents missed.**
WebFetch (second pass on ar5iv):
> "What it does establish: equations must be 'linear in the second order
> derivatives' to admit variational structure. For the Gauss-Bonnet case,
> the authors demonstrate that the truncated equations violate this
> linearity condition (Appendix A), proving non-variationality independent
> of Theorem 3.1."

This is a separate, structural condition: the source form must be **linear
in the highest-order derivatives**. PGT b5·R̃² *linearised* satisfies this
(linearisation is linear by construction), but as Agent E notes in their
Phase 2, beyond-quadratic PGT could violate it. The criterion is a *necessary*
condition for variationality, distinct from the homogeneity-of-fibres
condition that drives the integral divergence. The agents conflate the two
under "VT convergence", but they are independent: Voicu's discussion in
Section 4 / Appendix A (linearity in highest derivatives) and Section 3.2
(homogeneity for integral) are two separate pathologies.

**The Round 2 Agent E phase-2 argument** ("no analog of Voicu's pathology in
TIDAL") rests on:
1. Polynomial source forms (degree +1) — correct.
2. No conformal-rescaling trade-off — paper is quoted for this but it's
   actually a property of the linearised regime, not of the homotopy itself.
3. No topological obstruction — true at linearised order, but Voicu's
   *primary* argument in Appendix A is about loss of variationality from
   *truncating* a topological invariant (Lanczos-Lovelock dimensional
   continuation). TIDAL's PS reduction is also a truncation (drops higher
   orders in b5), and an analogous failure mode could in principle appear
   at order b5² or higher. The agents say "TIDAL never truncates" —
   strictly false, since PS reduction is itself a series truncation.

**Round-4 verdict:** Agent E's outcome (VT integral converges for the
linearised quadratic-Lagrangian case) is correct, but the framing
oversimplifies Voicu 2020. Two substantive issues:
- The convergence criterion is "homogeneity degree > −1", not "polynomial
  vs rational".
- Linearity-in-highest-derivatives is a *separate* condition. The agents
  treat "convergence" as a unitary property; it actually has at least two
  independent failure modes.

For a publishable paper, both criteria need to be stated explicitly and
checked separately for the actual PGT source form (not just for T2, T3, T4
toys).

---

## 3. Lyakhovich 2021 (arXiv:2102.10579)

**Agent claim** (Round 1, Agent A):
> "Lyakhovich's existence theorem is constructive (homological perturbation
> theory) but proves only existence, not finite termination."

**What the paper actually says** (from WebSearch / abstract / WebFetch attempts):
> "The Batalin-Vilkovisky form of inclusion the Stueckelberg fields is worked
> out and existence theorem for the Stueckelberg action is proven."

The paper explicitly proves an existence theorem. The procedure is iterative
and based on involutive closure. WebSearch confirms:
> "With the most general closure generators, the consistent Stueckelberg
> gauge invariant theory is iteratively constructed, without obstructions
> at any stage."

So: existence theorem ✓, iterative construction ✓, no obstruction at any
stage. **The paper does NOT prove finite termination** — the iteration is
"without obstructions at any stage" but may continue indefinitely (consistent
with general homological perturbation theory).

**Agent A's claim about constraint Poisson matrix det ∝ b5:** This is a
sympy-verified result from Agent A's specific lift, not a property of
Lyakhovich's theorem itself. The structural claim in Round 1 — that *any*
Lyakhovich-recipe lift would inherit the b5 dependence — is not a logical
consequence of Lyakhovich's existence theorem. It is a separate observation
that Agent D then promotes to a no-go via a parameter-counting argument.

**Round-4 verdict:** Agent A's interpretation is essentially correct: the
Lyakhovich recipe gives existence but not finite termination, and the
specific lift Agent A computed does have det(M) ∝ b5. The structural
"no Lyakhovich lift can bridge b5 = 0" claim is *not* in Lyakhovich's
theorem; it follows from Agent A's parameter-counting and Agent D's
extension. Agent A is appropriately careful with this distinction in their
writeup.

**Caveat the agents missed:** Lyakhovich's existence theorem assumes the
involutive closure of the Lagrangian equations is well-defined and finite-
dimensional in some appropriate sense (since it begins with "the original
Lagrangian equations are complemented by all the lower order consequences").
For higher-derivative theories with infinite-dimensional jet-space
consequences, the "all lower order consequences" closure may not be a
finite operation. This is a subtle issue that the agents don't engage with
but which could matter for the b5 → 0 limit (where the constraint structure
itself rank-jumps).

---

## 4. Abakumova-Lyakhovich 2021 (arXiv:2106.09355)

**Agent claim** (Round 1, Agent D):
> "Reducible Stückelberg generators are b5-independent by construction;
> cannot bridge b5 discontinuity."

**What the paper actually says** (WebFetch on abstract):
> "We propose a general procedure for iterative inclusion of Stueckelberg
> fields to convert the theory into gauge system being equivalent to the
> original one. … In so doing, we admit reducibility of the Stueckelberg
> gauge symmetry. In this case, no pairing exists between Stueckelberg fields
> and gauge parameters, unlike the irreducible Stueckelberg symmetry."

**Agent D's structural argument** (paraphrased): the reducibility null-vectors
`Z^a` are built from order-0 EOM consequence-generator structure, which is
b5-independent; therefore the null-vectors are b5-independent; therefore
they cannot transform a b5-dependent Poisson bracket into a b5-independent
one.

This argument is sound *if* the consequence-generator structure is genuinely
b5-independent at order 0. But Agent A's own toy involved
`L = ½(∂_tφ)² − ½m²φ² − λφh − ½Mh² + ½b5(∂_t²h)²`. The order-0 consequence
structure is `λφ − Mh = 0` (variation w.r.t. h, b5 → 0 limit), which is
indeed b5-independent. But the *full* consequence structure (i.e. consequences
of variations of the b5-dependent EOMs) does involve b5. Agent D's argument
implicitly equates "consequence-generator structure" with "b5 = 0 consequence
structure", which is a particular reading of Lyakhovich's recipe.

**Round-4 verdict:** Agent D's structural argument is *defensible* but
*not airtight*. The claim "any reducible Stückelberg lift inherits
b5-independence of the null vectors and so cannot bridge the discontinuity"
is contingent on a specific choice in the consequence-generator construction
(the b5 = 0 starting point). A paper-quality version of this no-go would
need to either:
- prove that Lyakhovich's procedure *forces* the b5 = 0 starting point, or
- enumerate the alternative starting points and show each one fails.

Agent D does neither, and this is a genuine gap in the no-go theorem as
currently stated.

The Abakumova-Lyakhovich paper itself does NOT discuss parameter-discontinuity
cases — only fixed-parameter equivalences. So Agent D is correct that the
*paper* does not provide a positive recipe for bridging b5 = 0. The negative
claim "no application of this paper's framework can bridge b5 = 0" is
plausible but not conclusively proven.

---

## 5. Hinterbichler-Saravani 2015 (arXiv:1508.02401)

**Agent claim** (Round 1, Agent C):
> "HS algebraic obstruction: HS Eq. 5.5 (conformal identity for Einstein
> tensor under g→e^{-2M²π}g) absorbs all π-dependence. No parity-odd analog
> for R̃ — R̃ transforms anomalously under conformal rescalings (topological
> term), giving □π·ε·R couplings with no clean reabsorption."

**What the paper actually says:**

WebFetch confirms HS Eq. 5.5 is the conformal-rescaling identity for the
Einstein tensor:
> "Gᵤᵥ[e²M²ᵖgᵤᵥ] = Gᵤᵥ − 2M²(∇ᵤ∇ᵥπ − gᵤᵥ□π) + 2M⁴(∇ᵤπ∇ᵥπ + ½gᵤᵥ(∂π)²)"

And the paper explicitly works only with parity-even quadratic curvature:
> "(b) Parity Considerations: The paper handles only parity-even terms.
> The authors work exclusively with R² terms, Weyl tensor squared
> (Cᵤᵥᵨₛ Cᵘᵛᵨˢ), Ricci squared combinations. No parity-odd invariants like
> the Pontryagin class or Chern-Simons-type terms appear in the analysis."

**Round-4 verdict:** Agent C's interpretation is correct. The HS construction
*genuinely* does not extend to parity-odd `R̃²` because (a) the paper itself
restricts to parity-even, and (b) the conformal identity (5.5) is for the
Einstein tensor, not for the parity-odd Pontryagin density. Agent C is
appropriately cautious about this.

**Caveat:** Agent C's claim about "□π·ε·R couplings with no clean
reabsorption" is a derived structural argument, not a quote from HS. This is
fine as Agent C's own analysis, but it needs to be flagged in any publication
as a derivation rather than a literature claim.

---

## 6. Aoki-Mukohyama 2020 (arXiv:2009.11739)

**Agent claim** (Round 2, Agent F):
> "Their dynamical spin-0⁻ scalar IS the φ auxiliary at linearised order —
> strong cross-validation."

**What the paper actually says** (read directly from local TeX):

The relevant section is §V (Inclusion of additional dof). From line 497:
> "The dynamical spin-0⁻ mode of the spin connection shows up around the
> flat background when the term `d⁴x|e|·𝒳²` is added where 𝒳 is so-called
> the Holst scalar"

And immediately after (line 502):
> "However, the naive inclusion of the 𝒳² term should break the special
> structure of the ghost-free theory. The YN ghost(s) must reappear, in
> general."

This is critical: **AM explicitly state that the naive `𝒳²` term BREAKS the
ghost-free structure of the surrounding theory**. The construction in §V
is not "add `𝒳²` and get a clean spin-0⁻ scalar"; it is "add `𝒳²` plus a
specific dRGT-mass-term tuning involving arbitrary functions `c_i(φ)` of
the parity-odd scalar φ" (Eqs. 513-516, 555-557).

The section culminates in Eq. 568 — the resulting "ghost-free" theory has
a very specific form:
```
S = (M_pl²/4) ∫ ε_abcd [ e∧e∧R^cd − α e∧T^b∧K^cd
                       − (1/M_*²) e∧S^b∧R^cd ]
    + (α_{0⁻} M_pl²)/(12 M_*²) ∫ d⁴x |e| 𝒳²
    + O(M_*^{-4})
```
This is `𝒳²` (i.e. `R̃²`) appearing as the *leading* term in an *infinite
tower* of corrections, with specific tuning.

**Agent F's claim "AM's φ is our φ":** Looking at Eq. 526:
> "− (3 M_pl²(1-α)/4) · |f^a_μ| · (∂_α φ)²/(1+φ²)"

φ here is a parity-odd scalar with a *non-canonical* kinetic term
`(∂φ)²/(1+φ²)` and a sinh-canonicalisation `φ = sinh θ` (Eqs. 530-543) to
get the canonical `(∂θ)²` kinetic. This is a *very specific* dRGT-induced
scalar with a tower of self-interactions, not a generic Lagrange-multiplier
auxiliary.

In Agent F's Bopp-Podolsky construction:
```
L_aux = -¼F² - ½m_A²A² + φ(∂·A) - (1/(2b))φ²
```
the auxiliary φ is a generic Lagrange multiplier with quadratic potential
1/(2b)·φ². This is **not** AM's φ, which has a sinh-tower of self-
interactions and lives inside a specific dRGT-tuned bigravity.

Furthermore, AM's construction requires the **infinite tower**
`R(1+R/(α m²))^{-1} R` (paper line 57: "the absence of ghost at non-linear
orders requires an infinite number of higher curvature terms"). They do
*not* solve the b5·R̃² problem with a single auxiliary — they solve it by
embedding it in a dRGT-bigravity structure with infinite curvature
corrections.

**Round-4 verdict:** Agent F's "cross-validation against AM" is **substantially
weaker than claimed**. AM's φ is not the same as Agent F's φ:
- AM's φ has non-canonical kinetic `(∂φ)²/(1+φ²)` with sinh-canonicalisation.
- AM's construction requires *infinite* curvature corrections; the bare
  `𝒳²` term BREAKS ghost-freedom.
- AM's spin-0⁻ scalar lives in a dRGT-bigravity tuned mass sector, not in a
  generic single-auxiliary lift.

The most that can honestly be said is: **at strictly linearised order**, both
constructions feature *some* parity-odd scalar mode in the spectrum.
That is consistent with both, but does not constitute "strong cross-
validation" — AM's construction is a non-linearly tuned ghost-free
extension, while Agent F's is a single auxiliary that loses ghost-freedom
beyond linearised order (which Agent F themselves note: the construction
fails at O(h²) on curved background).

The cross-validation claim should be downgraded to "qualitative consistency
at linearised order; AM's full ghost-free completion requires an infinite
tower of curvature corrections that this single-auxiliary lift does not
provide".

---

## 7. Chatzistavrakidis-Ranjbar-Zekoč 2024 (arXiv:2411.16928)

**Agent claim** (Round 3, Agent J):
> "Stückelberg construction with three auxiliaries (h, b, a) gives gauge
> invariance + rank uniformity for the Curtright tensor field."

**What the paper actually says** (WebFetch on ar5iv):

The construction is for **massive (2,1) Curtright fields with FREE kinetic
structure**. The auxiliaries are graviton `h_μν`, Kalb-Ramond `b_μν`, and
1-form `a_μ` (Section 5.1, Eq. 57). The gauge-invariant field strength is
Eq. 58:
```
F̊_{μν|ρ} = T_{μν|ρ} − 2∂_{[μ}h_{ν]ρ} − 2∂_{[μ}b_{ν]ρ} + 2∂_ρ b_{μν}
                    − 2∂_ρ ∂_{[μ}a_{ν]}
```

WebFetch on the ar5iv-rendered paper:
> "The analysis is restricted to **free, massless and massive field
> theories without parity-breaking interactions**. The paper:
>  - Works exclusively with free field Lagrangians (kinetic and mass terms only)
>  - Does not discuss parity-odd couplings or CP violation
>  - Focuses on Abelian gauge theories (U(1) symmetries)
>  - Develops conserved Noether currents for minimal coupling to background fields
>  - Excludes interaction terms beyond the coupling to classical background sources"

**Agent J's verification of "gauge invariance + rank uniformity"**: The
gauge-invariance check `δF̊ = 0` for the (2,1) Stückelberg combination
follows directly from Eq. 58 by inspection — it is an algebraic identity
in the *free* case. The rank-uniformity check (`det(H_kin) = 1 - λ_a²`,
b5-independent) is a property of Agent J's specific 4-field toy, not of
the published Curtright construction itself.

**Critical: the paper's m → 0 limit is the GOLDSTONE limit**, not the
parameter-discontinuity limit:
> "(c) Massless Limit Discussion: The paper frames the discussion in terms
> of Goldstone boson interpretation rather than a classical m→0 limit.
> Section 5.2 interprets these fields as 'Nambu-Goldstone bosons for
> spontaneously broken tensor global symmetry.'"

So Agent J's translation `b5 → 0 ⇔ m_q → ∞` is the *opposite* limit from
what the paper handles. Agent J acknowledges this directly:
> "PGT mapping reversal: the constraint promotion barrier is at b5 → 0,
> which corresponds to m_q² = M/b5 → ∞. This is the 'infinite-mass
> decoupling' limit, not the massless Goldstone limit. The Stückelberg
> field strength F̊ makes BOTH limits smooth in the auxiliary phase space."

The claim "F̊ makes both limits smooth" is Agent J's own derived assertion,
not in the paper. The paper proves only the m → 0 case. Agent J's sympy
check that `det(H_kin) = 1 - λ_a²` is b5-independent is genuine, but it is
a 4-field 1+1D toy verification, not a proof for the full 4D Curtright
field.

**Caveats Agent J already notes:**
> "C1: Parity-odd ε·DT·DT cross-terms (38 terms in q-projection) NOT
> covered by published paper. A parity-odd extension would be open
> research."
> "C3: Validity at linear-flat order; nonlinear extension needs non-abelian
> Stückelberg work."

These are correct caveats. But the framing in the Round 3 synthesis —
"Path B-tensor-q ✅ Sympy-verified gauge invariance + rank uniformity" — is
overstated. What is verified is gauge invariance of `δF̊` in 1+1D plus
rank uniformity of a 4×4 Hessian in 1+1D. The paper itself does not handle:
- Parity-odd contractions (the main TIDAL use case for `b5·R̃²`).
- The b5 → 0 limit (only m → 0 Goldstone).
- Interactions / non-abelian Stückelberg.
- 4D verification of the full 36×36 Hessian.

**Round-4 verdict:** The paper is correctly identified, but Agent J's
"applies with caveats" verdict is more honestly "applies in 1+1D toy
under conditions C1-C5 that are NOT met by the actual TIDAL use case".
For TIDAL `b5·R̃²` — which is *parity-odd* by construction (the `R̃` is
the Holst dual `(1/2)ε^{abcd}R_{abcd}`) — the paper's published construction
does not apply. The "first published Stückelberg recipe applicable to PGT
tensor-torsion" framing in Round 3 synthesis is wrong: it's the first
*candidate* recipe; whether it applies is conditional on the parity-odd
extension being publishable, which Agent J flags as open research.

---

## 8. Blagojević-Cvetković 2018 (arXiv:1804.05556)

**Agent claim** (Round 1):
> "Appendix D contains the verbatim quote 'the diagonal matrix D in (D.2)
> has no valid limit for b̄ → 0'."

**What I found:**

I could not access Appendix D in full via WebFetch — both the abstract page
and the ar5iv rendering truncated before the appendix. However, the
ar5iv page does indicate Appendix D's *purpose*:
> "Extension of the formalism to include vanishing critical parameters is
> outlined in Appendix D"

This phrasing — "extension to include vanishing critical parameters" —
suggests Appendix D *enables* the vanishing-critical-parameter analysis,
not that it identifies a no-go. This is the *opposite* of Agent's claimed
verbatim quote.

Furthermore, the agent claim references "the diagonal matrix D in (D.2)",
which is awkward terminology: Appendix *D*, equation D.2, contains a *matrix
D*. This could be a real coincidence, but it could also indicate the agent
was confabulating from the appendix label.

**Round-4 verdict:** **The verbatim quote is unverified and likely
hallucinated**. I cannot definitively prove the quote does not appear in
Appendix D — I lack full access — but:
1. The paper's stated purpose for Appendix D is to *enable* analysis of
   vanishing critical parameters, not to declare them inadmissible.
2. The "diagonal matrix D in (D.2) … no valid limit" phrasing has the
   structural texture of a confabulation (labels reused from the appendix
   designation, opaque quantitative claim).

**Recommendation:** any publication citing this quote MUST be re-verified
by reading the actual paper PDF. As currently written, the quote is
unsourced.

The agents' broader claim — that BC's analysis exhibits a constraint
discontinuity at b̄ → 0 — is plausible from the paper structure but the
specific verbatim attribution is suspect and should not be reproduced
without independent verification.

---

## 9. Barker et al. 2024 (arXiv:2406.12826v3)

**Agent claim** (Round 1, Agent C):
> "Trace torsion S_μ is a Stückelberg/Goldstone field for broken Weyl
> invariance."

**What the paper actually says** (read directly from local TeX):

The paper's main thesis is that **vector torsion `T_μ/3` is identified with
the Yang-Mills-type gauge field `V_μ` of extended Weyl gauge theory (eWGT)**
when expressed in scale-invariant variables. From line 168:
> "We will now connect these three observations by showing that eWGT is
> the unique scale-invariant embedding of PGT. This will identify `T_μ/3`
> with the vector `V_μ` when expressed in scale-invariant variables, and
> thereby reveal `∂_{[μ}T_{ν]}∂^{[μ}T^{ν]}` to be a Yang-Mills-type term."

So Barker identifies the trace torsion with the Weyl *gauge* field — a
genuine propagating vector with a Maxwell-type kinetic — *not* a Goldstone
of broken Weyl invariance. The compensator field (the would-be Goldstone)
is the *scalar* `φ`, not the vector `T_μ`. Quoting line 297:
> "The compensator is purely gauge, so that the embedding theory is
> completely indistinguishable from PGT after gauge-fixing."

The compensator is the Goldstone; trace torsion is the gauge field that
*it* is the Goldstone *of*. These are dual roles.

**Critical: parity-odd terms are explicitly excluded.** From line 109,
footnote:
> "In this letter we omit the parity-odd invariants only out of simplicity;
> there are no convincing theoretical grounds for excluding them"

And line 300:
> "The Yang-Mills-type actions … are restricted to parity-even terms for
> simplicity: the parity-odd extensions should be considered. The mass
> spectrum of parity-violating PGT was found in [Karananas:2014pxa], and
> confirmed in [Blagojevic:2018dpz], however the massless spectra and
> unitarity conditions of the various critical cases of the theory have
> not been thoroughly explored, nor have more than a handful of such
> cases been identified to date."

So Barker explicitly leaves the **parity-odd extension** as future work.
And `b5·R̃²` is parity-odd by construction.

**Round-4 verdict:** Agent C's claim is **substantially mistaken**. Barker
2024 does NOT support the proposition that "trace torsion is a Stückelberg/
Goldstone for broken Weyl invariance". The paper says:
- Trace torsion is the Weyl *gauge* field of eWGT (a propagating vector),
  not a Goldstone.
- The Goldstone (compensator) is the *scalar* φ.
- The construction is for parity-EVEN PGT only; parity-odd extensions
  (i.e. the actual TIDAL `b5·R̃²` case) are explicitly future work.

Path B-trace as "✅ Established literature (Barker et al. 2024)" in the
Round 2 / Round 3 synthesis is **wrong**. The construction is not
established for parity-odd, and the cited paper does not even attempt the
parity-odd case. The trace-sector handling for `b5·R̃²` should be re-marked
as "open / future work; awaits parity-odd extension of Barker 2024".

---

## 10. arXiv:2602.12114 (alleged "Cabo Bizet-Bartocci 2026")

**Agent claim** (Round 1):
> "Modern FJ Schur-complement view confirms FJ termination requires
> non-degeneracy."

**What the paper actually is** (WebFetch verified):

- **Authors**: E. Chan-López, A. Martín-Ruiz, Jaime Manuel Cabrera, and
  Jorge Mauricio Paulin Fuentes.
- **Title**: "Matrix bordering structure of the Faddeev-Jackiw algorithm:
  Schur complement regularization and symbolic automation".
- **Submitted**: 12 Feb 2026 to arXiv.
- **Peer review status**: Preprint only, not peer-reviewed (as of WebFetch
  retrieval).

**There is no author named Cabo Bizet or Bartocci on this paper.** The
attribution as "Cabo Bizet-Bartocci 2026" is **fabricated**.

I found no other arXiv:26XX preprints by Cabo Bizet & Bartocci on this
topic, and I have low confidence such a paper exists. (Cabo Bizet and
Bartocci are real physicists, but the citation as given does not match a
real paper that I can verify.)

**Round-4 verdict:** Citation is wrong. Whoever produced this citation
either confused authors or hallucinated the attribution. Any downstream
publication citing this paper MUST update the author list to:
"Chan-López, Martín-Ruiz, Cabrera, Paulin Fuentes (2026), arXiv:2602.12114".

The technical content of the paper (FJ Schur complement) is real and
relevant, but it is a 2026 preprint, not peer-reviewed, and so any claim
"Modern FJ Schur-complement view confirms …" should be hedged: it is a
recent preprint making a preliminary case, not an established result.

---

## Cross-cutting issues identified

### A. Conflation of three distinct conditions for "Path A works"

The agents repeatedly bundle three logically independent properties under
"Path A preflights pass":

1. **Helmholtz residue δE = 0** (Agent B): tests whether the source form
   is variational at all. This is the *necessary* condition for *any*
   Lagrangian to exist for the source form.
2. **VT integral convergence** (Agent E): tests whether the *specific* KV
   homotopy formula (Eq. 11) gives a finite result. Voicu 2020 identifies
   this as a separate gate.
3. **Linearity in highest derivatives** (Voicu 2020 Appendix A): a *third*
   condition for the result to be a well-defined Lagrangian. Agents do not
   discuss this separately.

For the linearised PGT pipeline, all three happen to hold, but conflating
them obscures the genuine physics: it is *because* the Lagrangian is quadratic
that the source forms are linear, polynomial-in-fibres, and degree-+1.
Beyond linearisation, all three could fail independently, and the "full
PGT b5·R̃² is unobstructed" claim relies on all three holding *and* on the
source-form structure being preserved under going beyond toy models. None
of the three has been verified for the actual 38-component PGT theory.

### B. "Cross-validation" inflation

Several agent claims of "cross-validation" between paths are genuine
qualitative consistency checks but not rigorous equivalences:

- Agent F vs AM 2020: linearised parity-odd scalar in spectrum ≈ AM's
  spin-0⁻ in a dRGT-tuned bigravity. Qualitatively consistent; not
  equivalent constructions.
- Agent J vs Agent F: scalar-auxiliary lift (axial) vs three-field
  Stückelberg (tensor-q). Mechanism *similar* (heavy-mode decoupling) but
  different constructions on different sectors.
- Agent I joint VT consistency check on 6-field parent: confirms the
  parent L is variational (which we already know — it was derived from a
  Lagrangian). This is a sanity check, not new content.

The cumulative effect is to make Path A + Path B sound like 2-3 independent
verifications when they are actually one verification + 1-2 consistency
checks.

### C. The metric h₄, h₇, h₉ blocker is repeatedly framed as "the genuine
remaining problem". This framing is correct, but the agents have **not
explicitly verified** that the BC h₄, h₇, h₉ subspace exhibits the same
PU structure they assume. Agent I's T4 toy is *modelled on* the BC
structure but is not the BC structure — it has generic K_{ab} cross-
couplings rather than the actual mass matrix from BC Appendix D (which
the agent explicitly couldn't access). The "phase-space dim 6 → 30" jump
is for the toy, not for actual PGT.

### D. None of the agents independently verified the Recipe 1 preflight
claim — that `b5·R̃² → q-irreducible` produces *standard-kinetic* `(∂q)²`
not Pais-Uhlenbeck `(∂²q)²`. Agent H (Round 3) ran sympy on a q-projection
toy, but it's a sympy verification of a *projection schema*, not of the
actual `b5·R̃²` decomposition in xAct. The Round 3 verdict "PASS" is
therefore conditional on a sympy implementation of the projection
matching the actual PGT pipeline — which has not been independently
verified.

---

## Summary of overstated / wrong / unverified claims

| Claim | Status | Severity |
|-------|--------|----------|
| "Krupka-Voicu Theorem 1 verified symbolically" | Misnomer — paper has Definition 1, not Theorem 1 | LOW (cosmetic) |
| "VT diverges for negative fibre homogeneity" | Threshold is degree ≤ −1, not "negative" generally | LOW |
| "Linearity in highest derivatives" gate not separately verified | Voicu 2020 condition skipped | MEDIUM |
| Lyakhovich 2021 existence theorem applies | True | OK |
| Abakumova-Lyakhovich 2021 b5-independence | Argument depends on starting-point choice; not airtight | MEDIUM |
| Hinterbichler-Saravani parity-even only | Correct | OK |
| AM 2020 cross-validation with axial BP | Substantially weaker than claimed (different φ, infinite tower vs single aux) | MEDIUM |
| Curtright Stückelberg "applies to PGT tensor-torsion" | Conditional on parity-odd extension that is open research | HIGH |
| BC 2018 "Appendix D verbatim quote" | UNVERIFIED, possibly hallucinated | HIGH |
| Barker 2024 "trace torsion is Goldstone" | WRONG — Barker says trace torsion is gauge field, parity-odd excluded | HIGH |
| arXiv:2602.12114 "Cabo Bizet-Bartocci 2026" | Author attribution FABRICATED | HIGH |

---

## Headline-claim revision

Round 3's executive summary states the constraint-promotion barrier is
"substantially resolved" with two complementary publishable paths. After
this audit:

**Path A (Vainberg-Tonti):**
- VT integral converges for linearised quadratic Lagrangians ✓
- Helmholtz residue vanishes ✓
- "Krupka-Voicu Theorem 1" should be re-titled "Krupka-Voicu Definition 1"
- Linearity-in-highest-derivatives needs separate explicit verification
- Caveat: only verified at toy level (T2, T3, T4); 38-component PGT not
  done.
- The metric h₄, h₇, h₉ subspace inherits PU at the Lagrangian level
  (Round 3 Agent I).

**Path B-axial (Bopp-Podolsky):**
- Construction is self-consistent at linear-flat order ✓
- "AM cross-validation" should be downgraded to "qualitative consistency,
  not equivalence" (different φ, AM requires infinite tower for ghost-
  freedom).
- Curved-background extension fails at O(h²) (Agent F flagged this).

**Path B-trace (Barker conformal embedding):**
- Status: **NOT established**. Barker 2024 explicitly excludes parity-odd
  PGT, which is the actual TIDAL use case. Path B-trace for `b5·R̃²` is
  open / future work, not "established literature".

**Path B-tensor-q (Curtright Stückelberg):**
- Status: **conditional on a parity-odd extension that does NOT exist in
  the published literature**. The Chatzistavrakidis-Ranjbar-Zekoč paper
  handles only parity-even free fields. Agent J's translation to PGT
  parity-odd `b5·R̃²` is a hypothesis, not a verified construction.

**Three convergent no-go theorems for metric Pais-Uhlenbeck subspace:**
- Round 1 Agent A and D are the strongest (sympy-verified det = b5^N).
- Round 2 Agent G's "dual no-go" is convincing but specific to 2-form
  auxiliaries; broader auxiliary classes are not excluded.
- The BC 2018 Appendix D quote is unverified.

**Net effect on publishability of the three proposed publications:**

1. **Publication A (Path A applied to PGT b5·R̃²)** — still publishable in
   principle, but should: (a) drop "Theorem 1" language, (b) explicitly
   verify all three Voicu gates (homogeneity, linearity in highest
   derivatives, integral convergence), (c) acknowledge that the verification
   is at toy level and full-PGT verification is still future work.

2. **Publication B (Sectoral Stückelberg recipe)** — significantly weaker
   than presented. Two of the four sectors (trace, tensor-q) are
   *conditional on parity-odd extensions of published works*, which makes
   the "complete sectoral recipe modulo metric" framing premature. Honest
   framing: only the axial sector has a verified construction; trace and
   tensor-q await published parity-odd extensions.

3. **Publication C (Three convergent no-go theorems)** — strongest of the
   three, but the BC 2018 Appendix D verbatim quote MUST be verified by
   reading the actual paper before any citation. If the quote is not in
   the paper, the historical-anchor argument needs to be re-grounded on
   what BC actually wrote.

---

## Recommendations

1. **Verify the BC 2018 Appendix D quote**: download the actual PDF (not
   ar5iv) and grep for "diagonal matrix D" and "b̄ → 0". If absent, remove
   the verbatim attribution from all downstream artefacts.

2. **Fix arXiv:2602.12114 attribution**: the authors are Chan-López,
   Martín-Ruiz, Cabrera, Paulin Fuentes — not Cabo Bizet-Bartocci. Update
   any TeX or Markdown that cites this paper.

3. **Downgrade Path B-trace from "established" to "open / parity-odd
   extension required"**: Barker 2024 does not handle parity-odd, and the
   "trace torsion is Goldstone" framing is wrong.

4. **Downgrade Path B-tensor-q from "applies with caveats" to "conditional
   on parity-odd extension"**: Chatzistavrakidis-Ranjbar-Zekoč 2024 handles
   only parity-even free fields; the parity-odd extension is open research.

5. **Re-title "Krupka-Voicu Theorem 1" verifications**: there is no
   Theorem 1; use "Krupka-Voicu Definition 1" or "Helmholtz consistency
   check" depending on what is actually being verified.

6. **Add Voicu 2020 linearity-in-highest-derivatives gate** as a separate
   preflight, distinct from VT integral convergence.

7. **Revise AM cross-validation language**: it's qualitative consistency at
   linearised order, not full cross-validation. AM's φ has non-canonical
   kinetic and lives in a dRGT-bigravity tower, not a single-auxiliary
   lift.

8. **Independently verify Recipe 1**: Agent H's sympy-projection result
   should be cross-checked against the actual xAct `b5·R̃²` decomposition,
   not just a hand-coded sympy schema.

The investigation has identified a real and interesting mathematical
problem; the technical sympy work is genuine; but the literature attributions
and the "publishable now" framing of Round 3 need substantial tightening
before any artefact goes out. The serious literature errors (Barker, BC,
the 2602.12114 attribution) would not survive peer review and need to be
fixed before submission.
