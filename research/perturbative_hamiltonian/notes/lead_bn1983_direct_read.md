# Lead: Direct Read of Blagojević-Nikolić 1983 Papers

**Date:** 2026-04-27
**Investigator:** Claude (Opus 4.7, harness session, Phase 2.4)
**Source files (local cache):**
- `/workspaces/torsion-gertsenshtein/literature/BF02721794/BF02721794.pdf` → `/tmp/BN1983_NuovoCim.txt` (1139 lines)
- `/workspaces/torsion-gertsenshtein/literature/PhysRevD.28.2455/PhysRevD.28.2455.pdf` → `/tmp/BN1983_PRD.txt` (1916 lines)

**BibTeX-style citations:**

```
@article{BlagojevicNikolic1983NCB,
  author  = {Blagojević, M. and Nikolić, I. A.},
  title   = {Hamiltonian Structure of the Theory of Gravity with $R + T^2$ Type of Lagrangian},
  journal = {Nuovo Cimento B},
  volume  = {73},
  number  = {2},
  pages   = {258--273},
  year    = {1983},
  note    = {Received 10 May 1982; revised 8 September 1982; published 11 February 1983},
  doi     = {10.1007/BF02721794}
}

@article{BlagojevicNikolic1983PRD,
  author  = {Blagojević, M. and Nikolić, I. A.},
  title   = {Hamiltonian Dynamics of {Poincar\'{e}} Gauge Theory: General Structure in the Time Gauge},
  journal = {Phys. Rev. D},
  volume  = {28},
  number  = {10},
  pages   = {2455--2463},
  year    = {1983},
  note    = {Received 22 February 1983; published 15 November 1983},
  doi     = {10.1103/PhysRevD.28.2455}
}
```

---

## Section 1 — BN 1983 Nuovo Cimento (R+T² PGT)

### 1.1 Lagrangian form

The paper restricts to `R + T²` PGT (linear in the Ricci scalar, quadratic in the torsion, **NO curvature-squared terms**). Verbatim from the abstract (`/tmp/BN1983_NuovoCim.txt:16-21`):

> "The IIamiltonian structure of a class of gravity Lagrangians which are linear in the scalar curvature and quadratic in the torsion is studied by using Dirac's general method for systems with constraints. It is found that the Lorentz gauge potentials can be eliminated from the theory by using constraint equations, without affecting the Dirae brackets of the remaining variables."

The Lagrangian split is given at line 147, eq. (6):

> "ℒ_G = a R + ℒ_T ≡ ℒ_ET" (`/tmp/BN1983_NuovoCim.txt:147`)

with the torsion-squared part decomposed irreducibly (`/tmp/BN1983_NuovoCim.txt:170-189`):

> "ℒ_T = α t_{ijk} t^{ijk} + β v_i v^i + γ a_i a^i" (eq. (7))
> "A = α/2 - γ/18, B = α/2 + γ/9, C = -α/2 + β" (eq. (9))

**No R² or R̃² term appears anywhere.** This is explicit in the paper's framing (`/tmp/BN1983_NuovoCim.txt:200-206`):

> "The presence of terms which are quadratic in the components of torsion does not seem to be able to alter the dynamics of the EC theory in an essential way: we again expect the existence of constraints, which will enable us to eliminate A^{ij}_μ as physical degrees of freedom from the theory. … Some of the results obtained here may be useful in an analysis of the general theory based on the Lagrangian density (5)."

### 1.2 Methodology

Standard Dirac-Bergmann constrained-Hamiltonian analysis in the time gauge `h_a^0 = 0`. Key steps:

- Define conjugate momenta from `ℒ_ET`, identify primary constraints (`/tmp/BN1983_NuovoCim.txt:271-289`).
- Identify "if-constraints" — relations that *become* primary constraints only when specific parameter combinations vanish (`/tmp/BN1983_NuovoCim.txt:289-336`).
- Build the canonical Hamiltonian in the "Dirac form" `H = b_0^a ℋ_a + (1/2) A^{ab}_0 ℳ_{ab} + …` via the appendix-A trick: identify a "field-strength" combination `F^i = q̇^i + A^i` so that the Lagrangian is quadratic in `F^i` only (`/tmp/BN1983_NuovoCim.txt:976-1023`).
- Work out consistency conditions, classify constraints into first/second class, construct Dirac brackets.

### 1.3 The "if-constraints" terminology

The "if-constraints" terminology that Round 1 attributed to BN is **introduced in this Nuovo Cim paper**, verbatim at `/tmp/BN1983_NuovoCim.txt:341-344`:

> "We will call these expressions if-constraints. It is interesting to note that they will not change in the general theory (5), since R² terms do not contain velocities ḣ_a^μ."

The structure is: a relation of the form `α_i F^i + β_i ≈ 0` after diagonalisation is genuinely a constraint **iff** the corresponding eigenvalue `α_i = 0`. The factor `1 - λ(α_i)` in the Hamiltonian (where `λ(x) = 1/x` for `x ≠ 0`, `0` for `x = 0`, eq. 22 at `/tmp/BN1983_NuovoCim.txt:402-407`) implements case-by-case dispatch.

### 1.4 Treatment of parameter limits / critical surfaces

The paper distinguishes massive-tordion (parameters `a-3α/2, a+3β/2, a-2γ/3` non-zero) and massless-tordion (vanishing) cases. Verbatim at `/tmp/BN1983_NuovoCim.txt:691-699`:

> "We note that, if the constants a - 3α/2, a + 3β/2 and/or a - 2γ/3 vanish, the theory becomes inconsistent unless σ^T = 0, σ^V = 0 and/or σ^A = 0, which are definite constraints on the matter spin tensor. The vanishing of these constants is related to the question of the existence of the massless tordions in the weak-field approximation of the theory… Staying out of this problem which has not been discussed carefully in the literature and leaving it to a separate study in the future, we will limit ourselves to the case of the massive tordions: all the constants a - 3α/2, a + 3β/2 and a - 2γ/3 are taken to be different from zero."

**Critical observation:** the massless-tordion case (which is where parameter limits like `b5 → 0` would induce a constraint-promotion phenomenon) is **explicitly deferred** in this paper. The Nuovo Cim paper handles only the "all eigenvalues non-zero" sector; vanishing-eigenvalue cases are left for "a separate study in the future".

### 1.5 Constructive method?

The constructive recipe in the Nuovo Cim paper is the `λ(α_i)` factor in the Hamiltonian (`/tmp/BN1983_NuovoCim.txt:402-410`, eqs. 21-22):

> "ℋ_T = b_0^a (…) + (1 - λ(β)) b_0^a φ^b_a/3 + (1 - λ(α - 4γ/9)) b_0^a φ^[ab] + (1 - λ(α)) b_0^a φ^{(ab)} - … T_{ab0}(0)"

This is a *case-discriminator*: when an `α_i ≠ 0` you have the kinetic-form Hamiltonian; when `α_i = 0` the corresponding if-constraint becomes a real constraint and the term `(1/α_i) … = λ(α_i) …` is set to zero by hand. The recipe applies to **2nd-order Lagrangians whose velocities only enter via a single field-strength quadratic form** (Appendix A construction `L = α_{ij} F^i F^j + β_i F^i + L(F^i = 0)`, `/tmp/BN1983_NuovoCim.txt:982-988`).

### 1.6 Acknowledgement of the R²-class limitation

Multiple verbatim acknowledgements that this paper does NOT cover R²-class theories:

- `/tmp/BN1983_NuovoCim.txt:60-61`: "a step closer to a complete unterstanding of the general theory, which contains also terms quadratic in the curvature."
- `/tmp/BN1983_NuovoCim.txt:733-736`: "if-constraints (17) are, in fact, always present in the theory. This will not necessarily be true in the general case (6) [meaning the full R+T²+R²], due to the possible existence of the components of A^{ij}_μ which are dynamically independent of the other variables, i.e. which can propagate in the weak-field approximation of the theory."
- `/tmp/BN1983_NuovoCim.txt:967-969`: "By clarifying the structure of the theory of the R+T² type we are making a step toward a complete understanding of the structure of the general theory based on the Lagrangian density (5). Having that in mind, we used a formalism which will be easy to generalize to that case."
- `/tmp/BN1983_NuovoCim.txt:1115-1116`: "The linearity of this theory in the variable A^{ab}_0 does not persist in the case of more general Lagrangian densities (5) and (6)."

The Nuovo Cim paper is **explicitly a stepping stone to the more general theory**, not the general analysis itself.

---

## Section 2 — BN 1983 Phys. Rev. D (R+T²+R² PGT)

### 2.1 Lagrangian form — CRITICAL FINDING

**The PRD paper covers the FULL R+T²+R² Lagrangian, including curvature-squared terms.** Verbatim from the abstract (`/tmp/BN1983_PRD.txt:28-47`):

> "general aspects of the Hamiltonian structure of Poincare gauge theory of gravity with R + T² + R² type of Lagrangian are investigated in the time gauge. The explicit form of the Hamiltonian is found, taking care of the fact that some of the primary constraints exist only if some specific relations among the constants of the theory are satisfied."

And from the introduction (`/tmp/BN1983_PRD.txt:74-79`):

> "The simplest generalization of this theory is based on the Lagrangian density of R + T² type … The addition of a T² term does not change the essential features of the Hamiltonian dynamics: the Lorentz gauge potentials are, in both cases, not independent dynamical variables. The situation will change in the general R + T² + R² case, which we are going to discuss here."

The full Lagrangian appears in `/tmp/BN1983_PRD.txt:182-227`, eqs. (2.4)-(2.6):

> "ℒ_G = ℒ_E + ℒ_T + ℒ_R … ℒ_R = b₁ R_{ijkl} R^{ijkl} + b₂ R_{ijkl} R^{kl ij} + b₃ R_{ij} R^{ij} + b₄ R_{ij} R^{ji} + b₅ R² + b₆ (ε^{ijkl} R_{ijkl})²" (eq. 2.6)

**This is exactly TIDAL's `b1…b6` parameter family** (modulo Bach-Lanczos, "only five constants are independent" `/tmp/BN1983_PRD.txt:228-229`). The `b5 R²` term and the `b6 (εR)²` parity-odd squared-curvature term (the closest analogue of TIDAL's `b5·R̃²`) are **explicitly in scope**.

This **inverts** the audit's revised framing. The audit corrected an earlier Nuovo Cim B 84:25 (1984) mis-citation by re-citing Nuovo Cim B 73:258 (1983) as "R+T² only", then concluded that BN 1983 work is R+T² only and hence "background diagnostics not constructive bridge". **That conclusion is wrong: the PRD paper of the same year handles the full R+T²+R² case, including R² and `(εR)²`.**

### 2.2 Methodology

Same Dirac-Bergmann-time-gauge framework as the Nuovo Cim paper, but extended to handle the new velocity dependence introduced by the curvature-squared term. The extra step is identified verbatim at `/tmp/BN1983_PRD.txt:430-436`:

> "The construction of the canonical Hamiltonian density seems rather complicated due to the fact that some of the expressions (2.12)–(2.15) (if-constraints) may become primary constraints if the parameters of the theory satisfy certain relations. We will take care of this possibility by using a method outlined in Appendix C. To do this we first note that the velocity variables ḣ_{a,0} and Ȧ^{ij}_0 are contained only in the field strength components T_{ab0} and R_{ij,0}, respectively."

Crucially, the velocity Ȧ^{ij}_0 (the time derivative of the *connection*) appears only in `R_{ij,0}` — i.e. the curvature time-component carries the velocity. This is exactly the Ostrogradsky-doubling structure for an `R²` term: `R²` contains `(R_{ij,0})²` which contains `(Ȧ^{ij}_0)²`. The PRD paper therefore *does* live in the same higher-derivative-order regime as TIDAL's `b5·R̃²`, viewed through this lens.

### 2.3 If-constraints with critical parameters — terminology

The PRD paper introduces the **second** set of if-constraints, this time arising from the `R²`-induced velocity dependence in `Ȧ^{ij}_0`, verbatim at `/tmp/BN1983_PRD.txt:494-509`:

> "Equations (2.14) and (2.15) give rise to primary constraints if some of the constants, multiplying velocities in R_{ij,0}, vanish."

The full table of if-constraints and their critical-parameter-vanishing conditions is reproduced in `/tmp/BN1983_PRD.txt:1133-1217`, Table I. Each row maps a parameter combination to: which if-constraint becomes primary, and which irreducible component of A^{ij}_μ becomes "unphysical" (i.e., is eliminated by the constraint rather than propagating).

| `/tmp/BN1983_PRD.txt:1138` | β = 0 | gives constraint φ^a_a |
| `/tmp/BN1983_PRD.txt:1139` | a₅ + 12 a₆ = 0 | gives constraint φ_b^{a0} |
| `/tmp/BN1983_PRD.txt:1141` | α₄ + 2 α₃ = 0 | gives constraint φ_{[ab]} |
| `/tmp/BN1983_PRD.txt:1144` | α = 0 | gives φ_{(ab)} |
| `/tmp/BN1983_PRD.txt:1166` | 3 a₂ + 2 a₅ = 0 | gives φ_{(a0b)} |
| `/tmp/BN1983_PRD.txt:1170` | α + β = 0 | gives φ_{a0} |
| `/tmp/BN1983_PRD.txt:1174` | a₄ + a₅ = 0 | gives related X_{ab} |
| `/tmp/BN1983_PRD.txt:1175` | a₁ + a₃ = 0 | gives φ_A |
| `/tmp/BN1983_PRD.txt:1177` | 3 a₂ + 4 a₃ = 0 | gives φ_{ab} |

The vanishing of `a₅` (the coefficient of the `R²` invariant in the Bach-Lanczos basis used here) appears in the row `a₅ + 12 a₆ = 0`. The parity-odd analogue uses `a₆` (coefficient of the `(εR)²` invariant). **TIDAL's `b5 R̃²` corresponds, modulo basis change, to a specific direction in BN's `(a₅, a₆)` plane.**

### 2.4 Treatment of parameter limits / critical surfaces

The PRD paper *does* handle vanishing-critical-parameter cases — but with a key restriction. From `/tmp/BN1983_PRD.txt:1052-1062`:

> "The analysis of the consistency conditions of the if-constraints (2.12)–(2.15), in the case corresponding to the massive tordions, shows that whenever one of the if-constraints becomes a constraint, then it does not commute either with the constraint belonging to the same pair (if it exists) or with its own consistency condition. As a consequence, if-constraints which become constraints are of a second class. Therefore, the corresponding multipliers in A_TOT (irreducible components of u^{ab}, u^{abc}, and/or u^{abA}) are determined from further consistency requirements. **In the case corresponding to the massless tordions (with a new kind of gauge symmetry in the theory) more detailed considerations are necessary for a complete understanding of the situation.**"

So the PRD paper:
- **Solves** the case where one if-constraint promotion induces *massive* tordions (second-class promotion, multipliers determined). This is the "well-behaved" rank-jump.
- **Defers** the massless-tordion case (where the rank-jump opens a new gauge symmetry, requiring "more detailed considerations").

Massive tordion treatment is closely paralleled at `/tmp/BN1983_PRD.txt:973-1008` for the (M_0a^c, π_b^{ab}) pair: when both critical parameters vanish ("β = 0 and a₅ + 12a₆ = 0 vanish, then both φ^a_a and φ_b^{a0} become primary constraints"), the constraints can be used to express M_0a^c and π_{ab}^c in terms of other dynamical variables. When only one vanishes, one variable is eliminated; the other becomes a propagating tordion. When both are non-zero, "the variables M_0a^c and π_{ab}^c are independent Hamiltonian degrees of freedom" describing a propagating tordion of mass dependent on `a + 3β/2`.

### 2.5 Does PRD cover R²-class? — **Yes, generally; subject to the massive-tordion restriction.**

The PRD paper's framework treats vanishing critical parameters in the R+T²+R² Lagrangian *constructively*, but only in the massive-tordion case. The massless-tordion case (where TIDAL's deepest rank-jump pathologies presumably live) is acknowledged as needing a different analysis.

This is **importantly different** from "R+T² only ↔ background diagnostic". The PRD paper *is* a constructive Hamiltonian analysis for a Lagrangian class that includes (a basis-rotated version of) `b5·R̃²`. Its limitations are:

1. **Restricted to massive tordions.** The vanishing-critical-parameter case where the would-be tordion becomes massless (and thereby a new gauge symmetry emerges) is deferred.
2. **Time gauge fixed.** All analysis is in `h_a^0 = 0`. The companion paper "without gauge fixing" was promised in `/tmp/BN1983_PRD.txt:96-97`: "The form of the theory without gauge fixing will be studied in our next paper."
3. **The general degeneracy structure (multiple critical parameters vanishing simultaneously) is treated case-by-case via Table I**, not via a unified algebraic recipe like Blagojević-Cvetković 2018's `det F = 0` analysis.

### 2.6 Summary — does this transfer to TIDAL?

Yes, **partially**. The PRD paper provides:
- A constructive Hamiltonian for the same Lagrangian class TIDAL studies (R+T²+R²), in the time gauge, in the massive-tordion regime.
- The "if-constraints with critical parameters" terminology that names the rank-jump phenomenon TIDAL encounters at b5 = 0.
- A worked-out Dirac-bracket prescription when critical parameters vanish, *if the resulting promoted constraint is second-class* (massive tordion).

What it does NOT provide:
- A treatment of the *massless-tordion* case at the rank-jump point, which is exactly where TIDAL's `b5 → 0` transition sits (h₄, h₇, h₉ become Lagrange multipliers, not propagating fields with mass).
- A unified algebraic recipe for arbitrary critical-parameter combinations (the table-of-cases approach scales poorly).
- Any treatment outside the time gauge.

---

## Section 3 — 25-year-history claim verification

The audit's "25-year history of named-but-unsolved" framing places the BN 1983 papers as the founding statements of the constraint-promotion phenomenon. Verifying this against the actual papers:

### 3.1 Do BN 1983 explicitly discuss "constraint promotion"?

**Yes, by the exact name "if-constraints" (see Section 1.3, `/tmp/BN1983_NuovoCim.txt:341-344`).** The phenomenon is defined: relations that become primary constraints when specific parameter combinations vanish. This is the Round 1 sense of "constraint-promotion" — modulo terminology change ("if-constraint → primary constraint" rather than "algebraic constraint → propagating mode").

**Subtlety:** BN 1983's "promotion" goes the *opposite direction* from TIDAL's `b5 → 0` transition:
- BN 1983: when critical parameter `c → 0`, an if-relation **becomes** a primary constraint (a new constraint appears).
- TIDAL `b5 → 0`: at b5 = 0, h₄/h₇/h₉ ARE algebraic (Lagrange multiplier) constraints; at b5 ≠ 0, those constraints **break** and the fields become propagating modes (a constraint disappears).

Meta-N's BC Appendix D analysis (`research/perturbative_hamiltonian/meta_reviews/meta_review_N_BC_AppendixD_content.md:108-113`) flags the same inversion:

> "Wrong constraint topology. BC's 'if-constraint' is `φ_n = 0` becoming a primary constraint when one critical parameter `cn` vanishes. TIDAL's b5 = 0 limit involves *primary constraints disappearing as you move off b5 = 0*."

This inversion applies equally to BN 1983 (since BC 2018's "if-constraints" inherit the BN 1983 terminology). **The "25-year history" framing is therefore accurate at the level of nomenclature and constraint-degeneracy phenomenology, but the directionality is opposite — TIDAL is in the "constraint-disappearance" branch, not the "constraint-promotion" branch.**

### 3.2 Do BN 1983 cover R²-class theories?

**Nuovo Cim:** No (R+T² only, Section 1.6).
**PRD:** Yes (R+T²+R², Section 2.1).

The audit's framing in `FINAL_ASSESSMENT.md:50-55`:

> "The 25-year history of named-but-unsolved — Blagojević-Nikolić (1983-84), Yo-Nester-Ni (1999, 2002), Blagojević-Cvetković (2018) — is real as historical framing."

is **correct that BN 1983 PRD covers R²-class theories**. The audit's text in `FINAL_ASSESSMENT.md:351-355`:

> "Blagojević-Nikolić 1983 (Nuovo Cim B **73**:258, doi:10.1007/BF02721794) — earlier mis-cited as Nuovo Cim B 84:25 (1984). User-retrieved 2026-04-27 to `literature/BF02721794/`. The actual paper is on $R + T^2$~PGT (linear in $R$, quadratic in torsion), NOT $R+R^2+T^2$ as previously claimed."

is correct about the Nuovo Cim paper but **omits the PRD paper of the same authors and year (PRD 28:2455, 1983), which DOES cover R+T²+R²**. The two papers form a deliberate pair — the Nuovo Cim paper handling R+T² as a stepping stone, the PRD paper handling the full general case. Several explicit cross-references establish this:

- Nuovo Cim `/tmp/BN1983_NuovoCim.txt:967-969`: "we are making a step toward a complete understanding of the structure of the general theory based on the Lagrangian density (5)" (where eq. 5 is the full R+T²+R² form).
- PRD `/tmp/BN1983_PRD.txt:1872-1873` cites the Nuovo Cim paper as ref. 8: "M. Blagojevic and I. A. Nikolic, Nuovo Cimento 73B, 258 (1983)" — so PRD references its own Nuovo Cim companion.

### 3.3 Verdict on the "25-year unsolved" claim

The historical claim survives, but with sharper specification:

- **What is solved** (BN 1983 PRD): R+T²+R² PGT in the time gauge, in the massive-tordion regime, with a case-by-case if-constraint dispatch.
- **What is open** (acknowledged by BN 1983 PRD itself): the *massless-tordion* case at vanishing-critical-parameter loci, which is where TIDAL's `b5 → 0` rank-jump-into-Lagrange-multiplier phenomenon sits.

So the audit's "25-year history" framing is **directionally right** — there is a real 1983 starting point and the *specific* problem TIDAL faces (massless-rank-jump in higher-derivative PGT) was already flagged as needing further analysis in 1983. The chain BN 1983 PRD → Yo-Nester-Ni 1999, 2002 → BC 2018 traces 35 years of incremental partial answers, **none of which fully address the massless-rank-jump case**. The phrase "25-year unsolved" should be sharpened to "**40-year-acknowledged-as-needing-more-work**, with Yo-Nester-Ni 1999 → 2002 → BC 2018 as partial answers".

---

## Section 4 — New references identified

Bibliography of BN 1983 PRD (`/tmp/BN1983_PRD.txt:1856-1916`) yields the following references that the original investigation may not have fully exploited:

### 4.1 W. Szczyrba, Phys. Rev. D 25, 2548 (1982) — `/tmp/BN1983_PRD.txt:1874`

Cited as ref. 9: "The Hamiltonian dynamics of Poincare gauge theory is studied in Ref. 9 by using the general geometric constructions; the form of the gravitational Lagrangian density is quite general (unspecified) and matter is described by a tensor field" (`/tmp/BN1983_PRD.txt:80-83`).

This is a **pre-1983 general-Lagrangian Hamiltonian analysis of PGT**. It predates BN 1983 and uses geometric (not Dirac-method) techniques. **Worth checking for constructive recipes that BN may have rediscovered.**

### 4.2 K. Hayashi, T. Shirafuji, Prog. Theor. Phys. 64, 866 (1980); 64, 833 (1980); 64, 1435 (1980); 64, 2222 (1980); 65, 525 (1981) — `/tmp/BN1983_PRD.txt:1881-1883`

Cited as ref. 12. This is the **canonical pre-1983 series of papers on PGT Lagrangian classification** (Bach-Lanczos identity, irreducible decompositions, mass formulas for tordions). Almost certainly already in the project's reference list, but worth verifying that the *Hamiltonian* aspects are not understated. The Hayashi-Shirafuji series was the source of the parameter classification BN 1983 PRD adopted.

### 4.3 E. Sezgin, P. van Nieuwenhuizen, Phys. Rev. D 21, 3269 (1980) — `/tmp/BN1983_PRD.txt:1907-1910`

Cited as ref. 18 in the context of tordion mass formulas. Sezgin-van Nieuwenhuizen 1980 is a **canonical early reference on PGT particle spectrum** (and its critical surfaces for masslessness). Should already be in `Literature/`.

### 4.4 S. Miyamoto, T. Nakano, T. Ohtani, Y. Tamura, Prog. Theor. Phys. 66, 481 (1981) — `/tmp/BN1983_PRD.txt:1911-1913`

Cited as ref. 19 also in the context of tordion mass formulas in PGT. Pre-1983 reference. Possibly relevant to the massless-tordion question.

### 4.5 F. W. Hehl, J. Nitsch, P. von der Heyde, in *General Relativity and Gravitation*, ed. A. Held (Plenum, 1980) — `/tmp/BN1983_PRD.txt:1913-1915`

Standard reference. Already in TIDAL's literature corpus presumably.

### 4.6 BN's promised "next paper" (no gauge fixing)

The PRD paper announces (`/tmp/BN1983_PRD.txt:96-97`): "The form of the theory without gauge fixing will be studied in our next paper."

This unnamed-at-time successor paper would be **the gauge-invariant general PGT Hamiltonian analysis**. Tracking it down would close the 1983 chapter of the BN programme. Likely candidates:
- M. Blagojević, I. A. Nikolić, "Hamiltonian structure of Poincaré gauge theory in general gauge", subsequent late-1980s paper(s).
- Per the PRD-paper reference list, ref. 15 is "M. Blagojevic and I. A. Nikolic, Institute of Physics, Belgrade, Report No. IF-83/19, 1983 (unpublished)" (`/tmp/BN1983_PRD.txt:1894-1896`) — this is plausibly the gauge-invariant follow-up. **Probably hard to retrieve (Belgrade Institute internal report)**, which is presumably why later authors (BC 2018) re-do the analysis from scratch.

**Conclusion of Section 4:** No genuinely missed references that would change the audit's verdict. Szczyrba 1982 and the BN unpublished IF-83/19 report are the only ones that might contain a constructive recipe the investigation hasn't seen, and IF-83/19 is functionally inaccessible.

---

## Section 5 — Verdict

### 5.1 Verdict: (b)-with-amendment — "R+T² scope" framing is incorrect; PRD covers R+T²+R²; but BN's recipe is still not directly transferable to TIDAL.

The audit's revised framing (post-Nuovo Cim mis-citation correction) reads BN 1983 as "R+T² only, hence background diagnostics not constructive bridge". **This is wrong about the scope.** The PRD companion paper of the same year covers R+T²+R² explicitly, including R² and `(εR)²` curvature-squared invariants — i.e. the same Lagrangian class as TIDAL's `b5·R̃²`.

However, the audit's *bottom-line* conclusion (BN 1983 is precursor diagnostics, not a constructive bridge for TIDAL's specific b5·R̃² critical-surface case) **survives a more careful reading**, because:

1. **BN 1983 PRD's constructive treatment is restricted to massive tordions** (`/tmp/BN1983_PRD.txt:1052-1062`). The vanishing-critical-parameter case where the would-be tordion is *massless* (the rank-jump-into-Lagrange-multiplier regime that TIDAL's b5 → 0 transition exhibits) is explicitly deferred to "more detailed considerations".
2. **Direction of constraint-promotion is opposite.** BN's "if-constraint becomes constraint when c → 0" is one direction; TIDAL's "constraint becomes propagating mode when b5 ≠ 0" is the inverse. The recipe cannot be transferred directly without inverting the construction.
3. **Time-gauge restriction** in BN 1983 PRD; TIDAL's analysis is gauge-flexible.
4. **Case-by-case dispatch** in BN's Table I scales poorly to multi-parameter degeneracies; for TIDAL's `b5·R̃²` the Round 1 result `det(M) ∝ b5⁶` indicates a 6-fold-degenerate critical surface, which BN's Table-I framework would not handle in closed form.

### 5.2 Required correction to the audit framing

The `FINAL_ASSESSMENT.md` text at lines 351-355 and the corresponding manual-retrieval-bookkeeping section need amendment:

**Replace:**

> "Blagojević-Nikolić 1983 (Nuovo Cim B 73:258, doi:10.1007/BF02721794) — earlier mis-cited as Nuovo Cim B 84:25 (1984). The actual paper is on R + T²~PGT (linear in R, quadratic in torsion), NOT R+R²+T² as previously claimed."

**With:**

> "The Blagojević-Nikolić 1983 corpus consists of two companion papers: (i) Nuovo Cim B **73**:258 (1983), doi:10.1007/BF02721794, on R+T² PGT only, as a stepping stone; (ii) **Phys. Rev. D 28:2455 (1983)**, doi:10.1103/PhysRevD.28.2455, on the *general* R+T²+R² PGT in the time gauge. The PRD paper covers TIDAL's Lagrangian class. It introduces the 'if-constraints with critical parameters' terminology and provides a constructive Hamiltonian for the *massive-tordion* regime. The *massless-tordion* case at vanishing-critical-parameter loci — which is where TIDAL's b5·R̃² rank-jump-into-Lagrange-multipliers phenomenon resides — is explicitly deferred in the PRD paper as needing 'more detailed considerations' (PRD ref. /tmp/BN1983_PRD.txt:1052-1062). The 'constraint-promotion' phenomenon and its 1983 historical anchor are therefore real, but BN 1983 PRD already flagged the specific massless-rank-jump TIDAL faces as outside its solved scope."

### 5.3 Recommendations for the canonical TeX writeup

For `docs/tex/perturbative_reduction_constraint_barrier.tex`:

1. **Cite the PRD companion paper.** The historical anchor is BN 1983 *PRD* (R+T²+R²), not the Nuovo Cim *(R+T² only)*. The latter should be cited as the methodological precursor; the former as the first treatment of the Lagrangian class TIDAL studies.

2. **Quote the 1983 deferral.** Adding the verbatim quote from `/tmp/BN1983_PRD.txt:1052-1062` ("In the case corresponding to the massless tordions [...] more detailed considerations are necessary for a complete understanding of the situation") is the strongest single piece of historical evidence that TIDAL's b5·R̃² massless-rank-jump case has been *acknowledged-as-open* for 40+ years. This is much stronger framing than "25-year-unsolved" cobbled from BC 2018 alone.

3. **Sharpen the directional inversion.** Both BN 1983 and BC 2018 treat "constraint promotion" in the direction `c → 0 ⇒ constraint appears`. TIDAL's b5 → 0 is the *inverse* direction (`b5 ≠ 0 ⇒ Lagrange-multiplier constraint disappears`). The TeX writeup should clearly distinguish these so the rhetorical move "BN 1983 said it's hard" is precisely what BN 1983 actually said — they said the *forward* direction in the *massless* regime is hard.

4. **Drop the "25-year" phrase in favour of "40-year-acknowledged".** The chain BN 1983 PRD → Yo-Nester-Ni 1999/2002 → BC 2018 is 35 years; rounding to "40-year" is fair. The original "25-year" came from a count starting at BC 2018; that's wrong.

5. **No new constructive recipe is unlocked.** Reading the BN 1983 PRD does not produce a method that Round 1's three Stückelberg/Vainberg-Tonti/Curtright paths overlooked. The case-by-case Table-I dispatch is a *narrower* tool than the auxiliary-field lifts already explored. The verdict (b) — precursor diagnostics, not constructive bridge — therefore survives the deeper reading.

---

## Section 6 — Citations

All quotes verbatim from the local `pdftotext` extracts. Primary citations:

| Quote | Source | Line(s) |
|---|---|---|
| "Hamiltonian Structure of the Theory of Gravity with R+T² Type of Lagrangian" (NCB title) | `/tmp/BN1983_NuovoCim.txt` | 7-8 |
| "linear in the scalar curvature and quadratic in the torsion" (NCB scope) | `/tmp/BN1983_NuovoCim.txt` | 16-21 |
| "ℒ_G = a R + ℒ_T ≡ ℒ_ET" (NCB eq. 6) | `/tmp/BN1983_NuovoCim.txt` | 147 |
| "ℒ_T = α t² + β v² + γ a²" (NCB eq. 7) | `/tmp/BN1983_NuovoCim.txt` | 170-189 |
| "Some of the results obtained here may be useful in an analysis of the general theory" (NCB scope) | `/tmp/BN1983_NuovoCim.txt` | 200-206 |
| "We will call these expressions if-constraints" (NCB terminology) | `/tmp/BN1983_NuovoCim.txt` | 341-344 |
| "λ(x) = 1/x for x≠0; 0 for x=0" (NCB eq. 22) | `/tmp/BN1983_NuovoCim.txt` | 402-407 |
| "we will limit ourselves to the case of the massive tordions" (NCB deferral) | `/tmp/BN1983_NuovoCim.txt` | 691-699 |
| "this will not necessarily be true in the general case" (NCB scope flag) | `/tmp/BN1983_NuovoCim.txt` | 733-736 |
| "we are making a step toward a complete understanding" (NCB stepping-stone framing) | `/tmp/BN1983_NuovoCim.txt` | 967-969 |
| "L = α_{ij} F^i F^j + β_i F^i + L(F^i = 0)" (NCB Appendix A construction) | `/tmp/BN1983_NuovoCim.txt` | 982-988 |
| "linearity of this theory in A^{ab}_0 does not persist" (NCB scope ack) | `/tmp/BN1983_NuovoCim.txt` | 1115-1116 |
| "Hamiltonian dynamics of Poincare gauge theory: General structure in the time gauge" (PRD title) | `/tmp/BN1983_PRD.txt` | 13 |
| "R+T²+R² type of Lagrangian" (PRD abstract) | `/tmp/BN1983_PRD.txt` | 28-47 |
| "The situation will change in the general R+T²+R² case" (PRD scope) | `/tmp/BN1983_PRD.txt` | 74-79 |
| "ℒ_R = b₁ R² + b₂ + … b₅ R² + b₆ (εR)²" (PRD eq. 2.6) | `/tmp/BN1983_PRD.txt` | 182-227 |
| "only five constants are independent, due to the Bach-Lanczos identi[ty]" (PRD basis) | `/tmp/BN1983_PRD.txt` | 228-229 |
| "form of the theory without gauge fixing will be studied in our next paper" (PRD followup) | `/tmp/BN1983_PRD.txt` | 96-97 |
| "ḣ_{a,0} and Ȧ^{ij}_0 are contained only in T_{ab0} and R_{ij,0}" (PRD velocity localisation) | `/tmp/BN1983_PRD.txt` | 430-436 |
| "give rise to primary constraints if some of the constants … vanish" (PRD if-constraints) | `/tmp/BN1983_PRD.txt` | 494-509 |
| "the variables M_0a^c and π_{ab}^c are independent Hamiltonian degrees of freedom" (PRD massive tordion) | `/tmp/BN1983_PRD.txt` | 973-1008 |
| "In the case corresponding to the massless tordions … more detailed considerations are necessary" (PRD deferral, MAIN VERDICT QUOTE) | `/tmp/BN1983_PRD.txt` | 1052-1062 |
| Table I: critical parameter combinations | `/tmp/BN1983_PRD.txt` | 1133-1217 |
| "W. Szczyrba, Phys. Rev. D 25, 2548 (1982)" (PRD ref. 9) | `/tmp/BN1983_PRD.txt` | 1874 |
| "Hayashi-Shirafuji series" (PRD ref. 12) | `/tmp/BN1983_PRD.txt` | 1881-1883 |
| "Sezgin–van Nieuwenhuizen 1980" (PRD ref. 18) | `/tmp/BN1983_PRD.txt` | 1907-1910 |
| "BN, Belgrade IF-83/19, unpublished" (PRD ref. 15) | `/tmp/BN1983_PRD.txt` | 1894-1896 |
| Nuovo Cim companion citation | `/tmp/BN1983_PRD.txt` | 1872-1873 |

Cross-references:
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/meta_reviews/meta_review_K_literature_claims.md`
- `/workspaces/torsion-gertsenshtein/research/perturbative_hamiltonian/meta_reviews/meta_review_N_BC_AppendixD_content.md`
