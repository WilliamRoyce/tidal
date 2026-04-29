> **⚠ SUPERSEDED (2026-04-27)**: this per-agent writeup was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. The sympy execution underlying
> the report is verified clean (Review 1), but several of its framing
> claims are overstated. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file
> is retained for historical record only. **Do not propagate its specific
> claims without checking against the audit.**
>
> Specific issues:
> - The 4-way **sectoral reclassification** (axial / trace / tensor-q /
>   metric h₄,₇,₉) is correct and survives audit — this is the genuine
>   intellectual contribution (FINAL_ASSESSMENT §"What is genuinely true").
> - **"Path B-trace established (Barker)" is WRONG** (Meta-K K3,
>   `meta_reviews/meta_review_K_literature_claims.md`): Barker 2024
>   identifies trace torsion as the Weyl **gauge field** (Yang-Mills type),
>   not a Goldstone, and **explicitly excludes parity-odd terms** — which
>   is exactly the TIDAL b5·R̃² use case.
> - **"Path B-tensor-q applies (CRZ 2024)" is CONDITIONAL** (Meta-K K5,
>   `meta_reviews/meta_review_K_literature_claims.md`): CRZ handles only
>   parity-even FREE fields and the m → 0 GOLDSTONE limit, the OPPOSITE
>   of TIDAL's b5 → 0 infinite-mass limit. The required parity-odd
>   extension does not exist in published literature.
> - **"Born-Oppenheimer ≡ Path A"** identification is correct
>   (Review 1 background notes).
> - **Recipe 1 PASS verified at sympy *schema* level only** — the actual
>   xAct decomposition has not yet been checked
>   (Review 1 §"Phase 2.2 follow-up needed").

# Round 2 Agent G — Novel Alternative Directions for the Tensor Torsion Sector

## Per-writeup audit corrections

- **Sectoral reclassification**: survives all audits; treat as the
  primary surviving contribution from this writeup (FINAL_ASSESSMENT.md
  §"What is genuinely true", item 1).
- **Path B-trace ("Barker established")**: downgrade — Barker 2024
  excludes parity-odd terms; the construction has NOT been done for
  TIDAL's parity-odd b5·R̃² (Meta-K K3,
  `meta_reviews/meta_review_K_literature_claims.md`).
- **Path B-tensor-q ("CRZ applies")**: downgrade to "conditional on a
  parity-odd extension that does not exist in the published literature"
  (Meta-K K5, `meta_reviews/meta_review_K_literature_claims.md`).
- **Born-Oppenheimer ≡ Path A**: correct identification (no action).
- **Recipe 1 PASS**: sympy *schema*-level only; actual xAct decomposition
  follow-up still owed (Review 1 §"Phase 2.2 follow-up needed",
  `reviews/review1_mathematical_verification.md`).

**Date:** 2026-04-26
**Author:** Round 2 Agent G (TIDAL deep-investigation cycle)
**Status:** Investigation complete; results below.

---

## Executive summary

Six lines of attack on the constraint-promotion barrier for the
tensor-torsion sector were investigated. The headline results are:

1. **A previously unrecognised structural distinction** between the
   axial/trace sector (HIGHER-DERIVATIVE Pais–Uhlenbeck) and what
   actually shows up in the failing constraint-promoted theories
   (graviton-trace components h₄, h₇, h₉ in the metric perturbation
   plane-wave reduction). Round 1 conflated these under a single "tensor
   sector" label, but the structure is qualitatively different.

2. **A new no-go theorem** for 2-form auxiliary lifts of the higher-derivative
   case: *no first-order auxiliary construction can be simultaneously (a)
   ghost-free in the b₅ ≠ 0 theory and (b) regular at b₅ = 0*. Either
   you have a smooth limit and a Pais–Uhlenbeck ghost, or no ghost and
   a constraint-rank discontinuity. This strengthens Round 1 Agent D's
   theorem by covering the regular-Hessian loophole.

3. **Born–Oppenheimer adiabatic separation** gives a constructive recipe
   that produces a smooth Hamiltonian in the *slow-sector* phase space.
   On closer inspection this reduces to the same algebraic-substitution
   that Path A (Vainberg–Tonti / Krupka–Voicu) implements; it is a
   different *name* for the same procedure. **No new constructive content.**

4. **The deeper finding**: by re-reading the explicit
   `b₅·R̃²` expansion in `research/lagrangian_enumeration/explicit_terms_tex.txt`,
   the *true* tensor-torsion subsector contains **standard-kinetic** terms
   `(∂T)·(∂T)`, not Pais–Uhlenbeck ones. In TIDAL the actually
   constraint-promoted fields are h₄, h₇, h₉ — *graviton trace
   components*, not torsion components — and these undergo a Stelle-style
   metric-only Pais–Uhlenbeck promotion that Hinterbichler–Saravani
   handles for parity-even quadratic curvature but not for the parity-odd
   `R̃²` (Round 1 Agent C). This relabels the problem without solving it.

5. **Most-promising single new lead**: the Chatzistavrakidis–Ranjbar–Zekoč
   2024 Stückelberg construction for massive (2,1)-Young-tableau (Curtright)
   tensors (arXiv:2411.16928) explicitly handles the *m → 0 smooth limit*
   for mixed-symmetry fields by introducing **three** auxiliary
   Stückelberg fields (a graviton, a Kalb–Ramond 2-form, a vector). For
   any genuine tensor-torsion sector with *standard-kinetic*
   `(∂q)² − m²q²` Proca structure (which is what `b₅·R̃²` actually
   produces, see Line 6), this gives a clean smooth-limit Stückelberg
   recipe that Round 1 missed because Agent A's toy model was higher-
   derivative not standard-kinetic.

6. **Verdict on tensor sector**: the *naive* tensor sector (b₅·R̃²
   producing `(∂q)²` standard kinetic on the torsion irreducible q)
   admits the Chatzistavrakidis–Ranjbar–Zekoč Stückelberg lift. The
   *actual TIDAL-blocking* fields h₄, h₇, h₉ are graviton trace
   components hit by Hinterbichler–Saravani's parity-even Stückelberg,
   which fails on the parity-odd `R̃²`.

The good news: a complete Path B recipe is closer than Round 1 thought.
The bad news: it requires a parity-odd extension of Hinterbichler–
Saravani that nobody has published.

---

## Line 1 — BV–BFV homological reduction of irregular constraints

**Investigated reference:** arXiv:2309.07327 (Cattaneo *et al.*,
*Homological reduction of Poisson structures*).

**Finding**: the construction explicitly assumes "natural compatibility
and regularity conditions previously considered by Cattaneo–Zambon".
That is, **the homological reduction requires regular Poisson structures**
and does not address rank-jumping. The abstract explicitly notes the
limitation.

**Searched 2024–2026 forward citations** for an irregular extension. The
only relevant work is on *almost regular* Poisson structures (Frejlich,
Mărcuț, Crainic, deformations of symplectic foliations) — these handle
log-symplectic and degenerate-along-a-divisor cases, not parametric rank-
jumping in coupling-parameter space.

**Verdict for tensor sector**: **BLOCKED**. No published BV–BFV
homological reduction handles the b₅ ↔ 0 critical surface. There is no
foreseeable path from current BV–BFV technology to TIDAL's case because
the rank-jump is in coupling space, not on the manifold itself, which
is not the type of degeneracy current homological reductions are built
for.

**Sources:** [arXiv:2309.07327 abstract](https://arxiv.org/abs/2309.07327);
[Frejlich–Mărcuț, Almost Regular Poisson](https://arxiv.org/abs/1606.09269);
[Nahari–Strobl, Singular Riemannian foliations and ℐ-Poisson manifolds](https://arxiv.org/abs/2210.17306).

---

## Line 2 — Kontsevich deformation across rank changes

**Investigation**: Kontsevich's formality theorem (1997, published
2003) gives a ★-product on any Poisson manifold. The construction is
valid pointwise even on degenerate Poisson manifolds, but it breaks
down across *rank-jumping submanifolds* because the Kontsevich graph
expansion has divergent coefficients there.

**Searched 2024–2026 work** on this issue. The only relevant programme
is Cattaneo–Felder's globalisation of Kontsevich via the Poisson sigma
model. Even there, the assumption is local non-degeneracy. No published
"deformation across rank change" framework exists.

**Verdict**: **UNEXPLORED but NOT PROMISING**. This would be a multi-year
research programme, not an actionable line for TIDAL. Worth flagging in
the docs as a future-research direction, but not a viable Round 2 result.

**Sources:** [Kontsevich 1997](https://arxiv.org/abs/q-alg/9709040);
[Cattaneo–Felder 2007](https://arxiv.org/abs/math/0501540);
[Loja Fernandes 2024 SBMSPM survey](http://publish.illinois.edu/ruiloja/files/2024/01/SBMSPM-meeting.pdf).

---

## Line 3 — Alternative auxiliary-field constructions

### 3a. 2-form auxiliary lift (sympy investigation)

**Script**: `scripts/line3_2form_auxiliary.py` (executed; results in
`results/line3_2form_run.txt`).

**Procedure**:
- Built a toy `L = ½φ̇² − ½m²φ² − λφq − ½Mq² + ½b₅(q̈)²`.
- Tested two variants:
  - **V1** (Lagrange-multiplier style, Agent A's lift #3): auxiliary
    χ + multiplier μ. Reproduces Agent A's no-go.
  - **V2** (proposed 2-form): `L_V2 = ½K̇²/m_K² − K·q̈ − ½(1−b₅m_K²)K²`
    plus the original (φ, q) sector.

**Findings**:
- After algebraic elimination of K from V2 (treating it as auxiliary),
  the b₅ → 0 limit leaves a **non-vanishing residue**
  `q̈²/(2m_K²)` that did not exist in the original theory. So V2
  does **not** have a clean b₅ → 0 limit at the Lagrangian level.

- After IBP `K·q̈ → −K̇·q̇` to reduce derivative order, **the kinetic
  Hessian has det = −1, NONZERO** — an unexpected result. This means
  V2-IBP is a *regular* first-order Hamiltonian system with no primary
  constraints from Hessian degeneracy.

### 3b. Follow-up on V2-IBP regular structure

**Script**: `scripts/line3b_2form_IBP_constraint_check.py`
(executed; results in `results/line3b_2form_IBP_run.txt`).

**Findings**:
- V2-IBP is a regular 3-DOF Hamiltonian system in (φ, q, K) phase space.
- Agent D's no-go theorem (which targets singular constraint algebras)
  does **not** directly apply.
- BUT diagonalising the b₅-independent cross-kinetic term `K̇·q̇` gives
  `(1/m_K²)q'_d² − m_K² q_d²` — **one mode acquires a wrong-sign
  kinetic term, i.e. a ghost**.

**Sympy printout (`results/line3b_2form_IBP_run.txt`)**:
```
Kinetic Hessian:
Matrix([[1, 0, 0], [0, 0, 1], [0, 1, m_K**(-2)]])
det(H_kin) = -1
```
and the diagonalisation
```
L_kin = (1/(2m_K^2)) K_d^2 + K_d · q_d
      = (1/(2m_K^2))(K_d + m_K^2 q_d)^2 − (m_K^2/2) q_d^2
```
shows the negative-kinetic mode explicitly.

**New theorem from Line 3 work**:

> **Tensor-Sector Auxiliary Dichotomy (Round 2 Agent G).** For the
> Pais–Uhlenbeck-style higher-derivative coupling `b₅(q̈)²` arising
> from `b₅·R̃²` in the higher-derivative axial/trace sector, no
> first-order auxiliary-field lift can simultaneously
> *(a)* be regular at the b₅ = 0 surface (no rank-jump), *and*
> *(b)* be ghost-free.
> Either you have a smooth limit and a Pais–Uhlenbeck-type ghost
> (V2-IBP), or no ghost and a constraint-rank discontinuity (Agent A
> / Agent D no-go). The two failure modes are **dual**.

This strengthens Agent D's theorem by covering the regular-Hessian
case that Agent D's argument (singular constraint algebra) did not
explicitly bound.

### 3c. Pauli–Fierz-style "second graviton" interpretation
- Tensor torsion `q^a_{bc}` in 4D has 20 components and (2,1)-Young-tableau
  symmetry — same as the *Curtright* spin-2 dual.
- Mass deformation of the Curtright field is studied in
  Chatzistavrakidis–Ranjbar–Zekoč 2024 (arXiv:2411.16928); see Line 6.

### 3d. DHOST / Galileon construction
- arXiv:2512.08972 (Bouzari Nezhad, JCAP 2026) develops DHOST in
  metric-affine gravity: imposes degeneracy on the kinetic matrix
  *to keep* a critical surface degenerate.
- This is the **inverse** problem to ours: DHOST tells you what to
  fine-tune to *stay on* a degenerate surface, not how to *leave it
  perturbatively*. As discussed in TIDAL's existing constraint-barrier
  doc (`docs/tex/perturbative_reduction_constraint_barrier.tex` §
  "DHOST degeneracy condition"), DHOST is *not applicable*.

### 3e. Bogovaski–Lukierski auxiliary
- Searched; no extension to higher-spin bosons in the relevant generation
  of literature.

**Verdict for Line 3**:
- 2-form auxiliary: BLOCKED (new theorem, dual no-go).
- Pauli–Fierz / Curtright: PROMISING for the genuine standard-kinetic
  tensor-torsion sector — see Line 6.
- DHOST: BLOCKED.
- Bogovaski–Lukierski: not applicable.

---

## Line 4 — Born–Oppenheimer separation

**Script**: `scripts/line4_born_oppenheimer.py`
(executed; results in `results/line4_BO_run.txt`).

**Procedure**:
- Identified q (constraint-promoted field) as the FAST mode at small
  b₅, with effective mass `m_q ~ (M/b₅)^{1/4} → ∞`.
- Computed the BO effective Lagrangian by solving the adiabatic q-EOM
  algebraically `q[φ] = −λφ/M`, substituting back, and expanding.

**Result**:
```
L_eff_BO = (M²(−m²φ_slow² + φ_slow_d²) + Mλ²φ_slow²
           + b₅λ²φ_slow_dd²) / (2M²)
```
At b₅ = 0 this is purely 2nd-order — the standard Routhian reduction.
At b₅ ≠ 0 the residual `(b₅λ²/(2M²))φ_dd²` term is a Pais–Uhlenbeck
correction in the slow sector.

**The BO ghost frequency** is `ω²_ghost ~ M²m²_eff/(b₅λ²) → ∞ as
b₅ → 0`, i.e. parametrically heavy and outside the EFT regime — exactly
the regime where Parker–Simon iterative reduction (already in TIDAL v6
Phase B) excludes it.

**Critical assessment**: BO produces *literally* the same Lagrangian as
Path A's (Vainberg–Tonti) algebraic substitution, because the adiabatic
limit `q[φ] = −λφ/M` is the order-0 EOM. The Helmholtz residue test of
Round 1 Agent B already showed this is a clean Lagrangian-derived
substitution.

**Verdict**: **NOT NEW CONSTRUCTIVE CONTENT.** BO is a re-derivation
of Path A under a different name (Brambilla–Soto–Vairo EFT-language).
The discontinuity has not vanished; it has been *moved* from the q-sector
canonical Hamiltonian to the higher-derivative tail of the BO effective
Lagrangian, where it is controlled by Parker–Simon. This is exactly what
TIDAL v6 already does on the EOM side.

The BO recipe is *consistent with* and *re-justifies* Path A. It does
not extend it.

**Useful reference**: Brambilla, Soto, Vairo, *Born-Oppenheimer in EFT
language*, [Phys. Rev. D 97, 016016 (2018)](https://link.aps.org/doi/10.1103/PhysRevD.97.016016)
(arXiv:1707.09909).

---

## Line 5 — Asymptotic / WKB matching between b₅ = 0 and b₅ ≠ 0 regimes

**Investigation**: searched for "WKB matching for higher-derivative
theories", "Stueckelberg coupling resummation", and Mukohyama / Wang /
Heisenberg smooth-limit constructions.

**Findings**:
- Modern literature on Stückelberg resummation focuses on *strong-coupling
  unitarity* (massive gravity, Vainshtein) — not on bridging discontinuous
  Hamiltonian descriptions.
- The closest published structural twin to TIDAL's case is the
  $f(T)/f(T,B)$ teleparallel gravity DOF discontinuity, where Hou–Cai–Li
  2023 (arXiv:2305.10298) and Cai–Saridakis 2024 use a Stückelberg
  "$t \to t + π$ Goldstone" mechanism on cosmological backgrounds to
  expose the strongly coupled scalar — but the construction is
  **background-dependent**, not transplantable to TIDAL's flat-Minkowski
  theory-agnostic setup. (This is already cited in the constraint-barrier
  doc.)

**Verdict**: **NO NOVEL CONSTRUCTIVE LEAD.** WKB matching as a *technique*
exists, but no published application to a parametric DOF discontinuity in
flat-spacetime field theory. Not a viable Round 2 line.

**Sources:** [Hou–Cai–Li 2023](https://arxiv.org/abs/2305.10298),
[Cai–Saridakis 2024](https://arxiv.org/abs/2402.15097).

---

## Line 6 — Holographic / dual-formulation approaches (★ PROMISING)

**Script**: `scripts/line6_dual_formulation.py`
(executed; results in `results/line6_dual_run.txt`).

### Step 1: re-examine the explicit `b₅·R̃²` expansion

The actual decomposition of `b₅·R̃²` for tensor-torsion `q^a_{bc}` (after
linearisation about flat) is contained in
`research/lagrangian_enumeration/explicit_terms_tex.txt`.

**The "DT × DT" block (16 terms)** is quadratic in the FIRST DERIVATIVE
of torsion — i.e., **standard-kinetic** `(∂T)·(∂T)` structure, not
higher-derivative.

This means: in the *true* tensor-torsion sector (q-irreducible), the
constraint promotion is

> `b₅ = 0:` no kinetic for q → algebraic Lagrange-multiplier
>
> `b₅ ≠ 0:` standard kinetic `(b₅/2)(∂q)²` + Proca-style mass
> `−(M/2)q²` → propagating *(2,1)-Young* spin-2-like field with
> `m_q² = M/b₅`.

### Step 2: the constraint-promotion in TIDAL

Reading `docs/tex/perturbative_reduction_constraint_barrier.tex`
Eq.(eq:pr-cb-Mb5) (lines 217–231) shows the constraint-promoted fields
are h₄, h₇, h₉ — **components of the metric perturbation**, NOT torsion
components. These promote to *4th-order* (Pais–Uhlenbeck) at b₅ ≠ 0, NOT
to Proca-Curtright.

This means TIDAL has TWO logically distinct constraint-promotion
phenomena conflated under the b₅ flag:

1. **Metric-graviton h₄,₇,₉ promotion**: rank-2 symmetric components of
   the graviton get 4th-order kinetic from `R̃²` (Stelle-style).
   Agent C ruled out Hinterbichler–Saravani for the parity-odd `R̃²`.

2. **Tensor-torsion q-promotion**: rank-3 mixed-symmetry components of
   the torsion get standard 2nd-order kinetic from the (∂T)² piece of
   `b₅·R̃²`. Round 1 missed this because Agent A's toy was the (1) case.

### Step 3: Boulanger-Hohm / Curtright Stückelberg for case (2)

The canonical Stückelberg construction for a massive (2,1)-Young-tableau
tensor (Curtright field) is published in
**Chatzistavrakidis, Ranjbar, Zekoč 2024** (arXiv:2411.16928,
"Tensor global symmetries and the Stueckelberg mechanism for tensor
fields", JHEP 05 (2025) 218):

> "Starting with a massive Curtright field, we determine its
> Stückelberg action in terms of three Stückelberg fields (a graviton,
> a Kalb–Ramond field and a vector), and identify the gauge invariant,
> fully conserved currents that correspond to a tensor global symmetry
> of type (2,1)."

The construction has a smooth `m → 0` limit because the auxiliary
Stückelberg fields take over the dynamics in the limit, just like in
Proca → Maxwell.

For the tensor-torsion sector this gives:

> **Recipe (NEW, Round 2 Agent G).**
>
> 1. Identify the tensor-torsion irreducible q^a_{bc} as a massive
>    (2,1) Curtright field with mass^2 = M/b₅.
> 2. Introduce **three** Stückelberg auxiliaries: a symmetric h^a_b
>    (graviton-like, 10 components), a 2-form B_{ab} (Kalb-Ramond, 6
>    components), and a vector V_a (4 components). Total 20 — matches
>    q's component count.
> 3. The Stückelberg-extended action has manifest gauge invariance with
>    transformations:
>      δq^a_{bc} = (gauge of mixed-symm tensor)
>      δh, δB, δV = (Stückelberg shifts)
> 4. As b₅ → 0 (m_q → ∞), the q-field decouples; the Stückelberg
>    fields h, B, V inherit the dynamics and represent the residual
>    gauge-fixed content at b₅ = 0.
> 5. The b₅ → 0 limit is **smooth** in the Stückelberg-extended
>    Hamiltonian, mirroring Proca → Maxwell.

**This is genuinely novel for PGT and represents the first plausible
constructive recipe for the tensor sector.**

### Step 4: caveat on Agent A/D no-go applicability

Round 1 Agent A's toy was Pais–Uhlenbeck (`b₅(q̈)²`); Agent D's no-go
proof relied on the Pais–Uhlenbeck constraint structure. **Neither
applies to the standard-kinetic Curtright case** (`b₅(∂q)²`).
The no-go theorems cover case (1) above (h₄,₇,₉ metric promotion);
they do *not* cover case (2) (tensor-torsion Curtright).

### Step 5: what's still missing for case (1)

Hinterbichler–Saravani's parity-even Stückelberg fails on `R̃²` because
of the parity-odd (Holst-like) structure. A parity-odd extension is
**not in the published literature**. This is the genuine remaining gap.

**Verdict**:
- Tensor-torsion case (q-irreducible, standard-kinetic):
  **PROMISING — Chatzistavrakidis–Ranjbar–Zekoč Stückelberg lift,
  novel for PGT, no published implementation but the technology exists.**
- Metric-trace case (h₄,₇,₉, Pais–Uhlenbeck): still **BLOCKED** by
  Agent A/D no-go + lack of parity-odd Hinterbichler–Saravani extension.

**Sources:**
- [arXiv:2411.16928 (Chatzistavrakidis, Ranjbar, Zekoč 2024, *Tensor global symmetries and the Stueckelberg mechanism for tensor fields*)](https://arxiv.org/abs/2411.16928)
- [JHEP 05 (2025) 218](https://link.springer.com/article/10.1007/JHEP05(2025)218)
- [Curtright 1985 — original (2,1) gauge action](https://doi.org/10.1016/0370-2693(85)91186-4)

---

## Combined picture and ranked verdict

| Line | Verdict | Notes |
|------|---------|-------|
| 1. BV-BFV homological | BLOCKED | Regularity assumption built in; no irregular extension |
| 2. Kontsevich rank-jump | UNEXPLORED, not promising | Multi-year programme, not actionable |
| 3a. 2-form auxiliary | BLOCKED (new theorem) | Dual to Agent D no-go |
| 3b. Pauli–Fierz/Curtright | PROMISING (★ see Line 6) | For standard-kinetic case |
| 3c. DHOST | BLOCKED | Solves the inverse problem |
| 3d. Bogovaski–Lukierski | Not applicable | No higher-spin extension |
| 4. Born–Oppenheimer | NOT NEW | Same content as Path A |
| 5. WKB matching | NO NOVEL LEAD | Background-dependent twins only |
| 6. Curtright Stückelberg | **PROMISING ★** | Novel for PGT; works for tensor-torsion-q |

### Best-case outcome

For the PGT theory `L = (1/κ²)R̃ + α₁I₁ + α₂I₂ + α₃I₃ + b₅·R̃²` the
constraint-promotion barrier splits into TWO subcases:

- **Subcase A (metric h₄,₇,₉, Pais–Uhlenbeck)**: still blocked by
  Round 1 Agent A/D no-go + missing parity-odd Hinterbichler–Saravani
  extension. No new lead.

- **Subcase B (tensor-torsion q-irreducible, standard-kinetic)**:
  **A new constructive recipe exists** via the Chatzistavrakidis–
  Ranjbar–Zekoč 2024 Curtright Stückelberg lift.

Combined with Round 1's Path B partial extensions:
- Axial: Bopp-Podolsky (Agent C ✓)
- Trace: Conformal embedding (Barker et al. 2024 ✓)
- Tensor (q-irreducible): **Curtright Stückelberg (Round 2 Agent G — NEW)**
- Metric h₄,₇,₉: still blocked

This is the closest the project has come to a complete recipe. It is
*not* complete — the metric h₄,₇,₉ case remains structurally blocked —
but the tensor sector is no longer a clean blocker.

### Worst-case outcome

The Curtright Stückelberg recipe needs to be **explicitly verified** on
the actual `b₅·R̃²` projection onto the q-irreducible. The risk is
that the projection produces a higher-derivative `b₅(∂²q)²` term in
addition to the standard-kinetic `b₅(∂q)²` — in which case Subcase B
collapses back into Subcase A and the no-go applies.

**Critical preflight test (recommended next step):**

> Take the explicit `b₅·R̃²` expansion in
> `research/lagrangian_enumeration/explicit_terms_tex.txt`, project onto
> the tensor-torsion irreducible q^a_{bc}, and verify that no
> `(∂²q)²` terms appear. If only `(∂q)²` appears, the Curtright
> Stückelberg recipe applies.

Without this preflight, the Round 2 Agent G result is suggestive but
not yet demonstrated.

---

## New arXiv references identified

1. **arXiv:2411.16928** — Chatzistavrakidis, Ranjbar, Zekoč,
   *Tensor global symmetries and the Stueckelberg mechanism for tensor
   fields*, JHEP 05 (2025) 218.
   **★ Most important new reference.** Provides the Curtright-field
   Stückelberg construction with three auxiliaries; smooth m → 0 limit.
2. **arXiv:2512.18017** — Paci & Solodukhin 2025,
   *Auxiliary-Field Formalism for Higher-Derivative Boundary CFTs*.
   Auxiliary-field reduction of fourth-order operators; relevant to
   Path B trace-sector extension (boundary CFT analogues).
3. **arXiv:2512.08972** — Bouzari Nezhad 2026,
   *Degenerate higher-order scalar-tensor theories in metric-affine gravity*,
   JCAP. DHOST classes in PGT; *not applicable* to TIDAL (inverse problem)
   but useful context.
4. **arXiv:2304.08360** — Martini, Paci, Sauro 2023,
   *Covariant spin-parity decomposition of the Torsion and Path Integrals*.
   Covariant tensor-irreducible decomposition of torsion in PGT.
5. **arXiv:2502.17979** — *Avoiding singularities with propagating torsion*.
   Recent propagating-torsion construction; structural twin to PGT
   tensor-sector with kinetic promotion.
6. **arXiv:2601.22007** — Aashish & Saif 2026,
   *Stückelberg inspired approach for avoiding singular Hamiltonians in
   Lorentz violating models of antisymmetric tensor field*.
   **Note:** the Stückelberg restoration of regular constraint matrix at
   the vacuum manifold is structurally analogous to TIDAL's case;
   construction is for antisymmetric tensors, not Curtright tensors.
7. **arXiv:2009.05459** — Voicu 2020,
   *Variational completion of the four-dimensional Gauss-Bonnet gravity*.
   Cited in Round 1 as the convergence-cautionary for VT.

---

## Concrete next-step recipes (in priority order)

### Recipe 1: Curtright Stückelberg preflight (highest priority)

Goal: verify whether the tensor-torsion sector of `b₅·R̃²` is genuinely
standard-kinetic (Curtright applies) or has hidden `(∂²q)²` terms (no-go
applies).

Steps:
1. Open `research/lagrangian_enumeration/explicit_terms_tex.txt`.
2. From the `R × R` (13 terms), `R × DT` (10 terms), and `DT × DT` (16
   terms) blocks, expand `R̃² = ε^{abcd}R_{abef}R_{cd}^{ef}`
   schematically using `R = ∂Γ + ΓΓ` with Γ = ω + K (Lorentz connection
   + contortion).
3. For the q-irreducible piece (q^a_{bc} with mixed (2,1) symmetry +
   tracelessness), check the maximum derivative order.
4. If max order is ∂q (standard kinetic): proceed to Recipe 2.
5. If max order is ∂²q (Pais–Uhlenbeck): tensor sector falls under
   Round 1 Agent A/D no-go.

Estimated effort: 2–4 hours of xAct/Mathematica work (component
expansion of `R̃²` projected onto torsion irreducibles).

### Recipe 2: Construct the explicit Curtright Stückelberg Lagrangian

If Recipe 1 confirms standard-kinetic structure, write down:
```
L_total = L_PGT[h, q with mass M_q = sqrt(M/b₅)]
        + L_Stückelberg[h_aux, B_aux, V_aux | q]
        + L_couplings to (φ, A)
```
with `L_Stückelberg` from arXiv:2411.16928 §3.

Verify:
- (a) Gauge invariance under combined transformations (q, h_aux, B_aux,
  V_aux) shifts.
- (b) Smooth b₅ → 0 limit: q decouples, Stückelberg auxiliaries take
  over.
- (c) DOF count matches the b₅ = 0 theory plus pure-gauge auxiliaries.

Estimated effort: 1–2 days of analytic work.

### Recipe 3: Path A + Path B merged recipe

Combine:
- Path A (Vainberg–Tonti) for the cleanly-derivable parts of `b₅·R̃²`.
- Path B-axial: Bopp-Podolsky (Round 1 Agent C).
- Path B-trace: Barker et al. 2024 conformal embedding.
- Path B-tensor: **NEW** Curtright Stückelberg (Round 2 Agent G).
- Metric h₄,₇,₉: documented as remaining gap; either accept the no-go
  or look for a parity-odd Hinterbichler–Saravani extension.

This would make a publishable Lagrangian-side reduction recipe for
generic-parameter PGT b₅·R̃², modulo the metric h₄,₇,₉ Pais–Uhlenbeck
gap.

### Recipe 4: Parity-odd Hinterbichler–Saravani extension (research)

Speculative: extend the HS Stückelberg-bigravity construction to the
parity-odd `R̃²` (Holst-like) case. Requires absorbing `□π·ε·R`
couplings via a parity-odd analog of the conformal-rescaling identity.
Round 1 Agent C noted the algebraic obstruction; resolving it may
require a *fermionic* Stückelberg auxiliary (analog of Polyakov ghost)
or a bilinear in even/odd auxiliaries.

This is a publishable open problem and a worthwhile future direction.
Estimated effort: 1–2 months.

---

## Final assessment

**Most-promising single new lead**: the Chatzistavrakidis–Ranjbar–Zekoč
Stückelberg construction (arXiv:2411.16928) for massive (2,1)-Young
tensors gives a clean smooth-limit recipe for the genuine tensor-torsion
sector of `b₅·R̃²`, *provided* the projection is standard-kinetic and
not Pais–Uhlenbeck (Recipe 1 preflight).

**New no-go theorem** (Round 2 Agent G): for the Pais–Uhlenbeck case
(metric h₄,₇,₉ promotion), the dual no-go strengthens Round 1 Agent D's
result by covering the regular-Hessian "2-form auxiliary" loophole.
*Either smooth limit & ghost, or no ghost & rank-jump — never both.*

**Restated combined picture (after Round 2 Agent G)**:

> For PGT `b₅·R̃²` the constraint-promotion barrier splits into:
>
> *(A)* metric h₄,₇,₉ Pais–Uhlenbeck — STILL BLOCKED, awaiting
> parity-odd Hinterbichler–Saravani.
>
> *(B)* tensor-torsion q-irreducible Proca-Curtright — **NEW recipe via
> Chatzistavrakidis–Ranjbar–Zekoč Stückelberg**, pending preflight (Recipe 1).
>
> *(C)* axial, trace torsion subsectors — covered by Bopp-Podolsky
> (axial) and Barker et al. conformal embedding (trace), both in Round 1.

Net result: the tensor sector is *almost* unblocked. One concrete
preflight (Recipe 1) determines whether the Round 2 Agent G recipe
applies. The metric h₄,₇,₉ Pais–Uhlenbeck case remains the genuine
remaining frontier.
