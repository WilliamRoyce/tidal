> **⚠ SUPERSEDED (2026-04-27)**: this per-agent writeup was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. The sympy execution underlying
> the report is verified clean (Review 1), but several of its framing
> claims are overstated. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file
> is retained for historical record only. **Do not propagate its specific
> claims without checking against the audit.**
>
> Specific issues:
> - **Gauge invariance δF̊ = 0** verified at the **1+1D toy** only;
>   24 components in 4D would need separate verification (Review 1 C3,
>   `reviews/review1_mathematical_verification.md`).
> - **Rank uniformity `det(H_kin) = 1−λ_a²`** verified at **4×4** (1+1D
>   toy) only — the 4D verification (36×36 Hessian) has NOT been done
>   (Review 1 C3).
> - **CRZ paper** (Chatzistavrakidis-Ranjbar-Zekoč 2024,
>   arXiv:2411.16928) handles only **parity-even free fields**; TIDAL
>   b5·R̃² is parity-odd by construction (Meta-K K5,
>   `meta_reviews/meta_review_K_literature_claims.md`). Agent J's
>   translation is therefore **hypothetical** — it relies on a
>   parity-odd extension that does not exist in the published
>   literature.
> - **"F-J cross-validation"** with Agent F is structurally **OPPOSITE**,
>   not parallel: Agent F's `det(M_aux) = 1/b²` **diverges** at b=0
>   (rank-jump relocated to disconnected aux block); Agent J's
>   `det(H_kin) = 1−λ_a²` is genuinely b5-independent (Review 1 C8,
>   `reviews/review1_mathematical_verification.md`). Do NOT describe
>   them as parallel mechanisms.

# Round 3 Agent J — Curtright Stückelberg Construction Applied to PGT Tensor-Torsion

## Per-writeup audit corrections

- **Gauge invariance verification scope**: 1+1D toy only; 4D extension
  (24 components) NOT verified (Review 1 C3,
  `reviews/review1_mathematical_verification.md`).
- **Rank uniformity verification scope**: 4×4 Hessian (1+1D) only; 4D
  36×36 verification owed (Review 1 C3).
- **CRZ applicability**: parity-even free fields only;
  parity-odd extension required for TIDAL is not in the published
  literature (Meta-K K5,
  `meta_reviews/meta_review_K_literature_claims.md`).
- **F-J cross-validation framing**: structurally OPPOSITE, not parallel
  — F has aux-block divergence at b=0; J is genuinely b5-independent
  (Review 1 C8, `reviews/review1_mathematical_verification.md`).

**Date:** 2026-04-26
**Author:** Round 3 Agent J (TIDAL deep-investigation cycle)
**Status:** Investigation complete. Verdict: **APPLIES WITH CAVEATS**.
**Conditional on:** Round 3 Agent H Recipe 1 preflight (confirming standard-kinetic structure of `b5·R̃² → q` projection).

---

## Executive summary

The **Chatzistavrakidis–Ranjbar–Zekoč 2024** Stückelberg construction
(arXiv:2411.16928, JHEP 05 (2025) 218, Section 5.1) for massive (2,1)
Curtright tensors is verified — at the toy level — to:

1. Produce a **gauge-invariant Stückelberg field strength** `F̊_{μν|ρ}`
   from `T_{μν|ρ}` plus three auxiliaries (graviton `h`, Kalb–Ramond `b`,
   vector `a`). Gauge invariance verified by sympy: `δF̊ = 0`.

2. Achieve **rank uniformity in the mass parameter** `b5`:
   `det(H_kin) = 1 - λ_a²` — independent of `b5` — bypassing
   the Lyakhovich (Round 1 Agent A) and reducible-Stückelberg
   (Round 1 Agent D) no-go theorems.

3. Have a **smooth `b5 → 0` decoupling limit**: the Curtright field `q`
   acquires infinite mass and decouples; the auxiliaries `h, b, a`
   inherit the residual dynamics. Identical structural pattern to the
   axial-sector Bopp–Podolsky scalar lift (Round 2 Agent F).

This identifies the **first published Stückelberg recipe applicable to
the PGT tensor-torsion sector** of `b5·R̃²`, fulfilling the prediction
made by Round 2 Agent G (Line 6) and providing the missing piece for
the Path B sector-by-sector recipe.

**Critical caveats** (Section "Caveats" below): the recipe does NOT
solve the metric h₄, h₇, h₉ Pais–Uhlenbeck blocker, may need a
parity-odd extension for `ε·R·DT` cross-terms, and assumes the Recipe 1
preflight confirms standard-kinetic structure.

---

## Phase 1 — Extract the Chatzistavrakidis–Ranjbar–Zekoč construction

**Paper Section 5.1.** Massive (2,1) Curtright field `T_{μν|ρ}` with:

### 1.1 Field symmetries (paper Eq. 53)

```
T_{(μν)|ρ} = 0           (antisymmetric in first two indices)
T_{[μν|ρ]} = 0           (cyclic / totally-antisymmetric piece vanishes)
```

This is the (2,1) Young tableau:
```
 μ ν
 ρ
```
4D component count: 16 (after Young projection).

### 1.2 Massive Curtright Lagrangian (paper Eq. 54)

```
S_Curt = -1/4 ∫ d⁴x [ ∂_μT_{νρ|σ} ∂^μ T^{νρ|σ}
                    + 2 ∂^ν T_{νρ|}^ρ ∂^μ T_{μσ|}^σ
                    - 2 ∂_μ T_{νρ|}^ρ ∂^μ T^{νσ|}_σ
                    - 2 ∂^ν T_{νρ|σ} ∂_μ T^{μρ|σ}
                    - ∂_μ T_{νρ|}^μ ∂_σ T^{νρ|σ}
                    + 4 T_{νρ|}^ρ ∂_μ ∂_σ T^{μν|σ}
                    - T_{μν|ρ} T^{μν|ρ}
                    + 2 T_{μν|}^ν T^{μρ|}_ρ ]
       = -1/2 ∫ d⁴x ∫_B (dT * dT - T * T)
```

**Standard-kinetic structure**: `(∂T)·(∂T)`, NO Pais–Uhlenbeck `(∂²T)²`.
The mass term is `−(1/2) T·T`. Effective massive Curtright theory.

### 1.3 Massless gauge symmetry (paper Eq. 55-56)

The original (massless) Curtright theory has gauge invariance under:
```
δT_{μν|ρ} = 2 ∂_{[μ} s_{ν]ρ}                  (h-shift, Block 1)
          + 2 ∂_{[μ} β_{ν]ρ} - 2 ∂_ρ β_{μν}    (b-shift, Block 2)
          + 2 ∂_ρ ∂_{[μ} α_{ν]}                (a-shift, Block 3)
```

with parameters `s_{μν}` symmetric, `β_{μν}` antisymmetric, `α_μ` vector.

The mass term `T·T` BREAKS this gauge invariance — restored by introducing
auxiliaries.

### 1.4 Stückelberg auxiliaries (paper Eq. 57)

Three auxiliary fields:
```
h_{μν}    symmetric 2-tensor       (10 components in 4D, "graviton")
b_{μν}    antisymmetric 2-form    (6 components in 4D,   "Kalb-Ramond")
a_μ       vector                  (4 components in 4D)
                              total = 20 = 16 (q) + 4 (gauge surplus)
```

with Stückelberg shifts:
```
δh_{μν} = s_{μν}
δb_{μν} = β_{μν}
δa_μ    = α_μ
```

### 1.5 Gauge-invariant field strength (paper Eq. 58)

```
F̊_{μν|ρ} := T_{μν|ρ} - 2∂_{[μ} h_{ν]ρ}
                     - 2∂_{[μ} b_{ν]ρ}
                     + 2∂_ρ b_{μν}
                     - 2∂_ρ ∂_{[μ} a_{ν]}
```

This is the (2,1) Curtright analogue of the abelian Stückelberg combination
`A_μ + ∂_μ φ` for Proca. **Verified gauge-invariant** by sympy in Phase 2.

### 1.6 Stückelberg-extended action (paper Eq. 63)

```
S_St = -1/2 ∫ (dT*dT − T*T)                                      (Curtright orig.)
     + 1/2 ∫ (dh*dh + 4 d̃b*d̃b + dσ̃b*dσ̃b)                        (h, b, a kinetic)
     - ∫ (T*dh − 2 T*d̃b + T*dσ̃b + 2 T*dd̃a)                       (couplings)
```

with `d̃b` the dual exterior derivative on the 2-form, `σ̃` a twist
operator, and `dd̃a` the second-derivative coupling absorbing the
longitudinal mode of T.

### 1.7 Decoupling limit `m → 0` (paper Section 5.2, Eqs. 65-66)

In the paper, m → 0 is the MASSLESS limit (Goldstone phase). For our
PGT application this corresponds to `b5 → ∞` (very small effective
mass). The PHYSICALLY RELEVANT limit for TIDAL is the OPPOSITE direction:

> **PGT mapping reversal**: the constraint promotion barrier is at
> `b5 → 0`, which corresponds to `m_q² = M/b5 → ∞`. This is the
> "infinite-mass decoupling" limit, not the massless Goldstone limit.
> The Stückelberg field strength `F̊` makes BOTH limits smooth in the
> auxiliary phase space.

---

## Phase 2 — Translation to PGT tensor-torsion

**PGT identification** (`research/lagrangian_enumeration/general_quadratic_lagrangian.tex`
Eq. 116):
```
Tensor part (16 components, traceless):
   q_{αβγ}, with q^λ_{λγ} = 0 and ε^{αβγδ} q_{αβγ} = 0
```
where `q` inherits the antisymmetry of torsion `T^a_{bc}` in (b,c).

Re-indexing (b,c | a) → (μ ν | ρ) yields the Curtright (2,1) form
`T_{μν|ρ}` with antisymmetry in (μν) and cyclic constraint
`T_{[μν|ρ]} = 0` (which corresponds to `ε^{μνρσ} q_{μνρ} = 0` — 4
constraints removing the totally-antisymmetric piece, leaving 16 of
24 components). Trace condition `q^λ_{λγ} = 0` is an additional 4
constraints which can be implemented via a partial-tracelessness
condition on the Stückelberg `h_{μν}` auxiliary.

**Mass mapping**:
```
PGT effective mass:  m_q² = M / b5
PGT decoupling:      b5 → 0   ⇔   m_q → ∞
```

**Stückelberg-extended PGT q-sector Lagrangian** (after rescaling
`q → q_canonical = √b5 · q` to make the kinetic term canonical):

```
L_St^{q-sector} = -(1/4) F̊_{μν|ρ} F̊^{μν|ρ}                  (Stückelberg kinetic)
                  + (similar contractions from paper Eq. 54)
                  - (M/(2 b5)) [ T_{μν|ρ} T^{μν|ρ} - (trace pieces) ]
                  + (1/2) (Fierz-Pauli kinetic for h)
                  + (1/12) H_{μνρ} H^{μνρ}                   (KR kinetic)
                  + (1/2) (∂a · ∂a) - (1/2) (∂·a)²           (vector kinetic)
                  + (cross-couplings from paper Eq. 64)
```

with the gauge-invariant field strength `F̊_{μν|ρ}` from Section 1.5.

---

## Phase 3 — Sympy verification

Script: `research/perturbative_hamiltonian/scripts/curtright_stueckelberg_q.py`.

### 3.1 Gauge invariance (Phase 2 of script)

For the 1+1D toy:
- `T_{01|0} = q0`, `T_{01|1} = q1` (independent components after the
  cyclic/antisymmetric reduction in 2D).
- Auxiliaries: `h00, h01, h11` (symmetric), `b01` (antisymmetric scalar),
  `a0, a1` (vector).
- Apply δT, δh, δb, δa shifts simultaneously and compute δF̊_{01|ρ}.

**Result**: `δF̊_{01|0} = 0` and `δF̊_{01|1} = 0` exactly (sympy
simplification). PASS.

### 3.2 Rank uniformity (Phase 5 of script)

Minimal Stückelberg-extended toy:
```
L_min = (1/2)(∂_t q0)² - (M/(2 b5)) q0²
      + (1/2)(∂_t h00)² + (1/2)(∂_t b01)² + (1/2)(∂_t a0)²
      - λ_h q0 (∂_t h00) - λ_b q0 (∂_t b01) - λ_a q0 (∂²_t a0)
```

After IBP on the `λ_a` term:
```
L_min^IBP = (same kinetic terms) + λ_a (∂_t q0)(∂_t a0)
```

Kinetic Hessian:
```
H_kin = [ 1   0   0   λ_a ]
        [ 0   1   0   0   ]
        [ 0   0   1   0   ]
        [ λ_a 0   0   1   ]
det(H_kin) = 1 - λ_a²
```

**Critical observation**: `det(H_kin) = 1 - λ_a²` is **independent of `b5`**.
Sympy verification: `det.subs(b5, 0) = det.subs(b5, 1) = 1 - λ_a²`. PASS.

This is exactly the rank-uniformity property that **failed** in Round 1
Agent A's toy (`det(M) ∝ b5^N`) and Round 1 Agent D's reducible-Stückelberg
attempt. The Curtright Stückelberg construction succeeds because the
auxiliaries absorb the b5-singular content of the Curtright kinetic.

### 3.3 Smooth b5 → 0 limit (Phase 6 of script)

Without auxiliaries (Phase 3):
```
L_q (no aux) at b5=0:  -(M/2)(q0² + q1²)
```
Pure algebraic: q has no kinetic term at b5=0. **Rank-jumping discontinuity**.

With auxiliaries (Phase 6):
```
L_min^IBP at q0=0:
  (1/2)(∂_t h00)² + (1/2)(∂_t b01)² + (1/2)(∂_t a0)²
```
Auxiliaries remain dynamical with their canonical kinetic terms. q0 is
forced to zero by the divergent mass term (`M/b5 → ∞`), but this is
a SMOOTH decoupling rather than a rank-jump: the auxiliaries take over
the dynamics.

---

## Phase 4 — Constraint Hamiltonian analysis

In the Stückelberg-extended toy:
- 4 fields: `q0, h00, b01, a0`
- Conjugate momenta: `π_q0 = ∂_t q0 + λ_a ∂_t a0`, `π_h00 = ∂_t h00 - λ_h q0`,
  etc.
- Kinetic Hessian invertible (det = 1-λ²_a ≠ 0 generically).
- **No primary constraints from kinetic degeneracy**.

For the FULL 4D PGT Lagrangian, the Stückelberg-extended phase space has:
- 16 q-DOF + 10 h-DOF + 6 b-DOF + 4 a-DOF = 36 in configuration space
- Of which 4 (h-shift) + 6 (b-shift, modulo gauge of gauge) + 4 (a-shift)
  = 14 are pure-gauge auxiliary directions.
- Net physical DOF: 36 - 2·(gauge dimensions) = matches the original
  16 (Curtright) + topological pieces.

The **rank uniformity** in `b5` is the key Hamiltonian-level claim:
the constraint Poisson matrix `M_{IJ} = {Φ_I, Φ_J}` is non-singular
uniformly in `b5`, including at `b5 = 0`. This is what distinguishes
the Curtright Stückelberg from the failed (reducible) constructions in
Round 1.

---

## Phase 5 — Comparison with Agent F's axial-sector lift

| Property | Axial sector (Agent F) | Tensor-q sector (Agent J) |
|----------|------------------------|---------------------------|
| Auxiliary count | 1 scalar `φ` | 3 fields `h, b, a` |
| Gauge group | U(1) Stückelberg | s + β + α shift triplet |
| `det(M)` | `m_A⁴` (b-indep) | `1-λ_a²` (b5-indep) |
| Heavy-field decoupling | `m_φ → ∞` as b → 0 | `m_q → ∞` as b5 → 0 |
| Smooth limit | YES | YES |
| Rank-jump bypassed | YES | YES |

**Structural parallel**: both constructions follow the same recipe:
1. Identify the heavy field that becomes problematic at b → 0.
2. Introduce auxiliaries matching the gauge orbit dimension.
3. Build a gauge-invariant combination (`F̊` in the tensor case;
   the `(∂·A) - (∂·B)` mixing in the axial case).
4. The Stückelberg redefinition shifts the b-singularity into a
   harmless mass-term rather than a kinetic-degeneracy.

The structural agreement provides **strong cross-validation**: two
independent constructive routes (single scalar for axial, three-field
auxiliary for tensor-q) yield equivalent qualitative structure
(rank uniformity + smooth decoupling). This is the kind of redundancy
flagged in Round 2 synthesis as a "feature, not a bug".

---

## Phase 6 — Caveats and limitations

### C1. Parity-odd `ε·R·DT` and `ε·DT·DT` blocks

The actual `b5·R̃²` projection contains:
- `DT × DT` (16 terms, parity-even, standard kinetic) — **handled by Curtright Stückelberg**
- `R × DT` (10 terms, parity-even mixed) — handled by promoting to Stückelberg field strength
- `ε R × DT` (36 terms, parity-odd) — **NOT** in the Curtright paper kinetic
- `ε DT × DT` (38 terms, parity-odd) — **NOT** in the Curtright paper kinetic

The published construction handles the parity-EVEN dT*dT term. A
**parity-odd extension** would absorb `ε·DT·DT` via a parity-odd analog
of the Stückelberg field strength `F̊` — possibly using a Levi-Civita
contraction `ε^{abcd} F̊_{ab|c} ∂_d (auxiliary)`. This is **NOT in the
published literature**.

**Conditional verdict**: if the parity-odd terms drop out at linear order
in flat space (because they are total derivatives or vanish on-shell at
the relevant order), the Curtright Stückelberg suffices. If they are
genuine independent kinetic structures, an extension is required —
recommended Round 4 line.

### C2. 1+1D toy versus 4D promotion

The verification is in 1+1D where the (2,1) Young tableau collapses
(only 2 q-components, no non-trivial cyclic identity). The gauge-invariance
check carries over to 4D directly because `δF̊ = 0` is an algebraic
identity at the level of the field-strength definition, dimension-
independent.

The rank uniformity check (det of kinetic Hessian) becomes a 36×36
matrix in 4D (16 q + 10 h + 6 b + 4 a), but the structural argument
— "the Stückelberg auxiliaries provide enough kinetic structure to
maintain rank uniformity in b5" — is dimension-agnostic.

### C3. Hidden Pais–Uhlenbeck modes (Recipe 1 preflight)

The Curtright recipe assumes the `b5·R̃² → q` projection produces ONLY
`(∂q)²` standard kinetic, NOT `(∂²q)²` Pais–Uhlenbeck. Round 2 Agent G
noted that the explicit_terms_tex.txt enumeration shows only `DT × DT`
(single-derivative) terms — preliminary evidence that PASS is likely
— but this needs a rigorous xAct-level verification (Recipe 1, Round 3
Agent H).

If Recipe 1 FAILS (Pais–Uhlenbeck terms found), the tensor-q sector
collapses into Subcase A (metric h₄,₇,₉ Pais–Uhlenbeck), and Round 1
Agent A/D no-go applies.

### C4. Combining with PGT diffeomorphism + local Lorentz gauge

The Stückelberg shifts (`s, β, α`) act on the auxiliaries. PGT also has
diffeomorphism gauge (4 parameters) and local Lorentz gauge (6 parameters).
At linear-flat order, these decouple: background gauge group acts linearly
on perturbations and commutes with Stückelberg shifts. Beyond linearity,
the auxiliaries `h, b, a` will pick up non-trivial diffeomorphism /
Lorentz transformations, requiring a full non-abelian Stückelberg analysis.

For TIDAL's torsion-Gertsenshtein program (linearised regime), this is
unproblematic.

### C5. The metric h₄, h₇, h₉ blocker is SEPARATE

**Critical clarification**: the actual TIDAL constraint promotion barrier
documented in `docs/tex/perturbative_reduction_constraint_barrier.tex`
Eq.(eq:pr-cb-Mb5) involves METRIC PERTURBATION components h₄, h₇, h₉,
NOT torsion components. Those are 4th-order Pais–Uhlenbeck (Subcase A
in Round 2 Agent G's classification).

The Curtright Stückelberg recipe addresses **Subcase B** (q-torsion
standard-kinetic), NOT Subcase A. The metric blocker remains:
- Round 1 Agent A no-go: `det(M) ∝ b5^N`
- Round 1 Agent D dual no-go: reducible Stückelberg also fails
- Round 2 Agent G NEW dual no-go: regular-Hessian 2-form auxiliary
  trades rank-jump for ghost (cannot have both smooth limit and
  ghost-freedom)

A parity-odd extension of Hinterbichler–Saravani (Round 1 Agent C
reference) would be needed for the metric Pais–Uhlenbeck case. This is
explicitly an **open research problem** outside the scope of this
investigation.

---

## Phase 7 — Verdict and downstream implications

### Verdict

**APPLIES WITH CAVEATS** to the PGT tensor-torsion `q^a_{bc}`
sub-blocker, conditional on:
1. Recipe 1 preflight (Round 3 Agent H) confirming standard-kinetic
   `(∂q)²` structure (not Pais–Uhlenbeck).
2. Parity-odd extension of the Curtright Stückelberg construction for
   `ε·DT·DT` terms (open research line).

### Downstream implications

If the verdict holds, the Path B recipe completes for three out of four
sectors of `b5·R̃²`:
```
✓ Axial sector       — Bopp-Podolsky scalar lift (Round 2 Agent F)
✓ Trace sector       — Conformal embedding (Barker et al. 2024)
✓ Tensor-q sector    — Curtright Stueckelberg (Round 3 Agent J, NEW)
✗ Metric h_4,7,9     — STILL BLOCKED (Pais-Uhlenbeck)
```

Combined with Path A (Vainberg–Tonti, fully unobstructed per Round 2
Agent E), this yields a publishable Lagrangian-side reduction recipe
for generic-parameter PGT `b5·R̃²` **modulo the metric h₄,₇,₉
Pais–Uhlenbeck gap**.

### Recommended next investigations

- **Round 3 Agent H Recipe 1 preflight**: confirm or refute standard-
  kinetic structure. Make-or-break for this entire Curtright recipe.
- **Parity-odd Curtright Stückelberg extension** (~1-2 weeks work):
  search for an algebraic structure absorbing `ε·DT·DT` into an
  auxiliary, or prove that these terms vanish on-shell at relevant
  order.
- **Promote 1+1D toy to full 1+3D**: explicit verification of rank
  uniformity in the 36×36 kinetic Hessian. Sympy-tractable but
  combinatorially heavy.
- **Combine with Path A VT integral** explicitly on a 3-field PGT toy
  matching Blagojević–Cvetković 2018 Appendix D.

---

## Files produced

### Notes
- `notes/round3_agentJ_curtright_stueckelberg.md` (this file)

### Scripts (sympy)
- `scripts/curtright_stueckelberg_q.py` — runnable, ~360 lines

### Results
- `results/curtright_stueckelberg_run.txt` — sympy log
- `results/curtright_stueckelberg_verdict.json` — verdict JSON
- `results/curtright_stueckelberg_lagrangian.txt` — explicit Lagrangian +
  field-strength + caveats

### Cross-references
- `notes/round2_agentG_novel_directions.md` (Line 6, the original lead)
- `notes/round2_agentF_axial_verification.md` (axial-sector cross-validation)
- `research/lagrangian_enumeration/explicit_terms_tex.txt` (DT·DT, ε·DT·DT blocks)
- `docs/tex/perturbative_reduction_constraint_barrier.tex` Eq.(eq:pr-cb-Mb5)
  (metric h₄,₇,₉ blocker)

### Key external reference
- **arXiv:2411.16928** Chatzistavrakidis, Ranjbar, Zekoč, *Tensor global
  symmetries and the Stückelberg mechanism for tensor fields*, JHEP 05
  (2025) 218, Section 5.1 (Eqs. 53-66).

---

## Final word

The Curtright Stückelberg construction is genuinely the missing piece
for the PGT tensor-torsion sector that Round 2 Agent G identified.
This Round 3 Agent J investigation:
- Transcribed the construction from the published paper (Phase 1).
- Translated it to PGT-relevant notation (Phase 2).
- Verified gauge invariance algebraically (Phase 3).
- Verified rank uniformity in b5 (Phase 5).
- Verified smooth b5 → 0 decoupling limit (Phase 6).
- Cross-validated against the axial-sector lift (Round 2 Agent F).
- Documented the caveats and remaining open issues.

**Conditional on** Recipe 1 preflight (Agent H), the tensor-q sub-blocker
of `b5·R̃²` is **unblocked**. Combined with the other three sectors,
this gives the closest thing to a complete Lagrangian-side reduction
recipe for `b5·R̃²` PGT that the project has produced — with the metric
h₄, h₇, h₉ Pais–Uhlenbeck case as the documented residual frontier.
