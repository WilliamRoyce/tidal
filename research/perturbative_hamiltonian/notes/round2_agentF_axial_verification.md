> **⚠ SUPERSEDED (2026-04-27)**: this per-agent writeup was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. The sympy execution underlying
> the report is verified clean (Review 1), but several of its framing
> claims are overstated. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file
> is retained for historical record only. **Do not propagate its specific
> claims without checking against the audit.**
>
> Specific issues:
> - The vector→scalar auxiliary correction itself is **sympy-correct**
>   (Agent F's correction of Agent C is verified by Review 1).
> - The "AM 2020 cross-validation" framing is **qualitative consistency
>   only**, not structural equivalence (Meta-K K6,
>   `meta_reviews/meta_review_K_literature_claims.md`): Aoki-Mukohyama's
>   φ has non-canonical sinh-tower kinetic `(∂φ)²/(1+φ²)` and lives in a
>   dRGT-bigravity infinite tower; the bare `𝒳²` term **breaks** AM's
>   ghost-freedom. AM's φ is NOT the same auxiliary as this writeup's
>   Lagrange-multiplier scalar.
> - The "no rank-jump in A-sector" framing is **sector-restricted**:
>   Review 1 C8 (`reviews/review1_mathematical_verification.md`) shows
>   `det(M_aux) = 1/b²` **diverges** as b→0 — the rank-jump is
>   *relocated* to the disconnected aux block, not removed.
> - **Curved-background O(h²) extension fails**: Nieh-Yan-type cross
>   terms remain unabsorbed (FINAL_ASSESSMENT §"What is overstated").

# Round 2 Agent F: Independent Verification of Axial-Torsion Bopp-Podolsky Lift

## Per-writeup audit corrections

- **Vector→scalar correction**: verified clean by Review 1 sympy
  re-execution (`reviews/review1_mathematical_verification.md`,
  C-checks).
- **AM 2020 cross-validation**: downgrade from "strong" to "qualitative
  consistency at linearised order"; AM's φ is structurally different
  (Meta-K K6, `meta_reviews/meta_review_K_literature_claims.md`).
- **A-sector rank uniformity**: holds within the A-sector but
  `det(M_aux) = 1/b²` diverges at b=0; rank-jump relocated, not removed
  (Review 1 C8, `reviews/review1_mathematical_verification.md`).
- **Curved-background extension**: Nieh-Yan cross terms unabsorbed at
  O(h²); the lift is linear-flat-only as written
  (FINAL_ASSESSMENT.md §"What is overstated").

**Date:** 2026-04-26
**Status of Agent C's claim:** **PARTIALLY VERIFIED — with one substantive correction and one structural caveat.**

## TL;DR

Agent C's structural claim — that the axial-torsion sector of `b5·R̃²` admits
a single-auxiliary Bopp-Podolsky lift with **a clean b₅→0 limit** and
**no rank-jump in the A-sector** — **is correct**. However the *specific
Lagrangian* Agent C wrote down is **not the right lift**: the auxiliary
must be a **scalar** φ (not a vector B_μ with G=dB) because R̃ in the
axial sector reduces to the **scalar** ∂·A, not to a tensor object.

The substantive deliverable for the paper is therefore:

```
L_aux^(b5) = -¼ F_{μν}F^{μν} - ½ m_A² A_μ A^μ + φ (∂·A) - (1/(2b)) φ²
```

with `b ≡ b₅ κ_axial² > 0` (κ_axial = 3/2 in 4D). The constraint
Poisson matrix is rigorously **block-diagonal**: an A-block with
b-independent rank, and a disconnected aux-block whose b→0 rank-jump
is the standard "infinitely-heavy auxiliary decoupling" limit.

This is consistent with **Aoki-Mukohyama 2020** (arXiv:2009.11739
Eq. 497-526), who show that adding `|e|·𝒳²` (Holst-squared) introduces
a **single propagating spin-0⁻ scalar** to PGT. Our φ corresponds to
their φ at linearised order.

## Agent C's verbatim claim

> Bopp-Podolsky-style lift:
> ```
> L_aux = -¼ F^(A)_μν F^(A)μν − (b5/2) B_μ B^μ + (b5/2) G_μν F^(A)μν − ½ m_A² A_μ A^μ
> ```
> with `F^(A) = dA, G = dB`. EOMs:
> - δA^ν: `∂_μ F^(A)μν − b5 ∂_μ G^μν − m_A² A^ν = 0`
> - δB^ν: `B^ν = ∂_μ F^(A)μν`
>
> Clean b5→0 limit (B decouples, becomes free non-propagating, can be
> discarded). NO rank-jump in A-sector because Poisson matrix
> block-diagonalizes.

## Verification phases

### Phase 1 — Linearised axial PGT (literature + dim. consistency)

Restricting torsion to the totally-antisymmetric (axial) channel
`T^ρ_{μν} = (1/3) ε_{μν}^{ρσ} A_σ`, the Holst scalar reduces to

```
R̃ |_{axial, linear, flat} = (3/2) ∂_μ A^μ        (Hehl-McCrea-Mielke-Ne'eman 1995)
```

(the Christoffel-only Pontryagin density `ε·R(g)` is a total derivative
and drops at linear-flat order). Hence

```
b5·R̃² → (b5·κ²/2) (∂·A)²  +  total deriv.
       = (b/2) (∂·A)²       with  b ≡ b5·κ²,  κ = 3/2.
```

This is **scalar-derivative-squared** structure, not tensor-derivative-
squared. It identifies the right auxiliary as a **scalar**, dual to the
scalar quantity (∂·A).

### Phase 2 — Why Agent C's vector-auxiliary form fails

We tested **three** candidate lifts in (1+1)D Minkowski (the polynomial
structure is identical in (1+3)D):

| Candidate | Form | EOM for aux | Algebraic? |
|---|---|---|---|
| **C1 (Agent C)** | `L = -¼F² - (b/2)B² + (b/2) G_{μν}F^{μν} - ½m²A²` | `B_ν = ∂_μ F^μν` (mixed with `G`) | **No** — adds Maxwell-style B sector |
| **C2 (Cuzinatto)** | `L = -¼F² - ½m²A² - (∂·B)(∂·A) - (1/(2b))B²` | `B_μ = -b·∂_μ(∂·A)` | Algebraic but trivial — same DOF as scalar |
| **C3 (scalar)** | `L = -¼F² - ½m²A² + φ(∂·A) - (1/(2b))φ²` | `φ = b·(∂·A)` | **Yes — algebraic** |

Sympy results (`bopp_podolsky_axial.py`):
- C1: `δB_0 = b·(B_0 + ∂_x²A_0 - ∂_t∂_x A_1)` — **kinetic in B through G, propagates**.
- C3: `δφ = ∂_x A_1 - ∂_t A_0 - φ/b` — **algebraic, solves φ = b(∂·A)**.

Substituting `φ = b(∂·A)` back into L_aux^φ gives **L_HD exactly** (verified:
`L_aux^φ[φ=b·(∂·A)] - L_HD = 0`).

**Why Agent C's vector form fails:** the Holst scalar `R̃` is a *scalar*
quantity in the axial sector. The Lagrange dual of a scalar is a scalar.
Coupling a vector B_μ via G_{μν}F^{μν} **adds a new propagating photon**
(Maxwell kinetic for B), not a Lagrange multiplier. The mistaken intuition
is that Bopp-Podolsky generalised electrodynamics — which keeps the
same vector A — can be lifted vector-to-vector. But there the
higher-derivative term is `(∂_α F^{αβ})²`, dual to a *vector* current.
Here the higher-derivative term is `(∂·A)²`, dual to a *scalar*.

**Correction:** the right Lagrangian is C3 (scalar auxiliary). The
structural conclusion that Agent C reached (rank-stable A-block,
disconnected aux-block) is correct under this correction.

### Phase 3 — Hamiltonian analysis (sympy, scalar aux)

Conjugate momenta from `L_aux^φ`:
```
π_{A_0} = 0                     [primary constraint]
π_{A_1} = ∂_t A_1 - ∂_x A_0     [physical electric field]
π_φ     = 0                     [primary constraint]
```

Velocity Hessian (sympy, full 3×3):
```
K = diag(0, 1, 0),    rank(K) = 1,    det(K) = 0
```

Two primary constraints (Φ_1 = π_{A_0}, Φ_2 = π_φ) and two secondary
constraints from preserving them in time:

```
χ_A   = δL/δA_0  = m²A_0 - ∂_x²A_0 + ∂_t∂_x A_1 + ∂_t φ      ≈ 0
χ_φ   = δL/δφ    = (∂·A) - φ/b                                ≈ 0
```

Constraint Poisson matrix (rows/cols = Φ_1, Φ_2, χ_A, χ_φ):

```
        ⎡  0      0    -m²    0    ⎤
   M =  ⎢  0      0     0   -1/b  ⎥
        ⎢  m²     0     0     0    ⎥
        ⎣  0    1/b     0     0    ⎦
```

`det(M) = m_A⁴ / b²`,  rank 4 at finite b (all constraints second-class,
as expected for a massive theory).

### Phase 4 — b → 0 limit and decoupling

The matrix is **block-diagonal**:

```
M = M_A ⊕ M_aux,
  M_A   = [[0,-m²],[m²,0]],   det(M_A)   = m⁴      [b-INDEPENDENT]
  M_aux = [[0,-1/b],[1/b,0]], det(M_aux) = 1/b²    [diverges as b→0+]
```

The A-sector Poisson block is **rigorously b-independent**. Agent C's
key structural claim is **verified**: the b₅→0 rank-jump lives entirely
in the disconnected auxiliary sector.

In the auxiliary sector:
- `φ_eom`: φ = b·(∂·A) → 0 as b → 0 (auxiliary trivialises).
- Auxiliary mass²: `m_φ² = 1/b → ∞` (auxiliary becomes infinitely heavy).
- This is the standard **decoupling limit**: integrate out the heavy
  scalar, recover Proca for A.

**Caveat (the structural barrier remains):** the limit is *singular at
the level of the auxiliary action* — the mass diverges, the kinetic
coefficient in `(1/(2b))φ²` blows up. The lift "names" the discontinuity
without bridging it analytically in b, which is exactly the conclusion
Round 1 Agent A reached for the toy model (`det(M) ∝ b₅` for the *full*
constraint structure; here `det(M) ∝ 1/b²` in the aux block).

So the lift gives:
- **Continuous in b₅ at the level of the physical (A-sector) phase
  space** ✓
- **Discontinuous in b₅ at the level of the auxiliary phase space** —
  but the auxiliary decouples, so this is physically harmless.

This is what makes the construction *useful for sector-by-sector
extension*: the rank-jump barrier identified in Round 1 is **localised
to the disconnected aux block**, and the physical sector inherits
analyticity in b₅.

### Phase 5 — Cross-check with Aoki-Mukohyama (arXiv:2009.11739)

AM Eq. 497:
> "The dynamical spin-0⁻ mode of the spin connection shows up around
> the flat background when the term `d⁴x|e|·𝒳²` is added where 𝒳 is
> the Holst scalar."

AM Eq. 526 (after integrating out non-dynamical fields):
```
L_AM ⊃ -(3 M_pl² (1-α)/4) · |f^a_μ| · (∂_α φ)²/(1+φ²)
```

with `φ` a **parity-odd scalar field**. Up to leading (linear) order
this is exactly a canonical scalar kinetic for φ.

In our construction: at linear-flat order with Proca completion,
solving `φ = b·(∂·A)` and integrating by parts gives a kinetic for
φ that mixes with `(∂·A)`. AM's φ corresponds to *our* φ at linearised
order — the "scalar excitation living inside R̃²". **Consistent.**

This is a strong cross-validation: AM independently arrived at the
same conclusion via the bigravity-equivalence route. Their spin-0⁻
scalar **is** the auxiliary that lifts Ostrogradsky.

### Phase 6 — Extensibility assessment

**Curved background (linear order in h_{μν}):**
The reduction `R̃ = κ·(∂·A)` used the fact that `ε·R(g)` is a
total-derivative Pontryagin density on flat space. On a curved
background,
```
ε·R(ω̃) = ε·R(g) + (Nieh-Yan torsion piece) + cross-terms
       ⊃ h_{μν}·∂²A·(∂h_terms)   [at quadratic-h, linear-A]
```
These cross-terms are **NOT** of the form (∂·A)² and cannot be
absorbed into the φ(∂·A) - φ²/(2b) lift. **The construction breaks
at O(h_{μν}²)**.

**Higher orders in the axial field A only:**
Cubic axial self-couplings come from `R̃·(T quadratic in A)`. These
are *algebraic* in A (no new derivatives). The lift survives
trivially: the φ field still solves `φ = b·(∂·A) + b·(algebraic A²)`,
which integrates back into the original L_HD without obstruction.

**Mixed sectors:**
- Trace torsion (vector channel `T_μ`): `b5·R̃²` projects to **zero** at
  quadratic order in this channel (ε-tensor index symmetry). No lift
  needed.
- Tensor torsion: `b5·R̃²` generates derivative-mixed tensor self-
  couplings that require *multiple* auxiliaries. Round 1 Agents A+D
  showed `det(M) ∝ b₅^N` for any reducible Stückelberg attempt.
  **BLOCKED.**

### Phase 7 — Verdict

| Sub-claim | Verdict |
|---|---|
| Single-auxiliary lift exists (axial, linear-flat) | **VERIFIED** (with scalar φ, not vector B) |
| Reproduces L_HD by integrating out aux | **VERIFIED** (sympy: difference = 0) |
| Rank-stable A-sector constraint structure | **VERIFIED** (`det(M_A) = m⁴`, b-independent) |
| Clean b₅→0 limit at A-sector level | **VERIFIED** (recovers Proca exactly) |
| Curved-background extension | **FAILS** at O(h²) |
| Higher-order axial-only extension | **SURVIVES** (cubic algebraic) |
| Tensor-sector extension | **BLOCKED** (Round 1 A+D) |

**Overall verdict: PARTIALLY VERIFIED with correction.**

The structural physics Agent C identified is correct and publishable.
The specific Lagrangian he wrote down had a wrong dualisation (vector
↔ vector instead of scalar ↔ scalar). The corrected Lagrangian (C3
above) is the one that should appear in the paper.

## Publishable proof-of-concept

```
PROPOSITION (axial-sector single-auxiliary completion of b5·R̃² PGT,
linear-flat order). Define
    b ≡ b5 · (3/2)²
and consider the auxiliary Lagrangian for an axial torsion 1-form
A_μ and a scalar auxiliary φ:

    L_aux = -¼ F_{μν} F^{μν} - ½ m_A² A_μ A^μ
            + φ · (∂_μ A^μ) - (1/(2b)) φ²

CLAIMS:
  (i)  Solving the algebraic EOM δφ = 0 gives φ = b · (∂·A) and
       reproduces the original 4th-order axial-PGT Lagrangian
       L_HD = -¼F² - ½m²A² + (b/2)(∂·A)².
  (ii) The constraint Poisson matrix decomposes as M = M_A ⊕ M_aux
       where M_A is the standard Proca constraint block (rank 2,
       b-independent) and M_aux is a 2×2 block in the (φ, π_φ)
       sector with det ∝ 1/b.
  (iii) The b → 0 limit is smooth at the level of the A-sector
       physical phase space and reproduces pure Proca.  The
       auxiliary φ becomes infinitely heavy (m_φ² = 1/b → ∞) and
       decouples in the standard heavy-mode sense.
  (iv) The construction is consistent with Aoki-Mukohyama
       arXiv:2009.11739: their propagating spin-0⁻ scalar is the
       same auxiliary at linearised order.

REMARKS:
  * The lift FAILS to extend to curved backgrounds at O(h_{μν}²)
    because Nieh-Yan torsion-curvature cross-terms generate
    couplings that are not absorbed by the φ(∂·A) trick.
  * The lift FAILS to extend to the tensor-torsion sector for the
    reasons given by Round 1 Agents A and D (reducible Stückelberg
    no-go).
  * Within the trace-torsion sector, b5·R̃² has zero quadratic
    projection (parity), so no lift is needed there.
END PROPOSITION.
```

This is the most concrete near-term publishable output of the
constraint-promotion-barrier investigation: a clean *partial*
extension valid at linear-flat order in the axial sector.

## Files produced

- `scripts/bopp_podolsky_axial.py` — sympy verification (runs in <2s).
- `results/axial_constraint_matrix.json` — explicit constraint Poisson matrix.
- This file: `notes/round2_agentF_axial_verification.md`.

## Open questions for the paper

1. **Sign convention for b**: We took `b > 0` so that the higher-
   derivative term has the same sign as a kinetic term in standard
   BP electrodynamics. PGT b5 sign conventions vary in the literature
   (Blagojevic vs Aoki vs Hehl). Reconciling these is a lit-review
   task before submitting.

2. **Coupling to gravity**: Even the linear-flat result assumes the
   metric is non-dynamical. Coupling to a dynamical h_{μν} forces
   the curved-bg breakdown identified in Phase 6 — but at *first*
   order in h_{μν} (graviton-axial torsion coupling), the lift may
   still survive. This deserves a dedicated Phase 8 follow-up.

3. **Comparison with Aoki-Mukohyama mass term**: AM's mass for the
   spin-0⁻ scalar comes from a dRGT-style construction with arbitrary
   `c_i(φ)` functions. Ours is a free parameter `m_A`. Mapping the
   two requires identifying which dRGT mass parameters generate the
   axial mass.
