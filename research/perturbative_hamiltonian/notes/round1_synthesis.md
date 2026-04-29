> **⚠ SUPERSEDED (2026-04-27)**: this synthesis was audited by Reviews 1-3 and
> Meta-Reviews K/L/M/N. Several of its claims are overstated, miscited, or
> wrong. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file is
> retained for historical record only. **Do not propagate its specific claims
> without checking against the audit.**
>
> Specific issues: see Reviews 1-3 in `reviews/` and Meta-Reviews K-N in
> `meta_reviews/`. Headline overstatements: KV "Theorem 1" doesn't exist;
> "Path B-trace established" is wrong (Barker excludes parity-odd); F-J
> cross-validation is qualitative not strong; L_VT inherits PU structure;
> NEW M_c²→0 pathology not flagged; phase-space jump factor 3 vs 5 is a
> convention drift.

# Round 1 Synthesis (2026-04-26)

Four parallel agents investigated the constraint-promotion barrier in
PGT b5·R̃². Verbatim conclusions below.

## Agent A — Lyakhovich's general Stückelberg recipe

**Toy model**:
```
L = ½(∂_tφ)² − ½m²φ² − λφh − ½Mh² + ½b5(∂_t²h)²
```

**Three Stückelberg lifts attempted**:
1. Naive `(∂_t²h − ψ)²` with `δh=ε, δψ=∂_t²ε` → gauge non-invariant; fix needs non-local `∂_t⁻²` (locality failure)
2. Auxiliary `ψ` with `½b5·ψ²` → degenerates to original 4th-order theory, no Ostrogradsky reduction
3. Lagrange-multiplier Ostrogradsky `L̃ = L|_{b5=0} + ½b5χ² − μ(χ − ∂_t²h)` → succeeds at Lagrangian level, terminates at order ξ², b5→0 decoupling clean

**CRITICAL NEGATIVE**: Even with lift #3, the 4×4 constraint Poisson matrix in extended phase space has **det(M) ∝ b5**. The Lyakhovich lift NAMES the discontinuity but does NOT bridge it.

**Recommendation**: investigate Lyakhovich-Sharapov (reducible generators) → Agent D.

## Agent B — Helmholtz residue / Vainberg-Tonti

**T2 falsifier toy** (with derivative-mixing cross-coupling):
```
L_T2 = ½(∂_tφ)² + ½(∂_tχ)² − ½m_φ²φ² − ½m_χ²χ² − λφh − λ_χ χ ∂_t h − ½h² + ½b5(∂_t²h)²
```

**Result**: Helmholtz residues `H^k_{φχ}, H^k_{χφ}, H^k_{φφ}, H^k_{χχ}` for k = 0..6 — **ALL ZERO**.

**Why T2 isn't a falsifier**: h is purely auxiliary at order 0; substitution of `h⁽⁰⁾[φ, χ̇]` into L is exact Routhian reduction; the constraint is holonomic.

**True falsifiers require non-Lagrangian sources** (Case D: γχ̇ injected directly into E_h with no Lagrangian preimage gives nonzero residue at k=1, 5). **TIDAL's EOMs are NEVER non-Lagrangian** — they come from a Wolfram-derived Lagrangian by construction.

**Verdict**: **TIDAL PGT direction (d) preflight PASSES on Helmholtz grounds.** The remaining gate is **VT integral convergence** (Voicu 2020 4D-GB cautionary).

## Agent C — Hinterbichler-Saravani extension to PGT torsion

**HS algebraic obstruction**: HS Eq. 5.5 (conformal identity for Einstein tensor under `g→e^{-2M²π}g`) absorbs all π-dependence. **No parity-odd analog for R̃** — R̃ transforms anomalously under conformal rescalings (topological term), giving `□π·ε·R` couplings with no clean reabsorption.

**Aoki-Mukohyama 4D obstruction** (line 385 of arXiv:2009.11739):
> "the Einstein-Hilbert is non-linear in the vielbein in four dimensions"

In 3D: ξ∧R is linear → clean f∧R for second vielbein. In 4D: requires quadratic ξ → infinite tower R(1+R/αm²)⁻¹R. Plus: generic PGT has FOUR independent mass scales; AM handles only ONE.

**Constructive partial extension — AXIAL torsion sector**:

The parity-odd Holst contracts cleanly with axial vector `A_μ = (1/6)ε^νρσ_μ T_νρσ`:
```
R̃ ⊃ ∂·A,  R̃² ⊃ (∂·A)²
```

Bopp-Podolsky-style lift:
```
L_aux = -¼F^(A)_μν F^(A)μν − (b5/2)B_μ B^μ + (b5/2)G_μν F^(A)μν − ½m_A² A_μ A^μ
```

with `F^(A) = dA, G = dB`. EOMs:
- δA^ν: `∂_μ F^(A)μν − b5 ∂_μ G^μν − m_A² A^ν = 0`
- δB^ν: `B^ν = ∂_μ F^(A)μν`

**Clean b5→0 limit** (B decouples, becomes free non-propagating, can be discarded). **NO rank-jump in A-sector** because Poisson matrix block-diagonalizes: A-block has constant rank, B-sector is degenerate but DISCONNECTED.

**Sector verdicts**:
- Axial: ✅ Bopp-Podolsky single-auxiliary lift works at linear-flat order
- Trace: ✅ Conformal embedding via Barker et al. 2024 (parity-even)
- Tensor: ❌ blocked, needs reducible Stückelberg

**Critical correction**: Bopp-Podolsky structure applies to AXIAL (via R̃⊃∂·A), NOT trace-only b5·R̃² (which vanishes by ε-tensor index symmetry in trace-only sector at linear-flat order).

## Agent D — Abakumova-Lyakhovich reducible Stückelberg

**Reference correction**: arXiv:2106.09355 (not 2107.08240).

**Three convergent arguments → DEFINITIVE NO-GO**:

1. sympy on 1-promoted-field toy: `det(M) = b5²`. Reproduces Agent A.
2. sympy on 2-promoted-field toy with mixed coupling `b5·ḧ₁·ḧ₂`: `det(M) = b5⁴`.
3. Structural impossibility theorem: reducibility null-vectors `Z^a` are b5-independent by construction (built from order-0 EOM consequence-generator structure). A b5-independent gauge structure cannot transform a b5-dependent Poisson bracket into a b5-independent one. QED.

**What Abakumova-Lyakhovich actually does**: duality-construction tool. Equivalence at fixed parameter values, NOT interpolation across critical surfaces. The 2022 Nordström follow-up confirms same usage pattern.

**For PGT b5·R̃²**: with three constraint-promoted fields, the (χ_a, μ_a) sector contributes a 6×6 block-diagonal constraint matrix with **det ∝ b5⁶**. Lorentz/diffeo gauge generators do not enter the (χ_a, μ_a) bracket.

**Verdict**: Stückelberg lifting in any form (irreducible or reducible) is genuinely BLOCKED for generic PGT b5·R̃² as a way to bridge the b5=0 critical surface.

## Combined picture: TWO complementary viable paths

### Path A — Krupka-Voicu / Vainberg-Tonti direct on full PGT
- Agent B: δE = 0 generically when constraint and EOM derive from the same Lagrangian (TIDAL's pipeline ✓)
- Remaining test: VT integral convergence (Voicu 2020 4D-GB pathology check)
- If converges: publishable Lagrangian-side reduction recipe

### Path B — Sector-by-sector partial extension
- Axial: Bopp-Podolsky single-auxiliary (Agent C explicit Lagrangian)
- Trace: Barker et al. conformal embedding (parity-even)
- Tensor: blocked

These paths are NOT mutually exclusive. The most likely outcome is a hybrid: VT for parts where it converges, sector-by-sector for the rest.

## Most important pending investigations

1. **VT integral convergence on a minimal PGT toy** (the make-or-break for Path A)
2. **Independent verification of Agent C's Bopp-Podolsky axial-sector construction** (derive EOMs from scratch, check b5→0 limit, constraint structure)
3. **Novel alternative directions** for the tensor-torsion sector specifically (BV-BFV homological, Kontsevich deformation, alternative auxiliary fields)
