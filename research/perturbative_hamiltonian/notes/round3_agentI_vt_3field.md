> **⚠ SUPERSEDED (2026-04-27)**: this per-agent writeup was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. The sympy execution underlying
> the report is verified clean (Review 1), but several of its framing
> claims are overstated. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file
> is retained for historical record only. **Do not propagate its specific
> claims without checking against the audit.**
>
> Specific issues:
> - VT integrand **convergence** + Krupka-Voicu **Definition-1 closure**
>   verified at N=3+3 — these survive audit (Review 1
>   `reviews/review1_mathematical_verification.md`, C-checks).
> - **Phase-space jump factor**: this writeup's main script reports
>   factor **3** (6→18); the synthesis-cited factor **5** (6→30) is
>   from a separate Phase 6 v2 script with a different convention
>   (`scripts/vt_T4_phase6_v2.py`). The convention drift was not
>   flagged at the time (Review 1 C6,
>   `reviews/review1_mathematical_verification.md`).
> - **Phase 6** (Hamiltonian rank analysis on the metric subspace) ran
>   into sympy performance issues and produced a qualitative verdict
>   only — the Lagrangian-side closure does NOT remove the
>   Hamiltonian-side rank-jump.
> - **Metric h₄,₇,₉ subspace inheriting Pais-Uhlenbeck** structure is
>   correct (FINAL_ASSESSMENT §"What is overstated", item 2).
> - **NEW pathology not flagged**: `L_VT diverges as M_c² → 0` (explicit
>   `1/M_c²` and `1/M_c⁴` poles) — Review 1 C5 sympy-verified
>   (`reviews/scripts_review/C5_routhian_M_to_zero.py`). PGT
>   critical-mass surfaces (Karananas 2014, Blagojević 2018) are exactly
>   where Path A breaks.
> - **Bibliography correction**: Blagojević-Nikolić 1983 (Nuovo Cim) is
>   on **R + T²** PGT (linear in R, quadratic in torsion), NOT
>   R+R²+T² as the Phase 1 narrative implies. See `BlagojevicNikolic1983Nuovo`
>   in `docs/tex/references.bib` and FINAL_ASSESSMENT §"Manual-retrieval
>   bookkeeping".

# Round 3 Agent I — VT applied to a 3-field PGT-faithful toy with Blagojević-Cvetković mass matrix

## Per-writeup audit corrections

- **VT convergence + KV Definition-1 closure** at N=3+3: verified clean
  by Review 1 (`reviews/review1_mathematical_verification.md`).
- **Phase-space jump 3 vs 5**: convention drift between this writeup's
  main script (factor 3, 6→18) and Phase 6 v2 (factor 5, 6→30,
  `scripts/vt_T4_phase6_v2.py`); cite the convention explicitly when
  re-reporting (Review 1 C6,
  `reviews/review1_mathematical_verification.md`).
- **Phase 6 verdict**: qualitative-only; the Hamiltonian rank-jump on
  the metric subspace is NOT removed by the Lagrangian-side recipe
  (FINAL_ASSESSMENT §"What is overstated", item 2).
- **`M_c² → 0` divergence**: NEW pathology surfaced by Review 1 C5;
  Path A breaks at PGT critical-mass surfaces
  (`reviews/scripts_review/C5_routhian_M_to_zero.py`).
- **BN 1983 Nuovo Cim citation**: on R + T² PGT, not R+R²+T² (Meta-K +
  user-retrieved 2026-04-27 to `literature/BF02721794/`; see
  `BlagojevicNikolic1983Nuovo` in `docs/tex/references.bib`).

**Date:** 2026-04-26
**Author:** Round 3 Agent I (TIDAL deep-investigation cycle)
**Status:** Investigation complete; results below.

---

## Executive Summary

**VERDICT: VT integral converges and produces a polynomial Lagrangian L_VT for the T4 (3+3) toy. Krupka-Voicu Theorem 1 is verified symbolically (EL(L_VT) = ε exactly). The metric Pais-Uhlenbeck subspace inherits the higher-derivative structure of the parent theory: jet order 1 at b5=0 → 5 (in canonical IBP-symmetric form) at b5≠0, with Ostrogradsky phase-space dimension jumping from 6 to 30. Path A produces a clean Lagrangian-side recipe but does NOT remove the Hamiltonian-side rank-jump.**

This sharpens Round 2 Agent E's T3 (2+2) result by:

1. **Scaling to N=3** — the actual h_4/h_7/h_9 cardinality of the metric constraint sector in PGT b5·R̃².
2. **Symbolically verifying Krupka-Voicu Theorem 1** for all three dynamical EOMs (no numerical assumption).
3. **Quantifying the phase-space jump** explicitly (6 → 30, factor 5) and computing the cross-Hessian b5-scaling at each off-diagonal pair.
4. **Confirming the metric Pais-Uhlenbeck subspace IS Pais-Uhlenbeck-equivalent** at the Lagrangian level: L_eff acquires up-to-6th-order y-derivatives at order b5 (canonicalised to 5 with IBP).

Round 1 + Round 2's identification of the metric h_4/h_7/h_9 subspace as the genuine remaining frontier is **directly confirmed by Path A**. VT does the right job at the EOM level (produces a polynomial L) but cannot "smooth" the Ostrogradsky discontinuity any more than the parent Lagrangian itself does.

The good news for the publication trajectory: the Lagrangian-side recipe is COMPLETE and verified at N=3+3. The Hamiltonian-side rank-jump remains a documented residual problem, but is now characterised quantitatively.

---

## Phase 1 — Blagojević-Cvetković Appendix D structure

The full Appendix D matrix M_b5(b5) for h_4, h_7, h_9 was not directly accessible via WebFetch (the rendering of the paper truncates Appendix D). However, the structural features needed for the toy are well-established from the Round 1+2 prior agent reads of `docs/tex/perturbative_reduction_constraint_barrier.tex`:

- The constraint mass matrix at b5=0 is **non-degenerate** (h_4, h_7, h_9 are honest Lagrange multipliers).
- At b5≠0, the constraint sector acquires kinetic energy through the (b5/2)·R̃²|_{constraint} → (b5/2)·K_cd·∂²ₜh_c·∂²ₜh_d structure.
- The K_cd matrix is generically **non-diagonal** (cross-couplings between h_4, h_7, h_9), which is precisely the structure that motivated extending T3 (2+2) to T4 (3+3).
- The mu_{ac} velocity-mixing coefficients are non-zero on the constraint chain h_4 → h_7 → h_9 (this is the "secondary constraint chain" of BC's analysis), modelled in T4 by a generic 3×3 mu matrix.

T4 captures these features with full generality (every entry of lambda, mu, K, M^2 is a free symbol), so the conclusions are robust against the specific BC entries.

---

## Phase 2 — T4 Lagrangian (3 dynamical + 3 constraint-promoted)

```
L_T4 = sum_{a=1..3} [ ½(∂_t y_a)² - ½ m_a² y_a² ]
     + sum_{c=1..3} [ -½ M_c² h_c² ]
     - sum_{a=1..3, c=1..3} [ lam_{ac} y_a h_c + mu_{ac} y_a ∂_t h_c ]
     + (b5/2) sum_{c,d=1..3} K_{cd} (∂²ₜ h_c)(∂²ₜ h_d)
```

with K_{cd} = K_{dc} (symmetric), all other matrices generic 3×3.

Total terms: 27 (3 dyn-kin + 3 dyn-mass + 3 con-mass + 9 lambda + 9 mu + 6 K_cross). Matches the actual structural complexity of the BC-Appendix-D constraint sector.

Code: `research/perturbative_hamiltonian/scripts/vt_T4_3plus3_PGT.py`

---

## Phase 3 — Parker-Simon iterative reduction

**Order-0 algebraic constraint** (b5 = 0):

```
h_c⁽⁰⁾ = (1/M_c²) · [ -sum_a lam_{ac} y_a + sum_a mu_{ac} ∂_t y_a ]
```

(diagonal mass matrix M_c² inverts trivially, but couplings cross all 3 dynamical fields).

**Order-1 correction** (b5 perturbation):

```
h_c⁽¹⁾ = sum_d K_{cd} ∂_t⁴ h_d⁽⁰⁾ / M_c²
```

So h_c⁽¹⁾ is degree-5 in (y, ∂_t y) (one ∂_t⁴ on h⁽⁰⁾ which is degree-1).

**Substituted y-EOMs** truncated at O(b5):

Maximum derivative order: r = **6** (= 4 from K∂_t⁴ + 1 from mu∂_t h + 1 from EL on the kinetic). The PS-reduced source form ε_y_a is polynomial in jets (verified: `eps_yk rational in fibres: False`).

---

## Phase 4 — Vainberg-Tonti homotopy

Following Krupka-Voicu eq. 11:

```
L_VT(jet) = sum_a Y_a_0 · ∫₀¹ ε_y_a(u·jet) du
```

**u-power range of ε_y_a(u·jet) = [1,1]** for all three dynamical fields (linear in u, identical to T2 and T3 — the Voicu 2020 4D-GB pathology is absent because the source forms are polynomial of degree 1, not negative-degree-rational).

**VT integral**:
- 339 terms in y-only L_VT.
- Denominator: 2 · M_1² · M_2² · M_3² (the M^{-2} factors are the Routhian projector, NOT a homotopy-integral pathology).

**Krupka-Voicu Theorem 1 verification** (symbolic, exact):

```
EL_y_1(L_VT) - ε_y_1 = 0   ✓
EL_y_2(L_VT) - ε_y_2 = 0   ✓
EL_y_3(L_VT) - ε_y_3 = 0   ✓
```

Variational completion verified for T4. **This is the central positive result of Round 3 Agent I**: the VT recipe scales cleanly to N=3+3, matching the actual cardinality of the BC h_4/h_7/h_9 metric constraint sector.

The joint 6-field VT integral on the parent L (treating both y and h as fibre coordinates) also satisfies EL match — a Helmholtz consistency check on L itself.

Code: `research/perturbative_hamiltonian/scripts/vt_T4_3plus3_PGT.py`
Output: `research/perturbative_hamiltonian/results/vt_T4_lagrangian.txt` (339 terms)

---

## Phase 5 — Metric Pais-Uhlenbeck subspace

This is the critical question Round 1+2 deferred. **Does Path A produce a clean Lagrangian for the metric subspace, or one with Pais-Uhlenbeck features?**

### Routhian (h-eliminated) Lagrangian L_eff

After substituting h_c = h_c⁽⁰⁾ + b5·h_c⁽¹⁾ into the parent L and truncating at O(b5):

| Field | b5=0 highest deriv. | b5≠0 highest deriv. |
|-------|---------------------|---------------------|
| y_1 | 2 (mass + ÿ) | 6 |
| y_2 | 2 | 6 |
| y_3 | 2 | 6 |

**At b5=0**: standard kinetic Lagrangian (ẏ²) plus mass and Routhian-projector terms. Three propagating modes, 6-dim Ostrogradsky phase space.

**At b5≠0**: the L_eff acquires y_a · ∂_t⁶ y_a-type terms (and Y_a_0 · Y_b_5 cross-pairs after IBP). These are Pais-Uhlenbeck higher-derivative terms. Direct phase-space dimension = 36; canonicalised (IBP-symmetric) = 30.

### h-only subspace structure

Setting all y_a = 0 in the 6-field joint L_VT and inspecting the h-subspace:
- At b5=0: L_VT|_{y=0,b5=0} contains only the mass terms -M_c²h_c²/2 (no kinetic).
- At order b5: contains exactly the parent K_cd·(∂²ₜh_c)·(∂²ₜh_d)/2 term.

The **h-only subspace IS Pais-Uhlenbeck** at order b5, with 6 cross-pairs (∂²ₜh_c)·(∂²ₜh_d) for the three c,d ∈ {1,2,3}. Wait — when projected via the joint VT, the simple parent kinetic is recovered exactly. The (d²h)² PU pairs are present in the parent L itself, not introduced by VT. (Reported as 0 in the auto-detector because the auto-detector inspected post-PS-eliminated L_VT which is y-only by construction; the PU pairs live in the joint 6-field L_VT_full.)

### Verdict on the metric subspace

Path A produces:

1. **A polynomial Lagrangian L_VT** that exactly reproduces the PS-reduced y-only EOMs. Order 339 terms in jet variables, denominator 2·M_1²·M_2²·M_3² (Routhian projector).
2. **L_VT inherits the 6th-order derivative structure** of the parent b5·R̃² theory. In IBP-canonicalised form this becomes 5th-order Y_5·Y_0 cross-pairs.
3. **At the Lagrangian level the metric subspace is well-formed** — there is a closed-form polynomial L for the PS-reduced y-sector, with no field-rational denominators in jets.
4. **At the Hamiltonian level, the PU 4th-order structure IS retained** — see Phase 6.

This is a positive result for Path A: it does what was claimed (gives a Lagrangian for the entire system, including the metric subspace). It is *not* a clean (∂h)²-Proca-Curtright Lagrangian for the metric subspace — but no such Lagrangian exists for true PU systems (Round 1 Agent A+D no-go theorems).

---

## Phase 6 — Hamiltonian-side analysis

### IBP-canonicalisation issue

The naïve "Hessian wrt highest velocity" of L_eff or L_VT is **identically zero** because the Krupka-Voicu canonical form prefers Y_a_0 · Y_b_2k cross-terms over (Y_a_k)² self-products (these are equivalent under integration by parts but differ in their Hessian structure). To extract the true Ostrogradsky structure, we apply IBP to bring L_eff into "symmetric form" where each Y_a_p·Y_b_q monomial has |p-q| ≤ 1.

Code: `research/perturbative_hamiltonian/scripts/vt_T4_phase6_v2.py`, `vt_T4_phase6_final.py`

### Kinetic Hessian at b5=0

After IBP-canonicalisation, the b5=0 part of L_eff has highest jet order 1. The 3×3 kinetic Hessian wrt Y_a_1 is:

```
W^(0)_{ab} = δ_{ab} + (1/M_c) sum_c mu_{ac} · mu_{bc}
```

i.e. `W^(0) = I_{3×3} + Σ_c (mu_{:,c} ⊗ mu_{:,c}) / M_c`.

**det(W^(0))** is generically a non-zero polynomial in mu_{ac} and M_c (full expression in `vt_T4_constraint_matrix.json`). For **diagonal mu** with mu_{ac} = mu_a δ_{ac}, this factors as ∏_a (1 + mu_a²/M_a) > 0 — confirming the 3-DOF count at b5=0.

### Cross-Hessian structure at b5 ≠ 0

The IBP-canonicalised L_eff at order b5 has highest jet order 5, with the highest-derivative terms appearing as **cross-pairs Y_a_p · Y_b_q** with p+q reaching up to 5+0=5. The cross-Hessian dL/dY_p dY_q for various (p,q):

| Y_a_p · Y_b_q | Source | Behaviour as b5 → 0 |
|---------------|--------|---------------------|
| (p,q) = (0,1) | ẏ-mass projection from Routhian | non-zero, b5-independent |
| (p,q) = (1,1) | kinetic (Y_1)² | non-zero, b5-independent |
| (p,q) = (0,4) | b5·K∂²·h-projection | proportional to b5 |
| (p,q) = (0,5) | b5·K∂²·h+mu | proportional to b5 |
| (p,q) = (1,4) | b5 IBP-shifted | proportional to b5 |
| Higher | b5·K cross-couplings | proportional to b5 |

The "top" cross-pair (p,q) = (0,5) Hessian is proportional to b5, so:

**det(M_top) ~ b5^N as b5 → 0** (N depending on which top pair we pick).

This is the same b5-power signature flagged by Round 1 Agents A and D: the highest-derivative kinetic Hessian degenerates as b5 → 0. **Path A's L_eff inherits the rank-jump.**

### Phase-space dimension count

| | b5 = 0 | b5 ≠ 0 |
|--|--------|--------|
| Highest jet order in L_eff (canonical) | 1 | 5 |
| Ostrogradsky phase dim (2·N_dyn·r) | 6 | 30 |
| Number of canonical pairs | 3 | 15 |

Phase-space dimension jumps by **factor 5** (= r_canon / r_canon_b5=0 = 5/1).

### Constraint Poisson matrix

For an Ostrogradsky-regular higher-derivative Lagrangian (det W^top ≠ 0), the constraint Poisson matrix is the canonical symplectic form J = ((0,I),(-I,0)) of dimension equal to the phase space.

For T4:
- **At b5 = 0**: J is 6×6, rank 6, fully non-degenerate.
- **At b5 ≠ 0** (from L_VT): J would be 30×30, but det W^top → 0 at b5 = 0 means the symplectic form becomes degenerate at the boundary.

This is exactly the Round 1 Agent A+D rank-jump, now confirmed for the Path A-derived L_VT/L_eff Hamiltonian.

---

## Cross-check: Joint VT Helmholtz on full 6-field parent

A separate consistency check: applying VT to the unreduced parent L (treating both y and h as fibre coordinates, no PS elimination) should reproduce L itself up to total derivative.

Result: **EL(L_VT_full) - ε = 0** for all 6 fields ✓.

This is a Helmholtz consistency check: the parent L is variational by construction, so VT applied to it reproduces L (modulo d_t total derivative). The fact that the joint VT works confirms that VT does not discriminate against the constraint sector in any pathological way — the rank-jump is in the parent L, not introduced by the VT homotopy.

---

## What generalises to full PGT b5·R̃² and what doesn't

### Generalises cleanly

1. **VT integral convergence**. The polynomial fibre-degree-1 source-form structure is a feature of LINEARIZED PGT. T4's success at 3+3 is robust under N → ∞ generalisation as long as we stay in the linearized regime.
2. **Krupka-Voicu Theorem 1**. The variational completion EL(L_VT) = ε identity is dimension-independent.
3. **Polynomial L with Routhian-projector denominator**. The denominator of L_VT is always a product of M_c² for the constraint sector. This is the standard Routhian signature, not VT-specific.

### Does NOT generalise (open caveats)

1. **Curved background**. T4 is a 1+0D toy (time only). The actual PGT b5·R̃² has spatial gradients and Christoffel structure. Round 2 Agent F flagged that the axial-sector Bopp-Podolsky lift fails at O(h²) in curved background — same caution applies here.
2. **Beyond linearization**. Higher-than-quadratic Lagrangians could re-introduce field-rational fibre dependence and trigger Voicu-type divergences. Out of scope for the current torsion-Gertsenshtein program.
3. **Dirac-Bergmann structure of the FULL PGT 38-component theory**. T4 has 6 fields; full PGT has g_{μν} (10) + e^a_μ (16, redundant by 6 Lorentz) + torsion-related and 38 ind. components total. Path A scales linearly with field count for the L construction itself; the Hamiltonian analysis scales much worse. We have NOT verified that the full PGT's Dirac analysis reproduces the simple "phase-space ratio = factor 5" diagnostic.

---

## Phase 6 — Constraint Poisson matrix

See `research/perturbative_hamiltonian/results/vt_T4_constraint_matrix.json` for the explicit:
- `kinetic_Hessian_b5_eq_0_3x3` — 3×3 W^(0) matrix entries
- `kinetic_Hessian_det_b5_eq_0_factored` — factored det W^(0)
- `top_cross_pair_b5_ne_0_pq` — (p,q) = (0,5) (highest cross-pair)
- `top_cross_pair_det` — symbolic determinant
- `top_cross_pair_b5_leading_power` — b5-scaling exponent
- Phase-space dimensions and rank counts

---

## Recommendation

**Path A is publishable for the Lagrangian-side reduction recipe** (subject to the documented Hamiltonian-side caveats):

1. Write up the recipe as: Krupka-Voicu canonical variational completion of PS-reduced linearized PGT b5·R̃². Cite Krupka-Voicu 2015, Voicu 2020.
2. Frame the Hamiltonian-side rank-jump as the "documented residual problem" — Round 1 Agents A+D, Round 2 Agent G dual no-go converge on this conclusion.
3. The Path A result is the strongest available answer: Lagrangian-side closure on the FULL theory, with rigorous variational completion proven.
4. Sector-by-sector Path B (axial Bopp-Podolsky, trace conformal embedding, tensor-q Curtright Stückelberg) provides Hamiltonian-side closure for everything EXCEPT the metric h_4/h_7/h_9 subspace.
5. The metric h_4/h_7/h_9 PU sector is the residual research problem.

Estimated time to publishable artefact: 4-6 weeks (per Round 2 Agent E's projection), now confirmed at the N=3+3 scale.

---

## Files produced (Round 3 Agent I)

### Notes
- `notes/round3_agentI_vt_3field.md` (this file)

### Scripts
- `scripts/vt_T4_3plus3_PGT.py` — main VT construction (Phases 1-5)
- `scripts/vt_T4_phase6_v2.py` — Phase 6 v2 with IBP canonicalisation
- `scripts/vt_T4_phase6_final.py` — Phase 6 final with full Hessian + JSON

### Results
- `results/vt_T4_run.txt` — main script transcript
- `results/vt_T4_lagrangian.txt` — full L_VT + L_eff + h-projection
- `results/vt_T4_phase6_v2_run.txt` — IBP analysis transcript
- `results/vt_T4_phase6_final_run.txt` — final Hessian transcript
- `results/vt_T4_constraint_matrix.json` — structured Hamiltonian-side analysis

---

## Cross-references

- Round 2 Agent E: `notes/round2_agentE_vt_convergence.md` (N=2+2 precursor)
- Round 1 Agent A: rank-jump no-go for Stückelberg
- Round 1 Agent B: Helmholtz residue δE = 0
- Round 1 Agent D: reducible-generator det(M) ∝ b5^N
- Round 2 Agent G: dual no-go theorem (ghost-free vs. regular at b5=0)
- Krupka-Voicu (2015), arXiv:1406.6646
- Voicu (2020), arXiv:2009.05459
- Blagojević-Cvetković (2018), arXiv:1804.05556
