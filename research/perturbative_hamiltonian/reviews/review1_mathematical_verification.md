# Review 1: Mathematical Verification of the Perturbative-Hamiltonian Investigation

**Reviewer**: Critical mathematical audit (sympy-level verification)
**Date**: 2026-04-27
**Scope**: Re-execute and audit Round 2/3 sympy scripts; verify that headline claims survive critical scrutiny.

---

## Executive Summary

All 18 Round 2 / Round 3 sympy scripts re-execute cleanly. The core mathematical computations are correct as far as they go. **However, the Round 3 synthesis framing overstates several claims and understates several caveats.**

Headline findings:

1. **NEW PATHOLOGY (C5)**: Path A's L_VT contains explicit `1/M_c²` and `1/M_c⁴` poles (Routhian projector denominator). It diverges as any constraint mass `M_c² → 0`. Round 3 Agent I did not flag this; the synthesis does not warn that Path A is undefined on critical-mass surfaces (which are exactly where the most physically interesting PGT critical cases live).
2. **CONVENTION DRIFT (C6)**: The "phase-space jump 6 → 30, factor 5" cited by Round 3 synthesis uses a non-standard PS-reduced-jet count. The standard Ostrogradsky count gives factor 3 (matching `vt_T4_3plus3_PGT.py`'s own output of `6 → 18`). The two conventions are reported inconsistently across the Round 3 results.
3. **CROSS-VALIDATION OVERSTATED (C8)**: Agent F's axial-sector lift and Agent J's Curtright lift do NOT have parallel mathematical structure. Agent F has `det(M) ~ 1/b²` (rank-jump moved into aux block); Agent J has b5-independent `det(H_kin) = 1 − λ_a²` (genuine kinetic uniformization). The "structural parallel" claim in Round 3 synthesis is qualitatively true (both decouple a heavy field) but quantitatively misleading.
4. **VOICU SECOND CONDITION CHECKED (C4)**: Linearity in highest-order derivatives (a separate Voicu gate, missed by Agent E) does pass for the toys at O(b5) and O(b5²). This is a useful new verification but the gate becomes nontrivial for beyond-quadratic PGT.
5. **CURTRIGHT GAUGE INVARIANCE EXTENDS (C3)**: δF̊ = 0 verified in 4D for all 24 Curtright components. **But** rank uniformity and cyclic-constraint preservation are NOT verified in 4D; Agent J's verdict relied on the 1+1D toy.
6. **NO-GO THEOREM (C7)**: Agent G's "dual no-go" extends to 3-form, tensor-scalar, and derivative-coupled scalar auxiliaries. Local first-order auxiliary lifts are exhausted. Non-local lifts (Aoki-Mukohyama infinite tower) and field-redefinition approaches are NOT covered.
7. **HELMHOLTZ ROBUST (C2)**: Closure holds to k=16 jet orders. The non-trivial-ness of the check was demonstrated by injecting spurious non-variational terms — Helmholtz correctly detects them.
8. **1+1D EXTENSION OF T2 (C1)**: Spatial-gradient generalization of Agent E's 0+1D toy preserves convergence and KV closure. The toy → PGT spatial-gradient gap is closed at one cardinality.

This review complements:
- Review 2 (literature) — found author fabrications, missing theorems, and parity-mismatch in the Chatzistavrakidis-Ranjbar-Zekoč (CRZ) result.
- Review 3 (project relevance) — found jet-order incompatibility with TIDAL's JSON schema, denominator pathologies in L_VT, and that Path B excludes the metric subspace that matters for Gertsenshtein.

The combined verdict from all three reviews: the mathematical work is genuine, but Round 3 synthesis overclaims on three independent fronts (literature attribution, project applicability, internal mathematical framing). Significant edits are required before any of the three proposed publications goes out.

---

## Per-script audit (re-execution results)

### Round 2 — Agent E

| Script | Re-execution | Verdict |
|---|---|---|
| `vt_convergence_T2.py` | clean (sympy verifies `EL(L_VT) − ε = 0`) | **PASS — but limited**: 0+1D (no spatial gradient), Helmholtz only checked to `k = 6` |
| `vt_convergence_T3.py` | clean (`EL(L_VT)−ε = 0`, denom `2·M2_1²·M2_2²`) | **PASS — but limited**: 0+1D, denom rational in mass parameters |

### Round 2 — Agent F

| Script | Re-execution | Verdict |
|---|---|---|
| `bopp_podolsky_axial.py` | clean; sympy reports `det M = m_A^4 / b^2`, `det(b·M) = b² m_A^4` (vanishes at `b=0`) | **PASS at linear-flat**, but the constraint-Poisson **det vanishes at b=0** (rank-jump at the parameter limit). Author reads this as "rank-jump pushed into the auxiliary block where it is harmless" — see C8 below for why this is structurally identical to the metric blocker, not parallel to it. |

### Round 2 — Agent G

| Script | Re-execution | Verdict |
|---|---|---|
| `line3_2form_auxiliary.py` | clean | **PASS** — confirms 2-form aux has `det ∝ (1 − b5·m_K²)`, fails for the same Agent A/D structural reason |
| `line3b_2form_IBP_constraint_check.py` | clean | **PASS** — establishes the dual no-go (smooth limit ↔ ghost) for 2-form auxiliaries only |
| `line4_born_oppenheimer.py` | clean | **PASS** — qualitative; no novel sympy claim |
| `line6_dual_formulation.py` | clean | **PASS** — qualitative |

### Round 3 — Agent H

| Script | Re-execution | Verdict |
|---|---|---|
| `recipe1_preflight_q_projection.py` | clean | **PASS** — projection schema gives k-degree ≤ 2, but this is verified on a hand-coded sympy schema, NOT against the actual xAct `b5·R̃²` decomposition |
| `recipe1_explicit_q_substitution.py` | clean | **PASS** — same caveat as above |
| `recipe1_all_eps_RR_contractions.py` | clean | **PASS** — confirms d_5, d_13 nonzero with k-deg=2; d_1, d_9 vanish |
| `recipe1_check_all_R_components.py` | clean | **PASS** — 144/256 nonzero R̃ components |
| `recipe1_parity_even_RR.py` | clean | **PASS** — α_4 gives k-deg 2 |

### Round 3 — Agent I

| Script | Re-execution | Verdict |
|---|---|---|
| `vt_T4_3plus3_PGT.py` | clean; reports phase-space jump 6 → **18** (factor 3) | **PASS but with discrepancy**: Round 3 synthesis claims **6 → 30 (factor 5)**. Script itself reports **6 → 18**. The 6 → 30 figure comes from a *separate* Phase 6 script (`vt_T4_phase6_v2.py`) using `2·N·r` with `N=3, r=5`. The two counts disagree because they use different conventions: the parent script's "velocity-jet dim" stops at `(N, jet+1) = (3, 6) = 18` for jet ≤ 5, while phase6_v2's `2·N·r` is `2·3·5 = 30`. **The synthesis cites the higher number without flagging this convention dependence.** |
| `vt_T4_phase6_v2.py` | clean; 6 → 30 | **PASS** under the `2·N·r` convention. |

### Round 3 — Agent J

| Script | Re-execution | Verdict |
|---|---|---|
| `curtright_stueckelberg_q.py` | clean; reports `det(H_kin) = 1 − λ_a²`, b5-independent | **PASS in 1+1D toy.** The verdict explicitly notes: parity-odd extension needed (open research), `metric_h4_h7_h9_solved: False`. |

---

## Critical concerns C1–C8 (independent verification)

### C1 — 1+1D extension of T2 (toy → PGT spatial-gradient gap)

**Script**: `reviews/scripts_review/C1_T2_1plus1D_extension.py`

Built a 1+1D toy with `L = ½(∂_t φ)² − ½(∂_x φ)² − ½m²φ² − λφh − μφ∂_t h − ν φ∂_x h − ½M²h² + ½b5(∂_t² h)²`.

**Result**:
- PS-reduced ε_φ contains both `P02` (= ∂_x² φ) and `P40` (= ∂_t⁴ φ).
- VT integrand u-degree min = 1, max = 1 → polynomial → integral converges.
- `EL(L_VT) − ε_φ_PS = 0` exactly, including the cross-coupled spatial term.

**Verdict**: **PASS.** The 0+1D → 1+1D extension preserves both convergence and KV closure. The convergence argument carries over to spatial gradients.

**Caveat**: This is still a toy-level verification with one dynamical and one constraint field. The actual PGT `b5·R̃²` has 38 fields and cross-sector couplings; PS-reduced source forms may contain mixed terms `(∂_x^j ∂_t^k φ)·(∂_x^l ∂_t^m χ)` with cross-jet indices that this toy does not exhibit.

### C2 — Helmholtz residue at higher orders (k = 8, 10, 12)

**Script**: `reviews/scripts_review/C2_helmholtz_higher_orders.py`

The PS-reduced ε_phi at O(b5) has highest derivative `P6 = ∂_t^6 φ`. The Helmholtz operator for an order-r Lagrangian has support on jet orders up to 2r. With r = 5 (truncated PS-reduced order), Helmholtz needs support to k = 10. Round 1 Agent B checked only k = 0..6.

**Result**:
- Built the source form, ran the VT homotopy, computed `EL(L_VT) − ε_PS` with `MAX = KMAX + 4 = 16` jet orders.
- Residue is identically **zero** at this truncation.
- Robustness check: deliberately added a non-variational term `b5·P[7]·P[3]`. Helmholtz correctly **detects** the violation: residue ≠ 0.
- Note: adding `b5·P[8]` alone (linear) gave residue 0 — that term IS variational (it's a total derivative of a linear-in-P[8] Lagrangian). This is a legitimate counter-example demonstrating that the Helmholtz check is non-trivial.

**Verdict**: **PASS.** Helmholtz holds at all jet orders the PS-reduced source form can probe.

### C3 — Curtright Stückelberg in higher D

**Script**: `reviews/scripts_review/C3_curtright_2plus1D.py`

Agent J's verification was 1+1D where (2,1)-Young has 1 component after antisymmetry. Full PGT has 16 components in 4D after cyclic-constraint enforcement.

**Result**:
- Built plane-wave amplitudes `T_{μν|ρ}` (D²(D−1)/2 components, antisymmetric in [μν]).
- Auxiliaries: `h_{μν}` (symmetric), `b_{μν}` (antisymmetric), `a_μ` (vector).
- Stückelberg shifts: `s_{μν}` (symm), `β_{μν}` (antisym), `α_μ` (vector).
- **D=3 (2+1D)**: All 9 components of `δF̊` vanish identically. ✓
- **D=4 (3+1D)**: All 24 components of `δF̊` vanish identically. ✓

**Verdict**: **PASS — gauge invariance only.** δF̊ = 0 holds in 4D.

**Critical caveats this script does NOT address**:
1. **Cyclic constraint q_{[μν|ρ]} = 0**: in 4D this removes 8 of the 24 components, leaving 16. The script does NOT verify that the Stückelberg shifts preserve the cyclic constraint. This is non-trivial; the paper claims it holds, but it's not what δF̊ = 0 alone proves.
2. **Rank uniformity in 4D**: Agent J verified `det(H_kin) = 1 − λ_a²` only in the 1+1D toy. The 4D version requires a 36×36 (or after-cyclic-quotient 28×28) Hessian — Agent J flagged this as "sympy-tractable but heavy" but did not run it.
3. **Parity-odd extension**: Review 2 confirmed the published paper handles only parity-even free fields. δF̊ = 0 holds independent of parity, but the actual `b5·R̃²` parity-odd contractions (38 ε·DT·DT terms per Recipe 1) need their own gauge-invariance check, which has NOT been performed.

So the 4D δF̊ = 0 result is genuine but only a small piece of what Agent J's "applies with caveats" verdict requires.

### C4 — Voicu linearity-in-highest-derivative gate

**Script**: `reviews/scripts_review/C4_voicu_linearity.py`

Voicu 2020's *second* failure mode (separate from "homogeneity ≤ −1"): the source form must be linear in highest-order derivatives. Review 2 noted Agent E only addressed homogeneity.

**Result**:
- O(b5): highest derivative is `P6`. `d²ε/dP6² = 0` ✓
- O(b5²): highest derivative is `P10`. `d²ε/dP10² = 0` ✓

**Verdict**: **PASS.** The PS-reduced source form is linear in highest derivatives at both orders.

**Mechanism**: PS reduction by iterative substitution of polynomial-in-jets `h0 + b5·h1 + b5²·h2 + ...` produces linear-in-each-jet expressions because:
- Each PS step is `h_{n+1} = (b5/M²)·D_t^4 h_n`, which is a linear differential operator.
- Linear operators acting on linear-in-jet expressions stay linear-in-jet.

**Caveat**: This is a property of the **toy at PS order** and survives because the parent Lagrangian is quadratic. For the full PGT `b5·R̃²` Lagrangian — which is bilinear in linearised R̃ — the source form is also linear in jets. So the linearity gate passes for the actual PGT linearised problem. **For beyond-quadratic PGT (cubic and higher self-interactions), this is not automatic** — Voicu's gate becomes a real obstruction.

### C5 — Routhian denominator at M_c² → 0 — **NEW PATHOLOGY FOUND**

**Script**: `reviews/scripts_review/C5_routhian_M_to_zero.py`

The L_VT for T4 has denominator `2·M_1²·M_2²·M_3²`. What happens at `M_c² → 0`?

**Result** (with simplified 2+2 toy):
- L_VT denominator: `2·M_1⁴·M_2⁴`
- 16 terms with `1/M_1²` poles, 6 terms with `1/M_1⁴` poles
- `limit(L_VT, M_1 → 0) = ∞·sign(...)` — **L_VT genuinely diverges as a constraint mass goes to zero**.
- The `O(b5)` terms have HIGHER negative powers of M_c than the `O(b5⁰)` terms.

**Verdict**: **NEW PATHOLOGY — Round 3 Agent I did not flag this.**

For PGT, constraint masses M_c are functions of the Lagrangian coupling constants. At certain critical surfaces (notably ghost-free PGT critical cases — see Karananas 2014, Blagojević 2018), some M_c² → 0. **Path A's L_VT is undefined on those surfaces.** This is independent of the b5 → 0 critical surface and adds a *second* class of singular limits where Path A fails.

For the *Gertsenshtein* application specifically: the actual `torsion_gertsenshtein.json` parent Lagrangian's M_c values for h_4/h_7/h_9 should be checked. If any equals zero or scales with the conversion-relevant parameter, Path A fails there too.

**Recommendation**: Round 3 synthesis should add this caveat to the Path A verdict.

### C6 — Phase 6 simpler version (rank-jump quantitative)

**Script**: `reviews/scripts_review/C6_simpler_phase6.py`

Built N=1+1 toy. Computed the kinetic Hessian explicitly in (∂_t Y, ∂_t z = ∂_t² H) basis: `diag(1, b5)`.

**Result**:
- At b5=0: rank 1, primary constraint π_z = 0.
- At b5≠0: rank 2, no primary constraint.
- Phase-space dimension: 2 → 6 (factor 3).

**Verdict**: **PASS — rank-jump robust at N=1+1.** The qualitative claim from Agent I generalizes structurally.

**Discrepancy uncovered in Round 3 synthesis claim**: The "phase-space jump 6 → 30 (factor 5)" comes from `vt_T4_phase6_v2.py`'s convention of `2·N·r` with `r = 5` (PS-reduced jet order). The parent script `vt_T4_3plus3_PGT.py` reports `6 → 18` (factor 3) using `2·N·r_parent` with `r_parent = 2` (parent Lagrangian order). For N=1, my C6 calc gives factor 3 (matching parent-jet count). **The factor-5 number cited in Round 3 synthesis uses the PS-reduced jet count, which is a non-standard counting; the physically correct Ostrogradsky count is factor 3.** The synthesis should clarify which convention it uses, or report both.

### C7 — Dual no-go theorem generality

**Script**: `reviews/scripts_review/C7_3form_aux_counterexample.py`

Tried 3-form, tensor-scalar, derivative-coupled scalar auxiliaries. None evades Agent G's dichotomy:
- Algebraic auxiliary → rank-jump (Agent A/D no-go applies)
- Derivative-kinetic auxiliary → ghost OR wrong-sign OR nonlocal

**Verdict**: **PASS** for the class of *local first-order auxiliary lifts*. Agent G's no-go genuinely extends beyond 2-forms.

**Caveats** (not addressed by C7):
- Non-local lifts (Pauli-Villars heavy-mode towers) — not a counter-example to the no-go because they introduce infinitely many DOF; out-of-scope for "perturbative Hamiltonian reduction".
- Higher-order Stückelberg (jet-2 auxiliaries) — not enumerated.
- Field-redefinition-only approaches (no new fields, only nonlinear redefinitions) — would be a genuine alternative; not yet investigated.

So Agent G's no-go is generic for first-order lifts, but the literature has examples of *non-perturbative* recipes (e.g., Aoki-Mukohyama 2020's infinite tower) that Agent G's argument doesn't cover. This was already flagged by Review 2 in the AM cross-validation discussion.

### C8 — Cross-validation independence (F vs J)

**Script**: `reviews/scripts_review/C8_F_vs_J_structural_diff.py`

Compared Agent F's `axial_constraint_matrix.json` to Agent J's `curtright_stueckelberg_verdict.json`.

**Findings**:

| Property | Agent F (axial) | Agent J (Curtright-q) |
|---|---|---|
| Object | Constraint Poisson M | Kinetic Hessian H |
| det at b/b5 → 0 | DIVERGES (`1/b²`) | UNCHANGED (`1 − λ_a²`) |
| Rank-jump? | YES (`b·M` has vanishing det at `b=0`) | NO |
| Mechanism | Aux mass → ∞, rank-jump moved into aux block | q mass → ∞, kinetic uniform |

**Verdict**: **Agent F's "rank uniformity" claim applies ONLY to the A-sector block.** The auxiliary block carries the same Lyakhovich rank-jump as the parent theory; it has been re-located, not removed. Agent J's lift, by contrast, has b5-independent kinetic Hessian.

The Round 3 synthesis claim that the two lifts are "structurally parallel" via heavy-mode decoupling is *qualitatively* true but **quantitatively misleading**. Agent F's mechanism is a *Lyakhovich rank-shuffling*; Agent J's is a *Stückelberg kinetic uniformization*. The latter is what the metric h₄/h₇/h₉ sector actually needs (and Agent J explicitly flags that this lift does NOT apply to the metric sector).

This corroborates Review 2's "cross-validation inflation" finding from the mathematical side. The cross-validation between the three Round 3 paths (E, F, J) is *one verification (Agent J's) plus two qualitative consistency checks (E and F) of different things*, not three independent verifications of the same thing.

---

## Headline-claim survival

| Claim | Survives critical review? | Caveat |
|---|---|---|
| VT integral converges for PS-reduced linearised PGT | **YES** at toy level (T2, T3, T4); 1+1D extended (C1) | Only verified for parent Lagrangian quadratic; full 38-component PGT not done |
| Helmholtz residue δE = 0 | **YES** to k=16 (C2) | — |
| KV closure (`EL(L_VT) = ε`) | **YES** | Should be re-titled "KV Definition 1 closure", not "Theorem 1" (per Review 2) |
| Voicu linearity-in-highest-deriv gate | **YES** (C4) — was missed by Round 2 | Linear at toy level; nonlinear PGT extension not verified |
| Routhian denominator well-behaved | **NO** (C5) — diverges as M_c² → 0 | NEW PATHOLOGY: Path A is undefined on critical-mass surfaces |
| Curtright Stückelberg gauge invariance δF̊ = 0 in 1+1D | **YES** in 4D too (C3) | Cyclic constraint not checked; rank uniformity in 4D not verified |
| Curtright rank uniformity `det(H_kin) = 1 − λ_a²` | **YES** at 4×4 1+1D toy | 4D 36×36 Hessian NOT computed; parity-odd extension is open |
| Recipe 1 PASS (q-projection standard-kinetic) | **YES** at projection-schema level | Sympy schema NOT cross-checked against actual xAct `b5·R̃²` decomposition |
| Three convergent no-go theorems for metric subspace | **YES** (C6 confirms rank-jump robust; C7 confirms Agent G generality for local first-order auxiliaries) | Does NOT exclude non-local, higher-order Stückelberg, or field-redefinition approaches |
| Phase-space jump factor 5 (6 → 30) for T4 | **MISLEADING** — direct Ostrogradsky gives factor 3 (6 → 18); factor 5 uses PS-reduced jet count | Round 3 synthesis should specify the convention |
| F-J cross-validation as "structurally parallel" | **OVERSTATED** (C8) — F is rank-shuffling, J is uniformization | These are different mechanisms with different mathematical content |

---

## Recommendations

### Tier 1 — must fix before any publication

1. **Drop "factor 5" without context**: in Round 3 synthesis and any downstream artefact, either report the convention (`2·N·r_PS` with PS-reduced jet `r=5`) or use the standard Ostrogradsky count (factor 3 for the toy, generally `2 + 2·r_parent`).
2. **Add the M_c² → 0 caveat to Path A**: the L_VT denominator divergence is a genuine new pathology not flagged by Round 3 Agent I. Either (a) prove the denominator factor cancels in the actual PGT projection, or (b) add a precondition `M_c² ≠ 0 for all c` to Path A's applicability statement.
3. **Re-title KV verifications**: per Review 2, "Theorem 1" should be "Definition 1 closure check" or "Helmholtz consistency".
4. **Add Voicu linearity gate as a separate preflight** (Review 2 + C4): currently bundled into "VT convergence"; should be its own bullet.
5. **Tone down F-J "cross-validation"**: per C8, these are not parallel; they are two different mechanisms applied to two different sectors. The Round 3 synthesis "Cross-validation matrix" should be revised to reflect: F gives Lagrangian closure on axial via rank-shuffling; J gives kinetic uniformization on tensor-q.

### Tier 2 — should investigate before claiming "publishable now"

6. **Verify Recipe 1 against the actual xAct `b5·R̃²` decomposition**: Agent H's preflight verified a hand-coded sympy projection schema. The actual TIDAL pipeline runs in xAct. A direct check on the explicit_terms_tex.txt or analogous output would close this gap.
7. **Compute the 4D Curtright kinetic Hessian** (36×36 or after-cyclic-quotient 28×28) to verify rank uniformity in 4D. Agent J flagged this as "sympy-tractable but heavy"; it should be done before claiming the construction "applies".
8. **Check whether the cyclic constraint `q_{[μν|ρ]} = 0` is preserved** by the Stückelberg shifts in 4D. Gauge invariance of δF̊ does not imply preservation of the cyclic constraint; this needs an independent verification.
9. **Cross-reference the actual M_c values** for h_4/h_7/h_9 from `torsion_gertsenshtein.json` against the C5 divergence. If `M_c² > 0` everywhere on the parameter space of interest, the C5 caveat may be harmless in practice.

### Tier 3 — could improve the academic publications

10. **Beyond-quadratic PGT extension**: Voicu's linearity gate is automatic for quadratic Lagrangians. For cubic PGT, it becomes a real test. If the publication claims theory-agnosticism, this needs a worked example.
11. **Field-redefinition-only approaches**: not enumerated by any of the no-go theorems. A literature survey for `b5(∂²h)² → polynomial-Lagrangian` field redefinitions would close a real gap.
12. **Verify Agent F's M_aux divergence is "pushed into harmless block"**: Agent F's verdict says the rank-jump is in the aux block, which is "decoupled at the level of the A-sector". The full Hamiltonian observable depends on both blocks; "decoupling at the A-sector" needs a precise statement of what observables remain regular and which diverge.

### Tier 4 — do NOT pursue

13. Do not extend TIDAL's JSON schema to support jet-5 operators (per Review 3 Tier 4).
14. Do not rebrand Path A as "we have a working Hamiltonian for PGT b5·R̃²" — at minimum, the PU structure persists and the M_c² → 0 caveat is a new condition.
15. Do not cite the F-J cross-validation as evidence that the *metric sector* is solvable — it is explicitly NOT, and the cross-validation matrix in the synthesis needs revision.

---

## Cross-cutting issues (mathematical side)

### Issue A — Toy faithfulness boundary

The toys (T2, T3, T4, axial Bopp-Podolsky, 1+1D Curtright) faithfully capture the b5(∂²h)² Pais-Uhlenbeck barrier. They do NOT capture:
- Spatial-gradient cross-couplings to constraint fields (C1 is a partial extension; mixed-jet terms not exhibited).
- Cyclic Young constraints (C3 verifies δF̊ = 0 only).
- Cross-sector couplings (Review 3 already flagged).
- Constraint-mass critical surfaces (C5 is the new find).

The boundary between "what the toys verify" and "what the full PGT requires" is wider than Round 3 synthesis presents.

### Issue B — Verification convention drift

Two distinct phase-space dimension counts are reported in Round 3 results: `6 → 18` (parent script, parent-jet basis) vs `6 → 30` (Phase 6 v2, PS-reduced-jet basis). The synthesis cites the larger figure without flagging the convention. This is sloppy: it makes the rank-jump look more dramatic than the standard Ostrogradsky count gives.

### Issue C — Recipe 1 verification depth

Agent H's sympy verification is a *projection schema* — a hand-coded model of how `b5·R̃²` would project onto the q-irreducible sector. It is **not** a direct verification on the xAct output. This is a genuine "trust me" gap. The verdict "Recipe 1 PASS" is conditional on the schema correctly reflecting the actual xAct pipeline.

### Issue D — Helmholtz check robustness

C2's robustness check (deliberately adding spurious terms) demonstrated that the Helmholtz test correctly *detects* non-variational source forms. This is a useful check that Round 1 Agent B did not perform. It strengthens confidence in the Round 1 result (which only checked that the closure holds, not that the closure check is non-trivial).

---

## Bottom line for the user

**Path A (Vainberg-Tonti)**: closure verified at toy level, Helmholtz holds to k=16, Voicu linearity gate (newly checked) PASSES. **NEW concern**: L_VT diverges at constraint-mass critical surfaces (C5). The "factor 5 phase-space jump" claim should be tightened to "factor 3 by direct Ostrogradsky" or have the convention specified.

**Path B-axial (Bopp-Podolsky)**: verified at linear-flat. Constraint Poisson divergence at b=0 is real (det ~ 1/b²). The "rank-uniform" claim is sector-restricted to the A-block; the full lift inherits the parent rank-jump in the auxiliary block.

**Path B-tensor-q (Curtright)**: gauge invariance δF̊ = 0 verified in 4D (newly extended). Rank uniformity, cyclic-constraint preservation, and parity-odd extension all remain unverified in 4D.

**Three convergent no-go theorems**: robust against the alternative auxiliary types tested (3-form, tensor-scalar, derivative-coupled scalar). Local first-order auxiliary lifts are exhausted. Non-local, higher-order, and field-redefinition approaches not excluded.

**Cross-validation framing**: weakened. The three "independent paths" are not independent verifications of the same fact; they are different recipes for different sectors with different mathematical structures. Review 2's "cross-validation inflation" finding is corroborated on the math side.

The core conclusion of Round 3 — *the metric h₄/h₇/h₉ subspace is the genuine residual frontier* — survives critical review. What does NOT fully survive is the framing that "Path A + Path B-axial + Path B-tensor-q" together provide a complete sectoral solution; in fact, only Path B-tensor-q has a genuinely uniform construction, only at toy level, and only modulo open caveats.

**The investigation has produced real mathematical content. The framing in the Round 3 synthesis overstates the strength of the cross-validation and understates the new caveats (M_c² → 0 pathology, factor-3-vs-5 convention drift, Agent F's rank-jump-in-aux-block). These should be tightened before any publication.**
