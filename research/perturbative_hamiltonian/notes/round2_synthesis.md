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

# Round 2 Synthesis (2026-04-26)

Three deep theoretical investigations completed, building on Round 1's
four-agent literature + structural survey. Round 2's agents wrote actual
sympy code (verified, transcript-saved) and produced concrete computational
results.

## Headline results

1. **Path A (Vainberg-Tonti direct on full PGT b5·R̃²) is fully unobstructed.**
   Both preflight tests pass:
   - Helmholtz residue δE = 0 (Round 1 Agent B)
   - VT integral converges to a finite polynomial (Round 2 Agent E)
   The Voicu 2020 4D-Gauss-Bonnet pathology requires *negative* fibre
   homogeneity; TIDAL's PS-reduced ε_σ are degree-1 polynomials in jets
   (because the parent Lagrangian is quadratic — linearised regime).
   The integrand is `u^(+1)·polynomial(jet, params)`, manifestly convergent
   over [0,1]. **Verified on T2 toy and T3 (PGT-faithful 2+2) toy with
   sympy.**

2. **Path B-axial (Bopp-Podolsky single-auxiliary lift) is verified with
   one substantive correction.** The auxiliary is a SCALAR φ, not a
   vector B_μ as Agent C originally proposed. The Lagrangian dual of
   `(∂·A)²` is a scalar. Agent F's corrected Lagrangian:
   ```
   L_aux = -¼ F_{μν}F^{μν} - ½m_A² A_μ A^μ + φ (∂·A) - (1/(2b)) φ²
   ```
   with `b ≡ b5·κ_axial² = (9/4)·b5`. Constraint Poisson matrix is
   **rigorously block-diagonal** with `det(M_A) = m_A⁴` (b-independent)
   and `det(M_aux) = 1/b²` (decoupled aux block diverges as b→0, but
   the divergence is an infinite mass for φ — standard decoupling-limit
   behaviour).
   **Cross-validated against Aoki-Mukohyama 2020 arXiv:2009.11739** —
   their dynamical spin-0⁻ scalar (their varphi in Eqs. 497-526) IS
   this auxiliary at linearised order.

3. **Sectoral reclassification — major correction to Round 1's framing.**
   Round 1 conflated two structurally distinct constraint-promotion
   subcases under "tensor sector":
   - **Subcase A — metric h₄, h₇, h₉ promotion**: genuine Pais-Uhlenbeck
     4th-order from `b5·R̃²` acting on metric perturbations. Round 1
     Agents A and D's no-go theorems cover this case correctly.
   - **Subcase B — tensor-torsion q-irreducible promotion**: standard-
     kinetic Proca-Curtright field from `b5·R̃²`'s `(∂T)²` projection
     on the tensor torsion sector. **Not actually Pais-Uhlenbeck.**
     Accessible via published Stückelberg construction
     (Chatzistavrakidis-Ranjbar-Zekoč 2024, arXiv:2411.16928, JHEP 05
     (2025) 218).
   The actual TIDAL constraint-promotion barrier identified in
   `docs/tex/perturbative_reduction_constraint_barrier.tex` Eq.(eq:pr-cb-Mb5)
   is Subcase A — the metric Pais-Uhlenbeck — which **remains genuinely
   blocked**.

4. **NEW no-go theorem (Round 2 Agent G, dual to Agent D)**: no
   first-order auxiliary lift of `b5·(q̈)²` can be **simultaneously**
   ghost-free at b5≠0 AND regular at b5=0. Either ghost-free + rank-jump,
   or smooth limit + Pais-Uhlenbeck ghost. Sympy-verified for the
   2-form auxiliary candidate. This **strengthens Agent D's theorem** by
   covering the regular-Hessian loophole: even constructions that pass
   the Hessian regularity check can still fail the rank-uniformity test.

5. **Born-Oppenheimer ≡ Path A** (Agent G). The BO adiabatic separation
   of fast/slow degrees of freedom produces the same Lagrangian as
   Path A's Vainberg-Tonti algebraic substitution (`q[φ] = -λφ/M` is
   just the order-0 EOM). Useful as alternative justification but not
   new constructive content.

## Refined sectoral picture

| Sector | Type | Recipe | Status |
|--------|------|--------|--------|
| **Axial torsion `A_μ`** | Higher-derivative `(∂·A)²` | Bopp-Podolsky scalar-auxiliary lift `+φ(∂·A) − φ²/(2b)` | ✅ verified (Agent F) |
| **Trace torsion `S_μ`** | Goldstone-Weyl | Conformal embedding (Barker et al. 2024) | ✅ established |
| **Tensor torsion `q_{μνρ}`** | Standard-kinetic Proca-Curtright | Chatzistavrakidis-Ranjbar-Zekoč Stückelberg (3 auxiliaries) | 🟡 pending preflight |
| **Metric h₄, h₇, h₉** | Pais-Uhlenbeck (true 4th-order) | None — genuinely blocked | ❌ Agent A+D no-go + Agent G dual no-go |

PLUS the orthogonal:

| Path | Approach | Status |
|------|----------|--------|
| **Path A — Direct VT on full PGT** | Krupka-Voicu canonical variational completion | ✅ Both preflights passed (Agent B + Agent E) |

**The Path A result MAY subsume the sector-by-sector Path B if the metric
Pais-Uhlenbeck case is excluded** — VT would directly produce a
perturbative Lagrangian for the entire b5·R̃² PGT minus the metric h₄,₇,₉
contribution. The metric Pais-Uhlenbeck blocker remains.

## Critical caveats

1. **Recipe 1 preflight** (Agent G, ~2-4 hours xAct work): Verify that
   the b5·R̃² → q-irreducible projection produces only `(∂q)²`
   (standard-kinetic), not `(∂²q)²` (Pais-Uhlenbeck). If the latter,
   Subcase B collapses into Subcase A and the Curtright Stückelberg
   recipe doesn't apply.

2. **Path A non-linear extension**: VT integral convergence relies on
   polynomial source forms. Quadratic linearised PGT is fine, but
   beyond-quadratic could re-introduce field-rational structure. Out
   of scope for current torsion-Gertsenshtein program but worth flagging.

3. **Curved-background extension** (Agent F): the axial-sector
   Bopp-Podolsky lift fails at O(h²) because Nieh-Yan torsion-curvature
   cross-terms generate `h_μν · ∂²A · (∂h)` couplings not absorbed by
   the φ(∂·A) trick. TIDAL's flat-Minkowski regime is unaffected, but
   cosmological applications would need extension work.

4. **The metric Pais-Uhlenbeck blocker (Subcase A) is the genuine
   remaining frontier.** Three independent investigations (Round 1
   Agent A, Round 1 Agent D, Round 2 Agent G dual no-go) converge on
   the same conclusion: no Stückelberg-type lift bridges the b5=0
   discontinuity for true Pais-Uhlenbeck higher-derivative metric
   couplings. Path A handles this via VT directly **at the EOM
   level**, but the resulting "Lagrangian" may have non-standard
   features (worth checking explicitly for the metric h₄,₇,₉
   subspace specifically).

## Best-case publication trajectory

If Recipe 1 preflight confirms standard-kinetic structure for the
tensor-torsion-q sector:

```
Path A (Vainberg-Tonti)        — full theory, modulo metric h₄,₇,₉
+ Path B-axial (Bopp-Podolsky) — explicit closed-form for A_μ sector
+ Path B-trace (conformal)     — Barker et al. 2024 closed-form for S_μ
+ Path B-tensor-q (Curtright)  — Chatzistavrakidis-Ranjbar-Zekoč Stückelberg
                                 closed-form for q_μνρ sector
```

Yields a publishable Lagrangian-side reduction recipe for generic
`(α₁, α₂, α₃, b5)` PGT b5·R̃² **modulo the metric h₄,₇,₉ Pais-Uhlenbeck
gap**. The metric-sector gap remains as the documented residual research
problem — a clean, well-defined open question rather than a vague obstruction.

## Worst-case outcome

Even if every Path B sector verification fails, Path A alone produces a
publishable result for the FULL theory (subject to Path A's caveats —
the resulting Lagrangian may have unusual features for the metric
sector specifically, which would need separate investigation).

## Files produced in Round 2

### Notes
- `notes/round2_agentE_vt_convergence.md` — VT integral convergence analysis
- `notes/round2_agentF_axial_verification.md` — Bopp-Podolsky axial verification
- `notes/round2_agentG_novel_directions.md` — Six lines of investigation
- `notes/round2_synthesis.md` (this file)

### Scripts (sympy, all executable)
- `scripts/vt_convergence_T2.py` (456 lines) — VT integral on Agent B's T2 toy
- `scripts/vt_convergence_T3.py` (373 lines) — VT integral on PGT-faithful 2+2 toy
- `scripts/bopp_podolsky_axial.py` — full axial-sector verification incl. Hamiltonian + constraint matrix
- `scripts/line3_2form_auxiliary.py` — 2-form auxiliary attempt (failed)
- `scripts/line3b_2form_IBP_constraint_check.py` — dual no-go theorem verification
- `scripts/line4_born_oppenheimer.py` — BO ≡ Path A demonstration
- `scripts/line6_dual_formulation.py` — tensor-torsion dual analysis

### Results
- `results/vt_T2_run.txt`, `results/vt_T3_run.txt` — VT sympy output
- `results/axial_constraint_matrix.json` — explicit Poisson matrix decomposition
- `results/line3_2form_run.txt`, `results/line3b_2form_IBP_run.txt` — 2-form lift transcripts
- `results/line4_BO_run.txt`, `results/line6_dual_run.txt` — BO + dual transcripts
- `results/novel_directions_assessment.md` — Agent G's ranked-list assessment

## Round 3 priorities (next investigation)

1. **Recipe 1 preflight**: Verify the b5·R̃² → q-irreducible projection
   structure. xAct work; ~2-4 hours. Make-or-break for Path B-tensor-q.
2. **Apply Path A (VT) explicitly to a 3-field PGT toy** matching the
   actual mass-matrix from Blagojević-Cvetković 2018 Appendix D.
   Verify the resulting Lagrangian is well-formed for the
   constraint-promoted h₄,₇,₉ subspace (or document precisely how it
   fails there).
3. **Attempt to apply the Curtright Stückelberg recipe (Agent G's
   arXiv:2411.16928) to a tensor-torsion toy** to see whether the
   construction produces a viable PGT-applicable Lagrangian.
4. **Investigate metric Pais-Uhlenbeck (Subcase A) extreme cases** —
   is there a regime (e.g. specific gauge fixings, particular limits)
   where the PU 4th-order can be tamed despite the no-go theorems?
   Or document the no-go more rigorously by constructing a complete
   proof from the three convergent results.

## Cross-cutting observations

- **Compounding redundancy is a feature, not a bug.** Path A and Path
  B-axial give equivalent results for the axial sector via different
  routes — useful cross-validation. The fact that Agent F's corrected
  scalar-auxiliary matches Aoki-Mukohyama's spin-0⁻ mode is independent
  external validation.
- **The 25-year history of named-but-unsolved status (Round 1 Finding
  1) is consistent with where we landed**: the metric Pais-Uhlenbeck
  blocker IS the long-standing barrier; the tensor-torsion-q and
  axial sectors were not the actual obstruction — they are tractable
  with techniques that postdate or were missed by the prior literature.
  Our genuine contribution is the **clarified sectoral decomposition
  + Path A proof-of-viability** more than any fundamentally new
  technique for the metric blocker.
