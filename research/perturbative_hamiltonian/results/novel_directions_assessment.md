> **⚠ SUPERSEDED (2026-04-27)**: this ranked assessment was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. Several headline claims about
> the "literature support" for Lines 1-6 are overstated. The verified
> picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This
> file is retained for historical record only. **Do not propagate its
> specific claims without checking against the audit.**
>
> Specific corrections:
> - **Path B-trace** "established literature" framing is **WRONG**
>   (Meta-K K3, `meta_reviews/meta_review_K_literature_claims.md`):
>   Barker 2024 identifies trace torsion as Yang-Mills GAUGE FIELD,
>   not Goldstone, and explicitly excludes parity-odd terms.
> - **Path B-tensor-q** is **CONDITIONAL** on a parity-odd extension
>   that does not exist in the published literature (Meta-K K5,
>   `meta_reviews/meta_review_K_literature_claims.md`): CRZ 2024
>   handles only parity-even free fields and the m → 0 Goldstone
>   limit, opposite to TIDAL's b5 → 0 infinite-mass limit.

# Round 2 Agent G — Novel Directions: Ranked Assessment

## Final ranking (most-to-least promising)

| Rank | Line | Verdict | Concrete outcome | Action item |
|------|------|---------|------------------|-------------|
| 1 ★★★ | 6 (dual / Curtright Stückelberg) | **PROMISING + NOVEL** | Identified arXiv:2411.16928 Chatzistavrakidis–Ranjbar–Zekoč Stückelberg lift for massive (2,1)-Young tensors. Covers tensor-torsion q-irreducible *if* Curtright projection is standard-kinetic. | Run Recipe 1 preflight (sympy/xAct expansion of b5·R̃² projected onto q-irreducible). |
| 2 ★★ | 3a/b (2-form auxiliary) | NEW NO-GO THEOREM | Constructed sympy-verified dual to Agent D's theorem: regular Hessian + smooth limit ⇒ ghost. Extends Agent D's no-go to cover the regular-Hessian loophole. | Document the dual no-go in the constraint-barrier doc. |
| 3 ★ | 4 (Born–Oppenheimer) | RE-DERIVATION OF PATH A | BO adiabatic reduction reproduces Path A's Vainberg-Tonti recipe under a different name. Useful as auxiliary justification but not new content. | Cite BO as alternative justification for Path A's algebraic substitution. |
| 4 ◯ | 1 (BV-BFV homological) | BLOCKED | arXiv:2309.07327 explicitly assumes regularity; no irregular extension in 2024–2026 literature. | Mark as future-research direction in docs; not actionable. |
| 5 ◯ | 5 (WKB matching) | NO NOVEL LEAD | Closest twin (f(T) Stückelberg in cosmology) is background-dependent, not transplantable to TIDAL's flat Minkowski theory-agnostic setup. | Cite Hou-Cai-Li 2023 as twin in docs; no constructive recipe. |
| 6 ◯ | 2 (Kontsevich rank-jump) | UNEXPLORED | Multi-year research programme; no published "deformation across rank change" framework. | Mark as multi-year programme; not actionable. |

Legend: ★★★ = promising + novel; ★★ = new theorem; ★ = derivative-but-useful; ◯ = no-go or not-applicable.

## Sectoral status of the constraint-promotion barrier

Refined picture after Round 2 Agent G:

| Sector | Constraint-promotion type | Recipe | Source |
|--------|--------------------------|--------|--------|
| Axial torsion | Higher-derivative `(∂·A)²` | Bopp-Podolsky single-aux | Round 1 Agent C |
| Trace torsion | Goldstone of broken Weyl | Conformal embedding | Barker et al. 2024 (arXiv:2406.12826) |
| Tensor torsion (q-irreducible) | **Standard-kinetic Proca-Curtright (NEW)** | **Chatzistavrakidis–Ranjbar–Zekoč Stückelberg** | **Round 2 Agent G + arXiv:2411.16928** |
| Metric h₄,₇,₉ trace promotion | Pais–Uhlenbeck (genuine 4th-order) | **Still blocked** | Round 1 Agents A+D no-go + Round 2 Agent G dual no-go |

**Net result**: 3 of 4 constraint-promoted sectors now have plausible
constructive recipes; only the metric Pais–Uhlenbeck case remains.

## Concrete next steps (in priority order)

1. **Recipe 1 preflight (highest priority)**: verify whether b5·R̃² projected
   onto the tensor-torsion irreducible q^a_{bc} produces standard-kinetic
   `(∂q)²` structure or higher-derivative `(∂²q)²` structure. Estimated
   effort: 2–4 hours of xAct work. Without this preflight, the Round 2
   Agent G result is suggestive but not yet demonstrated.

2. **Recipe 2 (if Recipe 1 confirms standard-kinetic)**: construct the
   explicit Curtright Stückelberg Lagrangian using the auxiliary triple
   (graviton, Kalb-Ramond, vector) from arXiv:2411.16928 §3, and verify
   the smooth b5 → 0 limit at the level of the canonical Hamiltonian
   (Dirac analysis, 1–2 days).

3. **Document the dual no-go theorem** (Round 2 Agent G) in
   `docs/tex/perturbative_reduction_constraint_barrier.tex` as a
   strengthening of the Agent D result.

4. **Cite the unified Path A-via-BO interpretation** in the same doc as
   alternative justification for the algebraic-substitution recipe.

5. **Speculative future direction**: parity-odd Hinterbichler–Saravani
   extension for the metric h₄,₇,₉ case. This is a multi-month research
   programme.

## Files produced

Scripts (in `research/perturbative_hamiltonian/scripts/`):
- `line3_2form_auxiliary.py` — 2-form auxiliary lift sympy investigation
- `line3b_2form_IBP_constraint_check.py` — V2-IBP constraint-rank check
  (key result: det(Hessian) = -1, ghost via diagonalisation)
- `line4_born_oppenheimer.py` — Born–Oppenheimer adiabatic reduction
- `line6_dual_formulation.py` — Curtright dual analysis + sectoral
  re-classification

Result transcripts (in `research/perturbative_hamiltonian/results/`):
- `line3_2form_run.txt`
- `line3b_2form_IBP_run.txt`
- `line4_BO_run.txt`
- `line6_dual_run.txt`

Notes (in `research/perturbative_hamiltonian/notes/`):
- `round2_agentG_novel_directions.md` — comprehensive writeup
- (this file) — `results/novel_directions_assessment.md`

## New arXiv references identified (not in Round 1)

- **arXiv:2411.16928** (Chatzistavrakidis–Ranjbar–Zekoč 2024) — *Tensor
  global symmetries and the Stueckelberg mechanism for tensor fields*,
  JHEP 05 (2025) 218. **★ critical reference** for tensor-torsion
  Curtright Stückelberg.
- **arXiv:2512.18017** (Paci–Solodukhin 2025) — *Auxiliary-Field
  Formalism for Higher-Derivative Boundary CFTs*. Higher-derivative
  auxiliary-field reduction; relevant context.
- **arXiv:2512.08972** (Bouzari Nezhad 2026) — *DHOST in metric-affine
  gravity*, JCAP. Inverse problem to TIDAL; useful comparison.
- **arXiv:2304.08360** (Martini–Paci–Sauro 2023) — *Covariant
  spin-parity decomposition of Torsion*. Covariant framework for
  q-irreducible projection.
- **arXiv:2502.17979** — *Avoiding singularities with propagating
  torsion*. Recent propagating-torsion construction.
- **arXiv:2601.22007** (Aashish–Saif 2026) — *Stückelberg approach for
  antisymmetric tensor field in Lorentz violation*. Structurally
  analogous to TIDAL's Stückelberg restoration of regular constraint
  matrix.
- **arXiv:2305.10298** (Hou–Cai–Li 2023) — *f(T) Stückelberg*.
  Background-dependent twin.
- **arXiv:2402.16938** — Mixed symmetry gauge fields in flat
  background. Generic mixed-symm Stückelberg literature.
