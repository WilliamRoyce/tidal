> **⚠ SUPERSEDED (2026-04-27)**: this per-agent writeup was audited by
> Reviews 1-3 and Meta-Reviews K/L/M/N. The sympy execution underlying
> the report is verified clean (Review 1), but several of its framing
> claims are overstated. The verified picture is in
> `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`. This file
> is retained for historical record only. **Do not propagate its specific
> claims without checking against the audit.**
>
> Specific issues: KV "Theorem 1" does not exist (paper has Definition 1
> only — Meta-K K2); only ONE of Voicu 2020's two failure modes is checked
> (Review 2 §2: linearity-in-highest-derivatives gate omitted); NEW
> pathology `L_VT diverges as M_c² → 0` was not flagged at the time
> (Review 1 C5).

# Round 2 Agent E — VT Integral Convergence Test

**Date**: 2026-04-26
**Scope**: Make-or-break gate for Path A (Krupka-Voicu canonical
variational completion) of TIDAL's perturbative-Hamiltonian reduction
for `b5·R̃²`-class theories.

## Executive Summary

**VERDICT: VAINBERG-TONTI INTEGRAL CONVERGES.** Path A is *unobstructed*
on integral-convergence grounds for TIDAL's PGT `b5·R̃²` linearised
sector.

Concretely:

1. **T2 (Agent B's two-field falsifier)** — VT integrand
   `y^σ · ε_σ(u·y)` is a polynomial in `u` of degree exactly 1
   (`min u-power = +1`). The integral over `[0,1]` is therefore a
   finite, regular polynomial in jet variables. The Euler-Lagrange of
   the resulting `L_VT` reproduces the original PS-reduced EOMs `E^(1)_φ`,
   `E^(1)_χ` *exactly* (Krupka-Voicu canonical-completion theorem 1).
2. **T3 (a more PGT-faithful 2+2 toy with `b5·K_ab·ḧ_a·ḧ_b` cross-couplings)**
   — same structure: polynomial integrand in `u` of degree 1, finite
   `L_VT`, EL match verified.
3. **No analog of Voicu 2020's 4D-Gauss-Bonnet pathology exists** in
   the TIDAL setup. The 4D-GB divergence is driven by the *negative*
   homogeneity degree of the truncated Lanczos-Lovelock tensor `A_μν`
   under `g_μν → t·g_μν` (`A_μν` scales as `t^{-1}`, giving a `(D-4)^{-1}`
   pole in the regulated integral). TIDAL's PS-reduced source forms are
   *polynomial* in fibre variables of degree ≥ 1, so the integrand
   carries `min u-power = +1`, not `−1`.

The remaining `M^{-2}` factors in `L_VT` (T3) come from inverting the
algebraic constraint mass matrix during PS reduction — they are *not*
VT homotopy pathologies but the standard Routhian-projector signature.
The `M^2 → 0` limit is the massless-constraint singular surface, a
distinct physical limit from `b5 → 0`.

**Recommendation**: proceed with Path A (full PGT `b5·R̃²` Lagrangian-side
canonical variational completion). Estimated 4–6 weeks to publishable
recipe.

---

## Phase 1 — VT integral on T2

### Setup recap (Agent B)

Lagrangian:
```
L_T2 = ½(∂_t φ)² + ½(∂_t χ)²
     − ½ m_φ² φ² − ½ m_χ² χ²
     − λ φ h − λ_χ χ ∂_t h
     − ½ h² + ½ b5 (∂²_t h)²
```

Order-0 algebraic constraint: `h⁽⁰⁾ = −λ φ + λ_χ ∂_t χ`.
Iterative correction: `h = h⁽⁰⁾ + b5 · ∂⁴_t h⁽⁰⁾ + O(b5²)`.

PS-reduced order-1 EOMs (from `vt_convergence_T2.py`):
```
ε_φ(jet) =  λ² Φ₀ − m_φ² Φ₀ − Φ₂ − λ λ_χ Χ₁
          + b5 [ λ² Φ₄ − λ λ_χ Χ₅ ]

ε_χ(jet) =  −m_χ² Χ₀ − Χ₂ − λ_χ² Χ₂ + λ λ_χ Φ₁
          + b5 [ −λ_χ² Χ₆ + λ λ_χ Φ₅ ]
```
(where `Φ_k = ∂_t^k φ`, `Χ_k = ∂_t^k χ`.)

### VT integrand

Krupka-Voicu eq. 11:
```
L_VT(jet) = Φ₀ · ∫₀¹ ε_φ(u·jet) du + Χ₀ · ∫₀¹ ε_χ(u·jet) du
```

Under the homothety `χ_u : jet → u·jet` (uniform scaling of every fibre
coordinate by `u`, parameters left alone), each ε is *linear* in fibre
variables (because the Lagrangian is *quadratic* in fields → EL is
linear). So
```
ε_φ(u·jet) = u · ε_φ(jet)
ε_χ(u·jet) = u · ε_χ(jet)
```
and the integrand `Φ₀ · u · ε_φ(jet)` has `min u-power = +1`. The integral
is trivially finite:
```
∫₀¹ u · ε_φ(jet) du = ½ ε_φ(jet)
```

### VT Lagrangian (T2 explicit)

```
L_VT = ½ λ² Φ₀² − ½ m_φ² Φ₀² − ½ Φ₀ Φ₂
     + b5 · ½ λ² Φ₀ Φ₄
     − ½ λ λ_χ (Φ₀ Χ₁ + Χ₀ Φ̇)        — derivative-mixing
     + b5 · (−½ λ λ_χ) (Φ₀ Χ₅ − Χ₀ Φ₅)
     − ½ m_χ² Χ₀² − ½ (1 + λ_χ²) Χ₀ Χ₂
     + b5 · (−½ λ_χ²) Χ₀ Χ₆
```

### Consistency check (Krupka-Voicu Theorem 1)

```
EL_φ(L_VT) − ε_φ = 0   ✓
EL_χ(L_VT) − ε_χ = 0   ✓
```

Variational completion verified.

### Comparison with the on-shell Routhian

Substituting `h = h_sol(b5)` directly into `L_T2` gives
`L_eff(jet)` (Agent B's "direct" Lagrangian). It contains *rational*
dependence on jet variables (a `(Χ₁ λ_χ − Φ₀ λ)^{-1}` denominator from
the on-shell substitution). However, its *Euler-Lagrange equations*
match `L_VT`'s EL exactly:
```
EL_φ(L_eff) − EL_φ(L_VT) = 0
EL_χ(L_eff) − EL_χ(L_VT) = 0
```
So `L_eff` and `L_VT` agree on-shell (modulo total derivatives + EOM
contractions); `L_VT` is the *polynomial canonical representative*.
This is exactly the Krupka-Voicu point: VT picks out the canonical
member of the equivalence class.

### Parameter-limit pathology check

Each of `b5 → 0`, `m_φ → 0`, `m_χ → 0`, `λ → 0`, `λ_χ → 0` gives a
finite polynomial `L_VT`. No divergence in any limit. The only
denominator factor is `2` (from `∫₀¹ u du = ½`).

---

## Phase 2 — Comparison with Voicu 2020 4D-GB

### What goes wrong in 4D-GB

From Voicu 2020 §3.2 / §4 (extracted via `pdftotext` from
arXiv:2009.05459):

The classical VT Lagrangian (their eq. 12, equivalent to Krupka-Voicu eq. 11)
applied to the truncated Gauss-Bonnet field equations:
```
L = − g_μν ∫₀¹ t^{D/2} √(−g) [ t^{−2} M_P² G^{μν} + t^{−1} Λ_0 g^{μν}
                              + 2t^{−3} α A^{μν} + W^{μν}/(D−4) ] dt
```
(their eq. 24). Under uniform metric rescaling `g_μν → t·g_μν`:
- `G^μν` is degree 0 (Einstein tensor is conformally minimal in this sense)
- `A^μν` is degree `−3`
- `W^μν` is degree `−3`

When raised to upper-upper and densitised by `√(−g)` (degree `D/2`),
we get an integrand with `t^{D/2 − 3} = t^{(D−6)/2}` for the GB pieces.
Combined with the leading factor of `t` from `g_μν` in the formula,
the GB sector has integrand `t^{(D−4)/2 − 1}`. The integral
`∫₀¹ t^{(D−4)/2 − 1} dt = 2/(D−4)` exists only for `D > 4` and
**diverges as `(D−4)^{-1}` at D = 4**. This is a logarithmic
divergence in the original integral, regulated as a pole when
analytically continued in D.

**Driver**: negative homogeneity of the truncated Lanczos-Lovelock
tensor in the metric.

**Voicu's modified construction (their Theorem 1, eq. 16)**: replace
the lower endpoint `t=0` with a generic `t=a` such that
`lim_{t→a} t · E_A(t·y) = 0`. This rescues *some* divergent cases but
*not* 4D-GB: the regulated Lagrangian still has a `(D−4)^{-1}` pole
because the GB invariant is intrinsically topological in 4D.

### Why TIDAL's PS-reduced ε_σ avoids this

Three structural differences:

1. **Polynomial source forms.** TIDAL's PS-reduced EOMs `ε_σ(jet)` are
   polynomials of degree 1 in fibre variables (because the parent
   Lagrangian is *quadratic* — linearised regime). Under
   `jet → u·jet` they scale as `u^{+1}`, not `u^{-1}`. The integrand
   carries `u^{+1}`, the integral is `½ y^σ ε_σ`.

2. **No conformal-rescaling trade-off.** Voicu's homothety acts on the
   metric `g_μν` — a *physical* field. TIDAL's homothety acts on
   linearised perturbation jets `(δg, δA, δh, ...)` around a fixed
   background. There is no `√(−g)` density factor that brings in
   `t^{D/2}`: the Lagrangian is already a scalar density in the *full*
   non-linear pipeline, but at linear order around a fixed background
   the relevant volume measure is fixed.

3. **No topological obstruction at the linear order.** Gauss-Bonnet is
   a topological term in 4D *only after non-linear assembly*. The
   linearised `b5·R̃²` is *not* a topological derivative — its
   linearised EOMs are honest 4th-order kinetic operators acting on
   the linearised torsion fields, not boundary-only contributions.

The combination means TIDAL's PS-reduced source form is a **degree-1
polynomial in fibres**, automatically convergent under the standard VT
homotopy.

### Could the analog ever appear in TIDAL?

The only routes that could re-introduce a Voicu-style divergence are:

- **Beyond linearisation**: if one extends the `b5·R̃²` analysis to the
  full non-linear PGT theory, the source form may become
  *non-polynomial* in the jets (`R̃²` carries rational factors of
  `det e_μ^a` after vielbein decomposition). At quadratic order this
  is fine; at cubic and higher, non-polynomial structure is possible.
  *Linear-order Path A does not encounter this.*
- **External truncations**: if one *drops* terms from the linearised
  EOM that would otherwise restore variationality (Voicu's "truncation
  of Lanczos-Lovelock"). TIDAL never does this — the Wolfram pipeline
  derives EOMs from the *full* perturbative Lagrangian and only *then*
  PS-reduces.

Both of these failure modes are absent from TIDAL's current
linear-order PGT pipeline. They become potential issues *only* if
one tries to extend Path A to non-linear PGT — which is well beyond
the make-or-break gate for the current torsion-Gertsenshtein program.

---

## Phase 3 — T3: PGT-faithful toy with constraint-sector cross-couplings

### Motivation

T2 has only one constraint field. PGT `b5·R̃²` linearised has *three*
constraint-promoted scalars (`h_4, h_7, h_9`) with nontrivial
cross-couplings via the kinetic `(∂² h_a)(∂² h_b) K_{ab}` matrix from
`R̃²`. If T2's clean polynomial structure were a feature of having only
one constraint field, T3 would expose it.

### Setup

```
L_T3 = Σ_a [ ½(∂_t q_a)² − ½ m_a² q_a² ]   (dynamical sector, N fields)
     + Σ_a [ −½ M_a² h_a² ]                 (constraint masses)
     − Σ_{a,b} [ λ_{ab} q_a h_b + μ_{ab} q_a ∂_t h_b ]  (mixing)
     + (b5/2) Σ_{a,b} K_{ab} (∂²_t h_a)(∂²_t h_b)       (PGT cross-coupling)
```

with `K_{ab} = K_{ba}` symmetric, all matrices generic.

Implemented for `N = 2` (analytically tractable, exercises full matrix
structure). Order-0:
```
h_a^(0) = (M²)^{-1}_{ab} ( −λ_{cb} q_c + μ_{cb} ∂_t q_c )
```

### Result

Both `ε_q1(u·jet)` and `ε_q2(u·jet)` have **u-power range [1, 1]** (pure
linear in u, identical to T2). The integral gives a finite polynomial
in jet variables. EL match verified for both fields.

`L_VT` denominator: `2 · M²_1² · M²_2²`. The `M^{-2}` factors are the
algebraic-constraint Routhian projector, not a VT homotopy artefact.
`L_VT` is finite (and polynomial in jets) for any non-zero `M²_a`.

### Why the matrix structure doesn't introduce divergence

The PS reduction is a **finite-step Neumann series** in `b5`:
```
h = (M²)^{-1} (Λ + μ ∂_t) q   [order 0]
  + b5 · (M²)^{-1} K ∂_t^4 [(M²)^{-1} (Λ + μ ∂_t) q]   [order 1]
  + O(b5²)
```
Each application of `(M²)^{-1}` adds an algebraic factor; each
application of `K ∂_t^4` adds a homogeneous degree-1 jet operator. The
result is a polynomial in jets times rational functions of the
*parameter* `M²`. The VT homotopy acts only on the *jets*, not on
parameters, so the integrand is `u^{+1} · polynomial(jet, params)` —
manifestly convergent.

This generalises immediately to `N = 3` (the actual PGT `h_4, h_7, h_9`
case): the matrix algebra is identical, just larger.

---

## Comparison table

| Feature | TIDAL PGT b5·R̃² (PS-reduced) | Voicu 2020 4D-GB (truncated) |
|---|---|---|
| Field-space homogeneity of source form | degree +1 (linear) | degree −1 (rational in metric) |
| Integrand u-power | `u^{+1}` | `u^{(D−4)/2 − 1}`, divergent at D=4 |
| Driver of divergence | none | non-linear topological term + truncation |
| `(D−4)^{-1}` pole | absent | manifest |
| EL of L_VT reproduces source | yes (Helmholtz residue δE = 0, Agent B) | no — only after adding canonical correction |
| Routhian / VT relationship | `L_eff` rational, `L_VT` polynomial; EL agree | not applicable (full PDE system) |

---

## Recommendation: proceed with Path A on the full PGT linearised theory

The make-or-break gate is cleared. Concrete next steps (4–6 week
timeline):

1. **Generalise T3 to N = 3** (3 dynamical TT modes + 3 constraint-promoted
   `h_4, h_7, h_9`) using the actual PGT mass-matrix structure from
   Blagojević-Cvetković (2018) Appendix D. ~1 day of sympy.
2. **Apply VT to the full linearised PGT `b5·R̃²` source form** in
   Wolfram (read from existing `tidal_pgt.json`). Construct `L_VT`
   symbolically. ~1 week.
3. **Verify EL match**: compute Euler-Lagrange of `L_VT` and confirm it
   reproduces `tidal_pgt.json`'s `equations[]` array. ~3 days
   (canonical Helmholtz check; should be automatic given Agent B's δE=0).
4. **Document and benchmark**: write up in `docs/tex/perturbative_reduction_vt.tex`,
   compare runtime against the current Hamiltonian-modal pipeline.
   ~3 days.
5. **Publishable artefact**: a short note titled "Canonical variational
   completion of perturbatively-reduced higher-derivative Poincaré
   gauge theory" with explicit `L_VT` for `b5·R̃²`, the Krupka-Voicu
   theorem application, and a direct comparison with the failed
   Stückelberg attempts (Round 1 Agents A, D). ~2 weeks.

### Risk register

- **Background-field dependence**: TIDAL operates around a fixed
  Minkowski background. A curved background introduces position-dependent
  coefficient operators in the source form. The VT homotopy still
  acts only on fibres, not on `x`, so this should preserve convergence,
  but worth checking with one curved test (`sphere_kg`-class) before
  committing.
- **Higher-than-quadratic Lagrangians**: this entire convergence
  argument assumes linear source forms. Cubic-and-higher PGT
  Lagrangians could re-introduce field-rational structure. Out of
  scope for the current campaign but worth flagging.
- **Equivalence with Path B**: For the axial sector, Agent C's
  Bopp-Podolsky construction is *also* viable. Path A and Path B may
  produce different but equivalent reduced Lagrangians (related by a
  field redefinition). A comparison would be a clean cross-check.

---

## Files produced

- `research/perturbative_hamiltonian/scripts/vt_convergence_T2.py` — sympy
  driver for T2.
- `research/perturbative_hamiltonian/scripts/vt_convergence_T3.py` — sympy
  driver for T3 (N=2+2 PGT-faithful toy).
- `research/perturbative_hamiltonian/results/vt_T2_run.txt` — full T2
  run output with explicit `L_VT`.
- `research/perturbative_hamiltonian/results/vt_T3_run.txt` — full T3
  run output.

## References

- Krupka, Voicu (2015), *Canonical variational completion of differential
  equations*, J. Math. Phys. **56**, 043507; arXiv:1406.6646.
  Eq. 11 (VT Lagrangian), Eq. 12 (canonical correction).
- Voicu, N. et al. (2020), *Canonical variational completion and 4D
  Gauss-Bonnet*, Eur. Phys. J. Plus; arXiv:2009.05459. Eq. 16 (extended
  VT with arbitrary endpoint), Eq. 24 (4D-GB integrand and divergence).
- Round 1 Agent B (2026-04-26), Helmholtz residue computation: residues
  vanish identically for T2 → confirmed Path A preflight passes on
  Helmholtz grounds. (Reproduced inside `vt_convergence_T2.py`.)

---

## Per-writeup audit corrections (appended 2026-04-27)

The following corrections supersede framing claims in the body above. They
do **not** invalidate the sympy execution, which Review 1 re-ran cleanly.

- **"Krupka-Voicu Theorem 1 verified" is wrong.** The Krupka-Voicu paper
  arXiv:1406.6646 has *Definition 1* (canonical variational completion),
  not a Theorem 1. What was verified is the tautology that an already-
  variational source form admits a canonical completion. See
  `meta_reviews/meta_review_K_literature_claims.md` §K2 (Krupka-Voicu
  closure verification).
- **VT integrand convergence is sympy-correct but only one of two
  failure modes is checked.** Voicu 2020 (arXiv:2009.05459) §4 +
  Appendix A list TWO independent gates: (i) homogeneity / integral
  convergence; (ii) linearity-in-highest-derivatives. Agent E checked
  only (i). See `reviews/review2_literature_interpretation.md` §2 for
  the missing gate. Voicu's linearity gate must be added as a separate
  preflight.
- **NEW pathology not flagged at the time: `L_VT` diverges as
  `M_c² → 0`.** Path A's L_VT carries explicit `1/M_c²` and `1/M_c⁴`
  poles. PGT critical-mass surfaces (Karananas 2014, Blagojević 2018)
  are exactly where Path A breaks. See
  `reviews/review1_mathematical_verification.md` C5
  (`scripts_review/C5_routhian_M_to_zero.py`) — sympy-verified.
- **Path A produces a Lagrangian, not a Hamiltonian for the phase
  space.** Round 3 Agent I's Phase 6 confirms L_VT inherits Pais-Uhlenbeck
  structure for the metric subspace; the Legendre transform retains the
  rank-jump. See `notes/FINAL_ASSESSMENT.md` "What is overstated or
  wrong" §2.

For the consolidated verdict see
`research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`.
