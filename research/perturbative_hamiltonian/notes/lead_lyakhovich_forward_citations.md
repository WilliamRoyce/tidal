# Lead: Lyakhovich-Sharapov forward citations 2022-2026

**Date**: 2026-04-27
**Investigator**: Phase 2.3 forward-citation audit
**Scope**: Find papers in 2022-2026 that cite arXiv:2106.09355 (Abakumova-Lyakhovich)
or arXiv:2102.10579 (Lyakhovich) and assess whether they extend the recipe in a
way that bridges TIDAL's `b5=0` critical surface in PGT `b5·R̃²`.

## Search queries used

1. WebSearch: `"2106.09355" Abakumova Lyakhovich citations 2023 2024 2025`
2. WebSearch: `"2102.10579" Lyakhovich Stückelberg general gauge symmetry citations 2024 2025`
3. Inspire-HEP REST API: `refersto:recid:1869064` (cites 2106.09355) and
   `refersto:recid:1847795` (cites 2102.10579) and the OR-combination — confirmed
   record numbers via `q=arxiv:<id>&fields=control_number`.
4. WebSearch: `"Abakumova" Stueckelberg 2025 OR 2026 Poincare gauge torsion higher derivative`
5. WebSearch: `"reducible Stueckelberg" "rank jump" OR "constraint promotion" OR "irregular constraint" 2024 2025`
6. Direct grep on local arXiv TeX sources after download.

The Inspire-HEP query is the authoritative listing — the WebSearch results were
used to triangulate that no major paper is missing from Inspire-HEP and that no
2025/2026 PGT-Stückelberg paper has surfaced under different keywords.

---

## Section 1 — Forward-citation list with relevance triage

The Inspire-HEP combined query
(https://inspirehep.net/api/literature?q=refersto%3Arecid%3A1869064%20OR%20refersto%3Arecid%3A1847795)
returned **9 total hits**.  After removing
(a) arXiv:2106.09355 itself (it is in the result set because it cites 2102.10579
internally) and
(b) arXiv:2008.02112 (Pandey 2022) which is a 2020 BRST-on-Riemann-manifolds
preprint with only a tangential bibliographic citation,
the genuine forward-citation set is **7 papers**.

| # | arXiv | Year | First authors | Title | Self-cite? | Higher-deriv? | Constraint-rank-jump? | Triage |
|---|-------|------|---------------|-------|------------|---------------|----------------------|--------|
| F1 | 2206.07891 | 2022 | Abakumova, Frolovsky, Herbig, Lyakhovich | Gauge symmetry of linearised Nordström gravity and the dual spin two field theory | YES | NO (linearised Nordström) | NO | duality application |
| F2 | 2207.10634 | 2023 | Guleryuz | (Super)universal attractors and de Sitter vacua in string landscape | NO | NO (nilpotent superfield/string) | NO | tangential — string-landscape attractors, not constraints |
| F3 | 2303.02616 | 2023 | Abakumova, Lyakhovich | Dualisation of free fields | YES | NO (free spin-1, spin-2) | NO | duality generalization |
| F4 | "Bull.Lebedev Phys.Inst. 50(3):108" | 2023 | Abakumova, Lyakhovich | Dual Formulation for the Massless Spin Two Theory | YES | NO | NO | proceedings note, subsumed by F3 |
| F5 | 2402.12437 | 2024 | Casadio, Chataignier | Relaxation of first-class constraints and the quantization of gauge theories | NO | NO (relativistic particle, EM, gravity) | NO (deliberately *avoids* rank changes) | conceptual/historical Fock-Stueckelberg |
| F6 | 2405.13706 | 2024 | Delplanque, Skvortsov | Symmetric vs. chiral approaches to massive fields with spin | NO | NO (Singh-Hagen / Zinoviev / chiral HS) | NO (second-class → first-class is via well-known Zinoviev recipe) | massive higher-spin auxiliaries |
| F7 | "Imperial Coll. London thesis" (Erickson) | 2022 | Erickson | Localizations with noncompact transverse spaces and covert symmetry breaking | NO | NO (localization in SUSY field theory) | NO | unrelated context |

**Total non-self-citation papers**: 3 (F2, F5, F6). F7 is a thesis tangential
to gauge theory and not an arXiv preprint; F2 is a string-landscape paper that
cites Lyakhovich's machinery only as background.

**Direct-read candidates** (Task 2 filter — applies recipe to higher-derivative
or to parameter-dependent constraint rank): **none of F2/F5/F6 satisfies the
filter on its face**. Going through them anyway because (a) the audit demands
direct reads not abstract scans and (b) Round 3 Agent J already direct-read
2411.16928 (Curtright-Stückelberg) which is even more remote, so consistency
demands the same of the genuine forward citations.

---

## Section 2 — Direct-read findings

### F1: arXiv:2206.07891 (Abakumova-Frolovsky-Herbig-Lyakhovich 2022)

Local cache: `/tmp/forward_cit/2206.07891/main.tex` (16 kB, downloaded).

Citation context — the paper applies the 2102.10579 recipe to the linearised
Nordström equation `R = Λ`, casting it into Lagrangian form by introducing a
hook-tensor Stückelberg field. Direct quote (line 265):

> "Equations (LELGH) for spin two are obviously non-Lagrangian. Making use of
> the general Stueckelberg scheme of the article [Lyakhovich:2021lzy], these
> equations can be cast into Lagrangian setup."

And on iteration termination (line 266):

> "If the original system is linearised, the inclusion of Stueckelberg fields
> terminates at squares. In the non-linear case it continues to the higher
> orders."

**Higher-derivative content?** No. The third-order consequence with hook Young
diagram is a *kinematical* derivative consequence of the second-order Einstein
equations, not a 4th-order Lagrangian.

**Parameter-dependent rank-jump?** No. The construction is at fixed Λ and does
not study limits Λ → 0 or any other critical-surface limit.

**TIDAL relevance**: NONE. Same duality-construction usage pattern Round 1
Agent D already flagged. Confirms Agent D's reading.

### F3: arXiv:2303.02616 (Abakumova-Lyakhovich 2023, "Dualisation of free fields")

Local cache: `/tmp/forward_cit/2303.02616/main.tex` (24 kB, downloaded).

This is the natural follow-up to 2106.09355 — a *general* procedure for
constructing parent actions for dual formulations of free field theories.
Direct quote on the recipe scope (line 99, the methodology paragraph):

> "There is a systematic way to cast the involutive closure of the Lagrangian
> system back into Lagrangian framework by inclusion of Stueckelberg fields
> [Lyakhovich:2021lzy]. The method of the article [Lyakhovich:2021lzy] can be
> viewed as manifestly covariant counterpart of the procedure of converting
> the Hamiltonian second class constraints into first class ones by
> introducing 'conversion variable' for every second class constraint
> [Faddeev:1986pc]–[Batalin:2005df]."

And on the irreducibility of the iteration (line 190):

> "The procedure does not have obstructions [Lyakhovich:2021lzy], and it
> allows one to iteratively construct both the action and its gauge symmetry.
> For the free field theories, the iterative procedure have to terminate at
> quadratic terms in Stueckelberg fields in the action, while the gauge
> generators are field-independent, so they are immediately defined by
> zero-order boundary conditions."

**Higher-derivative content?** No. Worked examples are Proca (spin-1, 2nd-order),
massive and massless symmetric-tensor spin-2 (2nd-order). Hook Young-diagram
Stückelbergs sit on top of 2nd-order EoMs.

**Parameter-dependent rank-jump?** No. Section 4 examples are at fixed mass
parameters — no limit-of-coupling analysis at all. The "topological subsystem"
identified for dualisation is the differential-cohomological subsystem, not a
locus in parameter space where the principal-symbol rank drops.

**Constructive scheme worth noting**:  the paper's generalisation introduces a
*second* layer of Stückelbergs `ω^A` for "extra consequences being gauge
variations of the original action w.r.t. gauge transformations of the
topological subsystem" (line 38). This is structurally interesting because it
*expands* the auxiliary-field ansatz from Round 1 Agent D's irreducible single-
auxiliary lift to a reducible multi-layer lift.

**Does the multi-layer lift bridge `b5 → 0` for `b5·R̃²`?** No, on the same
structural grounds Round 1 Agent D's argument 3 establishes:

> "Reducibility null-vectors `Z^a` are b5-independent by construction (built
> from order-0 EoM consequence-generator structure). A b5-independent gauge
> structure cannot transform a b5-dependent Poisson bracket into a
> b5-independent one. QED."

The 2303.02616 procedure adds *more* b5-independent structure (the second-layer
`ω^A` and the topological subsystem of `τ_a`), which strengthens rather than
weakens Agent D's argument: each new Stückelberg `ω^A` introduces a new
auxiliary block whose Poisson bracket with the original phase-space variables
inherits the b5-dependence of the parent EoM, and the multi-layer reducibility
structure is built before any limit is taken.

**TIDAL relevance**: NONE. Confirms and strengthens Agent D's no-go.

### F5: arXiv:2402.12437 (Casadio-Chataignier 2024)

Local cache: `/tmp/forward_cit/2402.12437/main.tex` (51 kB, downloaded).

This paper cites 2102.10579 only in the bibliography (`StueckTrick2`, line 1227),
and the citation appears in a footnote distinguishing the "Fock-Stueckelberg
proper-time" mechanism (their subject) from the "Stueckelberg field-shifting"
mechanism (Lyakhovich's subject) — line 547:

> "The works [Stueck,Stueck1,Stueck2] of Stueckelberg have a conceptual
> similarity with but are distinct from the other famous Stueckelberg method
> [StueckTrick,StueckTrick1,StueckTrick2]: that of introducing a scalar field
> to preserve gauge invariance in a massive Abelian gauge theory of the
> Yang-Mills type."

So the citation is *explicitly* contrastive — the paper does not use Lyakhovich's
recipe. The paper itself is about treating first-class constraints as
"relaxation parameters" via Fock-Stueckelberg time, applied to GR and the problem
of time. It does not address constraint *rank* changes at all, let alone the
b5=0 surface.

**Higher-derivative content?** No.
**Parameter-dependent rank-jump?** No.
**TIDAL relevance**: NONE.

### F6: arXiv:2405.13706 (Delplanque-Skvortsov 2024)

Local cache: `/tmp/forward_cit/2405.13706/main.tex` (1254 lines, downloaded).

The citation to `Abakumova:2021evc` (= 2106.09355) is in a list of "covariant
ideas" for higher-spin interactions on line 124:

> "Some other approaches to the problem of interactions of massive higher spin
> fields include covariant ideas [...,Abakumova:2021evc,Abakumova:2023wve,
> Skvortsov:2023jbn] and the light-cone gauge [...]."

This is a one-shot bibliographic mention — Lyakhovich's recipe is not invoked
in the body of the paper. The paper itself works with massive higher-spin fields
in the Singh-Hagen / Zinoviev framework, with second-class constraints
converted to first-class by adding a tower of auxiliary fields. The conversion
is well-known (Zinoviev 2001) and predates Lyakhovich's general recipe by 20
years.

**Higher-derivative content?** No (Lagrangian is 2nd-order in fields).
**Parameter-dependent rank-jump?** No — the entire paper assumes generic mass
`m ≠ 0`; the `m → 0` partially-massless limit is considered (footnote 7), but
that limit is the *opposite* topology to TIDAL's b5 → 0 (it removes Stückelbergs
rather than freeing constraint-promoted modes).

**TIDAL relevance**: NONE.

---

## Section 3 — Aashish-Saif 2026 status (`literature/2601.22007/main.tex`)

The user's task brief flagged this paper as "the FJ-related reference Round 2
Agent G actually cited" and asked whether it has constructive content for
TIDAL's case.

**Result**: arXiv:2601.22007 does **NOT** cite Lyakhovich (2102.10579) or
Abakumova-Lyakhovich (2106.09355). I grep'd the entire `main.tex` and `Ref.bib`
for "Lyakhovich" and "Abakumova" — zero hits.

The paper's Stückelberg references are
(`main.tex:50`):

> "The quantization of such models involves the use of Stückelberg mechanism
> [Stueckelberg_38-1, Stueckelberg_38-2, RueggRuizAltaba_04, Buchbinder_08,
> Hinterbichler_12, Govindarajan_25] wherein a so-called Stückelberg field
> (a vector field $C_\mu$ in case of antisymmetric tensor field model) is
> introduced in the Lagrangian to restore the gauge symmetry so as to apply
> standard quantization procedures like the Faddeev-Popov method
> [Faddeev_67, Peskin_95]."

This is the *original* 1938 Stückelberg-Higgs mechanism + Hinterbichler 2012
review — entirely different from Lyakhovich's "general method".

What Aashish-Saif actually does: apply the *standard* Stückelberg trick (single
auxiliary vector $C_\mu$) to a Kalb-Ramond/antisymmetric-tensor model with
spontaneous Lorentz violation, where the *bumblebee-vacuum-induced* singularity
of the Hamiltonian on the vacuum manifold is the problem (Seifert 2019,
arXiv:1810.09584). Direct quote from abstract:

> "Spontaneous Lorentz violation models of antisymmetric tensor field are
> known to possess singular Hamiltonian on the vacuum manifold, leading to
> unresolvable pathologies that render such theories unfit for cosmological
> studies. In this work, we show that by introducing an auxiliary vector field
> inspired by the Stückelberg mechanism to restore the gauge symmetry of the
> Lagrangian, it is possible to resolve such pathologies on vacuum manifold.
> The constraint analysis using Dirac-Bergmann method leads to a constraint
> matrix that acquires dependence on gradients and conjugate momentum of the
> Stückelberg field and therefore remains non-singular on the vacuum manifold."

Structural mapping to TIDAL:
- Their "singularity": `det(M_constraint)` vanishes on a *field-configuration
  submanifold* (vacuum manifold of the antisymmetric tensor — `B_{ab} = b_{ab}`).
- TIDAL's "singularity": `det(M_constraint) ∝ b5` vanishes on a *parameter
  submanifold* (`b5 = 0`).

These are different categories. Aashish-Saif's mechanism — the Stückelberg
field shifts the constraint matrix by terms in `∇C` and `π_C` so that the
matrix is no longer evaluable purely on `B = b` — works because it changes the
field-configuration where the determinant vanishes. **It cannot work for
TIDAL's case** because the parameter `b5` does not appear in the auxiliary
field `C_\mu`; introducing a Stückelberg field cannot make `b5` not appear in
the Poisson matrix.

This is the same structural reason Round 1 Agent D's argument 3 gives —
b5-independent gauge structure cannot bridge a b5-dependent rank jump.
Aashish-Saif doesn't help, doesn't even attempt to help, and (correctly) does
not cite Lyakhovich's machinery for help.

**TIDAL relevance**: NONE. Different singularity category. Not a forward
citation of either Lyakhovich paper despite being in the same year and topic.
Round 2 Agent G's "FJ-related reference" pointer was misleading — Aashish-Saif
is not Faddeev-Jackiw, it is Dirac-Bergmann; and the singularity-resolution
mechanism is non-transferable to parameter-rank-jump case.

---

## Section 4 — Verdict

**Verdict: (b) NO NEW LEAD.**

The forward-citation set is small (3 unique non-self-citing papers in 2022-2026)
and uniformly applies the Lyakhovich-Sharapov recipe in the way Round 1 Agent
D characterised: as a duality-construction tool at fixed parameter values.
None of the forward-citation papers attempts to extend the recipe to handle
parameter-dependent rank changes in the constraint Hessian, and none applies
it to a higher-derivative theory.

Specifically:
- **F1, F3, F4** (Abakumova-Lyakhovich 2022/2023) all stay within the
  duality / hook-Young-diagram framework on 2nd-order EoMs. F3 generalises the
  recipe with a multi-layer reducible structure but the new layers are still
  built from b5-independent topological subsystems, so Agent D's structural
  argument applies a fortiori.
- **F2** (Guleryuz 2023) is a string-landscape paper with only background
  citation.
- **F5** (Casadio-Chataignier 2024) explicitly *contrasts* itself with the
  Lyakhovich approach.
- **F6** (Delplanque-Skvortsov 2024) cites only as background in a list of
  covariant approaches.
- **Aashish-Saif 2026** (arXiv:2601.22007) does not cite Lyakhovich at all
  and addresses a structurally different singularity category
  (field-configuration vacuum manifold, not parameter critical surface).

This confirms Round 1 Agent D's no-go is **robust under the 2022-2026
literature**.  The gap stays open: there is no published recipe that bridges
the b5=0 critical surface for a 4th-order PGT Lagrangian.

### Concrete next steps

This investigation is closed. No new download is needed for the literature
cache (the relevant forward-citation papers are not load-bearing — they should
be flagged in `MANUAL_RETRIEVAL_NEEDED.md` as low-priority for completeness
but no claim in any TIDAL writeup depends on them).

Implications for the TIDAL roadmap:
1. **Do not chase Lyakhovich extensions.**  The 2022-2026 literature has
   produced no constructive technique for the parameter-rank-jump case.
2. **Aashish-Saif 2026 should not be cited as a Stückelberg lead.** It
   addresses a different singularity category and does not actually use
   Lyakhovich's machinery.
3. **The FINAL_ASSESSMENT.md "Optional Tier-2 work" framing is unchanged**:
   if the project later needs the metric h_4/h_7/h_9 sector specifically (it
   doesn't, per Meta-L), Glavan-Zlosnik-Lin 2024 (`literature/2311.17459/`)
   remains the methodologically closest published work, not any Lyakhovich
   forward citation.

### Bookkeeping

The 4 forward-citation papers have been downloaded to `/tmp/forward_cit/`
during this investigation but are **not** copied to `literature/` because no
TIDAL claim depends on them.  They are catalogued in
`MANUAL_RETRIEVAL_NEEDED.md` (updated below) as "tier-3 — non-load-bearing".

---

## Section 5 — Citations

### Primary papers

- **Abakumova-Lyakhovich 2021**, "Reducible Stueckelberg symmetry and
  dualities", Phys. Lett. B 820:136552, arXiv:2106.09355.
  Local cache: `literature/2106.09355/main.tex`.
- **Lyakhovich 2021**, "General method for including Stueckelberg fields",
  Eur. Phys. J. C 81:472, arXiv:2102.10579.
  Local cache: `literature/2102.10579/main.tex`.

### Forward-citation set (Inspire-HEP REST API, 2026-04-27)

Authoritative listing:
- https://inspirehep.net/api/literature?q=refersto%3Arecid%3A1869064%20OR%20refersto%3Arecid%3A1847795 (9 hits, 7 unique non-source non-tangential)
- Inspire record IDs verified via `q=arxiv:<id>&fields=control_number`:
  2106.09355 → 1869064, 2102.10579 → 1847795.

Forward-citation papers (in publication order):
- arXiv:2206.07891 (Abakumova, Frolovsky, Herbig, Lyakhovich) Eur. Phys. J. C 82:780 (2022). https://arxiv.org/abs/2206.07891. Self-citation; linearised-Nordström / hook-tensor duality. Local cache: `/tmp/forward_cit/2206.07891/main.tex`.
- arXiv:2207.10634 (Guleryuz) JCAP 05:039 (2023). https://arxiv.org/abs/2207.10634. String-landscape attractors. Background citation only.
- arXiv:2303.02616 (Abakumova, Lyakhovich) Annals Phys. 453:169322 (2023). https://arxiv.org/abs/2303.02616. Self-citation; multi-layer reducible Stückelberg / parent-action dualisation. Local cache: `/tmp/forward_cit/2303.02616/main.tex`.
- Bull. Lebedev Phys. Inst. 50(3):108 (2023), Abakumova-Lyakhovich, "Dual Formulation for the Massless Spin Two Theory". Self-citation, proceedings note subsumed by 2303.02616.
- arXiv:2402.12437 (Casadio, Chataignier) Annals Phys. 470:169783 (2024). https://arxiv.org/abs/2402.12437. Fock-Stueckelberg time relaxation; explicit contrast with Lyakhovich method. Local cache: `/tmp/forward_cit/2402.12437/main.tex`.
- arXiv:2405.13706 (Delplanque, Skvortsov) Class. Quant. Grav. 41(24):245018 (2024). https://arxiv.org/abs/2405.13706. Massive higher-spin chiral vs symmetric. Bibliographic mention only. Local cache: `/tmp/forward_cit/2405.13706/main.tex`.
- Erickson PhD thesis, Imperial 2022. Off-topic (SUSY localization).
- arXiv:2008.02112 (Pandey 2022). Excluded — earlier preprint with tangential mention.

### Tangential — Aashish-Saif 2026

- arXiv:2601.22007 (Aashish, Saif) "Stückelberg inspired approach for
  avoiding singular Hamiltonians in Lorentz violating models of antisymmetric
  tensor field", 2026. https://arxiv.org/abs/2601.22007.
  Local cache: `literature/2601.22007/main.tex` (+ `Ref.bib`).
  **Does not cite either Lyakhovich paper.** Different singularity category
  (field-configuration vacuum manifold, Seifert 2019).

### Verbatim quote pointers

All verbatim quotes in this lead file are tagged with their source file and
line number. Reproduced here for completeness:

- 2206.07891 line 265: "Equations (LELGH) for spin two are obviously
  non-Lagrangian. Making use of the general Stueckelberg scheme of the article
  [Lyakhovich:2021lzy], these equations can be cast into Lagrangian setup."
- 2206.07891 line 266: "If the original system is linearised, the inclusion of
  Stueckelberg fields terminates at squares. In the non-linear case it
  continues to the higher orders."
- 2303.02616 line 99: "There is a systematic way to cast the involutive closure
  of the Lagrangian system back into Lagrangian framework by inclusion of
  Stueckelberg fields [Lyakhovich:2021lzy]. The method [...] can be viewed as
  manifestly covariant counterpart of the procedure of converting the
  Hamiltonian second class constraints into first class ones by introducing
  'conversion variable' for every second class constraint."
- 2303.02616 line 190: "The procedure does not have obstructions
  [Lyakhovich:2021lzy], and it allows one to iteratively construct both the
  action and its gauge symmetry. For the free field theories, the iterative
  procedure have to terminate at quadratic terms in Stueckelberg fields in the
  action [...]."
- 2402.12437 line 547: "The works [Stueck,Stueck1,Stueck2] of Stueckelberg have
  a conceptual similarity with but are distinct from the other famous
  Stueckelberg method [StueckTrick,StueckTrick1,StueckTrick2]: that of
  introducing a scalar field to preserve gauge invariance in a massive Abelian
  gauge theory of the Yang-Mills type."
- 2405.13706 line 124: "Some other approaches to the problem of interactions
  of massive higher spin fields include covariant ideas
  [...,Abakumova:2021evc,Abakumova:2023wve,Skvortsov:2023jbn] and the
  light-cone gauge."
- 2601.22007 abstract: "Spontaneous Lorentz violation models of antisymmetric
  tensor field are known to possess singular Hamiltonian on the vacuum
  manifold [...]. In this work, we show that by introducing an auxiliary
  vector field inspired by the Stückelberg mechanism to restore the gauge
  symmetry of the Lagrangian, it is possible to resolve such pathologies on
  vacuum manifold."
- Round 1 Agent D's argument 3 (from `notes/round1_synthesis.md` line 93):
  "Reducibility null-vectors `Z^a` are b5-independent by construction (built
  from order-0 EoM consequence-generator structure). A b5-independent gauge
  structure cannot transform a b5-dependent Poisson bracket into a
  b5-independent one. QED."
