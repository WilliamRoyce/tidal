# Lead — Beltrán Jiménez / Heisenberg / Koivisto metric-affine corpus

**Date:** 2026-04-27
**Investigator:** Claude (Opus 4.7, harness session)
**Time spent:** ~75 min (Phase 2.5 of post-FINAL_ASSESSMENT investigation)
**Context:** Round 1+2 agents flagged BJHK metric-affine + Stückelberg work
without local TeX access. Three papers were claimed to be locally cached;
audit instructed to read them directly to test transferability claims.

---

## Papers

| # | Cite | Local file | Status |
|---|------|------------|--------|
| 1 | de la Cruz-Dombriz, Maldonado Torralba, Mazumdar, "Ghost-free higher-order theories of gravity with torsion", arXiv:1911.08846v1 (2019) | (not in `literature/`) `/tmp/lit_fetch/1911.08846.txt` (pdftotext from arxiv.org) | Read in full |
| 2 | Heisenberg, Hohmann, "Gauge-invariant cosmological perturbations in general teleparallel gravity", arXiv:2311.05597v2 (2023) | (not in `literature/`) `/tmp/lit_fetch/2311.05597.txt` (pdftotext from arxiv.org) | Read in full |
| 3 | Saha, Sanyal, "Vulnerability of f(Q) gravity theory and a possible resolution", arXiv:2503.18972 (2025) | `/workspaces/torsion-gertsenshtein/literature/2503.18972/f_R,Q_.tex` (local) | Read in full |

> **NOTE.** The audit prompt asserted that papers 1 and 2 are cached at
> `literature/1911.08846/` and `literature/2311.05597/`. They are **not**.
> `ls literature/` returns only `2503.18972/` from this trio. I retrieved the
> two missing papers via `curl https://arxiv.org/pdf/<id>` + `pdftotext -layout`
> to `/tmp/lit_fetch/`. The local-paper claim in the audit prompt is wrong;
> readers of this note should not rely on `literature/{1911.08846,2311.05597}/`
> existing. (Optional follow-up: stash these two PDFs into `literature/` per
> the project's literature workflow.)

---

## Section 1 — 1911.08846 (de la Cruz-Dombriz / Maldonado Torralba / Mazumdar)

### 1.1 Authors and provenance

This is **NOT a BJHK paper**, despite the audit prompt's framing. The
author list is `Álvaro de la Cruz-Dombriz, Francisco José Maldonado Torralba,
Anupam Mazumdar` (`/tmp/lit_fetch/1911.08846.txt:2`). It builds on Biswas-
Gerwick-Koivisto-Mazumdar IDG (`[12]` = arXiv:1110.5249, *Phys. Rev. Lett.*
**108**:031101). Koivisto is a co-author of the IDG foundation paper but
not of this 2019 manuscript; Heisenberg and Beltrán Jiménez are not authors.
Reference `[43]`:

> "J. Beltrán Jiménez and F. J. Maldonado Torralba, Revisiting the Stability
> of Quadratic Poincaré Gauge Gravity, 1910.07506."
> (`/tmp/lit_fetch/1911.08846.txt:1101`)

is cited only for the IR-ghost-freedom condition on `∇μ K^μνρ ∇σ K^σνρ`
terms. The paper is more accurately characterised as **"infinite-derivative
gravity with torsion"** (IDG-T), in the de la Cruz-Dombriz/Mazumdar lineage.

### 1.2 Starting Lagrangian

The full action (eq. 4) sums **47** quadratic-in-curvature/contortion
invariants R̃F̃₁(□)R̃ + R̃F̃₂(□)∂μ∂νR̃^μν + ... + KμνρF̃ᵢ(□)∂...K^σλ + R̃μνλσF̃₄₇(□)R̃^μλνσ
(`/tmp/lit_fetch/1911.08846.txt:138-180`), where each `F̃ᵢ(□)` is a
*Krasnikov-Tomboulis non-local form factor*

> "F̃ᵢ(□) = Σ_{n=0}^N f̃ᵢ,ₙ (□/M_S²)ⁿ"
> (`/tmp/lit_fetch/1911.08846.txt:225-231`)

with `n` allowed to be **infinite** ("infinite derivative theories of
gravity, IDG"). The local-limit Lagrangian (`/tmp/lit_fetch/1911.08846.txt:308-311`):

> "L_GPG = R̃ + b₁R̃² + b₂R̃μνρσR̃^μνρσ + b₃R̃μνρσR̃^ρσμν + 2(b₁−b₂−b₃)R̃μνρσR̃^μρνσ + b₅R̃μνR̃^μν − (4b₁+b₅)R̃μνR̃^νμ + a₁KμνρK^μνρ + a₂KμνρK^μρν + a₃Kν^μ_μK^νρ_ρ + c₁K^μνρ∇μ∇σK^σνρ + c₂K^μνρ∇μ∇σK^σρν + c₃Kμρ^ν∇ρ∇σK^μνσ + c₄Kμρ^ν∇ρ∇σK^μσν"

is **purely parity-even**. The closing remark of §IV is decisive
(`/tmp/lit_fetch/1911.08846.txt:572-573`):

> "This is indeed possible since parity breaking terms in the action are
> not considered, so there are no mixed trace-axial terms."

There is **no Holst-squared b̄₅·R̃² term** anywhere in the action. TIDAL's
critical-surface coupling `b5·R̃²` (with R̃ = the parity-odd Holst dual,
indices contracted via Levi-Civita ε) lies outside the universe of theories
this paper considers.

### 1.3 Ghost-freedom mechanism

The technique is the **Krasnikov-Tomboulis exponential-of-entire-function**
form factor. The five "ghost-freedom conditions" (`/tmp/lit_fetch/1911.08846.txt:537-551`) read:

> "p₂(□) − p₁(□) = C₁ e^{−□/M_S²}
>  q₁(□) − q₂(□) − q₃(□) + q₄(□) = C₂ e^{−□/M_S²}
>  2p₁(□) + p₂(□) + 3p₃(□) + ½ s(□)□² = C₃ e^{−□/M_S²}
>  q₁(□) + q₃(□) + 2q₄(□) − 3q₆(□) = C₄ e^{−□/M_S²}
>  3q₁(□) + 3q₂(□) + 9q₅(□) − s(□)□ = C₅ e^{−□/M_S²}"

Three additional sign conditions (`C₅ = 0`, `C₂ > 0`, `C₄ < 0`,
`/tmp/lit_fetch/1911.08846.txt:566-569`) eliminate the longitudinal
ghost in the trace vector, give the axial vector positive kinetic energy,
and give the trace vector negative-of-kinetic-energy with the right sign
to be a Proca field.

This is **not** an auxiliary-field lift, **not** a Stückelberg construction,
**not** a Lagrange-multiplier promotion, and does **not** treat any
critical surface. It is a class-restriction: pick the form factors so the
propagator is `e^{p²/M_S²}/p²`, with no extra zeros/poles, by construction.

### 1.4 Treatment of small-parameter / coupling limits

There is a single small-parameter limit in the paper, the IR / local limit
`M_S → ∞`. In that limit (Appendix C), the form factors collapse to their
zeroth Taylor coefficients and one recovers a local PGT Lagrangian. There
is **no treatment of any coupling parameter going to zero** (no `b5 → 0`,
no `M_S → 0`, no rank-jump analysis). The "stability" enforcement
`v₂(□) = −3u(□)` in eq. (29-30) (`/tmp/lit_fetch/1911.08846.txt:443`)
is a *condition imposed* on the form factors to suppress the trace-vector-
curvature non-minimal coupling, not a study of what happens when it fails.

### 1.5 Constraint promotion / rank jump

Not addressed. The whole construction is engineered so the propagator
has *no* extra poles to begin with — there is no constraint structure to
promote. The local-limit field equations (Appendix C) reduce to local PGT
with eq. (C7) parameter map, but the analysis stops there.

### 1.6 Verdict for paper 1

**Not transferable to TIDAL's b5·R̃² critical surface.** Three structural
gaps:

(i) The ghost-freedom mechanism is *non-locality from infinite-derivative
form factors*, not anything that talks to TIDAL's local b5·R̃² Ostrogradsky
rank-jump.
(ii) The Lagrangian space considered is *parity-even* (eq. 16; explicit
disclaimer at line 572-573). TIDAL's b5 multiplies the parity-odd R̃²
(Holst-squared); this paper says nothing about that sector.
(iii) There is no parameter limit `M_S → 0` or `b₅ → 0` analysed; the
authors restrict to entire-function form factors and never study the
boundary of that class.

A constructive transfer would require: (a) extending the 47-term Lagrangian
ansatz to include parity-odd invariants, (b) deriving ghost-freedom
conditions for the parity-odd sector, (c) analysing the local-limit map
(C7) when one of the local couplings (e.g. `b5`) approaches zero. Each
of (a-c) is a separate research project of comparable difficulty to the
original paper.

---

## Section 2 — 2311.05597 (Heisenberg, Hohmann)

### 2.1 Authors and provenance

This **is** a BJHK-adjacent paper — Lavinia Heisenberg is the H of BJHK.
Co-author Manuel Hohmann (Tartu) is a long-standing teleparallel-gravity
collaborator; Beltrán Jiménez and Koivisto are not authors but appear as
references `[20]` (CANTATA review) and `[21]` (Geometrical Trinity).

### 2.2 Starting framework

**General teleparallel gravity** = flat affine connection (curvature
`R^ρ_σμν = 0`, `/tmp/lit_fetch/2311.05597.txt:97-98`) with both torsion
and nonmetricity allowed. The dynamical fields are `g_μν` and `Γ^μ_νρ`
constrained by flatness. This is fundamentally **teleparallel** (zero
curvature) — the opposite of Riemann-Cartan / PGT, where curvature is
non-zero by construction.

### 2.3 Gauge-invariant perturbation framework

The technique is **standard Bardeen-style gauge-invariant decomposition**
adapted to teleparallel geometry. Method:

> "In particular, we have decomposed the perturbations of the metric and
> the teleparallel connection, the matter variables and the general form
> of the gravitational field equations into irreducible components, and
> studied their behavior under gauge transformations. Making use of the
> latter, we have shown how gauge-invariant quantities can be constructed
> from these components."
> (`/tmp/lit_fetch/2311.05597.txt`, conclusion)

Concretely: scalar/vector/tensor split (eqs. 9-15 + later); diffeomorphism
gauge transformations parameterised by ξ^μ; combinations of
perturbation variables that are invariant under those transformations.

### 2.4 Background and parameter dependence

**FLRW only.** "linear cosmological perturbations" / "around a homogeneous
and isotropic background geometry" (`/tmp/lit_fetch/2311.05597.txt:9-12`).
The classification of cosmologically symmetric backgrounds is taken from
refs. [23-26] (Hohmann et al. 2019-2022).

The paper applies the framework to f(T), f(Q), f(G) (eqs. 156, 162, 173)
and derives modified tensor-perturbation equations. Notable findings:

> "the tensor perturbations in the f(G) class of gravity theories are
> formally described by the same set of equations as in general relativity
> and in f(T) gravity, and observational differences can appear only
> through a modified background evolution influencing the cosmological
> background value of the function f and its derivatives in the tensor
> field equations."
> (`/tmp/lit_fetch/2311.05597.txt`, end of §VII.C)

There is **no analysis of parameter-dependent constraint structure**, no
rank jump, no critical-surface treatment. A keyword scan returns zero
hits for "strong coupling", "critical", "Stückelberg", "Holst", "rank",
"b5". The only "Lagrange multiplier" mention is the standard variational-
principle device for imposing flatness or vanishing-torsion (lines 213,
250).

### 2.5 Higher-derivative / Holst content

None. The actions f(T), f(Q), f(G) are arbitrary functions of *quadratic*
torsion / nonmetricity / general scalars, but the perturbation equations
(157-181) are second-order in time. There is no R²/R̃²/Holst-squared
structure and no Ostrogradsky-type analysis.

### 2.6 Verdict for paper 2

**Cosmological-background-specific; not transferable to TIDAL.** Four
structural reasons:

(i) **Wrong geometry.** Teleparallel = flat connection; PGT/Riemann-Cartan
= non-flat connection. The Bardeen-style classification of perturbations
in this paper does not extend to non-flat-connection backgrounds without
additional work.
(ii) **Wrong background.** FLRW only; TIDAL's Gertsenshtein conversion
runs on flat-Minkowski with localised wave-packet perturbations.
(iii) **Wrong order in time.** Second-order field equations throughout;
TIDAL's b5·R̃² is fourth-order in time for the constraint-promoted fields.
(iv) **No parameter-dependent constraint structure.** The paper does not
encounter any singular-limit phenomenon analogous to b5 → 0. f(T), f(Q),
f(G) classes are treated *generically*; the strong-coupling literature
on f(Q) (which the paper references via [22, 31] etc.) is not engaged.

---

## Section 3 — 2503.18972 (Saha, Sanyal)

### 3.1 Authors and provenance

Dalia Saha, Abhik Kumar Sanyal — Indian-school constraint-analysis lineage.
Their predecessor work (cited as [DA1, DA2]) is on Hamiltonian formulation
of GMTG (f(T)) and GSTG (f(Q)). Not BJHK-adjacent — they cite Heisenberg
only via H3 (D'Ambrosio-Heisenberg-Zentarra arXiv:2308.02250 "Failure of
Dirac-Bergmann Algorithm for Teleparallel Theories"; line 554 of f_R,Q_.tex).

### 3.2 The "vulnerability"

GSTG (= f(Q) with `f,QQ ≠ 0`) suffers (from abstract,
`literature/2503.18972/f_R,Q_.tex:20`):

> "the strong coupling issue and the ghost degrees of freedom. It has also
> been cognised that GSTG does not admit diffeomorphic invariance in
> general."

Specifically (line 27):

> "the Hamiltonian formulation of f(Q) gravity not only brings about the
> ghosts but also fails to establish diffeomorphic invariance [H1,H2] or
> even the 'Dirac-Bergmann algorithm' itself, due to the fact that the
> inhomogeneous system of partial differential equations required to
> check consistency condition cannot be solved exactly [H3]. These
> problems stem from the fact that the shift vector (N_i) and the
> auxiliary variable (Φ) act as non-propagating dynamical variables and
> the auxiliary variable turns out to be the ghost."

### 3.3 The "possible resolution"

**Add a Ricci-scalar term.** The proposed modified action (line 191):

> "f(R,Q) = α₁R + βR² + α₂Q + γQ²"

with R the Levi-Civita Ricci scalar and Q the standard f(Q) non-metricity
scalar. The mechanism (line 27):

> "the higher-order terms appearing in the function f(R) take control of
> higher degree terms [KNRA] associated with f(Q) theory, and hence the
> 'ghost degree of freedom' does not appear."

The Hamiltonian construction (§4.1, "Modified Dirac-Bergmann constraint
analysis") proceeds:
- introduce auxiliary variable `x = ż/N` via Lagrange multiplier `u` (eq.
  DA2, line 202);
- compute primary Hamiltonian (eq. AHp, line 208);
- enforce two second-class constraints `φ₁ = N p_z − u ≈ 0`, `φ₂ = p_u ≈ 0`
  (line 227);
- consistency conditions `φ̇ᵢ = {φᵢ, H_p1} ≈ 0` give Lagrange multipliers
  `u₁ = x`, `u₂ = −N ∂H_p1/∂z` (eq. Apc, line 234);
- final reduced Hamiltonian (eq. AHF, line 248) is `H = N[xp_z + √z·p_x²/(36β) + 3αx²/(2√z) − 9γx⁴/(4z^{5/2})]`,
  manifestly free of branched/divergent terms.

### 3.4 Does this generalise to PGT b5=0?

**No** — the resolution is f(Q)-specific in a precise sense. Three reasons:

(i) **Wrong sector.** The "strong coupling" diagnosed in f(Q) is the
inability to invert `p_z = ∂L/∂ż` into `ż(p_z, z, ...)` because the lapse
N enters the kinetic term non-perturbatively. Adding a `βR²` term gives a
new propagating scalar (the standard f(R) scalaron) which carries the
extra DOF cleanly. This is the **f(R) trick for fixing ill-defined
canonical structure in f(Q)**. TIDAL's b5·R̃² problem is opposite: the
problem there is that *too many* fields become dynamical when you turn
on b5 — h_4, h_7, h_9 stop being Lagrange multipliers and acquire fourth-
order kinetic terms. Adding more curvature invariants makes that *worse*,
not better.

(ii) **Wrong constraint topology.** Saha-Sanyal's f(Q) issue is that the
Dirac-Bergmann algorithm halts because the consistency conditions cannot
be solved (line 27 quote above). Their fix is to introduce *more*
propagating modes (R² scalaron) so the consistency conditions become
solvable. TIDAL's b5·R̃² issue is that the consistency conditions
*change qualitative form* across b5 = 0: at b5 = 0 the algorithm produces
algebraic constraints; at b5 ≠ 0 those constraints turn into 4th-order
kinetic equations. Saha-Sanyal do not analyse a parameter-dependent change
in constraint topology.

(iii) **Wrong background.** Saha-Sanyal restrict to spatially flat /
spatially curved RW minisuperspace throughout. The shift vector N_i is
explicitly dropped to side-step the diffeomorphism-invariance issue
(line 158). TIDAL's PGT analysis is on flat Minkowski with local
perturbations; the shift-vector N_i sector is intrinsic to the rank-jump,
not bypassable.

### 3.5 Verdict for paper 3

**f(Q)-specific, not transferable.** The resolution adds a propagating
scalaron via `βR²` to fix an *under-constrained* f(Q) Hamiltonian; TIDAL's
b5·R̃² produces an *over-constrained-becoming-Ostrogradsky* PGT Hamiltonian.
The two pathologies are dual. Adding R² to TIDAL's Lagrangian would make
the rank-jump worse, not better, because R² also produces propagating
4th-order modes.

---

## Section 4 — New references identified

### 4.1 From paper 1's bibliography

- **Ref [38]**: Koivisto-Tsimperis 1810.11847, "The spectrum of teleparallel
  gravity". Also extends IDG to teleparallel; same lineage. **Not transferable**
  by the same teleparallel-vs-PGT argument as paper 2.
- **Ref [39]**: Conroy-Koivisto 1710.05708, "The spectrum of symmetric
  teleparallel gravity". Same comment.
- **Ref [43]**: Beltrán Jiménez-Maldonado Torralba 1910.07506, "Revisiting
  the Stability of Quadratic Poincaré Gauge Gravity". This is a true BJHK-
  adjacent paper. Cited only for the IR ghost-freedom condition on
  `∇μ K^μνρ ∇σ K^σνρ` terms — i.e., for *qualitative consistency* with
  paper 1's local limit, not for transferable techniques. **Worth a
  separate read** if the project later needs the parity-even / vector-mode
  sector. Out of scope for this Phase 2.5.
- **Refs [12, 13]**: Biswas-Gerwick-Koivisto-Mazumdar 1110.5249 and
  Biswas-Koshelev-Mazumdar 1602.08475, the IDG foundational papers. Same
  framework as paper 1.

### 4.2 From paper 2's bibliography

- **Ref [21]**: Beltrán Jiménez-Heisenberg-Koivisto 1903.06830, "The
  Geometrical Trinity of Gravity". Foundational review of the metric /
  metric-teleparallel / symmetric-teleparallel equivalence. Theoretical
  framework; **not a constructive Hamiltonian-recipe paper**.
- **Ref [37]**: Beltrán Jiménez-Heisenberg-Koivisto 1710.03116, "Coincident
  General Relativity". Original f(Q) ↔ STEGR paper. Not Hamiltonian.
- **Ref [38]**: Heisenberg-Hohmann-Kuhn 2311.05495, "Cosmological
  teleparallel perturbations". Companion paper to 2311.05597. Same scope
  / same limitations (FLRW, generic f, no parameter-dependent constraint
  structure). **Not transferable** by the same argument as paper 2.

### 4.3 From paper 3's bibliography (line numbers in `f_R,Q_.tex`)

- **[H3]** D'Ambrosio-Heisenberg-Zentarra, "Hamiltonian Analysis of f(Q)
  gravity and the Failure of Dirac-Bergmann Algorithm, for Teleparallel
  Theories of Gravity", *Fortschr. Phys.* **71**, 2300185 (2023),
  arXiv:**2308.02250** (cited at `f_R,Q_.tex:554`). This **is** a BJHK-
  adjacent paper (Heisenberg co-author) and **is** about a constraint-
  algorithm pathology in a metric-affine setting. **Worth checking**
  whether it cites or extends TIDAL-relevant techniques. Was not flagged
  in `MANUAL_RETRIEVAL_NEEDED.md`. *Recommended retrieval.*
- **[H1, H2]** unspecified at this scan; line 553 has Tomonari-Bahamonde
  "Dirac-Bergmann analysis and degrees of freedom of coincident f(Q)-
  gravity", *Eur. Phys. J. C* **84**, 349 (2024). Same DBCA-on-f(Q)
  lineage. Probably not transferable but worth a fast check.

### 4.4 Recommendation for `MANUAL_RETRIEVAL_NEEDED.md`

Add to a new "Tier 2 — BJHK-adjacent leads" section:

| Paper | arXiv | Why |
|---|---|---|
| Beltrán Jiménez-Maldonado Torralba 2019, "Revisiting Stability of Quadratic PGT" | 1910.07506 | Cited by paper 1 for IR ghost-freedom; closest BJHK-style paper to TIDAL's PGT setting found in this lead |
| D'Ambrosio-Heisenberg-Zentarra 2023, "Failure of Dirac-Bergmann Algorithm for Teleparallel Theories" | 2308.02250 | Hamiltonian-pathology paper; BJHK-adjacent; not previously flagged |

Both are forward-citation candidates if the project later needs to
revisit the constraint-promotion barrier with fresh eyes. Neither is on
the critical path for the headline `h_5 ↔ a_1` channel.

---

## Section 5 — Verdict

### **(b) NO TRANSFER.**

All three papers are scope-restricted in ways that prevent transfer to
TIDAL's flat-Minkowski PGT b5·R̃² critical surface:

| Paper | Scope restriction | What it would take to transfer |
|---|---|---|
| 1911.08846 | Parity-even local Lagrangian + non-local form-factor ghost-freedom; no `b5` limit | Extend to parity-odd Lagrangian + new ghost-freedom conditions in the parity-odd sector + analyse local limit at `b5 → 0` |
| 2311.05597 | FLRW background + flat-connection (teleparallel) + 2nd-order EOMs | Replace teleparallel with non-flat PGT + extend gauge-invariant decomposition to flat-Minkowski with localised wave-packet perturbations + handle 4th-order time derivatives |
| 2503.18972 | f(Q) DBCA pathology of *too few* propagating modes; minisuperspace only | Saha-Sanyal's mechanism (add R² scalaron) goes the wrong direction — adding kinetic invariants makes TIDAL's *over*-propagating Ostrogradsky structure worse |

Each transfer is, in effect, a separate full research programme. None of
these three papers gives a constructive recipe one can apply to TIDAL
within a sitting. **The constraint-promotion barrier in PGT b5·R̃² stays
open**, consistent with Meta-N's verdict on BC 2018 Appendix D.

### What this confirms

- The Round 1+2 agents' framing "BJHK metric-affine + Stückelberg work"
  was overstated. Of the three papers nominated as load-bearing, only
  paper 2 (2311.05597) is genuinely BJHK-adjacent (Heisenberg author).
  Paper 1 is a Maldonado-Torralba/Mazumdar IDG-T paper that cites BJ-MT
  1910.07506; paper 3 is a Saha-Sanyal f(Q) paper that does not cite any
  BJHK author directly.
- None of the three contains a Stückelberg lift, an auxiliary-field
  promotion, or any constraint-promotion recipe that addresses TIDAL's
  case. The "BJHK Stückelberg" framing in earlier rounds appears to have
  been generated from titles + abstracts without source-code reading.
- The audit-prompt assertion that papers 1 and 2 are cached locally is
  factually incorrect: only 2503.18972 is in `literature/`. Future leads
  should verify path existence before assuming local access.

### Recommended next steps

1. **Stash the two retrieved PDFs** to `literature/1911.08846/` and
   `literature/2311.05597/` per the project's literature workflow, so the
   next agent doesn't have to re-fetch.
2. **Add the two Tier-2 leads** (1910.07506 BJ-MT 2019, 2308.02250 DHZ 2023)
   to `MANUAL_RETRIEVAL_NEEDED.md` if the project later needs to extend
   the BJHK lead. Skip if the FINAL_ASSESSMENT operational pivot stays in
   force.
3. **No documentation change needed in `docs/tex/`.** This Phase 2.5
   investigation does not produce new claims to add or retract — it
   confirms the FINAL_ASSESSMENT's recommendation that the Stückelberg /
   metric-affine lead does not advance the headline observable.
4. **Do not reopen Path B-trace or Path B-tensor-q.** Paper 1's parity-
   even-only restriction is now a verified literature gap, reinforcing
   Reviews 2 + Meta-K's finding that Barker 2024 and Chatzistavrakidis-
   Ranjbar-Zekoč 2024 do not have published parity-odd extensions.

---

## Section 6 — Citations

All quotes verified by direct read. Line numbers refer to the file paths
listed in the table at the top of this document.

### Paper 1 (1911.08846)

- Title and authors: `/tmp/lit_fetch/1911.08846.txt:1-7`
- Abstract: `/tmp/lit_fetch/1911.08846.txt:8-11`
- 47-term Lagrangian (eq. 4): `/tmp/lit_fetch/1911.08846.txt:108-220`
- Form-factor definition (eq. 5): `/tmp/lit_fetch/1911.08846.txt:225-237`
- Local-limit Lagrangian (eq. 16): `/tmp/lit_fetch/1911.08846.txt:308-311`
- "parity breaking terms ... not considered": `/tmp/lit_fetch/1911.08846.txt:572-573`
- Ghost-freedom conditions (eq. 37): `/tmp/lit_fetch/1911.08846.txt:537-551`
- Sign conditions on C₂, C₄, C₅: `/tmp/lit_fetch/1911.08846.txt:566-569`
- Local-limit map (eq. C7): `/tmp/lit_fetch/1911.08846.txt:1015-1018` (in
  Appendix C; shown via grep, not a direct re-quote in this note)
- Reference [43] BJ-MT 1910.07506: `/tmp/lit_fetch/1911.08846.txt:1101`

### Paper 2 (2311.05597)

- Title and authors: `/tmp/lit_fetch/2311.05597.txt:1-8`
- Abstract: `/tmp/lit_fetch/2311.05597.txt:9-25`
- General teleparallel framework (eqs. 1-8): `/tmp/lit_fetch/2311.05597.txt:96-138`
- Field-content notation (eqs. 9-15): `/tmp/lit_fetch/2311.05597.txt:144-178`
- f(T) tensor-perturbation eq. (eq. 159): `/tmp/lit_fetch/2311.05597.txt:1325-1332`
- f(Q) tensor-perturbation eq. (eq. 164): `/tmp/lit_fetch/2311.05597.txt:1366-1373`
- Conclusion ("formally same as GR"): paper §VIII (lines around 1490-1510)
- Reference [37] BJHK 1710.03116: bibliography
- Reference [38] HHK 2311.05495: bibliography

### Paper 3 (2503.18972)

- Title and authors: `literature/2503.18972/f_R,Q_.tex:9-14`
- Abstract: `literature/2503.18972/f_R,Q_.tex:20`
- f(Q) pathologies summary: `literature/2503.18972/f_R,Q_.tex:27`
- f(R,Q) action: `literature/2503.18972/f_R,Q_.tex:191`
- MDBA Hamiltonian construction: `literature/2503.18972/f_R,Q_.tex:194-251`
- Reduced Hamiltonian (eq. AHF): `literature/2503.18972/f_R,Q_.tex:248`
- D'Ambrosio-Heisenberg-Zentarra 2308.02250 reference: `literature/2503.18972/f_R,Q_.tex:554`
- Tomonari-Bahamonde reference: `literature/2503.18972/f_R,Q_.tex:553`

### Cross-references

- FINAL_ASSESSMENT: `research/perturbative_hamiltonian/notes/FINAL_ASSESSMENT.md`
- Meta-N (BC 2018 Appendix D): `research/perturbative_hamiltonian/meta_reviews/meta_review_N_BC_AppendixD_content.md`
- MANUAL_RETRIEVAL_NEEDED: `research/perturbative_hamiltonian/MANUAL_RETRIEVAL_NEEDED.md`
