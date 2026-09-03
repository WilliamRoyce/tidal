# Spectrum module design — numerical polology for the cosmology program (H6 / WS6)

**Status:** design settled 2026-09-03 (H6). Tracking: #495 (WS6), #360 (canonical PSALTer
tracker), umbrella #488, decision D6.
**Scope:** design only — no code in this document's commit. Module implementation is WS6
proper; the Wolfram-side Stage-1 engineering is proposed as its own handoff (§13).

## TL;DR

**Problem.** The cosmology pipeline must screen out vacuum-sick theories (ghosts, tachyons)
before spending compute on them, at arbitrary coupling values, in milliseconds per point —
the regime where symbolic computation dies of expression swell (arXiv:2606.30785).

**Design.** Two stages. **Stage 1** (once per Lagrangian structure, Wolfram, unmodified
symbolic PSALTer): our auto-generated export turns the theory's vacuum Lagrangian into
PSALTer's input, and PSALTer's run *discovers the structure* — the wave-operator blocks, which
modes are massive/massless, gauge ranks, source constraints, and closed-form health conditions
where they exist. **Stage 2** (per coupling sample, JAX): a numerical evaluator consumes that
structure and returns a health verdict in order milliseconds. The primary ghost/tachyon
algorithm is the **Schur-complement kinetic-matrix criterion** of arXiv:2506.02111 — published,
worked, and implemented nowhere — which fits our theories (massless graviton/photon kinetically
mixed with massive torsion), avoids radicals and residues, and carries parity-violation support
with it. The residue route from the released code is implemented alongside it as a permanent
cross-check. The massless sector — which the released code guards out, and which our program
cannot ignore — is handled by counting (free), compiled Stage-1 closed forms, and a numerical
massless-residue fallback.

**Next step.** Dispatch the proposed Stage-1 handoff (§13): install PSALTer, stand up the
spectrum derivation branch and exporter. WS6 proper then implements Stage 2 against §12's
validation set.

**Done when** the seven validation groups of §12 pass — headline gates: the three
Lin–Hobson–Lasenby inequalities reproduced (the WS6 gate in `docs/COSMOLOGY_PROGRAM.md`),
agreement with the published implementations, the TorC and Einstein–Cartan oracles, GR+Maxwell
healthy, and a median per-sample verdict ≤ 1 ms.

---

## 1. The two-stage architecture

Numerical polology does **not** replace symbolic PSALTer; it consumes its output.
arXiv:2606.30785 §implementation-appendix ¶Platforms, verbatim:

> "An initial application of computer algebra is required to extract the symbolic
> `M_n^{J^P}` from `S(θ)`, and this problem was previously solved for any theory of the form
> given in [eq. EFTLag] by the Wolfram Language implementation in [Barker:2024juc,
> Barker:2025qmw]. The remainder of the pipeline is implemented in Python, specifically within
> the JAX framework."

Confirmed in the code: the published JAX entry point is `sample(wxf_path, …)`, and
`_extracting/wxf.py` reads `xAct`PSALTer`WaveOperator` out of a PSALTer-emitted WXF file. No
path in the numerical implementation accepts a Lagrangian.

| | **Stage 1** — once per Lagrangian *structure* | **Stage 2** — per coupling *sample* |
| --- | --- | --- |
| what | Lagrangian → wave-operator blocks, **plus structure discovery**: gauge ranks, the massive/massless mode partition, source-constraint matrix, massless content, closed-form conditions where obtainable | blocks + coupling vector → masses, ghost/tachyon verdict |
| where | Wolfram, symbolic PSALTer (unmodified — we generate its input and harvest its output, both from outside) | JAX, our module — adapted from the released code **and extended** with the unimplemented Schur criterion (§5) |
| cost | minutes–hours, amortized; one-wolframscript rule applies | **order milliseconds** (Barker's stated target) |

Stage 1 is not a necessary evil. It is where the structural facts that make Stage 2 both
correct and fast are discovered — which sectors are gauge-degenerate and by what rank, which
modes are massive vs massless (the Schur reordering), the constraint matrix, how many genuinely
massless polarizations exist. Those facts are what let Stage 2 use fixed shapes (§9) and what
the massless verdict asserts against (§6).

## 2. Purpose and scope

The module is a **screening gate**: theories/points that fail never reach the solver. Its
verdicts are **Minkowski-vacuum** statements, and that is sufficient by the same split TorC
used (spectrum screens in vacuo; dynamics evolve on FRW). Three honesty clauses, stated here
rather than footnoted:

- **One-way filter.** A vacuum-healthy verdict says nothing about behavior on a background
  (`B ≠ 0`, FRW); vacuum-sick is a reason to discard, vacuum-healthy is not a certificate.
- **Necessary, not sufficient.** Linear (tree-level, free-spectrum) health only; strong
  coupling and non-linear pathologies are out of scope.
- **Parity-odd status.** Symbolic PSALTer supports parity-violating theories today
  (arXiv:2506.02111 upgraded it); the published *numerical* work restricts to parity-preserving
  ("may be applied in future work", arXiv:2606.30785 §implementation fn. 9). Our chosen
  algorithm (§5) is where parity violation already lives, so support arrives as a consequence
  rather than as extra work — in scope, and what gates O4/#499's Chern–Simons couplings.

**Term-roster scope (user decision, 2026-09-03):** curvature-squared terms (`R̃²`-type, the
old `b5`) are **excluded** — their perturbation theory is problematic and we will not sample
them. Direct `(∇T)²` torsion kinetic terms (as in the thesis's T-series) are **in**, and are
degree 2 in momentum — inside the wave-operator contract. The `k⁴` degree assertion in the
code becomes a misconfiguration guard rather than a live limitation.

## 3. Reference hierarchy and provenance policy

**Citable sources** (everything load-bearing in this document stands on these):

- arXiv:2606.30785 — "Numerical polology" (Barker, Handley, Hobson, Lasenby, Marzo, Santoni,
  Torcellini). Local: `literature/2606.30785/`.
- arXiv:2406.09500 — PSALTer v1. Local: `literature/2406.09500/`.
- arXiv:2506.02111 — PSALTer v2, parity-violating; the Schur-complement criterion. Local:
  `literature/2506.02111/`.
- `github.com/wevbarker/SupplementalMaterials-2607` — the companion repository (verified live
  2026-09-03): `Julia/`, `JAX/`, `WolframLanguage/` including committed
  `ParticleSpectrograph<Name>.wxf` exports and `*_unitarity_NS` chains.
- arXiv:2507.09228 — TorC. Local: `literature/2507.09228/`.

`psalter.tar.gz` at the repo root is a WIP snapshot the author sent us; it is **cited nowhere**
here, but it is not off-limits: it carries pieces absent from the release (a JAX draft of the
`U`-likelihood, a per-sample classifier, cubed-sphere tile geometry, recursive sampling) that
the implementation may adopt where they are the best available, with provenance recorded.
Where it *disagrees* with the release, investigate the difference — it can be a fix, an
experiment, or a regression — with the release as the presumption, not an automatic winner.

**Provenance policy** (author permission is explicit — D6): heavy, deliberate inspiration and
copying from Barker's code is authorized; every borrowed function records in its docstring
what came from where, in the form:

```text
Adapted from SupplementalMaterials-2607 Julia/src/run/evaluator/unitarity.jl
(find_gauge_modes), arXiv:2606.30785 eqs. (BlockMatrix)-(RegularizedWaveOperator).
Changes: fixed-shape projectors decided per structure (see spectrum_design.md §9).
```

Formal attribution is settled with the supervisor at publication time. License: TIDAL is MIT
(`pyproject.toml:7`), the released `psalter` is GPL-3.0-or-later — recorded on #495 as a
**release** blocker (the trigger is distribution), with options ranked there: (a) relicense
our side GPL-3.0-or-later, (b) separately-licensed subpackage behind a process boundary,
(c) written permission. Design consequence carried regardless: the spectrum module's
dependencies are one-directional — the package calls into it, it calls back into nothing.

### 3.1 The two published implementations are the paper's two halves

They are not two ports of one thing; learn from both:

| | **Julia** (`Julia/`) | **JAX** (`JAX/src/psalter/`) |
| --- | --- | --- |
| implements | the paper's **tuned** case — survey coupling space for the unitary sub-volume of a theory whose particle content is known | the paper's **untuned** case — *discover* theories on measure-zero decoupling hypersurfaces |
| driver | `run.jl` + YAML; nested sampling on the `U`/compression unitarity likelihood | `sample(wxf_path, …)`: NS on the residual `Likelihood_Det` → LM refine → measure → reweight |
| evidence | examples `FierzPauliTheory`, `FierzPauliProcaTheory`, `TDiffTheory`, `VectorTheory` with `*_unitarity_NS` outputs — i.e. the paper's figs. `fp`/`fpp`/`tdiff`/`mv`; configs carry `nlive: 10000`, `scale: 10.0`, matching the stated `10⁴` live points and `λ = 10` | tests `A23Theory`, `S123Theory`, `TensorTheory`, `VectorTheory` with `*_spin{zero,one,two}_direct_sphere` outputs |
| carries | **gauge regularization** (`evaluator/unitarity.jl`), the **soft unitarity likelihood**, **mass dimensions**, **emergent symmetries**, phenomenology reweighters | the fast per-sample kernels: WXF ingest, Vandermonde det coefficients, companion-matrix roots, Kato residues, compile/jit/vmap idioms |

**"Tuned" and why we are squarely in it.** Tuned = the couplings already satisfy relations
producing a *known* particle content (Fierz–Pauli, Maxwell, Proca are the paper's examples);
the only question is where in coupling space that content is healthy. Our theories are exactly
this — PGT+EM's `J^P` content is settled literature (Sezgin–van Nieuwenhuizen, Karananas,
Blagojević, Lin–Hobson–Lasenby) — with one caveat and one resolution. Caveat: for a *new*
combination of torsion kinetic terms the content is not always known a priori (this project's
own history: the T-series kinetic sectors; the implicit-dynamical-sector findings of GH
issues 457/468). Resolution: **the symbolic PSALTer run tells us what propagates, and that puts us
in the tuned case** — Stage 1 computes the content per structure, the discovered content
becomes the declared baseline the per-sample verdict asserts against, and the survey asks
where parameter choices keep it healthy. The untuned machinery (decoupling residuals, LM
refinement) is never needed: we always hold one specific Lagrangian, not a search over theory
space. Content can still change *within* coupling space (a kinetic term degenerates, an
accidental symmetry appears, a mass hits zero) — those are flagged findings (§6, §8), never
scores; the Julia emergent-symmetries machinery is the purpose-built diagnostic if a theory
needs such loci mapped.

There is a second sense in which our theories are tuned, and it explains the code split: tuned
loci are *where gauge symmetries live*. That is why the tuned (Julia) implementation carries
`find_gauge_modes` and the regularized wave operator, and the untuned (JAX) one — whose
appendix assumes "no gauge symmetries" — does not.

**So `Julia/` is the reference implementation of our exact use case**, and the *shape* of its
per-sample flow is the one to adopt (`evaluator.jl`): `wave_operator(c)` → sector evaluation →
masses → **early return if not healthy** → and only then derived quantities. The verdict needs
only the front half: v1 ships the gate, and the tail — **mass dimensions** (needed the moment
anyone wants masses on a physical scale; the bare verdict never does) and **emergent
symmetries** — is *explicitly deferred*, with the extension point exactly where `evaluator.jl`
continues after its early return.

**References are algorithm sources, not language commitments.** The entire Stage-2 algorithm
is implemented in **JAX** (jit + vmap, fp64) — that is where the ms budget is won — with NumPy
only off the hot path. Nothing obliges a Julia runtime; algorithms carry over.

### 3.2 Reconciliation with H4's package design (`docs/cosmology/repo_reshape.md`)

H4 landed before this study and its `tidalcosmo/spectrum/` README anticipates revision by H6.
What H6 **confirms**: the module lives at `tidalcosmo/spectrum/`, on the sampling path,
feeding `validity/` and the Cobaya prior surface (H4 §6 — "WS6 is not independent"); it is
built at **M∥** in the migration order, wired in before the first observable rung; the verdict
is a flag with a reason through the shared flagged-rejection mechanism (H6's payload in §10
realizes that); per-sample speed is a requirement (D5). And **#513** — adopt CAMB and PSALTer
conventions natively, no conversion layer — is confirmed independently here with a hard
technical reason H4 did not have: the parity factor in every ghost verdict is *defined by* the
signature (§4.3). Adopted from H4 §2.8 without change: PSALTer's vocabulary (`Theory`,
`Sector`, `couplings`, `tiles`) wherever it means the same thing, and the frozen-dataclass
settings idiom for `tidalcosmo/config/`.

What H6 **revises** in H4 §2.8's PSALTer picture, which was written from the WIP tarball
before this study established the published record: the released JAX package's API differs
from the tarball's (no `evaluate`/`classify`/`tile` surface; the untuned driver instead), the
tuned pipeline — our use case — is the **Julia** implementation (§3.1), and the "emit WXF or a
schema PSALTer's `extract` can consume" alignment item is **superseded**: the WXF interchange
as it exists is insufficient for our theories (unlabeled `J`-blocks, placeholder slots — §6.1),
so our Stage-1 exporter emits the richer §6.1 contract, and Stage 2 grows the Schur criterion
that exists in no released code (§5). None of this re-litigates H4 — its own README invites
exactly this revision.

## 4. Stage 1: conventions and input contract — native, not converted

**Design principle:** the spectrum-facing symbolic branch is written *in the Barker ecosystem's
conventions from the start*, so its quadratic Lagrangian is PSALTer-ready by construction and
no conversion layer exists to maintain or get wrong. Legacy TIDAL conventions are not
inherited here (WS1: the new package is designed whole, from the goal backwards).

### 4.1 The input form is the ecosystem's definition of theory space

arXiv:2606.30785 eq. `EFTLag`: `S(θ) = ∫d⁴x Σᵢ θᵢ Oᵢ` — every operator carries its own
coupling, with mass dimensions *derived*, not assumed (eq. `MassRescaling`). arXiv:2506.02111
§TheoreticalDevelopment restates it as a wave-operator requirement: "a homogeneous, linear
parametrization … its every term must be linear in these couplings."

So "a symbolic coupling on every term" is not a PSALTer quirk to patch in — it is the
definition of the theory space the method operates on. Contract rules:

- **Reject, never auto-assign.** A Lagrangian term with a bare numeric coefficient (`−¼F²`)
  is a config error with a hint naming the term. Auto-assigning would silently invent a
  parameter the user never declared and then report verdicts and posteriors in terms of it.
  This matches PSALTer's own posture (it validates — `ValidateLagrangian.m`,
  `EnsureLinearInCouplings.m` — and the paper leaves even `α = ½` free "because PSALTer
  requires that all terms … be linearly parameterised"); confirm the exact failure behavior on
  a live install (§13).
- **Fields**: real tensors, rank ≤ 3, from the enumerated symmetry list (identical in v1 and
  v2): scalar; vector; rank-2 general/`[μν]`/`(μν)`; rank-3 general plus eight symmetrizations
  including `ζ_{μ[νσ]}`. **Torsion `T^λ{}_{[μν]}` is exactly that class** (PSALTer's
  `RegisterFieldRank3Antisymmetric23`). Validate at config load, `error_with_hint` otherwise.
- **Sources are PSALTer's business.** `DefField` defines each field's conjugate source and
  `ParticleSpectrum` appends the source coupling itself; we supply only the quadratic field
  part and assert nothing linear in the fields survives. Sources reach us only as *output we
  read*: `ComponentSourceConstraints` is the constraint matrix, and the massless residue
  eigenvalues live on the reduced source space (§6).
- **No distinguished `κ`.** Neither reference paper has a κ convention to match — polology's
  couplings are generic `θᵢ`; **TorC has no Einstein–Hilbert term at all**
  (`literature/2507.09228/paper_Qtorsion.tex:171-180`: the EH term is emergent from
  condensation; its conditions `λ ≥ 0`, `μ < 0`, `(ν+2μ)(ν−μ) > 0` involve no κ). Give the EH
  operator its own generic coupling like every other operator; the mass-dimension machinery
  then works as designed. Nothing in the numerical code assumes a κ convention.

### 4.2 Field normalization: what matters and what does not

- The **sign** of `h` (`g − η` vs `η − g`) is verdict-irrelevant: `h → −h` is an invertible
  field redefinition; the quadratic action transforms by congruence, preserving pole positions
  and residue signs.
- **Positive rescalings** (`h → λh`, a `κ`/`M_p` absorbed into the field) are also
  verdict-irrelevant (residues scale by `λ² > 0`).
- Both **matter for oracle comparisons** at the matrix level (the overall factor in §5's
  published `K₀` is convention-laden) and for mapping couplings onto literature conditions.

**Rule: bare `h_{μν} ≡ g_{μν} − η_{μν}`, plus sign, no rescaling, throughout the spectrum
branch** — matching TorC's post-Riemannian form (and its tetrad analogue, the Kronecker gauge
`e_i{}^μ ≡ δ_i{}^μ + f_i{}^μ`, `paper_Qtorsion.tex:588`). xPert's `DefMetricPerturbation`
already produces exactly this bare `h`. PSALTer itself imposes nothing — `h` is whatever field
the quadratic Lagrangian is written in. The CAMB side has its own `h_{ij}` conventions
(synchronous gauge, Ma–Bertschinger), and the solver branch is free to adopt them — per-branch
native applies to field conventions too.

### 4.3 Signature: two ecosystems, opposite conventions, one answer per branch

| | signature | `ε` | source |
| --- | --- | --- | --- |
| Barker ecosystem (PSALTer v1/v2, polology) | `(+,−,−,−)` | `ε₀₁₂₃ = +1` | 2606.30785:364 ("the particle physics signature"); 2506.02111:100; 2406.09500 ¶Conventions |
| CAMB/CLASS ecosystem (Ma–Bertschinger) | `(−,+,+,+)` | — | `literature/astro-ph_9506072/9506072.tex:295`: `ds² = a²(τ){−dτ² + (δᵢⱼ + hᵢⱼ)dx^i dx^j}` |
| legacy TIDAL | `(−,+,+,+)` | `ε₀₁₂… = −1` | `tidal/wolfram/CommonUtilities.wl:31-33` |

"Adopt the other package's standard" therefore has two different answers, and that is fine:
**spectrum branch → Barker-native; solver branch → CAMB-native** (where legacy TIDAL already
sits). Signature and `ε` are declared, machine-checked properties of each derivation, never a
global project constant. The two branches never exchange tensors (§4.5), so nothing is
converted anywhere.

Where the signature actually enters the spectrum branch (the content of "native"): the branch
never evaluates components — PSALTer does its own SPO decomposition — but three things are
signature-tied even at abstract-index level: (1) the `ε·ε` contraction identities (they carry
the sign of `det g`), so parity-odd terms expand differently; (2) the sign conventions of the
curvature/torsion definitions and of the EH-type operator; (3) any
`ToCanonical`/`ContractMetric` simplification, which uses the session's declared metric. Being
native costs nothing: the emitted script is standalone and PSALTer pre-defines its own
`M4`/`G`/`epsilonG`/`CD` (reserved indices `a`–`z`) and forbids tampering — the branch simply
works in that geometry from its first line.

**The hard reason a conversion layer is forbidden:** arXiv:2506.02111 (SPO normalization
appendix) fixes the parity of a state by `P ≡ sgn(ζ*·ζ)` **"due to our choice of signature"**,
and the appendix's Hermitization step is exactly multiplication by `diag(1, −1)` on the parity
blocks. The parity factor that multiplies every residue — `parity_sign` in the released code —
is therefore *defined by* the signature; a conversion layer is precisely where a ghost verdict
would silently flip sign. `docs/tex/pgt_stability_priors.tex` §"Literature cross-check" (the
worked signature/Proca reconciliation for the torsion sector) changes role accordingly: from
conversion recipe to the one-time reconciliation against legacy results (validation (e)).

### 4.4 Covariant derivatives are expanded away before handoff

PSALTer's `CD` is the flat derivative — `CD[-m]@` **is** `∂_μ` (2406.09500 ¶"Loading the
software"). Our Lagrangians use `CDT`, the torsionful Riemann–Cartan connection. The spectrum
branch therefore reduces everything to `η`, `∂` and the perturbation fields: post-Riemannian
decomposition (connection → Levi-Civita + contortion, contortion → torsion via
`K^a{}_{bc} = ½(T^a{}_{bc} + T_b{}^a{}_c − T_{bc}{}^a)`), then quadratic expansion about
`(η, T̄ = 0)`. Christoffel and torsion couplings appearing explicitly is the point, not a side
effect. The design of this stage exists in legacy TIDAL
(`tidal/cli/_derive.py::_wls_torsion_curvature_decomposition`, line 2088 — `ChangeCurvature[L,
CDT, CD]` plus the contortion identity) and is ported, not imported (WS1's no-adapters rule).

This lands the analysis in the **post-Riemannian / Einstein–Cartan formulation** — fields
`(h_{μν}, T^λ{}_{μν}, a_μ)` — matching TorC's `ParticleSpectrographTorCECT` (15 gauge
generators), not the tetrad + spin-connection `…TorCPGT` (21). That choice fixes which
published spectrograph is the primary oracle (§12c).

### 4.5 Two derivation branches, not one derivation with a fork

| | **spectrum branch** | **solver branch** |
| --- | --- | --- |
| background | Minkowski vacuum | FRW with `a(η)` |
| background fields | absent from its input (separate solver-only config block) | kept and central (the `B` field is the physics) |
| product | quadratic **Lagrangian** | linearized **EOM** |
| representation | abstract index, never evaluated to components | components |
| conventions | Barker `(+,−,−,−)`, `ε₀₁₂₃ = +1` | CAMB `(−,+,+,+)`, conformal `τ` |

They differ in every row: **two derivations from one input**. The solver branch does not
detour through the PSALTer-shaped quadratic Lagrangian — it goes straight to the EOM, which is
faster and what it needs.

**Shared: everything convention-free.** Config parse and validation; field declarations and
symmetry-class classification; the coupling roster and its bookkeeping; the term-level
Lagrangian structure; gauge-symmetry declarations; scope guards. That shared input model is
decided once and consumed by both branches. **Duplicated: only the convention-laden symbolic
work** — the post-Riemannian rewrite (it carries curvature-sign conventions) and everything
downstream. Run it inside each branch's own session; measure later whether the duplication is
worth optimizing.

**Input architecture (fits H4's config design).** The config separates into a **Lagrangian
block** in ecosystem-aligned vacuum covariant form — the shared input model — and a
**background block** consumed only by the solver branch. This slots into H4's `tidalcosmo/config/`
typed-config layout (`docs/cosmology/repo_reshape.md` §2.8, §2.11); designing the background
block's content is WS2's job within that layout, and it is genuinely richer on FRW than the
legacy uniform-`B₀`-on-Minkowski setup (time-dependent, possibly stochastic fields). Related
legacy loss: the plane-wave 1-D reduction speedup (`∂_x = ∂_y = 0`) is unlikely to survive FRW
with structured backgrounds — assume full derivations.

**Gauge treatment is necessarily different per branch — by requirement.** The spectrum branch
hands PSALTer the **gauge-unfixed** quadratic Lagrangian: finding the gauge symmetries and
imposing source constraints is PSALTer's own job. The solver branch applies its own gauge
choices (CAMB-side conventions; the evolved `[[gauge]]` machinery). What must agree is the
**gauge-invariant content**: spectrum verdicts are gauge-invariant by construction, and
cross-branch comparisons compare gauge-invariant quantities (effective signs, invariant
combinations — the existing comparison rules in `.claude/rules/wolfram.md`); the solver side's
`gauge_certificate` culture already flags gauge-dependent numbers there.

**Spectrum-branch order:** (1) consume the vacuum Lagrangian block only, asserting no
background quantity leaks in; (2) `ChangeCurvature` + contortion→torsion; (3)
`Perturbation[·,2]` + `ExpandPerturbation` about `(η, T̄ = 0)` → the quadratic Lagrangian, in
Barker-native conventions from the start; (4) emit `DefField`s from the validated field specs
and `DefConstantSymbol` per coupling, pre-expand derived fields (`F = dA`), assert nothing
linear in the fields; (5) emit a **standalone** `.wls` running `ParticleSpectrum`
(`Method→"Hard"`, `MaxLaurentDepth→1`) plus our own export (§6.1). Model the generator on
`tidal/cli/_derive.py::{_wls_fields,_wls_lagrangian}` + `tidal/cli/_wls_helpers.py`; model the
emitted script on the published `WolframLanguage/ParticleSpectroscopy/{FieldKinematics.m,
Models/*.m}`.

## 5. The algorithm choice: Schur-complement kinetic matrix, residues as cross-check

arXiv:2506.02111 §MassiveSpectrum proposes a no-ghost criterion that "does not involve any
inversion of the wave operator nor the computation of residues of the propagator at massive
poles". It is fully laid out — derivation in its `NoGhostCriterionAppendix`, two worked
examples — and **implemented nowhere**: "our new method is not yet implemented in PSALTer, but
we illustrate its use" (§ECSpectroscopy). Per `J` sector:

1. **Reorder** rows/columns of the sector coefficient matrix into
   `Õ_J = [[O_VV, O_VL], [O_LV, O_LL]]` — massive modes in `O_VV`, massless modes in `O_LL`,
   their kinetic mixings off-diagonal (eq. `canonical_block_structure`).
2. **Block-diagonalize** by the unit-triangular congruence (eq. `tilde_OJ`); the **Schur
   complement** `Ô_VV = O_VV − O_VL O_LL⁻¹ O_LV` "encodes all the information about the
   massive spectrum" and decomposes as `Ô_VV = K_J k² + M_J`, with `K_J`, `M_J` functions of
   the couplings only (eq. `DiagonalizedWaveOperator`).
3. **No ghost ⟺ `K_J` negative-definite** (eq. `kin_matrix_negative`).
4. **Masses** = eigenvalues of `−K_J⁻¹ M_J`; no-tachyon = all real-positive, testable for a
   2×2 via `tr² − 4 det > 0 ∧ tr > 0 ∧ det > 0` (eq. `mass_constraints`).

**Why it is the primary algorithm here:**

- **It targets our exact difficulty.** The paper: the residue criterion "becomes extremely
  difficult … the complications are due to the presence of massless modes and — more
  specifically — because of their kinetic mixings with the massive ones", with
  **Einstein–Cartan gravity as the worked example**. PGT+EM — massless graviton and photon
  kinetically mixed with massive torsion — is that hard case as our default case.
- **It is cheaper and better conditioned** — the ms budget's best friend: per sector, one
  small inverse, two products, a definiteness test (Cholesky), one small eigenproblem. It
  replaces det-coefficient recovery → polynomial roots → per-root SVD → per-root residue, and
  it **eliminates the Abel–Ruffini problem outright** (masses are eigenvalues, never
  closed-form roots — the polology paper flags exactly this failure of closed forms).
- **It is where parity violation lives.** The blocks are chequer-Hermitian
  (2506.02111:234: Hermitian diagonal blocks, skew-Hermitian off-diagonal, states ordered
  `P=+1` then `P=−1` within each `J`), and the generalized residue criterion
  `NoGhostCriterionMixed` — `Res(tr Ô⁺_{J,+} − tr Ô⁺_{J,−}) > 0` — is derived in the same
  appendix. One implementation delivers the ghost verdict, the massive–massless mixing
  treatment, and parity-odd support together.
- **It is genuinely unimplemented** — a real contribution, plausibly the "finishing the
  numerical package so it can handle these theories" the supervisor described. **Planned and
  attempted, not proposed-and-gated**; tell the supervisor (it is a contribution to his
  package, and he may have work in flight).

**Caveats, all checkable:** `O_LL` must be invertible, so gauge modes are removed first
(PSALTer's Moore–Penrose gauge fixing does this) — the method *relocates* the gauge problem
into Stage 1 rather than removing it (§8), which is precisely what shrinks its per-sample
burden. The massive/massless reordering is structural — a Stage-1 discovery, fixed-shape per
structure (§9). Chequer-Hermitian bookkeeping is a genuine trap: PSALTer's convention for the
parity-violating SPOs differs from Karananas's by a factor of `i` (`app:PGT_comparison`;
either is valid), and getting it wrong flips ghost verdicts — the published `K₀` below must be
reproduced including its `i`'s.

**Cross-check discipline (a design principle, not just a validation item):** implement the
residue route alongside, never as a second production path — the shipped verdict and soft
likelihood always come from the Schur quantities. The two criteria are proven equivalent, so
any disagreement localizes a bug immediately; and our Schur implementation can additionally be
checked against the *published* residue implementations (JAX kernels; Julia
`evaluate_sectors_with_regularization`), giving an external oracle to an algorithm that
otherwise has none.

### 5.1 Published worked examples — the validation oracles

**`ScalarParityViolatingPGT`** (2506.02111 §ECSpectroscopy): the most general parity-indefinite
Einstein–Cartan action propagating the graviton plus two scalars,

```text
S = ∫d⁴x e [ c₁R + c₂H + c₃R² + c₄RH + c₅H² + c₆T_{μνρ}T^{μνρ} + c₇T_{μνρ}T^{νρμ}
             + c₈T_μT^μ + c₉ε^{μνρσ}T_{λμν}T^λ_{ρσ} + c₁₀ε^{μνρσ}T_{μνλ}T_{ρσ}^λ ]
```

(`H` the Holst pseudoscalar). Published targets, in increasing order of what they test:

| target | value | tests |
| --- | --- | --- |
| massless graviton no-ghost | `c₁ < 0` (eq. `massless_noghost`) | the massless path on EC field content |
| the four `J = 0` blocks | printed in eq. `tilde_OJ_actual` | reordering and block identification — a checkable *intermediate* step |
| kinetic matrix | `K₀ = 6·[[−4c₅, ic₄],[ic₄, c₃]]` | Schur complement + chequer `i` convention |
| mass matrix | `M₀` in eq. `eq:kin_mass_matrices` (Υᵢ abbreviations from the spectrograph; `Υ₁+Υ₃ = 2c₁`) | same, plus abbreviation handling |
| no-ghost verdict | `c₅ > 0 ∧ 4c₃c₅ − c₄² > 0` (⟹ `c₃ > 0`) | definiteness test end to end |
| no-tachyon verdict | eq. `no_tachyon_cond_3`, **fully explicit in bare couplings**: `(2c₆−c₇+3c₈)[(2c₁−2c₆+c₇−3c₈)(c₁−4(c₆+c₇)) + 8(c₂+c₁₀−2c₉)²] > 0 ∧ (2c₆−c₇+3c₈)(2c₁−2c₆+c₇−3c₈) < 0` — "the well-known constraints [Blagojević]" | the `−K⁻¹M` eigenvalue test end to end |
| gauge generators | **10** (Poincaré redundancy; no accidental symmetries) | Stage-1 gauge-rank discovery |

Intermediate matrices being published — not just final inequalities — is what makes this a
usable oracle for an algorithm with no reference implementation.

Two more facts from the same figure caption that shape this design: PSALTer **halts** its
symbolic analysis "once it has been determined that the square masses are not rational
functions of the Lagrangian coupling coefficients", so its exported closed-form conditions
resolve *only the graviton* for this theory — i.e. **for exactly our theory class, L0 (§6) is
partial and the numerical massive path is not optional**; and the improved algorithm "not yet
implemented in PSALTer" recovers the rest — which is what we implement.

**Second oracle:** `app:PGT_comparison` — general parity-violating PGT; masses "identical to
[Karananas:2014pxa]", no-tachyon conditions matching Blagojević, no-ghost "almost trivially by
inspection" via the same method; and the spin-one sector contains **no massless propagating
particles** (gauge redundancies) — a counting check for §6.

## 6. What Stage 1 hands over, and the massless sector

### 6.1 The export

The committed `ParticleSpectrograph<Name>.wxf` files (decoded 2026-09-03 with a minimal WXF
reader) show what PSALTer can emit and what our own exporter must add:

- `…TDiffTheory.wxf` (gauge theory) has **populated**: `WaveOperator` (the sector blocks);
  `PseudoDeterminant`; **`UnitarityConditions`** in closed form
  (`Theta1<0 ∧ Theta3<0 ∧ …`); `ComponentSourceConstraints` — the constraint matrix of
  2406.09500 eq. `MatrixEquation` in explicit components with `En`, `Mo` and
  `SourceRank0/1/2Symmetric` heads, **coupling-independent** (structural gauge symmetries);
  `ComponentSaturatedPropagator` (a 646 KB rational expression).
- `…TensorTheory.wxf`: `UnitarityConditions = Text["(Demonstrably impossible)"]`;
  `ComponentSourceConstraints` empty. Both cases must be handled.
- `…VectorTheory.wxf` / `…A23Theory.wxf`: several slots are *unevaluated placeholder symbols*,
  and key names vary (plural vs singular) between exports → read the association keys and
  normalize; never hardcode; treat a bare symbol as absent. The two small files are ready-made
  unit-test fixtures.

**Block labeling — settled by decoding `…A23Theory.wxf`.** The A23 field (rank-3, antisym in
the last pair — torsion's class, 24 components, `J^P` content `0±,1±(×2),2±`) exports
**three** `WaveOperator` entries of dims 2, 4, 2: `2·1 + 4·3 + 2·5 = 24` ✓. These are
**`J`-blocks containing both parities**, with the `J = 1` block visibly `2+2` block-diagonal
(parity decoupled, as it must be for a parity-preserving theory) — 2506.02111's ordering
convention, observed in the wild. The released JAX `SPIN_LABELS = ["0+","1-","2+","3-",…]`
positional convention is therefore wrong even for the release's own test theory once a field
carries both parities. **Our export carries explicit per-state `J^P` labels; position is never
a label.** (Parity feeds the Hermitization `diag(1,−1)` and hence every ghost verdict — §4.3.)

**Exported form: coupling-linear coefficient tensors, not symbolic matrices.** PSALTer
guarantees linearity in the couplings, so each block is `A[n,0] + Σᵢ cᵢ·A[n,i]` with `A`
numeric — no sympy/lambdify on the per-sample path, exact and free `∂M/∂c`, expression swell
kept Wolfram-side. The Schur blocks are fixed permutations of the same tensors; one export
serves both routes. `z_degree` (residue path) is exported as Stage-1 metadata or probed at
random coupling points, with the non-generic-point failure mode documented.

### 6.2 The massless sector

The paper's license for skipping it (arXiv:2606.30785 §implementation ¶"Mass spectra"): the
massless case needs "a full component decomposition" (little-group change), and "additional
radiative degrees of freedom are excluded by the thermal history" — footnoted as holding only
"for populated species in thermal equilibrium". Sound for their models (their targets are new
*massive* particles; `TDiff` is massive gravity with no massless graviton at all). **Inverted
for us**: the massless graviton and photon are the carriers of every observable this program
targets, their health is coupling-dependent (TorC's `(ν+2μ)(ν−μ) > 0` *is* a massless residue
condition), and the released health guards (`|Re m| > 0.01` in JAX, `abs(k) > TOL` in Julia)
would reject GR+Maxwell — they are scope guards implementing the license, not health criteria,
and no massless pole may ever reach them.

Four layers, by cost:

- **L1 — massive spectrum, per sample**: primary = Schur (§5); residue route as cross-check.
- **L2 — massless counting, essentially free, per sample**: under Schur, the structurally
  massless modes are the `O_LL` block by construction, and a massive mode *becoming* massless
  is a zero eigenvalue of `−K_J⁻¹M_J`; under the residue route, the multiplicity of `z = 0`
  sits in the det coefficients (`c₀ = det M₀`). Compared against the Stage-1-declared content.
- **L3 — massless health, compiled closed form**: Stage 1's massless analysis
  (`MaxLaurentDepth→1`) exports `UnitarityConditions` including the massless residue
  conditions (TorC's third condition; EC's `c₁ < 0`); compile to a per-sample boolean. Partial
  for our theory class (§5.1) — present for the graviton-type conditions, absent where PSALTer
  halts.
- **L4 — numerical massless residue, fallback and cross-check** (2406.09500 §MasslessSpectrum,
  eqs. `MatrixEquation`→`TotalNoGhost`, implemented numerically nowhere): constraint matrix
  `C(En, Mo)` → right null vectors with the `(E−p)^{μ_b}` spurious-pole guard → reduced
  Hermitian `Π_red = ξ†O⁺ξ` (Moore–Penrose by SVD — expensive symbolically, cheap
  numerically) → residue at `E→p` by extrapolation in `E` at fixed `p` (the direct analogue of
  the released Vandermonde det-coefficient trick) → Hermitian eigendecomposition →
  `λ_c > 0 ∀c`, polarization count = number of non-zero `λ_c`. The paper notes the `λ_c` are
  rational in `p²`, which obfuscates *symbolic* conditions; numerically it does not — evaluate
  at 2–3 values of `p` and require sign agreement.

L3 when present; L4 is what keeps the module general (closed forms provably do not always
exist — Abel–Ruffini, and `TensorTheory`'s "(Demonstrably impossible)"). L3 ↔ L4 agreement is
a validation gate. **Extra** massless poles beyond the declared content are *flagged and
routed to the program's `ΔN_eff ≲ 0.1` validity check* (`docs/COSMOLOGY_PROGRAM.md` §"Validity
enforcement"; Domcke & Garcia-Cely arXiv:2006.01161, local), not auto-rejected — footnote 8's
thermalization caveat is precisely the loophole a weakly-coupled spectator sector lives in. A
count *below* declared content means the structure changed and compiled L3 conditions no
longer apply: a loud finding, never a score.

## 7. Stage 2, step by step, with sources

`pub` = `SupplementalMaterials-2607/JAX/src/psalter/`; `jl` = `Julia/src/run/evaluator/`.
Rows 2–8 are the **residue route** — the cross-check path; the primary path (reorder → Schur →
`K_J` definiteness → `−K_J⁻¹M_J` eigenvalues) is §5 and has no published code to cite. Rows 1
and 9–12 apply to both routes.

| # | Step | Paper (arXiv:2606.30785 unless noted) | Code |
| --- | --- | --- | --- |
| 1 | `M(k) = M0 + k·M1 + k²·M2` per sector, Hermitian (chequer-Hermitian if parity-odd), linear in couplings | §implementation ¶Symmetries; 2406.09500 eq. `BuildWaveOperator` | Stage-1 output |
| 2 | Split `Def` powers; assert degree ≤ 2 | — | `pub _extracting/wxf.py::split_def_powers` |
| 3 | `z ≡ k²` reduction; `det M(k)` even in `k` | eq. `DetExpansion` | `pub wxf.py::_z_degree` |
| 4 | `det` at nodes `z_i = i`; precomputed Vandermonde inverse | eqs. `DetNodes`, `Vandermonde` | `pub _extracting/{kernels,coeffs}.py` |
| 5 | Roots in `z`: closed form deg ≤ 2, Frobenius companion ≥ 3 | ¶Measuring | `pub _extracting/masses.py` |
| 6 | Residue `Z = P·Re(v†(M1 + 2\|m\|M2)v)`, `v` from SVD of `M(k=m)` | eqs. `ResidueDefinition`, `NoGhost` | `pub _extracting/residues.py` |
| 7 | Health: `m` real ∧ `Z > 0` (∧ the `\|Re m\| > 0.01` scope guard) | eqs. `NoTachyon`, `MassiveNoGhost` | `pub _measuring/likelihood.py`; `pub measuring.py` |
| 8 | Smooth unitarity log-likelihood | eqs. `Likelihood`, `Compression`, `RelevanceKernel` | `jl unitarity.jl` |
| 9 | Gauge handling: Toeplitz block SVD → polynomial null modes `V(k)`; regularized `M̃ = M + VV†`; projection; spurious-root stripping | eqs. `PolyExpansions`, `BlockMatrix`, `RegularizedWaveOperator`, `GaugeProjection` | `jl unitarity.jl::{find_gauge_modes, regularize_wave_operator, project_out_gauge, is_gauge_root, solve_pep}` |
| 10 | Mass dimensions | eq. `MassRescaling` | `jl mass_dimensions.jl` — deferred (§3.1) |
| 11 | Gauge-surface proximity / emergent symmetries | eq. `SVHellmannFeynman` | `jl broken_symmetries.jl` — optional diagnostic, never on the verdict path |
| 12 | Massless spectrum | 2406.09500 §MasslessSpectrum | implemented numerically nowhere — L4, §6.2 |

## 8. Why the released code drops gauge sectors — worked from the published files

`_z_degree` computes `det M(k)` symbolically over the coupling symbols; returning `-1` means
the determinant vanishes identically in `k` *and* the couplings — a gauge symmetry acting
inside the sector — and the evaluator skips the sector. Decoding the published
`ParticleSpectrographVectorTheory.wxf` (general vector theory
`L = −½Θ₁∂_aA_b∂^aA^b + ½Θ₂(∂·A)² − ½Θ₃A²`) gives two 1×1 blocks, verbatim:

```text
0+ :  ½[ k²(Θ₂ − Θ₁) − Θ₃ ]
1- :  ½[ −k²Θ₁ − Θ₃ ]
```

Maxwell is the locus `Θ₁ = Θ₂, Θ₃ = 0` (that combination is `−¼F²` up to a total derivative):
the `0⁺` block vanishes identically → sector skipped; the `1⁻` root moves to `z = 0` → a
massless pole the scope guard discards. **Maxwell hits both failure modes at once.** The paper
itself analyzes this point — as a *destination*: its technical-naturalness trajectories
"terminate in the massless limit Θ₃ = 0 … the emergent symmetry is that of Maxwell theory",
where "two out of three vector polarizations survive" (2606.30785:515 and footnote). For them
a measure-zero boundary reached by steepest descent; **we live on that surface**.

Answers to the natural objections:

- *Is skipping harmful for pure EM?* No — Maxwell's `0⁺` sector is entirely pure gauge;
  skipping it is correct. The damage is confined to **mixed** sectors holding a gauge
  direction *and* physical states — and in PGT+EM the `0⁺` block is indexed by all fields'
  `0⁺` states at once, so one U(1) null direction throws away every torsion scalar in the
  block. TorC counts 15 (ECT) to 21 (PGT) gauge generators; scalar-EC counts 10.
- *Do new couplings fix it?* Not if they are gauge-invariant — built from `F`, curvature and
  torsion covariants, they enforce the degeneracy structurally at every coupling value. Only a
  symmetry-*breaking* term (a bare `A_μ`, a photon mass) lifts it — the deferred `cs1`–`cs3`
  case, where the newly visible mode's health is exactly what wants checking.
- *Both pipelines?* Split by issue. **Gauge degeneracy**: a problem only for the untuned JAX
  code; the tuned Julia code handles it (row 9), and the Schur route relocates it to Stage 1.
  **Massless health**: computed numerically by *neither* — both guard massless poles out; that
  is the genuinely missing piece §6.2 supplies.
- *So the problem is inherently that our theories are tuned?* **Yes — that is the one-sentence
  version.** The released fast path assumes genericity; gauge-invariant theories are maximally
  non-generic; §§4–6 are the response to that single fact.

Two failure modes, both handled: **structural** degeneracy (identical in couplings — fixed
rank, fixed-shape projectors, §9) and **accidental** degeneracy at special coupling points
(numerically singular only — detect and flag, never silently absorb).

### 8.1 Deforming off the symmetric locus — a capability, not the default

Since §4.1 gives every operator its own coupling, the symmetric locus is one point of our
parameter space, and the framework explicitly endorses expanding a symmetric model into the
surrounding non-symmetric space — one-way, with small, technically natural breaking couplings
(2606.30785:343, 345). A deformed theory (`Θ₃ ≠ 0`) is generic: the released residue code
works as-is and the massless problem softens. It cannot *replace* proper gauge handling, for
two concrete reasons: the symmetric limit is **discontinuous and branch-dependent** (the
paper's own footnote: the transverse branch's `Θ₃ = Θ₁ = 0` limit "propagates *nothing*",
unlike the Proca branch's surviving two polarizations; for gravity this is vDVZ), so a
deformed verdict is not automatically a verdict about the gauge-invariant theory; and it
trades a structural zero for an **ill-conditioned near-zero** — tolerance-based rank
identification, the ill-posed regime documented in this project's #473/#477 arc. Use
deformation as: a continuity diagnostic (sweep the breaking coupling → 0 against the
gauge-handled answer), a discontinuity detector (a vDVZ-type finding, not a bug), and a
bootstrap (early verdicts flagged `deformed`, never presented as verdicts about the symmetric
theory). Deformation changes the *theory*; Schur changes the *algorithm* — independent
mitigations, and deformation keeps its diagnostic value under either.

## 9. Performance

Anchors: the polology paper's tuned surveys — the corner-plot figures — ran at "a walltime of
minutes on a single CPU core at double precision" with `10⁴` live points
(§implementation-appendix ¶Platforms footnote), putting per-evaluation at µs–ms; the released
`sample.py` reports refine cost in ms/row. Those are small theories — the figure anchors the
scale, not a promise for PGT+EM-sized blocks.

- **Compile-once, per structure**: coefficient tensors; the reordering permutation and block
  shapes (Schur); det kernels and root finders (residue route); gauge null-space structure;
  constraint null vectors; JIT compile.
- **Per sample, primary (Schur)**: tensordot for the blocks; per sector one small inverse, two
  products, a Cholesky definiteness test, one small eigenproblem.
- **Per sample, residue cross-check**: det evaluations + Vandermonde solve + roots + one SVD
  per root — on demand, not in the survey hot loop.
- **Per sample, L4**: one pseudoinverse of the component-basis operator (~38–44² for PGT+EM)
  at a few `E` nodes.

**The rule that protects the budget: every rank and shape decision is made once per structure,
never per sample.** The Julia reference must not be copied verbatim here — it calls
`find_gauge_modes` per sample with tolerance-based rank determination, data-dependent shapes
that `jax.jit`/`vmap` cannot take. Structural gauge symmetries are coupling-independent
(verified: TDiff's exported source constraints contain no couplings), so the rank is fixed at
a generic coupling point at compile time; accidental degeneracies are detected and flagged.
The Julia `evaluator.jl` early-return shape (nothing computed for sick samples beyond the
verdict) is adopted; its internals are not.

**Threshold, written down before coding:** median per-sample verdict **≤ 1 ms** on the
`torsion_gertsenshtein`-family sector matrices (roster per §2) at fp64 on one CPU core;
compile time reported separately.

## 10. Design decisions (consolidated)

1. **Two-stage split** (§1); Stage 1 unmodified PSALTer driven from outside; Stage 2 all-JAX.
2. **Ecosystem-native conventions per branch**; no conversion layer anywhere (§4).
3. **Coupling-per-operator input contract; reject, never auto-assign** (§4.1).
4. **Bare `h = g − η`**, no rescaling, spectrum branch (§4.2).
5. **Primary ghost/tachyon algorithm: Schur-complement kinetic matrix; residue route as
   permanent cross-check** (§5). Redundant computation is a design principle wherever two
   published algorithms exist for one quantity (Schur/residue, L3/L4, deformed/gauge-handled)
   — the only affordable substitute for a reference implementation that mostly does not exist.
6. **One evaluator, two outputs**: hard verdict and smooth `log L` from the same Schur
   quantities — the manuscript's `U`/compression applied to the eigenvalues of `−K_J` (ghost
   margin) and `−K_J⁻¹M_J` (tachyon margin) — never `U` wired to residues while the verdict
   reads `K_J` (`feedback_shared_stage_not_copied_stage`).
7. **Explicit `J^P` labels per state**; position is never a label (§6.1).
8. **Massless layers L2/L3/L4** with declared-content assertion; extra massless poles routed
   to the `ΔN_eff` flag (§6.2).
9. **Fixed shapes per structure** (§9).
10. **Gate-only v1**: mass dimensions and emergent symmetries explicitly deferred, extension
    point named (§3.1).
11. **Verdict payload** follows the `gauge_certificate` honest-flag convention: per-state
    `J^P`, masses, `K_J` spectra, health, massless count vs declared content,
    `massless_sector` provenance (`compiled`/`numerical`/`not-assessed`),
    `extra_massless_poles` (→ `ΔN_eff` flag), `ghost_method` (`schur-kinetic`/`residue`),
    parity-odd support state, `k4_guard`, tolerances.
12. **One-directional dependency**; provenance docstrings; license handled as a release
    blocker on #495 (§3).

**Rework list against the released code:** (1) gauge sectors silently dropped (§8) — large on
the residue route, relocated to Stage 1 on the Schur route; (2) massless poles scored as sick
(§6.2); (3) positional spin labels (§6.1) — now demonstrated wrong on the release's own
`A23Theory` fixture; (4) the `k⁴` degree assertion — kept as a guard; the roster excludes
curvature-squared terms (§2).

## 11. Integration points

- **Screen before compute**: sick points never reach the solver.
- **`run_status` tagging in inference**: this is #360's Approach D, recorded in
  `docs/V3_2_DESIGN_INVESTIGATION.md:36` as "the ultimate goal", deferred only because
  "PSALTer integration is non-trivial". This workstream expires that deferral.
- **Soft prior designed in, not wired now**: the smooth `log L` (decision 6) composes with the
  joint-prior surface (`tidal/inference/_prior.py::RadialAngularPrior` / `parse_joint_prior`),
  which D9 schedules for re-implementation as a Cobaya component. Under a hard gate,
  positivity is imposed by fiat and invisible in the posterior; under the soft/tagged form it
  emerges as a learned constraint, and `D_KL` quantifies how much the data vs the unitarity
  likelihood contributes — a result, not plumbing.
- **Validity flags**: `extra_massless_poles` feeds the program's spectator-validity
  enforcement (`ΔN_eff`); the verdict payload joins `gauge_certificate` in the honest-flags
  family.

## 12. Validation set

| | check | tolerance / criterion |
| --- | --- | --- |
| (a) | **Lin–Hobson–Lasenby / TIDAL inequalities** — the WS6 gate. `docs/tex/pgt_stability_priors.tex` §"TIDAL conventions", `examples/torsion_gertsenshtein/theory_nonminimal.toml:91-97`: tensor `q²` (2⁻) `β₁+β₂ < −½/κ²`; axial `S²` (1⁺) `β₁−β₂ > +¼/κ²`; trace `T²` (1⁻) `β₃+(2β₁+β₂)/3 < 4/(3κ²)`. Working point `(0, −0.6, 0.5)`, recorded margins `−0.1, −0.058, −1.033` | reproduce each margin's **sign** at the working point; locate each boundary by 1-D scan to `10⁻⁶` relative. Sectors `2⁻/1⁺/1⁻` exercise the `J^P` labeling |
| (b) | **Against the published implementations, one per half**: kernel-level vs JAX on the committed WXF fixtures (`VectorTheory`, `TensorTheory`, `S123Theory`, `A23Theory`, `TDiffTheory`, `FierzPauliTheory`); pipeline-level vs the committed Julia `*_unitarity_NS` chains (`FierzPauli`, `FierzPauliProca`, `TDiff`, `Vector`) — no Wolfram or Julia run needed | masses/residues ≤ `10⁻¹⁰` relative (fp64, shared non-degenerate sectors); healthy-region agreement statistical vs chains |
| (c) | **TorC** — primary `ParticleSpectrographTorCECT` (our formulation; 15 gauge generators), `…TorCPGT` (21) secondary; caption at `paper_Qtorsion.tex:656` | pseudoscalar `0⁻`, `m = √λ`; `λ ≥ 0` and `μ < 0` (L1); massless `(ν+2μ)(ν−μ) > 0` (L3, L4 cross-check); 2 massless polarizations (L2); gauge counts from Stage 1 |
| (d) | **GR + Maxwell healthy** — declared massless content `2+2`, no spurious sick verdict | exact; catches rework item 2 and exercises the gauge path where the answer is known by hand |
| (e) | **Cross-branch reconciliation** (restricted limit): solver branch on Minkowski, background block off; vary the spectrum branch's quadratic Lagrangian and match the solver EOM up to the §4.3 dictionary, comparing gauge-invariant quantities (effective signs per `.claude/rules/wolfram.md`) | exact symbolic match of effective signs |
| (f) | **Performance** | median verdict ≤ 1 ms (§9) |
| (g) | **Schur route vs its oracles** (§5.1): `ScalarParityViolatingPGT` intermediate blocks, `K₀`, `M₀`, `c₅ > 0 ∧ 4c₃c₅ − c₄² > 0`, eq. `no_tachyon_cond_3`, `c₁ < 0`, 10 gauge generators; `app:PGT_comparison` (Karananas masses, Blagojević no-tachyon, empty massless spin-one); **and** Schur ↔ residue agreement everywhere both are valid — proven equivalent, so disagreement localizes a bug | matrices to `10⁻¹²` after convention alignment; verdicts exact |

Plus, from §8.1: the deformed → symmetric Maxwell limit reproduces the paper's own statement
(Proca branch: two of three polarizations survive; transverse branch at `Θ₃ = Θ₁ = 0`
propagates nothing).

## 13. Proposed follow-up handoff — Stage 1 is its own workstream (H8)

This document specifies the Stage-1 → Stage-2 **contract**: labeled `J^P` blocks as
coupling-linear coefficient tensors; gauge ranks and the massive/massless partition; the
source-constraint matrix; declared massless content; closed-form conditions when PSALTer
resolves them; scope-guard metadata. The Wolfram engineering behind that contract deserves its
own session:

- **Install PSALTer** (Mathematica 14+/xAct/Inkscape; the SupplementalMaterials repo carries
  `third_party/WolframEngine_14.3.0_LIN.sh`). Only xAct is installed today. The install also
  settles the open items needing a live session: the exact `WaveOperator` label metadata, the
  missing-coupling failure behavior, `Method→"Hard"` cost on a PGT+EM theory (currently
  unmeasured — record, do not guess).
- The **spectrum derivation branch** (§4.5) and the auto-generated model/export script (§4.4,
  §4.5), with the enumerated-symmetry validation and the coupling-per-operator contract.
- **Our own Wolfram-side exporter** reaching into PSALTer's result objects to emit the
  contract above. The published `JuliaExport.m` shows the pattern but emits `.jl` source; the
  WXF slots come from PSALTer itself; the coefficient-tensor and explicit-label export exists
  nowhere and must be written.

**We do not modify PSALTer** — integration means generating its input and harvesting its
output from outside. The division: H8 owns Stage-1 Wolfram engineering (coordinated with M3's
`derive/` port, whose "port and substantial extension" scope in H4 §5.4 the spectrum branch
extends); WS6 proper owns the Stage-2 numerics at `tidalcosmo/spectrum/`, built at M∥ per H4's
migration order (Schur, residue cross-check, verdict/likelihood surface); both build against
this document.

## 14. Open items and risks

1. **`Method→"Hard"` wall time on PGT+EM** — unmeasured; TorC's spectrographs prove
   feasibility. (→ H8.)
2. **Exact `WaveOperator` labeling metadata** — the `J`-block structure is established
   (§6.1), but the per-state ordering convention should be confirmed against a live PSALTer
   and the field-kinematics tables before the exporter hardcodes it. (→ H8.)
3. **Chequer-Hermitian sign conventions** — pinned only up to the `i`-vs-skew-symmetric
   freedom noted in `app:PGT_comparison`; must be fixed against the `K₀` oracle before the
   definiteness test ships. (→ WS6, validation (g).)
4. **`O_LL` conditioning away from strict masslessness** — the Schur complement needs an
   invertible `O_LL`; near-degenerate cases fall into the accidental-degeneracy flag path, but
   the threshold policy needs setting during implementation.
5. **PSALTer's closed-form coverage (L3) is partial for our class** — established (§5.1);
   L4's implementation priority follows from how often L3 comes back unresolved on the actual
   roster.
6. **License resolution before any public release** — #495, options ranked in §3.
7. **The `ΔN_eff` routing** for extra massless poles assumes the program's validity-flag
   machinery (WS2+) exists to receive it; until then the flag is recorded but nothing consumes
   it.
