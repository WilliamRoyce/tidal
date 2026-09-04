> **ARCHIVE — verbatim record of the 2026-08-29 planning session plan file.**
>
> Preserved here because it lived only in `~/.claude/plans/` (not version-controlled;
> its backup directory is gitignored) and carries the full reasoning trail, including
> several rounds of user correction on the physics that shaped the program.
>
> **Not authoritative.** The living documents supersede it wherever they differ:
> `docs/COSMOLOGY_PROGRAM.md` (operations), `docs/cosmology/primer.md` (physics
> orientation), `docs/cosmology/birefringence_notes.md` (O4), and
> `docs/cosmology/spectator_route.md` (spectator scope). Known points where this
> archive is now stale: O1's validation (narrowed by H1 — see #488), the CAMB
> inherit-vs-reapply question (settled: re-apply, #498), the `ν⁰`-vs-`ν²`
> framing (it is per-operator, now #503), **"WS6 independent"** (it is not — it gates
> O2/O3 on the sampling path; H2's dependency graph governs), the single-solver framing
> (O2 and O3 are different numerical problems — H2 §0.1), and the rung order (settled
> 2026-09-04: `O0 → O1 → O2 → O4a → O3 → O4b/V`).

# TIDAL → cosmology: perturbations over cosmic history, CMB observables, real inference

## Context

Supervisor meeting (2026-08-28) set a new direction. TIDAL has so far surveyed *parameter
space* — sweeping PGT couplings on flat Minkowski, scoring against a synthetic objective.
Nothing in it has ever touched real data. The new goal: integrate a candidate Lagrangian's
perturbations over the history of an expanding universe, produce CMB observables, and do
genuine inference against measured data — packaged as a Cobaya extension so others can
test their own Lagrangians against real likelihoods.

**The brief decoded.** The direct predecessor is already in the repo:
`literature/2507.09228/` — **Legner, Handley & Barker, "Alleviating the Hubble tension with
Torsion Condensation (TorC)"**. Its stack is *PSALTer particle spectrum → modified CAMB →
Cobaya → PolyChord*, against Planck 2018 TT/TE/EE+lowE+lensing + SH0ES. It states exactly
the limitation this work lifts:

> "this analysis modifies the background expansion in CAMB … **while the perturbation
> equations remain those of standard ΛCDM**" … "cosmological perturbation theory for TorC
> will be developed in future work."

Public assets from that work, all verified to exist this session:

- `ModifiedCAMB` → **<https://github.com/slegner/CAMB>** (fork of `cmbant/CAMB`, patched to
  read tabulated `ρ_Λ(a)`, `P_Λ(a)`; the `w(a)` interface develops poles when `ρ_Λ`
  changes sign)
- `ModifiedCobaya` → **<https://github.com/slegner/cobaya>** (forked via
  `AdamOrmondroyd/cobaya`, also Handley-group — H1 must diff against both that fork and
  upstream `CobayaSampler/cobaya`)
- Chains/supplementary → **doi:10.5281/zenodo.15866507**

TorC also uses **PSALTer explicitly** (§*Particle spectroscopy*, `paper_Qtorsion.tex:632`:
"we use the *PSALTer* software"), deriving the ghost/tachyon parameter constraints in two
independent formulations as a cross-check. Checked against the bibliography: TorC cites the
**symbolic Wolfram PSALTer** (`Barker:2024juc` = arXiv:2406.09500 and the parity-violating
upgrade `Barker:2025qmw`), *not* the numerical polology code — symbolic was sufficient
because they analyzed **one** theory once. We want the same check for *arbitrary*
Lagrangians and sampled coupling values, which is exactly the regime where the polology
paper says symbolic computation fails ("scales poorly … due to expression swell; the only
avenue is numerical"). That is why WS6 builds the numerical version.

### Background material: what a CMB pipeline computes, and where we plug in

*(Written out because this is a new field for us. H1 confirms it against the TorC text.)*

#### The pipeline, in stages

A Boltzmann code (CAMB, CLASS) turns cosmological parameters into a predicted CMB. Stages:

1. **Background.** Solve for the scale factor `a(η)` — how the universe expands. Driven by
   the total energy density `ρ(a)` and pressure `P(a)` of everything in it.
2. **Thermal history.** When do electrons and protons combine into neutral hydrogen so
   photons stop scattering? Gives `x_e(η)` and the visibility function `g(η)` — the
   probability a photon observed today last scattered at time `η`.
3. **Perturbations.** For each comoving wavenumber `k`, evolve the *lumpiness* — density
   contrasts `δρ`, velocities, anisotropic stress — of photons, baryons, CDM and neutrinos,
   coupled to each other through gravity, from the early universe to today.
4. **Projection.** Convert that evolution into what we actually see on the sky, via the
   line-of-sight integral `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`.
5. **Spectra.** `C_ℓ^{XY} = 4π ∫(dk/k) Δ²_R(k) Δ_ℓ^X Δ_ℓ^Y` — the angular power spectra.

**`C_ℓ` is the observable we compare to data.** Planck/ACT/SPT measure `C_ℓ^{TT}`,
`C_ℓ^{TE}`, `C_ℓ^{EE}` (and lensing); the likelihood scores our predicted `C_ℓ` against
those measurements. That comparison *is* the inference.

#### The two levels, and what TorC actually did

**Level 1 — background (stage 1).** TorC (the model *is* "Torsion Condensation" — a
condensed background torsion field) contributes an *effective* dark-energy density and
pressure `ρ_Λ^eff(a)`, `P_Λ^eff(a)`. "Feeding those to CAMB" means: CAMB needs the
universe's total energy content to solve for `a(η)`. Normally you hand it `Ω_b, Ω_c, Ω_Λ`
and it assumes `ρ_Λ` is constant. TorC's `ρ_Λ^eff` *varies with `a`*. CAMB's standard hook
for that is the equation of state `w(a) = P/ρ` — **but that develops poles when `ρ_Λ`
changes sign**, which TorC's does. So they patched CAMB to read tabulated `ρ_Λ(a)` and
`P_Λ(a)` separately. **That patch is `slegner/CAMB` — the "modified CAMB" is the feeding
mechanism *for background modification only*.** Our integration policy (user decision):
**the package targets the latest upstream CAMB by default** — the spectator spine needs
only standard CAMB products. The `slegner/CAMB` fork is an *optional hook*, used only when
the tabulated-background feature is active (TorC-class analyses, O1 reproduction). H1
decides whether to inherit its patch or re-apply it cleanly against current upstream.

Changing `a(η)` shifts the sound horizon `r_*` and the angular diameter distance `D_A^*`,
hence the acoustic-peak positions `θ_s = r_*/D_A^*`, and the early ISW — all real,
measurable effects on `C_ℓ`. That is what TorC computed and compared to Planck.

**How a general Lagrangian becomes `ρ^eff(a)`, `P^eff(a)` — the effective-fluid method**
*(context for the optional background feature only; the spectator spine never does this).*
"Field equations" here means the **gravitational** field equations — the theory's
generalization of Einstein's equations, obtained by varying the full action with respect
to the metric (and, in PGT, the connection); *not* the matter/torsion-field EOMs. (Yes —
these are the same equations whose *linearization* gives the `h`-component equations TIDAL
already derives; here they are used un-linearized, evaluated on the FRW ansatz instead.) Evaluate
them on the homogeneous-isotropic (FRW) ansatz; the **time-time (`00`) component of that
tensor equation** is the energy-constraint (Friedmann) equation and always takes the form
`3 M_p² H² = ρ_standard + (everything else)`. One *defines* `ρ^eff` as everything beyond
the GR terms, and `P^eff` similarly from the spatial-trace component. This is the standard
modified-gravity move — the new sector is disguised as a fluid CAMB already knows how to
carry.

**What TorC's "bi-scalar-tensor equivalent" is, and why it is NOT our method.** TorC
hand-derived a mapping specific to their theory: on homogeneous-isotropic backgrounds, the
24 torsion components reduce to two scalars (`ϖ`, `φ`), and the PGT Lagrangian maps to an
equivalent scalar-tensor theory whose Friedmann equations are then read off. That is
precisely the **per-theory hand-crafting this project exists to eliminate** — it required
theory-specific insight, holds only for that model class, and does not generalize. Its
only role for us: *if* the optional minisuperspace feature is ever built (generic FRW
reduction of any TOML Lagrangian, no equivalent-theory tricks), TorC's published result is
its validation case. It is on no critical path. The spectator spine needs no background
derivation at all.

**Level 2 — perturbations (stage 3).** Here my earlier phrasing was wrong, and the paper is
explicit. ΛCDM contains no torsion, so there is no such thing as "ΛCDM's torsion
perturbations." What TorC did is **omit the torsion perturbations entirely and use GR's
perturbation equations as a stand-in**:

> "This work focuses on the background evolution … and does not yet incorporate
> cosmological perturbations. For the purposes of this initial study, the perturbation
> equations remain those of standard ΛCDM." … §*Fiducial perturbation theory*: "we study
> the background dynamics in isolation, **adopting GR perturbation theory as a fiducial
> proxy**."

So the torsion was allowed to change *how fast the universe expands*, but not to have its
own ripples that gravitate, cluster, free-stream, or carry anisotropic stress.

**Including those torsion perturbations is precisely this project's work — it is O2, the
first genuinely new result.** TorC names it as its own future work; we are building the
engine that does it, for any Lagrangian.

**Promoting to Level 2** means deriving the torsion perturbations' own equations of motion
from the Lagrangian and evolving `δT(k,η)` alongside the standard species. **And yes — the
linear perturbation equations TIDAL already derives symbolically are exactly the object
that gets routed into the Boltzmann code**: re-derived about an FRW background instead of
Minkowski, they become the new evolution equations added to stage 3, coupled to the
standard species through the Einstein equations. TIDAL is the front end that generates
what the Boltzmann code evolves. The new-sector perturbations then:

- contribute `δρ, δP` to the Einstein constraints → change how gravitational potentials evolve;
- contribute **anisotropic stress** `σ`, which alters ISW **and lensing**. (Notation, since
  these are standard objects, not new fields of ours: `φ` and `ψ` are the two scalar
  *gravitational potentials* of the perturbed FRW metric,
  `ds² = a²[−(1+2ψ)dη² + (1−2φ)δᵢⱼdxⁱdxʲ]` — `ψ` is the Newtonian-potential-like
  time-time perturbation, `φ` the spatial-curvature one. In GR with no anisotropic stress
  they are equal; a sector carrying `σ` splits them via
  `k²(φ − ψ) = 12πG a²(ρ̄+P̄)σ`, and that split is what shifts ISW/lensing.);
- propagate with their own sound speed → change the driving of the acoustic oscillations;
- if parity-odd, **mix E and B polarization** — "polarizations converting into each other."

Running that through stages 4–5 is **"integrating the perturbations over the history of the
universe."**

**Why it matters — and TorC invites it directly.** The paper argues TorC is a research
*programme* with several realisations (bi-scalar-tensor, "tensor bypass", full PGT) that
**share background dynamics but differ in their perturbations**: *"Different realisations
may yield distinct perturbation theories."* So the perturbations are precisely what
discriminates between them — background-only cannot. Additionally, background-only is an
approximation of unquantified validity: if torsion perturbations are not negligible, using
GR's biases the inferred parameters, and the Hubble-tension claim rests on those parameters.

#### What we reuse versus what we build

**Yes — your reading is right, and it is the key architectural decision.** The ΛCDM sector
is enormously well solved: photons, baryons, CDM, neutrinos, recombination, the full
Boltzmann hierarchy, validated across 12 independent codes to sub-percent agreement
(arXiv:1709.09135). **We do not reimplement any of it.** We take it from CAMB and add only
our new sector, coupled in through the Einstein equations. That is the entire reason to
chain onto CAMB rather than write a Boltzmann code.

**Do we have to evolve the photons?** Depends on the theory, and this is what orders the
observable ladder:

- If the new physics is **gravity-sector only** (O1, O2): photon transport is unmodified
  ΛCDM — take it from CAMB untouched.
- If the Lagrangian **couples torsion to electromagnetism** (O3, O4 — the thesis's own
  physics): photon propagation *is* modified and we must intervene — either as a modified
  line-of-sight source (Tier 3a below) or by evolving a modified photon sector (Tier 3b).
  **This is a core goal of the programme, not an optional extension**: the thesis is built
  on the Gertsenshtein coupling, and seeing those couplings' cosmological effects is a key
  use case the supervisors want. The engine must handle modified photon sectors as a
  first-class case.

#### Where the new physics enters — and why the tier table matters

Stages 1–5 are a fixed pipeline. The practical question is **at which stage our new physics
gets inserted**, because that determines how much code we write and whether we can stay in
Python. The tier table in WS0 answers exactly that: Tier 1 touches only stage 1 (a tabulated
background — what TorC did, and our O1); Tier 3a adds a term to the stage-4 source function
and is pure Python via `camb.symbolic`; Tier 3b adds a *new evolving field* to stage 3 and
requires going into the Fortran/C, or using a code built to be extended. **TIDAL's torsion
is Tier 3b** — that is the hard part, and it is exactly what TorC deferred.

#### Division of labor: what CAMB computes, what our solver computes, what nobody computes

The clean split, per likelihood call:

| Computed by | What | When it recomputes |
|---|---|---|
| **CAMB (unchanged)** | ΛCDM background `a(η)`, `H(η)`; thermal history `x_e(η)`, `g(η)`; all standard-sector perturbations (photons, baryons, CDM, neutrinos); standard `C_ℓ`, lensing, matter power | only when ΛCDM parameters move (slow block; cached otherwise) |
| **Our solver (WS3)** | **the coupled block**: the new sector's modes `δT(k,η)` *together with the standard modes they directly couple to* (e.g. the tensor `h` in a graviton–torsion system; the photon polarization modes in a torsion–EM system) — a small linear ODE system per `k`, coefficients involving the CAMB-supplied `a(η)`, `H(η)`. We evolve the coupled subsystem; we do **not** re-solve the full Boltzmann hierarchy — everything not coupled to the new sector stays CAMB's | every call that moves a new-sector coupling (fast block) |
| **Nobody** | nothing ΛCDM is re-derived, re-implemented, or re-solved by us — ever | — |

**Why an in-house solver at all, if we're plumbing into existing packages?** Because no
existing package can evolve *our* equations: CAMB/CLASS hard-code the ΛCDM species; they
have no notion of a user-derived field with TIDAL-generated equations of motion. The one
genuinely new numerical task — integrating a small, k-parametrized, time-dependent-
coefficient linear system — is exactly TIDAL's existing competence (the modal engine),
generalized to `M = M(η)`. That generalization is WS3.

**How it slots back into the packages — and what a "Theory class" is.** Cobaya organizes
an inference run into pluggable components: *Theory* components compute predictions,
*Likelihood* components score them against data, and a *sampler* drives the loop. A
"Theory class" is simply a Python class implementing Cobaya's small interface
(`get_requirements()` — what it needs from other components; `calculate()` — what it
computes; `get_X()` — what it offers). **Yes, we write one — that IS the Cobaya extension.**
Ours declares that it requires CAMB's products, runs our per-k solver for the coupled
block, and returns the modified observables — either (a) a **replaced transfer function**
(the tensor mode `h(k,η)` is now the solution of our coupled graviton–torsion system
instead of CAMB's plain `h'' + 2ℋh' + k²h = 0`), or (b) a **conversion/rotation applied to
CAMB's output** (spectral-distortion amplitudes, polarization rotation of `C_ℓ`).
Downstream likelihoods see standard products and never know the difference.

**A worked chain, concretely (tensor channel), to make "evolution → observable" explicit:**

1. Inflation supplies the primordial tensor amplitude per mode `k` (two numbers: `A_t`,
   `n_t` — sampled parameters).
2. For each `k` on a grid: our solver integrates the coupled `(h, torsion)` system from
   deep radiation era to today in `η`, with `a(η)`, `H(η)` read from CAMB's table. The
   output is the **transfer function** — how much of the initial `h` survives (or is
   amplified/converted) by recombination and after.
3. CAMB's line-of-sight machinery projects that source onto the sky:
   `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`, then `C_ℓ^BB = 4π∫(dk/k) Δ²_t(k) |Δ_ℓ(k)|²`.
4. A likelihood (BICEP/Keck, LiteBIRD) scores that `C_ℓ^BB` against measured B-mode data.

Same shape for the photon channel: the torsion background modifies polarization
propagation along the line of sight → rotation/conversion accumulated between
recombination and today → applied to CAMB's `C_ℓ^{EE}` → scored against `EB`/`V` data.
**And yes — the CMB is literally those photon quanta**, free-streaming since
recombination; anything that alters them en route is imprinted on what we measure. That
is why propagation effects are so constraining.

**Resolving the `a(η)` analytic-vs-tabulated tension (blocker 5, precisely).** Two
different stages want `a(η)` in two different forms, and the current pipeline conflates
them: the **Wolfram derivation** stage today requires the metric as an *analytic string*
(that is how `de_sitter.toml` works), but the real ΛCDM `a(η)` exists only as CAMB's
*numeric table*. The resolution is to separate the stages: derive the equations
**symbolically with `a(η)`, `H(η)` left as unspecified background functions** (the
equations' coefficients then *contain* symbolic `a`, `H` — and yes, this is exactly how
the expansion details enter: Hubble-friction terms `∝ H`, mass terms `∝ a²`, etc.), and at
**solve time** evaluate those coefficients numerically from the CAMB table. Symbolic
derivation once per theory; numeric evaluation per call. The analytic-metric mode stays
for validation cases (de Sitter) where exact solutions exist.

#### Where TIDAL fits, and why Cobaya

TIDAL already derives linear perturbation equations from a Lagrangian symbolically — that
is its entire purpose. **Yes: that is the tool.** TorC would otherwise hand-derive them per
theory. Automating it *for any Lagrangian* is what makes this a general extension rather
than a one-off.

And Cobaya is the payoff. Once our component supplies `C_ℓ` (and lensing potential, matter
power spectrum, …) through Cobaya's standard interface, **every likelihood in that
ecosystem becomes available for free** — Planck, ACT, SPT, DESI/BAO, supernovae, weak
lensing. We build the theory tool once; testing it against many different datasets is then
configuration, not code. Equally important for impact: it makes the tool **drop into
people's existing setups** — a Cobaya component is something the community already knows
how to install, configure and cite, which is what turns a research code into an adopted
(and cited) one. That is the argument for chaining into Cobaya rather than building a
bespoke pipeline.

**Scope.** Programme plan spanning many sessions. Nothing implemented in this session.
This session acts as **orchestrator**: it holds the plan and issues handoff prompts to
other sessions (as with the campaigns). See *Session handoffs*.

## Decisions taken (user)

| # | Decision |
|---|---|
| D1 | Target **CMB observables via genuine perturbation evolution**. The point is the **general engine**: once built, many effects become investigable (growth/ISW/lensing, Gertsenshtein mixing, birefringence, …). Birefringence (E→B) is *one* supervisor-stated end goal among several, not the destination the engine is shaped around — see the observable ladder. |
| D2 | **No paper target fixed.** Both a physics paper and a software/method paper expected; order workstreams by dependency. |
| D3 | **Reshape the repo around the Cobaya-extension goal.** Where restructuring is the right answer it **should be done**, not merely tolerated — no attachment to the current layout. It grew without a target; the target now exists, and the design follows from it. Git history preserves the old shape, so nothing is lost by moving. WS1/H4 is explicitly licensed to propose deleting and re-laying-out packages, not just relocating files. |
| D4 | **No HPC without explicit permission.** Local only. |
| D5 | **Optimization is a standing first-class concern**, not a later phase — adoption depends on it. |
| D6 | **PSALTer-numerical: build our own version, taking heavy, deliberate inspiration from Barker's numerical code** (arXiv:2606.30785 + `psalter.tar.gz` on disk + the public SupplementalMaterials release). The author has **explicitly permitted** using their code to write ours — reading and copying is fine (user). Requirements: provenance recorded in docstrings (what came from where); formal attribution settled with the supervisor at publication time. License concerns are dropped — not a real constraint here (user). Not a clean-room reimplementation; the released versions serve as cross-checks, not as adversarial benchmarks. |
| D7 | **The solver targets the general case.** Cancellations are *fast paths the engine detects automatically* (as `can_use_modal` already does), never grounds to specialize to one observable. |
| D8 | **Don't overspecialize early.** The Chern–Simons/birefringence channel is one instance; the engine and the plan stay general. |
| D9 | **Migrate inference into the Cobaya ecosystem.** We currently drive `pypolychord` ourselves; Cobaya ships a PolyChord sampler, so this is a simplification — but much of `tidal/inference/` is then superseded and the redesign must be deliberate, not bolted on. |

## Parked

The **#477 arc** on `feat/ws2-localized-path-audit` is **already halted** (user): issues
arose that made it not possible, and those issues are likely unrelated to the new
direction. Park it as-is — record its stopping state on the issue/branch, no further work.
Also parked: Phase E (T6/T8), Phase A-γ, `NEXT_PHASES` G/H/I, Wolfram CI (#69). Nothing
gates on any of it; the cosmology workstreams start immediately.

## Observable ladder

Ordered by **pipeline-validation value**, not physics novelty. Each rung has a known or
checkable answer before the next adds unknowns. This replaces committing to a single
observable up front (D8).

| # | Observable | New capability exercised | How it's validated |
|---|---|---|---|
| **O0** | ΛCDM `C_ℓ` with new physics off | plumbing only: derive → background → CAMB → Cobaya → sampler | sub-percent vs CAMB/`nanoCMB`, `2 ≤ ℓ ≤ 2500` |
| **O1** | **Recreate the TorC analysis** through our pipeline, using TorC's *published* `ρ_Λ^eff(a)`, `P_Λ^eff(a)` as input tables. TorC is a **validation case, not a foundation** — the tabulated-background input is a general optional feature; the package depends on no particular model | end-to-end plumbing: Cobaya wiring, tabulated-background input, sampler, likelihoods | **against published chains** (Zenodo `10.5281/zenodo.15866507`). **Best first target.** |
| **O2** | **Spectator perturbations on ΛCDM** — the new sector's `δT(k,η)` evolved on the CAMB background (test-field mode); first observable: **modified tensor/GW propagation** (friction `ν`, mass `μ`, speed `c_T` derived from the Lagrangian) → tensor transfer function → B-modes | the spine: TIDAL-derived FRW perturbation equations + the WS3 solver, in the strict spectator regime | the thing TorC explicitly deferred; first genuinely new result. Validity flags reported per run |
| **O3** | **Gertsenshtein mixing on FRW** — graviton↔photon with Hubble damping | the time-dependent solver, on the thesis's own physics | **already scoped as GH #209**, citing Cembranos arXiv:2302.08186 (`h'' + 2H(1+ν)h' + (c_T²k² + μ²)h = 0`), local. Also `literature/2312.17636/` (inverse Gertsenshtein → CMB spectral distortion). Reduces to the known flat-space result as `H → 0`. |
| **O4** | **Cosmic birefringence (E→B)** | parity-odd sector + polarization observables | supervisor's stated end goal; currently ~4.8σ (`β = 0.277° ± 0.057°`, joint ACT+Planck, arXiv:2608.06480) |

**O0 first, then O1 (user).** O0 is not optional throat-clearing — it is the validation
gate. Precisely: our Cobaya component in **pass-through mode** (new physics off) must hand
CAMB's products through unmodified, reproduce the standard ΛCDM `C_ℓ` and a standard
ΛCDM posterior on a reference likelihood. Where any piece of the chain is ours rather than
CAMB's, that piece must agree with CAMB/`nanoCMB` sub-percent for `2 ≤ ℓ ≤ 2500`. Only
then O1.

**The central architecture — the spectator (test-field) route (user + supervisor).**
This is the regime the supervisor described in the meeting, and it is the programme's
spine, not one mode among several:

> **The background is always the established ΛCDM one, supplied by CAMB. The new sector's
> perturbations are evolved on top of it as test fields — "perturbations about the
> established background" — and their observable imprints (on gravitons, photons,
> polarization, lensing) are what we compare to data.**

Why this is the right spine:

- It is fully within TIDAL's linearized-only design — no nonlinear background derivation
  on the critical path at all.
- It is the regime the thesis worked in anyway (small couplings, linear mixing).
- It makes the per-likelihood-call cost small: CAMB's ΛCDM work is standard and cacheable;
  the only new computation per call is the new sector's mode evolution.
- Its validity condition (`ρ_new ≪ ρ_total`, plus small conversion probabilities) is
  checkable and **must be computed and reported per run** — the honest-flags philosophy
  the project already follows (gauge certificates).
- The same separation TorC uses for spectra (vacuum screening vs background dynamics)
  extends here: screen with PSALTer in vacuo, evolve as a spectator on ΛCDM.

**H7 (spectator-route scope investigation) has reported. Key findings, now load-bearing:**

- **The niche is empty — verified by enumerating the near-misses.** Symbolic
  Lagrangian→perturbation-equations exists (xPand/xAct — no numerics, no observables);
  automatic Lagrangian→solver→observables exists for inflation only
  (CppTransport/PyTransport, arXiv:1609.00380/1); covariant-Lagrangian→Boltzmann→likelihood
  exists for the Horndeski *gravity* sector with *full backreaction* (hi_class,
  arXiv:1909.01828); symbolic-model→differentiable-Boltzmann exists but takes equations,
  not actions (SymBoltz.jl, arXiv:2509.24740 — watch it; possible future backend). **No
  framework combines them for a general extra sector in the spectator limit, and nobody
  has treated torsion this way at all** — 2506.17017 and TorC both signpost cosmological
  perturbation theory as future work. The defensible novelty claim is *the combination
  plus the spectator restriction* (which is exactly what makes a general Lagrangian
  tractable where hi_class-style full self-consistency would not be).
- **The consistency of coupling-without-backreaction is a clean double expansion**, stated
  verbatim in Cembranos arXiv:2302.08186: order 0 in perturbations = background Einstein
  equations (discarded, replaced by CAMB's solution — the *only* place the spectator
  assumption enters); order 1 = tadpoles (vanish on-shell); order 2 = the full quadratic
  action *including all mixing terms* (kept whole). Coupling and backreaction sit at
  different orders; no inconsistency.
- **Quantitative validity criteria to enforce per run — a supervisor-flagged concern,
  hence a first-class requirement, not a nicety.** The supervisor raised that checking
  this approximation's validity is genuinely harder in this regime: we take the spectator
  approach, but must continuously verify the new perturbations are *not* feeding back into
  the overall growth. The literature asserts validity in one sentence and never enforces
  it — enforcing it numerically is both the safeguard and a cheap methodological
  contribution (the project's existing honest-flags style). Monitors: (a) `ρ_new/ρ_γ`
  against the `ΔN_eff ≲ 0.1` bound (Domcke & Garcia-Cely arXiv:2006.01161, local);
  (b) conversion probability `P_max ≪ 1`; (c) amplitudes `|h|, |f| ≪ 1`; and (d) a
  **growth-impact monitor**: per mode, the ratio of the new sector's source terms in the
  Einstein constraints to the standard sectors' — directly answering "would these
  perturbations have affected the growth we froze?"
- **⚠ The one silent-failure risk:** the CAMB background solves *Einstein's* equations,
  but our quadratic action comes from a PGT Lagrangian — consistency requires the PGT
  background field equations be satisfied by (FRW, `T̄ = 0` or tracking torsion) to the
  needed accuracy. `literature/2003.02690/` (the group's PGT background-cosmology paper,
  local) is where the tracking/frozen solutions live. The needed check is a
  **background-EOM residual** evaluated on the CAMB background. The *concept* is proven
  useful (the #477 arc measured exactly this quantity in another context), but per the
  design rule the implementation is **designed fresh for the new architecture** — the
  legacy machinery is not carried over.
- **Strategic constraint:** CMB bounds from conversion channels are typically weaker than
  the `N_eff` bound *unless* the signal is spectrally narrow or structurally distinctive —
  which is precisely what V-modes (arXiv:2502.12517, recent and thin literature) and
  birefringence are. This justifies the polarization-channel emphasis on physics grounds,
  not just because the supervisor named it.
- Useful negative result for O3 (Cembranos, local): the graviton mass `μ ≤ 10⁻³³ eV`
  changes nothing, but the friction `ν` matters strongly (amplification for `ν < 0`,
  vanishing at `ν = −1`) — the interesting torsion target is the friction/running term.

**Background modifications are optional extensions, not the spine.** Two explicitly
optional capabilities, kept because they widen the tool and enable validation:

- **Tabulated-background input** (`ρ^eff(a)`, `P^eff(a)` as tables): lets the pipeline
  *recreate* TorC-class background analyses. **Validation/compatibility feature only —
  the package must not depend on TorC or on any particular model.** TorC is a test case,
  not a foundation.
- **Minisuperspace derivation** (FRW-reduce the full TOML Lagrangian, vary in 0+1D):
  automates producing such tables from a Lagrangian. Optional extension, later.

*Ladder ↔ workstream interplay:* O1's critical path is H1's inheritance decision + a thin
slice of WS5 (Cobaya wiring) only. It waits on **neither** WS3 (the time-dependent solver
— first exercised by O2/O3) nor the optional background capabilities.

**O1 is confirmed as the first science target (user).** It validates the *entire* stack
against a published result from the same group, needs no new perturbation theory, and
automating the Lagrangian→`ρ_Λ(a),P_Λ(a)` step is publishable on its own. It also fails
loudly and unambiguously if the plumbing is wrong. H2 confirms the ordering of O2–O4 behind
it.

**The spectator boundary — a correction from H7 to the earlier O2 framing.** First, what
"spectator limit" means, plainly: **(a)** the new sector does not affect the background
expansion (its energy density is negligible in the Friedmann equation); **(b)** it does
not gravitationally disturb the standard perturbations (its energy lumps are negligible in
the Einstein constraint equations); **(c)** everything stays linear (all perturbations
small). Yes — "doesn't affect the background, and perturbations assumed to remain small"
is exactly it, with (b) as the extra subtlety.

Two clarifications on how this coheres:

- **"Standard perturbations don't feed back on the background" is not our assumption — it
  is standard cosmological perturbation theory, used by everyone.** Linear perturbations
  never back-react on the background at first order, and `δρ/ρ ~ 10⁻⁵` at the CMB makes
  the truncation superb. Our *additional* assumption is only about the **new** sector:
  that its *background* energy density is also negligible. Two separate assumptions; only
  the second is ours to monitor.
- **"But we modify the propagation of standard perturbations (photons) — doesn't that
  violate (b)?" No — (b) is about gravity; our modification is a direct coupling.**
  Condition (b) forbids the new sector influencing standard perturbations *through its
  stress-energy in the Einstein equations* (gravitational sourcing). What we do instead
  is influence them through the **explicit coupling terms in the Lagrangian** (the
  torsion–photon and torsion–graviton mixing terms), kept fully at linear order. That
  channel is the entire point of the calculation and is untouched by the spectator
  truncation — coupling-mediated, not gravity-mediated.

Within that limit, three classes of effect, of which we reach two:

1. **Propagation (reachable).** Standard quanta — photons, gravitational waves — travel
   *through* the new sector's background and couplings, which alter *how they travel*:
   speed, phase, damping, polarization. No energy changes hands. Examples: polarization
   rotation (birefringence), modified GW friction/dispersion. Observable because we
   measure those quanta exquisitely well.
2. **Conversion (reachable).** Quanta *oscillate between* sectors — a graviton becomes a
   photon (Gertsenshtein) and vice versa — moving a *small* fraction of energy across.
   Observable as small additions to or distortions of standard signals: CMB spectral
   distortion, V-mode polarization, excess radio background. Stays consistent with the
   spectator limit precisely because the converted fraction is tiny — and that tiny
   fraction *is* the measurement.
3. **Gravitational sourcing (NOT reachable in the strict limit).** The new sector's own
   energy lumps `δρ, σ` pulling gravitationally on photons and matter, *creating* CMB
   anisotropy the way CDM does. That requires its stress-energy inside the Einstein
   constraints — which is exactly assumption (b) being dropped. Doing it is the
   axionCAMB-style full-component route: legitimate, but a deliberate later extension,
   not silently blended into O2. O2 as scoped stays strictly inside the limit.

**O2's scientific motivation comes from the TorC paper itself:** it presents TorC as a
research *programme* with several realisations (bi-scalar-tensor, "tensor bypass", full PGT)
that **share background dynamics but differ in their perturbations** — *"Different
realisations may yield distinct perturbation theories."* The perturbations are therefore
what discriminates between them, and background-only cannot. That is a direct invitation to
exactly this work.

**Reference material already local for the later rungs** (do not re-derive): for O3, GH
#209 + `literature/2302.08186/`; for O4, `research/lagrangian_enumeration/` §`sec:chern_simons`
(line 529) enumerates `L_CS = ε^{μνρσ}V_μA_νF_ρσ` — the Carroll–Field–Jackiw structure —
with `CS: eps Tor A F` (6 terms) marked `ghost_status: safe`, `topological: False`, and
`literature/gr-qc_0307063/` (Itin & Hehl, verified this session) derives a torsion-induced
axion giving frequency-independent optical activity.

## Existing assets (audited this session)

- **FRW is already derivable.** [de_sitter.toml](examples/curved_spacetime/de_sitter.toml)
  uses `diagonal = ["-Exp[2*dSH*t[]]", ...]`; `IsNonConstantMetric`
  ([CommonUtilities.wl:462](tidal/wolfram/CommonUtilities.wl#L462)) gates Christoffels; the
  `√|g|` E-L measure produces the **correct** FRW Hubble friction, ground-truthed against
  `(D−2)H` in [test_component_euler_lagrange.wls:84-112](tests/wolfram/test_component_euler_lagrange.wls#L84-L112)
  (GH #394). `time_dependent` flows Wolfram → JSON → the L3 cache tier.
- **Forward model:** `run_inference_step(...) -> SimulationData`
  ([_sweep.py:619](tidal/cli/_sweep.py#L619)) + `_measure_from_sim_data(...)`
  ([_sweep.py:687](tidal/cli/_sweep.py#L687)), in-memory, ~82 ms/call.
- **Derived parity-odd theory:** `examples/data/torsion_gertsenshtein_parity_odd.json`
  (T6, 22 params).
- **Anchor issues:** **#209** (cosmological Gertsenshtein — becomes O3) and **#43**
  (priority-high, "whether only Minkowski is sufficient" — this pivot answers it).

## Blockers

1. **The modal solver cannot survive this** — it refuses non-trivial `volume_element`
   ([modal.py:371-372](tidal/solver/modal.py#L371-L372)) *and* any `time_dependent` term
   ([modal.py:409-413](tidal/solver/modal.py#L409-L413)); the analytical Jacobian too
   ([analytical_jacobian.py:917](tidal/solver/analytical_jacobian.py#L917)). Replacing it
   is a **research workstream** (WS3), not a port.
2. **Library core lives inside the CLI package — layering inverted.** `tidal/cli/` is
   25,566 of 67k lines and holds `_simulate`, `run_inference_step`,
   `_measure_from_sim_data`. `tidal/inference/_likelihood.py` (5 sites),
   `measurement/_posthoc_audit.py`, `inference/_prior_stability.py`,
   `measurement/_sweep_results.py` all import *private* functions back out of it, each
   annotated `# pyright: ignore[reportPrivateUsage]` or `# noqa: PLC0415`. This is *why*
   the forward model takes an `argparse.Namespace`: the CLI **is** the config object.
3. **Measurements return scalars only.** `C_ℓ` is an array; needs a vector path.
4. **No likelihood compares to data.** `maximize`/`minimize`/`extremize` are optimization
   objectives built for the thesis campaigns; `gaussian` targets a CLI-typed scalar.
   Superseded by Cobaya (D9); removed once the new objectives are in place, per the WS1
   strangler policy — not before.
5. **`a(η)` must be an analytic Wolfram string.** Real ΛCDM `a(η)` is tabulated from a
   Boltzmann code — exactly the problem `ModifiedCAMB` solved by reading tabulated inputs.
6. **Time-dependent energy path is wrong** — [ExportJSON.wl:1638](tidal/wolfram/ExportJSON.wl#L1638)
   filters `t` out of Hamiltonian coordinate deps; [_energy.py:993-1019](tidal/measurement/_energy.py#L993-L1019)
   hardcodes `t=0.0`. *A correctness bug to fix (silently wrong results are never
   acceptable), but low priority — energy is no longer a target observable.*
7. **Scope conflict to reverse:** the 2026-08-17 PSALTer plan proposed making Minkowski
   first-class project-wide and deferring curved backgrounds. Superseded for the solver.
   **PSALTer itself stays Minkowski-only and that is fine** — per the supervisor, the
   particle spectrum is evaluated around flat vacuum and everything else is treated as
   perturbations. To be explicit about "TorC used exactly this split": the TorC paper
   computes its particle spectrum with PSALTer **around Minkowski with vanishing background
   torsion** (§Particle spectroscopy: the spectrum is "very well understood around
   Minkowski spacetime with *vanishing* background torsion") while running its cosmology on
   FRW, and explicitly defers spectra around the condensate as future work. Same
   separation, same justification — vacuum health screens the theory; the background is a
   separate computation.
8. **Doc drift:** `CLAUDE.md`/`README.md`/`docs/tex/inference.tex` document a dynesty
   backend absent from the code — historical residue from before the PolyChord integration
   replaced it (per user). Pure cleanup; fold into the doc-update pass.

## Workstreams

Order: **WS0 → WS1 → WS2 ∥ WS3 → WS4 → WS5.** WS6 independent.

### WS0 — Research & scoping *(no code)*

Delegated to other sessions via the handoff prompts below; this session integrates the
results into `docs/COSMOLOGY_PROGRAM.md` and files the issue tree.

Must settle: the observable-ladder ordering; the integration-target decision; and how much
of `slegner/CAMB` + `slegner/cobaya` we inherit versus rebuild for the optional
tabulated-background hook.

**The integration-target decision (H3), spelled out.** For the parts where standard-sector
evolution must be modified (Tier 3b), the candidates are: (i) patch CAMB's Fortran;
(ii) build on **DISCO-EB** — a from-scratch Einstein–Boltzmann solver in JAX
(arXiv:2311.03291): it does the same job as CAMB (background + perturbations + `C_ℓ`),
agrees with CAMB at the per-mille level, and was explicitly designed to be modular and
extensible with new physics, plus it is differentiable (gradients for free — enables HMC
samplers); (iii) evolve only our coupled block in our own solver and hand results back to
unmodified CAMB (the pure chaining route — the current default assumption). H3
investigates and recommends; this is delegated, not decided here. Note (user): **language
is not a barrier** — with Claude Code, Fortran/C/Julia are all writable; the decision
should be made on architecture and maintainability, not language familiarity, with
environment setup costed in.

**Time variable (user question, answered as a decision):** the cosmology package works in
**conformal time `η`** — the standard variable of cosmological perturbation theory (CAMB's
internal variable; the equations take their cleanest form, `ds² = a(η)²[−dη² + dx²]`). The
existing derivation pipeline uses coordinate time `t`; supporting `η` as the time
coordinate in the symbolic derivation is part of WS2's coordinate work. Conversions
`t ↔ η ↔ z` happen at the CAMB interface, which supplies all three.

**Where new physics enters the pipeline — this table sets WS2/WS5 scope.** The five pipeline
stages (background → thermal history → perturbations → projection → spectra) are fixed
machinery we inherit. What varies between theories is *which stage* the new physics touches,
and that determines how much we write and in what language. It is the single biggest driver
of cost, which is why it is settled in WS0 before any code:

| Tier | Enters via | Implementation | Precedent |
|---|---|---|---|
| 1 | background `H(a)` only | tabulated `ρ(a),P(a)` into patched CAMB | **TorC** (`literature/2507.09228`) → our **O1** |
| 2 | ionization history `x_e(z)` | patch RECFAST | `literature/2204.06302` |
| 3a | new **line-of-sight source** from existing ΛCDM variables | `camb.symbolic` + `set_custom_scalar_sources()` — pure Python, cross-spectra free | CAMB docs |
| 3b | **new propagating DOF**, `φ ≠ ψ`, tensor/parity-odd sector | Fortran/C, or DISCO-EB | arXiv:2311.03291 |

TIDAL's torsion is **3b** — what TorC deferred, and the heart of this programme.

### WS1 — Repo reshape *(D3, D9; unblocks WS5)*

**Strategy: strangler-fig, not big-bang (user direction) — with a hard design rule: the
existing framework is legacy, not a template.** The current structure is acknowledged as
naive and gradually patched-on; nothing new is modeled on it, and holding onto its
patterns is explicitly what we are *not* doing. We do *not* refactor the existing 67k
lines in place. Instead:

1. H4 designs the ideal target layout for a Cobaya extension **from the goal backwards,
   from first principles** — the legacy code is consulted only for physics content and
   test oracles, never as a design reference.
2. **New packages are written clean from the ground up in that layout**
   (`tidal/cosmology/`, `tidal/background/`, the new solver work) — the existing code is
   left untouched and keeps working, so nothing is destabilized and the thesis-era
   pipeline stays reproducible throughout.
3. **New code never reaches into old code — no adapters, no shims, no imports (user
   directive).** Where an existing tool is genuinely useful (an operator kernel, a
   measurement algorithm, the coefficient evaluator's ideas), it is **fully ported**:
   redesigned to fit the new architecture and moved over — **preserving its documentation
   in the move: docstrings, why-comments, and GitHub-issue references travel with the
   code**. The port is reviewed against the original as its test oracle. The old copy
   stays in place, functional, untouched.
4. **Two coexisting, fully separated packages/CLIs during the transition.** The new
   package gets its **own entry point** — its CLI (and Cobaya component) never mixes with
   the legacy `tidal` CLI. The old `tidal` toolchain remains available as-is for
   cross-checking new results against thesis-era behavior. H4 designs the exact
   separation (naming, package boundary, whether new lives under `tidal.*` or a fresh
   top-level namespace — user preference to be captured there).
5. Legacy paths are retired **deliberately**, once the new package demonstrably covers
   their role — thesis-objective code (`maximize`/`minimize`/`extremize` likelihoods,
   sweep panels, …) included. Never two live implementations *presented* for the same job;
   the legacy one is explicitly labeled legacy.

This converts WS1 from a blocking overhaul into building the correct package from scratch
beside a frozen reference implementation — WS5 depends only on items 1–3.

Retained specifics:

- `SimulationConfig` spec source: `_backfill_simulate_args`
  ([_posthoc_audit.py:278-318](tidal/measurement/_posthoc_audit.py#L278-L318)) enumerates
  the ~20 attributes `_simulate` direct-accesses, with rationale comments.
- **Fate of `tidal/inference/`** (D9): Cobaya supplies priors, samplers (including
  PolyChord), and output. The capabilities worth having (cubed-sphere joint prior, KL/BMD
  importance diagnostics) are **re-implemented in the new structure as Cobaya-compatible
  components, informed by the old code** — not lifted from it; the new package's cohesion
  comes from being designed whole for the end goal, not assembled from ad-hoc parts.
  Retire the rest as the new path covers it. Never two live inference stacks for the same
  job.
- Target layering: `symbolic`/`solver`/`measurement` (no CLI deps) → `driver` →
  `inference`, `cosmology`, `spectrum` → `cli`.

The suite must stay green at every step — old tests keep passing untouched until the code
they cover is retired deliberately.

### WS2 — General expanding background *(D7)*

- **Spectator mode is the engine's spine:** new-sector perturbations evolved on the
  standard ΛCDM background from CAMB, with the `ρ_new ≪ ρ_total` validity condition
  computed and reported per run. FRW-background derivation in TIDAL means: re-derive the
  new sector's linear EOM about the CAMB-tabulated `a(η)` instead of Minkowski.
- **Tabulated-background input** (`ρ^eff(a)`, `P^eff(a)` tables → patched CAMB): optional
  general feature enabling TorC-class background analyses and the O1 validation. No model
  dependence.
- **Minisuperspace derivation** (FRW-reduce the TOML Lagrangian, vary in 0+1D): optional
  later extension automating table production; reference = TorC's bi-scalar-tensor
  reduction (H1 extracts it).
- Conformal time + comoving coordinates as first-class. Time is currently hardcoded as
  `coords[0]` across `_derive.py`, `CommonUtilities.wl`, `ExportJSON.wl`, `json_loader.py`,
  `_eval_utils.py`, `coefficients.py`.
- **Tabulated background interface** — ingest `a(η), H(η), x_e(η), g(η)` from CAMB rather
  than requiring an analytic metric string (blocker 5). `ModifiedCAMB` already does the
  reverse direction; read it first.
- Fix blocker 6 Wolfram-side at `ExportJSON.wl:1638`, never by patching Python
  (`.claude/rules/wolfram.md`). Fix the `_validate.py:181-216` volume-element false positive.
- **Conformal-weight analysis as a solver-selection rule**, computed from the spec, general
  path always the fallback:

```text
conformal-weight analysis of derived spec
  ├─ sector conformally invariant in D=4? → a(η) provably drops out; evolve on the flat
  │                                         slice. EXACT, not an approximation.
  └─ otherwise                            → general FRW: a(η) as time-dependent coefficient
```

  A worked instance: `√−g·ε^{μνρσ}V_μA_νF_ρσ = ε̃^{μνρσ}V_μA_νF_ρσ` is metric-independent,
  and Maxwell is conformally invariant in exactly D=4. Does **not** apply to
  gravitons/tensor modes (`a''/a`), massive torsion modes, or derivative couplings. Those
  exercise the general path — which is why it must be first-class. Side-effect: the
  conformal case is a cheap exact test the general path must reproduce to machine precision.

Success: reproduce de Sitter analytics; reproduce the conformal case to machine precision
via the general path.

### WS3 — Solver research: time-dependent linear systems *(D5; the hard one)*

The modal solver is gone (blocker 1). Replace `expm(M·t)` when `M = M(η)`, for O(10²–10³)
k-modes, blocks of 5–10, huge dynamic range in η, inside an inference loop.

| # | Method | Buys | Cost |
|---|---|---|---|
| 1 | **Exponential midpoint** `expm(h·M(η+h/2))` | correct physics immediately, 2nd order; near one-line change | trivial |
| 2 | **4th-order Gauss–Legendre Magnus** `exp(½h(M₁+M₂) − (√3/12)h²[M₁,M₂])`, batched over k | structure-preserving, fixed-cost, vectorizes; reduces *exactly* to `expm(Mt)` for constant M — generalizes the current engine rather than discarding it | 2 M-evals + 1 commutator + 1 `expm`/step |
| 3 | **Adiabatic/WKB regime switching** for `kη ≫ 1` | where the real speedup is | the research problem |
| 4 | **Piecewise-analytic transfer matrices** | reduces the ODE to repeated matrix multiplication | highest strategic value |
| 5 | Emulation | only if 1–4 insufficient | last resort |

**Critical:** Magnus converges only for `∫‖M‖ds < π`, and `‖M‖ ~ k`, so **Magnus alone
still resolves every oscillation.** It fixes time-dependence, not oscillation. Step 3 is
not optional.

**Step 3 is a genuine gap.** `oscode` (Agocs, Handley, Lasenby, Hobson — arXiv:1906.01421)
and `riccati` (arXiv:2212.06924) solve **scalar** 2nd-order ODEs; ours is matrix-valued.
**No off-the-shelf matrix RKWKB solver exists.** Matrix generalization = instantaneous
eigenbasis for the adiabatic part + Magnus for the non-adiabatic residual; template is
neutrino oscillations in matter (arXiv:0803.1967). No published application of Magnus to
cosmological perturbations was found either. Either is publishable.
**Step 4 is Handley's own work:** Haddadin & Handley, arXiv:1809.11095.

**Batching:** stacking N k-modes into one adaptive ODE is an **anti-pattern** — the error
controller forces every mode onto the worst mode's step. Batching is correct only for
fixed-step structured methods (another argument for Magnus). `tidal/solver/modal_jax.py`
already exists. `scikit-sundae` exposes only SUNDIALS' **serial** N_Vector.

**Budget:** CAMB is <1 s per `C_ℓ`; inference needs 10⁴–10⁵ Boltzmann calls ⇒ ~7 CPU-hours.
10× slower is fine; **100× is fatal** — which is what happens without a tight-coupling
analogue (CLASS: 1069 s → 19.4 s with TCA). Budget for implementing one.

### WS4 — Observables

- **Vector observables** in the measurement dispatcher (blocker 3).
- Line-of-sight sources `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`, then
  `C_ℓ^{XY} = 4π ∫(dk/k) Δ²_R(k) Δ_ℓ^X Δ_ℓ^Y`.
- Per-rung observables from the ladder. For O4, note post-processing rotation of `C_ℓ` is
  exact **only** for constant, isotropic, frequency-independent `β`; otherwise the rotation
  goes inside the LOS integral (arXiv:2209.07804).
- Cross-validate against `nanoCMB`, CAMB, and (for O4) `class_rot` (arXiv:2111.14199).

### WS5 — Cobaya extension *(D3, D9)*

- A `Theory` subclass **chaining off CAMB** (`{"CAMBdata": None}` →
  `provider.get_CAMBdata()`, which exposes `a(η)`, `H(η)`, `x_e(η)`, `g(η)`,
  transfer/source functions), wrapped in our own `get_X` accessors so no likelihood
  touches `CAMBdata`. **Built against latest upstream CAMB + Cobaya** — consistent with
  the CAMB policy above; the `slegner/*` forks are consulted only for the optional
  tabulated-background hook (H1: inherit vs re-apply cleanly).
- `get_helper_theories()` to split slow η-integration from fast post-processing — exactly
  how Cobaya's own CAMB wrapper splits `CambTransfers` from `CAMB`.
- **Packaging: Cobaya uses plain dotted-path `importlib`, not setuptools entry points.**
  Users write `tidal.cosmology.TidalTheory:` in YAML; defaults in `TidalTheory.yaml` beside
  the module. Template: `simonsobs/cosmopower_cobaya`. Add a `cosmo` extra.
- Samplers come from Cobaya (PolyChord included) — D9. Given the large speed hierarchy
  (ΛCDM params trigger CAMB; new-sector params trigger only the cheap spectator ODE), use
  Cobaya's fast/slow blocking, and prefer **dragging** over oversampling (recommended for
  large hierarchies with fast/slow degeneracies, e.g. `r` vs `ν`).
- Likelihoods: standard Planck/ACT for O0–O2. For O4 there is **no drop-in Cobaya
  birefringence likelihood**; escalate (a) Gaussian prior on published `β`, (b) fork
  `LilleJohs/cosmic-birefringence-planck-act` (MIT, reusable `.npz` spectra+covariances),
  (c) SPT-3G BB lite for anisotropic (no `α` degeneracy). The `α`–`β` miscalibration
  degeneracy is the central experimental difficulty there.

### WS6 — Numerical polology *(independent; D6)*

Build our version in the new package, **adapting Barker's numerical code directly** —
`psalter.tar.gz` (on disk; a version of the numerical psalter) and the public
`SupplementalMaterials-2607` release, guided by the algorithm in arXiv:2606.30785. Author
permission is explicit; copy freely, record provenance in docstrings, settle attribution
with the supervisor at publication. Cross-check our adaptation against the supplementary-
materials version and reproduce the Lin–Hobson–Lasenby inequalities. GH #360 has a
complete unexecuted plan (`~/.claude/plans/the-future-of-this-polished-crane.md`) — start
there, noting blocker 7. Minkowski-only is **correct and sufficient** (supervisor).
Purpose: screen vacuum-sick theories before spending compute — exactly TorC's use of
PSALTer.

## Session handoffs

This session orchestrates. Each handoff is self-contained, read-only unless stated, and
returns a written artifact this session integrates.

| ID | Session task | Returns |
|---|---|---|
| **H1** | **Read TorC end to end** (`literature/2507.09228/paper_Qtorsion.tex`) plus `github.com/slegner/CAMB` and `github.com/slegner/cobaya` diffs — CAMB against upstream `cmbant/CAMB`; cobaya against both its parent fork `AdamOrmondroyd/cobaya` and upstream `CobayaSampler/cobaya`. What exactly was patched, how `ρ_Λ(a)`/`P_Λ(a)` are threaded, how Cobaya was modified, how PolyChord was driven, what the priors/likelihoods were. Record TorC's published `ρ_Λ^eff(a), P_Λ^eff(a)` formulas (needed as O1 input tables; their bi-scalar-tensor derivation is reference material only for the *optional* minisuperspace feature — not our method). Pull the Zenodo chains listing. | `docs/cosmology/torc_pipeline_audit.md` — the inheritance decision for O1 |
| **H2** | **Observable-ladder feasibility.** For each of O1–O4: what must exist, known-answer validation target, rough cost, and what could invalidate it. Include GH #209 + `literature/2302.08186/` for O3 and `literature/2312.17636/` for the spectral-distortion route. Recommend the ordering. | `docs/cosmology/observable_ladder.md` |
| **H3** | **Solver design study (WS3).** Deepest task. Matrix-valued WKB/adiabatic generalization; Magnus convergence vs oscillation; what CAMB/CLASS/DISCO-EB actually do; prototype-level pseudocode and a benchmark protocol. Decide: patch CAMB-Fortran vs DISCO-EB vs own. | `docs/cosmology/solver_design.md` + recommendation |
| **H4** | **Repo reshape design (WS1).** Design the ideal package layout for a Cobaya extension **from the target backwards**, then map the existing code onto it — explicitly not constrained to preserve the current structure (D3). Include the `SimulationConfig` spec, migration order, and — critically — what of `tidal/inference/` survives Cobaya adoption versus is deleted as redundant (D9). | `docs/cosmology/repo_reshape.md` |
| **H5** | **Literature acquisition.** Download to `literature/` per `.claude/rules/literature.md`: `astro-ph/9603033`, `astro-ph/9506072`, `2602.23466` (nanoCMB), `2311.03291` (DISCO-EB), `2606.30785`, `1809.11095`, `1906.01421`, `2212.06924`, `0810.5488` (`2302.08186` is already local — verified), plus O4 set (`2202.13919`, `0908.0629`, `2011.11254`, `2205.13962`, `2608.06480`, `2209.07804`, `2111.14199`). Update `docs/references.md`. **Writes files.** | populated `literature/`, updated references |
| **H6** | **WS6 polology**, from the existing #360 plan. Independent of the rest. Note: the supervisor personally confirmed the expression-swell/numerical-only argument in the meeting. | `tidal/spectrum/` design doc |
| **H7** | **Spectator-route scope investigation** — **DONE (this session)**. Findings integrated above: niche empty (near-misses enumerated: xPand, CppTransport, hi_class, SymBoltz.jl); torsion-as-spectator nonexistent; double-expansion consistency (Cembranos); `ΔN_eff` validity criterion; PGT-background consistency risk; strict-spectator scope boundary; `N_eff`-vs-distinctive-channel strategy. Remaining action: write the findings up as `docs/cosmology/spectator_route.md` when WS0 executes. | `docs/cosmology/spectator_route.md` (from this session's results) |

H1, H2, H5 first (H2 depends on H1). H3 and H4 can run in parallel. H6 anytime.

## Verification

- **WS1:** suite green, behavior preserved for what survives; no `tidal.cli` imports in
  `tidal/inference/` or `tidal/measurement/`; `uv run pyright` clean.
- **WS2:** de Sitter matches analytics; the conformal case agrees between fast and general
  paths to machine precision; `tidal validate` passes on an FRW spec.
- **WS3:** each rung agrees with the previous to a stated tolerance; constant-`M` limit
  reproduces `expm(Mt)` to machine precision; per-call timing recorded at every rung
  against the ~82 ms baseline and the <1 s CAMB reference.
- **WS4/O0:** pass-through Cobaya run reproduces CAMB's `C_ℓ` and a standard ΛCDM
  posterior; any TIDAL-computed piece agrees with CAMB/`nanoCMB` sub-percent for
  `2 ≤ ℓ ≤ 2500`.
- **WS5/O1:** **reproduce the published TorC posterior** from the Zenodo chains using
  TorC's published `ρ_Λ(a)`/`P_Λ(a)` as input tables through the general
  tabulated-background feature. This is the headline end-to-end plumbing gate. (If the
  optional minisuperspace extension is built, its derived functions must reproduce TorC's
  published formulas symbolically and leave this posterior unchanged.)
- **WS2 (spectator mode):** the new sector's FRW-derived linear EOM reduce to the
  Minkowski-derived ones as `a → const` (machine precision); every run artifact carries
  the validity flags (`ρ_new/ρ_γ` vs the `ΔN_eff` bound, `P_max ≪ 1`, amplitudes `≪ 1`)
  plus the growth-impact monitor and the PGT-background consistency check (background-EOM
  residual on the CAMB background — designed fresh; the #477 work proves the concept, its
  code is not carried over).
- **WS6:** reproduces the three Lin–Hobson–Lasenby inequalities exactly and
  `psalter.tar.gz` pole masses/residues numerically.

All local (D4).

## What this session actually produces (execution deliverables)

This is the **orchestrating session**. On approval, executing this plan means producing
these concrete artifacts — including the **actual, ready-to-launch handoff prompts**, in
the campaign pattern:

1. **`docs/COSMOLOGY_PROGRAM.md`** — the program document: this plan's content
   (architecture, ladder, workstreams, decisions D1–D9, H7 findings) as the durable
   in-repo record.
2. **`docs/cosmology/handoffs/H1.md` … H6.md** — one self-contained prompt file per
   handoff, written in full: context the target session needs, precise task, inputs
   (paths, URLs, issue numbers), expected artifact and its location, and success criteria.
   Each is launchable by pasting into a fresh session (or via the Agent tool) with no
   further context from this one.
3. **`docs/cosmology/spectator_route.md`** — H7's findings written up (research already
   done this session).
4. **GitHub issue tree** — one tracking issue per workstream (WS0–WS6), linked; updates
   to #209 (→ O3) and #43 (answered by this pivot); a comment on #360 (WS6 supersedes its
   scope decision); #477's halted state recorded on its issue/branch.
5. **Memory updates** — project memory for the new direction (this pivot survives session
   loss), MEMORY.md index line, and `bash .devcontainer/scripts/sync-claude-memory.sh
   backup`.

**This session does NOT launch the handoffs (user directive).** It produces the prompt
files and stops there; the user hands each prompt to a separate session. This session
remains the orchestrator: it receives the artifacts those sessions produce, integrates
them into the program doc, and keeps the plan/issue tree current.

## Immediate next step

**H7 is done** (completed during planning; findings integrated above). Produce deliverables
1–5. The user then dispatches **H1**, **H2**, **H5** to separate sessions; this session
integrates what comes back.
