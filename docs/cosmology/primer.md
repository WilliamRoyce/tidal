# Primer: what a CMB pipeline computes, and where TIDAL plugs in

**Purpose.** The physics orientation for the cosmology program (`docs/COSMOLOGY_PROGRAM.md`).
Written out deliberately, because cosmology is a new field for this project — the program
document says *what* we are building and *when*; this says *why the physics works that way*.
Every claim here was checked against the TorC paper (`literature/2507.09228/`) or the
Boltzmann-code literature during the 2026-08-29 planning session and refined over several
rounds of correction; H1's audit subsequently confirmed the TorC-specific parts.

Read this before `docs/COSMOLOGY_PROGRAM.md` if you are new to the program.

---

## 1. The pipeline, in stages

A Boltzmann code (CAMB, CLASS) turns cosmological parameters into a predicted CMB:

1. **Background.** Solve for the scale factor `a(η)` — how the universe expands. Driven by
   the total energy density `ρ(a)` and pressure `P(a)` of everything in it.
2. **Thermal history.** When do electrons and protons combine into neutral hydrogen, so
   photons stop scattering? Gives `x_e(η)` and the visibility function `g(η)` — the
   probability that a photon observed today last scattered at time `η`.
3. **Perturbations.** For each comoving wavenumber `k`, evolve the *lumpiness* — density
   contrasts `δρ`, velocities, anisotropic stress — of photons, baryons, CDM and neutrinos,
   coupled to one another through gravity, from the early universe to today.
4. **Projection.** Convert that evolution into what we see on the sky, via the
   line-of-sight integral `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`.
5. **Spectra.** `C_ℓ^{XY} = 4π ∫(dk/k) Δ²_R(k) Δ_ℓ^X Δ_ℓ^Y` — the angular power spectra.

**`C_ℓ` is the observable we compare to data.** Planck/ACT/SPT measure `C_ℓ^{TT}`,
`C_ℓ^{TE}`, `C_ℓ^{EE}` (and lensing); the likelihood scores our predicted `C_ℓ` against
those measurements. That comparison *is* the inference.

## 2. The two levels, and what TorC actually did

### Level 1 — background (stage 1)

TorC (the model *is* "Torsion Condensation" — a condensed background torsion field)
contributes an *effective* dark-energy density and pressure `ρ_Λ^eff(a)`, `P_Λ^eff(a)`.

"Feeding those to CAMB" means: CAMB needs the universe's total energy content to solve for
`a(η)`. Normally you hand it `Ω_b, Ω_c, Ω_Λ` and it assumes `ρ_Λ` is constant. TorC's
`ρ_Λ^eff` *varies with `a`*. CAMB's standard hook for that is the equation of state
`w(a) = P/ρ` — **but that develops poles when `ρ_Λ` changes sign**, which TorC's does. So
they patched CAMB to read tabulated `ρ_Λ(a)` and `P_Λ(a)` separately. **That patch is
`slegner/CAMB`; the "modified CAMB" is the feeding mechanism, for background modification
only.**

Changing `a(η)` shifts the sound horizon `r_*` and the angular diameter distance `D_A^*`,
hence the acoustic-peak positions `θ_s = r_*/D_A^*`, and the early ISW — all real,
measurable effects on `C_ℓ`. That is what TorC computed and compared to Planck.

### How a general Lagrangian becomes `ρ^eff(a)`, `P^eff(a)` — the effective-fluid method

*Context for the optional background feature only; the spectator spine never does this.*

"Field equations" here means the **gravitational** field equations — the theory's
generalization of Einstein's equations, obtained by varying the full action with respect to
the metric (and, in PGT, the connection); *not* the matter/torsion-field EOMs. These are the
same equations whose *linearization* gives the `h`-component equations TIDAL already
derives; here they are used un-linearized, evaluated on the FRW ansatz instead.

Evaluate them on the homogeneous-isotropic (FRW) ansatz. The **time-time (`00`) component
of that tensor equation** is the energy-constraint (Friedmann) equation and always takes the
form `3 M_p² H² = ρ_standard + (everything else)`. One *defines* `ρ^eff` as everything
beyond the GR terms, and `P^eff` similarly from the spatial-trace component. This is the
standard modified-gravity move — the new sector is disguised as a fluid CAMB already knows
how to carry.

**What TorC's "bi-scalar-tensor equivalent" is, and why it is NOT our method.** TorC
hand-derived a mapping specific to their theory: on homogeneous-isotropic backgrounds the 24
torsion components reduce to two scalars (`ϖ`, `φ`), and the PGT Lagrangian maps to an
equivalent scalar-tensor theory whose Friedmann equations are then read off. That is
precisely the **per-theory hand-crafting this project exists to eliminate** — it required
theory-specific insight, holds only for that model class, and does not generalize. Its only
role for us: *if* the optional minisuperspace feature is ever built (generic FRW reduction of
any TOML Lagrangian, no equivalent-theory tricks), TorC's published result is its validation
case. It is on no critical path. **The spectator spine needs no background derivation at all.**

### Level 2 — perturbations (stage 3)

ΛCDM contains no torsion, so there is no such thing as "ΛCDM's torsion perturbations." What
TorC did is **omit the torsion perturbations entirely and use GR's perturbation equations as
a stand-in**:

> "This work focuses on the background evolution … and does not yet incorporate cosmological
> perturbations. For the purposes of this initial study, the perturbation equations remain
> those of standard ΛCDM." … §*Fiducial perturbation theory*: "we study the background
> dynamics in isolation, **adopting GR perturbation theory as a fiducial proxy**."

So the torsion was allowed to change *how fast the universe expands*, but not to have its own
ripples that gravitate, cluster, free-stream, or carry anisotropic stress.

**Including those torsion perturbations is precisely this project's work — it is O2, the
first genuinely new result.** TorC names it as its own future work; we are building the
engine that does it, for any Lagrangian.

**Promoting to Level 2** means deriving the torsion perturbations' own equations of motion
from the Lagrangian and evolving `δT(k,η)` alongside the standard species. **The linear
perturbation equations TIDAL already derives symbolically are exactly the object that gets
routed into the Boltzmann code**: re-derived about an FRW background instead of Minkowski,
they become the new evolution equations added to stage 3, coupled to the standard species
through the Einstein equations. TIDAL is the front end that generates what the Boltzmann code
evolves. The new-sector perturbations then:

- contribute `δρ, δP` to the Einstein constraints → change how gravitational potentials evolve;
- contribute **anisotropic stress** `σ`, which alters ISW **and lensing**.
  *(Notation — these are standard objects, not new fields of ours: `φ` and `ψ` are the two
  scalar gravitational potentials of the perturbed FRW metric,
  `ds² = a²[−(1+2ψ)dη² + (1−2φ)δᵢⱼdxⁱdxʲ]`. `ψ` is the Newtonian-potential-like time-time
  perturbation, `φ` the spatial-curvature one. In GR with no anisotropic stress they are
  equal; a sector carrying `σ` splits them via `k²(φ − ψ) = 12πG a²(ρ̄+P̄)σ`, and that split
  is what shifts ISW/lensing.)*
- propagate with their own sound speed → change the driving of the acoustic oscillations;
- if parity-odd, **mix E and B polarization** — "polarizations converting into each other."

Running that through stages 4–5 is **"integrating the perturbations over the history of the
universe."**

**Why it matters — and TorC invites it directly.** The paper argues TorC is a research
*program* with several realizations (bi-scalar-tensor, "tensor bypass", full PGT) that
**share background dynamics but differ in their perturbations**: *"Different realisations may
yield distinct perturbation theories."* So the perturbations are precisely what discriminates
between them — background-only cannot. Additionally, background-only is an approximation of
unquantified validity: if torsion perturbations are not negligible, using GR's biases the
inferred parameters, and the Hubble-tension claim rests on those parameters.

## 3. What we reuse versus what we build

**The key architectural decision.** The ΛCDM sector is enormously well solved: photons,
baryons, CDM, neutrinos, recombination, the full Boltzmann hierarchy, validated across 12
independent codes to sub-percent agreement (arXiv:1709.09135). **We do not reimplement any of
it.** We take it from CAMB and add only our new sector, coupled in through the Einstein
equations. That is the entire reason to chain onto CAMB rather than write a Boltzmann code.

**Do we have to evolve the photons?** Depends on the theory, and this is what orders the
observable ladder:

- If the new physics is **gravity-sector only** (O1, O2): photon transport is unmodified
  ΛCDM — take it from CAMB untouched.
- If the Lagrangian **couples torsion to electromagnetism** (O3, O4 — the thesis's own
  physics): photon propagation *is* modified and we must intervene, either as a modified
  line-of-sight source (Tier 3a) or by evolving a modified photon sector (Tier 3b).
  **This is a core goal of the program, not an optional extension**: the thesis is built on
  the Gertsenshtein coupling, and seeing those couplings' cosmological effects is a key use
  case the supervisors want. The engine must handle modified photon sectors as first-class.

## 4. Where the new physics enters — and why the tier table matters

Stages 1–5 are a fixed pipeline we inherit. The practical question is **at which stage our
new physics gets inserted**, because that determines how much code we write and whether we
can stay in Python. It is the single biggest driver of cost, which is why it is settled
before any code. The tier table in `COSMOLOGY_PROGRAM.md` (WS0) answers it:

- **Tier 1** touches only stage 1 — a tabulated background. What TorC did, and our O1.
- **Tier 3a** adds a term to the stage-4 source function; pure Python via `camb.symbolic`.
- **Tier 3b** adds a *new evolving field* to stage 3, and requires going into the Fortran/C
  or using a code built to be extended.

**TIDAL's torsion is Tier 3b** — the hard part, and exactly what TorC deferred.

## 5. Why an in-house solver, and what a "Theory class" is

**Why an in-house solver at all, if we are plumbing into existing packages?** Because no
existing package can evolve *our* equations: CAMB/CLASS hard-code the ΛCDM species and have
no notion of a user-derived field with TIDAL-generated equations of motion. The one genuinely
new numerical task — integrating a small, k-parametrized, time-dependent-coefficient linear
system — is exactly TIDAL's existing competence (the modal engine), generalized to
`M = M(η)`. That generalization is WS3.

**What a "Theory class" is.** Cobaya organizes an inference run into pluggable components:
*Theory* components compute predictions, *Likelihood* components score them against data, and
a *sampler* drives the loop. A "Theory class" is simply a Python class implementing Cobaya's
small interface (`get_requirements()` — what it needs from other components; `calculate()` —
what it computes; `get_X()` — what it offers). **We write one; that IS the Cobaya extension.**
Ours declares that it requires CAMB's products, runs our per-k solver for the coupled block,
and returns the modified observables — either (a) a **replaced transfer function** (the tensor
mode `h(k,η)` is now the solution of our coupled graviton–torsion system instead of CAMB's
plain `h'' + 2ℋh' + k²h = 0`), or (b) a **conversion/rotation applied to CAMB's output**
(spectral-distortion amplitudes, polarization rotation of `C_ℓ`). Downstream likelihoods see
standard products and never know the difference.

### A worked chain (tensor channel), making "evolution → observable" explicit

1. Inflation supplies the primordial tensor amplitude per mode `k` (two numbers: `A_t`,
   `n_t` — sampled parameters).
2. For each `k` on a grid, our solver integrates the coupled `(h, torsion)` system from the
   deep radiation era to today in `η`, with `a(η)`, `H(η)` read from CAMB's table. The output
   is the **transfer function** — how much of the initial `h` survives (or is
   amplified/converted) by recombination and after.
3. CAMB's line-of-sight machinery projects that source onto the sky:
   `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))`, then `C_ℓ^BB = 4π∫(dk/k) Δ²_t(k)|Δ_ℓ(k)|²`.
4. A likelihood (BICEP/Keck, LiteBIRD) scores that `C_ℓ^BB` against measured B-mode data.

Same shape for the photon channel: the torsion background modifies polarization propagation
along the line of sight → rotation/conversion accumulated between recombination and today →
applied to CAMB's `C_ℓ^{EE}` → scored against `EB`/`V` data. **The CMB is literally those
photon quanta**, free-streaming since recombination; anything that alters them en route is
imprinted on what we measure. That is why propagation effects are so constraining.

## 6. The `a(η)` analytic-vs-tabulated tension, and its resolution

Two stages want `a(η)` in two different forms, and the current pipeline conflates them. The
**Wolfram derivation** stage today requires the metric as an *analytic string* (that is how
`examples/curved_spacetime/de_sitter.toml` works), but the real ΛCDM `a(η)` exists only as
CAMB's *numeric table*.

The resolution is to separate the stages: derive the equations **symbolically with `a(η)`,
`H(η)` left as unspecified background functions** — the coefficients then *contain* symbolic
`a`, `H`, which is exactly how the expansion details enter (Hubble-friction terms `∝ H`, mass
terms `∝ a²`) — and at **solve time** evaluate those coefficients numerically from the CAMB
table. Symbolic derivation once per theory; numeric evaluation per call. The analytic-metric
mode stays for validation cases (de Sitter) where exact solutions exist.

**Working time variable: conformal time `η`** — CAMB's internal variable, and the one in
which the perturbation equations take their cleanest form (`ds² = a(η)²[−dη² + dx²]`). The
existing derivation pipeline uses coordinate time `t`; supporting `η` is part of WS2.
Conversions `t ↔ η ↔ z` happen at the CAMB interface, which supplies all three.

## 7. Where TIDAL fits, and why Cobaya

TIDAL already derives linear perturbation equations from a Lagrangian symbolically — that is
its entire purpose. **That is the tool.** TorC would otherwise hand-derive them per theory;
automating it *for any Lagrangian* is what makes this a general extension rather than a
one-off.

And Cobaya is the payoff. Once our component supplies `C_ℓ` (and lensing potential, matter
power spectrum, …) through Cobaya's standard interface, **every likelihood in that ecosystem
becomes available for free** — Planck, ACT, SPT, DESI/BAO, supernovae, weak lensing. We build
the theory tool once; testing it against many different datasets is then configuration, not
code.

Equally important for impact: it makes the tool **drop into people's existing setups**. A
Cobaya component is something the community already knows how to install, configure and cite,
which is what turns a research code into an adopted (and cited) one. That is the argument for
chaining into Cobaya rather than building a bespoke pipeline.

---

## See also

- `docs/COSMOLOGY_PROGRAM.md` — the operational record: decisions, observable ladder,
  workstreams, verification gates.
- `docs/cosmology/spectator_route.md` — the spectator/test-field approximation: where it is
  used, its validity criteria, what it can and cannot reach, and the empty-niche argument.
- `docs/cosmology/torc_pipeline_audit.md` — H1's audit of the TorC paper, forks and archive.
- `docs/cosmology/handoffs/` — per-workstream session prompts.
