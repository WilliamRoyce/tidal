# TIDAL → Cosmology Program

**Status:** ACTIVE (program started 2026-08-29, after supervisor meeting 2026-08-28)
**Tracking:** umbrella #488 · WS0 #489 · WS1 #490 · WS2 #491 · WS3 #492 · WS4 #493 ·
WS5 #494 · WS6 #495 (anchors: #209 → O3, #43 answered, #360 scope updated, #477 halted)
**Orchestration:** one orchestrator session holds this document current; workstream
sessions are dispatched by the user from the prompts in `docs/cosmology/handoffs/`.

## Goal

Integrate a candidate Lagrangian's perturbations over the history of an expanding
universe, produce CMB observables, and do genuine Bayesian inference against measured
data — packaged as a **Cobaya extension** so others can test their own Lagrangians
against real likelihoods. This replaces the prior mode of work (coupling-space surveys on
flat Minkowski scored against synthetic objectives).

## The brief, decoded

The direct predecessor is the group's own paper, local at `literature/2507.09228/`:
**Legner, Handley & Barker, "Alleviating the Hubble tension with Torsion Condensation
(TorC)"** — stack: PSALTer particle spectrum → modified CAMB → Cobaya → PolyChord,
against Planck 2018 + SH0ES. It states exactly the limitation this program lifts:

> "this analysis modifies the background expansion in CAMB … while the perturbation
> equations remain those of standard ΛCDM" … "cosmological perturbation theory for TorC
> will be developed in future work."

Public assets from that work (all verified 2026-08-29):

- `ModifiedCAMB` → <https://github.com/slegner/CAMB> (fork of `cmbant/CAMB`; reads
  tabulated `ρ_Λ(a)`, `P_Λ(a)` — the `w(a)` interface develops poles when `ρ_Λ` changes
  sign)
- `ModifiedCobaya` → <https://github.com/slegner/cobaya> (forked via
  `AdamOrmondroyd/cobaya`)
- Chains/supplementary → doi:10.5281/zenodo.15866507

TorC used the **symbolic** Wolfram PSALTer (arXiv:2406.09500 + 2506.02111) — sufficient
for one theory analyzed once. Sampling arbitrary Lagrangians needs the **numerical**
polology route (arXiv:2606.30785): "symbolic computation scales poorly due to expression
swell; the only avenue is numerical" — confirmed verbatim by supervisor Barker in the
meeting.

## Central architecture: the spectator (test-field) route

**The background is always the established ΛCDM one, supplied by CAMB. The new sector's
perturbations are evolved on top of it as test fields, and their observable imprints (on
gravitons, photons, polarization, lensing) are what we compare to data.**

The spectator limit means: **(a)** the new sector does not affect the background
expansion (negligible in the Friedmann equation); **(b)** it does not gravitationally
disturb the standard perturbations (negligible in the Einstein constraints); **(c)**
everything stays linear. Two clarifications:

- "Standard perturbations don't back-react on the background" is standard first-order
  cosmological perturbation theory (`δρ/ρ ~ 10⁻⁵`), used by everyone — not our
  assumption. Our *additional* assumption is (a)+(b) for the **new** sector only.
- We *do* modify the propagation of standard quanta (photons, gravitons) — through the
  **explicit coupling terms in the Lagrangian** (torsion–photon, torsion–graviton
  mixing), kept at linear order. Condition (b) forbids only *gravitational* sourcing via
  the new sector's stress-energy. Coupling-mediated ≠ gravity-mediated; the former is the
  entire point.

The consistency of coupling-without-backreaction is a clean double expansion (stated
verbatim in Cembranos et al., arXiv:2302.08186, local): order 0 in perturbations →
background equations, discarded and replaced by CAMB's solution (the only place the
spectator assumption enters); order 1 → tadpoles, vanish on-shell; order 2 → the full
quadratic action *including all mixing terms*, kept whole.

### Reachable observables (and the boundary)

1. **Propagation (reachable):** standard quanta travel through the new sector's
   background/couplings, altering speed, phase, damping, polarization — no energy moves.
   Birefringence, modified GW friction/dispersion.
2. **Conversion (reachable):** quanta oscillate between sectors (Gertsenshtein
   graviton↔photon), moving a small energy fraction — CMB spectral distortion, V-modes,
   radio excess. The tiny converted fraction *is* the measurement.
3. **Gravitational sourcing (NOT reachable in the strict limit):** the new sector's
   `δρ, σ` in the Einstein constraints creating anisotropy the way CDM does — that is
   dropping assumption (b), the axionCAMB-style full-component route. A deliberate later
   extension, never silently blended in.

### Validity enforcement (supervisor-flagged; first-class requirement)

The literature asserts spectator validity in one sentence and never enforces it.
We enforce it numerically, per run (honest-flags style, like `gauge_certificate`):

- `ρ_new/ρ_γ` against the `ΔN_eff ≲ 0.1` bound (Domcke & Garcia-Cely arXiv:2006.01161,
  local);
- conversion probability `P_max ≪ 1`;
- amplitudes `|h|, |f| ≪ 1`;
- **growth-impact monitor**: per mode, the ratio of new-sector to standard-sector source
  terms in the Einstein constraints — "would these perturbations have affected the growth
  we froze?"

**Silent-failure risk to check:** the CAMB background solves *Einstein's* equations; our
quadratic action is PGT. Consistency requires the PGT background field equations be
satisfied by (FRW, `T̄ = 0` or tracking torsion) to needed accuracy — see
`literature/2003.02690/` for the tracking/frozen solutions. The check is a
**background-EOM residual** on the CAMB background, designed fresh for the new
architecture (the concept is proven by the #477 work; its code is not carried over).

### Division of labor per likelihood call

| Computed by | What | Recomputes when |
|---|---|---|
| **CAMB (unchanged)** | ΛCDM background `a(η)`, `H(η)`; thermal history `x_e(η)`, `g(η)`; all standard-sector perturbations; standard `C_ℓ`, lensing, matter power | ΛCDM params move (slow block; cached otherwise) |
| **Our solver** | **the coupled block**: new-sector modes `δT(k,η)` *plus the standard modes they directly couple to* (tensor `h`; photon polarization) — small per-`k` linear ODE system, coefficients involving CAMB's `a(η)`, `H(η)` | new-sector coupling moves (fast block) |
| **Nobody** | nothing ΛCDM is re-derived or re-solved by us — ever | — |

No existing package can evolve *our* equations — CAMB/CLASS hard-code the ΛCDM species.
The one genuinely new numerical task (small, k-parametrized, time-dependent-coefficient
linear systems) is TIDAL's existing modal competence generalized to `M = M(η)`.

**Worked chain (tensor channel):** primordial `A_t, n_t` → per-`k` integration of the
coupled `(h, torsion)` system over `η` with CAMB's `a(η)` → transfer function →
line-of-sight projection `Δ_ℓ(k) = ∫dη S(k,η) j_ℓ(k(η₀−η))` →
`C_ℓ^BB = 4π∫(dk/k) Δ²_t(k)|Δ_ℓ(k)|²` → BICEP/Keck/LiteBIRD likelihood. Photon channel
analogous: accumulated rotation/conversion applied to CAMB's `C_ℓ` → EB/V data.

**`a(η)` handling:** derive the equations **symbolically with `a(η)`, `H(η)` as
unspecified background functions** (once per theory — this is how expansion enters:
Hubble-friction `∝ H`, `a²`-dependent masses); evaluate coefficients numerically from
CAMB's table at solve time (per call). Analytic-metric mode retained for validation cases
(de Sitter). **Working time variable: conformal time `η`** (CAMB's internal variable);
`t ↔ η ↔ z` conversions at the CAMB interface.

**Cobaya integration:** our component is a `Theory` class (Cobaya's pluggable prediction
component: `get_requirements()` / `calculate()` / `get_X()`), chaining off CAMB
(`{"CAMBdata": None}` → `provider.get_CAMBdata()`), returning either replaced transfer
functions or conversions/rotations applied to CAMB output. Downstream likelihoods see
standard products. Packaging is plain dotted-path import (no entry points); template:
`simonsobs/cosmopower_cobaya`. Fast/slow blocking with **dragging** (large speed
hierarchy: ΛCDM params trigger CAMB, coupling params trigger only our ODE).

**CAMB policy:** target **latest upstream CAMB** by default. The tabulated-background
feature is an *optional hook*, off by default. **Settled by H1 (§7, decision R2, GH #498):
re-apply the patch ourselves on our own fork off the upstream `2.0.3` tag — do not inherit
`slegner/CAMB`.** That fork is 180 commits behind upstream, carries unrelated post-paper
work, a breaking public-signature change, a debug `print *`, and a stray settings file;
pinning it would also contradict the latest-upstream policy. Upstream's own
`set_w_a_table` is not an alternative — it is the pole-prone route the paper had to avoid.
The patch itself is small: one optional `P` output on `BackgroundDensityAndPressure` so
pressure comes from a spline instead of being rebuilt as `w·ρ`; the rest of the fork's
file count is signature churn from that one argument.

## Decisions (user, 2026-08-29)

| # | Decision |
|---|---|
| D1 | Target CMB observables via genuine perturbation evolution. The point is the **general engine**; birefringence is one supervisor-stated end goal among several. |
| D2 | No paper target fixed; order workstreams by dependency. |
| D3 | Reshape the repo around the Cobaya-extension goal — restructuring **should be done** where right; no attachment to the current layout. |
| D4 | No HPC without explicit permission. Local only. |
| D5 | Optimization is a standing first-class concern — adoption depends on it. |
| D6 | PSALTer-numerical: build our version **taking heavy, deliberate inspiration from Barker's code** (`psalter.tar.gz` + SupplementalMaterials-2607 + arXiv:2606.30785). Author permission is explicit; copy freely; provenance in docstrings; attribution settled with supervisor at publication. Not clean-room. |
| D7 | The solver targets the general case; cancellations are auto-detected fast paths (as `can_use_modal` does), never scoping decisions. |
| D8 | Don't overspecialize early — CS/birefringence is one instance; engine and plan stay general. |
| D9 | Migrate inference into the Cobaya ecosystem (Cobaya ships PolyChord); much of `tidal/inference/` is then superseded; redesign deliberate. |

## Observable ladder

Ordered by pipeline-validation value; each rung has a checkable answer before the next
adds unknowns.

| # | Observable | Exercises | Validation |
|---|---|---|---|
| **O0** | ΛCDM `C_ℓ`, new physics off (pass-through mode) | plumbing: Cobaya wiring, CAMB products, sampler, likelihood | reproduces CAMB `C_ℓ` + standard ΛCDM posterior; any TIDAL-computed piece sub-percent vs CAMB/nanoCMB, `2 ≤ ℓ ≤ 2500` |
| **O1** | **Fixed-table pass-through** (H1 §R1): one `(a, ρ, P)` table generated offline at a chosen `(Ω_Λ, ϖ_r)` from TorC's published formulas, fed through our optional tabulated-background hook. Table generation is one-off data prep, **not** a package feature; no TorC physics inside the package | the hook + Cobaya wiring + CAMB seam, end to end | drive CAMB to the same `H(a)` and `C_ℓ` as a reference build. **Not** a posterior reproduction — see the note below. **First target.** |
| **O2** | **Spectator perturbations on ΛCDM** — first: modified tensor/GW propagation (friction `ν`, mass `μ`, speed `c_T` derived from the Lagrangian) → tensor transfer → B-modes | the spine: TIDAL FRW-derived equations + new solver, strict spectator | TorC explicitly deferred this; first genuinely new result. Validity flags per run. Cembranos (local): `μ ≤ 10⁻³³ eV` does nothing; friction `ν` matters strongly |
| **O3** | **Gertsenshtein mixing on FRW** — graviton↔photon with Hubble damping | time-dependent solver on the thesis's own physics; **core goal** (thesis couplings) | GH #209; Cembranos arXiv:2302.08186 (local); reduces to flat-space result as `H → 0` |
| **O4** | **Cosmic birefringence (E→B)** + V-modes | parity-odd sector + polarization observables | `β = 0.277° ± 0.057°` (4.8σ, joint ACT+Planck, arXiv:2608.06480); distinctive channels beat the `N_eff` bound where broadband conversion does not |

**Two O4 prerequisites, found 2026-08-30 — neither is a lookup, both are derivation work:**

- **Frequency scaling has no single answer.** It is *per-operator*, and TIDAL's own theory
  documents do not derive it: `general_quadratic_lagrangian.tex` contains no
  frequency-scaling discussion, and its CS section states only that `χ^CS_2` sources
  magnetic helicity `A·B`. The Das et al. `ξ₁ ∇T×F̃` structure gives `α = 2ξ₁p²T₁t` (`ν²`);
  an induced-axion `θFF̃` gives `ν⁰`. These are *different sectors* in the enumeration's
  language (`ζ̃₁₋₆` vs `d₁₈`/`χ^CS`), so "the torsion coupling's scaling" is not a
  well-posed question. Scope it as a derivation task per operator (**GH #499**).
- **No example implements the Chern–Simons couplings.**
  `examples/torsion_gertsenshtein/theory_parity_odd.toml` explicitly *defers* `cs1`–`cs3`
  ("require bare `A_μ`, special handling"). O4 therefore needs new derivation work, and
  the bare-`A_μ` gauge handling is an unsolved problem, not a configuration step (**GH #499**).
  (That file's `description` field claimed CS *was* included, contradicting its own
  header; corrected 2026-08-30.)

O0 is the validation gate, then O1. Strategic note (H7): CMB bounds from broadband
conversion are typically weaker than `N_eff` unless the signal is spectrally narrow or
structurally distinctive — V-modes and birefringence are exactly that, justifying the
polarization emphasis on physics grounds.

## The niche (H7 investigation, 2026-08-29)

**"Arbitrary Lagrangian → spectator-sector perturbations on fixed ΛCDM → observables →
Cobaya likelihood" is empty**, verified by enumerating near-misses: xPand/xAct (symbolic
only), CppTransport/PyTransport (inflation only, arXiv:1609.00380/1), hi_class (Horndeski
gravity sector, full backreaction, arXiv:1909.01828), SymBoltz.jl (equations in, not
actions, arXiv:2509.24740 — possible future backend). **Nobody has treated torsion as a
spectator on FRW at all**; arXiv:2506.17017 and TorC both signpost cosmological
perturbation theory as future work. Full findings: `docs/cosmology/spectator_route.md`.

## Workstreams

Order: **WS0 → WS1 → WS2 ∥ WS3 → WS4 → WS5.** WS6 independent.

**FRW is blocked in three independent places, each owned by a different workstream** (found
2026-08-30). No single fix clears it, and they must be tracked separately:

| # | Blocker | Owner |
|---|---|---|
| a | Modal solver refuses non-trivial `volume_element` and any `time_dependent` term (`tidal/solver/modal.py`) | WS3 (#492) |
| b | Hamiltonian/energy export filters `t` (`ExportJSON.wl:1638`; `_energy.py` hardcodes `t=0.0`) | WS2 (#491) |
| c | **Conversion measurement is energy-ratio-based**, so it inherits (b)'s `t=0.0` bug *and* the `P_max`-vs-`P_final` distinction | WS4 (#493) |

(c) is the dangerous one: it would silently corrupt an O3 number rather than fail loudly.

**Cost basis:** use the **v0.33.9 measured table** in `docs/tex/derivation_performance.tex`.
The per-theory derivation-timing headers in the TOMLs are declared untrustworthy there
(§lines 108–114) — treat them as ceilings, not estimates.

- **WS0 — Research & scoping** (no code): handoffs H1–H5; integration-target decision
  (patch CAMB Fortran vs DISCO-EB [arXiv:2311.03291 — from-scratch differentiable
  Einstein–Boltzmann in JAX, per-mille vs CAMB, built to be extended] vs own coupled-block
  solver chained to unmodified CAMB, the default assumption). Language is not a barrier
  (user); decide on architecture, not language.
- **WS1 — New package, strangler-fig**: the existing framework is **legacy, not a
  template**. New package designed whole, from the goal backwards (H4); **new code never
  imports old code — no adapters**; useful capabilities are *fully ported* (redesigned,
  with docstrings/why-comments/issue references traveling along; original kept as test
  oracle). Two fully separated packages/CLIs during transition; old `tidal` stays for
  cross-checks; legacy retired deliberately once covered.
- **WS2 — General expanding background**: symbolic derivation with unspecified `a(η)`,
  `H(η)`; conformal time as a first-class coordinate; CAMB-table coefficient evaluation;
  conformal-weight fast-path detection (general path always the fallback; the conformal
  case doubles as a machine-precision test); fix the time-dependent Hamiltonian/energy
  export bug (`ExportJSON.wl:1638`, `_energy.py` `t=0.0`) Wolfram-side.
- **WS3 — Solver research** (the hard one): replace `expm(M·t)` for `M(η)`. Ladder:
  (1) exponential midpoint → (2) 4th-order Gauss–Legendre Magnus, batched over k
  (reduces exactly to `expm(Mt)` for constant M) → (3) adiabatic/WKB regime switching for
  `kη ≫ 1` — **the research gap**: oscode/riccati are scalar-only; no matrix RKWKB solver
  exists; template = neutrino oscillations in matter (arXiv:0803.1967) → (4)
  piecewise-analytic transfer matrices (Haddadin & Handley arXiv:1809.11095) → (5)
  emulation, last resort. Magnus alone still resolves every oscillation
  (`∫‖M‖ds < π`, `‖M‖ ~ k`) — step 3 is not optional. Batching only for fixed-step
  structured methods; never stack modes into one adaptive ODE. Budget: 10× slower than
  CAMB is fine, 100× is fatal (the tight-coupling lesson: CLASS 1069 s → 19.4 s).
- **WS4 — Observables**: vector observables; line-of-sight sources; per-rung observables;
  post-processing rotation exact only for constant/isotropic/frequency-independent
  effects, else inside the LOS integral (arXiv:2209.07804); cross-validate vs nanoCMB
  (arXiv:2602.23466), CAMB, class_rot (arXiv:2111.14199).
- **WS5 — Cobaya extension**: Theory class as above; latest upstream CAMB + Cobaya;
  likelihoods per rung (no drop-in birefringence likelihood exists — escalate: Gaussian
  prior on published `β` → fork `LilleJohs/cosmic-birefringence-planck-act` (MIT) →
  SPT-3G BB lite for anisotropic).
- **WS6 — Numerical polology**: adapt Barker's numerical code (D6) into the new package;
  cross-check vs SupplementalMaterials-2607; reproduce Lin–Hobson–Lasenby inequalities;
  Minkowski-only is correct and sufficient (spectrum screens in vacuo; TorC used the same
  split). Start from the #360 plan (`.claude-plans-backup/the-future-of-this-polished-crane.md`),
  noting its Minkowski-first *project-wide* scope decision is superseded.

## Parked

`feat/ws2-localized-path-audit` / #477 arc: **already halted** (issues arose making it
not possible; likely unrelated to the new direction). Parked as-is with state recorded on
the issue. Also parked: Phase E (T6/T8 rescue), Phase A-γ, NEXT_PHASES G/H/I, Wolfram CI
(#69). Nothing gates on any of it.

## Handoffs

Prompt files in `docs/cosmology/handoffs/` — self-contained, dispatched by the user to
separate sessions; **this orchestrator session does not launch them**. H7 was executed
during planning; its findings are in `docs/cosmology/spectator_route.md`.

| ID | Task | Artifact |
|---|---|---|
| H1 | TorC pipeline audit — **DONE** (2026-08-30). Settled R1 (O1 = fixed-table pass-through, above) and R2 (re-apply the CAMB patch, GH #498). Also for WS5: the provider→consumer wiring template (§3.1), the `planck_clik` NaN guard as a flagged rejection path (§3.3), and a caution that stock Cobaya's PolyChord does not give correct evidences under non-uniform priors without the Ormondroyd patch (§3.2). For the spectator flags: `ΔN_eff` (§1.7) is already parameterized for the torsion sector | `docs/cosmology/torc_pipeline_audit.md` ✅ |
| H2 | Observable-ladder feasibility (O1–O4; recommend ordering) | `docs/cosmology/observable_ladder.md` |
| H3 | Solver design study (matrix WKB; CAMB-Fortran vs DISCO-EB vs own; benchmark protocol) | `docs/cosmology/solver_design.md` |
| H4 | New-package design (from the goal backwards; config spec; CLI separation; inference fate) | `docs/cosmology/repo_reshape.md` |
| H5 | Literature acquisition — **DONE** (2026-08-30). 20/20 fetched and title-verified; `literature/README.md` now tracked and auto-generated by `scripts/bibaudit/index_literature.py` (GH #497) | populated `literature/` + `docs/references.md` ✅ |
| H6 | Numerical polology design (from #360 plan + Barker's code) | spectrum-package design doc |
| H7 | Spectator-route scope — **DONE** | `docs/cosmology/spectator_route.md` |

Dispatch order: **H1 ✅ and H5 ✅ done. H2 next** — its dependency (H1's audit) is now
satisfied, and H1 left it two specific inputs: weigh O1 purely as a gate (it produces no
new physics as scoped), and treat `ΔN_eff` as a ready-made spectator-validity input.
Then H3 ∥ H4. H6 anytime.

## Verification gates

- **WS1:** suite green; no old-code imports from new code; pyright clean.
- **WS2:** de Sitter analytics reproduced; conformal fast path ≡ general path to machine
  precision; FRW-derived EOM → Minkowski EOM as `a → const`; validity flags on every run
  artifact.
- **WS3:** each ladder rung agrees with the previous to stated tolerance; constant-`M`
  limit reproduces `expm(Mt)` to machine precision; per-call timing recorded per rung
  (~82 ms current baseline; <1 s CAMB reference).
- **WS4/O0:** pass-through reproduces CAMB `C_ℓ` and a standard ΛCDM posterior.
- **WS5/O1 (narrowed by H1 §R1):** a fixed `(a, ρ, P)` table driven through the
  tabulated-background hook reproduces a reference CAMB build's `H(a)` and `C_ℓ` at the
  same `(Ω_Λ, ϖ_r)`. **Second, free oracle:** at TorC's own fiducial point
  (`ϖ_r = 0.8`, `Ω_Λ = 0.685`) `ρ_Λ^eff` never changes sign and `|w| ≤ 0.998` — poles
  appear only at `ϖ_r > 1` — so the same table pushed through stock CAMB's
  `set_w_a_table` must agree with our patched path. Disagreement then localizes a
  re-apply bug rather than a physics one.
  **Posterior reproduction is NOT a gate**, but be precise about why (H1 correction).
  *Exact* reproduction is impossible: no Cobaya run configuration is archived, so
  `num_repeats`, the precision criterion and the seed are unrecoverable. *Statistical*
  comparison, however, **is** available — `_equal_weights.txt` is archived for all seven
  runs and `.paramnames` carries `H0`, `omegam`, `sigma8` and the `S_8` variants as
  derived parameters, so marginal posteriors are obtainable via `anesthetic`; the paper
  simply never printed a table. So a later O1b could compare contours and check `log Z`
  against the stated ±0.22 — but **not** lean on Δ`log Z`, which would carry an unknown
  systematic from the unknown `num_repeats`. The four evidences (H1 §5.3) are the only
  numbers the paper itself states.
- **WS6:** Lin–Hobson–Lasenby inequalities reproduced exactly; agreement with the
  supplementary-materials implementation on pole masses/residues.

All local (D4).
