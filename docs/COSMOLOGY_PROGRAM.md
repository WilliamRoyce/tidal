# TIDAL → Cosmology Program

**Status:** ACTIVE (program started 2026-08-29, after supervisor meeting 2026-08-28)
**Tracking:** umbrella #488 · WS0 #489 · WS1 #490 · WS2 #491 · WS3 #492 · WS4 #493 ·
WS5 #494 · WS6 #495 · O4 prerequisites #499 (anchors: #209 → O3; **#43 answered only in
its cosmological half** — its non-cosmological residue, `T̄ = 0` as the linearization's
load-bearing assumption and the unrebutted BGMF self-gravity criticism
(`literature/2510.17094/`), stays open; #360 scope updated; #477 halted)
**Orchestration:** one orchestrator session holds this document current; workstream
sessions are dispatched by the user from the prompts in `docs/cosmology/handoffs/`.

**Companion documents** — this file is the *operational* record (decisions, ladder,
workstreams, gates). The physics lives alongside it:

| Document | Covers |
|---|---|
| `docs/cosmology/primer.md` | How a CMB pipeline computes `C_ℓ`, the two levels (background vs perturbations), what we reuse vs build, why an in-house solver, what a Cobaya `Theory` class is, the `a(η)` symbolic-vs-tabulated resolution. **Read first if new to the program.** |
| `docs/cosmology/spectator_route.md` | The spectator/test-field approximation: validity criteria, reachable vs unreachable observables, the empty-niche argument (H7) |
| `docs/cosmology/birefringence_notes.md` | O4's foundation: the CS/CFJ mechanism, conformal invariance as a solver fast path, frequency scaling, the linear-vs-circular terminology trap, the miscalibration degeneracy, likelihood availability |
| `docs/cosmology/observable_ladder.md` | H2's per-rung feasibility study: capabilities, known-answer targets, costs, invalidators, likelihood availability, ordering |
| `docs/cosmology/magnetic_field_background.md` | O3's assumed primordial magnetic field: models, bounds, and why its spectator status must be computed rather than imported |
| `docs/cosmology/repo_reshape.md` | H4's new-package design: target layout, port inventory, the CAMB/PSALTer convention decision (#513), and why the `derive` gate is not a byte diff |
| `docs/cosmology/spectrum_design.md` | H6's spectrum-module design: two-stage architecture, the Schur-complement kinetic-matrix criterion as primary algorithm, massless-sector layers, Stage-1 export contract, validation oracles |
| `docs/cosmology/stage1_engineering_plan.md` | H8's Stage-1 engineering study: install route, three-tier verification, derivation-branch architecture, exporter design, cost protocol. **Corrects H6 on six points from live sources** (#521–#523) |
| `docs/cosmology/torc_pipeline_audit.md` | H1's audit of the TorC paper, forks and Zenodo archive |

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
- amplitudes `|h| ≪ 1` and `|f| ≪ |F̄|` *(amended 2026-09-06: this read `|h|, |f| ≪ 1`.
  `f` is a field-strength perturbation, so `|f| ≪ 1` is not dimensionally meaningful —
  the criterion in the source, `spectator_route.md` §2, is that it be small **relative
  to the background field** `F̄`)*;
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

## Decisions register

Every decision settled anywhere in the program, in one place. Today "what has been decided?"
otherwise means reading ~6,500 lines across eleven documents — nothing is wrong with any of
them, there was simply no index.

| Ref | What was settled | By / when | Detail |
|---|---|---|---|
| **D1–D9** | Program-level direction — see the table below | user, 2026-08-29 | this document |
| **Rung order** | `O0 → O1 → O2 → O4a → O3 → O4b/V-modes`. Technically neutral, so a pure scientific-priority call. **Gated on #503**: O4a is only cheap if the operator is `n = 0` and `β` is constant over recombination | user, 2026-09-04 (H2 recommended) | §Rung order, `observable_ladder.md` §5 |
| **Integration target** | **Option (iii)** — our own coupled-block solver chained to **unmodified** CAMB. Patching CAMB's Fortran and building on DISCO-EB both rejected; learn from DISCO-EB, never copy its code (#516 closed) | H3, 2026-08-31 | `solver_design.md` §6 |
| **Two solver kinds** | O2 needs an oscillation-resolving mode-equation solver (`~1–10³` oscillations); O3 needs an eikonal amplitude engine with coherence-patch averaging (**`~10²⁹`** — not steppable). **Two engines, one shared core** | H2 → H3, 2026-08-30/31 | `observable_ladder.md` §0.1, `solver_design.md` §7 |
| **Binding cost** | `η`-grid **assembly**, not the matrix exponential (#518). Optimize there first | H3, 2026-08-31 | `solver_design.md` §1 |
| **Stepper choice** | By **bake-off**, not a priori — every candidate implemented and trialed against the §10 protocol | H3, 2026-08-31 | `solver_design.md` §11 |
| **R1 — O1 scope** | O1 is a **fixed-table pass-through plumbing gate**, not a TorC reproduction. Full posterior reproduction (O1b) rejected for now — a nested-sampling campaign for a plumbing check, with the nonlinear TorC background inside the fast block | H1, 2026-08-30 | `torc_pipeline_audit.md` §R1 |
| **R2 — CAMB patch** | **Re-apply cleanly** off upstream `2.0.3`, own fork. Never inherit `slegner/CAMB` (180 commits behind, breaks a public signature); upstream's own `set_w_a_table` is the pole-prone route the paper abandoned | H1, 2026-08-30 | `torc_pipeline_audit.md` §7 |
| **O1 gate** | Constant table `P = −ρ` → ΛCDM `C_ℓ` to machine precision; `TorC_rhopa.py` table → reference `H(a)`/`C_ℓ`; plus the free `set_w_a_table` cross-check at `ϖ_r = 0.8` | H1, 2026-08-30 | #494 |
| **#513 — conventions** | The new package emits **CAMB and PSALTer conventions natively** — no translation layer. Conformal time already agrees with `camb.symbolic`. Gauge is an explicit named input, applied symbolically, carried as spec metadata, asserted at the CAMB seam | H4, 2026-08-31 | `repo_reshape.md` §2.8 |
| **`derive` port gate** | **Not a byte diff** — a consequence of #513. The frozen legacy spec is the oracle *for the physics*; equivalence is semantic, recorded as a written mapping | H4, 2026-08-31 | `repo_reshape.md` §5.2, #525 |
| **Namespace** | Fresh top-level **`tidalcosmo`**, in this repo, renaming to `tidal` after legacy deletion (M7) | user, 2026-08-31 | `repo_reshape.md` §1.3, #490 |
| **Migration shape** | Strangler fig: new code **never imports legacy**; capabilities *ported*, never adapted; two suites; oracle frozen as committed **data** at M0.5 *before* porting, which is what makes early deletion safe | H4, 2026-08-31 | `repo_reshape.md` §7–8, #525 |
| **WS6 algorithm** | **Schur-complement kinetic-matrix criterion** (arXiv:2506.02111) as primary, pole masses/residues as cross-check. Bypasses the radicals that make residue analysis intractable for parity-violating theories | H6, 2026-09-03 | `spectrum_design.md` §5 |
| **WS6 architecture** | **Two stages**: Wolfram/PSALTer once at derivation time exporting a coupling-linear tensor contract; numeric evaluation per sample, no Wolfram in any hot path | H6, 2026-09-03 | `spectrum_design.md` §1 |
| **WS6 scope** | **Minkowski-only is correct and sufficient** — the spectrum screens *in vacuo*; the same split TorC uses | H6, 2026-09-03 | `spectrum_design.md` §2 |
| **`J^P` labels** | Exported **explicitly**. The positional convention is disproven — A23Theory yields J-blocks mixing parities | H8, 2026-09-04 | `pub_wxf.py` warning |
| **Tier-1 vs Tier-2 gate** | Gate the **install** on Tier 1 (author's own input *and* committed `.mx` oracle, so a mismatch can only be the install); Tier 2 mixes install-correctness with our authoring and does not localize. **Requires pinned SHAs** — they default to `HEAD` no longer | H8 + orchestrator, 2026-09-04 | `stage1_engineering_plan.md`, #526 |
| **Term-roster exclusion** | Excluded from the Stage-1 contract | H6/H8, 2026-09-04 | `spectrum_design.md` §6.1 |
| **Config layer** | Build it provisionally, **mark it disposable** — WS1's real surface must *replace* it, not extend it. The failure mode is a provisional layer becoming the de facto config by accretion | orchestrator, 2026-09-04 | #527 |
| **`Method→"Hard"`** | A **dead option** on v2.0.2 — declared, never read. It never gated anything; the cost to measure is plain `ParticleSpectrum` wall time | H8, 2026-09-04 | #521 |
| **Linearity enforcement** | **Ours, not PSALTer's** — `NonLinearCouplings` is never thrown and a bare numeric coefficient passes silently. Our validator is load-bearing for correctness | H8, 2026-09-04 | #522, #527 |
| **WXF** | Something we **produce**, not an interchange we receive: PSALTer writes no WXF and populates two association keys; the exporter reads private `$Local*` globals | H8, 2026-09-04 | #523 |
| **WS6 dependency** | **Not** independent — contradicting an earlier claim in this document; H2's dependency graph governs | H4 correction, 2026-08-31 | `observable_ladder.md` |
| **WKB sequencing** | Matrix-WKB is **built alongside Magnus in the first WS3 handoff**, not gated behind a bake-off. The bake-off decides composition and thresholds on measured numbers; no candidate is discounted on paper estimates | user, 2026-09-05/06 | `solver_design.md` §8/§12, #519 |
| **Admissible theories** | Scoped **explicitly**: only theories admitting the assumed background are in scope, the user judges whether theirs qualifies, and the background-EOM residual is a **first-class per-theory gate** rather than a diagnostic | user, 2026-09-05 | `spectator_route.md` §3, #501, #531 |
| **Oracle re-run rule** | **If `tidal/` or `examples/data/` changes, re-run `scripts/oracles/` in the same commit.** Stated as an action, not a prohibition — "never edit legacy" is unenforceable and would forbid fixing an open bug, whereas this is checkable at review and keeps the frozen oracle honest by construction | I-525 + orchestrator, 2026-09-06 | `repo_reshape.md` §8, `tests_cosmo/data/oracles/README.md` |
| **Guardrails retire too** | `tests_cosmo/test_package_boundary.py` and the shell-out check are **deleted at M7, not adapted** — after `git mv tidalcosmo tidal` their own regexes match every legitimate import, so the guard inverts into a blocker | I-525, 2026-09-06 | #537, §7's M7 row |
| **`A[0]` owner** | Answered by **I-S1B**, which has a real export; I-S1A represents both forms and carries it as a named unknown. Documents cannot settle it | orchestrator, 2026-09-09 | #527, #522 |
| **Oracle CI, with its expiry** | `oracle.yml` runs `freeze_legacy_oracle.py --check`, path-filtered to `tidal/**`, `examples/**`, `scripts/oracles/**`, the fixtures, `uv.lock` and itself. Firing on `uv.lock` — every version bump, ~4 min — is accepted. **Expires with legacy `inspect`/`validate` (M6/M7)**, recorded in the workflow header and §7's M3 row | I-REM + orchestrator, 2026-09-09 | `.github/workflows/oracle.yml`, `scripts/oracles/README.md` |
| **Wolfram gates are lane-held, not CI jobs** | No engine in CI (single license), so the Tier-1 gate and `verify --require-psalter` run on the lane and commit their evidence; CI keeps the evidence readable (`summarize_diff.py`, no lane) and the Python suites green. No separate CI job is added for them | I-REM + orchestrator, 2026-09-09 | `docs/cosmology/evidence/`, `scripts/psalter/README.md` |
| **The certified PSALTer configuration, and what it is independent of** | Wolfram **14.3.0** × xAct **1.3.0** bundle × PSALTer v2.0.2 `bb45adb0` × local registration of `LinearlyIndependent` (the engine-bundled/userbase `ResourceFunctionHelpers` file, sha256 `7bc228a2…`) and `PolynomialDegree` (the committed definition notebook, sha256 `c233e226…`) under **fixed UUIDs** `d40a8dd6-…`, `2f89f2e6-…`. `verify --require-psalter` asserts *provenance* — those UUIDs on master and subkernel, the engine under its mount, the xAct fingerprint — and each assertion was watched fail. `scripts/psalter/ensure_registered.sh` re-creates the registration idempotently (rebuild simulated from an empty registry: identical UUIDs, exit 0). **Independent of the Wolfram-ID cloud login** (asserted logged in and logged out, stored credentials untouched) and of the engine minor version (14.2.1 cross-check MATCH, `evidence/tier1-20260911-engine-142/`, then retired). Nobody's login is logged out: it lives in the user's own home mounts, not the repo | orchestrator, 2026-09-11 | `docs/cosmology/evidence/tier1-20260911-pass/`, `scripts/psalter/README.md`, #559 |
| **`theory_radial` stays excluded, not retired** | Three of M0.5's four excluded theories derive in under a minute (9 s / 45 s / 49 s) and are now pairs; `theory_radial` aborts in `ParseMultiFieldRHS` with `a_0`'s equation arriving as literal `0`, measured twice. #547's own plan said retire it — **overruled**: retiring legacy is milestone work (§7's M5/M6/M7 schedule), so it is recorded as un-derivable with the measurement in the manifest instead. A generic exclusion reason became a *false* one the moment the lane was held, so reasons are now per-theory | orchestrator, 2026-09-12 | #547, `scripts/oracles/freeze_legacy_oracle.py` `MEASURED_EXCLUSIONS` |
| **Exit 0 is not evidence a derivation ran** | An uncaught `Throw` stops a `wolframscript` script and still returns status 0 (probed, `-code` and `-file`), and the Wolfram pipeline signals its own errors by throwing — so `derive` verifies the **artifact**: the output JSON must exist and its mtime must have advanced. The Wolfram-side fix (`Catch` + `Exit[1]`) is deferred to #513's own script emission, because it would change every generated script and so every `derivation_hash` | orchestrator, 2026-09-12 | #561, `tidal/cli/_derive.py` |
| **The engine test is a file, never `command -v wolframscript`** | The dev container image ships a **second**, cloud-only `wolframscript` at `/usr/bin` (#565), so `command -v` succeeds with no engine installed and every later call evaluates in the cloud — an installer reported "already installed" on a bare machine, a verifier passed against the cloud, a registration step burned `timeout 600` before failing. The test is `[ -x "$HOME/.local/wolfram/engine/<series>/Executables/WolframKernel" ]`. Enumerated and fixed at **nine** sites rather than the one that bit | orchestrator + I-ONB, 2026-09-12 | #559, #565, `0078b3c9`, `fdbfdcdc` and `dc664076` |
| **Reduce before you bisect** | A protocol rule, not advice: find the smallest input that still shows the defect and iterate there. #543 reproduces on one scalar field in ~30 s where the CTEG gate takes 7 min | I-543 + orchestrator, 2026-09-09 | delegation protocol, `scripts/psalter/repro_543.wl` |
| **Verification gates** | Made **able to fail** — `tidalcosmo/` had been outside pyright, coverage, `testpaths` and CI, and the never-import-legacy rule had no test | coherence pass, 2026-09-04 | `8b54fe6e`, #524 |

### Still open, deliberately

| | |
|---|---|
| **When to start the `tidalcosmo` version line** | Natural trigger: M0, or O0 first passing. No cost to deferring |
| **GPL/MIT on the Barker-derived code (#495)** | A **release** gate, not an implementation one. Implementation proceeds; distribution does not |
| **Whether non-sampled constants enter the coupling vector** | Decides whether `A[0]` is legitimate (recommended: yes, they stay out, so the linearity check becomes "no *sampled* coupling appears non-linearly"). **Documents cannot settle this — it needs a real export** (#522). Owner: **I-S1B** (the `A[0]` row above); I-S1A carries it as a named unknown |
| **O4a's precondition** | Whether the chosen operator is `n = 0` with `β` constant over recombination. Determines whether the settled rung order's cheap rung exists (#503) |

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
| **O3** | **Gertsenshtein mixing on FRW** — graviton↔photon with Hubble damping | time-dependent solver on the thesis's own physics; **core goal** (thesis couplings) | GH #209; Cembranos arXiv:2302.08186 (local); reduces to flat-space result as `H → 0` **Requires an assumed PMF model** (`magnetic_field_background.md`) — an extra, irreducible assumption no other rung carries; its own spectator status `r_B ≈ 10⁻⁷ B₋₉²` is computed per run, not imported. **Different solver kind from O2** (eikonal + patch averaging, see below). |
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

### Rung order after O1 — SETTLED (user, 2026-09-04)

**`O0 → O1 → O2 → O4a → O3 → O4b/V-modes`.** H2's recommendation
(`observable_ladder.md` §5), adopted. The numbering O1→O2→O3→O4 is **capability order**;
this is the execution order.

The choice was a pure scientific-priority call, because it is **technically neutral**: the
coherence audit confirmed H3's implementation order and H4's migration order both survive
either ordering, since O4a needs no solver at all. O2 is first either way and builds the
spine both routes reuse.

**Why O4a before O3.** O4a (isotropic `β`, post-processing rotation) is the cheapest new
observable in the ladder: the CS/photon sector is conformally invariant, so `a(η)` drops out
and `β` is a dispersion relation integrated against a zero-mode — **no FRW photon solver at
all**, and the observable is an array rotation of CAMB's `C_ℓ` (Tier 3a, against O2/O3's
Tier 3b). It reuses O1's wiring, it is the only rung targeting a **positive 4.8σ signal**
rather than an upper limit, and its derivation work (#503) sits on the thesis's parity-odd
critical path regardless. It also exercises the validity machinery on a cheap rung before
O3 needs it.

> **⚠ Precondition to carry forward.** O4a stays cheap **only if** the chosen operator is
> `n = 0` **and** `β` is effectively constant over recombination. Otherwise it is O4b, which
> is much more expensive, and the ordering advantage disappears. **#503's per-operator
> derivation is what determines whether the cheap route exists** — so #503 is not merely a
> prerequisite of O4a, it is the gate on this ordering.

O3 remains the core scientific goal and is unchanged in importance by being third; it
carries the most unknowns of any rung (a second solver kind, a mandatory PMF model, an
amplitude-based measurement, and a GR signal unobservable by 11+ orders, so its result is a
torsion-*enhancement bound* whose interpretation leans on exactly the background-consistency
and spectator flags O2 and O4a exercise first).

## The niche (H7 investigation, 2026-08-29)

**"Arbitrary Lagrangian → spectator-sector perturbations on fixed ΛCDM → observables →
Cobaya likelihood" is empty**, verified by enumerating near-misses: xPand/xAct (symbolic
only), CppTransport/PyTransport (inflation only, arXiv:1609.00380/1), hi_class (Horndeski
gravity sector, full backreaction, arXiv:1909.01828), SymBoltz.jl (equations in, not
actions, arXiv:2509.24740 — possible future backend). **Nobody has treated torsion as a
spectator on FRW at all**; arXiv:2506.17017 and TorC both signpost cosmological
perturbation theory as future work. Full findings: `docs/cosmology/spectator_route.md`.

## Workstreams

Order: **WS0 → WS1 → WS2 ∥ WS3 → WS4 → WS5.**

**WS6 is NOT independent** (corrected 2026-08-31, flagged by H4). It was listed that way
originally; H2's own dependency graph (`observable_ladder.md`) has
`WS6 -.gates.-> O2` and `WS6 -.gates.-> O3`, and PSALTer's `_extract/likelihood.py` shows
why: on the **sampling path** the vacuum-spectrum screen is what rejects sick coupling
points before they reach the solver. WS6 can be *built* in parallel, but O2 and O3 depend
on it once either is sampled rather than run at a single point.

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

- **WS0 — Research & scoping** (no code): handoffs H1–H5; integration-target decision —
  **SETTLED by H3 (2026-08-31, `docs/cosmology/solver_design.md` §6): option (iii), own
  solver chained to unmodified CAMB, for BOTH engines over one shared core.** O3 forces
  it regardless of the O2 verdict (Boltzmann codes are `k`-resolved but frequency-*integrated*
  — one photon fluid per `k`, not a photon state per `ν`; see `solver_design.md` §2.1's
  amendment — so there is no per-frequency photon
  propagation); DISCO-EB stays the fallback only if gradient sampling ever becomes a
  requirement (GPL decision recorded: learn freely, never copy literally — #516 closed).
- **WS1 — New package, strangler-fig**: the existing framework is **legacy, not a
  template**. New package designed whole, from the goal backwards (H4); **new code never
  imports old code — no adapters**; useful capabilities are *fully ported* (redesigned,
  with docstrings/why-comments/issue references traveling along; original kept as test
  oracle). Two fully separated packages/CLIs during transition; old `tidal` stays for
  cross-checks; legacy retired deliberately once covered.

**⚠ The `derive` port gate is NOT a byte diff** (H4 §5.2 — this supersedes earlier
orchestrator guidance that byte-equality on spec JSON was the strongest gate). That held
only while we assumed legacy notation would be preserved. **#513 decides the opposite**:
the rebuilt symbolic layer emits **CAMB and PSALTer conventions natively** — no conversion
layer — so variable names, formats and gauge conventions change *deliberately*, and
`tidal inspect OLD --diff NEW` would report every intended change as a failure.
The gate instead: the frozen legacy spec stays the oracle **for the physics**, and
equivalence is established by careful comparison — same equations, coefficients and signs
under the new naming — **recorded as a written mapping committed alongside the fixtures**.
Slower and more human than a diff; the honest gate for a deliberate reformatting. `--diff`
stays useful *within* the new package for unintended drift between re-derivations.

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
  **Design complete (H3, 2026-08-31): `docs/cosmology/solver_design.md`** — two engines
  (mode-equation ladder for O2; eikonal amplitude engine + patch averaging for O3) over
  one shared core; the binding cost is coefficient *assembly*, not the exponential;
  every candidate ships as a registry backend and the composition is decided by the
  bake-off (#492, #517–#520).
- **WS4 — Observables**: vector observables; line-of-sight sources; per-rung observables;
  post-processing rotation exact only for constant/isotropic/frequency-independent
  effects, else inside the LOS integral (arXiv:2209.07804); cross-validate vs nanoCMB
  (arXiv:2602.23466), CAMB, class_rot (arXiv:2111.14199).
- **WS5 — Cobaya extension**: Theory class as above; latest upstream CAMB + Cobaya;
  likelihoods per rung (no drop-in birefringence likelihood exists — escalate: Gaussian
  prior on published `β` → fork `LilleJohs/cosmic-birefringence-planck-act` (MIT) →
  SPT-3G BB lite for anisotropic).
- **WS6 — Numerical polology**: **design settled by H6** — `docs/cosmology/spectrum_design.md`
  is the authority (it supersedes the #360 plan as the starting point). Adapt Barker's
  published code (D6) into `tidalcosmo/spectrum/` **and implement the unimplemented
  Schur-complement criterion** (arXiv:2506.02111) as the primary algorithm, residues as
  cross-check; cross-check vs SupplementalMaterials-2607 (JAX kernels *and* the committed
  Julia `*_unitarity_NS` chains); reproduce Lin–Hobson–Lasenby inequalities; Minkowski-only is
  correct and sufficient (spectrum screens in vacuo; TorC used the same split). **H8's three open questions — answered (orchestrator, 2026-09-04):**

1. **Yes, amend the design.** `spectrum_design.md` §4.5/§14.1 name `Method→"Hard"`, and
   `handoffs/H8.md` went further and made its cost *the* go/no-go for the two-stage
   architecture. On v2.0.2 (`bb45adb0`) the option is **declared but never read** (#521), so
   both statements are wrong and one of them is a stated gate. A finding that invalidates an
   instruction belongs where the instruction lives, not only in a measurements record — the
   same rule applied to the `cond(V)` principle. Amend §4.5/§14.1, **pinned to the version**
   (a later PSALTer release could implement it), and the cost to measure is simply
   `ParticleSpectrum` wall time.
2. **Gate on Tier 1; Tier 2 is the physics gate.** H8's recommendation, confirmed, for the
   reason that an *install* gate must fail unambiguously: Tier 1's input and its committed
   `.mx` result are both the author's, so a mismatch can only be the install. Tier 2 mixes
   install-correctness with our own authoring, so a failure does not localize.
3. **Build the provisional config layer, but mark it disposable.** WS1 is next in the
   dependency order yet is *not dispatched*, so "imminent" is not knowable — and blocking
   Stage-1 on it could stall indefinitely. Build the minimal input model, name it so nobody
   mistakes it for the surface, and make it a hard requirement that WS1's real surface
   **replaces** it rather than extending it. The failure mode to avoid is a provisional layer
   quietly becoming the de facto config by accretion.

Stage-1 Wolfram engineering was split out as handoff **H8** — **done as a study**
(`docs/cosmology/stage1_engineering_plan.md`); the engineering itself is not yet dispatched.
Its sole blocker is the PSALTer install (#526), which also gates #521, #522 and #523.

## Parked

`feat/ws2-localized-path-audit` / #477 arc: **already halted** (issues arose making it
not possible; likely unrelated to the new direction). Parked as-is with state recorded on
the issue. Also parked: Phase E (T6/T8 rescue), Phase A-γ, NEXT_PHASES G/H/I, Wolfram CI
(#69). Nothing gates on any of it.

## Handoffs

Prompt files in `docs/cosmology/handoffs/` — self-contained, dispatched by the user to
separate sessions; **this orchestrator session does not launch them**. H7 was executed
during planning; its findings are in `docs/cosmology/spectator_route.md`.

**All eight are complete.** Each prompt file carries a status header giving its date,
artifact and outcome. They are kept, not deleted: a prompt records *what was asked*, which
is what lets a reader judge whether the deliverable answered it. **None is an open
assignment** — H2.md in particular opens "if it does not exist yet, ask before proceeding",
which was live instruction text for a task finished five days earlier.

| ID | Task | Artifact |
|---|---|---|
| H1 | TorC pipeline audit — **DONE** (2026-08-30). Settled R1 (O1 = fixed-table pass-through, above) and R2 (re-apply the CAMB patch, GH #498). Also for WS5: the provider→consumer wiring template (§3.1), the `planck_clik` NaN guard as a flagged rejection path (§3.3), and a caution that stock Cobaya's PolyChord does not give correct evidences under non-uniform priors without the Ormondroyd patch (§3.2). For the spectator flags: `ΔN_eff` (§1.7) is already parameterized for the torsion sector | `docs/cosmology/torc_pipeline_audit.md` ✅ |
| H2 | Observable-ladder feasibility — **DONE** (2026-08-30). Found O3 needs a *different solver kind* from O2 (10²⁹ vs 10³ oscillations → eikonal + patch averaging; H3.md extended accordingly), that O4 splits into three channels with different prerequisites, and that no paper enforces the PMF's own spectator status | `docs/cosmology/observable_ladder.md` + `magnetic_field_background.md` ✅ |
| H3 | Solver design study — **DONE** (2026-08-31). Two engines over one shared core; architecture (iii) settled; assembly (not expm) identified as the binding cost; six-stepper prototype registry with measured orders; matrix-WKB designed + prototyped (bake-off-gated, #519); 18-paper reference set acquired. Follow-ups #517–#520 | `docs/cosmology/solver_design.md` ✅ |
| H4 | New-package design — **DONE** (2026-08-31). Target layout + `tidalcosmo/` scaffold (READMEs only, no Python); port inventory; **#513** adopt CAMB/PSALTer conventions natively (gauge as spec metadata), which is why the `derive` gate cannot be a byte diff; #514–#516 | `docs/cosmology/repo_reshape.md` ✅ |
| H5 | Literature acquisition — **DONE** (2026-08-30). 20/20 fetched and title-verified; `literature/README.md` now tracked and auto-generated by `scripts/bibaudit/index_literature.py` (GH #497) | populated `literature/` + `docs/references.md` ✅ |
| H8 | Stage-1 engineering — **DONE as a STUDY** (2026-09-04; user pivoted it to "research + recommend, report to the orchestrator"). Install route, three-tier verification, derivation-branch architecture, exporter design, cost protocol, plus `scripts/research/psalter_stage1/` (rescues H6's WXF decoders from /tmp; pinned upstream fetch script — the route, not the GPL payload). **Six live-source findings correct H6** — chiefly that `Method→"Hard"` is a **dead option** (#521), so the "go/no-go cost" I wrote into the prompt was based on a false premise and is simply `ParticleSpectrum` wall time | `docs/cosmology/stage1_engineering_plan.md` ✅ |
| H6 | Numerical polology design — **DONE** (2026-09-03). Two-stage architecture (symbolic PSALTer once per structure → JAX per sample, order-ms); **primary algorithm = the Schur-complement kinetic-matrix criterion of arXiv:2506.02111** (published, worked, implemented nowhere — we implement it; it handles massive–massless kinetic mixing and carries parity-odd support), residue route kept as permanent cross-check; massless sector handled by counting + compiled Stage-1 closed forms + numerical fallback (the released code's `\|Re m\|>0.01` guard would reject GR+Maxwell); Stage-1 export contract defined (explicit `J^P` labels — the WXF `J`-blocks are unlabeled and mix parities, verified on `A23Theory`); seven-group validation set incl. the LHL gate and published EC/TorC oracles; reconciled with H4 (§3.2 — confirms #513 with a hard reason: the residue's parity factor is *defined by* the signature). License → release blocker on #495 | `docs/cosmology/spectrum_design.md` ✅ |
| H7 | Spectator-route scope — **DONE** | `docs/cosmology/spectator_route.md` |

Dispatch order: **H1–H8 all ✅ done.** Stage-1 *engineering* is not yet dispatched — H8 studied it and produced the route; see `stage1_engineering_plan.md`. (H3 absorbed the
post-dispatch scope extension — both solver kinds designed — plus the 2026-08-31 user
directions: Wolfram-at-derivation-only symbolic policy, experimentation-over-reading
bake-off, uniform+stochastic B̄ modes, soft classifier.)

## Implementation delegation protocol

How implementation work leaves this program document and comes back. **The wave board
below is the single source of truth for where things stand** — a session with no other
context should be able to read it and the queue, and know what to do next.

### The loop

1. **Orchestrator writes prompts** into `docs/cosmology/handoffs/` as an `I-<issue>.md`
   file (the `I-` series; `H-` was the research series, all complete). Every "the orchestrator
   will…" commitment inside a prompt is copied onto the wave-boundary checklist as a named item
   **in the same commit as the prompt** — at dispatch, not at the boundary, which is where two
   of Wave 0's were lost.
2. **The user dispatches** each prompt to a separate session. *This orchestrator session
   never launches them.*
3. A delegate works in **its own git worktree** off `feat/cosmology-program`
   (`git worktree add /tmp/tidal-i<issue> -b cosmo/i<issue>-<slug> feat/cosmology-program`),
   **never merges**, and reports its branch back.
4. **The user relays a short update**; the orchestrator runs the merge checklist, merges,
   and updates the board **in the same commit**.
5. At a wave boundary the orchestrator updates memory, and the next wave is **planned
   fresh** rather than executed from a pre-baked plan — later waves are shaped by what
   earlier ones found.

### Prompt template

Header (issue · milestone · wave · **Wolfram lane y/n** · dependencies · **owned paths**)
→ context and an *ordered reading list* naming sections, not whole documents → the
**interface contract** (what it consumes, what it produces, with exact references) →
**quantitative success criteria stated before any code**, plus the tests to add under
`tests_cosmo/` → a **scope fence** naming the adjacent temptations explicitly → working
rules → the **flaw protocol** → the report-back format.

**Working rules every prompt carries:**

- Worktree off `feat/cosmology-program`; branch `cosmo/i<issue>-<slug>`; **never merge**;
  never touch the shared working directory (it is the orchestrator's alone).
- **Never version-bump, tag, or edit the changelog.** This overrides `CLAUDE.md`'s
  default-to-bumping rule: parallel delegates share `.git` refs, so concurrent bumps
  guarantee a `pyproject.toml` conflict and a tag collision. The orchestrator bumps once
  per wave.
- **Open a draft PR into `feat/cosmology-program` as soon as you have a commit**, and
  report its number. This is not bookkeeping: `test.yml` fires on `pull_request` and on
  `push` to the integration branch, so **a branch with no PR is never seen by CI**, and the
  orchestrator's local checks provably cannot substitute (see the cspell note in the merge
  checklist). I-524 opened one and its first run caught a defect that had already been
  merged blind; I-525 and I-526 opened none, faced the gate only after merging, and reddened
  the trunk. A draft PR moves that discovery to before the merge, where it is free.
- **Stay inside your owned paths** (listed in the header). Two sessions writing the same
  directory is the one collision a worktree does not prevent.
- **Assertions get verified before they land; hypotheses do not have to be.** The
  distinction is not whether you checked before speaking — it is whether the unchecked
  thing was framed as provisional or as settled. A hypothesis with its test named beside it
  costs one command when wrong. An unverified claim written into an **instruction site**, or
  **acted on**, costs whoever follows it. Two entry points, one rule:
  1. **Writing it into a document another session is told to follow** makes it an assertion.
     Verify first. (Wave 0: an enumeration glob shipped in a prompt produced 45 pairs where
     the same paragraph claimed 46; a boundary test shipped into the repo inverts into a
     blocker at M7 — #537.)
  2. **Acting on it** makes it an assertion too. (Wave 0: a confirmation run was killed on
     the unchecked belief that a delegate was duplicating work; it was finishing an
     hour-old comparison.)

  Corollary already stated in `feedback_gate_not_real_until_it_fails` and worth repeating
  because it was missed anyway: **a rule stated in prose is not enforcement.** If a document
  says "skipped by default" or "enforced", find the mechanism — `--collect-only`, the test,
  the config — before reasoning from it. Three places asserted a publication lane that
  nothing implements (#540).
- **A verification command must not write into what it verifies.** Distinct from a vacuous
  gate, and it fails more quietly: the gate works correctly, it just writes its answer over
  the question. Found the hard way — the committed Tier-1 evidence shipped with a README
  telling the reader to re-verify with `--diff-only`, which *recomputes* the verdict and
  **rewrites** `tier1_diff.json`. With the `.mx` deliberately excluded from the commit, the
  documented command could never have succeeded, and it failed **destructively**: `mismatch`
  became `ours_unreadable`, still well-formed JSON, so nothing announced it. Anyone checking
  the recorded result before the supervisor meeting would have silently replaced it with an
  artifact reading "the gate could not run".
  **So:** the documented way to read back a recorded result is a **read-only** path
  (`summarize_diff.py`), and any recompute path **refuses** when its inputs are absent
  rather than producing a failure artifact in their place. Probe both — refusal *and*
  non-destruction — the same way a gate is probed failing.
- **Any guard you add carries its expiry.** When you write a test, assertion, lint rule or
  CI gate, ask in the same breath *what future change makes this wrong rather than merely
  unnecessary?* If there is one, record it **both** where the guard lives and in the row of
  §7's milestone table that will trigger it, and mark the guard's site with an
  `EXPIRES-WITH: #<issue>` comment so `git grep EXPIRES-WITH` finds every guard whose trigger
  has closed — that grep is a closing precondition for the issue and an item
  `scripts/wave_boundary_check.sh` prints. A guard's correctness has a lifetime, and the
  moment of writing is the only moment its expiry is obvious — the merge checklist cannot
  catch it, because such a failure bites at a milestone nothing is tested against yet.
  (Found the hard way: `test_package_boundary.py` inverts into a blocker at M7, when
  `git mv tidalcosmo tidal` makes its own regex match every legitimate import — #537.)
- **Reduce before you bisect.** Find the smallest input that still shows the defect, iterate
  there, and run the expensive case only to confirm. Wave 0's Tier-1 failure reproduces on a
  single scalar field in ~30 s where the CTEG gate takes 7 minutes — a 14× cheaper loop for
  the same knowledge, and the difference between an afternoon and a week when a hypothesis
  needs a dozen tests. This extends `feedback_minimal_tests` from choosing test theories to
  bisecting bugs.
- **An argument for why something cannot be checked is not evidence about how the check comes
  out.** The two travel together and must be separated. Answering I-ONB's xPerm question the
  orchestrator wrote both: *"this cannot be verified without running it"* — correct, and worth
  saying — and *"the wrapper solves a runtime loader problem that switching to `mcc` would not
  remove"* — a mechanism claim, offered with the same confidence, for a code path that **had
  never executed once**. The merge rehearsal ran it: `mcc`'s output needs no wrapper. Being
  unable to settle it was real and bought nothing; the reasoning smuggled in beside it was the part
  that misled. When you say a thing cannot be settled here, stop there — do not also say how it
  would come out.
- **A caveat that is a scope limit hides better than a doubt.** "I can't test this in my
  environment" feels discharged the moment you name who can, so it reaches the report and the PR
  body and stops there — while the *artifact* states the untested thing flatly. That is how
  `build-xperm.sh` became step 4b of the onboarding path instead of a fallback (I-ONB's own
  diagnosis, 2026-09-12). The test is the same as for any caveat: it must survive into the thing
  the next person follows, not into the thing you hand back.
- **Read the tool's own known-issues before hypothesizing.** PSALTer's README lists "gauge
  symmetries not identified" as known bug #1 — the symptom class #543 spent a wave attributing
  to the engine. I-REM found it by reading; nobody had.
- **A forbidden resource is only forbidden once everything in scope that acquires it is
  enumerated.** "No kernel" means: `wolframscript` directly, `tidal derive`, every
  `examples/*/run.sh`, the gate, the verify script. `wolfram-guard.sh` failed open on the
  wrappers until #555 enumerated them; a prompt that says "no lane" names the list.
- One `wolframscript` at a time, machine-wide — only a prompt flagged for the Wolfram lane
  may start a kernel.
- Conventional commits; no attribution trailers.

**Flaw protocol** — what to do on finding the design wrong:

- A **design-document contradiction or error**: amend it *at the instruction site*, pinned
  to the version, and report it. Do not work around it silently and do not only note it in
  the report — the next reader of that instruction must see the correction.
- An **architectural contradiction** (the design cannot be built as specified): **stop and
  report.** Do not redesign in a delegate session.

**Report back:** branch name · gates run, with their output · amendments made · anything
discovered that the orchestrator should route · suggested next step.

### Orchestrator merge checklist

Run for every reported branch — this is what keeps coherence continuous rather than
requiring another end-of-phase sweep:

1. Read the **full diff**, not a summary.
2. Run every gate: `ruff check`, `ruff format --check`, `pyright`, `cspell`, **both**
   suites, the boundary test.

   **Run `cspell` over the merge's full changed-file set, including your own edits to the
   board and prompt headers** — not over an ad-hoc selection:

   ```sh
   git diff --name-only <merge>~1 <merge> | grep -vE '^(uv\.lock|tests_cosmo/data/)' \
     | xargs npx cspell lint --no-progress
   ```

   *(Learned immediately: the `654b627a` merge turned the trunk red on `uncontended`, a
   word the orchestrator wrote in a status header during that very merge. Local cspell
   flags it — it simply was not run on the file that had just been edited. The checklist
   said "run cspell"; what it needed to say is "over what you changed".)*

   > **⚠ Local `cspell` is a pre-check, not the gate (found by I-524, 2026-09-06).**
   > `npx cspell lint` and the `cspell-action` **disagree**: the action flagged
   > `PYTHONIOENCODING` and `PYTHONNOUSERSITE` in `scripts/oracles/freeze_legacy_oracle.py`
   > where local reported zero — most likely `allowCompoundWords: true` splitting them
   > locally but not under the action's bundled version. Reproduced from both sides:
   > local-with-compounds reports **0**, the action **2**, and local with compounds
   > *disabled* reports **11**, so it over-reports in that mode.
   >
   > **Settled 2026-09-07: no local configuration reproduces the gate, and it is not the
   > compound setting.** A second divergence (`distro`, in the I-526 merge) is absent from
   > local output *even with compounds disabled* — local cspell 10.2.2 carries the word in a
   > bundled dictionary that the action's version does not. The disagreement is in
   > **dictionary contents**, so no flag closes it.
   >
   > **Therefore: CI is the gate, and the workflow is fix-forward.** Push, watch the trunk
   > run, and correct in a follow-up commit. **Run the changed-file check anyway** — it
   > catches typos and flagged words — but record it as a typo pass, never as evidence the
   > gate will hold.
   >
   > **Corrected 2026-09-07 — a premise of this rule was false.** It said "delegate branches
   > get their own CI before merge, so the only uncovered surface is orchestrator edits."
   > They do not: `test.yml` fires on `pull_request` and on `push` to this branch, so **a
   > delegate branch with no PR gets no run at all.** Measured: I-524 had 2 runs (it opened
   > PR #539); **I-525 and I-526 had zero** and were merged having never faced the gate. Both
   > then reddened the trunk. Hence the new working rule below: **every delegate opens a
   > draft PR into `feat/cosmology-program`.**
   >
   > That file reached trunk via `e310e125`, which no CI had ever observed. The new lane
   > caught it on its first execution — the criterion demonstrated rather than asserted.

   **Cap BLAS threads when running the suite** —
   `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.
   Several sessions share one 8-core machine, and an unrestricted run takes ~3 cores of
   BLAS on its own; that is enough to make the wall-clock assertions in
   `test_probe_perf` and `test_inference_eval_perf` fail spuriously (observed: 400 ms
   against an 80 ms budget, where the probe's known quiet cost is 13–27 ms). A timing
   failure under load is not a regression — check the load before reading it as one.
3. Verify the prompt's success criteria **positively, from artifacts** — never from the
   session's report alone.
4. Propagate any instruction-site amendment the session should have made and did not.
5. Merge; close issues with commit references; **update the wave board and the queue in
   the merge commit**; put a status header on the prompt file.

   **Before closing an issue that instruction sites are conditioned on, run the expiry
   grep** — `grep -rn 'EXPIRES-WITH: #<n>' docs/ scripts/` — and re-date or remove every
   hit. Amendments written while a blocker is open are correct only while it is open;
   #543 and #551 each carry an "On close, revisit" section naming this. A checklist, not
   a mechanism: the trigger is a human act, and a test that goes red *because someone
   fixed a bug* would teach everyone to ignore red.

   > **Never state a CI verdict without its run id and fetched conclusion.** Not another
   > reminder to check — a *format* requirement, because five rules did not fix this and a
   > sixth would not either. Write `CI <run-id>: <conclusion>`, which cannot be written
   > without having fetched it. Wave 0 produced **six** instances of reporting an outcome
   > read from expectation rather than from the artifact, every one caught by a delegate or
   > by CI rather than by the orchestrator — including "trunk CI green on `3ff90c2a`" when
   > that run was `completed/failure`, and a `git push --delete` of four branches that
   > exited 0 having deleted one. **The exit code is not the outcome; the run is not the
   > verdict.** Make the claim impossible to write without the evidence.
   >
   > **The rule validated itself within one commit.** A `gh run watch --exit-status` wrapper
   > reported "CI FAILED" for `3c1209fe`; the fetched conclusion was **`cancelled`** — the
   > concurrency group superseding an in-flight run because a later commit was pushed, which
   > is correct behavior and not a failure at all. Wrapping a tool's exit status in a binary
   > label discards exactly the distinction that matters. `gh run watch` for the wait,
   > `gh run view --json conclusion` for the verdict; never the former's exit code as the
   > latter.
6. At a wave boundary: memory + backup.

### The Wolfram lane

One `wolframscript` machine-wide, so **at most one lane-flagged session is ever active** —
the orchestrator included. The next occupant starts only once the previous lane work is
**merged and its gate independently re-run**: once a long run begins, the lane offers no
verification window. Runs that can span a day (the Stage-1 cost measurement, 24 h ceiling)
are **launch-and-collect** — a detached script, collected by a later session.

### Wave board

Status: `drafted → dispatched → reported → merged`.

| Wave | Prompt | Issue | Lane | Owned paths | Branch | Status |
| --- | --- | --- | --- | --- | --- | --- |
| 0 | I-524 — packaging, extras, CI lane | #524 | — | `pyproject.toml`, `.github/`, `tidalcosmo/__init__.py`, `tidalcosmo/cli/`, `cspell.json` | `cosmo/i524-packaging` | **merged** ✅ CI 34051208887 green on the merged SHA: 2885 passed, 39 skipped |
| 0 | I-525 — freeze the legacy oracle | #525 | — | `scripts/oracles/`, `tests_cosmo/data/` | `cosmo/i525-oracles` | **merged** ✅ gate re-run independently: 185 fixtures current, regeneration byte-identical |
| 0 | I-526 — install PSALTer, Tier-1 gate | #526 | **yes** | `scripts/{install-psalter.sh,verify-wolfram-setup.sh,psalter/}`, `tests_cosmo/fixtures/`, `.gitattributes` | `cosmo/i526-psalter` + `fix/i526-evidence-reverify` (#544) | **merged** ✅ install works, all three probes answered. Tier-1 reported MISMATCH on merge — **resolved 2026-09-11 by I-543** (#543: two undocumented Function Repository dependencies, not the engine), and the install is now **certified** at 14.3.0 × `bb45adb0` × local registration. Both evidence states committed: `evidence/tier1-20260907` (fail) and `evidence/tier1-20260911-pass` (`VERDICT: MATCH`) |
| 0c | I-REM — instruction sites, docs index, tooling, oracle CI | #545 #546 #540 | — | design docs, `docs/README.md`, `handoffs/README.md`, `tidalcosmo/**/README.md`, `scripts/`, skills, `Makefile`, `.github/`, config | `cosmo/irem-amendments` (#550) | **merged** ✅ CI 34393863566 + 34393863645 success on `1cde083d`. **`oracle.yml` proven in both directions**: CI 34393316535 failure on a corrupted fixture, CI 34392693815 success clean. Replaced the six-item exporter list with a rule + anchor + guard rather than a longer list. Found PSALTer's own README known-bug #1 (#543) |
| 0c | I-533 — retire the M0 drop rows | #533 | — | `tidal/`, `tests/`, `examples/**/run.sh`, legacy `scripts/`, `docs/tex/` | `cosmo/i533-retire-drop-rows` (#552) | **merged** ✅ CI 34401167660 + 34401167568 success on `f8d4001a`. `sweep`/`sample`/`analyze`/`plot` gone with their in-package plotting, +1182/−26727; each name exits 2 naming a `git show v0.53.0:` recovery path, handler probed. **`measure` deliberately kept** — a `drop` verdict is not a retire milestone; it is drop **+ M5**, amended at the site. Orphan inventory #553, spec drift #554 |
| 0c | I-543 — resolve the Tier-1 gate | #543 #542 #556 #551 | **yes** | `scripts/psalter/`, `docs/cosmology/psalter_543_*.md`, `evidence/tier1-20260911-pass/` | `cosmo/i543-psalter-gate` (#557) | **merged** ✅ CI 34619683909 success on `4a2c38fd`. **GATE PASSES — `VERDICT: MATCH`**, both keys identical. Cause: two undocumented PSALTer Function Repository dependencies, **not** the engine (14.2.1 behaves identically). Certified: 14.3.0 × `bb45adb0` × local registration. Upstream issue drafted, not filed |
| 0d | I-ONB — one onboarding path, fresh-host container, registry volume | #559 | — | `.devcontainer/**`, root `README.md` setup, `scripts/README.md`, `scripts/{install-wolfram-engine,install-xact-xcoba,activate-wolfram}.sh` | `cosmo/onb-one-path` (#564) | **merged** ✅ CI 34712322686 success on `1158235f`. +911/−1784 over 21 files: one six-step path in `WOLFRAM_GUIDE.md`, `initializeCommand` + the `wolfram-objects` volume, the 25-step chain's Wolfram wiring extracted to `setup-wolfram-links.sh` (always exits 0), two installers deleted and two redirected, `tests/test_devcontainer_creation.py`. Found a **second** fresh-host abort my audit missed (`install-lsp-wl.sh:17`), which had been silently skipping the Claude memory restore and session reindex |
| 1 | I-532 — CAMB seam, background protocol, flag schema | #532 | — | `tidalcosmo/{background,spectator,validity}/` | — | planned |
| 1 | I-503 — per-operator dispersion + zero-mode scope | #503 | — | `research/lagrangian_enumeration/`, `docs/` | — | planned |
| 1 | I-S1A — Stage-1 Python side | #527 | — | `tidalcosmo/{config,derive}/` (Python only), `tidalcosmo/spectrum/` | — | planned |
| 2 | I-S1B — Stage-1 Wolfram side + cost run | #495 | **yes** | `tidalcosmo/derive/wolfram/` | — | outline |
| 2 | M1b — Cobaya Theory + ΛCDM posterior | — | — | `tidalcosmo/{spectator,presets,likelihoods}/` | — | outline |
| 2 | M2/O1 — CAMB fork re-apply | #498 | — | fork repo + `tidalcosmo/background/` | — | outline |

### Decision on #543 — the Tier-1 mismatch (orchestrator, 2026-09-07)

> **RESOLVED 2026-09-11 — the gate passes and the install is certified** (#557, `7286ae94`).
> The cause was **two undocumented PSALTer Function Repository dependencies**, not the engine
> version this section spends its length on: `LinearlyIndependent` (`SymbolicNullSpace.m:31`
> master, `IsNullVectorOfSpace.m:6` subkernels) and `PolynomialDegree`. Unresolved, they are
> not Booleans, so no gauge symmetry is identified and every pseudo-determinant is zero — and
> the run still writes its `.mx`. 14.2.1 behaves identically to 14.3.0, so the engine
> hypothesis is refuted by execution. **Certified configuration: Wolfram 14.3.0 × PSALTer
> `bb45adb0` × local registration** (`scripts/psalter/register_resources.wl`, no edit to
> PSALTer). Evidence: `docs/cosmology/evidence/tier1-20260911-pass/` (`VERDICT: MATCH`), with
> the failing 2026-09-07 run kept beside it.
>
> **The reasoning below is left standing because it was wrong in an instructive way.** It
> reads as a careful elimination, and its two named hypotheses were both wrong: the engine
> version, and "the resources are ruled out by test" — the latter because the substitute that
> tested it reproduced the disabled behavior instead of restoring it (#556). An elimination is
> only as good as whether its controls discriminate.


**The gate stays as written and the install stays uncertified.** A gate whose whole
justification is *"a mismatch can only be the install"* is not relaxed the first time it
reports one.

**What the result certifies.** `WaveOperator` is **bit-exact** — 3 sectors, 303 leaves — so
field declaration, decomposition and wave-operator construction are demonstrably correct.
The `1/0` is inside `ConstructSaturatedPropagator`, which feeds the residue route.

> **⚠ Correction (I-526, 2026-09-07) — my narrowing was wrong, and the blast radius is not
> confined to the cross-check.** I argued here that because the Schur-complement criterion
> "does not involve any inversion of the wave operator nor the computation of residues of
> the propagator" (§5), the *primary* path was untouched and only the residue cross-check
> was degraded. **I verified that one link and generalized past it.** The Schur route has a
> **second input** I did not check: §5's own caveat that `O_LL` must be invertible, so gauge
> modes are removed first — and those come from `ConstructSourceConstraints`.
>
> Measured, by replaying the CTEG Lagrangian without its trailing `Quit[]`:
> **`$LocalSourceConstraints` is `{}` — zero rows**, where CTEG's 21-generator formulation
> requires 21 (§3 measures the generator count as exactly that row count). So gauge-mode
> removal is degraded too, and **#543 must be resolved before *either* criterion can be
> trusted on this install.** Evidence strength stated honestly by the session: there is no
> oracle for that key, so this is measured against a *published count* rather than a
> committed artifact — weaker than the Tier-1 diff, but 0-versus-21 is not a near miss.
>
> **A limit on Tier 1 itself, worth stating plainly:** the oracle holds only two keys, so
> even a *passing* Tier 1 would certify the wave operator and the pseudo-determinants and
> **nothing else** — not source constraints, spectrum, or unitarity conditions. That is why
> §3 calls Tier 2 the physics gate, and it means Tier 1 was never going to certify as much
> as its name suggests.

**Next step is the supervisor, not an engine install.** The leading hypothesis is Wolfram
14.2 (the oracle's `.mx` header) versus 14.3 (ours), and testing it needs a 14.2 engine.
Before spending that, ask **Barker** — it is his code, he generated the oracle under 14.2,
and he can likely settle it in a sentence. Bundled with the two undocumented Function
Repository dependencies that fail silently (#542) into one D6 conversation on **11
September**. Only if that is inconclusive is a 14.2 engine worth installing.

**Wave 1 is not blocked.** The Stage-1 Python handoff (I-S1A) tests its reader against the
**committed upstream `.wxf` fixtures**, which exist precisely so it "can be tested without
Wolfram installed" (`PROVENANCE.md`). #543 blocks **I-S1B** — the Wolfram-lane half, Wave 2
— which is where our *own* exports must match. Verified against the design rather than
assumed.

**Deferred, same class — legacy on a deletion path, record the decision rather than fix the
code:** #533 (retire the drop rows), #537 (delete the boundary test at M7), #540 (the
publication-lane markers document an exclusion nothing implements).

### Wave-0 boundary — completed 2026-09-07

Every row `merged`; full gate set green on the trunk (`d3330d23`, CI 34127903873);
delegate discoveries routed (#535, #536, #537, #538, #540, #541, #542, #543; probes #521,
#522, #523 closed); memory updated and backed up. **Wave 1 is planned fresh from here** —
the point of the boundary is that Wave 0's findings get to change it.

**Inputs Wave 1 has that the original outline did not:**

- ~~**#543 blocks both spectrum criteria** … Awaiting Barker, 11 September.~~ **Superseded 2026-09-11:** resolved by I-543 (#557); the current inherited-inputs list is in "What to implement next" below and governs.
- **The Stage-1 exporter must read all eight private globals**, not the six §5 specified —
  both of that section's justifications were refuted on the live install (#523).
- **§6's checkpoint machinery need not be built**: PSALTer already emits a per-function
  trace in CLI mode, so the gate only timestamps stdout.
- **A Tier-1 *pass* would certify two keys and nothing else** — not source constraints,
  spectrum, or unitarity conditions. Tier 2 is the physics gate, as §3 says.

**Wave 3+, named only:** WS2 FRW derive (lane), then the O2 fan-out — WS2 background +
residual ∥ WS3 solver ∥ WS4 line-of-sight ∥ WS6 Stage-2. WS3's first handoff carries
**Magnus and matrix-WKB together** (`solver_design.md` §12), written from the ground up for
a time-dependent background; nothing is ported from `modal.py` and it is not an oracle.

### Wave-0 completion + hardening boundary — completed 2026-09-12

The completion wave's three prompts merged 2026-09-09/11 (I-REM, I-533, I-543) and the
Tier-1 gate passes. **Foundation hardening then ran on top of it**, because four audits found
things that would have made Wave 1 wrong: the certified configuration was not robust to a
container rebuild *and failed green*; a fresh user could not reach a working container; the
"engine plays no role" claim had no artifact; the Wave-1 issue bodies were not self-sufficient;
and A7/#547 had never run.

**Closed here:** #547 (three of the four excluded oracle theories derived — 9 s, 45 s, 49 s;
`theory_radial` measured un-derivable), #554 (drift classes documented, `--staleness` detector
shipped, `run.sh` no longer derives in place), #543/#551's markers discharged. **Filed here:**
#559 (onboarding — the I-ONB prompt implements it), #560 (the #394 volume-element check rejects
every FRW measure), #561 (`derive` reported success for an aborted derivation — **fixed**),
#562 (a shipped spec now carries a Wolfram call form, invalidating #428's premise).

**Two defects came from execution, not review**, and both are the same shape: a check that
could not fail. `verify --require-psalter` asserted resource *behavior*, which a resource
fetched from the repository would also satisfy — so a rebuild could have certified an
uncertified leg. And `wolframscript` exits **0** after an uncaught `Throw`, so every
pipeline-detected derivation error looked like success. Each is now an assertion that was
watched failing.

**The boundary record** — `scripts/wave_boundary_check.sh --with-wolfram`, which exists
because prose had failed this orchestrator three times:

```
    == wave boundary check @ 7bf15d77 (2026-09-12T15:47:37Z) ==
      GREEN  working tree clean
      GREEN  nothing unpushed
      GREEN  no delegate worktrees
      GREEN  no remote cosmo/* branches
      GREEN  CI 34702806495 [Run Tests] completed/success
      GREEN  every open cosmology issue has a milestone
      GREEN  no EXPIRES-WITH marker for a closed issue left un-dated
      GREEN  every merged prompt carries a STATUS header
      GREEN  MEMORY.md 98 lines / 16992 B
      GREEN  memory index complete (every link resolves, no orphan)
      GREEN  tidal sweep exits 2 (retired, names its replacement)
      GREEN  tidal sample exits 2 (retired, names its replacement)
      GREEN  tidal analyze exits 2 (retired, names its replacement)
      GREEN  tidal plot exits 2 (retired, names its replacement)
      GREEN  oracle --check: 197 fixture(s) are current
      GREEN  oracle --staleness: 3 current, 40 stale, 6 without a hash, 0 errors (of 49). 'stale' = tidal derive 
      GREEN  version 0.54.1 is tagged (v0.54.1)
      GREEN  verify-wolfram-setup --require-psalter exit 0
    == ALL GREEN ==
```

`CI 34702806495: success` [Run Tests] on `7bf15d77`, and `CI 34701324911: success`
[Frozen Legacy Oracle] on `e9796be1`, which carried every oracle change. Version **v0.54.1**
tagged and pushed. The Tier-1 gate was re-run from scratch from an **empty** registry:
`VERDICT: MATCH`, `tier1_diff.json` identical to I-543's but for `generated_utc`
(`evidence/tier1-20260911-recert/`, with the run's own `.mx` committed).

**Kernel rehearsal at merge (2026-09-12), the orchestrator's half of I-ONB's criteria.**
Steps 4–6 run against an **empty** userbase, with `~/.WolframEngine` pointed at a scratch
tree — `$UserBaseDirectory` itself cannot be redirected (`-userbase` is accepted and then
refused: the symbol is `Protected`), and the symlink is the indirection the engine actually
follows, which is what a fresh user's points at. The engine stayed licensed throughout, via
the container-wide `/usr/share/WolframEngine/Licensing/mathpass` link. Result: xAct **25 s**,
PSALTer **26 s**, `verify --require-psalter` **51 s** → *All checks passed* on a userbase that
began with no xAct at all. Certified tree confirmed untouched afterwards (registry byte-tree
identical to its pre-rehearsal snapshot; PSALTer and xPerm timestamps unchanged; verify green).

It also **corrected an instruction in the merge it was verifying**: `install-xact-xcoba.sh`
alone produces a working xPerm — it recompiles with the engine's own `mcc` (the resolution
I-ONB repaired, which had never run before) and the result connects with no `LD_LIBRARY_PATH`
wrapper. `build-xperm.sh` is a fallback, not a step (`c75d6cc6`). That reverses the orchestrator's
own answer to I-ONB's xPerm question, which had said the wrapper's need would survive the switch
"unless mcc sets an RPATH, which cannot be verified without running it". It has now been run.

**What is open at this boundary: nothing in the completion wave.** I-ONB merged 2026-09-12
(#564, `CI 34712322686: success`), closing #559. Carried forward with owners, none of them a
Wave-1 dependency: #558 (β over recombination, with #503), #548 (two specs without a TOML, M3),
#560 (the FRW volume-element check, WS2), #562 (Wolfram call forms in `sign_algebra`, M3),
#563 (the undeclared `xPerm.m` patch, M0.5) and #565 (the image's second `wolframscript`,
M-parallel).

### Wave-boundary checklist

What the next planning session does *first*, before planning anything. **Run
`scripts/wave_boundary_check.sh` before reading on** — it prints the mechanical half of this
list (tree and unpushed commits, worktrees, remote `cosmo/*` branches, CI conclusion *fetched*
for `HEAD`, open cosmology issues without a milestone, `EXPIRES-WITH` markers whose issue is
closed, prompt files without a STATUS header, memory size and index, retired names exiting 2,
oracle `--check` and `--staleness`, version vs tag; `--with-wolfram` adds
`verify --require-psalter`). Three orchestrator omissions in the completion wave were all
items it prints red; prose had failed three times.

1. Confirm every row of the wave reads `merged`. **A row whose deliverable is a gate reads
   `merged` only when the gate *passed*** — I-526 was recorded `merged` with its gate
   reporting MISMATCH, which is how the wave was declared closed around an unresolved
   blocker.
2. Re-run the full gate set once on the merged trunk — **including any lane-held gate, from
   scratch, not read back.** A delegate's recorded verdict is not verification; only a second
   run, by another session, from a clean run directory, makes it independent of that session's
   procedure and environment. (Wave 0: the orchestrator substituted `--diff-only` — which only
   re-compares the delegate's own artifacts, and destroyed them, #544 — then a read-only
   summary. Attempting the real thing on 2026-09-09 found that the gate could not be run at
   all, from 66 seconds after its only run, #549.)
3. Read the delegates' reported discoveries and **route** each — amend at the site, open an
   issue, or fold into the next wave's scope.
4. **Prune stale worktrees.** `git worktree list`; for each entry whose branch is an
   ancestor of `HEAD` and whose tree is clean, `git worktree remove` it and delete the
   branch. **Never touch a worktree whose session is still active** — check for unmerged
   commits and uncommitted changes first. A stale registration holding a branch ref
   quietly confuses a later `git worktree list`.
5. **Bump the version once, for the wave, and tag it.** The rule lived only as a clause
   inside the *delegate* working rules ("the orchestrator bumps once per wave"), and neither
   checklist named it — so Wave 0 shipped three feature merges under `v0.53.0`, which had been
   cut *before* them and therefore looked current.
6. **Copy every "the orchestrator will…" commitment out of the wave's prompts and onto this
   checklist as a named item.** Both misses above share one cause: a commitment addressed to
   the orchestrator was recorded somewhere other than the checklist the orchestrator executes
   — I-526's criterion 2 said "the orchestrator will re-run this before Wave 1 is dispatched",
   and nothing carried it forward.
7. Update memory (`project_cosmology_program.md`, the MEMORY.md status line) and back up.
8. *Then* plan the next wave, in detail; the one after it in outline; nothing beyond.
   **The user initiates that planning session; the orchestrator never starts it on its own.**

## What to implement next

The program is **design-complete**, has passed the pre-implementation scientific review
(`docs/cosmology/scientific_review.md`), and is in implementation.

> **State, 2026-09-12. Wave 0, its completion wave and the hardening pass are COMPLETE, and
> nothing is open in either.** Seven prompts merged: #524 M0 packaging (`654b627a`), #525 M0.5
> with 185 frozen fixtures (`e310e125`, now **197 fixtures over 49 pairs** since #547 derived
> three of the four excluded theories), #526 PSALTer installed (`c8c57251`), I-REM instruction
> sites and the oracle CI gate (`32d21a3b`), I-533 retiring the drop rows (`df206443`), I-543
> resolving the Tier-1 gate (`7286ae94`), and I-ONB one onboarding path (`91e300af`). Trunk is
> **v0.54.2**, `CI 34714779268` and `CI 34714779285` both success on `456c1c15`.
>
> **The PSALTer install is certified and rebuild-proof**: Wolfram 14.3.0 × xAct 1.3.0 × PSALTer
> `bb45adb0` × local registration of two Function Repository resources the package depends on
> but never declares, under fixed UUIDs. `VERDICT: MATCH`, re-run from scratch from an **empty**
> registry with provenance asserted, independent of any cloud login. **Nothing in Wave 1 is
> blocked, and nothing in Wave 0 is left over.**
>
> Closed 2026-09-12: **#559** (one onboarding path), **#547** (three of the four excluded
> theories derived, the fourth measured un-derivable), **#554** (the oracle's drift classes and
> the `--staleness` detector). Carried forward with owners, none a Wave-1 dependency:
> **#558** (β over recombination — the second half of O4a's precondition; goes with the #503
> session), **#548** (two specs without a TOML, M3), **#560** (the #394 volume-element check
> rejects every FRW measure — **WS2 needs this before trusting `de_sitter` as its oracle**),
> **#561** (`derive` reported success for an aborted derivation — fixed, `bf9361c5`),
> **#562** (Wolfram call forms in `sign_algebra`, M3), **#563** (the undeclared `xPerm.m` patch
> in the reference userbase, M0.5), **#565** (the image's second, cloud-only `wolframscript`,
> M-parallel).
>
> **Wave 1 is planned by a fresh planning session, which the user initiates** — this
> orchestrator never starts it. Composition is settled and recorded so that session begins
> from committed state: **#532 (M1a) ∥ #503 ∥ I-S1A (#527)**. **No Wave-1 prompt exists yet;
> writing them is that session's first task.** Read the wave board above, then the
> wave-boundary checklist, then plan.
>
> **What Wave 1 inherits that the original outline did not:**
>
> - The Stage-1 exporter must **enumerate the private globals from
>   `ParticleSpectrum.m:74-81` at the pinned revision** and guard the set with a fixture test
>   — not hard-code a count, which is how the six-item list came to be wrong.
> - **Two Function Repository resources are load-bearing and undeclared.**
>   `scripts/psalter/ensure_registered.sh` registers them idempotently under fixed UUIDs
>   (`install-psalter.sh` calls it; the `wolfram-objects` volume from #559 keeps the registry
>   across rebuilds), and `verify-wolfram-setup.sh --require-psalter` asserts provenance, so a
>   lost or foreign registration is red rather than a silently wrong spectrum. Run the verify
>   before any gate.
> - **A Tier-1 pass certifies two keys and nothing else** — not source constraints, spectrum,
>   or unitarity conditions. Tier 2 is the physics gate.
> - `A[0]` is answered by **I-S1B** (it needs a real export); I-S1A represents both forms and
>   carries it as a named unknown.
> - The `§5.2` legacy↔new mapping is designed **with the symbolic stage (M3)**, not at M0.5.
> - Deferred with owners: #553 (M1b/M5 orphan inventory), #554 (spec drift, needs the lane),
>   #534 (port manifest at M3), #535/#536 (M3 `inspect`/`validate` requirements), #537 (M7),
>   #529 (WS3 O2 contract, blocks the O2 handoff not Wave 1), #530 (survey tags at WS3),
>   #548 (two specs without a TOML, M3), #558 (β over recombination, with #503),
>   #559 (onboarding, I-ONB, no lane), #560 (the FRW volume-element check, WS2 reads
>   `de_sitter` through it), #561's deferred half (Wolfram-side `Exit[1]` on abort, with #513).

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
