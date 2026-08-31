# Solver design study: time-dependent per-k systems (H3)

**Status:** H3 deliverable, 2026-08-31. Design + research record for WS3 (#492) of the
cosmology program (`docs/COSMOLOGY_PROGRAM.md`). Prototype evidence from three scratchpad
prototypes (steppers registry, O2 Magnus/WKB suite, O3 eikonal engine); numbers quoted
below are measured, not estimated, unless marked otherwise. Claims sourced from the
2026-08-30 web surveys and not yet re-verified against primary sources are tagged
**[survey]**; everything else is verified against local TeX/code (References appendix).

**TL;DR.** WS3 is **two solver problems, not one** (H2 §0.1): an **oscillation-resolving
mode-equation engine** for O2 (gravitational waves, ~1–10³ oscillations per mode — step
through them) and an **eikonal amplitude engine with patch averaging** for O3 (CMB-frequency
photons, ~10²⁶ oscillations — carrier removed symbolically at derivation time, only
slowly-varying amplitudes evolved). Recommended architecture: **our own solver chained to
unmodified CAMB (option iii) for both engines, sharing one numerical core** — O3 is only
possible this way, which settles the choice regardless of the O2 verdict. First
implementation: η-grid assembly + GL4/CF4 Magnus with Duhamel sources (O2), the generic
eikonal front-end (O3), with the matrix-WKB rung gated on bake-off evidence. Every
candidate method ships as a selectable backend and is trialed on the benchmark suite; the
final composition is decided by the experimentation matrix, not by reading.

## Symbol table

| Symbol | Meaning | Units (natural, ħ=c=1) | Computed from |
| --- | --- | --- | --- |
| `η` | conformal time | Mpc (O2) / s (O3 literature) | CAMB background |
| `a(η)`, `ℋ = a'/a` | scale factor, conformal Hubble | —, 1/Mpc | CAMB background |
| `k` | comoving wavenumber (O2 mode index) | 1/Mpc | mode grid |
| `ω` | photon frequency (O3 mode index) | eV | observation bands |
| `M(η)` | the per-mode system matrix | 1/Mpc or eV | WS2-derived coefficients on CAMB tables |
| `N_osc` | oscillation count `∫|Im λ| dη / 2π` | — | classifier |
| `ε_ad` | adiabaticity `max_{i≠m}|Γ_im|/|λ_i−λ_m|` | — | classifier |
| `Γ = V⁻¹V'` | non-adiabatic eigenbasis coupling | 1/Mpc | WKB stepper |
| `Δ_M` | graviton–photon mixing entry, `B_T/M_Pl` | eV | B model + couplings |
| `Δ_γ`, `Δ_g` | photon / graviton refractive entries `ω(n−1)` | eV | plasma, EH, ℋ, ν, μ |
| `ω_pl` | plasma frequency, `ω_pl² = e²n_e/m_e` | eV | CAMB `x_e(η)` |
| `Γ_γ` | photon absorption `σ_T n_e` | eV | CAMB `x_e(η)` |
| `l_osc` | mixing oscillation length `2π/|Δ_osc|` | Mpc | eikonal engine |
| `l_B`, `λ_B` | magnetic coherence length | Mpc | B model |
| `l_t` | field transition/smoothing scale | Mpc | B model |
| `l_free` | photon mean free path `1/Γ_γ` | Mpc | CAMB `x_e(η)` |
| `ν`, `μ`, `c_T` | tensor friction, mass, speed | —, eV, — | theory couplings |
| `κ` | `√(16πG)` (TIDAL convention, `docs/tex/gertsenshtein_formula.tex`) | 1/eV | — |
| `γ_LZ` | Landau–Zener adiabaticity parameter | — | crossing analysis |
| `P` | **conversion probability only** in this doc | — | engines |

`P` is overloaded in the surrounding literature (Mirizzi's `B_⊥` power spectrum; pressure
in the program docs); here it always means a conversion/transition probability.

## 1. Scope, baseline, and the assembly constraint

**Problem.** Evolve small coupled linear systems `y' = M(mode, η) y (+ s(η))`, n ≈ 4–10
components per polarization block, coefficients built from CAMB background tables
(`a, ℋ, x_e → n_e, ω_pl², B̄(a)`), over `η` spanning ~4 decades, inside an inference
loop. Two distinct regimes (H2, `observable_ladder.md` §0.1):

| | O2 (gravitational waves) | O3 (CMB photons) |
| --- | --- | --- |
| mode index | comoving `k ~ 10⁻⁴–1 Mpc⁻¹` | `ω ~ 100 GHz` ⇒ `k ≈ 2×10²² Mpc⁻¹` |
| oscillations over `η₀ ≈ 1.4×10⁴ Mpc` | ~1–10³ | ~10²⁶ |
| steppable? | yes — the mode-equation engine (§8) | no — eikonal amplitude engine (§4, §7) |
| per-call cost profile | budget-critical | `N_ω ~ 50 × N_η ~ 10³` small matrices, ≪ 1 s |

**Budget, derived.** Inference needs ~10⁴–10⁵ likelihood calls; at CAMB's <1 s/call that
is ~3–30 CPU-hours at parity. 10× slower stays feasible; 100× makes sampling impractical.
Gates: **≤1 s/call comfortable, ≤10 s acceptable, >10 s fail** (the 10×/100× rule of
`COSMOLOGY_PROGRAM.md`).

**The CAMB null hypothesis.** CAMB solves the O2-class problem by brute force: explicit
adaptive Dormand–Prince 5(4) (`RungeKuttaDP45`, `equations.f90:225` scalars, `:488`
tensors) resolving every oscillation, made affordable by approximation switching (TCA/RSA
remove oscillatory equations once they stop mattering) and by integrating sources only
where they contribute. "Do what CAMB does" — adaptive RK per mode, or batched fixed-step
Magnus, plus RSA-analogue freezing of dead modes — is therefore the **baseline every
fancier method must beat in the bake-off** (§11), not a strawman. The WKB rung is an
improvement opportunity, not a necessity.

**The assembly constraint.** In the legacy modal solver's measured 82.6 ms forward model
(`docs/MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md`), the matrix exponential is **~0.3 ms
(~0.4%)**; coefficient evaluation + matrix assembly is ~49 ms and the stability probe
~13 ms. Any time-dependent scheme that rebuilds `M(η)` through an assembly path of that
shape at every step multiplies the dominant cost by the step count and is dead on
arrival. The design consequence, shared by both engines: **coefficient functions are
evaluated once per likelihood call on the η-grid** (per segment, per mode bucket) into
batched arrays, and the steppers consume arrays, never evaluators. This is the single
strongest argument for a shared numerical core (§7).

Batched-`expm` throughput measured on this machine (scipy loop, float64): ~20–70 µs per
`expm` at n = 4–10. Naive full-resolution rung 2 at `k = 1 Mpc⁻¹` needs `kη₀/π ≈ 4500`
steps; 10³ modes × 4500 steps × ~30 µs ≈ **2 minutes/call — over budget**. With per-k
step counts (`N_steps(k) ∝ k`, log-spaced k), the total collapses by ~an order of
magnitude to the "acceptable" band; the WKB rung (§9) and mode freezing take it further.
Step counts scaling with k, not the exponential's unit cost, are the O2 cost driver.

**Relation to the legacy solver and H4.** The new solver is a **clean design from the end
goals** — nothing is ported from `tidal/solver/modal.py`, per H4's "not ported, not an
oracle" (`tidalcosmo/solver/README.md`) and the program's direction. Legacy anchors used
here are evidence and pitfalls only: the cost decomposition above; the **silent-freeze
bug** (`solve_modal` never re-checks `time_dependent` — called directly on a
time-dependent spec it silently freezes coefficients at t=0, `modal.py:444-448` gate is
CLI-only); the augmented-matrix Duhamel precedent (`modal.py:3775`); and the negative
result that JAX-on-CPU is 3× slower than the scipy loop for exactly this workload
(`docs/meetings/2026-05-22_supervisor.md`). The new implementation lands in
`tidalcosmo/solver/` behind H4's seams (`tidalcosmo/background/protocol.py` is the
background-table swap point) under the #367/#379 test contract (every dispatch path
consumes the kinetic matrix identically). `tidalcosmo/solver/README.md` expects revision
by this study; that revision is a follow-up issue, not part of this commit.

**Symbolic-manipulation policy (program decision, 2026-08-31).** Wolfram runs **once, at
derivation time** — EOM derivation plus PSALTer Lagrangian input; no Wolfram license is
required during sampling. Runtime code evaluates numerical coefficient tables. Where
runtime symbolic work is genuinely needed, sympy on proper expression trees is
permitted if fast; **string-based symbolic manipulation is banned** (the historical
fragility). Every solver path must accept *any* WS2-derived system — no hand-crafted
routes bespoke to one theory; literature closed forms are validation oracles only.

## 2. Prior art

Three surveys (2026-08-30) plus targeted verification. Verdicts are **provisional by
design**: each promising direction gets a real minimal implementation trialed on the §10
suite (§11) before the final architecture is composed — paper-reading ranks candidates,
experiments decide.

### 2.1 Boltzmann codes, and how sectors get added

| code | integrator | tensors/`C_ℓ` | extension mechanism | verdict for us |
| --- | --- | --- | --- | --- |
| CAMB | DP45 explicit adaptive + TCA/RSA switches (verified in `equations.f90`) | yes | `DarkEnergyInterface` classes (TorC used it, background-only); EFTCAMB fork pattern (`! EFTCAMB MOD` hooks) **[survey]**; MagCAMB adds PMF tensor stress into `rhopi` in `derivst` **[survey]** — the closest tensor-source precedent | the oracle and the hand-back target; not the host (§6) |
| CLASS | `ndf15` implicit NDF (verified, arXiv:1104.2933) | yes | `class_define_index` + `perturb_derivs` (hi_class, GDM) **[survey]** | approximation-scheme lessons (§8); not the host |
| DISCO-EB | Diffrax Kvaerno5, autodiff Jacobian | **scalars only**, no `C_ℓ`, no published timings (verified: the paper declines performance comparison) | JAX state-vector extension; **GPL-3.0** | see GPL decision below |
| SymBoltz.jl | Rodas5P via ModelingToolkit, approximation-free **[survey]** | scalars only, tensors future work **[survey]** | symbolic species composition — the closest analogue to TIDAL's Lagrangian→JSON pipeline | design inspiration for the WS2→solver interface |
| ABCMB | JAX, fluid classes **[survey]** | no B-modes **[survey]** | "add a fluid without opening a source file" — but two-way coupled by design | API inspiration; two-way coupling violates the spectator premise |
| Bolt.jl / PyCosmo / clax | KenCarp4 / LSODA / Kvaerno5 **[survey]** | clax has primordial BB (unreviewed) **[survey]** | various | context only |

**Two structural findings.** (1) **No existing code evolves graviton–photon or
axion–photon mixing inside a Boltzmann solver** — every CMB-facing conversion paper
post-processes along the line of sight. (2) The CAMB tensor equation, from the code
(`derivst`): `σ' = −2ℋσ + kH_T − ρπ/k`, `H_T' = −kσ`, i.e.
`H_T'' + 2ℋH_T' + k²H_T = κa² Σ_i ρ_i π_i` with `κ = 8πG` and the convention
`h_ij = 2H_T Q_ij` — the source normalization is taken from the code, not from a
textbook, because the `Q_ij` normalization is code-specific.

**DISCO-EB GPL decision (routed here by H4, `repo_reshape.md`).** DISCO-EB is GPL-3.0.
Decision: **read it, learn from it, and be influenced by its design freely — only literal
code copying triggers GPL, and that is what we avoid.** (D6's copy-freely provenance rule
applies to Barker's code, which has explicit permission — not to DISCO-EB.)

### 2.2 Oscillatory and adiabatic matrix ODE integrators

Beyond the ladder's original five sources (Magnus review 0810.5488; oscode 1906.01421;
riccati 2212.06924; Ioannisian–Smirnov 0803.1967; Haddadin–Handley 1809.11095):

| method | what it buys | limitation | verdict |
| --- | --- | --- | --- |
| Iserles modified Magnus (midpoint interaction picture, BIT 42 2002) | step set by variation of M, not `|λ|`, when the midpoint rotation removes most of the phase | fails when eigenvectors rotate fast | registry backend `ip-magnus` (prototyped) |
| Blanes–Moan CF4/CF6 (APNUM 56 2006); Alvermann–Fehske CFETs to order 8 (arXiv:1102.5071) | commutator-free: cheap for small n, structure-preserving, works for non-Hermitian M | fixed order; still `O(k)` step counts | registry backends `cf4` (prototyped), `cf6` (coefficients acquired) |
| **Lorenz–Jahnke–Lubich adiabatic Magnus (BIT 45 2005)** **[survey — journal source, bounds not re-derived]** | time-varying eigendecomposition; O(h²) uniform in the adiabaticity parameter for `h < √ε` under a gap condition; one diagonalization/step | loses accuracy in avoided-crossing zones; constants ∝ 1/gap² | **the §9 template** (prototyped as `adiabatic-magnus`) |
| Hu–Bremer cyclic-vector phase functions (arXiv:2309.13848) | the **only** published frequency-independent matrix method | conditioning collapses for n ≥ 5; needs distinct eigenvalues | reference oracle on n ≤ 4 blocks only |
| RCMS (Degani–Schiff) **[survey]** | frequency-independent error | large part must be *constant* — ours is `k ×` time-varying structure | rejected |
| Exponential integrators (Hochbruck–Ostermann) | — | reduce to Magnus/CF for linear non-autonomous systems | nothing new for us |
| Magnus for neutrinos: Casas–D'Olivo–Oteo (arXiv:1611.06814), QKE codes (arXiv:1608.01336, 2406.09504) | ~100× over generic integrators on a 3×3 oscillation problem; unitary at every truncation; scales in production | problem-specific tuning | the strongest published evidence the Magnus route pays off on our problem shape |
| nuSQuIDS (arXiv:2112.13804) | interaction picture w.r.t. the fast diagonal; non-Hermitian absorption + incoherent source coexisting on one density matrix | neutrino-specific | both patterns copied by the O3 engine (§7) |

**The gap is real, and it is the group's own gap.** The WS3 ladder's named authors were
checked directly: oscode is Agocs–Handley–Lasenby–Hobson; rung 4 (Haddadin & Handley) is
the Handley group's own line of work. Their subsequent output through 2026 contains **no
matrix/coupled follow-up**: Bamber–Handley (arXiv:1907.11638) state RKWKB is "essentially
one-dimensional", and the oscode paper itself records that the coupled extension proved
"more difficult than anticipated" (`main.tex:607`). **No published application of Magnus,
WKB-switching, or adiabatic integrators to cosmological perturbations was found** (nearest
neighbors: the neutrino QKE Magnus works above). Both gaps — a matrix RKWKB solver, and
Magnus for cosmological perturbations — remain open and publishable, and rung 4 being the
group's own method makes collaboration natural.

### 2.3 Mixing/conversion prior art and the averaging criteria

The eikonal reduction originates with Raffelt–Stodolsky (PRD 37, 1988): factor the wave
operator `(ω² + ∂_z²) = (ω + i∂_z)(ω − i∂_z)`, replace `(ω + i∂_z) → 2ω` on the forward
branch, giving `i∂_z A = (ω + M/2ω) A`; valid when every `|M_ij| ≪ ω²` and the envelope
varies slowly on the wavelength. Ejlli (arXiv:2004.02714, local) checks this "SVEA"
against the exact second-order solution. Cembranos et al. (arXiv:2302.08186, local) is the
FRW statement our O3 block uses (§4). gammaALPs (arXiv:2108.02061) supplies the
engineering pattern: ordered environment list, closed-form `expm` per constant cell,
density-matrix propagation, batched random-field realizations.

**The published regime chain** (each row a validity criterion, → §4's decision table):

| condition | method | source |
| --- | --- | --- |
| cell ≪ `l_osc`, given realization | coherent transfer-matrix product | gammaALPs |
| Born (`Δ_M·Δz ≪ 1`) ∧ `N·P ≪ 1` ∧ random patch orientations | incoherent rate sum; continuum `⅓[1−e^{−3Pz/2s}]` beyond | Mirizzi–Raffelt–Serpico arXiv:0704.3044 App. A |
| `l_osc ≲ l_t` (smooth field) | adiabatic propagation eigenstates / WKB in the *mixing* eigenvalue; hard cells wrong by ~10⁴ | Kartavtsev–Raffelt–Vogel arXiv:1611.04526; Carenza–Marsh arXiv:2302.02700 |
| `l_free < l_osc` | density matrix with anti-Hermitian Γ | Dolgov–Ejlli arXiv:1211.0500 §5 |

Graviton–photon prescriptions and traps: Domcke–Garcia-Cely (arXiv:2006.01161, the rate
integral O3 validates against); Ito–Kohri–Nakayama (arXiv:2309.14765, the conservative
per-interval integrand switch — drop the `N_G` enhancement once `l_osc < l_G`);
Fujita–Kamada–Nakai (arXiv:2002.07548, patch = photon mean free path in the early
universe); Chiba–Jinno–Nomura (arXiv:2505.10926, Born kernels giving expectation **and
variance** directly from `P_B(k)` — the route to a likelihood covariance); Addazi et al.
(arXiv:2401.15965, parametric resonance that the sin²→½ average misses when
`λ_B ≳ 10 l_osc` — a recorded blind spot); Pshirkov–Baskaran (arXiv:0903.4160, the
`N·P` error the flags must prevent).

## 3. The regime continuum and the classifier

**Regime is a property of an eigen-direction — not of a theory, and not even of a block.**
Mixing conserves comoving momentum, so partners share `k`, but *not* frequency: the
same-k photon partner of a CMB-scale GW carries an effective mass `ω_pl ≫ k` and
oscillates at ~`ω_pl`, not `k`. A block that mixes a slow graviton direction with a fast
heavy-photon direction is therefore not a collision between the two engines — it is
exactly what the eigenbasis split (§9) handles: fast directions propagate adiabatically,
slow ones are stepped, and `Γ = V⁻¹V'` carries the physics between them.

**Phase representability.** A direction whose accumulated phase `∫|Im λ| dη` exceeds
~10¹² cannot be carried in float64 even adiabatically (oscode's own documented limit —
accuracy is needed mod 2π). Such a direction must either (a) provably decouple below
tolerance and be dropped, with a flag, or (b) live in a rotating frame where only phase
*differences* appear — which is precisely what the WS2 eikonal reduction provides at
derivation time. **A spec whose classified `N_osc` exceeds float capacity and which
carries no `eikonal` section is refused with a hint pointing at the WS2 reduction — never
stepped silently.** (Prototype evidence: at `N_osc = 10⁶`, per-step GL3 quadrature of the
phase in float64 holds the accumulated phase error to ~10⁻⁹ rad mod 2π — float64 phase
integration is safe well past the O2 band; the 10¹² wall is real but far away.)

Three regimes on one axis:

| regime | criterion | engine | carrier handling | fallback |
| --- | --- | --- | --- | --- |
| full stepping | `N_osc` modest; near/below horizon crossing | Magnus ladder (§8) | resolved | rk-adaptive baseline |
| adiabatic/WKB split | `kη ≳ few`, `ε_ad ≪ 1`, gaps open | matrix WKB (§9) | fast directions: phase integral | Magnus segment on trigger |
| eikonal | all directions fast; `N_osc` beyond stepping | amplitude engine (§7) | removed symbolically (WS2) | refuse without `eikonal` section |

**Classifier** (`classify_block(spec, background_tables, mode)`), deliberately **soft at
this stage**: the solvers we end up building will have their own empirical regimes, so the
classifier ships as *logged heuristics with a `--solver` override*, its thresholds set and
revised by the bake-off table (§11) — hard gates are deferred until that table exists.
It computes per block/mode: `N_osc`, damping scale `max|Re λ|`, `ε_ad`, gap statistics,
presence of an `eikonal` spec section, and whether a stochastic background is declared.
Oscillation (`Im λ`) and damping (`Re λ`) are classified **separately** — `|λ|` ratios
are ill-defined where an eigenvalue crosses zero, which happens for *every* O2 mode
super-horizon. Its decision is logged in the run record beside `gauge_certificate`.
Proposed home in the H4 layout: `tidalcosmo/validity/`.

The few **solver-independent structural facts** it can rely on now:

- Every O2 companion block passes through an exact Jordan (defective) point at horizon
  crossing (`k = ℋ`; for the de Sitter block, `kη = −1`). **The Magnus rung always owns
  the frozen/crossing regime; WKB is entered only at `kη ≳ few`.** Prototype evidence:
  eigenbasis conditioning alone does *not* reliably detect the crossing (a straddling
  step returned an 11% error without tripping a `cond(V) > 10⁸` guard), while the
  `‖Γ‖·h` trigger hands the whole segment to Magnus and lands at 1.7×10⁻¹³.
- The state is nondimensionalized (`y = (h, h'/k)`) so `cond(V)` does not scale with k.
- The phase-representability rule above.

## 4. The coupled blocks, concretely — including sources

### O2 block (tensor channel)

`h'' + 2ℋh' + k²h = S_std + S_new`, coupled to the torsion 2⁺ modes; first-order state
`y = (h, h'/k, δT, δT'/k)`. Reduced-form oracle: Cembranos eq. 20,
`h'' + 2ℋ(1+ν)h' + (c_T²k² + μ²)h = 0` — the derived block must reduce to this with
`ν, μ, c_T` identified from the Lagrangian couplings (the O2 known-answer check in
`observable_ladder.md`).

**Sourced evolution is first-class.** CAMB's tensor equation carries the photon/neutrino
anisotropic-stress source (`ρπ/k` in `derivst`), which damps tensor amplitudes at the
percent-to-tens-of-percent level around horizon crossing; the O2 gate (sub-percent BB
agreement with CAMB) is **unreachable with a homogeneous-only formalism**. Design:
variation of constants `y(η) = U(η,η₀)y₀ + ∫ U(η,s) s(s) ds`, implemented with the
augmented-matrix exponential (`Z = [[M, s],[0, 0]]`; Al-Mohy–Higham; legacy precedent
`modal.py:3775`) which every Magnus-family stepper inherits without modification.
Prototype evidence (benchmark A8): a forced oscillator through the augmented matrix
reproduces the analytic particular solution to 2.7×10⁻¹¹ (magnus4, 2000 steps).

**Dense output.** WS4's line-of-sight integral needs `h, h'` on the whole source grid
convolved with the visibility function — not one end-to-end propagator. The stepper
interface therefore returns per-step propagators, with step boundaries aligned to (or
interpolated onto) the LOS η-grid.

**Initial conditions.** Deep-radiation adiabatic ICs (`kη ≪ 1`, `h → const`) sit in the
regime where WKB is invalid; **every mode starts on the Magnus rung**, with the
WKB handover at the classifier's `kη ≳ few` gate.

### O3 block (photon channel)

Eikonal amplitude system per polarization λ (Cembranos eqs. 25–27, local):
`(ω + i∂_η + M)ψ = 0`, `ψ = (A_λ, h_λ)`, with

```
M = [[Δ_γ, Δ_M], [Δ_M, Δ_g]],   Δ_M = B_T/M_Pl,   Δ_γ = ω(n_γ − 1)  (plasma + EH),
Δ_g = ω(n_g − 1) = iℋ(1+ν) − μ²/2ω − (A/2)ω^{α−1}   (complex in general)
```

Two convention findings from the prototype (benchmark B3), worth recording because the
surrounding literature is a documented normalization minefield
(`docs/tex/gertsenshtein_formula.tex`):

- **`M_*` in Cembranos' `Δ_M = B_T/M_*` is the non-reduced Planck mass** and the eq. 27
  normalization (no /2) is the one their numbers use: with `Δ_M = B_T/M_Pl`,
  `M_Pl = 1.22×10²⁸ eV`, their eq. 37 benchmark reproduces as `P = 8.7×10⁻¹²` vs their
  `8.6×10⁻¹²` (1.7%, inside their rounding); the reduced-mass and /2 variants miss by
  4–25×.
- The `n_e` exponent in their numeric `Δ_γ` formula (their eq. 36) is a typo: `ω_pl² ∝
  n_e` linearly, and the printed prefactor `1.05×10⁻⁶ (n_e/cm⁻³)(eV/ω) s⁻¹` is correct
  for the linear dependence (verified against `ω_pl² = e²n_e/m_e` directly).

**Density-matrix form**, written exactly (Cembranos §III):
`i dρ/dη = [M₊, ρ] + {M₋ − iΓ/2, ρ}` with `M = M₊ + M₋` the Hermitian/anti-Hermitian
split and `Γ_γ = σ_T n_e` the Thomson absorption — *not* the shorthand
`−i[M,ρ] − {Γ,ρ}`. Required whenever `l_free < l_osc` (§2.3). Prototype evidence
(benchmark B4): wavefunction and density-matrix evolution agree to 3.7×10⁻¹⁴ at Γ = 0 on
a time-dependent complex-symmetric M.

**Eikonal validity flags**, evaluated at runtime from WS2-exported dropped terms and the
background: `|M_ij|/ω² ≪ 1`, envelope slowness, `P ≪ 1`, `l_osc ≪ Δz`,
`P·N_patches ≪ 1`, `|dropped|/|kept|` per η.

### Patch averaging — an auto-selected decision chain

The B̄ background comes in two supported modes:

- **Uniform B̄** (used by parts of the literature, e.g. Cembranos' single-patch
  calculation): homogeneous, preserves exact k-decoupling — the validation and
  idealized-run mode.
- **Stochastic B̄** (the realistic PMF, `magnetic_field_background.md`): breaks strict
  homogeneity; rather than couple k-modes, the engine handles it statistically, selecting
  per η-interval from the published regime chain (the Ito–Kohri–Nakayama integrand-switch
  pattern):

| condition (evaluated on `l_osc(η), l_B(η), l_t, l_free`) | method | runtime flag | known blind spot |
| --- | --- | --- | --- |
| cell ≪ `l_osc` | coherent transfer-matrix product over cells (gammaALPs pattern: closed-form `expm` per cell) | — | cost ∝ N_cells |
| Born ∧ `N·P ≪ 1` ∧ random orientations | incoherent rate sum `𝓕 = ∫⟨P⟩/Δz`, `Δz = min[λ_EQ, λ_B⁰]/(1+z)` | `N·P ≪ 1` | resonance at `λ_B ≳ 10 l_osc` missed by sin²→½ (arXiv:2401.15965) |
| `l_osc ≲ l_t` | adiabatic propagation eigenstates / WKB in the mixing eigenvalue | `γ_LZ` regime | LZ invalid for `γ_LZ > 1` with near boundaries (arXiv:2302.02700) |
| `l_free < l_osc` | density-matrix form | — | — |

Optional: the Chiba–Jinno–Nomura kernel route (expectation *and variance* directly from
`P_B(k)` under Born) feeds the patch-realization variance into the likelihood covariance
rather than only shifting the mean. Every selection is logged.

### One-way vs two-way coupling, per engine

- **O2:** one-way (CAMB → us) is valid iff the new sector's back-reaction on the standard
  modes sits below CAMB's own accuracy — the same class of approximation CAMB itself
  makes when it drops the anisotropic-stress source where negligible (Seljak–Zaldarriaga
  `los.tex:355-362`). Gated by the program's growth-impact monitor; two-way would be
  needed only if `|S_new| ≳ |S_std|` at the source level, which contradicts the spectator
  premise.
- **O3:** structurally one-way — a per-frequency depletion/rotation operator applied to
  CAMB's output; two-way would require conversion to drain photon energy at
  recombination-relevant levels, excluded by `P ≪ 1`.

### The CAMB hand-back, concretely

- **O2 →** `h(k,η), h'(k,η)` on the LOS grid, handed to **WS4's own tensor line-of-sight
  code** (`tidalcosmo/observables/`); CAMB's Fortran is an oracle, never the integrator.
  Assembly: `C_ℓ^BB = C_ℓ^{lensed-scalar, CAMB} + C_ℓ^{tensor, ours}`.
- **O3 →** a band-averaged operator `P(ω, n̂)` (plus patch variance where computed)
  applied by WS4 to CAMB Stokes/`C_ℓ`.
- Both front-ends are **coupling-agnostic**: O4-mix (#511) shares the O3 engine, so the
  Gertsenshtein-specific pieces (the stochastic-B patch model) are cleanly separated from
  the amplitude-system evolver.

## 5. The WS2 symbolic interface for the eikonal path

The reduction is derivation-time work and lives in the Wolfram stage — consistent with
the once-at-derivation policy (§1). It must be **theory-agnostic**: one reduction stage
(candidate `EikonalReduce.wl`) operating on any derived second-order FRW system, never a
per-theory hand derivation. The symbolic stage may emit **multiple derived
representations** from one derivation (as PSALTer's pipeline returns several objects):
the second-order system, the first-order eikonal amplitude system, and any further
solver-facing forms found useful — each its own JSON section, so solver back-ends pick
what they need without re-running Wolfram.

Deliverables for WS2:

1. **The reduction stage**: the Raffelt–Stodolsky factorization applied symbolically —
   carrier ansatz `Φ = e^{−iωη}ψ(η)`, factor the wave operator, truncate the forward
   branch at first order in `1/ω`.
2. **A JSON `eikonal` section** alongside `equations[]`: carrier `ω`; state vector;
   `M` coefficient expressions in `(η, a, ℋ, n_e, B, ω)`; the Hermitian/anti-Hermitian
   split; the polarization basis — conforming to H4's conventions spec (#513).
3. **Dropped terms exported symbolically**, so the engine can evaluate
   `|dropped|/|kept|` per η as a runtime validity flag (soundness rule: flag, never
   assume).
4. **Acceptance test**: the reduction of the graviton–photon FRW system reproduces
   Cembranos eqs. 24–25 (including the `Δ_M = B_T/M_Pl` normalization pinned down in §4)
   on their setup.

Runtime consumption reuses the numeric coefficient-table path; no Wolfram at sampling
time, no string surgery.

## 6. The integration-target decision

Three architectures, five criteria, assessed against **both** engines:

| criterion | (i) patch CAMB Fortran | (ii) build on DISCO-EB | (iii) own solver + unmodified CAMB |
| --- | --- | --- | --- |
| correctness risk | re-derive TCA/RSA analogues per new stiff coupling; silent two-way coupling violates the spectator premise | scalars only today — we would add the tensor sector *and* LOS ourselves | risk concentrated in our own well-tested small-system code |
| per-call cost | native speed | GPU-oriented; our CPU-JAX measurement is 3× slower; no published timings | §1 numbers; within budget with k-scaled steps |
| extensibility to arbitrary TIDAL sectors | per-theory Fortran regeneration | JAX state-vector edits per theory | generic by construction (§1 policy) |
| maintenance | fork tracking (EFTCAMB shows it is possible, at real cost) | GPL-3.0 + young codebase | ours; CAMB stays stock (program's CAMB policy) |
| coupling handling | two-way "for free" — which is exactly wrong here | two-way by design | one-way by construction, monitor-gated |
| **O2 verdict** | plausible (MagCAMB precedent) but wrong trade | poor fit today | **fits** |
| **O3 verdict** | **structurally impossible** — CAMB has no per-frequency photon propagation | same | **the only option** |

**Decision: (iii) for both engines, sharing one numerical core.** The headline finding is
the asymmetry: **O3 forces (iii) regardless of the O2 verdict.** (ii) remains the fallback
if gradient-based sampling ever becomes a requirement (autodiff is DISCO-EB's genuine
advantage); (i)'s MagCAMB-style `rhopi` hook remains documented as the escape hatch if
two-way tensor coupling is ever *wanted*.

## 7. Shared core + two front-ends

**Shared core** (`tidalcosmo/solver/`):

- **η-grid coefficient assembly**: evaluate WS2 coefficient expressions once per
  likelihood call on the segment grids, into batched `(N_nodes, N_mode, n, n)` arrays.
  Memory estimate: a fully dense outer product at worst case is 1.6–2.6 GB. Chunking per
  η-segment and per k-bucket is the default (different k want different η-grids anyway),
  but **memory cost never disqualifies a route** — production runs land on HPC nodes with
  hundreds of GB, so the dense path stays available where simpler or faster.
- **Batched small-matrix kernels**: expm/eig/solve over a leading batch axis. Baseline is
  the scipy loop (measured §1); the JAX vmap pattern is reserved for a future GPU path
  (CPU measurement: 3× slower).
- **The stepper registry** (§11): all steppers consume node-sampled `M` arrays; adaptive
  control lives *outside* the steppers so batching stays legal. Magnus applies unchanged
  to complex non-Hermitian `M` (bounded finite-dimensional systems).
- **Validity-flag plumbing** shared with the classifier and `gauge_certificate`-style
  run-record logging.

**O2 front-end**: mode index k; real `M` plus source vector (augmented-matrix Duhamel);
ladder §8; per-step dense output onto the LOS grid; hand-back per §4.

**O3 front-end**: mode index ω (~50 values); complex symmetric `M`; wavefunction and
density-matrix forms; per-patch probabilities + the §4 averaging chain; hand-back per §4.

Why two front-ends over one engine: different mode loop (ω vs k), closed vs open system,
transfer function vs rate integral as the observable. Why one core: identical batched
small-matrix numerics, identical assembly constraint, identical flag plumbing.

## 8. The O2 ladder, rung by rung

1. **Exponential midpoint** `expm(h·M(η+h/2))` — 2nd order (measured 2.0); the `Ỹ₀` of
   the modified-Magnus picture. Debug rung; near drop-in.
2. **GL4 Magnus** `Ω = (h/2)(A₁+A₂) − (√3/12)h²[A₁,A₂]` — 4th order (measured 4.0);
   reduces exactly to `expm(Mh)` for constant M (measured: 7×10⁻¹⁶). Convergence needs
   `∫‖A‖ < π` (sufficient, not sharp), and Moan's bound `O(h^{p+1}λ^{p/2−1})` means a
   6th-order scheme degrades to 2nd at `h²λ ≈ 1` — **Magnus fixes time dependence, not
   oscillation**; `N_steps ≈ kη₀/π ≈ 4500` at `k = 1 Mpc⁻¹`. Variants: **CF4**
   (commutator-free, two exponentials; coefficients `(3 ± 2√3)/12` — the Magnus review's
   printed `½ ± √3/72` fails the sum-to-½ and 4th-order conditions and is a misprint;
   measured order 3.9 with the corrected values, and the working composition applies
   `exp(h(α₂A₁+α₁A₂))` first); **GL6** (measured order 5.9) for smooth early segments.
   **Validated: necessary, insufficient alone.**
3. **Matrix WKB** — §9. The research rung; gated on bake-off evidence (§12).
4. **Piecewise-analytic transfer matrices** (Haddadin–Handley) — requires a frictionless
   *scalar* form; for a matrix block it applies only per eigenmode after
   diagonalization, i.e. it is rung 3's adiabatic propagator with Airy/Bessel phase fits
   instead of quadrature. **Deferred**; revisit if rung 3's phase integrals dominate.
5. **Emulation** — last resort only.

**The tight-coupling-analogue question** (CLASS lessons, arXiv:1104.2933): their measured
tables say a stiff integrator buys more than higher-order TCA (rk 1069 s → 19.4 s with
first-order TCA; ndf15 16 s → 14.9 s), and the *large* win is RSA-style replacement of
oscillatory equations by non-oscillatory particular solutions once `kτ ≥ 100`. Our
analogues, both **auto-detected fast paths per D7, never scoping decisions**: strong
coupling `|M_ij| ≫ k` is absorbed by the eigenbasis split (fast eigenmodes propagate
exactly — no stiffness); negligible/decaying eigen-directions are frozen below tolerance
(the RSA analogue). The O3 eikonal reduction is the RSA idea taken to its limit — the
carrier removed *analytically at derivation time* — which is exactly why O3 is cheap.

## 9. Matrix-WKB design

Template: **Lorenz–Jahnke–Lubich adiabatic Magnus** lifted to first-order `y' = M(η)y`
with non-normal `M = V (Λ) V⁻¹` (general eigendecomposition instead of orthogonal Q;
`cond(V)` guarded; Λ complex), cross-checked against Ioannisian–Smirnov's 2×2 closed
forms. LJL's error model — O(h²) uniform in the adiabaticity parameter for `h < √ε`
under a gap condition — is the accuracy expectation **[survey — journal source]**;
`‖Γ‖ = ‖V⁻¹V'‖` spikes are the switch signal (their own observation at avoided
crossings, and our prototype's).

```text
per step [η, η+h], batched over the segment's mode bucket:
  M_j = M(η + c_j h) at GL3 nodes;  Λ_j, V_j = eig(M_j)
  continuity: match eigenvector columns to the previous node by max overlap
              (greedy assignment), fix phase so the overlap is real-positive
  Γ = V_mid⁻¹ V'   (V' by finite-difference stencil across the nodes)
  Φ_osc = diag ∫ i·Im Λ ds   (GL3 quadrature)      # oscillatory phase only
  damping ∫ Re Λ ds applied as a separate diagonal factor   # see below
  B(s) = −e^{−Φ_osc} Γ e^{Φ_osc}    # entries oscillate via Im-eigenvalue
                                    # differences only; ‖B‖ ~ ‖Γ‖, not ~k
  Ω_B = GL4 Magnus of B             # converges iff ∫‖Γ‖ < π (the adiabatic-
                                    # picture convergence statement)
  U_step = V(η+h) · e^{Φ} · e^{Ω_B} · V(η)⁻¹
```

Design decisions, each with its failure mode:

- **Complex-eigenvalue (damping) split.** With the full `Φ` in the similarity transform,
  `B_ij` would carry `exp(∫Re(λ_j − λ_i))` — exponential growth for unequal damping
  rates, destroying `∫‖B‖ < π`. Only the oscillatory phase enters the transform; damping
  is applied as a separate diagonal factor, with `h` capped by `|Re Δλ|·h ≲ 1` where
  damping-rate differences are large.
- **Switching criterion**: adiabaticity `ε_ad = max_{i≠m}|Γ_im|/|λ_i − λ_m|` (the matrix
  form of Ioannisian's `θ̇/Δ`); take a WKB step of size `h_ad = min_i|λ_i|/|λ_i'|`
  (riccati's `h_osc`) iff `h_ad > 5π/max|Im λ|` and `ε_ad·h_ad ≪ 1`; otherwise a Magnus
  step at `h ≤ π/‖M‖`. Error estimate: Ω_B at order 1 vs order 2 plus the GL quadrature
  difference (the oscode two-error pattern); accept iff the predicted next step grows.
  **Prototype evidence on the trigger threshold**: with `‖Γ‖·h ≤ 0.02` the switching
  driver reproduces dense-reference Landau–Zener transition probabilities to 0.1–0.5%
  with ≥95% of steps adiabatic; at a loose 0.15 threshold the error reaches 17%. The
  threshold **must** come from the bake-off, not from a guess — measured sensitivity is
  the proof.
- **Degeneracy and crossings.** The *generic* case is horizon crossing (§3) — owned by
  Magnus, entered by classifier gate, not by runtime detection (prototype: `cond(V)`
  guards alone miss a straddled Jordan point and silently return an 11% error; the `Γ`
  trigger hands over and lands at 1.7×10⁻¹³). Genuine mid-flight crossings (our systems
  can have them): detected by gap collapse or assignment ambiguity → cluster the
  near-degenerate eigenvalues and keep `Φ` block-diagonal (block-Schur; `Γ` is never
  divided by a small splitting), or drop the segment to Magnus. Loud, counted fallback —
  never silent.
- **Constant-M limit**: `Γ = 0`, `Ω_B = 0` ⇒ `V e^{Λh} V⁻¹ = expm(Mh)` to machine
  precision (measured: 8×10⁻¹⁵) — the WS3 gate.
- **Batching**: fixed-step-per-segment keeps batching legal; segment the η-range per mode
  a priori from `λ(η)`, bucket modes by step count.
- **k-independence** (the point of the rung): measured on the de Sitter adiabatic band at
  a fixed 60 steps, the adiabatic stepper's error is *identical* for k = 10, 100, 1000
  (the problem is self-similar in kη — cost per accuracy is k-independent), while GL4
  Magnus at the same step count is useless (error ~3). Absolute accuracy at few steps is
  second-order-limited; the LJL h² correction terms are the known upgrade path.

## 10. Benchmark protocol

Accuracy thresholds derive from the program gates, not assertion: sub-percent `C_ℓ^BB`
(2 ≤ ℓ ≤ 500) and LiteBIRD `δr < 0.001` imply a per-mode amplitude error budget of
~10⁻³ relative through the tensor band; solver-level targets sit two orders below that
budget so solver error never dominates. Timing gates: §1. Methodology:
`OMP_NUM_THREADS=1`; assembly cost reported separately from evolution; `N_k ∈ {100, 300,
1000}`; the CAMB reference (`get_tensor_cls`) measured on the same machine.

**O2 suite** (P = prototype evidence already in hand):

| # | system | reference | pass criterion | status |
| --- | --- | --- | --- | --- |
| A1 | constant M, n = 4, 6 | `scipy.linalg.expm` | rel L2 < 10⁻¹², every backend | **P: 10⁻¹⁵–10⁻¹⁶ all six** |
| A2 | de Sitter tensor, `ℋ = −1/η` | `h = (1−ikη)e^{ikη}` | measured order matches nominal; err within budget | **P: orders 2.0/3.9/4.0/3.9/5.9 (midpoint/cf4/magnus4/ip/magnus6)** |
| A3 | stiff graviton–torsion toy, `a²m²/k² → 10⁴`, coupling ε | converged GL6 | transition amplitude + symplectic drift < 10⁻⁶ | **P: magnus4 1.6×10⁻¹² @ 4000 steps** |
| A4 | `kη ≫ 1`: Cembranos eq. 20 with ν(η), μ(η) on a CAMB-like `a(η)` table | `oscode`/`riccati` (pip; fallback: tight-tol `solve_ivp`) | WKB cost k-independent; Magnus cost ∝ k | partial (see A2b) |
| A2b | de Sitter adiabatic band, fixed 60 steps, k ∈ {10, 10², 10³} | exact solution | error k-independent | **P: identical error across k; magnus4 useless at same steps** |
| A5 | Landau–Zener crossing | dense GL4 reference + analytic `exp(−πΔ²/λ)` | driver within 1% of dense ref; fallback provably triggers | **P: 0.1–0.5% at `‖Γ‖h ≤ 0.02`, ≥95% adiabatic steps; LZ formula itself 5% off at Δ=0.75 (its own asymptotic validity — the Carenza–Marsh point, reproduced)** |
| A6 | phase fidelity at `N_osc = 10⁶` | 40-digit mpmath quadrature | phase error mod 2π ≪ 1 | **P: ~10⁻⁹ rad at 50 GL3 steps** |
| A7 | Jordan point at horizon crossing | Magnus reference | classifier hands over; no NaN; loud | **P: Γ-trigger → 1.7×10⁻¹³; cond-only guard silently 11% wrong (finding)** |
| A8 | **sourced** GR tensor equation with CAMB's `π_ν` table | CAMB's own `h(k,η)` | sub-percent — *the* known-answer gate | augmented-matrix mechanism proven (2.7×10⁻¹¹ on the analytic toy); CAMB-table run is implementation work |

**O3 suite:**

| # | system | reference | pass criterion | status |
| --- | --- | --- | --- | --- |
| B1 | `H → 0`, no plasma, uniform B, Minkowski | TIDAL's own `P = sin²(κB₀D/2)` | < 0.4% | pending (needs WS2 FRW derivation) |
| B2 | detuned flat 2-level | `P = (μ/ω_m)² sin²(ω_m D)` | ~1% | **P: 10⁻¹⁶** |
| B3 | single-patch FRW, GR | Cembranos eq. 34 and `P = 8.6×10⁻¹²` at their eq. 37 parameters | ~10% (their rounding) | **P: 8.75×10⁻¹² (1.7%), with `Δ_M = B_T/M_Pl` pinned (§4)** |
| B4 | density matrix vs wavefunction at `Γ_γ → 0` | identity | < 10⁻¹⁰ | **P: 3.7×10⁻¹⁴** |
| B5 | friction sweep, ν at the GW170817 bounds | `2.09×10⁻¹³` / `1.67` | order of magnitude | pending |
| B6 | patch-averaged LOS, GR | Domcke–Garcia-Cely `𝓕`; He et al. (their `I(1100) ≈ 6.31×10⁶` as evaluated in `observable_ladder.md`, not printed in the paper) | factor ~2 (conventions) | pending |
| B7 | eikonal validity | full second-order solve with ω lowered to ~10³·|M| | error `O(|M|/ω)` | pending |

## 11. The bake-off — every candidate gets implemented and trialed

The final solution is **pieced together from experiments, not from reading**. Each
direction gets a real minimal implementation, is trialed on the §10 suite, and is
retained as a release option; retirement is the classifier's job, never deletion.

- **Interface**: `Propagator.step(M_nodes, source_nodes, h) → (U, particular)` over a
  registry `{rk-adaptive (the CAMB-style baseline), midpoint, magnus4, cf4, cf6,
  magnus6, ip-magnus, adiabatic-magnus, eikonal-cell, hu-bremer (n ≤ 4 reference
  only)}`. Adaptive control lives outside the steppers.
- **Experimentation matrix**: every backend × every §10 system it claims — error, wall
  time, step count, fallback count — produced by a benchmark harness (under
  `scripts/benchmarks/` in the new package) and printed here in a future revision as the
  evidence table. This artefact sets and revises the §3 classifier thresholds and decides
  the implementation order, replacing any a-priori ranking.
- **Prototype registry already run** (scratchpad, 2026-08-31): six steppers × A1–A3
  produced the measured numbers quoted in §8–§10 — including two findings reading alone
  would not have caught (the CF4 misprint + composition order; the cond-guard blindness
  at Jordan points).

## 12. Recommendation and open items

**Architecture: (iii)** — own solver chained to unmodified CAMB — **for both engines,
over one shared core** (§7). O3 forces it; O2 confirms it.

**Implementation order** (each stage gated by its benchmarks):

1. η-grid segmented assembly + the stepper registry with `rk-adaptive`, `midpoint`,
   `magnus4`, `cf4` (the shared prerequisite; A1/A2/A8 as gates).
2. O2 front-end: sourced evolution, dense output, LOS hand-off; **A8 against real CAMB
   tables is the gate.**
3. O3 front-end: generic eikonal engine (wavefunction + density-matrix forms) + the
   patch chain; gated on the WS2 reduction stage (§5) — its absence blocks O3, not O2;
   B1/B3/B6 as gates.
4. Matrix-WKB rung (§9): promoted from prototype when the bake-off shows the Magnus
   baseline exceeding budget in the `k ≳ 0.1 Mpc⁻¹` band (A4 evidence); publishable
   independently (§2.2).
5. Deferred: piecewise-analytic (rung 4); emulation (only on benchmark failure); GPU/JAX
   path (only with GPU allocations); CF6/CFET high-order variants (coefficients on hand).

**Open items** (issues filed alongside this doc):

- `solve_modal` silent t=0 freeze when `can_use_modal` is bypassed (legacy; documented
  hazard until the legacy path is retired).
- WS2: `EikonalReduce.wl` + the multi-representation JSON sections (§5), linked to #209
  and H4's conventions spec.
- η-grid segmented assembly as the WS3 prerequisite (#492).
- Matrix-WKB implementation tracker (bake-off-gated).
- O3 eikonal engine + patch chain tracker (O4-mix #511 as second client).
- `tidalcosmo/solver/README.md` revision to reflect this design.
- Re-verify remaining **[survey]** tags against primary sources as they are first used
  in implementation (LJL bounds; the Boltzmann-code extension-mechanism details).

## References appendix — local TeX anchors

All under `literature/` unless noted. Line anchors are to the named files as of
2026-08-31.

| source | file : lines | used for |
| --- | --- | --- |
| Magnus review (Blanes et al., arXiv:0810.5488) | `section5.tex:818-838` (GL4), `:797-878` (GL6), `:2141-2158` (CF4 — misprinted coefficients, see §8), `:1945-1962` (modified Magnus); `section2.tex:1108-1120` (convergence); `section6.tex:404-408` (Moan error); `section4.tex:512-520` (adiabatic-picture convergence) | §8, §9 |
| Ioannisian & Smirnov (arXiv:0803.1967) | `rmagnus.tex:351-425` (interaction picture), `:1067-1244` (adiabatic Magnus, error scaling) | §9 |
| CLASS II (arXiv:1104.2933) | `approximations2.tex:884-897` (TCA timing table), `:901-1000` (UFA), `:1190-1270` (RSA) | §8 |
| oscode (arXiv:1906.01421) | `main.tex:607` (coupled extension "more difficult than anticipated"); App. on error estimators | §2.2, §9 |
| riccati (arXiv:2212.06924) | `main-expanded.tex:890` (step-type proposition), `:788-920` (adaptive selection) | §9 |
| Haddadin & Handley (arXiv:1809.11095) | `Paper.tex:137-147` (frictionless scalar form required), `:371-372` (coupled case open) | §8 |
| Cembranos et al. (arXiv:2302.08186) | `conversionDef4arxiv.tex:270-274` (eq. 20), `:291-330` (eqs. 24-27), `:352-381` (eq. 34, eq. 37 parameters), `:443-485` (density matrix) | §4, §10 |
| Ejlli (arXiv:2004.02714) | SVEA vs exact | §2.3 |
| Dolgov & Ejlli (arXiv:1211.0500) | §5 (density-matrix necessity) | §2.3, §4 |
| Domcke & Garcia-Cely (arXiv:2006.01161) | `text.tex:207` (the `𝓕` rate) | §4, §10 |
| Seljak & Zaldarriaga (astro-ph/9603033) | `los.tex:355-362` (tensor equation, dropped stress source) | §4 |
| TorC (arXiv:2507.09228) | `paper_Qtorsion.tex:360, 389` (background-only CAMB modification) | §6 |
| DISCO-EB (arXiv:2311.03291) | `main.tex:251` (declines timings); architecture §2 | §2.1, §6 |
| CAMB source | `equations.f90:225, 488` (`RungeKuttaDP45`), `derivst` (tensor block; fetched copy, session scratchpad) | §1, §2.1 |
| H3 additions (18 papers, 2026-08-31) | see `docs/references.md` §"Solver design study" for the per-paper relevance lines | §2 |
| legacy evidence | `docs/MODAL_TEMPLATE_CACHE_RETROSPECTIVE.md:24-31, 242`; `docs/meetings/2026-05-22_supervisor.md:19-38`; `tidal/solver/modal.py:405-448, 1526-1549, 3612, 3750, 3775` | §1 |
| non-arXiv | Lorenz–Jahnke–Lubich BIT 45 (2005) 91; Jahnke–Lubich Numer. Math. 94 (2003) 289; Blanes–Moan APNUM 56 (2006) 1519; Iserles BIT 42 (2002) 561; Raffelt–Stodolsky PRD 37 (1988) 1237 | §2.2, §5, §9 |
