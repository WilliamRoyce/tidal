# O4 foundation: cosmic birefringence from parity-odd torsion–EM couplings

**Status:** research notes from the 2026-08-29 planning session, preserved here because they
were otherwise held only in a session plan file outside version control. O4 is the last rung
of the observable ladder (`docs/COSMOLOGY_PROGRAM.md`); these notes are its physics
foundation and its literature map. The open derivation work is **GH #499**.

> **Nothing here is a result.** It is a scoping map with the traps marked. Where a claim was
> reasoned rather than read out of a source, it says so.

---

## 1. The mechanism is already enumerated in this repo

`research/lagrangian_enumeration/general_quadratic_lagrangian.tex` §`sec:chern_simons`
(line 529) gives a 3-parameter Chern–Simons torsion–EM family:

```
L_CS = ε^{μνρσ} V_μ A_ν F_ρσ ,    V_μ ∈ { T_μ (trace torsion), S_μ (axial torsion), … }
```

This is the **Carroll–Field–Jackiw** structure — the canonical source of cosmic
birefringence. It is the same structure the axion case produces after integrating
`(g/4) χ F F̃` by parts, with `V_μ → (g/2) ∂_μ χ`.

`enumeration_classified.json` classifies the blocks:

| block | terms | `ghost_status` | `topological` | note |
|---|---|---|---|---|
| `CS: eps Tor A F` | 6 | safe | **False** | physical, not a total derivative |
| `CS: eps Tvec A F` | 1 | safe | **False** | |
| `eps F x F` | 1 | safe | **True** | the bare `F·F̃` θ-term — **no rotation on its own** |
| `eps T x DF` / `eps F x DT` | 6 each | borderline / safe | False | derivative couplings, see §3 |

**Why this fits FRW.** Homogeneity and isotropy permit a torsion background with only a
*temporal* vector component, `V_μ = (V(η), 0, 0, 0)`. Substituted into `L_CS` this is
precisely the birefringence term, and the observable is a line-of-sight integral over cosmic
history: `β = ∫_{η_LSS}^{η_0} V_0(η) dη` (up to normalization). A perturbed `δV(x, η)`
extends this to *anisotropic* birefringence — a rotation field `α(n̂)` with its own power
spectrum.

**⚠ No example implements these couplings.** `examples/torsion_gertsenshtein/theory_parity_odd.toml`
explicitly defers `cs1`–`cs3` ("require bare `A_μ`, special handling"). The bare-`A_μ` gauge
handling is an unsolved pipeline problem, not a configuration step. **GH #499.**

## 2. Conformal invariance — a solver fast path, NOT a scoping shortcut

Written with the Levi-Civita *tensor* `ε^{μνρσ} = ε̃^{μνρσ}/√−g`:

```
√−g · ε^{μνρσ} V_μ A_ν F_ρσ  =  ε̃^{μνρσ} V_μ A_ν F_ρσ
```

— all indices down, no metric anywhere. The term is *metric-independent*. And the Maxwell
action is conformally invariant in exactly D=4 (`√−g → Ω⁴√−g` cancels `g^{μρ}g^{νσ} → Ω⁻⁴`).
So on FRW in conformal time the photon + CS sector is **identical to flat space**; the scale
factor never enters. This is why the axion result `β = (g/2)Δχ` is cosmology-independent.

**How this is used (user directive D7):** as a **detection rule the solver computes from the
spec**, alongside the existing Modal → IDA → CVODE cascade — never as a reason to specialize
the engine to birefringence. The general FRW path is always the fallback:

```text
conformal-weight analysis of derived spec
  ├─ sector conformally invariant in D=4? → a(η) provably drops out; evolve on the flat
  │                                         slice. EXACT, not an approximation.
  └─ otherwise                            → general FRW: a(η) as time-dependent coefficient
```

Does **not** apply to gravitons/tensor modes (they pick up `a''/a`), massive torsion modes,
or the derivative-coupling blocks (`eps T x DF`, `eps F x DT`), which carry explicit metric
contractions. Those exercise the general path — which is why it must be first-class.
Useful side-effect: the conformal case is a cheap exact test the general path must reproduce
to machine precision.

*This is reasoned from standard field theory, not read from a repo source. WS0/WS2 should
verify it symbolically, and derive the detection rule from conformal weights rather than
hardcoding it.*

## 3. Frequency scaling — the decisive question, and it has no single answer

The scaling is fixed by the mass dimension of the effective photon-sector operator —
equivalently, by how many derivatives act on `F`.

| Lagrangian structure | eff. operator dim | `β(ν)` | effect |
|---|---|---|---|
| `(g/4) χ F F̃` (axion/CS) | 3 (CPT-odd) | `ν⁰` | pure rotation, E↔B |
| `(k_AF)_μ` constant (CFJ / SME `d=3`) | 3 | `ν⁰` | pure rotation |
| `θ(T) F F̃`, `θ` built from `T²` (Itin–Hehl) | 3 | `ν⁰` | pure rotation |
| `ξ T F ∂F̃` (one extra `∂`) | 5-equivalent | `ν²` | rotation, dispersive |
| `k_F`-type `χ^{abcd} F F` (`T²F²`, Preuss) | 4 (CPT-even) | — | **linear** birefringence → E→V, ellipticity |
| skewon part of `χ^{abcd}` | — | — | **dichroism**, absorption — not rotation |
| plasma / PMF Faraday rotation | — | `ν⁻²` | rotation |

Kostelecký & Mewes (arXiv:0905.0031): for SME operators of mass dimension `d`, CPT-odd
polarization change scales as `ω^{d−3}` — **only `d = 3` is frequency-independent**.

**Why it is decisive.** The measured birefringence is frequency-independent across 23–353 GHz
(Eskilt & Komatsu, arXiv:2205.13962). A `ν²` operator is therefore **already excluded as an
explanation** of the signal and can only be bounded; a `ν⁰` operator could explain it. This
also determines whether an off-the-shelf `β` posterior is usable or a multi-frequency
likelihood must be built.

**But there is no single answer for "the torsion coupling".** TIDAL's own theory documents
derive none of this: `general_quadratic_lagrangian.tex` has no frequency-scaling discussion,
and its CS section states only that `χ^CS_2` sources magnetic helicity `A·B`. The Das et al.
`ξ₁ ∇T×F̃` structure gives `α = 2ξ₁p²T₁t` (`ν²`); an induced-axion `θFF̃` gives `ν⁰`. These
are *different sectors* in the enumeration's language (`ζ̃₁₋₆` vs `d₁₈`/`χ^CS`). **Scope it as
a per-operator derivation task — GH #499.**

## 4. ⚠ Terminology trap — two different things are called "birefringence"

- **Hehl/Obukhov school:** "birefringence" means a **double light cone** — two phase
  velocities for two *linear* polarizations (the quartic Fresnel surface failing to
  factorize). That is *linear* birefringence → ellipticity, **E→V**.
- **CMB literature:** "cosmic birefringence" means **optical activity** — rotation of the
  polarization plane, `Q ± iU → e^{∓2iβ}(Q ± iU)` → **EB, TB**.

These are different observables. **`gr-qc/0305049` (Rubilar–Obukhov–Hehl, local) derives the
*linear* kind — citing it as support for a CMB `β` would be wrong.** Same caution for
Preuss et al. `gr-qc/0507071` (PGT gravity-induced birefringence, constrained with magnetic
white-dwarf polarimetry, not CMB).

**The premetric classification that settles which is which** (Hehl–Obukhov, local at
`literature/gr-qc_0001010/`): any `T²F²`-type coupling modifies the constitutive tensor
`χ^{abcd}` (36 components), which decomposes as

```
χ  =  χ_(principal, 20)  ⊕  χ_(skewon, 15)  ⊕  χ_(axion, 1)
```

- **axion part (1)** → optical activity → **rotation, EB/TB, frequency-independent** ✅ what we want
- **principal part (20)** → light-cone deformation → linear birefringence, E→V
- **skewon part (15)** → **dichroism / absorption**, damping one helicity — not rotation

## 5. Structural notes on torsion–photon couplings

- **`T^μ A_μ`-type (torsion vector × potential):** **not gauge invariant**; the price is a
  photon mass. Standard Riemann–Cartan lore. The known escape (arXiv:1701.07132) routes the
  torsion *trace* vector into a Weyl connection, where it cancels identically out of `F_μν` —
  so **the trace vector cannot give photon birefringence at all** without breaking gauge
  invariance.
- **Linear-in-`F` couplings** (`~F T²`) *modify Maxwell's equations* — they source `F` rather
  than rotating it. The bilinear `~F²T²` class leaves Maxwell intact and only modifies `χ`,
  which is where the axion piece lives. The gauge-invariant `F^{μν} R̃_{[μν]}`
  (arXiv:2507.02362, local) is a *mixing* term → photon↔torsion oscillation
  (Gertsenshtein-like), not a helicity-dependent phase.
- **Condition for frequency-independent rotation** *(reasoned, not a cited result — verify):*
  the effective term must be `(k_AF)_μ ε^{μνρσ} A_ν F_ρσ` with `(k_AF)_μ` a **closed
  one-form**, else gauge invariance fails. On homogeneous FRW with axial torsion
  `A^tors_μ = (A_0(η),0,0,0)`, `ε^{μνρσ}∂_ν A^tors_μ F_ρσ = 0` identically, so
  `(k_AF)_μ = c A^tors_μ` is admissible at background level, giving
  `β = (c/2)∫A^tors_0(η) a(η) dη`. For *perturbed* torsion this fails unless
  `A^tors_μ = ∂_μ θ`.

## 6. The observable

**Exact rotation formulae** (uniform rotation by `β`; ΛCDM has `C_ℓ^{EB} = C_ℓ^{TB} = 0`):

```
C_ℓ^{EE,o} = C_ℓ^{EE} cos²(2β) + C_ℓ^{BB} sin²(2β) − C_ℓ^{EB} sin(4β)
C_ℓ^{BB,o} = C_ℓ^{EE} sin²(2β) + C_ℓ^{BB} cos²(2β) + C_ℓ^{EB} sin(4β)
C_ℓ^{EB,o} = ½(C_ℓ^{EE} − C_ℓ^{BB}) sin(4β) + C_ℓ^{EB} cos(4β)
C_ℓ^{TB,o} = C_ℓ^{TE} sin(2β) + C_ℓ^{TB} cos(2β)
```

Small-`β`: `C_ℓ^{EB,o} ≃ 2β(C_ℓ^{EE} − C_ℓ^{BB})`, `C_ℓ^{TB,o} ≃ 2β C_ℓ^{TE}`. EB carries
essentially all the statistical weight.

**When post-processing rotation is exact — and when it is not.** Exact iff `β` is a single
constant, isotropic, and frequency-independent. It must go *inside* the line-of-sight
integral when `β` is time-dependent over recombination (Murai et al., arXiv:2209.07804),
tomographic (recombination vs reionization), anisotropic (arXiv:2404.13771), scale-dependent,
or **frequency-dependent** — in which case it must be applied per frequency channel *before*
component separation, interacting with foreground cleaning. **That last case is likely ours
if the operator is `ν²`; scope it early.** Public engine for the hard cases: `class_rot`
(arXiv:2111.14199).

**Anisotropic birefringence** — a rotation field `α(n̂)` with `C_L^{αα} = 2πA_CB/[L(L+1)]`,
reconstructed with EB quadratic estimators (Kamionkowski arXiv:0810.1286). **Crucially it
does NOT suffer the miscalibration degeneracy** (a constant instrument angle contributes only
to `L=0`), and public likelihoods exist — SPT-3G BB lite (arXiv:2510.07928),
combined limits arXiv:2504.13154.

## 7. ⚠ The miscalibration degeneracy — the central experimental difficulty

An instrument's absolute polarization-angle offset `α` rotates *every* photon it measures,
CMB and Galactic foreground alike. Cosmic birefringence rotates only the CMB (foreground
photons travel ≲ kpc). From CMB spectra alone `C_ℓ^{EB,o} ≃ 2(α+β)(C_ℓ^{EE} − C_ℓ^{BB})` —
**only the sum is measurable.** Planck's absolute-angle systematic (~0.28°) is the same size
as the claimed signal.

**Minami–Komatsu degeneracy-breaking** (arXiv:1904.12440, 2006.15982): `β` multiplies **only
the CMB** spectra while `α_i` multiplies the **total observed** (CMB + foreground). Because
the foreground amplitude is strongly frequency-dependent while the CMB's is not, and `β` is
common across channels while `α_i` is per-channel, fitting `{β, α_i, A_fg}` jointly is
well-posed. **Achilles heel:** intrinsic dust `EB` is not zero (Clark et al.,
arXiv:2105.00120) — the source of the ~0.28° systematic in map-space analyses.

**Practical implication:** an off-the-shelf `β` likelihood is only usable if our prediction is
a single constant `β`. Frequency- or ℓ-dependence changes the `α_i` nuisance structure, so a
published `P(β)` cannot simply be imported.

## 8. Measurements and likelihood availability

| Analysis | `β` | significance |
|---|---|---|
| Minami & Komatsu 2020 (arXiv:2011.11254), Planck PR3 | 0.35° ± 0.14° | 2.4σ |
| Diego-Palazuelos et al. 2022 (arXiv:2201.07682), PR4 | 0.30° ± 0.11° | ~2.7σ |
| Eskilt & Komatsu 2022 (arXiv:2205.13962), WMAP+Planck 23–353 GHz | — | 3.6σ, **no ν-dependence** |
| Diego-Palazuelos & Komatsu 2025, ACT DR6 | 0.215° ± 0.074° | 2.9σ |
| **Eskilt 2026 (arXiv:2608.06480), joint ACT+Planck** | **0.277° ± 0.057°** | **4.8σ** (3.5σ dust-robust) |

Authors caution that unresolved systematics must be understood before strong conclusions.
Note SPT and older BICEP **self-calibrate by forcing `C_ℓ^{EB}=0`**, which destroys the
isotropic signal by construction — they contribute anisotropic limits, not `β`.

**There is no drop-in Cobaya birefringence likelihood.** Escalation, cheapest first:

1. **Gaussian prior on a published `β`** — legitimate *only* for a constant,
   frequency-independent prediction.
2. **Fork `LilleJohs/cosmic-birefringence-planck-act` (MIT)** — its `.npz` binned spectra +
   covariances are exactly the `{C_b^{E_iB_j,o}, Cov}` needed; replace the model function
   with our `β(ν, ℓ)`. Best value for effort if we need frequency or ℓ dependence.
3. **Anisotropic route** — SPT-3G BB lite; cleanest statistically (no `α` degeneracy).
4. Full multi-frequency EB likelihood from NPIPE/ACT DR6 maps — months.

**Forecasts:** LiteBIRD `σ(β) ~ 0.02°` including systematics (arXiv:2503.22322); SO ~0.04°;
CMB-S4/PICO statistics-only ~0.0017° (arXiv:1904.07855) — all calibration-limited in practice.

## 9. Literature gaps this project could fill

1. **A modern update of Das, Mohanty & Prasanna (arXiv:0908.0629)** — the *only* existing
   torsion→CMB-birefringence paper (2009, `ν²`, pre-Minami–Komatsu data, no treatment of the
   `α` degeneracy). It appears **never to be cited by the modern CMB birefringence
   literature**.
2. **An explicit Barbero–Immirzi / Nieh–Yan `γ · F F̃` calculation.** Calcagni & Mercuri give
   the pseudoscalar; nobody appears to have coupled it to `F F̃` and computed `β`. Note the
   Nieh–Yan literature is otherwise uniformly about *gravitational-wave* birefringence, not
   photon.
3. **The general "PGT Lagrangian → `β(η, ν, ℓ)`" mapping** — no such framework exists in any
   Boltzmann ecosystem (EFTCAMB, MGCAMB, hi_class all modify the *gravity* sector; none
   touches the photon constitutive relation).

## 10. What is already local

`literature/gr-qc_0307063/` (**Itin & Hehl** — torsion-induced *axion* `θ(T) ∝ ε T T`, giving
frequency-independent optical activity; the best candidate for a signal-explaining channel;
verified present 2026-08-29), `gr-qc_0305049/` (Rubilar–Obukhov–Hehl — the *linear* kind, see
§4 trap), `gr-qc_0001010/` (Hehl–Obukhov premetric decomposition), `0804.4011/` (Kruglov,
axial-torsion–photon mixing), `2410.01355/` (Trukhanova–Obukhov, spin-torsion), `2507.02362/`
(Bahamonde et al., `F^{μν}R̃_{[μν]}`). `manuscript/references.bib` already carries these under
a "Torsion-EM Coupling Foundations" heading.

**Absent locally and needed:** the entire observational/statistical side (Minami, Komatsu,
Eskilt, Diego-Palazuelos, Planck/ACT/SPT/LiteBIRD). H5 fetched the O4 set — see
`docs/references.md` §Cosmology Program.

---

## See also

- `docs/COSMOLOGY_PROGRAM.md` — the operational record; O4 sits on the observable ladder.
- `docs/cosmology/primer.md` — how a CMB pipeline works and where TIDAL plugs in.
- `docs/cosmology/spectator_route.md` — why propagation/conversion effects are reachable in
  the spectator limit and gravitational sourcing is not.
- **GH #499** — the open derivation work (per-operator frequency scaling; CS implementation).
