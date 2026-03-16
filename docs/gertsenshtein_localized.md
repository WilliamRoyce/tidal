# Localized B-field Scattering: Boccaletti Formula and Phase Matching

## Overview

When a graviton wavepacket passes through a spatially localized magnetic field region, the
conversion probability is given by the **Boccaletti (1970) formula**:

```
P(graviton → photon) = sin²(κ/2 × ∫_{-∞}^{∞} B₀(z) dz)
```

For a Gaussian profile B₀(z) = Bpeak × exp(-z²/(2R²)), the integral gives
∫B₀ dz = Bpeak × R × √(2π), so:

```
P = sin²(κ × Bpeak × R × √(π/2))
```

**TIDAL numerical validation (48-point 2D sweep, κ=1, k=2, σ=5)**:

| R/σ | R | Mean err | Max err | Status |
|-----|---|----------|---------|--------|
| 0.40 | 2.0 | 0.00004 | 0.00008 | Boccaletti exact |
| 0.92 | 4.6 | 0.00029 | 0.00045 | Boccaletti exact |
| 1.44 | 7.2 | 0.00078 | 0.00172 | Boccaletti exact |
| 1.96 | 9.8 | 0.00113 | 0.00211 | Boccaletti exact |
| 2.48 | 12.4 | 0.00091 | 0.00144 | Boccaletti exact |
| 3.00 | 15.0 | 0.00058 | 0.00160 | Boccaletti exact |

The formula is confirmed to < 0.003 error across all 48 parameter combinations, including
R = 3σ (wavepacket width much narrower than field region). This is explained by perfect
phase matching in massless vacuum conversion (see [§ Phase matching](#phase-matching)).

---

## Boccaletti (1970) Formula: Derivation and Exact Validity

### Derivation

Boccaletti et al. (1970) extended Gertsenshtein's (1962) original result to finite-region
boundary conditions. For a localized background B₀(z), the conversion probability is:

```
P(graviton → photon) = sin²(κ/2 × ∫_{-∞}^{∞} B₀(z) dz)
```

For a Gaussian profile B₀(z) = Bpeak × exp(-z²/(2R²)):
```
∫_{-∞}^{∞} B₀(z) dz = Bpeak × R × √(2π)

P = sin²(κ × Bpeak × R × √(π/2))
```

**Confirmed numerically** (TIDAL, defaults κ=1, Bpeak=0.1, R=5.0, σ=5.0):
```
P_numerical = 0.3436  (at t=80.4)
P_Boccaletti = sin²(1.0 × 0.1 × 5.0 × √(π/2)) = sin²(0.627) = 0.3432
Agreement: 0.04%
```

### Connection to uniform-B formula

For a uniform field B₀ over length D, the Boccaletti integral gives:
```
∫₀ᴰ B₀ dz = B₀ D
P = sin²(κ B₀ D / 2)
```
This recovers the standard Gertsenshtein formula. Boccaletti's integral formula generalizes
this to arbitrary spatial profiles.

---

## Phase Matching: Why Boccaletti is Exact for Massless Vacuum Conversion {#phase-matching}

### Perfect phase matching

For massless graviton-photon conversion in vacuum (no plasma), both modes travel at the same
speed (c = 1): their dispersion relations are identical (k_h = k_γ = ω). This means:

```
Δk = k_γ - k_h = 0   (zero momentum transfer for forward scattering)
```

The conversion amplitude in the Born approximation is:
```
A ~ (κ/2) ∫ B₀(z) e^{i Δk z} dz = (κ/2) ∫ B₀(z) dz   (since Δk = 0)
```

Since Δk = 0 exactly, the integrand has **no oscillating phase factor** — all positions z
contribute coherently, regardless of whether R ≪ σ or R ≫ σ. The Boccaletti formula is
therefore **exact** for any spatial profile of B₀(z), not just a thin-lens approximation.

**Physical picture**: Unlike in quantum optics (where different wavelengths create phase
mismatch) or plasma-detuned conversion (where the photon has a mass), graviton and photon
in vacuum have identical dispersion. There is no coherence length limiting the effective
conversion region — the entire field extent contributes in phase.

### Born approximation validity

The formula is derived in the weak-coupling (Born) approximation. The small-angle condition
is satisfied when:
```
κ/2 × ∫B₀(z) dz ≪ 1   (weak conversion)
```
For our sweep range (Bpeak ∈ [0.02, 0.20], R ∈ [2, 15]), the argument ranges from 0.05 to
3.76 — this includes strong-field regimes near the first Rabi maximum, where the Born
approximation breaks down yet the formula still holds (it is the exact non-perturbative Rabi
formula sin²(θ/2), not just its small-angle limit).

### When does the Boccaletti formula fail?

The formula P = sin²(κ/2 × ∫B₀ dz) breaks down when:

1. **Plasma detuning**: If the photon acquires a plasma mass ω_p, the phase mismatch is
   Δk = ω_p²/(2ω), creating a coherence length L_coh = 2ω/ω_p². For R ≫ L_coh, contributions
   from different parts of the field cancel and the formula fails (Raffelt & Stodolsky 1988).

2. **Axion-photon mixing (non-zero axion mass)**: Same mechanism — the axion and photon travel
   at different speeds, and the effective integral picks up the Fourier component at q = Δk
   rather than q = 0.

3. **Resonant conversion**: When ω_p = m_axion (or analogously Δk = 0 locally at some point),
   the resonance condition creates a "thick-lens" Rabi behavior that is not captured by the
   Born approximation.

The TIDAL simulations in this document use massless fields without plasma, so none of these
apply and the Boccaletti formula holds exactly.

---

## Thin-Lens / Thick-Lens Distinction: When It Matters

The thin-lens / thick-lens terminology is used in the literature (and in the original
version of this document) to describe two regimes of localized scattering:

| Regime | Condition | Conversion formula | Context |
|--------|-----------|-------------------|---------|
| **Thin-lens** | R ≪ σ or R ≪ L_coh | P = sin²(κ/2 × ∫B dz) | Born approx., Boccaletti |
| **Thick-lens** | R ≫ L_coh | Local Rabi oscillation | WKB / adiabatic |
| Transitional | R ≈ L_coh | Neither formula exact | Numerical required |

where the relevant coherence length is **not the wavepacket width σ** but the phase-mismatch
coherence length L_coh = 1/|Δk|.

For **massless vacuum** conversion: Δk = 0, so L_coh → ∞. The "thin-lens" condition
R ≪ L_coh is always satisfied — the formula is valid for any R. This is why the TIDAL
sweep shows the Boccaletti formula holding even for R/σ = 3.

For **plasma-detuned** conversion (Raffelt & Stodolsky 1988): L_coh = 2ω/ω_p². When
R ≫ L_coh, the thick-lens regime applies. This is the relevant regime for axion dark matter
searches in magnetized plasmas (neutron stars, galactic clusters).

The MSW (Mikheyev-Smirnov-Wolfenstein) resonance in solar neutrino oscillations is the
extreme thick-lens case: a density gradient creates a slowly-varying phase mismatch that
passes through zero, enabling adiabatic level crossing. Raffelt & Stodolsky (1988) drew
this analogy explicitly for the graviton-photon case.

---

## Forward vs Backward Scattering

Domcke, Garcia-Cely & Lee (2025, arXiv:2507.16609) provide additional insight by treating
the problem as classical electromagnetic scattering:

**Transmitted (forward) wave** (Eq. 28):
```
P_{h→γ}^T = 4πG × B₀^{T2} × D²   (for uniform field, forward direction)
```
The amplitude is **proportional to domain size L** → coherent, resonant enhancement.
For a localized source: amplitude = ∫B₀(x)dx (the DC Fourier component = Boccaletti).

**Reflected (backward) wave** (Eq. 27):
```
amplitude ∝ ∫ B₀(x) e^{2iωx} dx   (Fourier transform of B₀ at wavenumber 2ω)
```
For large domains (ωL ≫ 1), reflected amplitude ∝ sin(ωL)/(ωL) → 0. Backward scattering
is suppressed relative to forward transmission for extended fields.

**Physical meaning**: In the forward direction, all parts of the field contribute in phase
(zero momentum transfer, Δk = 0). In the backward direction, different parts of the field
contribute with rapidly oscillating phases (momentum transfer 2ω), causing destructive
interference. This is why the Boccaletti integral formula captures only forward (transmitted)
conversion, and only the DC Fourier component q=0 enters.

---

## TIDAL Implementation Details

### theory_localized.toml

The localized B-field is implemented via the gauge potential:
```
Abar_y = -Bpeak * R * Sqrt[Pi/2] * Erf[z[] / (Sqrt[2]*R)]
∂_z(Abar_y) = -Bpeak * exp(-z²/(2R²)) = -B_x(z)  ✓
```

The derived equations of motion (from Euler-Lagrange via TIDAL pipeline) are:
```
d²t(h_7) = [-κ² B(x)²] h_7 + [-κ² B(x)] ∂_x(a_2) + ∂²_x(h_7)
d²t(a_2) = [B'(x)] h_7  + [B(x)] ∂_x(h_7) + ∂²_x(a_2)
```
where B(x) = Bpeak × exp(-x²/(2R²)) and B'(x) = -Bpeak × x × exp(-x²/(2R²))/R².

The photon source `B'(x) h_7 + B(x) ∂_x(h_7) = ∂_x(B(x) h_7)` is the **total derivative**
form, consistent with the Boccaletti integral formula (integration by parts gives ∫B'h + Bh'
= [Bh]_{-∞}^{+∞} - ∫B ∂_x(h))..

### Coefficient overflow fix

The Wolfram exporter sometimes serializes `Exp[-x²/R²]` as `1/E^(x²/R²)` in the
coefficient strings. After the `E^` → `exp` substitution, this becomes `1/exp(x²/R²)`,
which causes numpy to compute `exp(+x²/R²)` before dividing — triggering overflow warnings
for large |x| (e.g., at grid boundaries x = ±100 with small R = 2).

This is fixed in `tidal/symbolic/_eval_utils.py` (`_invert_exp_denominator`): `/exp(arg)` →
`*exp(-(arg))`, avoiding the positive-exponent computation entirely.

### Simulation setup

```bash
uv run tidal simulate examples/data/gertsenshtein_localized.json \
  --grid-shape 1024 "--bounds=-100:100" \
  --ic gaussian --ic-wavevector 2.0 --ic-amplitude 0.1 --ic-width 5.0 \
  "--ic-center=-50.0" --ic-component h_7 \
  --t-end 120.0 --param kappa=1.0 --param Bpeak=0.1 --param R=5.0
```

**Key choices**:
- Non-periodic domain (Neumann BCs): `Abar_y` is not periodic
- `t_end=120`: stops after one pass through B-field (wavepacket reaches x=0 at t≈50,
  exits field region at t≈65, hits right wall at t≈150)
- `--bounds="-100:100"` and `--ic-center=-50.0`: negative values require `=` syntax in CLI
- `--what conversion` (not `--what energy`): energy unreliable with position-dependent coefficients

### Validation scripts

- `examples/gertsenshtein/run_localized.sh` — single validation run (Bpeak=0.1, R=5)
- `examples/gertsenshtein/sweep_profile.sh` — 2D sweep (Bpeak ∈ [0.02,0.20] × R ∈ [2,15])

---

## Literature References

| Reference | Year | Key contribution |
|-----------|------|-----------------|
| [Gertsenshtein, JETP 14:84](https://www.scirp.org/reference/referencespapers?referenceid=1912195) | 1962 | Original graviton-photon conversion prediction (uniform B, plane wave) |
| [Boccaletti et al., Nuovo Cimento 70B:129](https://link.springer.com/article/10.1007/BF02710177) | 1970 | **Finite-region formula** P=sin²(κ/2 × ∫B dz); exact for massless vacuum conversion |
| [Raffelt & Stodolsky, PRD 37:1237](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.37.1237) | 1988 | 2×2 mixing matrix; MSW resonance analogy; **plasma-detuned thick-lens regime** |
| [Ejlli, JHEP 2020:029](https://arxiv.org/abs/2004.02714) | 2020 | **Exact (non-perturbative) solution** in uniform B — generalizes both limits |
| [Dandoy, Lella et al., arXiv:2406.17853](https://arxiv.org/abs/2406.17853) | 2024 | 4-component mixing; canonical normalization; confirms P=sin²(κB₀D/2) |
| [Berlin et al., arXiv:2405.08865](https://arxiv.org/abs/2405.08865) | 2024 | Axion-photon mixing numerics; WKB validity; free-space vs resonant regimes (analogous) |
| [Domcke, Garcia-Cely & Lee, arXiv:2507.16609](https://arxiv.org/abs/2507.16609) | 2025 | **GW scattering on 3D B-fields**; Born approx.; forward vs backward scattering; WKB for inhomogeneous B |

### Terminology mapping

| This document | Domcke et al. (2025) | Berlin et al. (2024) | Raffelt & Stodolsky (1988) |
|--------------|---------------------|---------------------|---------------------------|
| Exact (massless vacuum) | Born approx. = exact at Δk=0 | Free-space (non-resonant, q=0) | Diabatic / sudden |
| Plasma thick-lens | WKB / adiabatic | Resonant conversion | Adiabatic / MSW-like |
| Boccaletti integral | DC Fourier component ∫B₀dx | ∫B₀ e^{iq·r} dr at q=0 | ∫B · (phase factor) |
| Transition (plasma) | Born-adiabatic crossover | Coherence length condition | Adiabaticity parameter Ψ |

---

## Open Questions

1. **Plasma detuning + localized B**: For theory_plasma.toml (Proca mass for photon) combined
   with a Gaussian B(z), the thick-lens regime should manifest. The coherence length
   L_coh = 2ω/ω_p² provides the transition scale R ≈ L_coh. TIDAL can simulate this.

2. **Cross-regime applications**: Real magnetic field configurations (e.g., magnetar dipole
   B ∝ 1/r³, galactic B ∝ random) involve both thin-lens and thick-lens elements at different
   scales. The TIDAL framework can handle arbitrary B(z) profiles via the background field
   mechanism — the local equations of motion are exact regardless of regime.

3. **3D generalization**: The Domcke et al. (2025) framework treats full 3D B-field
   configurations. TIDAL currently handles 1+1D plane-wave reductions. Extending to 2+1D
   or 3+1D localized fields would require non-plane-wave IC and boundary conditions.
