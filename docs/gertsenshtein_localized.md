# Localized B-field Scattering: Thin-Lens vs Thick-Lens Physics

## Overview

When a graviton wavepacket passes through a spatially localized magnetic field region, the
conversion probability depends on the **ratio of the field region size R to the probe wavepacket
width σ**. Two distinct physical regimes emerge:

| Regime | Condition | Conversion formula | Analogy |
|--------|-----------|-------------------|---------|
| **Thin-lens** | R ≪ σ | P = sin²(κ/2 × ∫B dz) | Boccaletti (1970), Born approx. |
| **Thick-lens** | R ≫ σ | Local Rabi oscillation | WKB / adiabatic |
| Transitional | R ≈ σ | Neither formula exact | Numerical required |

TIDAL numerical validation (48-point 2D sweep, kappa=1, k=2, σ=5):

| R/σ | R | Mean err vs Boccaletti | Max err | Regime |
|-----|---|----------------------|---------|--------|
| ≤1 | 2–5 | < 0.001 | 0.001 | Thin-lens: **Boccaletti exact** |
| 1.4 | 7.2 | 0.004 | 0.028 | Transitional |
| ≥2 | 10–15 | 0.10–0.34 | 0.51–0.90 | Thick-lens: Boccaletti fails |

---

## Boccaletti (1970) Thin-Lens Formula

### Derivation and assumptions

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

**Physical picture**: In the thin-lens limit (R ≪ σ), the graviton wavepacket is much wider
than the field region. From the wavepacket's perspective, the entire B-field acts as an
instantaneous "phase kick." The conversion amplitude is:

```
a₂ ~ (κ/2 × ∫B dz) × h₇
```

The integral ∫B dz is the **DC Fourier component** of B₀(z), representing coherent,
phase-matched conversion across the entire field region. This is the Born approximation
(Domcke, Garcia-Cely & Lee 2025, arXiv:2507.16609, Section 4.2).

**Validity condition**: R ≪ σ (B-field region much narrower than probe coherence length).
Equivalently: the interaction is "sudden" — the field switches on and off before the
wavepacket can evolve within it.

**Confirmed numerically** (TIDAL, defaults κ=1, Bpeak=0.1, R=5=σ):
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

## Thick-Lens / Adiabatic Regime

### Physical picture

When R ≫ σ, the graviton wavepacket is much narrower than the magnetized region. The
wavepacket propagates through the field, experiencing a **locally quasi-uniform** environment.
At each position z, the local field amplitude B(z) drives local conversion according to the
Rabi formula for a uniform field.

The conversion is no longer determined by the global integral ∫B dz, but by the **local
field value** as the wavepacket propagates through it. This regime corresponds to the
**adiabatic / WKB approximation** (Raffelt & Stodolsky 1988; Domcke, Garcia-Cely & Lee 2025,
arXiv:2507.16609, Sections 5.1 and Appendix E).

**WKB validity condition** (from Domcke et al. 2025, arXiv:2507.16609):
```
|B₀'(z)/B₀(z)| ≪ ω     (field varies slowly on scale of wavelength)
|(μ²(z))'| ≪ ω³         (effective mass varies slowly)
```
For our Gaussian field B₀ = Bpeak × exp(-z²/(2R²)):
```
|B₀'/B₀| = |z|/R²  →  max at |z| = R:  1/R
```
WKB valid when 1/R ≪ ω = k = 2.0, i.e., R ≫ 0.5. Our domain always satisfies this.

**Note**: Both regimes (R ≪ σ and R ≫ σ) can satisfy the WKB condition — they are
distinguished by the probe coherence length σ, not the field variation scale.

### Local Rabi oscillation

For R ≫ σ with peak field Bpeak, the wavepacket of width σ sweeps through the peak field
region. Using the uniform-field formula with effective length ≈ σ:
```
P_local ≈ sin²(κ × B_eff × σ_eff / 2)
```

For Bpeak=0.174, R=15, σ=5 (R/σ = 3): the wavepacket sees B ≈ Bpeak near z=0 for a
duration ~σ. This gives P_local ≈ sin²(κ×Bpeak×σ/2) ≈ sin²(0.435) ≈ 0.18. In practice,
the wavepacket spends extended time inside the wide Gaussian field, accumulating conversion
over the full traversal — giving the numerically observed P ≈ 0.92.

A full treatment requires numerical simulation (as done here) or the adiabatic integral
formula of Domcke et al. (2025, arXiv:2507.16609):
```
A_h^T ~ ∫₀ᴸ dz' B₀(z') exp(i ∫₀^{z'} μ²(z'')/2ω dz'')
```

### Non-resonant thick-lens formula (Berlin et al. 2024)

For the analogous case of axion-photon mixing (Berlin et al. 2024, arXiv:2405.08865),
the free-space non-resonant conversion probability is:
```
P_{a→γ}^free = (ω g²/4k) |∫ dy' e^{iω|y-y'|} e^{ik_a y'} B₀(y')|²
```
This is the Fourier transform of B₀ at the momentum transfer — identical to the
graviton-photon Born approximation. In the forward scattering (k_a = k_γ) limit, the
momentum transfer vanishes and the formula reduces to the Boccaletti ∫B dz result.

**Key insight**: Resonant conversion (axion-photon) depends only on the **local** B at
the resonance point (where the plasma frequency equals the axion mass), not on the integral.
This is the extreme thick-lens limit where the probe coherence length → 0.

---

## Coherence Length and Regime Transition

The transition between regimes is governed by the probe **coherence length** L_coh, which
characterizes the spatial extent over which the wavepacket maintains phase coherence:

```
L_coh ~ σ  (Gaussian wavepacket width in position space)
       = 1/(Δk)  (inverse momentum width in Fourier space)
```

For our IC: L_coh = σ = 5 (in natural units where c=1).

**Regime boundaries**:
- Thin-lens (Boccaletti valid): R ≪ L_coh = σ
- Thick-lens (WKB/Rabi): R ≫ L_coh = σ
- Transitional: R ≈ L_coh

The condition R ≪ L_coh is equivalent to requiring that the B-field region be a "thin"
perturbation from the probe's perspective — the entire field region is traversed in a time
much shorter than the probe's oscillation period.

**Comparison with oscillation length**:
The oscillation length L_osc = 2π/(κ × Bpeak) characterizes the Rabi cycle. For our
parameters: L_osc = 2π/(1.0 × 0.1) ≈ 62.8. Since L_osc ≫ σ and L_osc ≫ R in all our
sweep cases, the thin-lens condition is independent of L_osc — it depends purely on R vs σ.

### Analogy with quantum optics

The thin/thick-lens distinction is the graviton-photon version of the **thin/thick target
approximation** in quantum optics and nuclear physics:

| Concept | Thin target | Thick target |
|---------|-------------|--------------|
| Field/barrier | Much narrower than wavepacket | Much wider than wavepacket |
| Conversion | Single integrated phase kick | Continuous accumulation |
| Formula | Born approximation / ∫V(z)dz | WKB / local coupling |
| Optical analogy | Thin lens (immediate refraction) | Thick lens (propagation inside) |
| Neutrino analogy | Thin matter slab (MSW inactive) | Solar interior (MSW resonance) |

The MSW (Mikheyev-Smirnov-Wolfenstein) resonance in solar neutrino oscillations is the
thick-lens analog: neutrinos propagate through a density gradient and undergo adiabatic
level crossing, which is only possible in the thick-target regime (density varies over many
oscillation lengths). Raffelt & Stodolsky (1988) explicitly drew this analogy for the
graviton-photon case.

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
(zero momentum transfer). In the backward direction, different parts of the field contribute
with rapidly oscillating phases (momentum transfer 2ω), causing destructive interference.
This is why the Boccaletti integral formula captures only forward (transmitted) conversion.

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
= [Bh]_{-∞}^{+∞} - ∫B ∂_x(h)).

### Simulation setup

For the thin-lens regime (R ≤ σ = 5):
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
| [Boccaletti et al., Nuovo Cimento 70B:129](https://link.springer.com/article/10.1007/BF02710177) | 1970 | **Finite-region formula** P=sin²(κ/2 × ∫B dz); thin-lens limit |
| [Raffelt & Stodolsky, PRD 37:1237](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.37.1237) | 1988 | 2×2 mixing matrix; MSW resonance analogy; **adiabatic / thick-lens framework** |
| [Ejlli, JHEP 2020:029](https://arxiv.org/abs/2004.02714) | 2020 | **Exact (non-perturbative) solution** in uniform B — generalizes both limits |
| [Dandoy, Lella et al., arXiv:2406.17853](https://arxiv.org/abs/2406.17853) | 2024 | 4-component mixing; canonical normalization; confirms P=sin²(κB₀D/2) |
| [Berlin et al., arXiv:2405.08865](https://arxiv.org/abs/2405.08865) | 2024 | Axion-photon mixing numerics; WKB validity; free-space vs resonant regimes (analogous) |
| [Domcke, Garcia-Cely & Lee, arXiv:2507.16609](https://arxiv.org/abs/2507.16609) | 2025 | **GW scattering on 3D B-fields**; Born approx.; forward vs backward scattering; WKB for inhomogeneous B |

### Terminology mapping

| This document | Domcke et al. (2025) | Berlin et al. (2024) | Raffelt & Stodolsky (1988) |
|--------------|---------------------|---------------------|---------------------------|
| Thin-lens | Born approx., far-field | Free-space (non-resonant) | Diabatic / sudden |
| Thick-lens | WKB / adiabatic | Resonant conversion | Adiabatic / MSW-like |
| Boccaletti integral | DC Fourier component ∫B₀dx | ∫B₀ e^{iq·r} dr at q=0 | ∫B · (phase factor) |
| Transition at R≈σ | Born-adiabatic crossover | Coherence length condition | Adiabaticity parameter Ψ |

---

## Open Questions

1. **Exact thick-lens formula**: What is the conversion probability for R ≫ σ with Gaussian
   B-field? The Domcke et al. (2025) adiabatic formula requires a resonance (μ=0 crossing);
   without plasma, B₀ variation alone doesn't create a resonance for massless fields.

2. **Wavepacket dispersion**: For R ≫ σ, the wavepacket disperses as it traverses the field.
   Position-dependent coupling can cause group velocity modification. Not captured by the
   simple local-Rabi picture.

3. **Cross-regime applications**: Real magnetic field configurations (e.g., magnetar dipole
   B ∝ 1/r³, galactic B ∝ random) involve both thin-lens and thick-lens elements at different
   scales. The TIDAL framework can handle arbitrary B(z) profiles via the background field
   mechanism — the local equations of motion are exact regardless of regime.

4. **σ dependence**: The thin-lens validation used σ=5. Verifying the formula for σ ≫ R at
   various σ values would strengthen the regime characterization. Expect the Boccaletti formula
   to hold for any σ ≫ R (the IC width should not appear in the conversion probability).
